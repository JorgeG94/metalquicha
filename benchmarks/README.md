# Benchmarks

Time this build on this machine, and say what changed since last time.

```bash
cmake --build build -j
python3 benchmarks/run_benchmarks.py --exe build/mqc --record   # first time
python3 benchmarks/run_benchmarks.py --exe build/mqc            # every time after
```

About twelve minutes on forty cores for a build in reasonable shape. It is a
performance regression suite and a "does this build perform the way it should
here" report, because those want the same measurements.

`--quick` skips the thread ladder and the fragmented cases, which is most of the
time. `--record` writes `baseline-<hostname>.json`; later runs compare against
it. Baselines are per host on purpose -- reference numbers checked into a
repository are meaningless on hardware unlike the machine that produced them.

## What it reports, and why each part is there

Every phase here exists because a plain wall-clock number was misleading during
a real investigation.

**The build configuration, first.** The BLAS choice alone moves a density
functional run by a factor of five, is set by one CMake cache variable, and
cannot be recovered from inside the program. A timing that does not say
`blas=Intel10_64lp_seq` is not reproducible. `--version` reports it and the
suite records it beside the numbers, warning if a baseline came from a
different build.

**Hartree-Fock beside every functional, as `xHF`.** A pure functional carries no
exact exchange and should cost *less* than Hartree-Fock. When this suite was
written it cost twenty-six times as much, and no absolute number would have said
so -- only the ratio to a reference sharing the same integrals. The ratio is
size-dependent, since the quadrature and the Fock build scale differently: the
same functional is 6.9x HF on ten waters and 4.2x on twenty. Compare like with
like.

**The same case at every thread count.** A stage that does not parallelise is
invisible at a single thread count; it just looks expensive. The
exchange-correlation quadrature was serial for the whole life of the code, and
this ladder is the report that would have shown it. Read the speedup column: a
run that stops improving has a serial stage in it somewhere.

**A fragmented case, serial and under MPI.** Fragment work pins itself to one
OpenMP thread and parallelises with MPI instead, so a change that helps a single
molecule can hurt the fragment path. One did: threaded BLAS is five times faster
on one molecule and thirty-one per cent slower on four ranks. Without this phase
that trade looks like a pure win.

Note that `fragmented/mpi4` being slower than `fragmented/serial` is not a fault.
The serial run puts every thread on one fragment at a time; four ranks put one
thread on each of four. They are different shapes of the same work, and which
wins depends on the system and the rank count. What the suite is watching is
whether either *changes*.

## How it decides something regressed

Against the measured spread, never a fixed percentage. A case that repeats to
0.3 per cent and moves by 3 has regressed; one that repeats to 16 per cent and
moves by 3 has not moved at all. One threshold for both is how a suite becomes
noise that people stop reading.

So every case is run more than once and the spread is stored with the number.
Where only one run was made, the suite says `1 run, no spread` rather than
`spread 0.0%`, and widens the band it judges by -- a single measurement is not a
measurement of reproducibility.

## Adding a case

Add it to `ENERGY_CASES` or `GRADIENT_CASES` in `run_benchmarks.py`. Decks are
generated into a temporary directory at run time rather than committed, since
they differ only in two fields. Geometries live in `geometries/`.

Keep this directory out of `validation/`: the CPU validation suite is generated,
and `tools/cpu_validation/gen_cpu_validation.py` deletes decks under
`validation/inputs/cpu/mqc/` that it did not write.
