# Performance

Measured facts to know before optimising or benchmarking anything. Read from `AGENTS.md`, which says when.

Four things worth knowing before optimising anything here, all of them measured.
`benchmarks/run_benchmarks.py` checks them on the machine it runs on; see `benchmarks/README.md`.

* **The exchange-correlation quadrature dominates a DFT run** -- 89 per cent of
  an LDA one and 98 of a meta-GGA one, against a Fock build under two per cent.
  It is threaded, and scales -- 32 s on one thread, 8 s on four, PBE on ten
  waters in 6-31G -- **provided the BLAS is sequential**. It calls `dgemm` from
  inside its own OpenMP threads, so a pthreads OpenBLAS left threaded (what a
  Debian or Ubuntu system library is) contends with them and the quadrature
  does not speed up at all: 33 s on four threads. `tools/run.sh` and
  `benchmarks/run_benchmarks.py` set `OPENBLAS_NUM_THREADS=1` for this; calling
  `mqc` directly does not. Held to one, PBE costs 3.4x Hartree-Fock on four
  threads.
* **Only xTB is pinned to one thread, and not for speed.** `needs_serial_execution`
  in `src/mqc_method_types.f90` names GFN1 and GFN2 and nothing else: threaded,
  tblite corrupts a result rather than failing. Every other method -- the whole
  cenzontle ab initio path -- keeps the threads the launcher gave it, on
  the unfragmented path and on the fragment workers alike, and a Fock build
  threads its own quartet loop. Clamping those was an old bug, not a policy:
  raising `OMP_NUM_THREADS` did nothing and the cause was nowhere near the
  command line. The measured "threaded MKL competes with the ranks" result --
  five times faster on one molecule, thirty-one per cent slower on four --
  was taken on that clamped xTB path; do not carry it over to ab initio.
* **Structure is a poor guide to where time goes.** Two separate optimisation
  attempts here were aimed at the wrong file on the strength of reading the
  code; one VTune hotspots collection settled it in minutes. `perf_event_paranoid`
  is 4 on the development box, so PMU counters are unavailable -- use
  `-collect hotspots -knob sampling-mode=sw`.
* **Repeat before believing a timing.** DFT energies repeat to under one per
  cent; gradients used to vary by sixteen, and a single-shot comparison across
  that produced a phantom twenty per cent regression during this work.
* **On Perlmutter, do not let Cray libsci be the BLAS.** The `ftn` wrapper
  links it whatever `BLAS_LIBRARIES` says; its `dgetrf` gets slower with every
  thread (3.8 s on one, 70 s on 64) and is not re-entrant, which segfaulted
  MakeFP at 128 threads. Use the `perlmutter-cpu` preset (sequential MKL, put
  ahead of libsci and checked after every link); see `CMAKE_STYLE.md`. The BLAS being
  sequential means a big GEMM on the unfragmented path has to be split over
  threads by hand -- `gemm_over_columns` and `gemm_over_inner` in
  `mqc_czt_gemm_threads`, and `getrf_threaded` for an LU -- or it runs on one
  core while 127 wait. `perf record -F 49` and binning the samples per second
  by symbol finds every such phase in minutes; `perf_event_paranoid` is 2 on
  Perlmutter compute nodes.
* **Memory decisions read the machine, or the deck.** `memory_budget` in
  `mqc_memory` returns `system.memory_gb` when the deck set it, else a share of
  MemAvailable from /proc/meminfo, else the caller's blind constant. The
  response solver plans on `RESPONSE_BUDGET_SHARE` of the machine's figure.
  The machine's figure is divided among the ranks sharing the node
  (`set_memory_sharers`, called once in `main` and in a session), because
  every fragment worker plans against the same free memory; the deck's figure
  is per rank and taken whole.
  Fixed laptop-sized limits sent a 545-function MakeFP on a 500 GB node down a
  fourteen-hour column build when the transform it refused needed a hundred
  gigabytes and forty seconds.
* **`tools/run.sh` runs `build/mqc` unless told otherwise.** `-b <exe>` or
  `MQC_BINARY` picks another. A LibSci-linked `build/mqc` segfaults in the
  fitted Coulomb build at 956 functions, because that build calls the BLAS from
  every thread; the same deck runs on the `perlmutter-cpu` build.
