"""Two-pass energy-screened MBE: cheap everywhere, expensive where it matters.

    mpirun -np 64 python energy_screened_mbe.py cluster.xyz

Pass 1 runs the whole expansion with GFN2, which is cheap enough to do every
term. Pass 2 keeps the terms whose n-body contribution came back above a
threshold and recomputes those with DFT. The point is that the criterion is
*energetic* rather than geometric: a distance cutoff throws away a strong
interaction that happens to be far apart and keeps a negligible one that
happens to be close, and neither mistake is visible in the answer.

The total is then

    E = E_low(everything) + sum over kept terms of [E_high - E_low]

which is why pass 1 has to cover every term, not just the ones that survive.

This is a script, not a library: read it, change the threshold, change the
functional. Everything it does is in the four calls below.
"""

import sys

import mqc

THRESHOLD = 1.0e-4  #: Hartree, on the n-body delta. See below -- calibrate it.
LEVEL = 3
LOW = dict(method="gfn2")
HIGH = dict(method="dft", functional="pbe0", basis="def2-svp")


def main(path):
    # Rank 0 only from here down. Every other rank is inside Fortran, waiting.
    system = mqc.System.from_xyz(path)
    system.auto_monomers()
    print(f"{system}  on {mqc.n_ranks()} rank(s)")

    # A cluster of intact molecules has no bonds crossing monomers. If this is
    # not zero the partition is cutting something, and the fragments would be
    # radicals -- `run` refuses rather than computing them, but knowing here
    # is better than finding out after the queue starts.
    cut = system.missing_bonds()
    if cut:
        raise SystemExit(
            f"the partition cuts {cut} bond(s). Declare the monomers and the "
            "bonds explicitly -- auto_monomers is for molecular clusters."
        )

    # ---- pass 1: everything, cheaply ------------------------------------
    low = mqc.MBE(system, level=LEVEL, **LOW)
    print(f"pass 1: {len(low.terms())} terms at {LOW['method']}")
    low_result = low.run("screen", write_to_file=True)
    print(f"pass 1: E = {low_result.energy:.10f} Ha")

    # ---- choose ---------------------------------------------------------
    # The breakdown is read back from the CSV pass 1 just wrote, which is also
    # what makes this restartable: point `terms` at an older run's file and
    # pass 2 goes ahead without repeating pass 1.
    breakdown = low_result.breakdown()
    strong = {t.monomers for t in breakdown if t.level > 1 and abs(t.delta) > THRESHOLD}
    print(
        f"screen: {len(strong)} of {sum(1 for t in breakdown if t.level > 1)} "
        f"interactions above {THRESHOLD:g} Ha"
    )
    report(breakdown, strong)

    if not strong:
        raise SystemExit("nothing above the threshold; lower it or stop here")

    # ---- pass 2: the survivors, expensively ------------------------------
    # `keep` closes the list under subsets afterwards, so this may run more
    # terms than the screen chose -- a trimer needs its dimers whether or not
    # they were interesting on their own.
    high = mqc.MBE(system, level=LEVEL, **HIGH)
    high.keep(lambda t: len(t) == 1 or t in strong)
    print(f"pass 2: {len(high.terms())} terms at {HIGH['functional']} (closed under subsets)")
    high_result = high.run("refined", write_to_file=True)

    # ---- the correction --------------------------------------------------
    # Only the kept terms are recomputed, so the difference between the two
    # passes over those terms is the correction to the cheap total.
    low_over_kept = sum(t.delta for t in breakdown if t.monomers in strong)
    high_over_kept = sum(
        t.delta for t in high_result.breakdown() if t.monomers in strong
    )
    total = low_result.energy - low_over_kept + high_over_kept

    print()
    print(f"  E[{LOW['method']}, all terms]        {low_result.energy:20.10f}")
    print(f"  correction over {len(strong)} terms  {high_over_kept - low_over_kept:+20.10f}")
    print(f"  E[two-level]                 {total:20.10f}")


def report(breakdown, strong):
    """Say what the threshold did, in the units the decision was made in."""
    deltas = sorted(
        (abs(t.delta) for t in breakdown if t.level > 1), reverse=True
    )
    if not deltas:
        return
    print(
        f"  deltas span {deltas[0]:.2e} to {deltas[-1]:.2e} Ha; "
        f"the cut discards {sum(deltas[len(strong):]):.2e} Ha of them"
    )
    # That last number is the error the screen introduces at the low level,
    # and is the thing to look at when choosing a threshold -- not the count.


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    # The session wraps everything, so a traceback closes it rather than
    # leaving the other ranks blocked until the job's wall clock runs out.
    with mqc.session():
        main(sys.argv[1])
