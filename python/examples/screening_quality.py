"""How good is the energy screen? Measure it against the answer it approximates.

    mpirun -np 64 python screening_quality.py cluster.xyz --level 2

Three calculations, in this order:

  1. the low-level pass over every term (GFN2) -- the screen itself;
  2. the *full* high-level expansion, every term, no screening -- the reference,
     and the expensive part of this script;
  3. the two-level energy-screened total, for a range of thresholds.

Step 3 costs nothing extra. Once step 2 has computed every high-level term, the
screened answer for *any* threshold is arithmetic over the two breakdowns:

    E(tau) = E_low(all) - sum_{t in S(tau)} delta_low(t) + sum_{t in S(tau)} delta_high(t)

where S(tau) is the kept set closed under subsets. So one reference run buys the
whole threshold curve, and the thing that would be paid for in production -- the
count of high-level terms -- is reported alongside the error it bought.

Two checks fall out of this and are worth more than the table:

  * tau = 0 keeps everything, so E(0) must equal the reference *exactly*. If it
    does not, the recombination arithmetic is wrong, not the screening.
  * --verify runs one threshold as a real screened calculation (pass 1 + pass 2
    over the survivors) and compares it to the reconstructed number. That is the
    check that the cheap curve describes the calculation you would actually run.

A distance screen is reported next to the energy screen at the *same term count*,
because "the energy criterion is better than a cutoff" is the claim this script
exists to test, and comparing them at different costs would not test it.
"""

import argparse
import sys

import mqc

HARTREE_TO_KCAL = 627.5094740631


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("xyz", help="cluster geometry; monomers are perceived from it")
    p.add_argument("--level", type=int, default=2, help="MBE truncation (default 2)")
    p.add_argument("--low-method", default="gfn2", help="screening method (default gfn2)")
    p.add_argument("--high-method", default="hf", help="reference method (default hf)")
    p.add_argument("--high-basis", default="6-31G", help="reference basis (default 6-31G)")
    p.add_argument("--high-functional", default=None, help="reference functional, if any")
    p.add_argument(
        "--high-df",
        action="store_true",
        help="density-fit the reference (never inferred from --high-aux-basis)",
    )
    p.add_argument(
        "--high-aux-basis",
        default=None,
        help="fitting basis; defaults to def2-universal-jkfit when --high-df is set",
    )
    p.add_argument(
        "--thresholds",
        default="0,1e-7,1e-6,1e-5,5e-5,1e-4,5e-4,1e-3",
        help="comma-separated tau values in Hartree on the n-body delta",
    )
    p.add_argument(
        "--verify",
        type=float,
        default=None,
        metavar="TAU",
        help="also run this threshold for real and compare against the reconstruction",
    )
    p.add_argument(
        "--checkpoint",
        action="store_true",
        help="checkpoint each pass, so an interrupted job resumes instead of repeating",
    )
    return p.parse_args(argv)


def main(args):
    system = mqc.System.from_xyz(args.xyz)
    system.auto_monomers()
    print(f"{system}  on {mqc.n_ranks()} rank(s)")

    # Same guard as energy_screened_mbe.py, same reason: auto_monomers is for
    # clusters of intact molecules, and a partition that cuts bonds would make
    # radicals out of every fragment.
    cut = system.missing_bonds()
    if cut:
        raise SystemExit(
            f"the partition cuts {cut} bond(s). Declare the monomers and the "
            "bonds explicitly -- auto_monomers is for molecular clusters."
        )

    low_kwargs = dict(method=args.low_method)
    high_kwargs = dict(method=args.high_method, basis=args.high_basis)
    if args.high_functional:
        high_kwargs["functional"] = args.high_functional
    if args.high_df:
        high_kwargs["density_fitting"] = True
    if args.high_aux_basis:
        # On its own this does nothing -- fitting is switched on by the flag
        # above, not by naming a fitting basis. Said out loud because the
        # combination that looks like it enables DF is aux_basis alone.
        high_kwargs["aux_basis"] = args.high_aux_basis
        if not args.high_df:
            print("  note: --high-aux-basis without --high-df does not fit anything")

    def build(label, **kwargs):
        return mqc.MBE(
            system,
            level=args.level,
            checkpoint=f"chk_{label}.h5" if args.checkpoint else None,
            **kwargs,
        )

    n_terms = len(build("probe", **low_kwargs).terms())
    print(f"level {args.level}: {n_terms} terms")

    # ---- 1: the screen ---------------------------------------------------
    print(f"\n[1/3] low level, every term ({args.low_method})")
    low = build("low", **low_kwargs)
    low_result = low.run("low", write_to_file=True)
    low_rows = {t.monomers: t for t in low_result.breakdown()}
    print(f"      E_low  = {low_result.energy:.10f} Ha")

    # ---- 2: the reference ------------------------------------------------
    # Every term, unscreened. This is what the screened numbers are errors
    # against, and the only reason a small cluster is the right test case.
    print(f"\n[2/3] high level, every term ({describe(high_kwargs)}) -- the reference")
    high = build("high", **high_kwargs)
    high_result = high.run("high", write_to_file=True)
    high_rows = {t.monomers: t for t in high_result.breakdown()}
    print(f"      E_high = {high_result.energy:.10f} Ha")

    missing = set(low_rows) ^ set(high_rows)
    if missing:
        # The reconstruction sums delta_high over the kept set, so a term the
        # reference does not have is a silent zero in the correction.
        raise SystemExit(
            f"the two passes did not compute the same {len(missing)} term(s); "
            "the reconstruction below would be comparing different expansions"
        )

    # ---- 3: the screened totals ------------------------------------------
    print(f"\n[3/3] energy screening, reconstructed from the two breakdowns")

    interactions = [t for t in low_rows.values() if t.level > 1]
    thresholds = sorted({float(x) for x in args.thresholds.split(",")})

    print()
    print("  tau (Ha)   high-level terms   E(tau) [Ha]        err vs ref [kcal/mol]")
    print("  " + "-" * 74)
    rows = []
    for tau in thresholds:
        kept = {t.monomers for t in interactions if abs(t.delta) > tau}
        closed = close(system, args.level, kept)
        total = recombine(low_result.energy, low_rows, high_rows, closed)
        err = (total - high_result.energy) * HARTREE_TO_KCAL
        rows.append((tau, closed, total, err))
        print(
            f"  {tau:8.1e}   {len(closed):>6d} / {n_terms:<8d}   "
            f"{total:18.10f}   {err:+14.6f}"
        )

    # The two ends of the table, named. E(0) is the invariant: keeping every
    # term must reproduce the reference bit for bit, because every low-level
    # delta is subtracted back out again.
    exact = next((r for r in rows if r[0] == 0.0), None)
    if exact and abs(exact[3]) > 1e-6:
        print(
            f"\n  WARNING: tau=0 differs from the reference by {exact[3]:+.3e} kcal/mol. "
            "That is a bug in the recombination, not a screening error."
        )
    monomers_only = recombine(
        low_result.energy, low_rows, high_rows,
        {t for t in low_rows if len(t) == 1},
    )
    print(
        f"\n  for scale: correcting monomers only (no interactions at all) gives "
        f"{(monomers_only - high_result.energy) * HARTREE_TO_KCAL:+.4f} kcal/mol.\n"
        f"  The bare low-level total is off by "
        f"{(low_result.energy - high_result.energy) * HARTREE_TO_KCAL:+.1f} kcal/mol, "
        "which is not an error in any useful sense -- two methods share no\n"
        "  absolute scale. It is there to show what the monomer correction is "
        "carrying, and why leaving it out is not an option."
    )

    # ---- the comparison the script exists for ----------------------------
    distance_table(system, args, low_rows, high_rows, low_result, high_result, rows)

    # ---- optional: is the curve real? ------------------------------------
    if args.verify is not None:
        verify(system, args, low_result, low_rows, high_rows, high_result, high_kwargs)


def recombine(low_energy, low_rows, high_rows, corrected):
    """E_low(all) with the corrected terms swapped out for their high-level ones.

    `corrected` must already be closed under subsets, and must include the
    monomers -- they carry the whole 1-body energy, and two methods share no
    absolute scale, so leaving them at the low level gives a number that is not
    a high-level answer at all.
    """
    delta_low = sum(low_rows[t].delta for t in corrected)
    delta_high = sum(high_rows[t].delta for t in corrected)
    return low_energy - delta_low + delta_high


def close(system, level, kept):
    """Close a set of terms under subsets, the way a real run would.

    Done through `MBE.keep` rather than by hand so the table counts the terms
    the calculation would actually run: a kept trimer drags its dimers back in
    whether or not they passed the screen.
    """
    scratch = mqc.MBE(system, level=level)
    scratch.keep(lambda t: len(t) == 1 or t in kept)
    return set(scratch.terms())


def distance_table(system, args, low_rows, high_rows, low_result, high_result, energy_rows):
    """The same budget, spent geometrically instead of energetically."""
    interactions = [t for t in low_rows.values() if t.level > 1]
    if any(t.distance is None for t in interactions):
        print("\n  (no distance column in the breakdown; skipping the cutoff comparison)")
        return
    by_distance = sorted(interactions, key=lambda t: t.distance)

    print("\n  same term counts, chosen by distance instead:")
    print()
    print("  high-level terms   energy screen [kcal/mol]   distance screen [kcal/mol]")
    print("  " + "-" * 74)
    for _tau, closed, _total, err in energy_rows:
        n_pairs = sum(1 for t in closed if len(t) > 1)
        if n_pairs == 0 or n_pairs == len(interactions):
            continue  # both screens agree trivially at the ends
        picked = {t.monomers for t in by_distance[:n_pairs]}
        d_closed = close(system, args.level, picked)
        d_total = recombine(low_result.energy, low_rows, high_rows, d_closed)
        d_err = (d_total - high_result.energy) * HARTREE_TO_KCAL
        print(f"  {len(closed):>6d}             {err:+14.6f}             {d_err:+14.6f}")
    print(
        "\n  Closure makes the two counts differ slightly above level 2; the "
        "distance rows keep the same number of *interactions*, not of terms."
    )


def verify(system, args, low_result, low_rows, high_rows, high_result, high_kwargs):
    """Run one threshold for real and check the reconstruction predicted it."""
    tau = args.verify
    kept = {t.monomers for t in low_rows.values() if t.level > 1 and abs(t.delta) > tau}
    if not kept:
        print(f"\n[verify] nothing survives tau={tau:g}; nothing to run")
        return

    print(f"\n[verify] running the screen at tau={tau:g} as a real calculation")
    run = mqc.MBE(system, level=args.level, **high_kwargs)
    run.keep(lambda t: len(t) == 1 or t in kept)
    corrected = set(run.terms())
    print(f"[verify] {len(corrected)} terms (closed under subsets)")
    run_result = run.run("verify", write_to_file=True)

    run_rows = {t.monomers: t for t in run_result.breakdown()}
    actual = (
        low_result.energy
        - sum(low_rows[t].delta for t in corrected)
        + sum(run_rows[t].delta for t in corrected)
    )
    predicted = recombine(low_result.energy, low_rows, high_rows, corrected)

    print(f"\n  reconstructed  {predicted:20.10f}")
    print(f"  really run     {actual:20.10f}")
    print(f"  difference     {(actual - predicted) * HARTREE_TO_KCAL:+20.10f} kcal/mol")
    print(
        "\n  These should agree to SCF convergence. A larger gap means the "
        "screened run's fragments are not the terms the reference computed --\n"
        "  the usual cause is a different geometry, partition, or threshold."
    )
    print(
        f"  error vs reference: "
        f"{(actual - high_result.energy) * HARTREE_TO_KCAL:+.6f} kcal/mol"
    )


def describe(kwargs):
    parts = [str(v) for k, v in kwargs.items() if not isinstance(v, bool)]
    parts += [k for k, v in kwargs.items() if v is True]
    return " ".join(parts)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    # The session wraps everything, so a traceback closes it rather than
    # leaving the other ranks blocked until the job's wall clock runs out.
    with mqc.session():
        main(args)
