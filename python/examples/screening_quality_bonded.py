"""How good is the energy screen when the partition cuts covalent bonds?

    mpirun -np 64 python3 screening_quality_bonded.py 6qm1_frag.json --level 2

The covalent twin of `screening_quality.py`, and the same three calculations in
the same order: the low-level pass over every term, the *full* high-level
expansion as the reference, and then the screened total for a range of
thresholds, reconstructed from those two breakdowns at no extra cost.

What changes is not the arithmetic -- it is what the term list means.

A cluster is partitioned by geometry, so `auto_monomers` finds the monomers and
every term of the expansion is optional: drop a distant dimer and you lose a
small interaction. A peptide is partitioned by *decision*: the fragments are
declared, the bonds between them are declared, and every bond that crosses a
fragment boundary is capped with hydrogen. That has three consequences this
script is built around.

  1. The partition is an input, not a perception. It can be wrong -- an atom in
     two fragments, an atom in none, fragment charges that do not add up to the
     molecule's -- in ways that still run and still converge. Those are checked
     before anything is computed.

  2. A bond you cut but never declared is a radical with a converged energy.
     `missing_bonds` is the check, and here it is the *only* thing that must be
     zero: cutting bonds is the point, cutting them silently is the bug.

  3. Most of the term list is not screenable. Two fragments joined by a cut bond
     share a pair of caps sitting on top of each other; their 2-body term is not
     an interaction that can be neglected, it is the repair of a bond the
     partition broke. On the example deck, four of six dimers are of that kind.

So the screen here is applied only to the terms that are through-space, and the
covalently connected ones are kept unconditionally. That is a policy, and a
policy in a screening study has to be defended rather than assumed, so:

  * the script prints the *separation* -- the smallest through-bond delta next
    to the largest through-space one. If through-bond terms sit orders of
    magnitude above everything else, the energy criterion already finds them and
    the protection never fires. That is the outcome that says the criterion is
    sound on covalent partitions, and it is a measurement, not a hope.
  * `--screen-bonded` removes the protection and screens everything on delta
    alone, which prices the policy in kcal/mol instead of arguing about it.

Connectivity, not adjacency in the file, decides this: a term is protected when
its fragments form a connected subgraph of the cut-bond graph. At level 2 that
is "joined by a cut bond"; at level 3 it is a run of three along the chain,
whose 3-body term carries the through-bond coupling a pairwise expansion over
capped fragments is worst at.

The two invariants from the cluster script survive unchanged and are worth more
than the table:

  * tau = 0 keeps everything, so E(0) must equal the reference *exactly*.
    Capping does not enter the recombination -- both passes cap identically --
    so a mismatch is broken arithmetic, not broken chemistry.
  * --verify runs one threshold as a real screened calculation and compares it
    against the reconstruction.

The distance comparison also survives, restricted to the through-space terms,
because a distance screen keeps the covalently bonded pairs for free (they are
1.5 Angstrom apart) and comparing the two screens on terms both of them keep
would flatter the cutoff for a reason that has nothing to do with screening.
"""

import argparse
import itertools
import json
import os
import sys

import mqc

HARTREE_TO_KCAL = 627.5094740631


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "deck",
        help="fragment JSON: the geometry, the partition, the cut bonds and the "
             "fragment charges. The same file mqc would be given.",
    )
    p.add_argument("--xyz", default=None, help="override the geometry the deck names")
    p.add_argument("--level", type=int, default=2, help="MBE truncation (default 2)")
    p.add_argument("--low-method", default="gfn2", help="screening method (default gfn2)")
    p.add_argument("--high-method", default=None, help="reference method (default: the deck's)")
    p.add_argument("--high-basis", default=None, help="reference basis (default: the deck's)")
    p.add_argument("--high-functional", default=None, help="reference functional, if any")
    p.add_argument(
        "--high-df",
        action="store_true",
        help="density-fit the reference (never inferred from an aux basis)",
    )
    p.add_argument(
        "--high-aux-basis",
        default=None,
        help="fitting basis; defaults to the deck's, and does nothing without --high-df",
    )
    p.add_argument(
        "--thresholds",
        default="0,1e-7,1e-6,1e-5,5e-5,1e-4,5e-4,1e-3",
        help="comma-separated tau values in Hartree on the n-body delta",
    )
    p.add_argument(
        "--screen-bonded",
        action="store_true",
        help="screen the through-bond terms too, on delta alone -- the experiment "
             "that prices the protection rather than assuming it",
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


# ---------------------------------------------------------------------------
#  the deck: geometry, partition, cut bonds
# ---------------------------------------------------------------------------


def load_deck(path):
    """Read the fragment JSON, or say why it is not readable.

    Read strictly. A deck that Python will not parse is a deck mqc will not
    parse either, and repairing it here would mean this study runs on a file
    that differs from the one the production job reads -- which is the one
    difference a screening study cannot afford.
    """
    with open(path) as handle:
        text = handle.read()
    try:
        return json.loads(text)
    except json.JSONDecodeError as error:
        raise SystemExit(
            f"{path}:{error.lineno}: {error.msg}.\n"
            "  A trailing comma before a closing brace is the usual cause; JSON "
            "does not allow one. Fix the deck rather than working around it here "
            "-- mqc reads this same file."
        )


def build_system(deck, deck_path, xyz_override):
    """The molecule, its declared partition, and its declared bonds."""
    if not deck.get("molecules"):
        raise SystemExit("the deck has no molecules block")
    if len(deck["molecules"]) > 1:
        raise SystemExit("this script expects one molecule; the deck has several")
    molecule = deck["molecules"][0]

    xyz = xyz_override or molecule.get("xyz")
    if not xyz:
        raise SystemExit(
            "the deck gives no xyz. Inline symbols/geometry are fine for mqc but "
            "not read here; point the deck at a file, or pass --xyz."
        )
    if not os.path.isabs(xyz):
        # Resolved against the deck, the way mqc resolves it -- not against the
        # working directory, so the script runs from anywhere.
        xyz = os.path.join(os.path.dirname(os.path.abspath(deck_path)), xyz)

    charge = int(molecule.get("molecular_charge", 0))
    multiplicity = int(molecule.get("molecular_multiplicity", 1))
    system = mqc.System.from_xyz(xyz, charge=charge, multiplicity=multiplicity)

    fragments = molecule.get("fragments")
    if not fragments:
        raise SystemExit(
            "the deck declares no fragments. That is the whole input of this "
            "script -- for a cluster whose monomers can be perceived, use "
            "screening_quality.py instead."
        )
    charges = molecule.get("fragment_charges") or [0] * len(fragments)
    multiplicities = molecule.get("fragment_multiplicities") or [1] * len(fragments)

    bonds = [tuple(entry[:2]) for entry in molecule.get("connectivity", [])]
    orders = [int(entry[2]) if len(entry) > 2 else 1 for entry in molecule.get("connectivity", [])]

    check_partition(fragments, charges, multiplicities, system.n_atoms, charge)

    system.set_monomers(fragments, charges=charges, multiplicities=multiplicities)
    if bonds:
        system.set_bonds(bonds, orders=orders)
    return system, fragments, bonds, orders


def check_partition(fragments, charges, multiplicities, n_atoms, molecular_charge):
    """Everything about a hand-written partition that runs anyway when wrong.

    None of these are caught by convergence. A duplicated atom gives a fragment
    energy that is too low and an MBE that double-counts it; a dropped atom
    gives a total that is missing a chunk and looks like a screening error;
    fragment charges that do not add up describe a different molecule than the
    one in the xyz, and every term of the expansion converges regardless.
    """
    if len(charges) != len(fragments) or len(multiplicities) != len(fragments):
        raise SystemExit(
            f"{len(fragments)} fragments but {len(charges)} charges and "
            f"{len(multiplicities)} multiplicities"
        )

    seen = {}
    for index, fragment in enumerate(fragments, start=1):
        for atom in fragment:
            if not 0 <= atom < n_atoms:
                raise SystemExit(
                    f"fragment {index} names atom {atom}, which is outside the "
                    f"{n_atoms}-atom geometry (indices are 0-based)"
                )
            if atom in seen:
                raise SystemExit(
                    f"atom {atom} is in fragments {seen[atom]} and {index}. "
                    "Overlapping fragments need GMBE, and this script's "
                    "recombination assumes a partition."
                )
            seen[atom] = index

    orphans = sorted(set(range(n_atoms)) - set(seen))
    if orphans:
        raise SystemExit(
            f"{len(orphans)} atom(s) are in no fragment: {orphans[:10]}"
            f"{' ...' if len(orphans) > 10 else ''}. The expansion would silently "
            "describe a smaller molecule than the geometry."
        )

    if sum(charges) != molecular_charge:
        raise SystemExit(
            f"the fragment charges sum to {sum(charges)} but the molecule is "
            f"{molecular_charge:+d}. Those are two different systems, and both "
            "of them converge."
        )


def cut_bond_graph(fragments, bonds):
    """Which monomers are joined by a cut bond, as an adjacency map.

    1-based, because that is how the expansion counts monomers and how the
    breakdown CSV reports them. A bond whose atoms land in the same fragment is
    not cut and does not appear.
    """
    owner = {atom: index for index, f in enumerate(fragments, start=1) for atom in f}
    adjacency = {index: set() for index in range(1, len(fragments) + 1)}
    cut = []
    for i, j in bonds:
        a, b = owner[i], owner[j]
        if a != b:
            adjacency[a].add(b)
            adjacency[b].add(a)
            cut.append(((i, j), (a, b)))
    return adjacency, cut


def connected(term, adjacency):
    """Is this term a connected piece of the cut-bond graph?

    The protected set. At level 2 it is a pair joined by a cut bond -- a term
    that is not an interaction to be neglected but the repair of a broken bond,
    caps and all. At level 3 it is three fragments in a row, whose 3-body term
    carries the through-bond coupling that a pairwise expansion over capped
    fragments approximates worst.

    A term that is *disconnected* -- two bonded fragments plus a third across
    the fold -- is not protected. Its through-bond part is already in the dimer;
    what the trimer adds is through-space, which is exactly what a delta
    measures well.
    """
    members = set(term)
    if len(members) < 2:
        return True  # monomers are never screened
    stack, seen = [next(iter(members))], set()
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack.extend(adjacency[node] & members - seen)
    return seen == members


# ---------------------------------------------------------------------------
#  the study
# ---------------------------------------------------------------------------


def main(args):
    deck = load_deck(args.deck)
    system, fragments, bonds, _orders = build_system(deck, args.deck, args.xyz)
    adjacency, cut = cut_bond_graph(fragments, bonds)
    print(f"{system}  on {mqc.n_ranks()} rank(s)")
    print(
        f"  {len(fragments)} declared fragments, sizes "
        f"{', '.join(str(len(f)) for f in fragments)}; "
        f"{len(cut)} of {len(bonds)} declared bonds are cut"
    )
    for (i, j), (a, b) in cut:
        print(f"    atoms {i}-{j}  cuts monomer {a} | monomer {b}  -> 2 H caps")

    # The one thing that must be zero. Not "the partition cuts nothing" -- it is
    # meant to cut -- but "it cuts nothing it did not declare". An undeclared cut
    # is an uncapped valence: a radical, with a converged energy, in a
    # closed-shell calculation that never complains.
    undeclared = system.missing_bonds()
    if undeclared:
        raise SystemExit(
            f"the geometry implies {undeclared} bond(s) across fragments that the "
            "deck does not declare. Those valences go uncapped -- add them to "
            "`connectivity`. Nothing below is meaningful until this is zero."
        )

    low_kwargs, high_kwargs = methods(deck, args)

    def build(label, **kwargs):
        return mqc.MBE(
            system,
            level=args.level,
            checkpoint=f"chk_{label}.h5" if args.checkpoint else None,
            **kwargs,
        )

    all_terms = build("probe", **low_kwargs).terms()
    n_terms = len(all_terms)
    bonded_terms = {t for t in all_terms if len(t) > 1 and connected(t, adjacency)}
    screenable = [t for t in all_terms if len(t) > 1 and t not in bonded_terms]
    print(
        f"\nlevel {args.level}: {n_terms} terms -- {len(fragments)} monomers, "
        f"{len(bonded_terms)} through-bond, {len(screenable)} through-space"
    )
    if not screenable:
        raise SystemExit(
            "every interaction in this expansion is through-bond, so there is "
            "nothing for a screen to decide. Raise --level, or fragment further."
        )
    if len(screenable) < 4:
        # Said plainly rather than left for the reader to infer from a
        # three-row table that all reads the same.
        print(
            f"  note: only {len(screenable)} screenable interaction(s). The "
            "invariants below still hold, but the threshold curve has too few\n"
            "  points to calibrate a production threshold from."
        )

    # ---- 1: the screen ---------------------------------------------------
    print(f"\n[1/3] low level, every term ({describe(low_kwargs)})")
    low = build("low", **low_kwargs)
    low_result = low.run("low", write_to_file=True)
    low_rows = {t.monomers: t for t in low_result.breakdown()}
    print(f"      E_low  = {low_result.energy:.10f} Ha")
    report_unconverged(low_result, "low")

    # ---- 2: the reference ------------------------------------------------
    print(f"\n[2/3] high level, every term ({describe(high_kwargs)}) -- the reference")
    high = build("high", **high_kwargs)
    high_result = high.run("high", write_to_file=True)
    high_rows = {t.monomers: t for t in high_result.breakdown()}
    print(f"      E_high = {high_result.energy:.10f} Ha")
    report_unconverged(high_result, "high")

    missing = set(low_rows) ^ set(high_rows)
    if missing:
        raise SystemExit(
            f"the two passes did not compute the same {len(missing)} term(s); "
            "the reconstruction below would be comparing different expansions"
        )

    # ---- the measurement the policy rests on -----------------------------
    separation(low_rows, bonded_terms)

    # ---- 3: the screened totals ------------------------------------------
    print("\n[3/3] energy screening, reconstructed from the two breakdowns")
    if args.screen_bonded:
        print("      --screen-bonded: through-bond terms are screened on delta too")
    else:
        print("      through-bond terms are kept at every threshold")

    interactions = [t for t in low_rows.values() if t.level > 1]
    thresholds = sorted({float(x) for x in args.thresholds.split(",")})

    print()
    print("  tau (Ha)   high-level terms   rescued   E(tau) [Ha]        err vs ref [kcal/mol]")
    print("  " + "-" * 84)
    rows = []
    for tau in thresholds:
        passed = {t.monomers for t in interactions if abs(t.delta) > tau}
        # How many through-bond terms the delta criterion would have dropped on
        # its own. Zero all the way down the column is the result that says the
        # protection is redundant -- which is the good news, not the boring news.
        rescued = len(bonded_terms - passed)
        kept = passed if args.screen_bonded else passed | bonded_terms
        closed = close(system, args.level, kept)
        total = recombine(low_result.energy, low_rows, high_rows, closed)
        err = (total - high_result.energy) * HARTREE_TO_KCAL
        rows.append((tau, closed, total, err))
        print(
            f"  {tau:8.1e}   {len(closed):>6d} / {n_terms:<8d}   {rescued:>7d}   "
            f"{total:18.10f}   {err:+14.6f}"
        )

    exact = next((r for r in rows if r[0] == 0.0), None)
    if exact and abs(exact[3]) > 1e-6:
        print(
            f"\n  WARNING: tau=0 differs from the reference by {exact[3]:+.3e} kcal/mol. "
            "That is a bug in the recombination, not a screening error, and\n"
            "  capping is not the suspect: both passes cap the same fragments the "
            "same way, so the caps cancel term by term."
        )

    monomers_only = recombine(
        low_result.energy, low_rows, high_rows, {t for t in low_rows if len(t) == 1},
    )
    print(
        f"\n  for scale: correcting the capped monomers only (no interactions at "
        f"all) gives {(monomers_only - high_result.energy) * HARTREE_TO_KCAL:+.4f} "
        "kcal/mol.\n"
        "  That number is larger here than it would be for a cluster, and the cut "
        "bonds are why: the caps have to be paid for by the very\n"
        "  through-bond terms this expansion cannot drop."
    )

    # ---- the comparison the script exists for ----------------------------
    distance_table(system, args, low_rows, high_rows, low_result, high_result,
                   rows, bonded_terms)

    # ---- optional: is the curve real? ------------------------------------
    if args.verify is not None:
        verify(system, args, low_result, low_rows, high_rows, high_result,
               high_kwargs, bonded_terms)


def separation(low_rows, bonded_terms):
    """Does the delta criterion find the cut bonds without being told?

    The whole defence of screening a covalent partition on energy. If the
    smallest through-bond delta is far above the largest through-space one, then
    any threshold that keeps the interactions worth keeping also keeps every
    broken bond, and the protection is a seatbelt that never engages. If the two
    overlap, a threshold set from the through-space terms will quietly drop a
    bond repair, and the protection is doing real work -- which `--screen-bonded`
    then measures in kcal/mol.
    """
    bonded = [abs(t.delta) for t in low_rows.values() if t.monomers in bonded_terms]
    space = [
        abs(t.delta) for t in low_rows.values()
        if t.level > 1 and t.monomers not in bonded_terms
    ]
    if not bonded or not space:
        return
    low_bonded, high_space = min(bonded), max(space)
    print(
        f"\n  through-bond deltas: {min(bonded):.2e} to {max(bonded):.2e} Ha\n"
        f"  through-space:       {min(space):.2e} to {max(space):.2e} Ha"
    )
    if low_bonded > high_space:
        print(
            f"  They do not overlap -- the smallest bond repair is "
            f"{low_bonded / high_space:.0f}x the largest through-space term, so a\n"
            "  threshold chosen for the latter cannot reach the former. The energy "
            "criterion finds the cut bonds by itself."
        )
    else:
        print(
            f"  They overlap: a through-space term at {high_space:.2e} Ha sits above "
            f"a bond repair at {low_bonded:.2e} Ha.\n"
            "  A threshold set from the through-space column would drop that bond "
            "repair. Run --screen-bonded to see what it costs."
        )


def recombine(low_energy, low_rows, high_rows, corrected):
    """E_low(all) with the corrected terms swapped out for their high-level ones.

    `corrected` must already be closed under subsets, and must include the
    monomers -- they carry the whole 1-body energy, caps included, and two
    methods share no absolute scale.
    """
    delta_low = sum(low_rows[t].delta for t in corrected)
    delta_high = sum(high_rows[t].delta for t in corrected)
    return low_energy - delta_low + delta_high


def close(system, level, kept):
    """Close a set of terms under subsets, the way a real run would.

    Done through `MBE.keep` rather than by hand so the table counts the terms
    the calculation would actually run: a kept trimer drags its dimers back in
    whether or not they passed the screen. On a covalent partition that closure
    does more work than on a cluster -- a through-space trimer typically pulls
    in the through-bond dimer underneath it.
    """
    scratch = mqc.MBE(system, level=level)
    scratch.keep(lambda t: len(t) == 1 or t in kept)
    return set(scratch.terms())


def distance_table(system, args, low_rows, high_rows, low_result, high_result,
                   energy_rows, bonded_terms):
    """The same budget, spent geometrically instead of energetically.

    Restricted to the through-space terms in both columns. A distance screen
    keeps the covalently bonded pairs for nothing -- they are 1.5 Angstrom
    apart, the shortest separations in the system -- so scoring it on terms both
    screens keep would credit the cutoff with the one thing the geometry cannot
    get wrong, and hide the question of what it does with the rest.
    """
    space = [
        t for t in low_rows.values()
        if t.level > 1 and t.monomers not in bonded_terms
    ]
    if not space:
        return
    if any(t.distance is None for t in space):
        print("\n  (no distance column in the breakdown; skipping the cutoff comparison)")
        return
    by_distance = sorted(space, key=lambda t: t.distance)

    print("\n  same through-space term counts, chosen by distance instead:")
    print()
    print("  through-space terms   energy screen [kcal/mol]   distance screen [kcal/mol]")
    print("  " + "-" * 78)
    shown = set()
    for _tau, closed, _total, err in energy_rows:
        n_space = sum(1 for t in closed if len(t) > 1 and t not in bonded_terms)
        if n_space in shown or n_space == 0 or n_space == len(space):
            continue  # both screens agree trivially at the ends
        shown.add(n_space)
        picked = {t.monomers for t in by_distance[:n_space]}
        if not args.screen_bonded:
            picked |= bonded_terms
        d_closed = close(system, args.level, picked)
        d_total = recombine(low_result.energy, low_rows, high_rows, d_closed)
        d_err = (d_total - high_result.energy) * HARTREE_TO_KCAL
        print(f"  {n_space:>9d}             {err:+14.6f}             {d_err:+14.6f}")
    if not shown:
        print("  (every threshold lands on the same term count; nothing to compare)")
    print(
        "\n  Closure makes the totals differ slightly above level 2; the rows keep "
        "the same number of through-space *interactions*, not of terms."
    )


def verify(system, args, low_result, low_rows, high_rows, high_result, high_kwargs,
           bonded_terms):
    """Run one threshold for real and check the reconstruction predicted it."""
    tau = args.verify
    passed = {t.monomers for t in low_rows.values() if t.level > 1 and abs(t.delta) > tau}
    kept = passed if args.screen_bonded else passed | bonded_terms
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
        "  on a covalent partition the usual cause is a changed partition or a "
        "changed bond list, either of which moves the caps and so\n"
        "  changes every fragment that touches them."
    )
    print(
        f"  error vs reference: "
        f"{(actual - high_result.energy) * HARTREE_TO_KCAL:+.6f} kcal/mol"
    )


def report_unconverged(result, label):
    """Capped fragments converge worse than whole molecules; say so when they do.

    A cap sits at a bond length from the atom it replaces and leaves a fragment
    that is not a molecule anyone would have drawn, so the pass that struggles
    is usually this one rather than the cluster equivalent. Every unconverged
    term below is an energy the screen and the reference both used.
    """
    bad = result.unconverged()
    if bad:
        print(f"      WARNING: {len(bad)} unconverged term(s) in the {label} pass: "
              f"{', '.join(str(t.monomers) for t in bad[:6])}"
              f"{' ...' if len(bad) > 6 else ''}")


def methods(deck, args):
    """The two levels. The high one defaults to whatever the deck asks for.

    Taking the reference from the deck is the point of reading the deck: the
    number this study calls "the answer" is then the number the production job
    would have produced, rather than one this script chose.
    """
    model = deck.get("model", {}) or {}
    low_kwargs = dict(method=args.low_method)

    method = args.high_method or model.get("method")
    if not method:
        raise SystemExit("no reference method: the deck has no model.method and "
                         "--high-method was not given")
    high_kwargs = dict(method=method)

    basis = args.high_basis or model.get("basis")
    if basis:
        high_kwargs["basis"] = basis
    functional = args.high_functional or model.get("functional")
    if functional:
        high_kwargs["functional"] = functional

    if args.high_df:
        high_kwargs["density_fitting"] = True
    aux = args.high_aux_basis or model.get("aux_basis")
    if aux:
        # On its own this does nothing -- fitting is switched on by --high-df,
        # not by naming a fitting basis. Worth saying out loud here because a
        # deck that names an aux basis and no `density_fitting` reads, to the
        # eye, as a density-fitted calculation.
        high_kwargs["aux_basis"] = aux
        if not args.high_df:
            print(f"  note: aux basis {aux} is set but nothing is fitted without --high-df")
    return low_kwargs, high_kwargs


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
