"""Check the --dump-terms output of the two screening-quality scripts.

    python3 test_screening_dump.py

Plain python3, no pytest, no MPI and no compiled library. The scripts it tests
need all three to run a calculation, which is exactly why the dump was factored
into a function that takes the two breakdowns and a kept set and nothing else:
that part is arithmetic over rows, and arithmetic is testable here.

Two stand-ins make that work, and both are deliberate:

  * `Term` is a copy of the class in `mqc/__init__.py` rather than an import of
    it, because importing `mqc` loads `libmqc.so` at module scope -- there is no
    way to reach the class without a build. It carries the attribute names the
    dump reads and nothing else. If those names ever change in mqc, this file
    keeps passing while the scripts break, which is the one thing it cannot
    check.
  * `sys.modules["mqc"]` is stubbed before the example scripts are imported, for
    the same reason: they import mqc at module scope, and none of the functions
    under test touch it.

The check worth having is the last one. The dump's error total and the scripts'
own recombination are the same quantity by two routes, so they must agree
exactly rather than closely; the test pins that at 1e-12 Hartree on numbers of
order 1e-2, which is float noise and no more.
"""

import contextlib
import csv
import importlib.util
import io
import os
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))


class Term:
    """The attributes of `mqc.Term` that a breakdown row is read through."""

    def __init__(self, monomers, energy, delta, distance=None):
        self.monomers = monomers
        self.energy = energy
        self.delta = delta
        self.distance = distance

    @property
    def level(self):
        return len(self.monomers)


def load(name):
    """Import one of the example scripts with a stubbed-out mqc."""
    sys.modules.setdefault("mqc", types.ModuleType("mqc"))
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
#  the hand-made expansion
# ---------------------------------------------------------------------------
#
#  Four monomers, all six dimers, and deltas chosen so that the ordering by
#  |delta_low| (what the screen sees) is not the ordering by |correction| (what
#  the error is made of). A dump that sorted on the wrong column would pass on
#  rows where the two agree.

MONOMERS = {
    (1,): (-76.0, -76.0),
    (2,): (-76.1, -76.1),
    (3,): (-76.2, -76.2),
    (4,): (-76.3, -76.3),
}

DIMERS = {
    #                delta_low   delta_high   distance
    (1, 2): (-2.0e-2, -2.05e-2, 2.8),
    (1, 3): (-4.0e-3, -3.0e-3, 4.1),
    (1, 4): (-9.0e-4, -1.4e-3, 6.2),   # small delta, large correction
    (2, 3): (-3.0e-3, -3.01e-3, 4.4),
    (2, 4): (-5.0e-4, -5.2e-4, 7.0),
    (3, 4): (-1.0e-4, -1.1e-4, 8.3),
}


def breakdowns(with_distance=True):
    low, high = {}, {}
    for key, (e_low, e_high) in MONOMERS.items():
        low[key] = Term(key, e_low, e_low)
        high[key] = Term(key, e_high, e_high)
    for key, (d_low, d_high, distance) in DIMERS.items():
        r = distance if with_distance else None
        low[key] = Term(key, sum(MONOMERS[(m,)][0] for m in key) + d_low, d_low, r)
        high[key] = Term(key, sum(MONOMERS[(m,)][1] for m in key) + d_high, d_high, r)
    return low, high


def low_total(low_rows):
    """E_low over the whole expansion: every delta, monomers included."""
    return sum(t.delta for t in low_rows.values())


def read_dump(path):
    """The dump back as (comment lines, header, rows)."""
    with open(path) as handle:
        lines = handle.readlines()
    comments = [ln for ln in lines if ln.startswith("#")]
    body = [ln for ln in lines if not ln.startswith("#")]
    reader = csv.reader(body)
    header = next(reader)
    return comments, header, [row for row in reader if row]


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def run_dump(module, *args):
    """Call `dump_terms` with its printing captured, and hand back both.

    The stdout half is part of what the option promises -- the two top-ten
    lists and the dropped total -- so it is checked rather than silenced,
    and capturing it keeps this test's own output to one line.
    """
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        path, dropped = module.dump_terms(*args)
    printed = buffer.getvalue()
    check("KEPT terms" in printed and "DROPPED terms" in printed,
          "the dump printed no kept/dropped listing")
    check("dropped corrections sum to" in printed,
          "the dump printed no total for the dropped corrections")
    return path, dropped


# ---------------------------------------------------------------------------
#  the tests
# ---------------------------------------------------------------------------


def check_cluster(module):
    low_rows, high_rows = breakdowns()
    tau = 1e-3
    kept = {k for k, t in low_rows.items() if t.level > 1 and abs(t.delta) > tau}
    closed = kept | {k for k in low_rows if len(k) == 1}

    with tempfile.TemporaryDirectory() as scratch:
        os.chdir(scratch)
        path, dropped = run_dump(module, low_rows, high_rows, closed, tau)
        check(path == "screening_terms_tau0.001.csv", f"tau formatting: {path}")
        comments, header, rows = read_dump(path)

    check(header == list(module.DUMP_COLUMNS), f"columns: {header}")
    check(header == ["monomers", "level", "distance", "kept", "delta_low_ha",
                     "delta_high_ha", "correction_ha", "correction_kcal"],
          f"columns: {header}")
    check(comments, "no header comment explaining the correction column")

    # One row per interaction, monomers excluded.
    check(len(rows) == len(DIMERS), f"{len(rows)} rows for {len(DIMERS)} dimers")
    check(all("-" in row[0] for row in rows), "a monomer reached the dump")

    corrections = [float(row[6]) for row in rows]
    check(corrections == sorted(corrections, key=lambda c: -abs(c)),
          f"not sorted by |correction| descending: {corrections}")
    # The sort is not the same as sorting on delta_low, or this proves nothing.
    deltas = [abs(float(row[4])) for row in rows]
    check(deltas != sorted(deltas, reverse=True),
          "the fixture no longer distinguishes the two orderings")

    for row in rows:
        key = tuple(int(m) for m in row[0].split("-"))
        check(int(row[3]) == (1 if key in closed else 0), f"kept flag on {key}")
        check(row[2] != "", f"distance dropped from {key}")
        check(abs(float(row[6]) - (high_rows[key].delta - low_rows[key].delta))
              < 1e-15, f"correction on {key}")
        check(abs(float(row[7]) - float(row[6]) * module.HARTREE_TO_KCAL)
              < 1e-6, f"kcal column on {key}")

    identity(module, low_rows, high_rows, closed, dropped)

    # No distance column in the breakdown: the field is empty, not "None".
    low_rows, high_rows = breakdowns(with_distance=False)
    with tempfile.TemporaryDirectory() as scratch:
        os.chdir(scratch)
        path, _ = run_dump(module, low_rows, high_rows, closed, tau)
        _, _, rows = read_dump(path)
    check(all(row[2] == "" for row in rows), "empty distance not written as empty")


def check_bonded(module):
    low_rows, high_rows = breakdowns()
    # Monomers 1-2 and 3-4 joined by cut bonds, as `connected` would find them.
    protected = {(1, 2), (3, 4)}
    tau = 1e-3
    passed = {k for k, t in low_rows.items() if t.level > 1 and abs(t.delta) > tau}
    kept = passed | protected
    closed = kept | {k for k in low_rows if len(k) == 1}

    with tempfile.TemporaryDirectory() as scratch:
        os.chdir(scratch)
        path, dropped = run_dump(module, low_rows, high_rows, closed, tau, protected)
        _comments, header, rows = read_dump(path)

    check(header == list(module.DUMP_COLUMNS), f"columns: {header}")
    check(header == ["monomers", "level", "distance", "kept", "protected",
                     "delta_low_ha", "delta_high_ha", "correction_ha",
                     "correction_kcal"], f"columns: {header}")

    corrections = [float(row[7]) for row in rows]
    check(corrections == sorted(corrections, key=lambda c: -abs(c)),
          f"not sorted by |correction| descending: {corrections}")

    for row in rows:
        key = tuple(int(m) for m in row[0].split("-"))
        check(int(row[3]) == (1 if key in closed else 0), f"kept flag on {key}")
        check(int(row[4]) == (1 if key in protected else 0), f"protected flag on {key}")

    # (3,4) has the smallest delta of any dimer and is kept anyway: that is the
    # protection, and it is what the extra column exists to make visible.
    by_key = {row[0]: row for row in rows}
    check(by_key["3-4"][3] == "1" and by_key["3-4"][4] == "1",
          "the protected dimer was screened out")

    identity(module, low_rows, high_rows, closed, dropped)


def identity(module, low_rows, high_rows, closed, dropped):
    """E(tau) - E_ref = -sum of the dropped corrections. Exactly, not nearly.

    The reference is the recombination over every term, which is what a full
    high-level pass would have produced; E(tau) is the recombination over the
    kept set. Both go through the script's own `recombine`, so this crosses the
    dump against the arithmetic the tables are built from.
    """
    energy = low_total(low_rows)
    reference = module.recombine(energy, low_rows, high_rows, set(low_rows))
    total = module.recombine(energy, low_rows, high_rows, closed)
    residual = (total - reference) + dropped
    check(abs(residual) < 1e-12, f"identity violated by {residual:.3e} Ha")
    # A dropped term is not absent from E(tau) -- it is there at its low-level
    # delta -- so the error carries the opposite sign to the corrections that
    # were dropped. Pinned here because getting that sign backwards produces a
    # dump that looks right and is off by a factor of two against the table.
    check(abs((total - reference) - -dropped) < 1e-12, "sign of the error")


def main():
    cluster = load("screening_quality")
    bonded = load("screening_quality_bonded")
    check_cluster(cluster)
    check_bonded(bonded)
    print("ok: dump columns, sort order, kept/protected flags and the error identity")


if __name__ == "__main__":
    cwd = os.getcwd()
    try:
        main()
    finally:
        os.chdir(cwd)
