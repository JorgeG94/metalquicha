"""metalquicha from Python.

Set a molecule up here; the MPI lives in Fortran.

    import mqc

    with mqc.session():
        water = mqc.System.from_xyz("cluster.xyz")
        water.auto_monomers()

        mbe = mqc.MBE(water, level=3, method="gfn2")
        low = mbe.run("screen")

**Every rank runs this file until `session()` is entered, and only rank 0
comes out of it.** Under `mpirun -np 64 python script.py` there are 64
interpreters; 63 of them enter `session()` and never return, spending the rest
of their lives inside Fortran waiting to be told what to compute. So:

  * Enter the session first. Anything above it -- reading a structure,
    parsing arguments -- happens 64 times, on a shared filesystem.
  * Anything after it happens once. Which is what you want, and is not
    obvious from reading the script.

The session must be closed or the other 63 ranks wait forever on a message
that never comes, which on a batch system burns the allocation rather than
failing. `with mqc.session()` handles that including on an exception; the bare
`begin()`/`end()` pair is there for a REPL.
"""

import contextlib
import ctypes
import json
import os

from . import _ffi

__all__ = [
    "session",
    "begin",
    "end",
    "n_ranks",
    "rank",
    "System",
    "MBE",
    "Result",
    "MQCError",
    "ELEMENTS",
]

ELEMENTS = (
    # H through Og, matching mqc_elements.f90 exactly. It is a second copy of
    # that table and there is no mechanism keeping the two in step, so if this
    # ever needs editing, the real fix is to stop parsing structures here --
    # see System.from_xyz.
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I "
    "Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt "
    "Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr "
    "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()

_NUMBER = {sym.upper(): i + 1 for i, sym in enumerate(ELEMENTS)}


class MQCError(RuntimeError):
    """Something Fortran refused to do, with the reason it gave."""


def _check(status, reader, reversed_args=False):
    if status != 0:
        raise MQCError(_ffi.message(reader, reversed_args) or f"status {status}")


# ---------------------------------------------------------------------------
#  Session
# ---------------------------------------------------------------------------


def begin():
    """Start MPI. Returns on rank 0; every other rank never comes back."""
    _check(_ffi.session_begin(), _ffi.session_last_error)


def end():
    """Release the workers and shut MPI down. Not optional."""
    _check(_ffi.session_end(), _ffi.session_last_error)


@contextlib.contextmanager
def session():
    """A session that is closed even if the body raises.

    An uncaught exception on rank 0 would otherwise leave every other rank
    blocked on a broadcast that never arrives -- a hang rather than a
    traceback, and on a batch system a hang that costs the whole allocation.
    """
    begin()
    try:
        yield
    finally:
        end()


def n_ranks():
    return int(_ffi.session_size())


def rank():
    return int(_ffi.session_rank())


# ---------------------------------------------------------------------------
#  System
# ---------------------------------------------------------------------------


class System:
    """A molecule, its partition into monomers, and its bonds.

    Coordinates are Angstrom here and Bohr inside; the boundary converts.
    Atom indices are 0-based, as Python expects. Monomer indices, where they
    appear in a term list, are 1-based -- that is the expansion's own
    convention and is not worth hiding, because the CSV a screen reads back
    uses it too.
    """

    def __init__(self, symbols=None, coords=None, numbers=None, charge=0, multiplicity=1):
        self._handle = _ffi.system_new()
        if not self._handle:
            raise MQCError("could not allocate a system")
        self._monomers = None
        if symbols is not None or numbers is not None:
            self.set_geometry(symbols, coords, numbers, charge, multiplicity)

    # -- construction -------------------------------------------------------

    @classmethod
    def from_xyz(cls, path, charge=0, multiplicity=1):
        """Read an .xyz file. Angstrom, as the format says."""
        with open(path) as handle:
            lines = handle.read().splitlines()
        n = int(lines[0].split()[0])
        symbols, coords = [], []
        for line in lines[2 : 2 + n]:
            parts = line.split()
            symbols.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
        return cls(symbols=symbols, coords=coords, charge=charge, multiplicity=multiplicity)

    def set_geometry(self, symbols=None, coords=None, numbers=None, charge=0, multiplicity=1):
        if numbers is None:
            if symbols is None:
                raise ValueError("give either symbols or numbers")
            numbers = [_NUMBER[str(s).strip().upper()] for s in symbols]
        flat = [float(x) for atom in coords for x in atom] if coords and hasattr(
            coords[0], "__len__"
        ) else [float(x) for x in coords]
        if len(flat) != 3 * len(numbers):
            raise ValueError(f"{len(numbers)} atoms but {len(flat)} coordinates")
        _check(
            _ffi.system_set_geometry(
                self._handle, len(numbers), _ffi.ints(numbers), _ffi.doubles(flat),
                int(charge), int(multiplicity),
            ),
            _ffi.system_last_error,
        )
        return self

    # -- partition ----------------------------------------------------------

    def auto_monomers(self, tolerance=1.2):
        """Work the monomers out from the geometry.

        Connected components under a covalent-radius criterion, so for a
        cluster of intact molecules this is exactly right and for anything
        held together by covalent bonds it returns one monomer and there is
        nothing to expand. Fine for water, wrong for a peptide -- there,
        say which atoms go where and which bonds you cut.
        """
        _check(_ffi.system_auto_monomers(self._handle, float(tolerance)), _ffi.system_last_error)
        self._monomers = None
        return self

    def set_monomers(self, monomers, charges=None, multiplicities=None):
        """The partition, as a list of lists of 0-based atom indices."""
        monomers = [list(m) for m in monomers]
        n = len(monomers)
        width = max(len(m) for m in monomers)
        # Column per monomer, which is what fragment_atoms already is.
        packed = []
        for m in monomers:
            packed.extend(list(m) + [0] * (width - len(m)))
        _check(
            _ffi.system_set_monomers(
                self._handle, n, width,
                _ffi.ints(len(m) for m in monomers),
                _ffi.ints(packed),
                _ffi.ints(charges if charges is not None else [0] * n),
                _ffi.ints(multiplicities if multiplicities is not None else [1] * n),
            ),
            _ffi.system_last_error,
        )
        self._monomers = monomers
        return self

    def set_bonds(self, bonds, orders=None):
        """Declare the connectivity, as (i, j) pairs of 0-based atoms.

        A bond whose atoms land in different monomers is marked broken and
        capped; that is derived, not declared, so there is no flag to get
        wrong. What cannot be derived is a bond you never mentioned -- see
        `missing_bonds`.
        """
        bonds = [tuple(b) for b in bonds]
        n = len(bonds)
        _check(
            _ffi.system_set_bonds(
                self._handle, n,
                _ffi.ints(b[0] for b in bonds),
                _ffi.ints(b[1] for b in bonds),
                _ffi.ints(orders if orders is not None else [1] * n),
                _ffi.ints([0] * n),
            ),
            _ffi.system_last_error,
        )
        return self

    def perceive_bonds(self, tolerance=1.2):
        """Fill the bond list in from the geometry, and declare it."""
        _check(_ffi.system_perceive_bonds(self._handle, float(tolerance)), _ffi.system_last_error)
        return self

    def missing_bonds(self, tolerance=1.2):
        """Bonds the geometry implies across monomers that you did not declare.

        Non-zero means fragments with uncapped valences -- radicals, and a
        converged energy for them. `run` refuses rather than reporting this,
        so this is for looking before you get there.
        """
        return int(_ffi.system_count_missing_bonds(self._handle, float(tolerance)))

    # -- state --------------------------------------------------------------

    @property
    def n_atoms(self):
        return int(_ffi.system_n_atoms(self._handle))

    @property
    def n_monomers(self):
        return int(_ffi.system_n_monomers(self._handle))

    @property
    def n_bonds(self):
        return int(_ffi.system_n_bonds(self._handle))

    @property
    def bonds_declared(self):
        """Whether bonds were stated. Not the same as having any."""
        return bool(_ffi.system_bonds_declared(self._handle))

    def __repr__(self):
        return (
            f"<mqc.System {self.n_atoms} atoms, {self.n_monomers} monomers, "
            f"{self.n_bonds} bonds{'' if self.bonds_declared else ' (undeclared)'}>"
        )

    def __del__(self):
        handle, self._handle = getattr(self, "_handle", None), None
        if handle:
            _ffi.system_free(handle)


# ---------------------------------------------------------------------------
#  Terms and results
# ---------------------------------------------------------------------------


class Result:
    """What came back, and where the rest of it was written."""

    def __init__(self, energy, label, wrote):
        self.energy = energy
        self.label = label
        self.wrote = wrote

    @property
    def breakdown_csv(self):
        """The per-fragment CSV, if files were written.

        This is where a screen reads from: one row per term, with the monomer
        indices, the term energy and its n-body delta.
        """
        return f"output_{self.label}_fragments.csv" if self.wrote else None

    def breakdown(self):
        """The per-term rows: monomers, energy, n-body delta, distance.

        Distinct from `MBE.terms`, which is only the monomer indices -- these
        rows carry what the calculation found, and are what a threshold is
        applied to.
        """
        path = self.breakdown_csv
        if not path or not os.path.exists(path):
            raise MQCError(
                f"no breakdown file for '{self.label}'; run with write_to_file=True"
            )
        rows = []
        with open(path) as handle:
            header = next(handle).strip().split(",")
            mcols = [i for i, name in enumerate(header) if name.startswith("m")]
            ie, idelta = header.index("energy"), header.index("delta_energy")
            idist = header.index("distance") if "distance" in header else None
            iscf = header.index("scf") if "scf" in header else None
            for line in handle:
                if not line.strip():
                    continue
                cells = line.split(",")
                rows.append(
                    Term(
                        tuple(int(cells[i]) for i in mcols if int(cells[i]) > 0),
                        float(cells[ie]),
                        float(cells[idelta]),
                        float(cells[idist]) if idist is not None else None,
                        _converged(cells[iscf]) if iscf is not None else None,
                    )
                )
        return rows

    def unconverged(self):
        """Terms whose SCF did not converge, from the breakdown.

        Empty is the answer you want. Non-empty means this energy was built
        from those terms anyway -- which only happens when the run asked for
        it with allow_crap_scf, and is exactly the list to screen out or
        recompute.
        """
        return [t for t in self.breakdown() if t.converged is False]

    def __repr__(self):
        return f"<mqc.Result {self.energy:.10f} Ha>"


class Term:
    """One term of the expansion, as read back from a breakdown."""

    __slots__ = ("monomers", "energy", "delta", "distance", "converged")

    def __init__(self, monomers, energy, delta, distance, converged=None):
        self.monomers = monomers  # 1-based, as the expansion counts them
        self.energy = energy
        self.delta = delta  #: n-body contribution -- what a threshold is on
        self.distance = distance
        self.converged = converged
            #: True, False, or None when the method did not report. None is
            #: not a claim that it converged.

    @property
    def level(self):
        return len(self.monomers)

    def __repr__(self):
        return f"<Term {self.monomers} delta={self.delta:+.3e}>"


# ---------------------------------------------------------------------------
#  MBE
# ---------------------------------------------------------------------------


class MBE:
    """A many-body expansion over a system.

    The settings become a JSON document -- the same one `mqc` reads from a
    deck, less the molecules -- which is what crosses to the other ranks.
    `keywords` is the escape hatch: anything the deck accepts can be put
    there, and it is merged over what the arguments produced.
    """

    def __init__(
        self,
        system,
        level=2,
        method="gfn2",
        basis=None,
        functional=None,
        driver="Energy",
        cutoffs=None,
        unchecked=False,
        allow_crap_scf=False,
        verbosity="info",
        keywords=None,
    ):
        self.system = system
        self.level = int(level)
        self.method = method
        self.basis = basis
        self.functional = functional
        self.driver = driver
        self.cutoffs = dict(cutoffs) if cutoffs else None
        self.unchecked = bool(unchecked)
        self.allow_crap_scf = bool(allow_crap_scf)
        self.verbosity = verbosity
        self.keywords = dict(keywords) if keywords else {}
        self._terms = None

    # -- the settings document ---------------------------------------------

    def settings(self):
        """The JSON that will be sent. Print it when a run surprises you."""
        model = {"method": self.method}
        if self.basis:
            model["basis"] = self.basis
        if self.functional:
            model["functional"] = self.functional

        scf = {}
        if self.allow_crap_scf:
            # Off unless asked for. A non-converged SCF yields a number of the
            # right magnitude, so the run must stop rather than quietly fold it
            # into the total; this says "I know, keep going, tell me at the end".
            scf["allow_crap_scf"] = True

        fragmentation = {"method": "MBE", "level": self.level}
        if self.cutoffs:
            fragmentation["cutoff_method"] = "distance"
            fragmentation["cutoffs"] = self.cutoffs

        document = {
            "schema": {"name": "mqc-python", "version": "1.0"},
            "model": model,
            "driver": self.driver,
            "keywords": {"fragmentation": fragmentation, **({"scf": scf} if scf else {})},
            "system": {
                "logger": {"level": self.verbosity},
                "unchecked_input": self.unchecked,
            },
        }
        _merge(document["keywords"], self.keywords)
        return document

    # -- the term list ------------------------------------------------------

    def terms(self):
        """The terms this expansion will compute: tuples of 1-based monomers.

        Indices only -- no energies, because none have been computed yet. For
        the rows a finished run produced, see `Result.breakdown`.

        Generated on demand and cached, so screening is read-modify-write on
        one list rather than a new one each time.
        """
        if self._terms is None:
            handle = _ffi.fraglist_new()
            try:
                _check(
                    _ffi.fraglist_generate(handle, self.system.n_monomers, self.level),
                    _ffi.fraglist_last_error, True,
                )
                _check(_ffi.fraglist_close_subsets(handle), _ffi.fraglist_last_error, True)
                self._terms = _read_terms(handle)
            finally:
                _ffi.fraglist_free(handle)
        return list(self._terms)

    def keep(self, predicate):
        """Screen the term list, then close it under subsets.

        The predicate is handed one term at a time as a tuple of 1-based
        monomer indices -- `len(t)` is its level -- so a screen against a
        previous run reads that run's breakdown into a set and tests
        membership.

        The closure is not optional and not a convenience. An n-body term's
        delta is its energy less the delta of every proper subset, so a
        screen that keeps a trimer and drops one of its dimers has not
        approximated the expansion -- it has made one that cannot be
        assembled. Dropped subsets are put back, so the list that runs may be
        longer than the one the predicate chose. `run` refuses a list that is
        not closed, so leaving this out fails loudly rather than quietly.
        """
        kept = [t for t in self.terms() if predicate(t)]
        if not kept:
            raise MQCError("the screen kept nothing")
        handle = _ffi.fraglist_new()
        try:
            width = max(len(t) for t in kept)
            packed = []
            for term in kept:
                packed.extend(list(term) + [0] * (width - len(term)))
            _check(
                _ffi.fraglist_set(handle, _ffi.ints(packed), len(kept), width),
                _ffi.fraglist_last_error, True,
            )
            _check(_ffi.fraglist_close_subsets(handle), _ffi.fraglist_last_error, True)
            self._terms = _read_terms(handle)
        finally:
            _ffi.fraglist_free(handle)
        return self

    def set_terms(self, terms):
        """Hand over a list outright -- a restart, or a screen done elsewhere."""
        self._terms = None
        return self.keep(lambda t: t in {tuple(x) for x in terms})

    # -- running ------------------------------------------------------------

    def run(self, label="mqc", write_to_file=True):
        """Compute, on every rank of the job.

        Returns the energy and, with `write_to_file`, leaves the JSON summary
        and the per-fragment CSV behind. The two are independent: a screening
        pass wants the number to decide with *and* the breakdown to screen on.
        """
        payload = json.dumps(self.settings()).encode()
        name = _check_label(label).encode()

        handle = None
        try:
            if self._terms is not None:
                handle = _ffi.fraglist_new()
                width = max(len(t) for t in self._terms)
                packed = []
                for term in self._terms:
                    packed.extend(list(term) + [0] * (width - len(term)))
                _check(
                    _ffi.fraglist_set(handle, _ffi.ints(packed), len(self._terms), width),
                    _ffi.fraglist_last_error, True,
                )

            energy = ctypes.c_double(0.0)
            status = _ffi.run(
                self.system._handle, handle or None,
                len(payload), payload, len(name), name,
                1 if write_to_file else 0, ctypes.byref(energy),
            )
            _check(status, _ffi.run_last_error)
            return Result(energy.value, str(label), bool(write_to_file))
        finally:
            if handle:
                _ffi.fraglist_free(handle)

    def __repr__(self):
        n = "?" if self._terms is None else len(self._terms)
        return f"<mqc.MBE level={self.level} {self.method} terms={n}>"


def _check_label(label):
    """Reject a label the output files will not come back under.

    The label names the output, and Fortran treats it as a *filename*: it
    strips any directory and everything after the last dot, the way it derives
    `output_water.json` from `water.json`. So a label like "tau0.001" writes
    `output_tau0_*`, and "tau0.002" writes over it -- two runs, one set of
    files, no complaint from anything. A later `breakdown()` then reads
    whichever run finished last while believing it has the one it asked for.

    Refused here rather than silently rewritten, because a label the caller
    chose and a file they cannot find is a better failure than a file they can
    find with the wrong numbers in it.
    """
    name = str(label)
    if not name:
        raise ValueError("a run needs a label; it names the output files")
    stem = name.rsplit("/", 1)[-1]
    if "." in stem:
        stem = stem.rsplit(".", 1)[0]
    if stem != name:
        raise ValueError(
            f"label {name!r} would write its output as {stem or '(nothing)'!r} -- "
            "the label names a file, so dots and slashes are eaten. Use "
            f"{name.replace('/', '_').replace('.', '_')!r} instead."
        )
    return name


def _converged(cell):
    """Map the table's scf column back to a tri-state."""
    text = cell.strip().lower()
    if text == "yes":
        return True
    if text == "no":
        return False
    return None


def _read_terms(handle):
    count = _ffi.fraglist_count(handle)
    width = _ffi.fraglist_max_level(handle)
    buf = (ctypes.c_int * (count * width))()
    _check(_ffi.fraglist_get(handle, buf, count, width), _ffi.fraglist_last_error, True)
    return [
        tuple(m for m in buf[i * width : (i + 1) * width] if m > 0) for i in range(count)
    ]


def _merge(into, extra):
    """Deep-merge `extra` over `into`, so a keywords dict can override one key."""
    for key, value in extra.items():
        if isinstance(value, dict) and isinstance(into.get(key), dict):
            _merge(into[key], value)
        else:
            into[key] = value
