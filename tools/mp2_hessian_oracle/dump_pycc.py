#!/usr/bin/env python3
"""Dump pycc's MP2 Hessian intermediates -- the oracle the metalquicha
correlation Hessian is gated against, unit by unit.

The design plan is ``mp2-hessian-plan.md``; the execution ladder that consumes
these dumps is ``mp2-hessian-phased-plan.md``, Unit 0.2. The point of dumping
*intermediates* rather than only the final Hessian is that a per-stage
comparison says "wrong here" where an end-of-phase comparison only says "wrong
somewhere".

Four gated configurations -- {all-electron, frozen-core} x {symmetric,
asymmetric geometry} -- plus two runs on psi4's own basis tables purely to
prove the harness drives pycc correctly.

**Basis provenance is not a detail.** Psi4's internal Pople tables carry eight
significant figures where the Basis Set Exchange bundle under ``basis_sets/``
carries ten (O exponent ``5484.6717`` against ``5484.67166``). Measured on the
pinned case, that difference alone moves the SCF energy by 7e-9 and the
correlation Hessian by up to 3e-9 -- over every 1e-10 gate by more than an
order, and it looks exactly like a persistent bug in our Fortran. So the dumps
our code gates against are generated with psi4 fed *the same file metalquicha
reads*, ``basis_sets/6-31g.json``, and the frozen ``HESS_COL`` columns are
regenerated on that basis too. The published pycc literals were computed on
psi4's internal tables, so they are checked on the internal tables and used for
nothing else.

Four of the quantities below (``fX``, ``wx``, the pass-1 skeleton scalar ``s``,
the pass-2 response ``resp``) are locals inside ``CorrelatedDerivs._hessian_blocks``
and cannot be pulled out. They are re-derived here from the reachable pieces by
transcribing pycc's own formulas -- and ``--check`` proves that transcription is
pycc's and not almost-pycc's by reassembling the full correlation block from the
re-derived parts and requiring it to equal ``_hessian_blocks()`` bit for bit.

Usage::

    source ~/dev/mqc_worktrees/mqc_env.sh
    tools/mp2_hessian_oracle/dump_pycc.py --out <dir>

Requires a configured build tree, because ``basis_sets/6-31g.json`` is extracted
from the BSE bundle at configure time.
"""

import argparse
import datetime
import hashlib
import json
import pathlib
import sys

import numpy as np

# --------------------------------------------------------------------------
# The pinned case. Bohr, frame pinned: a displacement must not move the centre
# of mass or re-orient, or it will not match analytic (bohr) integral
# derivatives. GEOM_SYM is pycc's own test_069 geometry, which is what the
# published HESS_COL literals were computed at. GEOM_ASYM breaks the residual
# planarity, because a symmetric geometry can hide a bug that an asymmetric one
# exposes -- pycc record a CCSD ket-swap bug masked exactly that way.
# --------------------------------------------------------------------------

SYMBOLS = ["O", "H", "H"]

GEOM_SYM = np.array(
    [
        [0.0, 0.000000, 0.000000],
        [0.0, 0.000000, 1.814137],
        [0.0, 1.756000, -0.454300],
    ]
)

GEOM_ASYM = np.array(
    [
        [0.000000, 0.000000, 0.000000],
        [0.031000, -0.024000, 1.814137],
        [0.000000, 1.756000, -0.454300],
    ]
)

BASIS_NAME = "6-31g"
ELEMENT_Z = {"H": 1, "O": 8}
ANGULAR = "SPDFGH"


def sha(array):
    """Content hash of an array, for the manifest."""
    return hashlib.sha256(np.ascontiguousarray(array, dtype=np.float64).tobytes()).hexdigest()[:16]


def geometry_block(coords):
    lines = ["units bohr", "symmetry c1", "no_com", "no_reorient"]
    for symbol, xyz in zip(SYMBOLS, coords):
        lines.append(f"{symbol} {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f}")
    return "\n".join(lines) + "\n"


def bse_gbs(basis_json):
    """Convert the repository's BSE JSON to a psi4 ``.gbs`` string.

    One gbs shell block per ``electron_shells`` entry, with one coefficient
    column per contraction -- which is how psi4 spells an SP shell, and is the
    only general-contraction form 6-31G needs. A basis whose shells carry more
    columns than angular momenta (cc-pVDZ's contracted S, say) does not survive
    this mapping, so it is refused rather than silently mistranscribed.
    """
    data = json.loads(pathlib.Path(basis_json).read_text())
    lines = ["cartesian", "****"]
    # Once per *element*, not once per atom. Emitting a block twice does not
    # fail: psi4 concatenates it, the extra functions are linearly dependent,
    # and the SCF drops them again -- so the span, the energy and every
    # MO-basis quantity come out right while nao is silently wrong and every
    # AO-basis array is a different size. Ask this file for unique symbols.
    for symbol in dict.fromkeys(SYMBOLS):
        z = str(ELEMENT_Z[symbol])
        lines.append(f"{symbol}     0")
        for shell in data["elements"][z]["electron_shells"]:
            momenta = shell["angular_momentum"]
            columns = shell["coefficients"]
            if len(momenta) not in (1, len(columns)):
                raise SystemExit(
                    f"{basis_json}: shell with {len(momenta)} angular momenta and "
                    f"{len(columns)} contractions has no gbs spelling"
                )
            if len(momenta) == 1 and len(columns) > 1:
                raise SystemExit(
                    f"{basis_json}: general contraction on l={momenta[0]} needs one "
                    "gbs shell per column; not implemented (6-31G does not need it)"
                )
            label = "".join(ANGULAR[a] for a in momenta)
            lines.append(f"{label}   {len(shell['exponents'])}   1.00")
            for i, exponent in enumerate(shell["exponents"]):
                row = f"      {float(exponent): .12E}"
                for column in columns:
                    row += f"  {float(column[i]): .12E}"
                lines.append(row)
        lines.append("****")
    return "\n".join(lines) + "\n"


def make_driver(psi4, pycc, coords, freeze_core, gbs=None):
    """An MPderiv driver at one geometry, basis source and core treatment."""
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.geometry(geometry_block(coords))
    if gbs is None:
        psi4.set_options({"basis": BASIS_NAME.upper()})
    else:
        psi4.basis_helper(gbs, name="mqcbse")
    psi4.set_options(
        {
            "scf_type": "pk",
            "freeze_core": freeze_core,
            "e_convergence": 1e-13,
            "d_convergence": 1e-13,
        }
    )
    energy, wfn = psi4.energy("scf", return_wfn=True)
    mp = pycc.MPwfn(wfn, orbital_basis="spatial")
    mp.compute_energy()
    return energy, pycc.MPderiv(mp)


def collect(drv, Perturbation):
    """Every intermediate on the ladder, in dependency order.

    Everything here except ``fX``/``wx`` is a direct pull from pycc. The stacked
    per-perturbation arrays are indexed ``[3*atom + cartesian]``.
    """
    wfn, c, d = drv.wfn, drv.contract, drv.wfn.derivatives
    v, ofull = wfn.v, slice(0, wfn.o.stop)
    ncore = wfn.o.stop - wfn.no
    natom = d.natom
    canonical = drv.perturbed_mo_gauge == "canonical"

    rec = drv._orbital_response()
    Drel, Gam = rec.Drel, rec.Gam
    out = {
        "t2": np.asarray(drv.mp.t2),
        "Drel": np.asarray(Drel),
        "Gam": np.asarray(Gam),
        "I": np.asarray(drv._lagrangian(Drel, Gam)),
        "GeffAO": np.asarray(drv._effective_2pdm_ao(Drel, Gam)),
        "C": np.asarray(wfn.C),
        "eps": np.diag(np.asarray(wfn.H.F)),
    }
    Doo, Dvv = drv._mp2_opdm()
    out["Doo"], out["Dvv"] = np.asarray(Doo), np.asarray(Dvv)

    pert = [Perturbation("nuclear", (a, ct)) for a in range(natom) for ct in range(3)]
    cphf = drv._full_occ_cphf()
    stacks = {k: [] for k in ("fX", "SX", "erix", "Ip", "Xx", "I2x", "U", "dDrel", "dI", "dGam")}
    pf_x = []
    for p in pert:
        atom, ct = p.comp
        erix = np.asarray(d.eri(atom)[ct])
        wx = 2.0 * erix - erix.transpose(0, 1, 3, 2)
        hx = np.asarray(d.core(atom)[ct])
        Sx = np.asarray(d.overlap(atom)[ct])
        # f^(X) = h^(X) + L^(X) traced over the FULL occupied space (core included).
        fx = hx + c("pmqm->pq", wx[:, ofull, :, ofull])
        Ip, xov, i2 = drv._skeleton_lagrangian(fx, Sx, wx, erix, Drel, Gam, out["I"])
        if canonical or ncore:
            xt, it, pf = drv._augment_with_canonical_pair_rotations(Ip, xov, i2)
        else:
            xt, it, pf = xov, i2, None
        response = drv._relaxed_response(p)
        stacks["fX"].append(fx)
        stacks["SX"].append(Sx)
        stacks["erix"].append(erix)
        stacks["Ip"].append(np.asarray(Ip))
        stacks["Xx"].append(np.asarray(xt))
        stacks["I2x"].append(np.asarray(it))
        stacks["U"].append(np.asarray(cphf.full_U(p, ncore, canonical=canonical)))
        stacks["dDrel"].append(np.asarray(response.dDrel))
        stacks["dI"].append(np.asarray(response.dI))
        stacks["dGam"].append(np.asarray(response.dGam))
        pf_x.append(pf)
    for key, value in stacks.items():
        out[key] = np.stack(value)
    if any(p is not None for p in pf_x):
        out["Pf_x"] = np.stack([np.asarray(p) for p in pf_x])
    out["dt2_pert2"] = np.asarray(drv._perturbed_t2(pert[2]))

    reference, correlation = drv._hessian_blocks()
    out["H_correlation"] = np.asarray(correlation)
    out["H_reference"] = np.asarray(reference)

    # The two pass-scalars that are locals in pycc, re-derived and kept per
    # component so Units 1.3 and 1.6 can gate against them directly.
    out["skel_s"], out["resp"] = _passes(drv, out, pert)
    return out, ncore, canonical


def _passes(drv, arr, pert):
    """Re-derive the pass-1 skeleton scalar and the pass-2 response, per (X, Y).

    Transcribed from ``_hessian_blocks``' ``'aod'`` route. ``--check`` is what
    licenses this: the two together must reassemble the correlation block
    exactly.
    """
    wfn, c, d = drv.wfn, drv.contract, drv.wfn.derivatives
    v, ofull = wfn.v, slice(0, wfn.o.stop)
    natom = d.natom
    nc = 3 * natom
    Drel, I, GeffAO = arr["Drel"], arr["I"], arr["GeffAO"]
    fX, SX, Xx, I2x, U = arr["fX"], arr["SX"], arr["Xx"], arr["I2x"], arr["U"]
    dDrel, dI, dGam = arr["dDrel"], arr["dI"], arr["dGam"]
    Pf_x = arr.get("Pf_x")

    skel = np.zeros((nc, nc))
    for a1 in range(natom):
        for a2 in range(a1, natom):
            core2s = [np.asarray(m) for m in d.core2(a1, a2)]
            ov2s = [np.asarray(m) for m in d.overlap2(a1, a2)]
            ao_eri = d.ao_eri2(a1, a2)
            for cx in range(3):
                for cy in range(3):
                    comp = cx * 3 + cy
                    s = (
                        c("mnls,mnls->", GeffAO, ao_eri[comp])
                        + c("pq,pq->", Drel, core2s[comp])
                        + c("pq,pq->", I, ov2s[comp])
                    )
                    ix, iy = a1 * 3 + cx, a2 * 3 + cy
                    skel[ix, iy] += s
                    if a1 != a2:
                        skel[iy, ix] += s

    resp = np.zeros((nc, nc))
    for atom in range(natom):
        erX = [np.asarray(m) for m in d.eri(atom)]
        for iy in range(nc):
            for cx in range(3):
                ix = atom * 3 + cx
                orb = 2.0 * c("ai,ai->", U[iy][v, ofull], Xx[ix][v, ofull]) + c(
                    "pq,pq->", SX[iy], I2x[ix]
                )
                if Pf_x is not None:
                    orb += c("pq,pq->", Pf_x[ix], fX[iy])
                resp[ix, iy] = (
                    c("pq,pq->", dDrel[iy], fX[ix])
                    + c("pqrs,pqrs->", dGam[iy], erX[cx])
                    + c("pq,pq->", dI[iy], SX[ix])
                    + orb
                )
    return skel, resp


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="oracle_dumps", help="output directory")
    parser.add_argument(
        "--repo",
        default=str(pathlib.Path(__file__).resolve().parents[2]),
        help="metalquicha worktree (for basis_sets/6-31g.json)",
    )
    parser.add_argument("--no-check", action="store_true", help="skip the Gate 0.2 checks")
    args = parser.parse_args()

    import psi4
    import pycc
    from pycc.cphf import Perturbation

    psi4.core.be_quiet()
    psi4.set_memory("2 GB")

    basis_json = pathlib.Path(args.repo) / "basis_sets" / f"{BASIS_NAME}.json"
    if not basis_json.exists():
        raise SystemExit(
            f"{basis_json} not found -- configure the build first; the BSE bundle is "
            "unpacked at configure time"
        )
    gbs = bse_gbs(basis_json)

    out_root = pathlib.Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)
    manifest = {
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "psi4": psi4.__version__,
        "pycc": pycc.__version__,
        "numpy": np.__version__,
        "basis": {
            "name": BASIS_NAME,
            "bse_json": str(basis_json),
            "bse_sha256": hashlib.sha256(basis_json.read_bytes()).hexdigest()[:16],
            "note": "dumps use the BSE file metalquicha reads; 'internal' runs exist "
            "only to check the published literals",
        },
        "configurations": {},
        "hess_col_bse": {},
        "gate_0_2": {},
    }

    failures = []
    for geometry_name, coords in (("sym", GEOM_SYM), ("asym", GEOM_ASYM)):
        for fc_name, fc in (("ae", "false"), ("fc", "true")):
            tag = f"{geometry_name}_{fc_name}"
            energy, drv = make_driver(psi4, pycc, coords, fc, gbs=gbs)
            arrays, ncore, canonical = collect(drv, Perturbation)

            target = out_root / tag
            target.mkdir(exist_ok=True)
            entries = {}
            for name, array in arrays.items():
                np.save(target / f"{name}.npy", array)
                entries[name] = {
                    "shape": list(array.shape),
                    "sha256": sha(array),
                }
            manifest["configurations"][tag] = {
                "geometry": coords.tolist(),
                "freeze_core": fc,
                "ncore": int(ncore),
                "gauge": drv.perturbed_mo_gauge,
                "canonical": bool(canonical),
                "scf_energy": float(energy),
                "basis_source": "bse",
                "arrays": entries,
            }

            H = arrays["H_correlation"]
            checks = {
                "reconstruction": float(
                    np.max(np.abs(arrays["skel_s"] + arrays["resp"] - H))
                ),
                "symmetry": float(np.max(np.abs(H - H.T))),
                "translational": float(
                    np.max(np.abs(H.reshape(3, 3, 3, 3).sum(axis=2)))
                ),
            }
            manifest["gate_0_2"][tag] = checks
            if checks["reconstruction"] != 0.0:
                failures.append(f"{tag}: reconstruction {checks['reconstruction']:.3e} != 0")
            print(
                f"[{tag:8s}] ncore={ncore} scf={energy:.12f}  "
                f"reconstruction={checks['reconstruction']:.1e}  "
                f"H=H^T {checks['symmetry']:.1e}  transl {checks['translational']:.1e}"
            )

            # Decision 2026-08-27: regenerate the frozen columns on OUR basis, so
            # no Phase 1 or Phase 2 gate ever reads psi4's internal tables.
            if geometry_name == "sym":
                manifest["hess_col_bse"][fc_name] = {
                    f"{atom},{cart}": H[:, atom * 3 + cart].tolist()
                    for atom in range(3)
                    for cart in range(3)
                }

    if not args.no_check:
        # The literals were computed on psi4's internal tables; check them there
        # and nowhere else.
        import pycc.tests.test_069_mp2_hessian as ref

        for fc_name, fc, literals in (
            ("ae", "false", ref.HESS_COL_631G),
            ("fc", "true", ref.HESS_COL_631G_FC),
        ):
            _, drv = make_driver(psi4, pycc, GEOM_SYM, fc, gbs=None)
            _, correlation = drv._hessian_blocks()
            correlation = np.asarray(correlation)
            worst = 0.0
            for (atom, cart), column in literals.items():
                worst = max(worst, float(np.max(np.abs(correlation[:, atom * 3 + cart] - column))))
            manifest["gate_0_2"][f"literals_{fc_name}_internal_basis"] = worst
            print(f"[literal {fc_name}] internal-basis max deviation {worst:.3e}")
            if worst > 1e-10:
                failures.append(f"literals {fc_name}: {worst:.3e} > 1e-10")

    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nwrote {out_root}/manifest.json")

    if failures:
        print("\nGATE 0.2 FAILED:")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print("GATE 0.2 GREEN")
    return 0


if __name__ == "__main__":
    sys.exit(main())
