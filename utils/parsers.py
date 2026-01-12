"""JSON and XYZ parsing for MQC input files."""

from __future__ import annotations
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

try:
    import numpy as np
except Exception:
    np = None

from .validation import die, require_only_keys, req_type, opt_type, check_index_0based
from .models import (
    SchemaTag, Model, XTB, SCF, Hessian, AIMD, FragCutoffs, Fragmentation,
    Logger, System, Molecule, Input,
    TBLITE_SOLVENTS, SOLVATION_MODELS,
)


def read_xyz_file(xyz_path: Union[str, Path]) -> Tuple[List[str], Any]:
    """Read an XYZ file and return (symbols, geometry)."""
    xyz_path = Path(xyz_path)
    lines = xyz_path.read_text(encoding="utf-8").strip().splitlines()
    if len(lines) < 3:
        die(f"XYZ file too short: {xyz_path}")

    try:
        natom = int(lines[0].split()[0])
    except Exception as e:
        die(f"XYZ first line must start with natom int: {xyz_path} ({e})")
    if natom <= 0:
        die(f"XYZ natom must be > 0: {xyz_path}")

    atom_lines = lines[2:]
    if len(atom_lines) < natom:
        die(f"XYZ expected {natom} atom lines, got {len(atom_lines)}: {xyz_path}")

    symbols: List[str] = []
    coords: List[float] = []
    for i in range(natom):
        parts = atom_lines[i].split()
        if len(parts) < 4:
            die(f"XYZ atom line malformed at line {i+3}: '{atom_lines[i]}' in {xyz_path}")
        sym = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except Exception as e:
            die(f"XYZ atom line has non-numeric coords at line {i+3}: '{atom_lines[i]}' ({e})")
        symbols.append(sym)
        coords.extend([x, y, z])

    geom = _reshape_geometry(coords, natom)
    return symbols, geom


def _reshape_geometry(flat: Any, natom: int) -> Any:
    """Reshape flat coordinate list to Nx3 array."""
    req_type(flat, list, "molecule.geometry")
    if len(flat) != 3 * natom:
        die(f"molecule.geometry: expected length {3*natom}, got {len(flat)}")

    if np is not None:
        return np.asarray(flat, dtype=float).reshape((natom, 3))

    out = []
    for i in range(natom):
        x = float(flat[3*i + 0])
        y = float(flat[3*i + 1])
        z = float(flat[3*i + 2])
        out.append([x, y, z])
    return out


def _parse_schema(obj: Any) -> SchemaTag:
    """Parse schema from string or dict."""
    if isinstance(obj, str):
        if "-" in obj:
            name, version = obj.rsplit("-", 1)
        else:
            name, version = obj, "unknown"
        return SchemaTag(name=name, version=version, index_base=0, units="angstrom")

    if isinstance(obj, dict):
        require_only_keys(obj, {"name", "version", "index_base", "units"}, "schema")
        name = req_type(obj.get("name"), str, "schema.name")
        version = req_type(obj.get("version"), str, "schema.version")
        index_base = obj.get("index_base", 0)
        units = obj.get("units", "angstrom")
        if index_base != 0:
            die("schema.index_base: only 0 supported in v1")
        units = req_type(units, str, "schema.units")
        return SchemaTag(name=name, version=version, index_base=0, units=units)

    die("schema must be string or object")


def _parse_model(d: Dict[str, Any]) -> Model:
    """Parse model section."""
    require_only_keys(d, {"method", "basis", "aux_basis"}, "model")
    method = req_type(d.get("method"), str, "model.method")
    basis = opt_type(d.get("basis"), str, "model.basis")
    aux = opt_type(d.get("aux_basis"), str, "model.aux_basis")

    return Model(method=method, basis=basis, aux_basis=aux)


def _parse_xtb_keywords(d: Dict[str, Any]) -> XTB:
    """Parse keywords.xtb section for XTB-specific settings."""
    require_only_keys(d, {"solvent", "solvation_model", "dielectric", "cpcm_nang", "cpcm_rscale"},
                      "keywords.xtb")

    solvent = opt_type(d.get("solvent"), str, "keywords.xtb.solvent")
    solvation_model = opt_type(d.get("solvation_model"), str, "keywords.xtb.solvation_model")

    # Parse CPCM-specific settings
    dielectric = d.get("dielectric")
    if dielectric is not None:
        if not isinstance(dielectric, (int, float)):
            die("keywords.xtb.dielectric must be a number")
        dielectric = float(dielectric)
        if dielectric <= 0:
            die("keywords.xtb.dielectric must be > 0")

    cpcm_nang = d.get("cpcm_nang")
    if cpcm_nang is not None:
        if not isinstance(cpcm_nang, int):
            die("keywords.xtb.cpcm_nang must be an integer")
        if cpcm_nang <= 0:
            die("keywords.xtb.cpcm_nang must be > 0")

    cpcm_rscale = d.get("cpcm_rscale")
    if cpcm_rscale is not None:
        if not isinstance(cpcm_rscale, (int, float)):
            die("keywords.xtb.cpcm_rscale must be a number")
        cpcm_rscale = float(cpcm_rscale)
        if cpcm_rscale <= 0 or cpcm_rscale > 2.0:
            die("keywords.xtb.cpcm_rscale must be in range (0, 2.0]")

    # Validate solvent if specified
    if solvent is not None:
        solvent_lower = solvent.lower()
        if solvent_lower not in TBLITE_SOLVENTS:
            die(f"keywords.xtb.solvent: unknown solvent '{solvent}'. "
                f"Supported solvents include: water, methanol, ethanol, acetone, acetonitrile, "
                f"benzene, toluene, chloroform, dichloromethane, dmso, dmf, thf, "
                f"diethylether, hexane, cyclohexane, and many more. "
                f"See documentation for full list.")
        solvent = solvent_lower

    # Validate solvation model if specified
    if solvation_model is not None:
        solvation_model_lower = solvation_model.lower()
        if solvation_model_lower not in SOLVATION_MODELS:
            die(f"keywords.xtb.solvation_model: unknown model '{solvation_model}'. "
                f"Supported models: {', '.join(sorted(SOLVATION_MODELS))}")
        solvation_model = solvation_model_lower

    # If dielectric is specified without solvation_model, default to cpcm
    if dielectric is not None and solvation_model is None:
        solvation_model = "cpcm"

    # If solvent is specified but model isn't, default to alpb
    if solvent is not None and solvation_model is None:
        solvation_model = "alpb"

    # Validation: need either solvent or dielectric for solvation
    if solvation_model is not None and solvent is None and dielectric is None:
        die("keywords.xtb.solvation_model: cannot specify solvation model without solvent or dielectric")

    # CPCM-specific validation
    if solvation_model == "cpcm":
        if solvent is None and dielectric is None:
            die("keywords.xtb: CPCM solvation requires either 'solvent' or 'dielectric'")
    else:
        if solvation_model is not None and solvent is None:
            die(f"keywords.xtb: {solvation_model.upper()} solvation requires a solvent name")

    return XTB(solvent=solvent, solvation_model=solvation_model,
               dielectric=dielectric, cpcm_nang=cpcm_nang, cpcm_rscale=cpcm_rscale)


def _parse_keywords(d: Dict[str, Any]) -> Tuple[Optional[SCF], Optional[Hessian], Optional[AIMD], Optional[Fragmentation], Optional[XTB]]:
    """Parse keywords section."""
    require_only_keys(d, {"scf", "hessian", "aimd", "fragmentation", "xtb"}, "keywords")

    scf = None
    if "scf" in d:
        sd = req_type(d["scf"], dict, "keywords.scf")
        require_only_keys(sd, {"maxiter", "tolerance"}, "keywords.scf")
        maxiter = req_type(sd.get("maxiter"), int, "keywords.scf.maxiter")
        tol = sd.get("tolerance")
        if not isinstance(tol, (int, float)):
            die("keywords.scf.tolerance must be number")
        if maxiter <= 0:
            die("keywords.scf.maxiter must be > 0")
        if float(tol) <= 0:
            die("keywords.scf.tolerance must be > 0")
        scf = SCF(maxiter=maxiter, tolerance=float(tol))

    hessian = None
    if "hessian" in d:
        hd = req_type(d["hessian"], dict, "keywords.hessian")
        require_only_keys(hd, {"finite_difference_displacement", "displacement", "temperature", "pressure"}, "keywords.hessian")
        disp = hd.get("finite_difference_displacement") or hd.get("displacement")
        if disp is None:
            disp = 0.001
        if not isinstance(disp, (int, float)):
            die("keywords.hessian.finite_difference_displacement must be number")
        if float(disp) <= 0:
            die("keywords.hessian.finite_difference_displacement must be > 0")

        temperature = hd.get("temperature", 298.15)
        if not isinstance(temperature, (int, float)):
            die("keywords.hessian.temperature must be number")
        if float(temperature) <= 0:
            die("keywords.hessian.temperature must be > 0")

        pressure = hd.get("pressure", 1.0)
        if not isinstance(pressure, (int, float)):
            die("keywords.hessian.pressure must be number")
        if float(pressure) <= 0:
            die("keywords.hessian.pressure must be > 0")

        hessian = Hessian(finite_difference_displacement=float(disp),
                          temperature=float(temperature), pressure=float(pressure))

    aimd_obj = None
    if "aimd" in d:
        ad = req_type(d["aimd"], dict, "keywords.aimd")
        require_only_keys(ad, {"dt", "timestep", "nsteps", "steps", "initial_temperature",
                               "temperature", "output_frequency", "output_freq"}, "keywords.aimd")
        dt = ad.get("dt") or ad.get("timestep")
        if dt is None:
            dt = 1.0
        nsteps = ad.get("nsteps") or ad.get("steps")
        if nsteps is None:
            nsteps = 0
        init_temp = ad.get("initial_temperature") or ad.get("temperature")
        if init_temp is None:
            init_temp = 300.0
        out_freq = ad.get("output_frequency") or ad.get("output_freq")
        if out_freq is None:
            out_freq = 1

        if not isinstance(dt, (int, float)):
            die("keywords.aimd.dt must be number")
        if not isinstance(nsteps, int):
            die("keywords.aimd.nsteps must be int")
        if not isinstance(init_temp, (int, float)):
            die("keywords.aimd.initial_temperature must be number")
        if not isinstance(out_freq, int):
            die("keywords.aimd.output_frequency must be int")
        if float(dt) <= 0:
            die("keywords.aimd.dt must be > 0")
        if nsteps < 0:
            die("keywords.aimd.nsteps must be >= 0")
        if float(init_temp) <= 0:
            die("keywords.aimd.initial_temperature must be > 0")
        if out_freq <= 0:
            die("keywords.aimd.output_frequency must be > 0")

        aimd_obj = AIMD(dt=float(dt), nsteps=nsteps,
                        initial_temperature=float(init_temp), output_frequency=out_freq)

    frag = None
    if "fragmentation" in d:
        fd = req_type(d["fragmentation"], dict, "keywords.fragmentation")
        require_only_keys(
            fd,
            {"method", "allow_overlapping_fragments", "level", "embedding",
             "cutoff_method", "distance_metric", "cutoffs"},
            "keywords.fragmentation",
        )

        method = req_type(fd.get("method"), str, "keywords.fragmentation.method")
        allow = req_type(fd.get("allow_overlapping_fragments"), bool,
                         "keywords.fragmentation.allow_overlapping_fragments")
        level = req_type(fd.get("level"), int, "keywords.fragmentation.level")
        embedding = req_type(fd.get("embedding"), str, "keywords.fragmentation.embedding")
        cutoff_method = req_type(fd.get("cutoff_method"), str, "keywords.fragmentation.cutoff_method")
        distance_metric = req_type(fd.get("distance_metric"), str, "keywords.fragmentation.distance_metric")

        if level <= 0:
            die("keywords.fragmentation.level must be > 0")

        cutoffs = None
        if "cutoffs" in fd:
            cd = req_type(fd["cutoffs"], dict, "keywords.fragmentation.cutoffs")
            cutoff_dict: Dict[int, float] = {}

            name_to_level = {
                "dimer": 2, "trimer": 3, "tetramer": 4, "pentamer": 5,
                "hexamer": 6, "heptamer": 7, "octamer": 8
            }

            for key, value in cd.items():
                if key in name_to_level:
                    nmer_level = name_to_level[key]
                elif isinstance(key, int):
                    nmer_level = key
                elif key.isdigit():
                    nmer_level = int(key)
                else:
                    die(f"keywords.fragmentation.cutoffs: unknown key '{key}'. "
                        f"Use 'dimer', 'trimer', etc. or numeric keys like '2', '3', ...")

                if nmer_level < 2:
                    die(f"keywords.fragmentation.cutoffs: n-mer level must be >= 2, got {nmer_level}")

                if not isinstance(value, (int, float)) or float(value) <= 0:
                    die(f"keywords.fragmentation.cutoffs.{key} must be a positive number")

                cutoff_dict[nmer_level] = float(value)

            if cutoff_dict:
                sorted_levels = sorted(cutoff_dict.keys())
                for i in range(len(sorted_levels) - 1):
                    level_low = sorted_levels[i]
                    level_high = sorted_levels[i + 1]
                    cutoff_low = cutoff_dict[level_low]
                    cutoff_high = cutoff_dict[level_high]

                    if cutoff_high > cutoff_low:
                        die(f"keywords.fragmentation.cutoffs: {level_high}-mer cutoff ({cutoff_high}) "
                            f"cannot be larger than {level_low}-mer cutoff ({cutoff_low}). "
                            f"Cutoffs must be monotonically decreasing.")

                cutoffs = FragCutoffs(cutoffs=cutoff_dict)

        frag = Fragmentation(
            method=method,
            allow_overlapping_fragments=allow,
            level=level,
            embedding=embedding,
            cutoff_method=cutoff_method,
            distance_metric=distance_metric,
            cutoffs=cutoffs,
        )

    xtb = None
    if "xtb" in d:
        xd = req_type(d["xtb"], dict, "keywords.xtb")
        xtb = _parse_xtb_keywords(xd)

    return scf, hessian, aimd_obj, frag, xtb


def _parse_system(d: Dict[str, Any]) -> System:
    """Parse system section."""
    require_only_keys(d, {"logger"}, "system")

    logger_data = req_type(d.get("logger"), dict, "system.logger")
    require_only_keys(logger_data, {"level"}, "system.logger")

    level = req_type(logger_data.get("level"), str, "system.logger.level")
    valid_levels = {"debug", "verbose", "info", "performance", "warning", "error", "knowledge"}
    if level.lower() not in valid_levels:
        die(f"system.logger.level must be one of: {', '.join(sorted(valid_levels))}")

    return System(logger=Logger(level=level))


def _parse_molecule(m: Dict[str, Any], base_dir: Path, mi: int) -> Molecule:
    """Parse a single molecule entry."""
    require_only_keys(
        m,
        {"symbols", "geometry", "xyz", "molecular_charge", "molecular_multiplicity",
         "connectivity", "fragments", "fragment_charges", "fragment_multiplicities"},
        f"molecules[{mi}]",
    )

    has_xyz = "xyz" in m
    has_sym = "symbols" in m
    has_geo = "geometry" in m

    if has_xyz:
        if has_sym or has_geo:
            die(f"molecules[{mi}]: provide either xyz OR (symbols+geometry), not both")
        xyz_rel = req_type(m.get("xyz"), str, f"molecules[{mi}].xyz")
        xyz_path = (base_dir / xyz_rel).resolve()
        symbols, geom = read_xyz_file(xyz_path)
        source_xyz = xyz_rel
    else:
        if not (has_sym and has_geo):
            die(f"molecules[{mi}]: must provide either xyz OR (symbols+geometry)")
        symbols = req_type(m.get("symbols"), list, f"molecules[{mi}].symbols")
        if not symbols or not all(isinstance(s, str) for s in symbols):
            die(f"molecules[{mi}].symbols must be a non-empty list of strings")
        natom = len(symbols)
        geom = _reshape_geometry(m.get("geometry"), natom)
        source_xyz = None

    natom = len(symbols)

    charge = req_type(m.get("molecular_charge"), int, f"molecules[{mi}].molecular_charge")
    multiplicity = req_type(m.get("molecular_multiplicity"), int, f"molecules[{mi}].molecular_multiplicity")
    if multiplicity <= 0:
        die(f"molecules[{mi}].molecular_multiplicity must be >= 1")

    frags_raw = m.get("fragments", None)

    if frags_raw is None:
        fragments = []
        charges = []
        mults = []
    else:
        frags = req_type(frags_raw, list, f"molecules[{mi}].fragments")
        if len(frags) == 0:
            die(f"molecules[{mi}].fragments must be non-empty if provided")

        fragments: List[List[int]] = []
        for fi, frag in enumerate(frags):
            req_type(frag, list, f"molecules[{mi}].fragments[{fi}]")
            if len(frag) == 0:
                die(f"molecules[{mi}].fragments[{fi}] empty")
            fragments.append([check_index_0based(a, natom, f"molecules[{mi}].fragments[{fi}]") for a in frag])

        nfrag = len(fragments)

        charges = m.get("fragment_charges", [0] * nfrag)
        mults = m.get("fragment_multiplicities", [1] * nfrag)

        req_type(charges, list, f"molecules[{mi}].fragment_charges")
        req_type(mults, list, f"molecules[{mi}].fragment_multiplicities")
        if len(charges) != nfrag:
            die(f"molecules[{mi}].fragment_charges length mismatch: expected {nfrag}, got {len(charges)}")
        if len(mults) != nfrag:
            die(f"molecules[{mi}].fragment_multiplicities length mismatch: expected {nfrag}, got {len(mults)}")
        if not all(isinstance(c, int) for c in charges):
            die(f"molecules[{mi}].fragment_charges must be ints")
        if not all(isinstance(mm, int) for mm in mults):
            die(f"molecules[{mi}].fragment_multiplicities must be ints")
        for i, mm in enumerate(mults):
            if mm <= 0:
                die(f"molecules[{mi}].fragment_multiplicities[{i}] must be >=1")

        total_frag_charge = sum(charges)
        if total_frag_charge != charge:
            die(f"molecules[{mi}]: fragment charges sum to {total_frag_charge}, "
                f"but molecular charge is {charge}. These must be equal.")

        if all(m == 1 for m in mults) and multiplicity != 1:
            die(f"molecules[{mi}]: all fragments have multiplicity=1 (closed-shell), "
                f"but molecular multiplicity is {multiplicity}. "
                f"This is inconsistent - molecular multiplicity should also be 1.")

        if multiplicity == 1 and not all(m == 1 for m in mults):
            import sys
            print(f"WARNING: molecules[{mi}]: molecular multiplicity=1 but some fragments "
                  f"have multiplicity>1. Ensure fragment spins are properly paired.",
                  file=sys.stderr)

    conn_raw = m.get("connectivity", None)
    connectivity = None
    if conn_raw is not None:
        req_type(conn_raw, list, f"molecules[{mi}].connectivity")
        connectivity = []
        for bi, triple in enumerate(conn_raw):
            req_type(triple, list, f"molecules[{mi}].connectivity[{bi}]")
            if len(triple) != 3:
                die(f"molecules[{mi}].connectivity[{bi}] must be [i,j,order]")
            i, j, order = triple
            i0 = check_index_0based(i, natom, f"molecules[{mi}].connectivity[{bi}].i")
            j0 = check_index_0based(j, natom, f"molecules[{mi}].connectivity[{bi}].j")
            if not isinstance(order, int) or order <= 0:
                die(f"molecules[{mi}].connectivity[{bi}].order must be positive int")
            connectivity.append((i0, j0, order))

    return Molecule(
        symbols=symbols,
        geom_xyz=geom,
        charge=charge,
        multiplicity=multiplicity,
        fragments=fragments,
        fragment_charges=charges,
        fragment_multiplicities=mults,
        connectivity=connectivity,
        source_xyz=source_xyz,
    )


def parse_json_to_input(json_path: Union[str, Path]) -> Input:
    """Parse a JSON input file into an Input object."""
    json_path = Path(json_path)
    base_dir = json_path.parent

    data = json.loads(json_path.read_text(encoding="utf-8"))
    req_type(data, dict, "top-level")

    require_only_keys(data, {"schema", "molecules", "model", "keywords", "driver", "title", "system"}, "top-level")

    schema = _parse_schema(data.get("schema"))
    model = _parse_model(req_type(data.get("model"), dict, "model"))
    driver = req_type(data.get("driver"), str, "driver")
    title = opt_type(data.get("title"), str, "title")

    kw = data.get("keywords", {})
    req_type(kw, dict, "keywords")
    scf, hessian, aimd, frag, xtb = _parse_keywords(kw)

    system = None
    sys_data = data.get("system")
    if sys_data is not None:
        system = _parse_system(req_type(sys_data, dict, "system"))

    mols = req_type(data.get("molecules"), list, "molecules")
    if len(mols) == 0:
        die("molecules must be non-empty list")
    molecules = []
    for mi, m in enumerate(mols):
        molecules.append(_parse_molecule(req_type(m, dict, f"molecules[{mi}]"), base_dir, mi))

    return Input(schema=schema, model=model, driver=driver, molecules=molecules,
                 title=title, scf=scf, hessian=hessian, aimd=aimd, fragmentation=frag, xtb=xtb, system=system)
