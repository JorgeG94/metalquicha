#!/usr/bin/env python3
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

try:
    import numpy as np
except Exception:
    np = None  # optional


# ----------------------------
# Strict validation helpers
# ----------------------------

def die(msg: str) -> None:
    raise ValueError(msg)

def require_only_keys(d: Dict[str, Any], allowed: set[str], ctx: str) -> None:
    extra = set(d.keys()) - allowed
    if extra:
        raise ValueError(f"{ctx}: unknown keys: {sorted(extra)}")

def req_type(x: Any, t, ctx: str) -> Any:
    if not isinstance(x, t):
        die(f"{ctx}: expected {t.__name__}, got {type(x).__name__}")
    return x

def opt_type(x: Any, t, ctx: str) -> Any:
    if x is None:
        return None
    if not isinstance(x, t):
        die(f"{ctx}: expected {t.__name__}, got {type(x).__name__}")
    return x

def check_index_0based(i: Any, natom: int, ctx: str) -> int:
    if not isinstance(i, int):
        die(f"{ctx}: expected int, got {type(i).__name__}")
    if i < 0 or i >= natom:
        die(f"{ctx}: index out of range {i} for natom={natom} (0-based expected)")
    return i


# ----------------------------
# Internal representation
# ----------------------------

@dataclass
class SchemaTag:
    name: str
    version: str
    index_base: int = 0
    units: str = "angstrom"

# ----------------------------
# Supported solvents from tblite
# ----------------------------
# From tblite's src/tblite/solvation/data.f90
TBLITE_SOLVENTS = {
    # Common solvents
    "water", "h2o",
    "methanol", "ch3oh",
    "ethanol", "c2h5oh",
    "acetone",
    "acetonitrile", "ch3cn",
    # Aromatics
    "benzene",
    "toluene",
    # Halogenated
    "chloroform", "chcl3",
    "dichloromethane", "ch2cl2", "dcm",
    "carbon tetrachloride", "ccl4",
    # Polar aprotic
    "dmso", "dimethylsulfoxide",
    "dmf", "dimethylformamide",
    "thf", "tetrahydrofuran",
    # Ethers
    "diethylether", "ether",
    "dioxane",
    # Alkanes
    "hexane", "n-hexane",
    "cyclohexane",
    "heptane", "n-heptane",
    "octane", "n-octane",
    # Other organics
    "pyridine",
    "aniline",
    "nitromethane",
    "nitrobenzene",
    "formamide",
    "phenol",
    "cs2", "carbondisulfide",
    # Alcohols
    "1-propanol", "propanol",
    "2-propanol", "isopropanol",
    "1-butanol", "butanol",
    "2-butanol",
    "1-octanol", "octanol",
    # Esters/acids
    "ethyl acetate", "ethylacetate",
    "acetic acid", "aceticacid",
    "formic acid", "formicacid",
    # Additional solvents
    "chlorobenzene",
    "furan",
    "hexadecane",
    "pentane",
    "decane",
    "decanol",
    "woctanol",  # wet octanol
    "inf",  # infinite dielectric (conductor)
}

# Supported solvation models
SOLVATION_MODELS = {"alpb", "gbsa", "cpcm"}


@dataclass
class Model:
    method: str
    basis: Optional[str] = None
    aux_basis: Optional[str] = None
    solvent: Optional[str] = None
    solvation_model: Optional[str] = None
    # CPCM-specific settings
    dielectric: Optional[float] = None      # Direct dielectric constant (for CPCM)
    cpcm_nang: Optional[int] = None         # Angular grid points for CPCM cavity
    cpcm_rscale: Optional[float] = None     # Radii scaling for CPCM cavity

@dataclass
class SCF:
    maxiter: int
    tolerance: float

@dataclass
class Hessian:
    finite_difference_displacement: float = 0.001  # Bohr

@dataclass
class AIMD:
    dt: float = 1.0  # femtoseconds
    nsteps: int = 0
    initial_temperature: float = 300.0  # Kelvin
    output_frequency: int = 1

@dataclass
class FragCutoffs:
    """
    Distance cutoffs for fragment generation.
    Maps n-mer level (2=dimer, 3=trimer, etc.) to cutoff distance in Angstrom.
    """
    cutoffs: Dict[int, float]  # {2: 5.0, 3: 4.0, ...} for dimer, trimer, etc.

@dataclass
class Fragmentation:
    method: str
    allow_overlapping_fragments: bool
    level: int
    embedding: str
    cutoff_method: str
    distance_metric: str
    cutoffs: Optional[FragCutoffs]

@dataclass
class Logger:
    level: str = "info"

@dataclass
class System:
    logger: Logger

@dataclass
class Molecule:
    symbols: List[str]
    geom_xyz: Any  # Nx3 list-of-lists or numpy array
    charge: int
    multiplicity: int
    fragments: List[List[int]]
    fragment_charges: List[int]
    fragment_multiplicities: List[int]
    connectivity: Optional[List[Tuple[int, int, int]]] = None  # allowed, not emitted v1 (yet)
    source_xyz: Optional[str] = None

@dataclass
class Input:
    schema: SchemaTag
    model: Model
    driver: str
    molecules: List[Molecule]
    title: Optional[str] = None
    scf: Optional[SCF] = None
    hessian: Optional[Hessian] = None
    aimd: Optional[AIMD] = None
    fragmentation: Optional[Fragmentation] = None
    system: Optional[System] = None


# ----------------------------
# XYZ reader
# ----------------------------

def read_xyz_file(xyz_path: Union[str, Path]) -> Tuple[List[str], Any]:
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

    geom = reshape_geometry(coords, natom)
    return symbols, geom


# ----------------------------
# Geometry reshape
# ----------------------------

def reshape_geometry(flat: Any, natom: int) -> Any:
    req_type(flat, list, "molecule.geometry")
    if len(flat) != 3 * natom:
        die(f"molecule.geometry: expected length {3*natom}, got {len(flat)}")

    if np is not None:
        return np.asarray(flat, dtype=float).reshape((natom, 3))

    out = []
    for i in range(natom):
        x = float(flat[3*i + 0]); y = float(flat[3*i + 1]); z = float(flat[3*i + 2])
        out.append([x, y, z])
    return out


# ----------------------------
# Strict JSON -> Input
# ----------------------------

def parse_schema(obj: Any) -> SchemaTag:
    # require schema present for strictness
    if isinstance(obj, str):
        # allow "mqc-frag-1.0" or "mqc-frag"
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

def parse_model(d: Dict[str, Any]) -> Model:
    require_only_keys(d, {"method", "basis", "aux_basis", "solvent", "solvation_model",
                          "dielectric", "cpcm_nang", "cpcm_rscale"}, "model")
    method = req_type(d.get("method"), str, "model.method")
    basis = opt_type(d.get("basis"), str, "model.basis")
    aux = opt_type(d.get("aux_basis"), str, "model.aux_basis")

    # Parse solvation settings
    solvent = opt_type(d.get("solvent"), str, "model.solvent")
    solvation_model = opt_type(d.get("solvation_model"), str, "model.solvation_model")

    # Parse CPCM-specific settings
    dielectric = d.get("dielectric")
    if dielectric is not None:
        if not isinstance(dielectric, (int, float)):
            die("model.dielectric must be a number")
        dielectric = float(dielectric)
        if dielectric <= 0:
            die("model.dielectric must be > 0")

    cpcm_nang = d.get("cpcm_nang")
    if cpcm_nang is not None:
        if not isinstance(cpcm_nang, int):
            die("model.cpcm_nang must be an integer")
        if cpcm_nang <= 0:
            die("model.cpcm_nang must be > 0")

    cpcm_rscale = d.get("cpcm_rscale")
    if cpcm_rscale is not None:
        if not isinstance(cpcm_rscale, (int, float)):
            die("model.cpcm_rscale must be a number")
        cpcm_rscale = float(cpcm_rscale)
        if cpcm_rscale <= 0 or cpcm_rscale > 2.0:
            die("model.cpcm_rscale must be in range (0, 2.0]")

    # Validate solvent if specified
    if solvent is not None:
        solvent_lower = solvent.lower()
        if solvent_lower not in TBLITE_SOLVENTS:
            # Provide helpful error with suggestions
            die(f"model.solvent: unknown solvent '{solvent}'. "
                f"Supported solvents include: water, methanol, ethanol, acetone, acetonitrile, "
                f"benzene, toluene, chloroform, dichloromethane, dmso, dmf, thf, "
                f"diethylether, hexane, cyclohexane, and many more. "
                f"See documentation for full list.")
        # Store normalized (lowercase) solvent name
        solvent = solvent_lower

    # Validate solvation model if specified
    if solvation_model is not None:
        solvation_model_lower = solvation_model.lower()
        if solvation_model_lower not in SOLVATION_MODELS:
            die(f"model.solvation_model: unknown model '{solvation_model}'. "
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
        die("model.solvation_model: cannot specify solvation model without solvent or dielectric")

    # CPCM-specific validation
    if solvation_model == "cpcm":
        # CPCM needs either solvent or dielectric
        if solvent is None and dielectric is None:
            die("model: CPCM solvation requires either 'solvent' or 'dielectric'")
    else:
        # ALPB/GBSA need solvent name
        if solvation_model is not None and solvent is None:
            die(f"model: {solvation_model.upper()} solvation requires a solvent name")

    return Model(method=method, basis=basis, aux_basis=aux,
                 solvent=solvent, solvation_model=solvation_model,
                 dielectric=dielectric, cpcm_nang=cpcm_nang, cpcm_rscale=cpcm_rscale)

def parse_keywords(d: Dict[str, Any]) -> Tuple[Optional[SCF], Optional[Hessian], Optional[AIMD], Optional[Fragmentation]]:
    require_only_keys(d, {"scf", "hessian", "aimd", "fragmentation"}, "keywords")

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
        require_only_keys(hd, {"finite_difference_displacement", "displacement"}, "keywords.hessian")
        disp = hd.get("finite_difference_displacement") or hd.get("displacement")
        if disp is None:
            disp = 0.001  # default
        if not isinstance(disp, (int, float)):
            die("keywords.hessian.finite_difference_displacement must be number")
        if float(disp) <= 0:
            die("keywords.hessian.finite_difference_displacement must be > 0")
        hessian = Hessian(finite_difference_displacement=float(disp))

    aimd_obj = None
    if "aimd" in d:
        ad = req_type(d["aimd"], dict, "keywords.aimd")
        require_only_keys(ad, {"dt", "timestep", "nsteps", "steps", "initial_temperature", "temperature", "output_frequency", "output_freq"}, "keywords.aimd")
        dt = ad.get("dt") or ad.get("timestep")
        if dt is None:
            dt = 1.0  # default
        nsteps = ad.get("nsteps") or ad.get("steps")
        if nsteps is None:
            nsteps = 0  # default (no AIMD)
        init_temp = ad.get("initial_temperature") or ad.get("temperature")
        if init_temp is None:
            init_temp = 300.0  # default
        out_freq = ad.get("output_frequency") or ad.get("output_freq")
        if out_freq is None:
            out_freq = 1  # default

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

        aimd_obj = AIMD(dt=float(dt), nsteps=nsteps, initial_temperature=float(init_temp), output_frequency=out_freq)

    frag = None
    if "fragmentation" in d:
        fd = req_type(d["fragmentation"], dict, "keywords.fragmentation")
        require_only_keys(
            fd,
            {
                "method",
                "allow_overlapping_fragments",
                "level",
                "embedding",
                "cutoff_method",
                "distance_metric",
                "cutoffs",
            },
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

            # Parse cutoffs: support both named (dimer, trimer) and numeric (2, 3, 4, ...) keys
            cutoff_dict: Dict[int, float] = {}

            # Mapping from names to n-mer levels
            name_to_level = {
                "dimer": 2, "trimer": 3, "tetramer": 4, "pentamer": 5,
                "hexamer": 6, "heptamer": 7, "octamer": 8
            }

            for key, value in cd.items():
                # Determine the n-mer level
                if key in name_to_level:
                    nmer_level = name_to_level[key]
                elif isinstance(key, int):
                    nmer_level = key
                elif key.isdigit():
                    nmer_level = int(key)
                else:
                    die(f"keywords.fragmentation.cutoffs: unknown key '{key}'. " +
                        f"Use 'dimer', 'trimer', etc. or numeric keys like '2', '3', ...")

                # Validate nmer_level
                if nmer_level < 2:
                    die(f"keywords.fragmentation.cutoffs: n-mer level must be >= 2, got {nmer_level}")

                # Validate cutoff value
                if not isinstance(value, (int, float)) or float(value) <= 0:
                    die(f"keywords.fragmentation.cutoffs.{key} must be a positive number")

                cutoff_dict[nmer_level] = float(value)

            # Validate that cutoffs are monotonically decreasing (or equal)
            if cutoff_dict:
                sorted_levels = sorted(cutoff_dict.keys())
                for i in range(len(sorted_levels) - 1):
                    level_low = sorted_levels[i]
                    level_high = sorted_levels[i + 1]
                    cutoff_low = cutoff_dict[level_low]
                    cutoff_high = cutoff_dict[level_high]

                    if cutoff_high > cutoff_low:
                        die(f"keywords.fragmentation.cutoffs: {level_high}-mer cutoff ({cutoff_high}) " +
                            f"cannot be larger than {level_low}-mer cutoff ({cutoff_low}). " +
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

    return scf, hessian, aimd_obj, frag

def parse_system(d: Dict[str, Any]) -> System:
    require_only_keys(d, {"logger"}, "system")

    logger_data = req_type(d.get("logger"), dict, "system.logger")
    require_only_keys(logger_data, {"level"}, "system.logger")

    level = req_type(logger_data.get("level"), str, "system.logger.level")
    # Validate log level (case-insensitive check)
    valid_levels = {"debug", "verbose", "info", "performance", "warning", "error", "knowledge"}
    if level.lower() not in valid_levels:
        die(f"system.logger.level must be one of: {', '.join(sorted(valid_levels))}")

    return System(logger=Logger(level=level))

def parse_molecule(m: Dict[str, Any], base_dir: Path, mi: int) -> Molecule:
    require_only_keys(
        m,
        {
            "symbols",
            "geometry",
            "xyz",
            "molecular_charge",
            "molecular_multiplicity",
            "connectivity",
            "fragments",
            "fragment_charges",
            "fragment_multiplicities",
        },
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
        geom = reshape_geometry(m.get("geometry"), natom)
        source_xyz = None

    natom = len(symbols)

    # Molecular charge and multiplicity (required)
    charge = req_type(m.get("molecular_charge"), int, f"molecules[{mi}].molecular_charge")
    multiplicity = req_type(m.get("molecular_multiplicity"), int, f"molecules[{mi}].molecular_multiplicity")
    if multiplicity <= 0:
        die(f"molecules[{mi}].molecular_multiplicity must be >= 1")

    # fragments (optional - not needed for unfragmented calculations)
    frags_raw = m.get("fragments", None)

    if frags_raw is None:
        # Unfragmented calculation - no fragments specified
        fragments = []
        charges = []
        mults = []
    else:
        # Fragmented calculation
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

        # fragment_charges/fragment_multiplicities: optional with defaults, but if present must match nfrag
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

        # Validate that fragment charges sum to molecular charge
        total_frag_charge = sum(charges)
        if total_frag_charge != charge:
            die(f"molecules[{mi}]: fragment charges sum to {total_frag_charge}, "
                f"but molecular charge is {charge}. These must be equal.")

        # Validate multiplicity consistency
        # If all fragments are closed-shell (mult=1), molecule should be too
        if all(m == 1 for m in mults) and multiplicity != 1:
            die(f"molecules[{mi}]: all fragments have multiplicity=1 (closed-shell), "
                f"but molecular multiplicity is {multiplicity}. "
                f"This is inconsistent - molecular multiplicity should also be 1.")

        # If molecule is closed-shell but fragments aren't, that's suspicious
        if multiplicity == 1 and not all(m == 1 for m in mults):
            import sys
            print(f"WARNING: molecules[{mi}]: molecular multiplicity=1 but some fragments "
                  f"have multiplicity>1. Ensure fragment spins are properly paired.",
                  file=sys.stderr)

    # connectivity allowed (strictly typed), but not emitted to v1 yet
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
    json_path = Path(json_path)
    base_dir = json_path.parent

    data = json.loads(json_path.read_text(encoding="utf-8"))
    req_type(data, dict, "top-level")

    # strict top-level keys
    require_only_keys(data, {"schema", "molecules", "model", "keywords", "driver", "title", "system"}, "top-level")

    schema = parse_schema(data.get("schema"))

    model = parse_model(req_type(data.get("model"), dict, "model"))

    driver = req_type(data.get("driver"), str, "driver")

    # Optional title for output filename
    title = opt_type(data.get("title"), str, "title")

    kw = data.get("keywords", {})
    req_type(kw, dict, "keywords")
    scf, hessian, aimd, frag = parse_keywords(kw)

    # Parse system section (optional)
    system = None
    sys_data = data.get("system")
    if sys_data is not None:
        system = parse_system(req_type(sys_data, dict, "system"))

    mols = req_type(data.get("molecules"), list, "molecules")
    if len(mols) == 0:
        die("molecules must be non-empty list")
    molecules = []
    for mi, m in enumerate(mols):
        molecules.append(parse_molecule(req_type(m, dict, f"molecules[{mi}]"), base_dir, mi))

    return Input(schema=schema, model=model, driver=driver, molecules=molecules, title=title, scf=scf, hessian=hessian, aimd=aimd, fragmentation=frag, system=system)


# ----------------------------
# Emit v1 Fortran-input format
# ----------------------------

def fmt_float(x: float) -> str:
    # stable, compact; good for machine files
    return f"{x:.12g}"

def write_indices_block(f, indices: List[int], per_line: int = 20) -> None:
    f.write("%indices\n")
    for i in range(0, len(indices), per_line):
        chunk = indices[i:i+per_line]
        f.write(" ".join(str(v) for v in chunk) + "\n")
    f.write("end  ! indices\n")

def find_atom_fragments(atom_idx: int, fragments: List[List[int]]) -> List[int]:
    """Find all fragments (0-indexed) that contain the given atom index"""
    result = []
    for fi, frag in enumerate(fragments):
        if atom_idx in frag:
            result.append(fi)
    return result

def write_connectivity(f, mol: Molecule) -> None:
    """Write connectivity section, marking broken bonds"""
    if mol.connectivity is None or len(mol.connectivity) == 0:
        return

    f.write("%connectivity\n")
    f.write(f"nbonds = {len(mol.connectivity)}\n\n")

    # Count broken bonds for summary
    broken_count = 0

    for bi, (i, j, order) in enumerate(mol.connectivity):
        # For overlapping fragments: find ALL fragments containing each atom
        frags_i = find_atom_fragments(i, mol.fragments)
        frags_j = find_atom_fragments(j, mol.fragments)

        # Bond is broken if there exists a fragment containing one atom but not the other
        # This happens when the sets of fragments are different
        frags_i_set = set(frags_i)
        frags_j_set = set(frags_j)

        # Check if atoms always appear together
        # Bond is preserved only if both atoms appear in exactly the same fragments
        is_broken = (frags_i_set != frags_j_set)

        if is_broken:
            broken_count += 1
            # Warn if breaking a non-single bond
            if order != 1:
                import sys
                print(f"WARNING: Breaking non-single bond between atoms {i} and {j} (bond order {order})",
                      file=sys.stderr)
            # Format: atom_i atom_j order broken
            f.write(f"{i} {j} {order} broken\n")
        else:
            # Format: atom_i atom_j order preserved
            f.write(f"{i} {j} {order} preserved\n")

    f.write(f"\nnbroken = {broken_count}\n")
    f.write("end  ! connectivity\n\n")

def write_structure(f, mol: Molecule) -> None:
    """Write %structure section with molecular properties"""
    f.write("%structure\n")
    f.write(f"charge = {mol.charge}\n")
    f.write(f"multiplicity = {mol.multiplicity}\n")
    f.write("end  ! structure\n\n")

def write_geometry_xyz(f, mol: Molecule) -> None:
    """Write %geometry section in XYZ format"""
    symbols = mol.symbols
    geom = mol.geom_xyz
    natom = len(symbols)

    f.write("%geometry\n")
    f.write(f"{natom}\n")
    f.write("\n")  # blank line (comment line in XYZ format)

    # geometry Nx3
    if np is not None and hasattr(geom, "shape"):
        for s, (x, y, z) in zip(symbols, geom):
            f.write(f"{s} {fmt_float(float(x))} {fmt_float(float(y))} {fmt_float(float(z))}\n")
    else:
        for i, s in enumerate(symbols):
            x, y, z = geom[i]
            f.write(f"{s} {fmt_float(float(x))} {fmt_float(float(y))} {fmt_float(float(z))}\n")

    f.write("end  ! geometry\n\n")

def write_molecule_sections(f, mol: Molecule) -> None:
    """Write all molecule-specific sections (%structure, %geometry, %fragments, %connectivity)"""
    # %structure (always)
    write_structure(f, mol)

    # %geometry (always)
    write_geometry_xyz(f, mol)

    # %fragments (only if fragments are specified)
    if len(mol.fragments) > 0:
        f.write("%fragments\n")
        f.write(f"nfrag = {len(mol.fragments)}\n\n")

        for fi, frag_atoms in enumerate(mol.fragments):
            charge = mol.fragment_charges[fi]
            mult = mol.fragment_multiplicities[fi]

            f.write("%fragment\n")
            f.write(f"charge = {charge}\n")
            f.write(f"multiplicity = {mult}\n")
            write_indices_block(f, frag_atoms, per_line=24)
            f.write("end  ! fragment\n\n")

        f.write("end  ! fragments\n\n")

    # Write connectivity information (identifies broken bonds)
    # This must be outside %fragments section
    write_connectivity(f, mol)

def emit_v1(inp: Input, json_path: Path) -> Tuple[str, Path]:
    """
    Generates v1 Fortran-input format and writes to file.

    Returns tuple of (text_content, output_path).
    Supports single or multiple molecules.

    Output filename is determined by:
    - If inp.title is set, use "{title}.mqc"
    - Otherwise use "{json_stem}.mqc" (e.g., input.json -> input.mqc)
    """
    if len(inp.molecules) == 0:
        die(f"v1 emitter requires at least 1 molecule; got {len(inp.molecules)}")

    # Determine output filename
    if inp.title:
        # Sanitize title for filesystem (remove problematic chars)
        safe_title = "".join(c if c.isalnum() or c in "._- " else "_" for c in inp.title)
        out_path = json_path.parent / f"{safe_title}.mqc"
    else:
        out_path = json_path.parent / f"{json_path.stem}.mqc"

    lines: List[str] = []
    from io import StringIO
    buf = StringIO()

    # %schema (always)
    buf.write("%schema\n")
    buf.write(f"name = {inp.schema.name}\n")
    buf.write(f"version = {inp.schema.version}\n")
    buf.write(f"index_base = {inp.schema.index_base}\n")
    buf.write(f"units = {inp.schema.units}\n")
    buf.write("end  ! schema\n\n")

    # %model (always)
    buf.write("%model\n")
    buf.write(f"method = {inp.model.method}\n")
    if inp.model.basis is not None:
        buf.write(f"basis = {inp.model.basis}\n")
    if inp.model.aux_basis is not None:
        buf.write(f"aux_basis = {inp.model.aux_basis}\n")
    if inp.model.solvent is not None:
        buf.write(f"solvent = {inp.model.solvent}\n")
    if inp.model.solvation_model is not None:
        buf.write(f"solvation_model = {inp.model.solvation_model}\n")
    # CPCM-specific settings
    if inp.model.dielectric is not None:
        buf.write(f"dielectric = {fmt_float(inp.model.dielectric)}\n")
    if inp.model.cpcm_nang is not None:
        buf.write(f"cpcm_nang = {inp.model.cpcm_nang}\n")
    if inp.model.cpcm_rscale is not None:
        buf.write(f"cpcm_rscale = {fmt_float(inp.model.cpcm_rscale)}\n")
    buf.write("end  ! model\n\n")

    # %driver (always)
    buf.write("%driver\n")
    buf.write(f"type = {inp.driver}\n")
    buf.write("end  ! driver\n\n")

    # %system (optional)
    if inp.system is not None:
        buf.write("%system\n")
        buf.write(f"log_level = {inp.system.logger.level}\n")
        buf.write("end  ! system\n\n")

    # Handle single molecule (backward compatible) vs multiple molecules
    if len(inp.molecules) == 1:
        # Single molecule - write at top level (backward compatible format)
        mol = inp.molecules[0]
        write_molecule_sections(buf, mol)
    else:
        # Multiple molecules - use %molecules section
        buf.write("%molecules\n")
        buf.write(f"nmol = {len(inp.molecules)}\n\n")

        for mi, mol in enumerate(inp.molecules):
            buf.write("%molecule\n")
            if hasattr(mol, 'name') and mol.name:
                buf.write(f"name = {mol.name}\n")
            write_molecule_sections(buf, mol)
            buf.write("end  ! molecule\n\n")

        buf.write("end  ! molecules\n\n")

    # %scf (optional)
    if inp.scf is not None:
        buf.write("%scf\n")
        buf.write(f"maxiter = {inp.scf.maxiter}\n")
        buf.write(f"tolerance = {fmt_float(inp.scf.tolerance)}\n")
        buf.write("end  ! scf\n\n")

    # %hessian (optional)
    if inp.hessian is not None:
        buf.write("%hessian\n")
        buf.write(f"finite_difference_displacement = {fmt_float(inp.hessian.finite_difference_displacement)}\n")
        buf.write("end  ! hessian\n\n")

    # %aimd (optional)
    if inp.aimd is not None:
        buf.write("%aimd\n")
        buf.write(f"dt = {fmt_float(inp.aimd.dt)}\n")
        buf.write(f"nsteps = {inp.aimd.nsteps}\n")
        buf.write(f"initial_temperature = {fmt_float(inp.aimd.initial_temperature)}\n")
        buf.write(f"output_frequency = {inp.aimd.output_frequency}\n")
        buf.write("end  ! aimd\n\n")

    # %fragmentation (optional)
    if inp.fragmentation is not None:
        fk = inp.fragmentation
        buf.write("%fragmentation\n")
        buf.write(f"method = {fk.method}\n")
        buf.write(f"allow_overlapping_fragments = {'true' if fk.allow_overlapping_fragments else 'false'}\n")
        buf.write(f"level = {fk.level}\n")
        buf.write(f"embedding = {fk.embedding}\n")
        buf.write(f"cutoff_method = {fk.cutoff_method}\n")
        buf.write(f"distance_metric = {fk.distance_metric}\n")

        if fk.cutoffs is not None:
            buf.write("\n%cutoffs\n")
            # Write cutoffs sorted by n-mer level for consistency
            for nmer_level in sorted(fk.cutoffs.cutoffs.keys()):
                cutoff_value = fk.cutoffs.cutoffs[nmer_level]
                # Use numeric notation: "2 = 5.0", "3 = 4.0", etc.
                buf.write(f"{nmer_level} = {fmt_float(cutoff_value)}\n")
            buf.write("end  ! cutoffs\n")

        buf.write("end  ! fragmentation\n\n")

    text = buf.getvalue()

    # Always write to file
    out_path.write_text(text, encoding="utf-8")

    return text, out_path


# ----------------------------
# CLI
# ----------------------------

def main() -> None:
    import argparse
    p = argparse.ArgumentParser(
        description="Strict MQC JSON -> v1 Fortran-input emitter",
        epilog="Output file is automatically named based on input JSON or 'title' field."
    )
    p.add_argument("json_file", help="Input JSON file")
    args = p.parse_args()

    json_path = Path(args.json_file)
    if not json_path.exists():
        die(f"Input file not found: {json_path}")

    inp = parse_json_to_input(json_path)
    text, out_path = emit_v1(inp, json_path)

    print(f"Generated: {out_path}")
    print(f"  Lines: {len(text.splitlines())}")
    print(f"  Size: {len(text)} bytes")

if __name__ == "__main__":
    main()
