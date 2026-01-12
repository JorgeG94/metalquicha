"""Data models for MQC input representation."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple


# Supported solvents from tblite (src/tblite/solvation/data.f90)
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
class SchemaTag:
    """Schema metadata for the input file."""
    name: str
    version: str
    index_base: int = 0
    units: str = "angstrom"


@dataclass
class Model:
    """Computational model specification."""
    method: str
    basis: Optional[str] = None
    aux_basis: Optional[str] = None


@dataclass
class XTB:
    """XTB-specific settings (solvation, etc.)."""
    solvent: Optional[str] = None
    solvation_model: Optional[str] = None  # "alpb", "gbsa", or "cpcm"
    # CPCM-specific settings
    dielectric: Optional[float] = None
    cpcm_nang: Optional[int] = None
    cpcm_rscale: Optional[float] = None


@dataclass
class SCF:
    """SCF convergence settings."""
    maxiter: int
    tolerance: float


@dataclass
class Hessian:
    """Hessian calculation settings."""
    finite_difference_displacement: float = 0.001  # Bohr
    temperature: float = 298.15  # K for thermochemistry
    pressure: float = 1.0  # atm for thermochemistry


@dataclass
class AIMD:
    """Ab initio molecular dynamics settings."""
    dt: float = 1.0  # femtoseconds
    nsteps: int = 0
    initial_temperature: float = 300.0  # Kelvin
    output_frequency: int = 1


@dataclass
class FragCutoffs:
    """Distance cutoffs for fragment generation by n-mer level."""
    cutoffs: Dict[int, float]  # {2: 5.0, 3: 4.0, ...} for dimer, trimer, etc.


@dataclass
class Fragmentation:
    """Fragmentation method settings."""
    method: str
    allow_overlapping_fragments: bool
    level: int
    embedding: str
    cutoff_method: str
    distance_metric: str
    cutoffs: Optional[FragCutoffs]


@dataclass
class Logger:
    """Logging configuration."""
    level: str = "info"


@dataclass
class System:
    """System-level settings."""
    logger: Logger


@dataclass
class Molecule:
    """Molecular structure and fragment information."""
    symbols: List[str]
    geom_xyz: Any  # Nx3 list-of-lists or numpy array
    charge: int
    multiplicity: int
    fragments: List[List[int]]
    fragment_charges: List[int]
    fragment_multiplicities: List[int]
    connectivity: Optional[List[Tuple[int, int, int]]] = None
    source_xyz: Optional[str] = None


@dataclass
class Input:
    """Complete MQC input specification."""
    schema: SchemaTag
    model: Model
    driver: str
    molecules: List[Molecule]
    title: Optional[str] = None
    scf: Optional[SCF] = None
    hessian: Optional[Hessian] = None
    aimd: Optional[AIMD] = None
    fragmentation: Optional[Fragmentation] = None
    xtb: Optional[XTB] = None
    system: Optional[System] = None
