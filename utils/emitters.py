"""V1 Fortran-input format emitter for MQC."""

from __future__ import annotations
from io import StringIO
from pathlib import Path
from typing import List, Tuple

try:
    import numpy as np
except Exception:
    np = None

from .validation import die
from .models import Input, Molecule


def _fmt_float(x: float) -> str:
    """Format float in compact form."""
    return f"{x:.12g}"


def _write_indices_block(f, indices: List[int], per_line: int = 20) -> None:
    """Write an indices block."""
    f.write("%indices\n")
    for i in range(0, len(indices), per_line):
        chunk = indices[i:i+per_line]
        f.write(" ".join(str(v) for v in chunk) + "\n")
    f.write("end  ! indices\n")


def _find_atom_fragments(atom_idx: int, fragments: List[List[int]]) -> List[int]:
    """Find all fragments (0-indexed) that contain the given atom index."""
    result = []
    for fi, frag in enumerate(fragments):
        if atom_idx in frag:
            result.append(fi)
    return result


def _write_connectivity(f, mol: Molecule) -> None:
    """Write connectivity section, marking broken bonds."""
    if mol.connectivity is None or len(mol.connectivity) == 0:
        return

    f.write("%connectivity\n")
    f.write(f"nbonds = {len(mol.connectivity)}\n\n")

    broken_count = 0

    for bi, (i, j, order) in enumerate(mol.connectivity):
        frags_i = _find_atom_fragments(i, mol.fragments)
        frags_j = _find_atom_fragments(j, mol.fragments)

        frags_i_set = set(frags_i)
        frags_j_set = set(frags_j)

        is_broken = (frags_i_set != frags_j_set)

        if is_broken:
            broken_count += 1
            if order != 1:
                import sys
                print(f"WARNING: Breaking non-single bond between atoms {i} and {j} (bond order {order})",
                      file=sys.stderr)
            f.write(f"{i} {j} {order} broken\n")
        else:
            f.write(f"{i} {j} {order} preserved\n")

    f.write(f"\nnbroken = {broken_count}\n")
    f.write("end  ! connectivity\n\n")


def _write_structure(f, mol: Molecule) -> None:
    """Write %structure section with molecular properties."""
    f.write("%structure\n")
    f.write(f"charge = {mol.charge}\n")
    f.write(f"multiplicity = {mol.multiplicity}\n")
    f.write("end  ! structure\n\n")


def _write_geometry_xyz(f, mol: Molecule) -> None:
    """Write %geometry section in XYZ format."""
    symbols = mol.symbols
    geom = mol.geom_xyz
    natom = len(symbols)

    f.write("%geometry\n")
    f.write(f"{natom}\n")
    f.write("\n")

    if np is not None and hasattr(geom, "shape"):
        for s, (x, y, z) in zip(symbols, geom):
            f.write(f"{s} {_fmt_float(float(x))} {_fmt_float(float(y))} {_fmt_float(float(z))}\n")
    else:
        for i, s in enumerate(symbols):
            x, y, z = geom[i]
            f.write(f"{s} {_fmt_float(float(x))} {_fmt_float(float(y))} {_fmt_float(float(z))}\n")

    f.write("end  ! geometry\n\n")


def _write_molecule_sections(f, mol: Molecule) -> None:
    """Write all molecule-specific sections."""
    _write_structure(f, mol)
    _write_geometry_xyz(f, mol)

    if len(mol.fragments) > 0:
        f.write("%fragments\n")
        f.write(f"nfrag = {len(mol.fragments)}\n\n")

        for fi, frag_atoms in enumerate(mol.fragments):
            charge = mol.fragment_charges[fi]
            mult = mol.fragment_multiplicities[fi]

            f.write("%fragment\n")
            f.write(f"charge = {charge}\n")
            f.write(f"multiplicity = {mult}\n")
            _write_indices_block(f, frag_atoms, per_line=24)
            f.write("end  ! fragment\n\n")

        f.write("end  ! fragments\n\n")

    _write_connectivity(f, mol)


def emit_v1(inp: Input, json_path: Path) -> Tuple[str, Path]:
    """
    Generate v1 Fortran-input format and write to file.

    Returns tuple of (text_content, output_path).
    """
    if len(inp.molecules) == 0:
        die(f"v1 emitter requires at least 1 molecule; got {len(inp.molecules)}")

    if inp.title:
        safe_title = "".join(c if c.isalnum() or c in "._- " else "_" for c in inp.title)
        out_path = json_path.parent / f"{safe_title}.mqc"
    else:
        out_path = json_path.parent / f"{json_path.stem}.mqc"

    buf = StringIO()

    # %schema
    buf.write("%schema\n")
    buf.write(f"name = {inp.schema.name}\n")
    buf.write(f"version = {inp.schema.version}\n")
    buf.write(f"index_base = {inp.schema.index_base}\n")
    buf.write(f"units = {inp.schema.units}\n")
    buf.write("end  ! schema\n\n")

    # %model
    buf.write("%model\n")
    buf.write(f"method = {inp.model.method}\n")
    if inp.model.basis is not None:
        buf.write(f"basis = {inp.model.basis}\n")
    if inp.model.aux_basis is not None:
        buf.write(f"aux_basis = {inp.model.aux_basis}\n")
    buf.write("end  ! model\n\n")

    # %driver
    buf.write("%driver\n")
    buf.write(f"type = {inp.driver}\n")
    buf.write("end  ! driver\n\n")

    # %system (optional)
    if inp.system is not None:
        buf.write("%system\n")
        buf.write(f"log_level = {inp.system.logger.level}\n")
        buf.write("end  ! system\n\n")

    # Molecules
    if len(inp.molecules) == 1:
        mol = inp.molecules[0]
        _write_molecule_sections(buf, mol)
    else:
        buf.write("%molecules\n")
        buf.write(f"nmol = {len(inp.molecules)}\n\n")

        for mi, mol in enumerate(inp.molecules):
            buf.write("%molecule\n")
            if hasattr(mol, 'name') and mol.name:
                buf.write(f"name = {mol.name}\n")
            _write_molecule_sections(buf, mol)
            buf.write("end  ! molecule\n\n")

        buf.write("end  ! molecules\n\n")

    # %scf (optional)
    if inp.scf is not None:
        buf.write("%scf\n")
        buf.write(f"maxiter = {inp.scf.maxiter}\n")
        buf.write(f"tolerance = {_fmt_float(inp.scf.tolerance)}\n")
        buf.write("end  ! scf\n\n")

    # %xtb (optional)
    if inp.xtb is not None:
        buf.write("%xtb\n")
        if inp.xtb.solvent is not None:
            buf.write(f"solvent = {inp.xtb.solvent}\n")
        if inp.xtb.solvation_model is not None:
            buf.write(f"solvation_model = {inp.xtb.solvation_model}\n")
        if inp.xtb.dielectric is not None:
            buf.write(f"dielectric = {_fmt_float(inp.xtb.dielectric)}\n")
        if inp.xtb.cpcm_nang is not None:
            buf.write(f"cpcm_nang = {inp.xtb.cpcm_nang}\n")
        if inp.xtb.cpcm_rscale is not None:
            buf.write(f"cpcm_rscale = {_fmt_float(inp.xtb.cpcm_rscale)}\n")
        buf.write("end  ! xtb\n\n")

    # %hessian (optional)
    if inp.hessian is not None:
        buf.write("%hessian\n")
        buf.write(f"finite_difference_displacement = {_fmt_float(inp.hessian.finite_difference_displacement)}\n")
        buf.write(f"temperature = {_fmt_float(inp.hessian.temperature)}\n")
        buf.write(f"pressure = {_fmt_float(inp.hessian.pressure)}\n")
        buf.write("end  ! hessian\n\n")

    # %aimd (optional)
    if inp.aimd is not None:
        buf.write("%aimd\n")
        buf.write(f"dt = {_fmt_float(inp.aimd.dt)}\n")
        buf.write(f"nsteps = {inp.aimd.nsteps}\n")
        buf.write(f"initial_temperature = {_fmt_float(inp.aimd.initial_temperature)}\n")
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
            for nmer_level in sorted(fk.cutoffs.cutoffs.keys()):
                cutoff_value = fk.cutoffs.cutoffs[nmer_level]
                buf.write(f"{nmer_level} = {_fmt_float(cutoff_value)}\n")
            buf.write("end  ! cutoffs\n")

        buf.write("end  ! fragmentation\n\n")

    text = buf.getvalue()
    out_path.write_text(text, encoding="utf-8")

    return text, out_path
