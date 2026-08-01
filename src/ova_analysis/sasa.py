"""Deterministic Shrake–Rupley solvent-accessible surface-area calculations."""

from __future__ import annotations

from collections import defaultdict
import math
from pathlib import Path

from .structure import Atom, read_atoms, read_poly_seq_scheme, validate_chain_sequence


# Bondi-style van der Waals radii in Å for elements present in the analysis context.
VDW_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


def sphere_points(count: int) -> list[tuple[float, float, float]]:
    """Return deterministic, approximately uniform Fibonacci-sphere points."""
    if count < 1:
        raise ValueError("SASA sphere-point count must be positive")
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    points = []
    for index in range(count):
        z = 1.0 - 2.0 * (index + 0.5) / count
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        angle = index * golden_angle
        points.append((radius * math.cos(angle), radius * math.sin(angle), z))
    return points


def _element_radius(atom: Atom) -> float:
    element = atom.element.upper()
    if element not in VDW_RADII:
        raise ValueError(f"Unsupported element for SASA: {element or '(blank)'}")
    return VDW_RADII[element]


def atom_sasa(
    atoms: list[Atom], probe_radius: float = 1.4, point_count: int = 960
) -> list[float]:
    """Calculate per-atom SASA with a grid-accelerated Shrake–Rupley algorithm."""
    if probe_radius <= 0:
        raise ValueError("SASA probe radius must be positive")
    if not atoms:
        raise ValueError("Cannot calculate SASA for an empty atom list")
    points = sphere_points(point_count)
    radii = [_element_radius(atom) + probe_radius for atom in atoms]
    cell_size = 2.0 * max(radii)
    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)

    def cell(atom: Atom) -> tuple[int, int, int]:
        return (
            math.floor(atom.x / cell_size),
            math.floor(atom.y / cell_size),
            math.floor(atom.z / cell_size),
        )

    for index, atom in enumerate(atoms):
        grid[cell(atom)].append(index)

    areas: list[float] = []
    for index, atom in enumerate(atoms):
        cx, cy, cz = cell(atom)
        neighbors: list[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for candidate in grid.get((cx + dx, cy + dy, cz + dz), []):
                        if candidate == index:
                            continue
                        other = atoms[candidate]
                        cutoff = radii[index] + radii[candidate]
                        center_distance_sq = (
                            (atom.x - other.x) ** 2
                            + (atom.y - other.y) ** 2
                            + (atom.z - other.z) ** 2
                        )
                        if center_distance_sq < cutoff * cutoff:
                            neighbors.append(candidate)

        accessible = 0
        expanded = radii[index]
        for ux, uy, uz in points:
            x = atom.x + expanded * ux
            y = atom.y + expanded * uy
            z = atom.z + expanded * uz
            blocked = False
            for candidate in neighbors:
                other = atoms[candidate]
                if (
                    (x - other.x) ** 2
                    + (y - other.y) ** 2
                    + (z - other.z) ** 2
                    < radii[candidate] ** 2
                ):
                    blocked = True
                    break
            if not blocked:
                accessible += 1
        areas.append(
            4.0 * math.pi * expanded * expanded * accessible / point_count
        )
    return areas


def calculate_residue_sasa(
    sequence: str,
    pdb_path: str | Path,
    cif_path: str | Path,
    chain: str = "A",
    probe_radius: float = 1.4,
    point_count: int = 960,
) -> list[dict[str, object]]:
    """Calculate residue SASA for 1OVA chain A with same-chain NAG occlusion."""
    all_atoms = read_atoms(pdb_path)
    scheme = read_poly_seq_scheme(cif_path)
    validate_chain_sequence(scheme, sequence, chain)
    context_atoms = [
        atom for atom in all_atoms
        if atom.chain == chain
        and (atom.record == "ATOM" or atom.residue_name == "NAG")
    ]
    areas = atom_sasa(context_atoms, probe_radius, point_count)
    area_by_residue: dict[tuple[int, str], dict[str, float | int]] = defaultdict(
        lambda: {"total": 0.0, "backbone": 0.0, "sidechain": 0.0, "atoms": 0}
    )
    for atom, area in zip(context_atoms, areas):
        if atom.record != "ATOM":
            continue
        key = (atom.residue_number, atom.insertion_code)
        values = area_by_residue[key]
        values["total"] = float(values["total"]) + area
        category = "backbone" if atom.atom_name in BACKBONE_ATOMS else "sidechain"
        values[category] = float(values[category]) + area
        values["atoms"] = int(values["atoms"]) + 1

    rows: list[dict[str, object]] = []
    for scheme_row in sorted(
        (row for row in scheme if row["pdb_strand_id"] == chain),
        key=lambda row: int(row["seq_id"]),
    ):
        if scheme_row["auth_seq_num"] is None:
            continue
        position = int(scheme_row["seq_id"])
        insertion = (
            "" if scheme_row["pdb_ins_code"] in {".", "?"}
            else str(scheme_row["pdb_ins_code"])
        )
        pdb_number = int(scheme_row["auth_seq_num"])
        values = area_by_residue.get((pdb_number, insertion))
        if values is None:
            continue
        rows.append({
            "uniprot_position": position,
            "residue_one_letter": sequence[position - 1],
            "pdb_id": "1OVA",
            "pdb_chain": chain,
            "pdb_residue_name": scheme_row["mon_id"],
            "pdb_residue_id": f"{pdb_number}{insertion}",
            "atom_count": values["atoms"],
            "total_sasa_angstrom2": round(float(values["total"]), 3),
            "backbone_sasa_angstrom2": round(float(values["backbone"]), 3),
            "sidechain_sasa_angstrom2": round(float(values["sidechain"]), 3),
            "probe_radius_angstrom": probe_radius,
            "sphere_point_count": point_count,
            "occlusion_context": "1OVA chain A protein plus chain A NAG",
        })
    return rows


def summarize_sasa_targets(
    residue_rows: list[dict[str, object]],
    structure_sites: list[dict[str, object]],
    epitopes: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Join residue SASA to modification sites and aggregate epitope intervals."""
    by_position = {int(row["uniprot_position"]): row for row in residue_rows}
    site_rows = []
    for site in structure_sites:
        position = int(site["uniprot_position"])
        sasa = by_position.get(position)
        if sasa is None:
            raise ValueError(f"No SASA value for annotated UniProt position {position}")
        site_rows.append({**site, **sasa})

    epitope_rows = []
    for epitope in epitopes:
        positions = range(
            int(epitope["uniprot_start"]), int(epitope["uniprot_end"]) + 1
        )
        residues = [by_position[position] for position in positions if position in by_position]
        expected_count = int(epitope["uniprot_end"]) - int(epitope["uniprot_start"]) + 1
        if len(residues) != expected_count:
            raise ValueError(
                f"Incomplete SASA coverage for IEDB epitope {epitope['iedb_epitope_id']}"
            )
        total = sum(float(row["total_sasa_angstrom2"]) for row in residues)
        sidechain = sum(float(row["sidechain_sasa_angstrom2"]) for row in residues)
        epitope_rows.append({
            **epitope,
            "resolved_residue_count": len(residues),
            "total_sasa_angstrom2": round(total, 3),
            "mean_residue_sasa_angstrom2": round(total / len(residues), 3),
            "sidechain_sasa_angstrom2": round(sidechain, 3),
            "mean_sidechain_sasa_angstrom2": round(sidechain / len(residues), 3),
            "probe_radius_angstrom": residues[0]["probe_radius_angstrom"],
            "sphere_point_count": residues[0]["sphere_point_count"],
            "occlusion_context": residues[0]["occlusion_context"],
        })
    return site_rows, epitope_rows
