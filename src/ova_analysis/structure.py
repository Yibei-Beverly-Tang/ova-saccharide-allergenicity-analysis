"""Focused, dependency-free parsing and analysis of the RCSB PDB 1OVA structure."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re
import shlex


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "ACE": "M",  # 1OVA represents the N-acetylated initiator as ACE.
    "SEP": "S",
}


@dataclass(frozen=True)
class Atom:
    record: str
    atom_name: str
    residue_name: str
    chain: str
    residue_number: int
    insertion_code: str
    x: float
    y: float
    z: float
    b_factor: float
    element: str


def read_atoms(path: str | Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        alternate = line[16].strip()
        if alternate not in {"", "A"}:
            continue
        atoms.append(Atom(
            record=line[0:6].strip(),
            atom_name=line[12:16].strip(),
            residue_name=line[17:20].strip(),
            chain=line[21].strip(),
            residue_number=int(line[22:26]),
            insertion_code=line[26].strip(),
            x=float(line[30:38]),
            y=float(line[38:46]),
            z=float(line[46:54]),
            b_factor=float(line[60:66]),
            element=line[76:78].strip(),
        ))
    if not atoms:
        raise ValueError(f"No coordinates found in {path}")
    return atoms


def read_structure_metadata(path: str | Path) -> dict[str, object]:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    title = " ".join(line[10:].strip() for line in lines if line.startswith("TITLE"))
    methods = " ".join(line[10:].strip() for line in lines if line.startswith("EXPDTA"))
    resolution = None
    for line in lines:
        if line.startswith("REMARK   2 RESOLUTION.") and "ANGSTROMS" in line:
            match = re.search(r"RESOLUTION\.\s+([0-9.]+)\s+ANGSTROMS", line)
            if match:
                resolution = float(match.group(1))
    return {
        "pdb_id": "1OVA",
        "title": title,
        "experimental_method": methods,
        "resolution_angstrom": resolution,
        "pdb_doi": "10.2210/pdb1OVA/pdb",
    }


def read_poly_seq_scheme(path: str | Path) -> list[dict[str, object]]:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    start = next(
        (i for i, line in enumerate(lines) if line == "_pdbx_poly_seq_scheme.asym_id "),
        None,
    )
    if start is None:
        raise ValueError("mmCIF poly-sequence scheme not found")
    fields: list[str] = []
    index = start
    while index < len(lines) and lines[index].startswith("_pdbx_poly_seq_scheme."):
        fields.append(lines[index].split(".", 1)[1].strip())
        index += 1
    rows: list[dict[str, object]] = []
    while index < len(lines) and lines[index].strip() != "#":
        tokens = shlex.split(lines[index])
        if tokens:
            if len(tokens) != len(fields):
                raise ValueError("Unexpected mmCIF poly-sequence row")
            row = dict(zip(fields, tokens))
            row["seq_id"] = int(str(row["seq_id"]))
            row["auth_seq_num"] = (
                None if row["auth_seq_num"] in {"?", "."}
                else int(str(row["auth_seq_num"]))
            )
            rows.append(row)
        index += 1
    if not rows:
        raise ValueError("Empty mmCIF poly-sequence scheme")
    return rows


def validate_chain_sequence(
    scheme: list[dict[str, object]], sequence: str, chain: str = "A"
) -> None:
    rows = sorted(
        (row for row in scheme if row["pdb_strand_id"] == chain),
        key=lambda row: int(row["seq_id"]),
    )
    mapped = "".join(THREE_TO_ONE[str(row["mon_id"])] for row in rows)
    if mapped != sequence:
        mismatch = next(
            (
                index + 1 for index, (observed, expected) in
                enumerate(zip(mapped, sequence)) if observed != expected
            ),
            None,
        )
        raise ValueError(
            f"Structure/reference sequence mismatch at position {mismatch}; "
            f"lengths {len(mapped)} and {len(sequence)}"
        )


def _distance(atom_a: Atom, atom_b: Atom) -> float:
    return math.sqrt(
        (atom_a.x - atom_b.x) ** 2
        + (atom_a.y - atom_b.y) ** 2
        + (atom_a.z - atom_b.z) ** 2
    )


def _residue_atoms(
    atoms: list[Atom], chain: str, number: int, insertion_code: str
) -> list[Atom]:
    return [
        atom for atom in atoms
        if atom.record == "ATOM"
        and atom.chain == chain
        and atom.residue_number == number
        and atom.insertion_code == insertion_code
    ]


def analyze_structure_sites(
    sites: list[dict[str, object]],
    sequence: str,
    pdb_path: str | Path,
    cif_path: str | Path,
    chain: str = "A",
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    atoms = read_atoms(pdb_path)
    scheme = read_poly_seq_scheme(cif_path)
    validate_chain_sequence(scheme, sequence, chain)
    scheme_by_position = {
        int(row["seq_id"]): row
        for row in scheme
        if row["pdb_strand_id"] == chain
    }
    nag_atoms = [
        atom for atom in atoms
        if atom.record == "HETATM"
        and atom.chain == chain
        and atom.residue_name == "NAG"
    ]
    if not nag_atoms:
        raise ValueError(f"No NAG atoms found for chain {chain}")

    chain_c_alpha = [
        atom for atom in atoms
        if atom.record == "ATOM" and atom.chain == chain and atom.atom_name == "CA"
    ]
    centroid = (
        sum(atom.x for atom in chain_c_alpha) / len(chain_c_alpha),
        sum(atom.y for atom in chain_c_alpha) / len(chain_c_alpha),
        sum(atom.z for atom in chain_c_alpha) / len(chain_c_alpha),
    )

    mapped: list[dict[str, object]] = []
    ca_by_study_site: dict[str, Atom] = {}
    for site in sites:
        position = int(site["uniprot_position"])
        row = scheme_by_position.get(position)
        if row is None:
            raise ValueError(f"UniProt position {position} is absent from chain {chain}")
        pdb_number = int(row["auth_seq_num"])
        insertion = "" if row["pdb_ins_code"] in {".", "?"} else str(row["pdb_ins_code"])
        residue_atoms = _residue_atoms(atoms, chain, pdb_number, insertion)
        expected = THREE_TO_ONE[str(row["mon_id"])]
        if expected != sequence[position - 1]:
            raise ValueError(f"Residue mismatch at UniProt position {position}")
        ca = next((atom for atom in residue_atoms if atom.atom_name == "CA"), None)
        resolved = ca is not None
        if resolved:
            min_nag = min(_distance(atom, nag) for atom in residue_atoms for nag in nag_atoms)
            radial = math.sqrt(
                (ca.x - centroid[0]) ** 2
                + (ca.y - centroid[1]) ** 2
                + (ca.z - centroid[2]) ** 2
            )
            mean_b = sum(atom.b_factor for atom in residue_atoms) / len(residue_atoms)
            ca_by_study_site[f"{site['study_id']}:{site['reported_position']}"] = ca
        else:
            min_nag = radial = mean_b = None
        mapped.append({
            **site,
            "pdb_id": "1OVA",
            "pdb_chain": chain,
            "pdb_residue_name": row["mon_id"],
            "pdb_residue_number": pdb_number,
            "pdb_insertion_code": insertion,
            "pdb_residue_id": f"{pdb_number}{insertion}",
            "coordinate_resolved": resolved,
            "atom_count": len(residue_atoms),
            "distance_to_nag_angstrom": round(min_nag, 3) if min_nag is not None else "",
            "ca_distance_from_chain_centroid_angstrom": (
                round(radial, 3) if radial is not None else ""
            ),
            "mean_b_factor": round(mean_b, 3) if mean_b is not None else "",
        })

    pairwise: list[dict[str, object]] = []
    keys = list(ca_by_study_site)
    for index, key_a in enumerate(keys):
        for key_b in keys[index + 1:]:
            pairwise.append({
                "site_a": key_a,
                "site_b": key_b,
                "ca_distance_angstrom": round(
                    _distance(ca_by_study_site[key_a], ca_by_study_site[key_b]), 3
                ),
            })
    return mapped, pairwise


def analyze_epitope_structure_relationships(
    epitopes: list[dict[str, object]],
    structure_sites: list[dict[str, object]],
    sequence: str,
    pdb_path: str | Path,
    cif_path: str | Path,
    chain: str = "A",
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Map epitope intervals to 1OVA and calculate minimum site-to-epitope Cα distances."""
    atoms = read_atoms(pdb_path)
    scheme = read_poly_seq_scheme(cif_path)
    validate_chain_sequence(scheme, sequence, chain)
    scheme_by_position = {
        int(row["seq_id"]): row
        for row in scheme
        if row["pdb_strand_id"] == chain
    }
    ca_by_position: dict[int, Atom] = {}
    pdb_id_by_position: dict[int, str] = {}
    for position, row in scheme_by_position.items():
        if row["auth_seq_num"] is None:
            continue
        pdb_number = int(row["auth_seq_num"])
        insertion = "" if row["pdb_ins_code"] in {".", "?"} else str(row["pdb_ins_code"])
        residue_atoms = _residue_atoms(atoms, chain, pdb_number, insertion)
        ca = next((atom for atom in residue_atoms if atom.atom_name == "CA"), None)
        if ca is not None:
            ca_by_position[position] = ca
            pdb_id_by_position[position] = f"{pdb_number}{insertion}"

    epitope_mappings: list[dict[str, object]] = []
    resolved_positions_by_epitope: dict[int, list[int]] = {}
    for epitope in epitopes:
        identifier = int(epitope["iedb_epitope_id"])
        start = int(epitope["uniprot_start"])
        end = int(epitope["uniprot_end"])
        peptide = str(epitope["linear_sequence"])
        if sequence[start - 1:end] != peptide:
            raise ValueError(f"IEDB epitope {identifier} no longer matches P01012")
        positions = list(range(start, end + 1))
        resolved = [position for position in positions if position in ca_by_position]
        if not resolved:
            pdb_start = pdb_end = ""
        else:
            pdb_start = pdb_id_by_position[resolved[0]]
            pdb_end = pdb_id_by_position[resolved[-1]]
        resolved_positions_by_epitope[identifier] = resolved
        epitope_mappings.append({
            **epitope,
            "pdb_id": "1OVA",
            "pdb_chain": chain,
            "total_residue_count": len(positions),
            "resolved_residue_count": len(resolved),
            "coordinate_coverage_percent": round(100 * len(resolved) / len(positions), 3),
            "pdb_start_residue_id": pdb_start,
            "pdb_end_residue_id": pdb_end,
        })

    relationships: list[dict[str, object]] = []
    for site in structure_sites:
        site_position = int(site["uniprot_position"])
        site_ca = ca_by_position.get(site_position)
        for epitope in epitopes:
            identifier = int(epitope["iedb_epitope_id"])
            start = int(epitope["uniprot_start"])
            end = int(epitope["uniprot_end"])
            candidates = resolved_positions_by_epitope[identifier]
            within = start <= site_position <= end
            sequence_separation = 0 if within else min(
                abs(site_position - start), abs(site_position - end)
            )
            if site_ca is None or not candidates:
                nearest_position = None
                nearest_distance = None
            else:
                nearest_position = min(
                    candidates,
                    key=lambda position: _distance(site_ca, ca_by_position[position]),
                )
                nearest_distance = _distance(site_ca, ca_by_position[nearest_position])
            relationships.append({
                "site_key": f"{site['study_id']}:{site['reported_position']}",
                "study_id": site["study_id"],
                "annotation_type": site["annotation_type"],
                "reported_site": f"{site['residue']}-{site['reported_position']}",
                "site_uniprot_position": site_position,
                "site_pdb_residue_id": site["pdb_residue_id"],
                "iedb_epitope_id": identifier,
                "epitope_sequence": epitope["linear_sequence"],
                "epitope_uniprot_start": start,
                "epitope_uniprot_end": end,
                "site_within_epitope": within,
                "sequence_separation_residues": sequence_separation,
                "nearest_epitope_uniprot_position": (
                    nearest_position if nearest_position is not None else ""
                ),
                "nearest_epitope_residue": (
                    sequence[nearest_position - 1] if nearest_position is not None else ""
                ),
                "nearest_epitope_pdb_residue_id": (
                    pdb_id_by_position[nearest_position]
                    if nearest_position is not None else ""
                ),
                "min_ca_distance_angstrom": (
                    round(nearest_distance, 3) if nearest_distance is not None else ""
                ),
                "pdb_id": "1OVA",
                "pdb_chain": chain,
            })
    return epitope_mappings, relationships
