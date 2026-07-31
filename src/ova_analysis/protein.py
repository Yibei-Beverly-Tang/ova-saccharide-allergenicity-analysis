"""Dependency-free protein sequence parsing and descriptive analysis."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
RESIDUE_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
    "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
    "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
    "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
PKA = {
    "n": 9.69, "c": 2.34, "C": 8.33, "D": 3.86, "E": 4.25,
    "H": 6.00, "K": 10.53, "R": 12.48, "Y": 10.07,
}


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


def read_fasta(path: str | Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    parts: list[str] = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(_record(header, parts))
            header, parts = line[1:].strip(), []
        else:
            if header is None:
                raise ValueError("FASTA sequence encountered before a header")
            parts.append(line)
    if header is not None:
        records.append(_record(header, parts))
    if not records:
        raise ValueError("No FASTA records found")
    return records


def _record(header: str, parts: list[str]) -> FastaRecord:
    sequence = "".join(parts).replace(" ", "").upper()
    invalid = sorted(set(sequence) - VALID_AA)
    if invalid:
        raise ValueError(f"Unsupported amino-acid symbols: {', '.join(invalid)}")
    identifier, _, description = header.partition(" ")
    return FastaRecord(identifier, description, sequence)


def molecular_weight(sequence: str) -> float:
    return sum(RESIDUE_MASS[aa] for aa in sequence) + 18.0153


def net_charge(sequence: str, ph: float) -> float:
    counts = Counter(sequence)
    positive = 1 / (1 + 10 ** (ph - PKA["n"]))
    positive += counts["K"] / (1 + 10 ** (ph - PKA["K"]))
    positive += counts["R"] / (1 + 10 ** (ph - PKA["R"]))
    positive += counts["H"] / (1 + 10 ** (ph - PKA["H"]))
    negative = 1 / (1 + 10 ** (PKA["c"] - ph))
    for aa in ("D", "E", "C", "Y"):
        negative += counts[aa] / (1 + 10 ** (PKA[aa] - ph))
    return positive - negative


def isoelectric_point(sequence: str) -> float:
    low, high = 0.0, 14.0
    for _ in range(80):
        mid = (low + high) / 2
        if net_charge(sequence, mid) > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def n_glycosylation_sequons(sequence: str) -> list[dict[str, int | str]]:
    return [
        {"position": index + 1, "motif": sequence[index:index + 3]}
        for index in range(len(sequence) - 2)
        if sequence[index] == "N"
        and sequence[index + 1] != "P"
        and sequence[index + 2] in {"S", "T"}
    ]


def analyze_record(record: FastaRecord) -> dict[str, object]:
    sequence = record.sequence
    counts = Counter(sequence)
    length = len(sequence)
    aromaticity = 100 * sum(counts[aa] for aa in "FWY") / length
    aliphatic_index = (
        100 * counts["A"] / length
        + 2.9 * 100 * counts["V"] / length
        + 3.9 * 100 * (counts["I"] + counts["L"]) / length
    )
    return {
        "identifier": record.identifier,
        "description": record.description,
        "length_aa": length,
        "molecular_weight_da": round(molecular_weight(sequence), 2),
        "theoretical_pi": round(isoelectric_point(sequence), 3),
        "gravy": round(sum(HYDROPATHY[aa] for aa in sequence) / length, 3),
        "aromaticity_percent": round(aromaticity, 3),
        "aliphatic_index": round(aliphatic_index, 3),
        "lysine_count": counts["K"],
        "arginine_count": counts["R"],
        "candidate_basic_residues": counts["K"] + counts["R"],
        "n_glycosylation_sequons": n_glycosylation_sequons(sequence),
    }

