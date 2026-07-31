"""Validated loading of publication-reported summary evidence."""

from __future__ import annotations

import csv
from pathlib import Path


EVIDENCE_COLUMNS = (
    "study_id", "publication_year", "doi", "pmid", "treatment", "metric",
    "value", "unit", "comparator", "evidence_level", "source_location", "notes",
)
SITE_COLUMNS = (
    "study_id", "publication_year", "doi", "pmid", "annotation_type", "residue",
    "reported_position", "reference_context", "evidence_level", "notes",
)


def _read(path: str | Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [name for name in required if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    return rows


def read_literature_values(path: str | Path) -> list[dict[str, object]]:
    parsed: list[dict[str, object]] = []
    for line, row in enumerate(_read(path, EVIDENCE_COLUMNS), start=2):
        if not row["study_id"] or not row["doi"] or not row["metric"]:
            raise ValueError(f"Missing study_id, doi, or metric on line {line}")
        if row["evidence_level"] != "abstract_reported_summary":
            raise ValueError(f"Unsupported evidence level on line {line}")
        try:
            year, value = int(row["publication_year"]), float(row["value"])
        except ValueError as exc:
            raise ValueError(f"Invalid year or value on line {line}") from exc
        if not 1900 <= year <= 2100:
            raise ValueError(f"Implausible year on line {line}")
        parsed.append({**row, "publication_year": year, "value": value})
    return parsed


def read_site_annotations(path: str | Path) -> list[dict[str, object]]:
    parsed: list[dict[str, object]] = []
    for line, row in enumerate(_read(path, SITE_COLUMNS), start=2):
        if not row["study_id"] or not row["doi"] or not row["residue"]:
            raise ValueError(f"Missing study_id, doi, or residue on line {line}")
        try:
            year = int(row["publication_year"])
            position = int(row["reported_position"])
        except ValueError as exc:
            raise ValueError(f"Invalid year or residue position on line {line}") from exc
        if position < 1:
            raise ValueError(f"Invalid residue position on line {line}")
        parsed.append({**row, "publication_year": year, "reported_position": position})
    return parsed


def map_sites_to_sequence(
    rows: list[dict[str, object]], sequence: str
) -> list[dict[str, object]]:
    """Map publication numbering to a reference sequence, allowing ±1 offsets."""
    residue_codes = {"Asn": "N", "Lys": "K", "Arg": "R"}
    mapped: list[dict[str, object]] = []
    for row in rows:
        expected = residue_codes.get(str(row["residue"]))
        if expected is None:
            raise ValueError(f"Unsupported residue name: {row['residue']}")
        reported = int(row["reported_position"])
        candidates = [reported, reported + 1, reported - 1]
        match = next(
            (
                position for position in candidates
                if 1 <= position <= len(sequence) and sequence[position - 1] == expected
            ),
            None,
        )
        if match is None:
            raise ValueError(
                f"{row['study_id']} {row['residue']}-{reported} does not match "
                "the reference sequence at offsets 0, +1, or -1"
            )
        mapped.append({
            **row,
            "uniprot_position": match,
            "numbering_offset": match - reported,
            "reference_residue": sequence[match - 1],
            "sequence_validation": "matched",
        })
    return mapped
