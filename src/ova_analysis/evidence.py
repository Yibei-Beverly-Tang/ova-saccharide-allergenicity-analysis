"""Validated loading of publication-reported summary evidence."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re


EVIDENCE_COLUMNS = (
    "study_id", "publication_year", "doi", "pmid", "treatment", "metric",
    "value", "unit", "comparator", "evidence_level", "source_location", "notes",
)
SITE_COLUMNS = (
    "study_id", "publication_year", "doi", "pmid", "annotation_type", "residue",
    "reported_position", "reference_context", "evidence_level", "notes",
)
DOI_PATTERN = re.compile(r"^10\.\d{4,9}/\S+$", re.IGNORECASE)
EVIDENCE_VALUE_LEVELS = {"abstract_reported_summary"}
SITE_EVIDENCE_LEVELS = {"abstract_reported_annotation"}
BOUNDED_PERCENT_METRICS = {"sequence_coverage", "modified_lysines"}


def _validate_source(row: dict[str, str], line: int) -> None:
    """Validate identifiers shared by numerical and site evidence."""
    if not row["study_id"] or not row["doi"]:
        raise ValueError(f"Missing study_id or doi on line {line}")
    if not DOI_PATTERN.fullmatch(row["doi"]):
        raise ValueError(f"Invalid DOI on line {line}: {row['doi']}")
    if row.get("pmid") and not row["pmid"].isdigit():
        raise ValueError(f"Invalid PMID on line {line}: {row['pmid']}")


def _parse_year(value: str, line: int) -> int:
    try:
        year = int(value)
    except ValueError as exc:
        raise ValueError(f"Invalid publication year on line {line}") from exc
    if not 1900 <= year <= 2100:
        raise ValueError(f"Implausible publication year on line {line}")
    return year


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
        _validate_source(row, line)
        if not row["metric"] or not row["unit"] or not row["source_location"]:
            raise ValueError(f"Missing metric, unit, or source location on line {line}")
        if row["evidence_level"] not in EVIDENCE_VALUE_LEVELS:
            raise ValueError(f"Unsupported evidence level on line {line}")
        try:
            value = float(row["value"])
        except ValueError as exc:
            raise ValueError(f"Invalid value on line {line}") from exc
        if not math.isfinite(value):
            raise ValueError(f"Non-finite value on line {line}")
        if "percent" in row["unit"].lower() and value < 0:
            raise ValueError(f"Negative percentage on line {line}")
        if row["metric"] in BOUNDED_PERCENT_METRICS and value > 100:
            raise ValueError(f"Bounded percentage above 100 on line {line}")
        year = _parse_year(row["publication_year"], line)
        parsed.append({**row, "publication_year": year, "value": value})
    return parsed


def read_site_annotations(path: str | Path) -> list[dict[str, object]]:
    parsed: list[dict[str, object]] = []
    for line, row in enumerate(_read(path, SITE_COLUMNS), start=2):
        _validate_source(row, line)
        if not row["residue"] or not row["annotation_type"]:
            raise ValueError(f"Missing residue or annotation type on line {line}")
        if row["evidence_level"] not in SITE_EVIDENCE_LEVELS:
            raise ValueError(f"Unsupported site evidence level on line {line}")
        try:
            position = int(row["reported_position"])
        except ValueError as exc:
            raise ValueError(f"Invalid residue position on line {line}") from exc
        if position < 1:
            raise ValueError(f"Invalid residue position on line {line}")
        year = _parse_year(row["publication_year"], line)
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
