"""Validation and sequence mapping for conservative IEDB epitope snapshots."""

from __future__ import annotations

import csv
from datetime import date
from pathlib import Path

from .protein import VALID_AA


EPITOPE_COLUMNS = (
    "iedb_epitope_id",
    "linear_sequence",
    "immune_response_type",
    "structure_type",
    "total_assay_count",
    "positive_assay_count",
    "positive_reference_count",
    "most_common_positive_host",
    "most_common_positive_mhc",
    "evidence_level",
    "source_url",
    "query_url",
    "retrieved_date",
    "notes",
)
SUPPORTED_EVIDENCE_LEVEL = "iedb_positive_assay_summary"


def read_iedb_epitopes(path: str | Path) -> list[dict[str, object]]:
    """Read a curated IEDB snapshot and reject incomplete provenance."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [
            column for column in EPITOPE_COLUMNS
            if column not in (reader.fieldnames or [])
        ]
        if missing:
            raise ValueError(f"Missing required epitope columns: {', '.join(missing)}")
        raw_rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]
    if not raw_rows:
        raise ValueError(f"No epitope rows found in {path}")

    parsed: list[dict[str, object]] = []
    seen_ids: set[int] = set()
    for line, row in enumerate(raw_rows, start=2):
        try:
            epitope_id = int(row["iedb_epitope_id"])
            total_count = int(row["total_assay_count"])
            assay_count = int(row["positive_assay_count"])
            reference_count = int(row["positive_reference_count"])
            date.fromisoformat(row["retrieved_date"])
        except ValueError as exc:
            raise ValueError(f"Invalid numeric or date value on line {line}") from exc
        sequence = row["linear_sequence"].upper()
        invalid = sorted(set(sequence) - VALID_AA)
        if not sequence or invalid:
            raise ValueError(f"Invalid linear peptide sequence on line {line}")
        if epitope_id < 1 or epitope_id in seen_ids:
            raise ValueError(f"Invalid or duplicate IEDB epitope ID on line {line}")
        if row["structure_type"] != "Linear peptide":
            raise ValueError(f"Unsupported epitope structure type on line {line}")
        if row["evidence_level"] != SUPPORTED_EVIDENCE_LEVEL:
            raise ValueError(f"Unsupported epitope evidence level on line {line}")
        if (
            total_count < 1
            or assay_count < 1
            or assay_count > total_count
            or reference_count < 1
            or reference_count > assay_count
        ):
            raise ValueError(f"Invalid positive evidence counts on line {line}")
        expected_url = f"https://www.iedb.org/epitope/{epitope_id}"
        if row["source_url"] != expected_url:
            raise ValueError(f"IEDB source URL does not match ID on line {line}")
        expected_query = (
            "https://query-api.iedb.org/api/v1/tcell_search?"
            f"structure_id=eq.{epitope_id}"
        )
        if row["query_url"] != expected_query:
            raise ValueError(f"IEDB query URL does not match ID on line {line}")
        seen_ids.add(epitope_id)
        parsed.append({
            **row,
            "iedb_epitope_id": epitope_id,
            "linear_sequence": sequence,
            "total_assay_count": total_count,
            "positive_assay_count": assay_count,
            "positive_reference_count": reference_count,
        })
    return parsed


def map_epitopes_to_sequence(
    rows: list[dict[str, object]], reference_sequence: str
) -> list[dict[str, object]]:
    """Map each linear epitope by an exact, unique P01012 sequence match."""
    mapped: list[dict[str, object]] = []
    for row in rows:
        peptide = str(row["linear_sequence"])
        starts: list[int] = []
        cursor = reference_sequence.find(peptide)
        while cursor >= 0:
            starts.append(cursor + 1)
            cursor = reference_sequence.find(peptide, cursor + 1)
        if not starts:
            raise ValueError(
                f"IEDB epitope {row['iedb_epitope_id']} does not match P01012"
            )
        if len(starts) > 1:
            raise ValueError(
                f"IEDB epitope {row['iedb_epitope_id']} maps ambiguously to P01012"
            )
        start = starts[0]
        mapped.append({
            **row,
            "uniprot_start": start,
            "uniprot_end": start + len(peptide) - 1,
            "sequence_validation": "exact_unique_match",
        })
    return mapped
