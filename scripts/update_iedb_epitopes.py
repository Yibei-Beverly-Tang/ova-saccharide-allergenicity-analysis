"""Regenerate the conservative OVA epitope snapshot from the official IEDB API."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import date
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


API = "https://query-api.iedb.org/api/v1/tcell_search"
EPITOPES = {
    58560: "SIINFEKL",
    28676: "ISQAVHAAHAEINEAGR",
}
SELECT = ",".join([
    "structure_id",
    "linear_sequence",
    "structure_type",
    "qualitative_measure",
    "reference_id",
    "host_organism_name",
    "mhc_allele_name",
    "parent_source_antigen_name",
    "parent_source_antigen_source_org_name",
])
FIELDS = [
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
]


def query_url(epitope_id: int, include_projection: bool = True) -> str:
    parameters = {"structure_id": f"eq.{epitope_id}"}
    if include_projection:
        parameters.update({"select": SELECT, "limit": "1000"})
    return f"{API}?{urlencode(parameters)}"


def fetch_records(epitope_id: int) -> list[dict[str, object]]:
    with urlopen(query_url(epitope_id), timeout=120) as response:
        rows = json.load(response)
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"IEDB returned no records for epitope {epitope_id}")
    return rows


def most_common(rows: list[dict[str, object]], field: str) -> str:
    counts = Counter(str(row[field]) for row in rows if row.get(field))
    if not counts:
        raise ValueError(f"IEDB records have no values for {field}")
    return counts.most_common(1)[0][0]


def summarize(epitope_id: int, expected_sequence: str) -> dict[str, object]:
    rows = fetch_records(epitope_id)
    sequences = {str(row.get("linear_sequence")) for row in rows}
    if sequences != {expected_sequence}:
        raise ValueError(
            f"Unexpected sequences for IEDB epitope {epitope_id}: {sequences}"
        )
    if any(row.get("structure_type") != "Linear peptide" for row in rows):
        raise ValueError(f"Non-linear record found for IEDB epitope {epitope_id}")
    if any(
        row.get("parent_source_antigen_name") != "Ovalbumin (UniProt:P01012)"
        for row in rows
    ):
        raise ValueError(f"Non-P01012 record found for IEDB epitope {epitope_id}")
    positive = [
        row for row in rows
        if str(row.get("qualitative_measure", "")).startswith("Positive")
    ]
    references = {
        int(row["reference_id"])
        for row in positive
        if row.get("reference_id") is not None
    }
    return {
        "iedb_epitope_id": epitope_id,
        "linear_sequence": expected_sequence,
        "immune_response_type": "T cell",
        "structure_type": "Linear peptide",
        "total_assay_count": len(rows),
        "positive_assay_count": len(positive),
        "positive_reference_count": len(references),
        "most_common_positive_host": most_common(positive, "host_organism_name"),
        "most_common_positive_mhc": most_common(positive, "mhc_allele_name"),
        "evidence_level": "iedb_positive_assay_summary",
        "source_url": f"https://www.iedb.org/epitope/{epitope_id}",
        "query_url": query_url(epitope_id, include_projection=False),
        "retrieved_date": date.today().isoformat(),
        "notes": (
            "Generated from the official IEDB query API by "
            "scripts/update_iedb_epitopes.py"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = [summarize(identifier, sequence) for identifier, sequence in EPITOPES.items()]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} verified IEDB epitope summaries to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
