"""Command-line interface for the public-data OVA workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from .evidence import (
    EVIDENCE_COLUMNS,
    map_sites_to_sequence,
    read_literature_values,
    read_site_annotations,
)
from .protein import analyze_record, read_fasta
from .report import write_csv, write_evidence_svg, write_protein_summary, write_report


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Analyze OVA sequence properties and traceable literature evidence."
    )
    result.add_argument("--fasta", required=True)
    result.add_argument("--evidence", required=True)
    result.add_argument("--sites", required=True)
    result.add_argument("--out-dir", default="outputs")
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    records = read_fasta(args.fasta)
    proteins = [analyze_record(record) for record in records]
    evidence = read_literature_values(args.evidence)
    if len(records) != 1:
        raise ValueError("Site mapping requires exactly one reference FASTA record")
    sites = map_sites_to_sequence(read_site_annotations(args.sites), records[0].sequence)

    write_protein_summary(output / "protein_summary.csv", proteins)
    write_csv(output / "literature_evidence.csv", evidence, list(EVIDENCE_COLUMNS))
    write_csv(
        output / "validated_site_annotations.csv",
        sites,
        [
            "study_id", "publication_year", "doi", "pmid", "annotation_type",
            "residue", "reported_position", "uniprot_position", "numbering_offset",
            "reference_residue", "sequence_validation", "reference_context",
            "evidence_level", "notes",
        ],
    )
    write_evidence_svg(
        output / "hwang_2014_relative_response.svg",
        [row for row in evidence if row["study_id"] == "Hwang_2014"],
        "Hwang 2014: relative immune response",
        "Percent of intact OVA; abstract-reported mouse-study summaries",
    )
    write_evidence_svg(
        output / "wang_2013_glycation_extent.svg",
        [
            row for row in evidence
            if row["study_id"] == "Wang_2013" and row["metric"] == "modified_lysines"
        ],
        "Wang 2013: modified lysines",
        "Percent of total lysines; approximate abstract-reported values",
    )
    write_report(output / "report.md", proteins, evidence, sites)
    print(f"Analysis complete: {output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
