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
from .report import (
    write_csv,
    write_evidence_svg,
    write_protein_summary,
    write_pymol_script,
    write_report,
    write_structure_report,
    write_structure_svg,
)
from .structure import (
    analyze_structure_sites,
    read_structure_metadata,
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Analyze OVA sequence properties and traceable literature evidence."
    )
    result.add_argument("--fasta", required=True)
    result.add_argument("--evidence", required=True)
    result.add_argument("--sites", required=True)
    result.add_argument("--structure-pdb", default="structures/1OVA.pdb")
    result.add_argument("--structure-cif", default="structures/1OVA.cif")
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
    structure_sites, pairwise = analyze_structure_sites(
        sites,
        records[0].sequence,
        args.structure_pdb,
        args.structure_cif,
    )
    structure_metadata = read_structure_metadata(args.structure_pdb)

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
    structure_fields = [
        "study_id", "doi", "annotation_type", "residue", "reported_position",
        "uniprot_position", "pdb_id", "pdb_chain", "pdb_residue_name",
        "pdb_residue_number", "pdb_insertion_code", "pdb_residue_id",
        "coordinate_resolved", "atom_count", "distance_to_nag_angstrom",
        "ca_distance_from_chain_centroid_angstrom", "mean_b_factor",
    ]
    write_csv(output / "structure_site_mapping.csv", structure_sites, structure_fields)
    write_csv(
        output / "structure_site_pairwise_distances.csv",
        pairwise,
        ["site_a", "site_b", "ca_distance_angstrom"],
    )
    write_structure_svg(output / "1ova_site_distance_to_nag.svg", structure_sites)
    write_structure_report(
        output / "structure_report.md",
        structure_metadata,
        structure_sites,
        pairwise,
    )
    scripts = Path("scripts")
    scripts.mkdir(exist_ok=True)
    write_pymol_script(scripts / "visualize_1ova_sites.pml", structure_sites)
    print(f"Analysis complete: {output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
