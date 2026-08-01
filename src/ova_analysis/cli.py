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
from .epitopes import EPITOPE_COLUMNS, map_epitopes_to_sequence, read_iedb_epitopes
from .protein import analyze_record, read_fasta
from .quality import QUALITY_FIELDS, classify_evidence
from .report import (
    write_csv,
    write_epitope_report,
    write_epitope_svg,
    write_evidence_svg,
    write_protein_summary,
    write_pymol_script,
    write_report,
    write_structure_report,
    write_structure_svg,
    write_site_epitope_report,
    write_site_epitope_svg,
    write_quality_report,
    write_sasa_report,
    write_sasa_svg,
)
from .sasa import calculate_residue_sasa, summarize_sasa_targets
from .structure import (
    analyze_structure_sites,
    analyze_epitope_structure_relationships,
    read_structure_metadata,
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Analyze OVA sequence properties and traceable literature evidence."
    )
    result.add_argument("--fasta", required=True)
    result.add_argument("--evidence", required=True)
    result.add_argument("--sites", required=True)
    result.add_argument("--epitopes", required=True)
    result.add_argument("--structure-pdb", default="structures/1OVA.pdb")
    result.add_argument("--structure-cif", default="structures/1OVA.cif")
    result.add_argument("--out-dir", default="outputs")
    result.add_argument(
        "--pymol-script",
        default="scripts/visualize_1ova_sites.pml",
        help="Path for the generated PyMOL visualization script.",
    )
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
    epitopes = map_epitopes_to_sequence(
        read_iedb_epitopes(args.epitopes), records[0].sequence
    )
    structure_sites, pairwise = analyze_structure_sites(
        sites,
        records[0].sequence,
        args.structure_pdb,
        args.structure_cif,
    )
    structure_metadata = read_structure_metadata(args.structure_pdb)
    epitope_structure, site_epitope_relationships = (
        analyze_epitope_structure_relationships(
            epitopes,
            structure_sites,
            records[0].sequence,
            args.structure_pdb,
            args.structure_cif,
        )
    )
    residue_sasa = calculate_residue_sasa(
        records[0].sequence,
        args.structure_pdb,
        args.structure_cif,
    )
    site_sasa, epitope_sasa = summarize_sasa_targets(
        residue_sasa, structure_sites, epitopes
    )
    quality_rows = classify_evidence(evidence, sites, epitopes)

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
    epitope_fields = list(EPITOPE_COLUMNS) + [
        "uniprot_start", "uniprot_end", "sequence_validation",
    ]
    write_csv(output / "validated_iedb_epitopes.csv", epitopes, epitope_fields)
    write_epitope_svg(
        output / "iedb_epitope_map.svg", epitopes, len(records[0].sequence)
    )
    write_epitope_report(output / "epitope_report.md", epitopes)
    epitope_structure_fields = epitope_fields + [
        "pdb_id", "pdb_chain", "total_residue_count", "resolved_residue_count",
        "coordinate_coverage_percent", "pdb_start_residue_id", "pdb_end_residue_id",
    ]
    write_csv(
        output / "epitope_structure_mapping.csv",
        epitope_structure,
        epitope_structure_fields,
    )
    site_epitope_fields = [
        "site_key", "study_id", "annotation_type", "reported_site",
        "site_uniprot_position", "site_pdb_residue_id", "iedb_epitope_id",
        "epitope_sequence", "epitope_uniprot_start", "epitope_uniprot_end",
        "site_within_epitope", "sequence_separation_residues",
        "nearest_epitope_uniprot_position", "nearest_epitope_residue",
        "nearest_epitope_pdb_residue_id", "min_ca_distance_angstrom",
        "pdb_id", "pdb_chain",
    ]
    write_csv(
        output / "modification_site_epitope_distances.csv",
        site_epitope_relationships,
        site_epitope_fields,
    )
    write_site_epitope_svg(
        output / "site_epitope_distances.svg", site_epitope_relationships
    )
    write_site_epitope_report(
        output / "site_epitope_report.md",
        epitope_structure,
        site_epitope_relationships,
    )
    residue_sasa_fields = [
        "uniprot_position", "residue_one_letter", "pdb_id", "pdb_chain",
        "pdb_residue_name", "pdb_residue_id", "atom_count",
        "total_sasa_angstrom2", "backbone_sasa_angstrom2",
        "sidechain_sasa_angstrom2", "probe_radius_angstrom",
        "sphere_point_count", "occlusion_context",
    ]
    write_csv(output / "residue_sasa.csv", residue_sasa, residue_sasa_fields)
    site_sasa_fields = [
        "study_id", "doi", "annotation_type", "residue", "reported_position",
        "uniprot_position", "pdb_residue_name", "pdb_residue_id",
        "total_sasa_angstrom2", "backbone_sasa_angstrom2",
        "sidechain_sasa_angstrom2", "probe_radius_angstrom",
        "sphere_point_count", "occlusion_context",
    ]
    write_csv(output / "annotated_site_sasa.csv", site_sasa, site_sasa_fields)
    epitope_sasa_fields = epitope_fields + [
        "resolved_residue_count", "total_sasa_angstrom2",
        "mean_residue_sasa_angstrom2", "sidechain_sasa_angstrom2",
        "mean_sidechain_sasa_angstrom2", "probe_radius_angstrom",
        "sphere_point_count", "occlusion_context",
    ]
    write_csv(output / "epitope_sasa_summary.csv", epitope_sasa, epitope_sasa_fields)
    write_sasa_svg(output / "annotated_site_sasa.svg", site_sasa)
    write_sasa_report(
        output / "sasa_report.md", residue_sasa, site_sasa, epitope_sasa
    )
    write_csv(
        output / "evidence_quality_classification.csv",
        quality_rows,
        QUALITY_FIELDS,
    )
    write_quality_report(output / "evidence_quality_report.md", quality_rows)
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
    pymol_script = Path(args.pymol_script)
    pymol_script.parent.mkdir(parents=True, exist_ok=True)
    write_pymol_script(pymol_script, structure_sites, epitope_structure)
    print(f"Analysis complete: {output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
