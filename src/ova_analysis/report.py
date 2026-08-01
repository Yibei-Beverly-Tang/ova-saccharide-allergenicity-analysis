"""Writers for machine-readable summaries, Markdown, and SVG charts."""

from __future__ import annotations

import csv
from collections import Counter
import html
from pathlib import Path


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_protein_summary(path: Path, proteins: list[dict[str, object]]) -> None:
    fields = [
        "identifier", "length_aa", "molecular_weight_da", "theoretical_pi",
        "gravy", "aromaticity_percent", "aliphatic_index", "lysine_count",
        "arginine_count", "candidate_basic_residues",
        "n_glycosylation_sequon_count", "n_glycosylation_sequons",
    ]
    rows = []
    for protein in proteins:
        motifs = protein["n_glycosylation_sequons"]
        rows.append({
            **protein,
            "n_glycosylation_sequon_count": len(motifs),
            "n_glycosylation_sequons": "; ".join(
                f"{hit['motif']}@{hit['position']}" for hit in motifs
            ),
        })
    write_csv(path, rows, fields)


def write_evidence_svg(
    path: Path,
    rows: list[dict[str, object]],
    title: str,
    subtitle: str,
) -> None:
    if not rows:
        raise ValueError("Cannot draw an empty evidence chart")
    width, height = 800, 130 + 58 * len(rows)
    left, chart_width = 270, 430
    scale = chart_width / max(100, max(float(row["value"]) for row in rows))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfcfe"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.title{font-size:20px;'
        'font-weight:700}.label{font-size:13px}.value{font-size:13px;fill:#43506a}</style>',
        f'<text x="28" y="38" class="title">{html.escape(title)}</text>',
        f'<text x="28" y="64" class="value">{html.escape(subtitle)}</text>',
    ]
    for index, row in enumerate(rows):
        y = 98 + 58 * index
        label = html.escape(
            str(row["metric"]).replace("relative_", "").replace("_", " ")
            + " — " + str(row["treatment"])
        )
        value = float(row["value"])
        bar_width = value * scale
        parts.extend([
            f'<text x="28" y="{y + 20}" class="label">{label}</text>',
            f'<rect x="{left}" y="{y}" width="{bar_width:.2f}" height="26" '
            'rx="5" fill="#457b9d"/>',
            f'<text x="{left + bar_width + 9:.2f}" y="{y + 19}" '
            f'class="value">{value:g}%</text>',
        ])
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_epitope_svg(
    path: Path,
    rows: list[dict[str, object]],
    protein_length: int,
) -> None:
    """Draw validated linear epitopes on the P01012 sequence axis."""
    if not rows or protein_length < 1:
        raise ValueError("Cannot draw an empty epitope map")
    width, height = 980, 150 + 62 * len(rows)
    left, axis_width = 95, 790
    scale = axis_width / protein_length
    axis_y = 92
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfcfe"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.title{font-size:20px;'
        'font-weight:700}.label{font-size:13px}.value{font-size:12px;fill:#43506a}</style>',
        '<text x="28" y="36" class="title">IEDB OVA linear epitopes mapped to P01012</text>',
        '<text x="28" y="61" class="value">Exact unique sequence matches; 386-aa reference</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{left + axis_width}" y2="{axis_y}" '
        'stroke="#7b879d" stroke-width="5" stroke-linecap="round"/>',
        f'<text x="{left - 4}" y="{axis_y - 12}" class="value">1</text>',
        f'<text x="{left + axis_width - 22}" y="{axis_y - 12}" class="value">'
        f'{protein_length}</text>',
    ]
    colors = ["#457b9d", "#e76f51", "#2a9d8f", "#8f5aa8"]
    for index, row in enumerate(rows):
        y = 124 + 62 * index
        start = int(row["uniprot_start"])
        end = int(row["uniprot_end"])
        x = left + (start - 1) * scale
        bar_width = max(5.0, (end - start + 1) * scale)
        label = html.escape(
            f"IEDB {row['iedb_epitope_id']}  {row['linear_sequence']}  "
            f"P01012 {start}-{end}"
        )
        evidence = html.escape(
            f"{row['positive_assay_count']} / {row['total_assay_count']} positive rows; "
            f"{row['positive_reference_count']} references"
        )
        parts.extend([
            f'<rect x="{x:.2f}" y="{y}" width="{bar_width:.2f}" height="18" '
            f'rx="4" fill="{colors[index % len(colors)]}"/>',
            f'<text x="{left}" y="{y + 38}" class="label">{label}</text>',
            f'<text x="{left + 480}" y="{y + 38}" class="value">{evidence}</text>',
        ])
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_epitope_report(path: Path, rows: list[dict[str, object]]) -> None:
    """Write a traceable summary of the conservative IEDB epitope snapshot."""
    retrieved_dates = sorted({str(row["retrieved_date"]) for row in rows})
    lines = [
        "# Experimentally Observed OVA Epitope Report",
        "",
        "> This is a conservative, non-exhaustive IEDB snapshot. Positive assay "
        "records do not by themselves establish human clinical allergenicity.",
        "",
        "## Inclusion criteria",
        "",
        "- Linear peptide assigned by IEDB to ovalbumin / UniProt P01012.",
        "- At least one positive T-cell assay record.",
        "- Peptide sequence maps exactly and uniquely to the repository's P01012 sequence.",
        "- Canonical reference epitopes only; overlapping peptide-scan variants are excluded.",
        "",
        f"IEDB data retrieved through the official query API: "
        f"{', '.join(retrieved_dates)}.",
        "",
        "## Validated reference epitopes",
        "",
        "| IEDB epitope | Sequence | P01012 interval | Response | Most common MHC | "
        "All assay rows | Positive rows | Positive references |",
        "|---|---|---:|---|---|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| [IEDB {row['iedb_epitope_id']}]({row['source_url']}) | "
            f"`{row['linear_sequence']}` | {row['uniprot_start']}-{row['uniprot_end']} | "
            f"{row['immune_response_type']} | {row['most_common_positive_mhc']} | "
            f"{row['total_assay_count']} | {row['positive_assay_count']} | "
            f"{row['positive_reference_count']} |"
        )
    lines.extend([
        "",
        "Counts summarize positive rows and distinct reference IDs returned by the IEDB "
        "query API for each epitope ID. They are database-record counts, not effect sizes "
        "or counts of independent biological replications.",
        "",
        "## Interpretation limits",
        "",
        "- Most records use mouse hosts and mouse MHC contexts.",
        "- A T-cell response is not equivalent to IgE binding or clinical food allergy.",
        "- Negative records also exist and are not erased by reporting positive counts.",
        "- The snapshot is intentionally not an exhaustive catalog of overlapping peptides.",
        "- Database content can change after the recorded retrieval date.",
        "",
        "![IEDB OVA epitope sequence map](iedb_epitope_map.svg)",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_site_epitope_svg(
    path: Path, relationships: list[dict[str, object]]
) -> None:
    """Draw minimum Cα distances for modification-site/epitope pairs."""
    resolved = [
        row for row in relationships if row["min_ca_distance_angstrom"] != ""
    ]
    if not resolved:
        raise ValueError("Cannot draw site-to-epitope distances without coordinates")
    resolved.sort(key=lambda row: float(row["min_ca_distance_angstrom"]))
    width, height = 980, 115 + 44 * len(resolved)
    left, chart_width = 430, 430
    maximum = max(float(row["min_ca_distance_angstrom"]) for row in resolved)
    scale = chart_width / max(1.0, maximum * 1.08)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfcfe"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.title{font-size:20px;'
        'font-weight:700}.label{font-size:12px}.value{font-size:12px;fill:#43506a}</style>',
        '<text x="28" y="36" class="title">1OVA modification-site to IEDB epitope distance</text>',
        '<text x="28" y="61" class="value">Minimum Cα distance; chain A; overlap gives 0 Å</text>',
    ]
    for index, row in enumerate(resolved):
        y = 84 + 44 * index
        distance = float(row["min_ca_distance_angstrom"])
        label = html.escape(
            f"{row['reported_site']} / P01012 {row['site_uniprot_position']}  ↔  "
            f"IEDB {row['iedb_epitope_id']} ({row['epitope_uniprot_start']}-"
            f"{row['epitope_uniprot_end']})"
        )
        color = "#e76f51" if row["site_within_epitope"] else "#457b9d"
        bar_width = max(3.0, distance * scale)
        parts.extend([
            f'<text x="28" y="{y + 17}" class="label">{label}</text>',
            f'<rect x="{left}" y="{y}" width="{bar_width:.2f}" height="22" '
            f'rx="4" fill="{color}"/>',
            f'<text x="{left + bar_width + 8:.2f}" y="{y + 16}" '
            f'class="value">{distance:.1f} Å</text>',
        ])
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_site_epitope_report(
    path: Path,
    epitope_mappings: list[dict[str, object]],
    relationships: list[dict[str, object]],
) -> None:
    """Write the 1OVA modification-site/epitope geometry report."""
    resolved_relationships = [
        row for row in relationships if row["min_ca_distance_angstrom"] != ""
    ]
    lines = [
        "# OVA Modification-Site and Epitope Spatial Report",
        "",
        "> Distances are geometric descriptors for chain A of experimental structure "
        "1OVA. They do not demonstrate biochemical interaction or altered allergenicity.",
        "",
        "## Epitope coordinate coverage",
        "",
        "| IEDB epitope | P01012 interval | PDB interval | Resolved residues | Coverage |",
        "|---|---:|---|---:|---:|",
    ]
    for row in epitope_mappings:
        lines.append(
            f"| [IEDB {row['iedb_epitope_id']}]({row['source_url']}) | "
            f"{row['uniprot_start']}-{row['uniprot_end']} | "
            f"{row['pdb_start_residue_id']}-{row['pdb_end_residue_id']} "
            f"(chain {row['pdb_chain']}) | {row['resolved_residue_count']}/"
            f"{row['total_residue_count']} | {row['coordinate_coverage_percent']}% |"
        )
    lines.extend([
        "",
        "## Site-to-epitope relationships",
        "",
        "| Reported site | P01012 site | IEDB epitope | Sequence separation | "
        "Within epitope | Nearest epitope residue | Minimum Cα distance |",
        "|---|---:|---|---:|---|---|---:|",
    ])
    for row in sorted(
        resolved_relationships,
        key=lambda item: float(item["min_ca_distance_angstrom"]),
    ):
        lines.append(
            f"| {row['reported_site']} | {row['site_uniprot_position']} | "
            f"[IEDB {row['iedb_epitope_id']}]"
            f"(https://www.iedb.org/epitope/{row['iedb_epitope_id']}) | "
            f"{row['sequence_separation_residues']} residues | "
            f"{row['site_within_epitope']} | "
            f"{row['nearest_epitope_residue']}{row['nearest_epitope_uniprot_position']} "
            f"(PDB {row['nearest_epitope_pdb_residue_id']}) | "
            f"{row['min_ca_distance_angstrom']} Å |"
        )
    lines.extend([
        "",
        "The minimum distance is measured between the modification site's Cα atom and "
        "all resolved Cα atoms in the epitope interval. A site inside an epitope has a "
        "trivial self-distance of 0 Å; this is sequence overlap, not evidence of an "
        "interaction.",
        "",
        "## Interpretation limits",
        "",
        "- 1OVA is one crystallographic state and may not represent solution dynamics.",
        "- Cα proximity does not measure side-chain contact, binding or accessibility.",
        "- The IEDB snapshot is a conservative mouse T-cell reference set.",
        "- These results do not establish human IgE binding or clinical food allergy.",
        "",
        "![Modification-site to epitope distances](site_epitope_distances.svg)",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_sasa_svg(path: Path, site_rows: list[dict[str, object]]) -> None:
    """Draw absolute residue SASA for annotated modification sites."""
    if not site_rows:
        raise ValueError("Cannot draw an empty SASA chart")
    ordered = sorted(site_rows, key=lambda row: float(row["total_sasa_angstrom2"]))
    width, height = 900, 115 + 52 * len(ordered)
    left, chart_width = 330, 430
    maximum = max(float(row["total_sasa_angstrom2"]) for row in ordered)
    scale = chart_width / max(1.0, maximum * 1.1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfcfe"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.title{font-size:20px;'
        'font-weight:700}.label{font-size:13px}.value{font-size:12px;fill:#43506a}</style>',
        '<text x="28" y="36" class="title">1OVA annotated-site SASA</text>',
        '<text x="28" y="61" class="value">Shrake–Rupley; 1.4 Å probe; 960 points; chain A + NAG</text>',
    ]
    for index, row in enumerate(ordered):
        y = 84 + 52 * index
        total = float(row["total_sasa_angstrom2"])
        sidechain = float(row["sidechain_sasa_angstrom2"])
        label = html.escape(
            f"{row['residue']}-{row['reported_position']} → "
            f"P01012 {row['uniprot_position']} → PDB {row['pdb_residue_id']}"
        )
        parts.extend([
            f'<text x="28" y="{y + 18}" class="label">{label}</text>',
            f'<rect x="{left}" y="{y}" width="{total * scale:.2f}" height="24" '
            'rx="4" fill="#8fb9cf"/>',
            f'<rect x="{left}" y="{y}" width="{sidechain * scale:.2f}" height="24" '
            'rx="4" fill="#457b9d"/>',
            f'<text x="{left + total * scale + 8:.2f}" y="{y + 18}" class="value">'
            f'{total:.1f} Å²</text>',
        ])
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_sasa_report(
    path: Path,
    residue_rows: list[dict[str, object]],
    site_rows: list[dict[str, object]],
    epitope_rows: list[dict[str, object]],
) -> None:
    """Write the reproducible residue-level SASA report."""
    settings = residue_rows[0]
    lines = [
        "# OVA Solvent-Accessible Surface-Area Report",
        "",
        "## Method",
        "",
        "SASA was calculated with a deterministic Shrake–Rupley point-sampling "
        "implementation using Bondi-style elemental van der Waals radii.",
        "",
        f"- Structure context: {settings['occlusion_context']}",
        f"- Solvent probe radius: {settings['probe_radius_angstrom']} Å",
        f"- Sphere points per atom: {settings['sphere_point_count']}",
        f"- Resolved protein residues reported: {len(residue_rows)}",
        "- Hydrogen atoms are absent from the deposited X-ray model.",
        "",
        "Algorithm reference: [Shrake and Rupley (1973)]"
        "(https://doi.org/10.1016/0006-3495(73)90011-9). Atomic-radius reference: "
        "[Bondi (1964)](https://doi.org/10.1021/j100785a001).",
        "",
        "## Annotated modification sites",
        "",
        "| Reported site | P01012 | PDB | Total SASA | Side-chain SASA |",
        "|---|---:|---|---:|---:|",
    ]
    for row in sorted(site_rows, key=lambda item: int(item["uniprot_position"])):
        lines.append(
            f"| {row['residue']}-{row['reported_position']} | "
            f"{row['uniprot_position']} | {row['pdb_residue_name']}"
            f"{row['pdb_residue_id']} | {row['total_sasa_angstrom2']} Å² | "
            f"{row['sidechain_sasa_angstrom2']} Å² |"
        )
    lines.extend([
        "",
        "## IEDB epitope intervals",
        "",
        "| IEDB epitope | P01012 interval | Residues | Total SASA | Mean per residue | "
        "Side-chain SASA |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for row in epitope_rows:
        lines.append(
            f"| [IEDB {row['iedb_epitope_id']}]({row['source_url']}) | "
            f"{row['uniprot_start']}-{row['uniprot_end']} | "
            f"{row['resolved_residue_count']} | {row['total_sasa_angstrom2']} Å² | "
            f"{row['mean_residue_sasa_angstrom2']} Å² | "
            f"{row['sidechain_sasa_angstrom2']} Å² |"
        )
    lines.extend([
        "",
        "## Interpretation limits",
        "",
        "- Values are absolute SASA for one crystal structure, not relative SASA.",
        "- No exposed/buried threshold is assigned because reference maxima vary by method.",
        "- Crystal packing, missing hydrogens and conformational dynamics can alter SASA.",
        "- Accessibility does not establish chemical reactivity or immune recognition.",
        "",
        "![Annotated-site SASA](annotated_site_sasa.svg)",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_quality_report(path: Path, rows: list[dict[str, object]]) -> None:
    """Write a non-ranked evidence provenance and analytical-scope report."""
    counts = Counter(str(row["quality_class"]) for row in rows)
    lines = [
        "# Evidence Quality and Analytical-Scope Report",
        "",
        "> Classes describe provenance, granularity and reanalysis limits. They are "
        "not numeric scores and do not rank biological truth.",
        "",
        "## Classification summary",
        "",
        "| Quality class | Records | Meaning |",
        "|---|---:|---|",
        f"| limited_abstract_summary | {counts['limited_abstract_summary']} | "
        "Traceable numerical summary without raw observations |",
        f"| limited_abstract_annotation | {counts['limited_abstract_annotation']} | "
        "Traceable site annotation without complete site-level measurements |",
        f"| curated_experimental_database_aggregate | "
        f"{counts['curated_experimental_database_aggregate']} | "
        "Curated experimental database counts with heterogeneous, non-independent rows |",
        "",
        "## Record-level classification",
        "",
        "| Record | Source | Class | Granularity | Reanalysis readiness |",
        "|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['source_record_id']} | {row['source_identifier']} | "
            f"{row['quality_class']} | {row['data_granularity']} | "
            f"{row['quantitative_reanalysis_readiness']} |"
        )
    lines.extend([
        "",
        "## Safeguards",
        "",
        "- No p-values, effect sizes, sample sizes or uncertainty estimates are invented.",
        "- IEDB assay rows are not treated as independent biological replicates.",
        "- A traceable source can still have limited granularity or model relevance.",
        "- Classification does not replace formal study-level risk-of-bias assessment.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_report(
    path: Path,
    proteins: list[dict[str, object]],
    evidence: list[dict[str, object]],
    sites: list[dict[str, object]],
) -> None:
    lines = [
        "# OVA–Saccharide Public-Data Report",
        "",
        "> Numerical evidence is transcribed from publication abstracts. It is not "
        "reconstructed raw data.",
        "",
        "## Protein sequence characterization",
        "",
        "| Record | Length | Mass (kDa) | pI | GRAVY | N-X-S/T motifs | Lys + Arg |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for protein in proteins:
        lines.append(
            f"| {protein['identifier']} | {protein['length_aa']} | "
            f"{float(protein['molecular_weight_da']) / 1000:.2f} | "
            f"{protein['theoretical_pi']} | {protein['gravy']} | "
            f"{len(protein['n_glycosylation_sequons'])} | "
            f"{protein['candidate_basic_residues']} |"
        )
    lines.extend([
        "",
        "Calculated properties are estimates for the unmodified sequence. Motifs are "
        "sequence candidates, not confirmed glycosylation sites.",
        "",
        "## Publication-reported numerical evidence",
        "",
        "| Study | Treatment | Metric | Value | Unit | Comparator | Source |",
        "|---|---|---|---:|---|---|---|",
    ])
    for row in evidence:
        lines.append(
            f"| {row['study_id']} | {row['treatment']} | {row['metric']} | "
            f"{float(row['value']):g} | {row['unit']} | {row['comparator']} | "
            f"[DOI](https://doi.org/{row['doi']}) |"
        )
    lines.extend([
        "",
        "These are author-reported summary values. They cannot support new p-values "
        "or pooled effect sizes without the original measurements.",
        "",
        "## Publication-reported residue annotations",
        "",
        "| Study | Annotation | Reported site | UniProt site | Offset | Validation | Source |",
        "|---|---|---|---|---:|---|---|",
    ])
    for row in sites:
        lines.append(
            f"| {row['study_id']} | {row['annotation_type']} | "
            f"{row['residue']}-{row['reported_position']} | "
            f"{row['reference_residue']}-{row['uniprot_position']} | "
            f"{int(row['numbering_offset']):+d} | {row['sequence_validation']} | "
            f"[DOI](https://doi.org/{row['doi']}) |"
        )
    lines.extend([
        "",
        "Residue numbering is retained as reported by each publication. A one-position "
        "difference can occur relative to the 386-aa UniProt sequence because of "
        "initiator-residue or mature-chain conventions.",
        "",
        "## Interpretation boundary",
        "",
        "The report is descriptive. It does not establish causality, clinical "
        "allergenicity, safety, or therapeutic efficacy.",
        "",
        "![Hwang 2014 relative response](hwang_2014_relative_response.svg)",
        "",
        "![Wang 2013 glycation extent](wang_2013_glycation_extent.svg)",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_structure_svg(path: Path, rows: list[dict[str, object]]) -> None:
    resolved = [row for row in rows if row["coordinate_resolved"]]
    width, height = 900, 120 + 54 * len(resolved)
    left, chart_width = 330, 460
    maximum = max(float(row["distance_to_nag_angstrom"]) for row in resolved)
    scale = chart_width / (maximum * 1.12)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfcfe"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.title{font-size:20px;'
        'font-weight:700}.label{font-size:13px}.value{font-size:13px;fill:#43506a}</style>',
        '<text x="28" y="36" class="title">1OVA modification-site distance to NAG</text>',
        '<text x="28" y="61" class="value">Minimum heavy-atom distance; chain A</text>',
    ]
    for index, row in enumerate(resolved):
        y = 88 + 54 * index
        distance = float(row["distance_to_nag_angstrom"])
        label = html.escape(
            f"{row['residue']}-{row['reported_position']} → "
            f"{row['reference_residue']}-{row['uniprot_position']} → "
            f"PDB {row['pdb_residue_name']}{row['pdb_residue_id']}"
        )
        color = "#e76f51" if row["annotation_type"] == "N-linked glycosylation" else "#457b9d"
        parts.extend([
            f'<text x="28" y="{y + 18}" class="label">{label}</text>',
            f'<rect x="{left}" y="{y}" width="{distance * scale:.2f}" height="24" '
            f'rx="5" fill="{color}"/>',
            f'<text x="{left + distance * scale + 9:.2f}" y="{y + 18}" '
            f'class="value">{distance:.1f} Å</text>',
        ])
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_pymol_script(
    path: Path,
    rows: list[dict[str, object]],
    epitope_mappings: list[dict[str, object]] | None = None,
) -> None:
    glycosylation_ids = [
        str(row["pdb_residue_id"])
        for row in rows if row["annotation_type"] == "N-linked glycosylation"
    ]
    glycation_ids = [
        str(row["pdb_residue_id"])
        for row in rows if row["annotation_type"] == "glycation-associated site"
    ]
    all_ids = glycosylation_ids + glycation_ids
    lines = [
        "# PyMOL visualization for the public-data OVA structure analysis",
        "# Run from the repository root: pymol scripts/visualize_1ova_sites.pml",
        "reinitialize",
        "load structures/1OVA.pdb, ova_1ova",
        "hide everything",
        "show cartoon, ova_1ova and chain A",
        "color gray70, ova_1ova and chain A",
        "set cartoon_transparency, 0.12",
        "show sticks, ova_1ova and chain A and resn NAG",
        "color tv_orange, ova_1ova and chain A and resn NAG",
    ]
    if glycation_ids:
        lines.extend([
            f"select glycation_sites, ova_1ova and chain A and resi {'+'.join(glycation_ids)}",
            "show sticks, glycation_sites",
            "color marine, glycation_sites",
        ])
    if glycosylation_ids:
        lines.extend([
            f"select n_glycosylation_site, ova_1ova and chain A and resi "
            f"{'+'.join(glycosylation_ids)}",
            "show sticks, n_glycosylation_site",
            "color firebrick, n_glycosylation_site",
        ])
    if all_ids:
        lines.extend([
            f"select annotated_sites, ova_1ova and chain A and resi {'+'.join(all_ids)}",
            'label annotated_sites and name CA, "%s%s" % (resn, resi)',
        ])
    for index, epitope in enumerate(epitope_mappings or [], start=1):
        lines.extend([
            f"select iedb_epitope_{epitope['iedb_epitope_id']}, ova_1ova and chain A "
            f"and resi {epitope['pdb_start_residue_id']}-{epitope['pdb_end_residue_id']}",
            f"color {'teal' if index % 2 else 'violet'}, "
            f"iedb_epitope_{epitope['iedb_epitope_id']}",
        ])
    lines.extend([
        "show surface, ova_1ova and chain A",
        "set transparency, 0.78, ova_1ova and chain A",
        "set label_size, 16",
        "set label_color, black",
        "orient ova_1ova and chain A",
        "bg_color white",
        "ray 1600, 1200",
        "png outputs/1ova_annotated_sites.png, dpi=200",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_structure_report(
    path: Path,
    metadata: dict[str, object],
    rows: list[dict[str, object]],
    pairwise: list[dict[str, object]],
) -> None:
    lines = [
        "# OVA Three-Dimensional Structure Report",
        "",
        "## Structure source",
        "",
        f"- PDB: [{metadata['pdb_id']}](https://www.rcsb.org/structure/{metadata['pdb_id']})",
        f"- PDB DOI: [{metadata['pdb_doi']}](https://doi.org/{metadata['pdb_doi']})",
        f"- Method: {metadata['experimental_method']}",
        f"- Resolution: {metadata['resolution_angstrom']} Å",
        f"- Model: {metadata['title']}",
        "",
        "Analysis below uses chain A. UniProt positions were mapped through the "
        "mmCIF polymer-sequence scheme rather than by assuming a constant numbering offset.",
        "",
        "## Modification-site mapping",
        "",
        "| Reported site | UniProt site | PDB site | Resolved | Distance to NAG | "
        "CA radius | Mean B-factor | Source |",
        "|---|---|---|---|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['residue']}-{row['reported_position']} | "
            f"{row['reference_residue']}-{row['uniprot_position']} | "
            f"{row['pdb_residue_name']}-{row['pdb_residue_id']} (chain {row['pdb_chain']}) | "
            f"{row['coordinate_resolved']} | {row['distance_to_nag_angstrom']} Å | "
            f"{row['ca_distance_from_chain_centroid_angstrom']} Å | "
            f"{row['mean_b_factor']} | [DOI](https://doi.org/{row['doi']}) |"
        )
    lines.extend([
        "",
        "The distance to NAG is the minimum heavy-atom distance to the crystallographic "
        "N-acetylglucosamine in the same chain. CA radius is distance from the chain-A "
        "Cα centroid; it is not a solvent-accessibility calculation.",
        "",
        "## Closest annotated-site pairs",
        "",
        "| Site A | Site B | Cα distance |",
        "|---|---|---:|",
    ])
    for row in sorted(pairwise, key=lambda item: float(item["ca_distance_angstrom"]))[:10]:
        lines.append(
            f"| {row['site_a']} | {row['site_b']} | "
            f"{row['ca_distance_angstrom']} Å |"
        )
    lines.extend([
        "",
        "## Interpretation limits",
        "",
        "- Crystal coordinates describe one experimental structural state.",
        "- A short geometric distance does not demonstrate biochemical interaction.",
        "- B-factors are structure- and refinement-dependent.",
        "- Epitope proximity is reported separately in `site_epitope_report.md`.",
        "- Solvent accessibility requires a separate analysis.",
        "",
        "![Distance of annotated sites to NAG](1ova_site_distance_to_nag.svg)",
        "",
        "A PyMOL script is provided at `scripts/visualize_1ova_sites.pml`.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")
