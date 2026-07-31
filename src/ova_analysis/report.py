"""Writers for machine-readable summaries, Markdown, and SVG charts."""

from __future__ import annotations

import csv
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


def write_pymol_script(path: Path, rows: list[dict[str, object]]) -> None:
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
        "- Epitope proximity and solvent accessibility require separate analyses.",
        "",
        "![Distance of annotated sites to NAG](1ova_site_distance_to_nag.svg)",
        "",
        "A PyMOL script is provided at `scripts/visualize_1ova_sites.pml`.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")
