"""Transparent evidence-scope classification without invented quality scores."""

from __future__ import annotations

from collections import Counter


QUALITY_FIELDS = [
    "source_record_id",
    "source_type",
    "source_identifier",
    "evidence_level",
    "quality_class",
    "traceability",
    "data_granularity",
    "replication_status",
    "quantitative_reanalysis_readiness",
    "primary_limitation",
]


def classify_evidence(
    literature_values: list[dict[str, object]],
    site_annotations: list[dict[str, object]],
    epitopes: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Classify provenance and analytical scope without ranking biological truth."""
    rows: list[dict[str, object]] = []
    for value in literature_values:
        rows.append({
            "source_record_id": (
                f"literature_value:{value['study_id']}:{value['metric']}:"
                f"{value['treatment']}"
            ),
            "source_type": "publication numerical evidence",
            "source_identifier": f"DOI:{value['doi']}",
            "evidence_level": value["evidence_level"],
            "quality_class": "limited_abstract_summary",
            "traceability": "DOI and source location retained",
            "data_granularity": "author-reported abstract summary",
            "replication_status": "individual replicates not available",
            "quantitative_reanalysis_readiness": "descriptive use only",
            "primary_limitation": (
                "No raw observations, sample size or dispersion estimate in source extract"
            ),
        })
    for site in site_annotations:
        rows.append({
            "source_record_id": (
                f"site_annotation:{site['study_id']}:{site['residue']}"
                f"{site['reported_position']}"
            ),
            "source_type": "publication residue annotation",
            "source_identifier": f"DOI:{site['doi']}",
            "evidence_level": site["evidence_level"],
            "quality_class": "limited_abstract_annotation",
            "traceability": "DOI and publication numbering retained",
            "data_granularity": "author-reported abstract annotation",
            "replication_status": "site-level replication not reported",
            "quantitative_reanalysis_readiness": "mapping and descriptive use only",
            "primary_limitation": (
                "Abstract annotation does not provide complete site-level measurements"
            ),
        })
    for epitope in epitopes:
        rows.append({
            "source_record_id": f"iedb_epitope:{epitope['iedb_epitope_id']}",
            "source_type": "IEDB experimental-record aggregate",
            "source_identifier": f"IEDB:{epitope['iedb_epitope_id']}",
            "evidence_level": epitope["evidence_level"],
            "quality_class": "curated_experimental_database_aggregate",
            "traceability": "IEDB epitope ID, query URL and retrieval date retained",
            "data_granularity": "positive assay-row and distinct-reference counts",
            "replication_status": (
                "reference count available; independence not established"
            ),
            "quantitative_reanalysis_readiness": "descriptive counts only",
            "primary_limitation": (
                "Database rows are heterogeneous and are not independent effect estimates"
            ),
        })
    identifiers = [str(row["source_record_id"]) for row in rows]
    duplicates = [key for key, count in Counter(identifiers).items() if count > 1]
    if duplicates:
        raise ValueError(f"Duplicate evidence quality record IDs: {duplicates}")
    return rows
