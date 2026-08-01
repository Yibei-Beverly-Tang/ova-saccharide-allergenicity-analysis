import unittest
from pathlib import Path

from ova_analysis.evidence import (
    map_sites_to_sequence,
    read_literature_values,
    read_site_annotations,
)
from ova_analysis.epitopes import map_epitopes_to_sequence, read_iedb_epitopes
from ova_analysis.protein import read_fasta
from ova_analysis.quality import classify_evidence


ROOT = Path(__file__).resolve().parents[1]


class QualityTests(unittest.TestCase):
    def test_all_public_evidence_is_classified_without_numeric_ranking(self):
        sequence = read_fasta(ROOT / "sequences/ova_uniprot_P01012.fasta")[0].sequence
        evidence = read_literature_values(ROOT / "data/public_literature_values.csv")
        sites = map_sites_to_sequence(
            read_site_annotations(ROOT / "data/public_site_annotations.csv"), sequence
        )
        epitopes = map_epitopes_to_sequence(
            read_iedb_epitopes(ROOT / "data/public_iedb_epitopes.csv"), sequence
        )
        rows = classify_evidence(evidence, sites, epitopes)
        self.assertEqual(len(rows), len(evidence) + len(sites) + len(epitopes))
        self.assertEqual(len({row["source_record_id"] for row in rows}), len(rows))
        self.assertEqual(
            {row["quality_class"] for row in rows},
            {
                "limited_abstract_summary",
                "limited_abstract_annotation",
                "curated_experimental_database_aggregate",
            },
        )
        self.assertTrue(all("score" not in row for row in rows))


if __name__ == "__main__":
    unittest.main()
