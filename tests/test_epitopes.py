import tempfile
import unittest
from pathlib import Path

from ova_analysis.epitopes import map_epitopes_to_sequence, read_iedb_epitopes
from ova_analysis.protein import read_fasta


ROOT = Path(__file__).resolve().parents[1]


class EpitopeTests(unittest.TestCase):
    def test_public_snapshot_maps_to_p01012(self):
        sequence = read_fasta(ROOT / "sequences/ova_uniprot_P01012.fasta")[0].sequence
        rows = read_iedb_epitopes(ROOT / "data/public_iedb_epitopes.csv")
        mapped = map_epitopes_to_sequence(rows, sequence)
        positions = {
            row["iedb_epitope_id"]: (row["uniprot_start"], row["uniprot_end"])
            for row in mapped
        }
        self.assertEqual(positions[58560], (258, 265))
        self.assertEqual(positions[28676], (324, 340))
        counts = {
            row["iedb_epitope_id"]: (
                row["total_assay_count"],
                row["positive_assay_count"],
                row["positive_reference_count"],
            )
            for row in mapped
        }
        self.assertEqual(counts[58560], (485, 471, 191))
        self.assertEqual(counts[28676], (188, 170, 58))
        self.assertTrue(
            all(row["sequence_validation"] == "exact_unique_match" for row in mapped)
        )

    def test_rejects_zero_positive_assays(self):
        content = (
            "iedb_epitope_id,linear_sequence,immune_response_type,structure_type,"
            "total_assay_count,positive_assay_count,positive_reference_count,"
            "most_common_positive_host,most_common_positive_mhc,evidence_level,"
            "source_url,query_url,"
            "retrieved_date,notes\n"
            "1,SIINFEKL,T cell,Linear peptide,1,0,1,mouse,H2-Kb,"
            "iedb_positive_assay_summary,https://www.iedb.org/epitope/1,"
            "https://query-api.iedb.org/api/v1/tcell_search?structure_id=eq.1,"
            "2026-08-02,note\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "epitopes.csv"
            path.write_text(content, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "evidence counts"):
                read_iedb_epitopes(path)

    def test_rejects_ambiguous_sequence_mapping(self):
        rows = [{"iedb_epitope_id": 1, "linear_sequence": "ABC"}]
        with self.assertRaisesRegex(ValueError, "ambiguously"):
            map_epitopes_to_sequence(rows, "ABCABC")


if __name__ == "__main__":
    unittest.main()
