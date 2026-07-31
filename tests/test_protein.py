import tempfile
import unittest
from pathlib import Path

from ova_analysis.protein import analyze_record, n_glycosylation_sequons, read_fasta


class ProteinTests(unittest.TestCase):
    def test_fasta_analysis(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "protein.fasta"
            path.write_text(">test\nACDEFGHIKLMNPQRSTVWY\n", encoding="utf-8")
            result = analyze_record(read_fasta(path)[0])
        self.assertEqual(result["length_aa"], 20)
        self.assertGreater(result["molecular_weight_da"], 2000)
        self.assertTrue(0 <= result["theoretical_pi"] <= 14)

    def test_motif_scan_excludes_proline(self):
        self.assertEqual(
            n_glycosylation_sequons("ANVTNPSNAT"),
            [{"position": 2, "motif": "NVT"}, {"position": 8, "motif": "NAT"}],
        )


if __name__ == "__main__":
    unittest.main()

