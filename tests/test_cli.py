import tempfile
import unittest
from pathlib import Path

from ova_analysis.cli import main


ROOT = Path(__file__).resolve().parents[1]


class CliTests(unittest.TestCase):
    def test_complete_public_data_workflow(self):
        with tempfile.TemporaryDirectory() as directory:
            temporary = Path(directory)
            output = temporary / "outputs"
            pymol_script = temporary / "visualize.pml"
            result = main([
                "--fasta", str(ROOT / "sequences/ova_uniprot_P01012.fasta"),
                "--evidence", str(ROOT / "data/public_literature_values.csv"),
                "--sites", str(ROOT / "data/public_site_annotations.csv"),
                "--epitopes", str(ROOT / "data/public_iedb_epitopes.csv"),
                "--structure-pdb", str(ROOT / "structures/1OVA.pdb"),
                "--structure-cif", str(ROOT / "structures/1OVA.cif"),
                "--out-dir", str(output),
                "--pymol-script", str(pymol_script),
            ])

            self.assertEqual(result, 0)
            expected = {
                "protein_summary.csv",
                "literature_evidence.csv",
                "validated_site_annotations.csv",
                "hwang_2014_relative_response.svg",
                "wang_2013_glycation_extent.svg",
                "report.md",
                "validated_iedb_epitopes.csv",
                "iedb_epitope_map.svg",
                "epitope_report.md",
                "structure_site_mapping.csv",
                "structure_site_pairwise_distances.csv",
                "1ova_site_distance_to_nag.svg",
                "structure_report.md",
            }
            self.assertEqual({path.name for path in output.iterdir()}, expected)
            self.assertTrue(pymol_script.is_file())


if __name__ == "__main__":
    unittest.main()
