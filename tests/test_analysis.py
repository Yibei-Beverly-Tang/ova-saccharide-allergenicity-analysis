import sys
import unittest
from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise unittest.SkipTest(
        "experimental-analysis tests require the project test extra"
    ) from exc


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analyze_experiment import compute_outcomes, validate_data  # noqa: E402


def row(sample_id, saccharide, free_amino=100.0, ige=80.0):
    return {
        "sample_id": sample_id,
        "branch": "Test branch",
        "treatment": saccharide,
        "saccharide": saccharide,
        "replicate": 1,
        "free_amino_groups_pct": free_amino,
        "ige_binding_pct": ige,
        "particle_size_nm": 100.0,
        "surface_hydrophobicity_au": 40.0,
    }


class ComputeOutcomesTest(unittest.TestCase):
    def test_control_relative_calculation(self):
        data = pd.DataFrame(
            [row("C1", "Control"), row("T1", "Lactose", 50.0, 40.0)]
        )
        treated = compute_outcomes(data).query("sample_id == 'T1'").iloc[0]
        self.assertAlmostEqual(treated["modification_extent_pct"], 50.0)
        self.assertAlmostEqual(treated["ige_reduction_vs_control_pct"], 50.0)

    def test_missing_control_raises_error(self):
        with self.assertRaises(ValueError):
            compute_outcomes(pd.DataFrame([row("T1", "Lactose")]))

    def test_duplicate_sample_id_raises_error(self):
        with self.assertRaises(ValueError):
            validate_data(
                pd.DataFrame([row("DUPLICATE", "Control")] * 2)
            )

    def test_missing_column_raises_error(self):
        with self.assertRaises(ValueError):
            validate_data(pd.DataFrame([{"sample_id": "C1"}]))


if __name__ == "__main__":
    unittest.main()
