import math
import unittest

from ova_analysis.sasa import atom_sasa
from ova_analysis.structure import Atom


class SasaTests(unittest.TestCase):
    def test_isolated_carbon_matches_analytic_sphere_area(self):
        atom = Atom(
            record="ATOM",
            atom_name="CA",
            residue_name="ALA",
            chain="A",
            residue_number=1,
            insertion_code="",
            x=0.0,
            y=0.0,
            z=0.0,
            b_factor=0.0,
            element="C",
        )
        observed = atom_sasa([atom], probe_radius=1.4, point_count=100)[0]
        expected = 4.0 * math.pi * (1.70 + 1.4) ** 2
        self.assertAlmostEqual(observed, expected, places=10)

    def test_overlapping_atoms_reduce_accessible_area(self):
        atoms = [
            Atom("ATOM", "CA", "ALA", "A", 1, "", 0, 0, 0, 0, "C"),
            Atom("ATOM", "CB", "ALA", "A", 1, "", 3, 0, 0, 0, "C"),
        ]
        areas = atom_sasa(atoms, probe_radius=1.4, point_count=960)
        isolated = 4.0 * math.pi * (1.70 + 1.4) ** 2
        self.assertTrue(all(0 < area < isolated for area in areas))


if __name__ == "__main__":
    unittest.main()
