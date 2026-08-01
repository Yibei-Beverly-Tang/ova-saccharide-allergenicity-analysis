import unittest
from pathlib import Path

from ova_analysis.evidence import map_sites_to_sequence, read_site_annotations
from ova_analysis.epitopes import map_epitopes_to_sequence, read_iedb_epitopes
from ova_analysis.protein import read_fasta
from ova_analysis.structure import (
    analyze_structure_sites,
    analyze_epitope_structure_relationships,
    read_atoms,
    read_poly_seq_scheme,
    read_structure_metadata,
    validate_chain_sequence,
)


ROOT = Path(__file__).resolve().parents[1]


class StructureTests(unittest.TestCase):
    def test_real_1ova_metadata(self):
        metadata = read_structure_metadata(ROOT / "structures/1OVA.pdb")
        self.assertEqual(metadata["pdb_id"], "1OVA")
        self.assertEqual(metadata["resolution_angstrom"], 1.95)

    def test_real_structure_matches_uniprot_sequence(self):
        record = read_fasta(ROOT / "sequences/ova_uniprot_P01012.fasta")[0]
        scheme = read_poly_seq_scheme(ROOT / "structures/1OVA.cif")
        validate_chain_sequence(scheme, record.sequence, "A")

    def test_real_sites_map_to_coordinates(self):
        record = read_fasta(ROOT / "sequences/ova_uniprot_P01012.fasta")[0]
        annotations = read_site_annotations(ROOT / "data/public_site_annotations.csv")
        sites = map_sites_to_sequence(annotations, record.sequence)
        mapped, pairwise = analyze_structure_sites(
            sites,
            record.sequence,
            ROOT / "structures/1OVA.pdb",
            ROOT / "structures/1OVA.cif",
        )
        self.assertEqual(len(mapped), 7)
        self.assertTrue(all(row["coordinate_resolved"] for row in mapped))
        self.assertEqual(len(pairwise), 21)
        glycosylation = next(row for row in mapped if row["study_id"] == "Ito_2007")
        self.assertEqual(glycosylation["pdb_residue_id"], "298")
        self.assertLess(float(glycosylation["distance_to_nag_angstrom"]), 2.0)

    def test_pdb_contains_four_nag_residues(self):
        atoms = read_atoms(ROOT / "structures/1OVA.pdb")
        residues = {
            (atom.chain, atom.residue_number)
            for atom in atoms
            if atom.residue_name == "NAG"
        }
        self.assertEqual(residues, {("A", 393), ("B", 393), ("C", 393), ("D", 393)})

    def test_epitope_site_relationships_use_resolved_coordinates(self):
        record = read_fasta(ROOT / "sequences/ova_uniprot_P01012.fasta")[0]
        sites = map_sites_to_sequence(
            read_site_annotations(ROOT / "data/public_site_annotations.csv"),
            record.sequence,
        )
        structure_sites, _ = analyze_structure_sites(
            sites,
            record.sequence,
            ROOT / "structures/1OVA.pdb",
            ROOT / "structures/1OVA.cif",
        )
        epitopes = map_epitopes_to_sequence(
            read_iedb_epitopes(ROOT / "data/public_iedb_epitopes.csv"),
            record.sequence,
        )
        mapped, relationships = analyze_epitope_structure_relationships(
            epitopes,
            structure_sites,
            record.sequence,
            ROOT / "structures/1OVA.pdb",
            ROOT / "structures/1OVA.cif",
        )
        self.assertEqual(len(mapped), 2)
        self.assertTrue(all(row["coordinate_coverage_percent"] == 100.0 for row in mapped))
        self.assertEqual(len(relationships), 14)
        overlap = [row for row in relationships if row["site_within_epitope"]]
        self.assertEqual(len(overlap), 1)
        self.assertEqual(overlap[0]["site_uniprot_position"], 264)
        self.assertEqual(overlap[0]["iedb_epitope_id"], 58560)
        self.assertEqual(overlap[0]["min_ca_distance_angstrom"], 0.0)
        glycan_to_class_ii = next(
            row for row in relationships
            if row["site_uniprot_position"] == 293
            and row["iedb_epitope_id"] == 28676
        )
        self.assertEqual(glycan_to_class_ii["nearest_epitope_uniprot_position"], 328)
        self.assertEqual(glycan_to_class_ii["min_ca_distance_angstrom"], 4.957)


if __name__ == "__main__":
    unittest.main()
