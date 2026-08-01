import tempfile
import unittest
from pathlib import Path

from ova_analysis.evidence import (
    map_sites_to_sequence,
    read_literature_values,
    read_site_annotations,
)


class EvidenceTests(unittest.TestCase):
    def test_reads_traceable_summary(self):
        content = (
            "study_id,publication_year,doi,pmid,treatment,metric,value,unit,"
            "comparator,evidence_level,source_location,notes\n"
            "S,2014,10.1000/example,,Treatment,metric,58.8,percent,control,"
            "abstract_reported_summary,Abstract,note\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "evidence.csv"
            path.write_text(content, encoding="utf-8")
            rows = read_literature_values(path)
        self.assertEqual(rows[0]["value"], 58.8)

    def test_rejects_unknown_evidence_level(self):
        content = (
            "study_id,publication_year,doi,pmid,treatment,metric,value,unit,"
            "comparator,evidence_level,source_location,notes\n"
            "S,2014,10.1000/example,,Treatment,metric,58.8,percent,control,"
            "unknown,Abstract,note\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "evidence.csv"
            path.write_text(content, encoding="utf-8")
            with self.assertRaises(ValueError):
                read_literature_values(path)

    def test_rejects_non_finite_percentage(self):
        content = (
            "study_id,publication_year,doi,pmid,treatment,metric,value,unit,"
            "comparator,evidence_level,source_location,notes\n"
            "S,2014,10.1000/example,,Treatment,metric,nan,percent,control,"
            "abstract_reported_summary,Abstract,note\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "evidence.csv"
            path.write_text(content, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Non-finite"):
                read_literature_values(path)

    def test_rejects_percentage_above_100(self):
        content = (
            "study_id,publication_year,doi,pmid,treatment,metric,value,unit,"
            "comparator,evidence_level,source_location,notes\n"
            "S,2014,10.1000/example,,Treatment,sequence_coverage,101,percent,control,"
            "abstract_reported_summary,Abstract,note\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "evidence.csv"
            path.write_text(content, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "above 100"):
                read_literature_values(path)

    def test_allows_relative_response_above_100(self):
        content = (
            "study_id,publication_year,doi,pmid,treatment,metric,value,unit,"
            "comparator,evidence_level,source_location,notes\n"
            "S,2014,10.1000/example,,Treatment,relative_IgE_production,125,"
            "percent of control,control,abstract_reported_summary,Abstract,note\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "evidence.csv"
            path.write_text(content, encoding="utf-8")
            rows = read_literature_values(path)
        self.assertEqual(rows[0]["value"], 125.0)

    def test_reads_site_annotation(self):
        content = (
            "study_id,publication_year,doi,pmid,annotation_type,residue,"
            "reported_position,reference_context,evidence_level,notes\n"
            "S,2007,10.1000/example,1,N-linked glycosylation,Asn,292,"
            "study numbering,abstract_reported_annotation,note\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sites.csv"
            path.write_text(content, encoding="utf-8")
            rows = read_site_annotations(path)
        self.assertEqual(rows[0]["reported_position"], 292)

    def test_rejects_invalid_site_doi(self):
        content = (
            "study_id,publication_year,doi,pmid,annotation_type,residue,"
            "reported_position,reference_context,evidence_level,notes\n"
            "S,2007,not-a-doi,1,N-linked glycosylation,Asn,292,"
            "study numbering,abstract_reported_annotation,note\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sites.csv"
            path.write_text(content, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Invalid DOI"):
                read_site_annotations(path)

    def test_maps_publication_numbering_offset(self):
        rows = [{
            "study_id": "S",
            "residue": "Asn",
            "reported_position": 2,
        }]
        mapped = map_sites_to_sequence(rows, "AAN")
        self.assertEqual(mapped[0]["uniprot_position"], 3)
        self.assertEqual(mapped[0]["numbering_offset"], 1)
        self.assertEqual(mapped[0]["sequence_validation"], "matched")


if __name__ == "__main__":
    unittest.main()
