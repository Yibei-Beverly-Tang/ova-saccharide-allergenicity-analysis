# Changelog

## [0.3.1] - 2026-08-02

### Added

- GitHub Actions validation on Python 3.10 and 3.12.
- End-to-end test coverage for the public-data command-line workflow.
- Optional dependency groups for experimental analysis and testing.

### Changed

- Added DOI, PMID, year, finite-value and percentage-range validation for
  literature evidence.
- Added validation of site-annotation evidence levels.
- Made the generated PyMOL script path configurable.
- Synchronized the package version with project metadata.
- Added ignore rules for caches, build artifacts and designated private-data
  directories.

## [0.3.0] - 2026-07-31

### Added

- Experimental 1OVA PDB and mmCIF structure data with checksums.
- UniProt-to-PDB residue mapping through the mmCIF polymer-sequence scheme.
- Coordinate-resolution checks for all publication-reported sites.
- Site-to-NAG and pairwise Cα distance calculations.
- Chain-centroid radius and mean B-factor descriptors.
- Structure CSV outputs, Markdown report and SVG visualization.
- PyMOL visualization script.
- Four structure-specific automated tests.

### Scientific safeguards

- Structure mapping does not assume a constant residue-number offset.
- Geometric proximity is not interpreted as biochemical interaction.
- Chain-centroid radius is not represented as solvent accessibility.

## [0.2.0] - 2026-07-31

### Added

- Real UniProtKB P01012 sequence analysis.
- Traceable publication-reported numerical evidence.
- Publication-to-UniProt residue validation.
- Privacy and evidence-provenance safeguards.
