# Changelog

## [1.0.0] - 2026-08-02

### Added

- Deterministic Shrake–Rupley solvent-accessible surface-area analysis with a
  1.4 Å probe, 960 sphere points per atom and documented elemental radii.
- Residue-level total, backbone and side-chain SASA for 383 resolved 1OVA
  chain-A residues with same-chain NAG occlusion.
- SASA tables and summaries for all annotated modification sites and both IEDB
  epitope intervals, plus a Markdown report and SVG chart.
- Record-level evidence classification by provenance, granularity, replication
  status, reanalysis readiness and primary limitation.
- Analytical sphere-area and atom-occlusion tests for the SASA implementation.

### Scientific safeguards

- Absolute SASA is not converted to relative SASA or exposed/buried labels
  without method-specific reference maxima.
- Accessibility is not interpreted as chemical reactivity or immune recognition.
- Evidence classes are descriptive and non-ranked; no numerical quality score
  or unavailable study statistic is invented.

### Milestone

- Completed every item in the original project roadmap.

## [0.5.0] - 2026-08-02

### Added

- Complete residue-level mapping of both validated IEDB epitope intervals to
  1OVA chain A with coordinate-coverage reporting.
- Fourteen modification-site/epitope relationships with sequence separation,
  nearest epitope residue and minimum Cα distance.
- Machine-readable structure and distance tables, Markdown report, SVG chart
  and PyMOL epitope selections.
- Regression checks for the exact K264 sequence overlap and the independently
  calculated 4.957 Å Asn293-to-IEDB-28676 minimum Cα distance.

### Scientific safeguards

- A 0 Å self-distance is explicitly labeled as sequence overlap rather than an
  interaction.
- Cα proximity is not interpreted as binding, accessibility, immune modulation
  or clinical allergenicity.

## [0.4.0] - 2026-08-02

### Added

- Conservative IEDB snapshot for canonical OVA T-cell epitopes 58560 and
  28676, including positive assay-row and distinct-reference counts.
- Exact, unique peptide-to-UniProt P01012 sequence mapping.
- Machine-readable validated epitope output, Markdown limitations report and
  SVG sequence map.
- Epitope validation and end-to-end workflow tests.
- Reproducible IEDB API updater with strict source-antigen and sequence checks.

### Scientific safeguards

- Overlapping peptide-scan variants are not represented as independent
  validated epitopes.
- Positive IEDB rows are not represented as effect sizes or independent
  biological replications.
- Mouse T-cell evidence is not interpreted as human IgE binding or clinical
  food allergy.

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
