# OVA–Saccharide Allergenicity Analysis

A reproducible Python project combining **ovalbumin protein characterization**
and **experimental three-dimensional structure analysis** with traceable
publication-reported evidence about glycation, N-glycan modification and
immune-response measurements.

This public version contains no résumé, manuscript draft, unpublished result or
private experimental dataset. It also contains no invented replicate values.

## Included real public data

- Reviewed chicken ovalbumin sequence:
  [UniProtKB P01012](https://www.uniprot.org/uniprotkb/P01012/entry), 386 aa.
- Wang et al. (2013), DOI
  [10.1016/j.foodchem.2013.04.045](https://doi.org/10.1016/j.foodchem.2013.04.045):
  abstract-reported sequence coverage and percentages of glucose-modified
  lysines after microwave or conventional heating.
- Hwang et al. (2014), DOI
  [10.1016/j.bbrc.2014.06.101](https://doi.org/10.1016/j.bbrc.2014.06.101):
  abstract-reported IgE, IL-4 and IL-5 production relative to intact OVA after
  terminal N-acetylglucosamine cleavage.
- Ito et al. (2007), DOI
  [10.1093/glycob/cwl077](https://doi.org/10.1093/glycob/cwl077):
  abstract-reported OVA N-glycosylation annotation.
- Mao et al. (2023), DOI
  [10.1016/j.ijbiomac.2023.123640](https://doi.org/10.1016/j.ijbiomac.2023.123640):
  abstract-reported glycation-associated residues.
- RCSB PDB [1OVA](https://www.rcsb.org/structure/1OVA), DOI
  [10.2210/pdb1OVA/pdb](https://doi.org/10.2210/pdb1OVA/pdb):
  X-ray structure of uncleaved ovalbumin at 1.95 Å resolution.
- IEDB reference T-cell epitopes
  [58560](https://www.iedb.org/epitope/58560) (`SIINFEKL`) and
  [28676](https://www.iedb.org/epitope/28676) (`ISQAVHAAHAEINEAGR`),
  retained as a conservative snapshot of positive experimental assay records.

Every numerical row retains its DOI, PMID when available, comparator, unit,
source location and evidence level.

## Analysis

Protein sequence module:

- FASTA validation;
- length and unmodified molecular-mass estimate;
- theoretical pI, GRAVY, aromaticity and aliphatic index;
- Lys/Arg abundance;
- canonical `N-X-S/T` sequon scanning (`X ≠ P`).

Evidence module:

- validates public literature-summary tables;
- rejects records without traceable study metadata;
- preserves publication residue numbering;
- creates source-specific charts without pooling incompatible experiments.

Epitope module:

- validates a dated, conservative snapshot of positive IEDB T-cell evidence;
- rejects missing provenance, duplicate IDs and invalid evidence counts;
- maps each peptide by an exact, unique match to UniProt P01012;
- preserves the distinction between assay rows, references and biological
  replication;
- generates a machine-readable table, sequence map and limitations report.

The committed IEDB snapshot can be regenerated from the official query API
with `scripts/update_iedb_epitopes.py`; no epitope counts are entered by
estimation or imputation.

Three-dimensional structure module:

- validates the 1OVA mmCIF polymer sequence against UniProt P01012;
- maps publication, UniProt and PDB residue numbering without assuming a
  constant offset;
- verifies whether each literature-reported site has resolved coordinates;
- calculates minimum heavy-atom distance from each annotated site to the
  crystallographic N-acetylglucosamine (NAG);
- calculates pairwise Cα distances, chain-centroid radius and mean B-factor;
- produces a structure report, SVG chart and PyMOL visualization script.

Current chain-A mapping includes:

| Publication site | UniProt site | PDB site | Distance to NAG |
|---|---|---|---:|
| Asn-292 | N293 | ASN298 | 1.47 Å |
| Arg-84 | R85 | ARG97 | 26.99 Å |
| Lys-92 | K93 | LYS105 | 29.84 Å |
| Lys-206 | K207 | LYS216 | 25.54 Å |
| Lys-263 | K264 | LYS271 | 29.85 Å |
| Lys-322 | K323 | LYS328 | 9.21 Å |
| Arg-381 | R382 | ARG387 | 19.05 Å |

## Run

No third-party Python packages are required.

```bash
PYTHONPATH=src python -m ova_analysis.cli \
  --fasta sequences/ova_uniprot_P01012.fasta \
  --evidence data/public_literature_values.csv \
  --sites data/public_site_annotations.csv \
  --epitopes data/public_iedb_epitopes.csv \
  --structure-pdb structures/1OVA.pdb \
  --structure-cif structures/1OVA.cif \
  --out-dir outputs
```

Generated locally:

- `outputs/report.md`
- `outputs/protein_summary.csv`
- `outputs/literature_evidence.csv`
- `outputs/validated_site_annotations.csv`
- `outputs/validated_iedb_epitopes.csv`
- `outputs/iedb_epitope_map.svg`
- `outputs/epitope_report.md`
- `outputs/hwang_2014_relative_response.svg`
- `outputs/wang_2013_glycation_extent.svg`
- `outputs/structure_report.md`
- `outputs/structure_site_mapping.csv`
- `outputs/structure_site_pairwise_distances.csv`
- `outputs/1ova_site_distance_to_nag.svg`
- `scripts/visualize_1ova_sites.pml`

Run tests:

```bash
python -m pip install -e ".[test]"
python -m unittest discover -s tests -v
```

The core `ova-analysis` workflow has no third-party runtime dependencies. The
test extra also installs the scientific Python packages used by the separate
simulated-experiment demonstration.

## Important evidence limitation

The cited publications do not provide individual replicate measurements in
their abstracts. This repository therefore uses the values strictly as
**publication-reported summaries**, not as raw observations. It does not invent
sample sizes, standard deviations, error bars or p-values.

Protein properties are sequence-based estimates and do not model
post-translational modifications. Motif scanning does not confirm occupancy.
The validation step maps publication residue numbering to the UniProt sequence
and records any 0, +1 or −1 numbering offset explicitly.
The project does not establish clinical allergenicity, causality, safety or
therapeutic efficacy.

IEDB positive assay rows are not effect sizes or counts of independent
biological replications. The current epitope snapshot is deliberately limited
to two canonical reference T-cell epitopes and is not an exhaustive catalog of
overlapping OVA peptides. Most supporting records use mouse models; T-cell
reactivity is not equivalent to human IgE binding or clinical food allergy.

Structure distances are geometric descriptors for one crystallographic state.
They do not demonstrate biochemical interaction or altered allergenicity.
Chain-centroid distance is not a substitute for solvent-accessibility
calculation, and B-factors are refinement-dependent.

## Roadmap

- [x] Replace synthetic examples with traceable public data
- [x] Add OVA sequence characterization
- [x] Validate literature residue numbering against UniProt P01012
- [x] Map public modification sites onto the experimental PDB 1OVA structure
- [x] Integrate a conservative reference set of experimentally observed OVA
  T-cell epitopes from IEDB
- [ ] Analyze spatial relationships between modification sites and epitopes
- [ ] Add solvent-accessible surface-area calculations
- [ ] Add evidence-quality classification
- [x] Add automated validation with GitHub Actions

## Privacy

See [PRIVACY.md](PRIVACY.md). Transient caches and designated private-data
directories are ignored to reduce accidental disclosure risk. The public,
reproducible example outputs in this repository remain version-controlled.

## License

MIT. Referenced sequences and publication data remain subject to their source
terms and attribution requirements.
