# OVA–Saccharide Allergenicity Analysis

A reproducible Python project combining **ovalbumin protein characterization**
with **traceable publication-reported evidence** about glycation,
N-glycan modification and immune-response measurements.

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

## Run

No third-party Python packages are required.

```bash
PYTHONPATH=src python -m ova_analysis.cli \
  --fasta sequences/ova_uniprot_P01012.fasta \
  --evidence data/public_literature_values.csv \
  --sites data/public_site_annotations.csv \
  --out-dir outputs
```

Generated locally:

- `outputs/report.md`
- `outputs/protein_summary.csv`
- `outputs/literature_evidence.csv`
- `outputs/validated_site_annotations.csv`
- `outputs/hwang_2014_relative_response.svg`
- `outputs/wang_2013_glycation_extent.svg`

Run tests:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

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

## Privacy

See [PRIVACY.md](PRIVACY.md). Generated outputs and private-data directories are
ignored to reduce accidental disclosure risk.

## License

MIT. Referenced sequences and publication data remain subject to their source
terms and attribution requirements.
