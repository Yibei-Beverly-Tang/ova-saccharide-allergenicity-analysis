# Data Schema

`example_simulated_ova_data.csv` is generated with a fixed random seed and is
provided only to demonstrate the workflow. It contains no laboratory results.

Required columns:

| Column | Meaning |
|---|---|
| `sample_id` | Unique sample identifier |
| `branch` | Proposed reaction branch |
| `treatment` | Human-readable treatment label |
| `saccharide` | Control, Lactose, Mannotriose, or Dextran |
| `replicate` | Demonstration replicate number |
| `free_amino_groups_pct` | Simulated free-amino-group value |
| `ige_binding_pct` | Simulated IgE-binding value |
| `particle_size_nm` | Simulated particle size |
| `surface_hydrophobicity_au` | Simulated hydrophobicity indicator |

Future experimental data should be stored separately and de-identified before
analysis or publication.

## Public IEDB epitope snapshot

`public_iedb_epitopes.csv` is a dated, conservative snapshot of two canonical
linear OVA T-cell epitopes. Records were queried from the official IEDB API by
IEDB epitope ID on 2026-08-02. Positive counts include qualitative measures
whose value begins with `Positive`; reference counts use distinct IEDB
reference IDs among those positive rows.

Inclusion is intentionally narrow: the peptide must be linear, assigned to
ovalbumin / UniProt P01012, supported by positive T-cell assay records, and map
exactly and uniquely to the repository's P01012 sequence. Overlapping
peptide-scan variants are not included in this reference snapshot.

The data do not imply human IgE binding, clinical allergy, an effect size, or a
count of independent biological replications. Database counts can change after
the recorded retrieval date.

Regenerate the snapshot directly from IEDB:

```bash
python scripts/update_iedb_epitopes.py \
  --output data/public_iedb_epitopes.csv
```

The updater rejects unexpected peptide sequences, non-linear records and
records whose parent source antigen is not `Ovalbumin (UniProt:P01012)`.
