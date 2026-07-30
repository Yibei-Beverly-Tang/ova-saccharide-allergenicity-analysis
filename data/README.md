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
