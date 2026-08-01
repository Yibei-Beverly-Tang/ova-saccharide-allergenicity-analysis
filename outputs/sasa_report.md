# OVA Solvent-Accessible Surface-Area Report

## Method

SASA was calculated with a deterministic Shrake–Rupley point-sampling implementation using Bondi-style elemental van der Waals radii.

- Structure context: 1OVA chain A protein plus chain A NAG
- Solvent probe radius: 1.4 Å
- Sphere points per atom: 960
- Resolved protein residues reported: 383
- Hydrogen atoms are absent from the deposited X-ray model.

Algorithm reference: [Shrake and Rupley (1973)](https://doi.org/10.1016/0006-3495(73)90011-9). Atomic-radius reference: [Bondi (1964)](https://doi.org/10.1021/j100785a001).

## Annotated modification sites

| Reported site | P01012 | PDB | Total SASA | Side-chain SASA |
|---|---:|---|---:|---:|
| Arg-84 | 85 | ARG97 | 31.607 Å² | 31.481 Å² |
| Lys-92 | 93 | LYS105 | 126.102 Å² | 117.526 Å² |
| Lys-206 | 207 | LYS216 | 72.174 Å² | 69.658 Å² |
| Lys-263 | 264 | LYS271 | 72.924 Å² | 72.319 Å² |
| Asn-292 | 293 | ASN298 | 37.078 Å² | 23.886 Å² |
| Lys-322 | 323 | LYS328 | 66.023 Å² | 66.023 Å² |
| Arg-381 | 382 | ARG387 | 12.621 Å² | 12.621 Å² |

## IEDB epitope intervals

| IEDB epitope | P01012 interval | Residues | Total SASA | Mean per residue | Side-chain SASA |
|---|---:|---:|---:|---:|---:|
| [IEDB 58560](https://www.iedb.org/epitope/58560) | 258-265 | 8 | 455.103 Å² | 56.888 Å² | 399.059 Å² |
| [IEDB 28676](https://www.iedb.org/epitope/28676) | 324-340 | 17 | 434.151 Å² | 25.538 Å² | 325.52 Å² |

## Interpretation limits

- Values are absolute SASA for one crystal structure, not relative SASA.
- No exposed/buried threshold is assigned because reference maxima vary by method.
- Crystal packing, missing hydrogens and conformational dynamics can alter SASA.
- Accessibility does not establish chemical reactivity or immune recognition.

![Annotated-site SASA](annotated_site_sasa.svg)
