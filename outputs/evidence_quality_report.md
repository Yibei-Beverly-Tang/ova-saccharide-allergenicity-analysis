# Evidence Quality and Analytical-Scope Report

> Classes describe provenance, granularity and reanalysis limits. They are not numeric scores and do not rank biological truth.

## Classification summary

| Quality class | Records | Meaning |
|---|---:|---|
| limited_abstract_summary | 7 | Traceable numerical summary without raw observations |
| limited_abstract_annotation | 7 | Traceable site annotation without complete site-level measurements |
| curated_experimental_database_aggregate | 2 | Curated experimental database counts with heterogeneous, non-independent rows |

## Record-level classification

| Record | Source | Class | Granularity | Reanalysis readiness |
|---|---|---|---|---|
| literature_value:Wang_2013:sequence_coverage:Native OVA | DOI:10.1016/j.foodchem.2013.04.045 | limited_abstract_summary | author-reported abstract summary | descriptive use only |
| literature_value:Wang_2013:sequence_coverage:Glycated OVA | DOI:10.1016/j.foodchem.2013.04.045 | limited_abstract_summary | author-reported abstract summary | descriptive use only |
| literature_value:Wang_2013:modified_lysines:Microwave-heated OVA-glucose | DOI:10.1016/j.foodchem.2013.04.045 | limited_abstract_summary | author-reported abstract summary | descriptive use only |
| literature_value:Wang_2013:modified_lysines:Conventionally heated OVA-glucose | DOI:10.1016/j.foodchem.2013.04.045 | limited_abstract_summary | author-reported abstract summary | descriptive use only |
| literature_value:Hwang_2014:relative_IgE_production:De-N-acetylglucosaminylated OVA | DOI:10.1016/j.bbrc.2014.06.101 | limited_abstract_summary | author-reported abstract summary | descriptive use only |
| literature_value:Hwang_2014:relative_IL4_production:De-N-acetylglucosaminylated OVA | DOI:10.1016/j.bbrc.2014.06.101 | limited_abstract_summary | author-reported abstract summary | descriptive use only |
| literature_value:Hwang_2014:relative_IL5_production:De-N-acetylglucosaminylated OVA | DOI:10.1016/j.bbrc.2014.06.101 | limited_abstract_summary | author-reported abstract summary | descriptive use only |
| site_annotation:Ito_2007:Asn292 | DOI:10.1093/glycob/cwl077 | limited_abstract_annotation | author-reported abstract annotation | mapping and descriptive use only |
| site_annotation:Mao_2023:Arg84 | DOI:10.1016/j.ijbiomac.2023.123640 | limited_abstract_annotation | author-reported abstract annotation | mapping and descriptive use only |
| site_annotation:Mao_2023:Lys92 | DOI:10.1016/j.ijbiomac.2023.123640 | limited_abstract_annotation | author-reported abstract annotation | mapping and descriptive use only |
| site_annotation:Mao_2023:Lys206 | DOI:10.1016/j.ijbiomac.2023.123640 | limited_abstract_annotation | author-reported abstract annotation | mapping and descriptive use only |
| site_annotation:Mao_2023:Lys263 | DOI:10.1016/j.ijbiomac.2023.123640 | limited_abstract_annotation | author-reported abstract annotation | mapping and descriptive use only |
| site_annotation:Mao_2023:Lys322 | DOI:10.1016/j.ijbiomac.2023.123640 | limited_abstract_annotation | author-reported abstract annotation | mapping and descriptive use only |
| site_annotation:Mao_2023:Arg381 | DOI:10.1016/j.ijbiomac.2023.123640 | limited_abstract_annotation | author-reported abstract annotation | mapping and descriptive use only |
| iedb_epitope:58560 | IEDB:58560 | curated_experimental_database_aggregate | positive assay-row and distinct-reference counts | descriptive counts only |
| iedb_epitope:28676 | IEDB:28676 | curated_experimental_database_aggregate | positive assay-row and distinct-reference counts | descriptive counts only |

## Safeguards

- No p-values, effect sizes, sample sizes or uncertainty estimates are invented.
- IEDB assay rows are not treated as independent biological replicates.
- A traceable source can still have limited granularity or model relevance.
- Classification does not replace formal study-level risk-of-bias assessment.
