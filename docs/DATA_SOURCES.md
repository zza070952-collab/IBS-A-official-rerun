# Data sources and input boundaries

The source URLs below are documented input locators from the original project's scripts and frozen protocol. This document was assembled from those local records; it does not claim that every URL was fetched again at repository-publication time.

| Input | Public source | Original expected location | Identity/checkpoint |
|---|---|---|---|
| IBS GWAS | GWAS Catalog accession **GCST90016564**, GRCh37 | `data/raw/gwas/GCST90016564_buildGRCh37.tsv` | Original full summary-statistic SHA256 recorded below. |
| Patient RNA-seq counts/RPKM | GEO **GSE166869**, `GSE166869_merged_gene_count_rpkm_by_gene.tsv.gz` | `data/raw/bulk/` | Use the raw integer count fields; the RPKM fields are not DESeq2 counts. |
| Patient sample metadata | GEO **GSE166869**, family SOFT | `data/raw/bulk/GSE166869_family.soft.gz` | Rebuild the 69-GSM/71-library mapping; preserve all exclusions and identity ambiguities. |
| Small-intestinal spatial data | Harnik author's **Data.zip**, Zenodo record **10715015** | `data/raw/spatial/Data.zip` | The frozen archive has MD5 `a855987d31f5d0f1c082a55dc4763c61`; full SHA256 below. |
| Spatial preprocessing reference source | [`yotamharnik/ST_Human_SI`](https://github.com/yotamharnik/ST_Human_SI) | `data/reference/ST_Human_SI` | Exact pinned commit `cbccfecca66b7896e0595905ad3c01f4e6484ecf`; individual source-file hashes are audited. |
| MAGMA software/reference | [Official MAGMA distribution](https://cncr.nl/research/magma/) | `env/tools/magma/`, `data/reference/gwas/` | MAGMA 1.10; GRCh37 gene locations and EUR 1000 Genomes Phase 3 reference. Acquire official distributions and respect their terms. |
| Gene identity | Official HGNC complete set | `data/reference/HGNC/hgnc_complete_set.txt` | Use the frozen source-verification receipt/hash; do not silently replace gene identities with a newer annotation. |
| RNA-seq gene sets | Hallmark and Reactome sets extracted by the recorded reference-preparation code | `data/reference/bulk/` | Preserve the actual package/gene-set versions and set membership hashes. |

## Original raw-file locators and hashes

These SHA256 values are preserved from the official-rerun bootstrap input records, not newly inferred from file names:

| File | SHA256 |
|---|---|
| `GCST90016564_buildGRCh37.tsv` | `06537a1b466084bd68f90b680055b0f1ed44b7b1edb2476d2736cb815516a8c9` |
| `GSE166869_merged_gene_count_rpkm_by_gene.tsv.gz` | `7ef6f5ac039e8ec6113ea398a80b95112f6a687d1d5dbf8a35547c50d194e44d` |
| `GSE166869_family.soft.gz` | `8ad4a1361a607fc12e88c5ab61293c4ceff958fd64596eb8bd0c47e1cf7b84a7` |
| Harnik `Data.zip` | `ff6730079b61485e4810eb3cbf78b26af76e2c912f99fc78150ea3af001798aa` |

Original locators:

- [GWAS Catalog summary file](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016564/GCST90016564_buildGRCh37.tsv)
- [GEO merged table](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166869/suppl/GSE166869_merged_gene_count_rpkm_by_gene.tsv.gz)
- [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166869/soft/GSE166869_family.soft.gz)
- [Harnik Data.zip](https://zenodo.org/records/10715015/files/Data.zip?download=1)

The historical bootstrap could reuse an individually SHA256-verified immutable **official raw-data cache**, then create independent copies in the new project. It did not make an old model, candidate ranking or calculated statistic an input to the new model fitting. Independent users should acquire the official source files or a legally accessible identical cache; an old analysis directory is not a biological input requirement.

## Sample and spatial scope

The RNA-seq QC is an independent, candidate-blind reanalysis. The project does not claim exact identity reproduction of the publication's 60 specimens without an authenticated specimen list. Non-rRNA library counts, detected genes, tissue/group identity and the public library mapping determine inclusion. Sensitivity thresholds are predefined, not selected to improve a target's P value. Multiple tissue samples must not be assumed to be independent patients or paired specimens unless source metadata establishes that relationship.

The eight spatial donors are P1, P5, P8, P9, P10, P11, P12 and P13. They provide surgical proximal small-intestinal reference tissue. Raw integer UMI, the background-corrected fractional matrix, array coordinates, tissue masks and manual exclusions are separate input layers. They must not be conflated, and no input layer should be replaced with a differently filtered object to force a published spot count.

The stricter spatial QC sensitivity in the preserved code compares barcode masks and descriptive gene detection/expression. It does **not** constitute a second BayesSpace or ligand–receptor model fit on the strict mask. Such unperformed analyses must remain untested.

## Files excluded from this code repository

Full raw GWAS SNP data, raw bulk/spatial matrices, 1000 Genomes LD files, native BayesSpace chain files, RDS/H5/H5AD models, analysis caches, large statistics/figures and private connection material are not needed for a source-code upload. Original third-party source and binaries are not implicitly relicensed by this repository. The source manifest and original data-reference records identify what must be acquired separately.
