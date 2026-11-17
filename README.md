# IBS A: official-source computational rerun

Core analysis code for the IBS A official-source rerun project, originally executed in `IBS_A_OFFICIAL_RERUN_20261105`. The workflow connects GWAS gene association, patient small-intestinal bulk RNA-seq and donor-specific spatial reference tissue. The analysis preserves complete statistics and distinguishes significant, non-significant, unavailable and untested evidence.

This repository is a source-code archive of the completed analysis. Uploading or inspecting it does not rerun the models. It does not contain raw matrices, full GWAS SNP data, 1000 Genomes LD files, fitted model objects, SSH credentials or the former project's computed results.

## Analysis modules

| Module | Main scripts | Role |
|---|---|---|
| GWAS and MAGMA | `scripts/gwas_phase0.py`, `scripts/gwas_run.py` | Validate the GRCh37 GWAS/reference inputs, audit variant identity and duplicates, run official MAGMA annotation and gene association, and retain all gene P values and BH-FDR values. |
| GWAS supporting checks | `scripts/gwas_effect_field_audit.py`, `scripts/hgnc_source_verify.py` | Audit effect fields and official gene identifiers. |
| Patient bulk RNA-seq | `scripts/bulk_phase0.py`, `scripts/bulk_run.R`, `scripts/bulk_sensitivity.R` | Reconstruct public sample/library identities and fixed QC cohorts; fit tissue-specific DESeq2 models; run full preranked fgsea analyses and predefined QC sensitivity analyses. |
| Author-source spatial preprocessing | `scripts/spatial_phase0.py`, `scripts/spatial_preprocess.py` | Validate the Harnik archive and pinned author source; reimplement the documented importer/background/QC logic in Python while keeping raw integer UMI and corrected fractional expression separate. |
| Spatial domains | `scripts/spatial_bayesspace.R` | Use official BayesSpace preprocessing, q tuning and clustering for each donor, with formal and independently seeded stability chains. |
| Spatial ligand–receptor association | `scripts/spatial_spatialdm.py` | Use official SpatialDM global and local permutation procedures with explicit eligibility, normalization and local-selection scopes. |
| Gene evidence integration | `scripts/integrate_continuous_evidence.py`, `scripts/rank_and_freeze_candidates.py` | Freeze continuous evidence before applying the recorded common exploratory ranking rule; keep statistical ranking separate from biological/model feasibility review. |
| Figures and audit support | `scripts/make_main_figures.py`, `scripts/plot_spatial_preprocessing_qc.py`, `scripts/audit_bayes_runtime_readonly.R`, `scripts/run_logged.py` | Build figures from completed outputs, inspect native chains and record real execution evidence. The `manuscript_figures/` directory contains the later R/grid manuscript renderers and their source-data preparation code, as a separate postprocessing layer. |

Historical raw-input continuation, coordinator, watcher and dependency-recovery helpers are not general-purpose starting commands and are not all included in this core source archive. They must not be started automatically on an existing completed project.

## Main methods and boundaries

- **MAGMA:** official MAGMA 1.10 gene-based association. Gene association does not establish a causal gene or the direction of expression regulation.
- **Bulk:** official DESeq2 and fgsea. The input is raw integer counts; tissue models and the 0.75M/1.25M library-QC sensitivity thresholds are fixed independently of candidate significance. The public sample mapping does not authenticate an otherwise unavailable published 60-sample identity list. Unknown patient pairing or technical batches are not invented.
- **Spatial preprocessing:** an explicit independent Python reimplementation of the pinned author's importer/background/QC logic. MATLAB was **not** executed. This does not reproduce the author's entire zonation workflow. The result is not tuned to the publication's spot count.
- **BayesSpace:** version 1.10.1; q tuning over q=4–8 with 2,000 iterations and burn-in 500; formal clustering 50,000 iterations/burn-in 1,000 and a separate independently seeded 50,000-iteration/burn-in 1,000 stability run. Low ARI or ESS remains a limitation, not a reason to silently change results.
- **SpatialDM:** version 0.3.1; donor-specific primary global and local analyses use 1,000 permutations. Local analysis is restricted to globally permutation-BH-significant pairs, capped at 100 per donor using global BH, global P and pair ID. Unselected pairs are untested locally. The project-specific non-log corrected CP10,000 local input is explicitly labelled and is not presented as raw integer UMI or the tutorial's raw-count default.
- **Interpretation:** the spatial specimens are surgical proximal small-intestinal reference tissue, not an IBS case–control spatial cohort. Spatial coexpression is not direct evidence of IBS-specific communication, directed signalling, a causal pain mechanism or therapeutic effect. No target, including NCAM1, is made causal by its rank.

## Inputs and environment

The raw public inputs and source checkpoints are described in [Data sources](docs/DATA_SOURCES.md). Actual package records and the portability limits are described in [Environment](docs/ENVIRONMENT.md).

These are preserved research scripts with original project-layout checks, paths, software gates and checkpoint semantics. They are **not** a tested one-command portable workflow. A GitHub clone alone is insufficient: raw inputs, licensed external software/reference resources, the appropriate environment, recorded parameters and validated gates are required. The original environment bootstrap cloned an existing local environment and therefore must not be presented as a standalone public installer.

Read [Reproducibility and command interfaces](docs/REPRODUCIBILITY.md) before attempting an independent reproduction. Do not run model commands against the completed original project. Use an isolated reproduction directory and record any path-only adaptation.

## Provenance and scope

The publication source is the **frozen delivered analysis/manuscript archives**, not an unverified current-server copy. The source manifest at [`provenance/SOURCE_CODE_MANIFEST.tsv`](provenance/SOURCE_CODE_MANIFEST.tsv) maps each file to its original server path, source archive/member, size, SHA256 and any lossless extraction. The archive identities are recorded in [`provenance/RELEASE_SOURCE_STATUS.json`](provenance/RELEASE_SOURCE_STATUS.json):

- `Supplementary_Source_Data.zip`: `ff1a390300cdb168d75b8cbf03a396e5bc3407c417e4c209234d484d1af947f9`.
- `IBS_A_OFFICIAL_RERUN_CORE_REVIEW.zip`: `5b3000722a2cc21a133579949e4cf92659c79e5c32cf704a6c2b4ad89ac9645c`.

These archives preserve the delivered computation and the 2026-11-08 manuscript-source snapshot. The single attempted SSH source check on 2026-11-17 failed before a session was established; this publication therefore makes **no claim** that the code is identical to any subsequent change currently present on the server. Stale local same-name or `.staged` files were not substituted for the archived bytes. Formal-model code, subsequent evidence exports and manuscript rendering can have different finalization dates; their layers are identified separately.

Execution evidence consists of real START/FINISH trace rows, native module traces, command lines, input/output hashes and exit statuses. File modification times do not establish process duration. Existing matching checkpoints, failed attempts and audited resumes remain distinguishable from new model fits.

No new biological experiment or computational result is generated by this repository publication. Third-party algorithms, author source and datasets remain subject to their own licenses and terms. Download their original distributions rather than assuming they are redistributed here.

## Repository integrity

Run `python tools/verify_source_manifest.py` to check preserved source byte counts and SHA256 without loading scientific packages or running models. `provenance/FILE_SHA256.tsv` lists the release files (excluding the checksum list itself).

`provenance/EXECUTED_COMMANDS.tsv` contains all 96 archived outer START/FINISH rows with original values. Its derivation record identifies omitted verbose hash/host columns; raw/reference input checksums are supplied separately. The source manifest and frozen parameter JSON records identify the actual delivered analysis snapshot.

The original scripts under `scripts/` and `manuscript_figures/` are preserved; this publication adds documentation and integrity tools only. See [Code-only release scope](docs/RELEASE_SCOPE.md).
