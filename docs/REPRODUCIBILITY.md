# Reproducibility and actual command interfaces

This is an archival guide to the preserved analysis interfaces. It is not a claim that an end-to-end reproduction was rerun during repository publication. The publication uses the verified frozen delivered computation/manuscript archives identified in `provenance/RELEASE_SOURCE_STATUS.json`. The current server could not be inspected through SSH on 2026-11-17, so later server changes are not asserted to be included. Exact historical commands, resume state and source hashes are established by the preserved execution evidence and source manifest.

Do **not** execute these commands in the completed original project. Model scripts create outputs and enforce checkpoint/gate rules. An independent reproduction needs an isolated project with the validated raw inputs, reference snapshots, software and original directory layout. Review path constraints first; changing only the variable below is insufficient for helpers with hardcoded historical paths.

## Historical working-directory conventions

The following values describe the original Linux project:

```bash
ROOT=/home/asset/projects/IBS_A_OFFICIAL_RERUN_20261105
PY="$ROOT/env/conda/bin/python"
R="$ROOT/env/r_exec.sh"
cd "$ROOT"
```

Use this only to interpret the interfaces and trace. It is not an instruction to start a new run in that directory.

## Stage interfaces

### 1. Public raw-input validation

These interfaces are present in the source:

```bash
"$PY" scripts/gwas_phase0.py --root "$ROOT"

"$PY" scripts/bulk_phase0.py \
  --root "$ROOT" \
  --merged "$ROOT/data/raw/bulk/GSE166869_merged_gene_count_rpkm_by_gene.tsv.gz" \
  --soft "$ROOT/data/raw/bulk/GSE166869_family.soft.gz"

"$PY" scripts/spatial_phase0.py --root "$ROOT"
```

The global Phase-0 gate is assembled only after all required component checks pass. Do not manufacture a PASS JSON to bypass input validation. Standalone resume helpers not present in the frozen executed source set are not substituted into this guide.

### 2. Author-source spatial preprocessing

```bash
"$PY" scripts/spatial_preprocess.py --root "$ROOT"
```

This runs the independent Python reimplementation of the documented author importer/background/QC logic on the validated raw cache. It does not execute MATLAB or all downstream author zonation analyses. Preserve raw integer counts and fractional corrected counts as different layers. The fixed stricter mask is a descriptive expression/detection sensitivity, not a new spatial-model fit.

The separate Phase-1 PASS gate must be supported by the actual all-donor preprocessing output; do not write it merely to permit the next command.

### 3. MAGMA

The model driver interface is:

```bash
"$PY" scripts/gwas_run.py \
  --root "$ROOT" \
  --phase1-gate manifest/PHASE1_GATE.json \
  --hgnc data/reference/HGNC/hgnc_complete_set.txt
```

The historical execution wrapped this driver in `run_logged.py`; its MAGMA child commands are themselves recorded by the shared logger. Inputs include the validated GRCh37 summary statistics, EUR LD reference and gene locations. Preserve official output and all-gene P/FDR; do not select candidate-associated SNPs to alter a gene result.

### 4. Bulk DESeq2, GSEA and fixed QC sensitivity

```bash
"$R" scripts/bulk_run.R --root "$ROOT" --mode prepare-reference

"$R" scripts/bulk_run.R \
  --root "$ROOT" --mode run \
  --global-gate "$ROOT/manifest/PHASE0_GATE.json"

"$R" scripts/bulk_sensitivity.R \
  --root "$ROOT" \
  --global-gate "$ROOT/manifest/PHASE0_GATE.json"
```

The main analysis uses tissue-specific raw-count DESeq2 models and all four predefined contrasts. GSEA uses the full eligible Wald-statistic ranking rather than a candidate-specific subset. Preserve native NA and its reason.

The sensitivity driver uses the same main implementation with the fixed 750k and 1250k library criteria. It explicitly links a sensitivity cohort to the already completed **new-project** primary fit if the sample set is identical; such reuse is not a new fit. A distinct cohort is fitted afresh. The code refuses to silently repeat an existing completed sensitivity branch.

### 5. BayesSpace

```bash
"$R" scripts/spatial_bayesspace.R \
  --root "$ROOT" \
  --donors P1,P5,P8,P9,P10,P11,P12,P13
```

The `--donors` R argument is comma-separated. Official BayesSpace 1.10.1 uses 2,000-iteration q tuning over q=4–8, followed by a 50,000-iteration formal chain and a separate independently seeded 50,000-iteration stability chain; both formal/stability burn-in values are 1,000. The original execution split donors between at most two single-threaded R workers. The single interface above shows supported donor selection, not a claim that all donors historically ran serially under one process.

Actual native chain lengths, retained samples, model objects, seed metadata, ARI and posterior ESS must be checked separately. A sentinel file is insufficient proof. Failed chain-inspection attempts and the audited continuation are part of the run history; do not shorten or refit already completed chains to simplify that history.

### 6. SpatialDM

```bash
"$PY" scripts/spatial_spatialdm.py \
  --root "$ROOT" \
  --donors P1 P5 P8 P9 P10 P11 P12 P13
```

The Python `--donors` argument uses space-separated donor IDs. Primary global analysis uses all eligible complete/detectable pairs, including required complex subunits, under official SpatialDM 0.3.1 and 1,000 permutations. Analytic z-score results and geometry sensitivity remain separately labelled.

Local selection is all globally permutation-BH-significant pairs, capped at 100 per donor by ascending global BH, global P and pair ID. The local selected set is not a gene-specific or pathway-specific shortlist. A pair outside it is **not tested locally**, even if globally significant. Local pair-by-spot BH inference is donor-specific; spots are not pooled across donors to create a biological replicate count.

The script's `.raw` local interface contains non-log CP10,000-normalized author-background-corrected expression, normalized using all corrected genes before selecting required LR genes. It is explicitly not raw integer UMI and not the raw-count tutorial default. Preserve this methodological distinction.

### 7. Continuous evidence, ranking and figures

```bash
"$PY" scripts/integrate_continuous_evidence.py --root "$ROOT"
```

This is the continuous-evidence stage; it does not rank or nominate candidates. All actual new source-model gates and hash freezes must already be complete.

The ranking interface is `rank_and_freeze_candidates.py --root ROOT --rule RULE --rank`. `--freeze-only` additionally needs the separate post-ranking `--review` record and does not recompute the ranking. Use the exact frozen rule and review from the documented run. Do not silently regenerate a rule, alter its thresholds or nominate NCAM1 because an earlier project used it.

`make_main_figures.py` accepts `--root`, `--freeze`, `--out`, `--formats` and `--validate-only`. It reads completed evidence. Manuscript-specific R figure renderers, if published, constitute a later postprocessing layer and must be identified as such, rather than mislabelled as the statistical model code. Existing result values, statistics and candidate ranks are not changed to make a figure fit a story.

## Timing, checkpoints and failure records

`run_logged.py` appends real START/FINISH evidence with commands, timestamps, elapsed duration, exit code and hashes. R/Python model stages also keep native call traces. File modification time is not used to infer runtime.

A matching, inspected checkpoint is reused transparently; an unverified existing file is not evidence that a fit completed. Preserve original failures and audit-specific fixes. The historical original coordinator and later recovery scripts do not form a tested general-purpose restart manager. Never launch a resume helper just because its filename contains “resume”.

An independent code review or file upload can validate syntax, source identity and provenance without refitting any model. That limited validation must not be reported as a successful new end-to-end biological analysis.

## Manuscript figure layer

The later manuscript renderer is preserved separately in `manuscript_figures/`. To reproduce that historical layer, copy its files into the independent manuscript OUTPUT directory's `code/` subdirectory and supply the original source-data tables under `OUTPUT/02_Main_Figures/source_data/`.

The preserved entry point is `env/r_exec.sh OUTPUT/code/render_manuscript_R.R OUTPUT final 1,2,3,4,5` after inspecting the installation-specific wrapper and source-data paths. The renderer sources `build_figure6_and_s8.R` even when rendering Figures 1–5, so that dependency is included. Figure 6 and S8 use later author experimental inputs, which are outside this code-only computational release; no Ct files, simulated results or images are uploaded.

The original shared ranking rule is included at `provenance/frozen_parameters/CANDIDATE_RANKING_RULE.json`. In a separate reproduction, preserve its content and original documented project location. The final biological/model-feasibility review is a separate input record, not created or silently replaced by these scripts.
