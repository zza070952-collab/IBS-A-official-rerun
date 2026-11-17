# Environment and portability

## Recorded methods

The original analysis used a project-private Python environment and an isolated R package library. The core algorithm versions enforced in the scripts are **MAGMA 1.10**, **BayesSpace 1.10.1** and **SpatialDM 0.3.1**. The archived server-exported records are included under `environment/`: `PYTHON_SOFTWARE_LOCK_FINAL.txt`, `CONDA_EXPLICIT_LOCK.txt`, `R_FINAL_PACKAGE_LOCK.tsv` and `SOFTWARE_GATE.json`. They describe the preserved executed environment; they are not a newly resolved environment specification or a check of the server's current state on 2026-11-17.

Key Python packages used by the analysis include NumPy, pandas, SciPy, AnnData, h5py, statsmodels, SpatialDM and Matplotlib. Some source-management/audit scripts also use Requests or subprocess calls to Git/curl. Native packages and executables must be installed separately when required.

Key R packages include DESeq2, fgsea, msigdbr, BayesSpace, SingleCellExperiment, SummarizedExperiment, scran, Matrix, BiocParallel, mclust, rhdf5, coda, data.table, jsonlite and digest. Compilation/system libraries for the recorded Bioconductor/R packages and HDF5 chain support must be compatible with the actual runtime.

The R final package lock contains native `installed.packages()` metadata, potentially including continuation lines and multiple library locations. Do not read it as an automatically installable two-column requirements file or select the last duplicate package row without checking the project-private library/session evidence.

## Preserved runtime layout

The historical project uses:

```text
IBS_A_OFFICIAL_RERUN_20261105/
  scripts/
  env/
    conda/bin/python
    r_exec.sh
    R_library/
    R_empty_site/
    tools/magma/magma
  data/
    raw/
    reference/
    metadata/
    processed/
  manifest/
  logs/
  results/
```

`env/r_exec.sh` selects the project's R library and invokes Rscript. Its publication copy is reconstructed literally from the archived bootstrap script's generated-wrapper heredoc and checked against the matching SHA256 in the archived file manifest. Its derivation is explicit in the source manifest. This establishes identity with that frozen archived wrapper, **not** a new live-server readback. The wrapper remains installation-specific. `run_logged.py` uses Unix process/file-lock facilities. The code is intended for Linux; it has not been validated as a native Windows workflow.

Several drivers require a root directory named `IBS_A_OFFICIAL_RERUN_20261105`; some setup/audit helpers retain original absolute installation paths. Merely setting a new `ROOT` variable does not adapt every helper. On a new machine, inspect those constraints in a separate source copy, document any path-only changes and retain the original source manifest. Do not modify the original completed project to make a reproduction run convenient.

The original `bootstrap_environment.sh` cloned an existing local Conda environment and copied an existing R library before installing/version-checking the required official packages. It is historical infrastructure, **not** a portable public installer. The preserved Conda/Python locks, R package metadata, software gate and session records describe the executed environment; their existence does not establish that the environment can be recreated unchanged on every host.

## Essential source-file relationships

- Keep `spatial_phase0.py`, `spatial_preprocess.py` and `spatial_spatialdm.py` together: the latter two import shared donor, hashing and preprocessing helpers.
- Keep `bulk_sensitivity.R` beside `bulk_run.R`: the sensitivity driver sources the main implementation.
- Preserve `run_logged.py` and the R wrapper expected by the model drivers.
- Download the original Harnik source at the pinned commit and the official MAGMA/reference distributions separately. A GitHub repository's location does not replace these data/software inputs.

For a new execution, first verify the package versions and input integrity, then preserve new command logs. Do not edit the version assertions to silence a mismatch or present a failed dependency installation as a completed model run.

The Python lock is an archival `pip freeze` record with 197 local `file://` build-source entries. It cannot be used unmodified as a portable `pip install -r` file. The Conda explicit lock records a Linux platform; it is not a cross-platform dependency solver. No environment was installed, resolved or executed as part of this upload.
