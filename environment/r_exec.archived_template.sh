#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/asset/projects/IBS_A_OFFICIAL_RERUN_20261105
export R_LIBS_USER="$ROOT/env/R_library" R_LIBS_SITE="$ROOT/env/R_empty_site"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /usr/bin/Rscript --vanilla "$@"
