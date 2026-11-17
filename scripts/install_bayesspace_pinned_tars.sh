#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/asset/projects/IBS_A_OFFICIAL_RERUN_20261105
cd "$ROOT"
env/r_exec.sh -e 'stopifnot(requireNamespace("digest")); pairs<-c("env/DirichletReg_0.7-2.tar.gz"="088b4c8cf4c87e661c741a87791f8bbbbe55175d76d4ed0a7b24e7fe76782d9b","env/RcppDist_0.1.1.tar.gz"="3dfcc6f34d781f055ecba64e21bd9bc9e2d504a8304f03b6a4acdde8783fe509","env/BayesSpace_1.10.1.tar.gz"="aee7f89dd2241852b6dea26e7d1a80ddedfe7b66827fcc68ef1b8193c2c4c1ed");for(f in names(pairs))stopifnot(digest::digest(file=f,algo="sha256")==pairs[[f]])'
env/r_exec.sh -e 'stopifnot(requireNamespace("Formula"),requireNamespace("maxLik"),requireNamespace("Rcpp"),requireNamespace("RcppArmadillo"));install.packages("env/DirichletReg_0.7-2.tar.gz",repos=NULL,type="source",lib="env/R_library");stopifnot(requireNamespace("DirichletReg"));install.packages("env/RcppDist_0.1.1.tar.gz",repos=NULL,type="source",lib="env/R_library");stopifnot(requireNamespace("RcppDist"))'
env/r_exec.sh -e 'install.packages("env/BayesSpace_1.10.1.tar.gz",repos=NULL,type="source",lib="env/R_library");stopifnot(requireNamespace("BayesSpace"),as.character(packageVersion("BayesSpace"))=="1.10.1",as.numeric(formals(BayesSpace::spatialCluster)$nrep)==50000);writeLines(capture.output(sessionInfo()),"env/R_FINAL_SESSIONINFO.txt");write.table(installed.packages(),"env/R_FINAL_PACKAGE_LOCK.tsv",sep="\t",row.names=FALSE,quote=FALSE)'
env/conda/bin/python -m pip freeze > env/PYTHON_SOFTWARE_LOCK_FINAL.txt
python3 -c 'import json,datetime;from pathlib import Path;p=Path("manifest/SOFTWARE_GATE.json");p.write_text(json.dumps({"status":"PASS","BayesSpace":"1.10.1","formal_default_nrep":50000,"R":"4.3.1","recorded_at":datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))'
echo SOFTWARE_GATE_PASS
