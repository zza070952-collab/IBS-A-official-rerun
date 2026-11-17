import json,subprocess,hashlib,datetime
from pathlib import Path
root=Path('/home/asset/projects/IBS_A_OFFICIAL_RERUN_20261105')
r='''p<-c("BayesSpace","DESeq2","fgsea","scran","DirichletReg","RcppDist");stopifnot(as.character(packageVersion("BayesSpace"))=="1.10.1");stopifnot(as.numeric(formals(BayesSpace::spatialCluster)$nrep)==50000);x<-lapply(p,function(q)list(package=q,version=as.character(packageVersion(q)),path=find.package(q)));cat(jsonlite::toJSON(x,auto_unbox=TRUE));'''
p=subprocess.run([str(root/'env/r_exec.sh'),'-e',r],capture_output=True,text=True,check=True)
packages=json.loads(p.stdout)
assert all(str(root/'env/R_library') in x['path'] for x in packages)
py=subprocess.run([str(root/'env/conda/bin/python'),'-c','import spatialdm,numpy,pandas,scipy,anndata,json,sys;print(json.dumps({"python":sys.version,"spatialdm":spatialdm.__version__,"numpy":numpy.__version__,"pandas":pandas.__version__,"scipy":scipy.__version__,"anndata":anndata.__version__}))'],capture_output=True,text=True,check=True)
python=json.loads(py.stdout);assert python['spatialdm']=='0.3.1'
result={'status':'PASS','recorded_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'R_packages':packages,'python_packages':python,'formal_default_nrep':50000,'software_only_no_models_run':True,
        'R_lock_sha256':hashlib.sha256((root/'env/R_FINAL_PACKAGE_LOCK.tsv').read_bytes()).hexdigest()}
(root/'manifest/SOFTWARE_GATE.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
