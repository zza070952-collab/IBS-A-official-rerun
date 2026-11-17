#!/usr/bin/env Rscript
# Independent native-output inspection only. Never sources the fitting program.
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(Matrix); library(rhdf5)
  library(jsonlite); library(digest)
})
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1")
argv <- commandArgs(trailingOnly=TRUE)
arg <- function(key) {i<-match(key,argv);if(is.na(i)||i==length(argv))stop("Required ",key);argv[i+1L]}
root <- normalizePath(arg("--root"),mustWork=TRUE)
donor <- arg("--donor")
donors <- c("P1","P5","P8","P9","P10","P11","P12","P13")
stopifnot(basename(root)=="IBS_A_OFFICIAL_RERUN_20261105",donor%in%donors)
out <- file.path(root,"reports/audits/bayes_runtime_QA")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
output <- file.path(out,paste0(donor,"_NATIVE_RUNTIME_QA.json"))
if(file.exists(output))stop("Preserve existing independent audit; no silent repeat: ",output)
iso <- function(x=Sys.time())format(x,"%Y-%m-%dT%H:%M:%OS6Z",tz="UTC")
begin <- iso();tick <- proc.time()["elapsed"]
sources <- list()
within <- function(value) {
  p <- if(startsWith(value,"/"))value else file.path(root,value)
  p <- normalizePath(p,mustWork=TRUE)
  if(!startsWith(p,paste0(root,"/")))stop("Source outside new root: ",p)
  p
}
relative <- function(p)substring(p,nchar(root)+2L)
register <- function(value) {
  p <- within(value);key <- relative(p)
  if(is.null(sources[[key]])) {
    start <- iso();h<-digest::digest(file=p,algo="sha256")
    sources[[key]] <<- list(path=p,sha256=h,bytes=unname(file.info(p)$size),
                           hash_read_started_at=start,hash_read_finished_at=iso())
  }
  p
}
sha <- function(p) {register(p);sources[[relative(within(p))]]$sha256}
js <- function(p)read_json(register(p),simplifyVector=TRUE)
check <- function(value,message)if(!isTRUE(value))stop(message)
atomic.json <- function(value,p) {
  temp<-paste0(p,".",Sys.getpid(),".tmp")
  write_json(value,temp,pretty=TRUE,auto_unbox=TRUE,na="null",digits=NA)
  if(!file.rename(temp,p))stop("Audit JSON rename failed")
}
record <- list(status="AUDIT_RUNNING_NOT_MODEL_COMPLETE",donor=donor,audit_started_at=begin,
               auditor_pid=Sys.getpid(),scope="independent_read_of_actual_native_outputs_not_certificate_only",
               model_fits_run=FALSE,native_labels_recalculated=FALSE,markers_or_PCA_recomputed=FALSE,
               large_HDF_parameter_arrays_fully_scanned=FALSE,
               software=list(R=R.version.string,BayesSpace=as.character(packageVersion("BayesSpace")),
                             SingleCellExperiment=as.character(packageVersion("SingleCellExperiment")),
                             rhdf5=as.character(packageVersion("rhdf5"))))
tryCatch({
  script.args<-commandArgs(trailingOnly=FALSE);script<-sub("^--file=","",script.args[startsWith(script.args,"--file=")])
  register(script)
  # The live control status legitimately changes while other donors run. Keep
  # the actual bytes used here, rather than treating a later control SHA as
  # immutable scientific evidence or choosing an archive by modification time.
  resume.snapshot<-file.path(out,paste0(donor,"_RESUME_STATUS_AT_AUDIT_START.json"))
  check(!file.exists(resume.snapshot)&&file.copy(within("manifest/CHAIN_AUDIT_RESUME_STATUS.json"),resume.snapshot,overwrite=FALSE),"Cannot preserve actual resume-control snapshot")
  resume<-js(resume.snapshot)
  record$resume_control_snapshot_of<-"manifest/CHAIN_AUDIT_RESUME_STATUS.json"
  check(!is.null(resume$archive_path),"Actual explicit resume archive path missing")
  original<-js(file.path(resume$archive_path,"ORIGINAL_NATIVE_FILE_HASHES.json"))
  directory<-file.path(root,"results/spatial/BayesSpace",donor)
  native.status.path<-file.path(directory,"BAYESSPACE_STATUS.json")
  status<-js(native.status.path)
  check(status$status=="BAYESSPACE_COMPLETE"&&status$donor==donor,"Donor is not natively complete; do not audit as complete")
  protocol.path<-"results/spatial/BayesSpace/FROZEN_MODEL_PROTOCOL.json"
  protocol<-js(protocol.path);protocol.sha<-sha(protocol.path)
  check(status$version=="1.10.1"&&as.character(packageVersion("BayesSpace"))=="1.10.1","Locked package version differs")
  check(status$d==15&&status$protocol_sha256==protocol.sha,"Native d/protocol identity mismatch")
  check(status$actual_chain_files_verified&&status$model_metadata_records_parameters&&status$is_new_fit,"Native completion flags missing")
  ip<-file.path(root,"data/processed/spatial",donor)
  prep<-js(file.path(ip,"PREPROCESS_STATUS.json"))
  check(prep$status=="PREPROCESS_COMPLETE"&&prep$donor==donor&&prep$manual_retained==0&&prep$n_retained==status$n_spots,"Raw-QC retained cohort identity mismatch")
  input.files<-file.path(ip,c("raw_integer_UMI_genes_by_spots.mtx","author_background_corrected_genes_by_spots.mtx","spots.tsv","genes.tsv"))
  input.shas<-unname(vapply(input.files,sha,""))
  check(identical(unname(status$input_sha256),input.shas),"Native status and actual input SHA differ")
  spots<-read.delim(register(input.files[3]),check.names=FALSE,stringsAsFactors=FALSE)
  genes<-read.delim(register(input.files[4]),check.names=FALSE,stringsAsFactors=FALSE)
  check(nrow(spots)==status$n_spots&&!anyDuplicated(spots$barcode),"Actual spatial spot identity differs")
  native.trace.path<-file.path(directory,"MODEL_EXECUTION_TRACE.tsv")
  trace<-read.delim(register(native.trace.path),check.names=FALSE,stringsAsFactors=FALSE,na.strings=character())
  check(all(c("module","source_command","input_hashes","start_iso8601","finish_iso8601","duration_seconds","exit_code","run_id")%in%names(trace)),"Actual native trace schema missing")
  base.expected<-list(donor=donor,input_hashes=input.shas,protocol_sha256=protocol.sha,version="1.10.1")
  results<-list();rows<-list()
  for(stage in c("formal","stability")) {
    stage.start<-iso()
    seed<-as.integer(20261105+100*(match(donor,donors)-1L)+if(stage=="stability")1L else 0L)
    check(status[[paste0("nrep_",stage)]]==50000&&status[[paste0("burn_in_",stage)]]==1000&&status[[paste0("seed_",stage)]]==seed,"Frozen native stage parameters changed")
    calls<-trace[trace$module==paste("BayesSpace",donor,stage,sep=":"),,drop=FALSE]
    success<-calls[as.character(calls$exit_code)=="0",,drop=FALSE]
    check(nrow(success)==1L,paste("Actual native API must have exactly one successful call:",donor,stage))
    check(is.finite(as.numeric(success$duration_seconds))&&as.numeric(success$duration_seconds)>=0&&nzchar(success$start_iso8601)&&nzchar(success$finish_iso8601),"Incomplete actual native API timing")
    rds.path<-file.path(directory,paste0(stage,"_50000.rds"))
    check(sha(rds.path)==status[[paste0(stage,"_rds_sha256")]],"Canonical native RDS hash differs from its own status")
    native<-readRDS(register(rds.path))
    metadata.actual<-metadata(native)$CLEAN_RERUN
    expected<-c(base.expected,list(stage=stage,seed=seed,q=as.integer(status$q),d=15L,nrep=50000L,burn.in=1000L))
    check(identical(metadata.actual,expected),"Independent actual RDS model-call metadata mismatch")
    raw<-assay(native,"counts")
    values<-if(inherits(raw,"sparseMatrix"))raw@x else as.vector(raw)
    check(all(is.finite(values))&&all(values>=0)&&all(values==floor(values)),"Actual RDS counts are not nonnegative integer raw UMI")
    check(ncol(native)==status$n_spots&&nrow(native)==nrow(genes)&&identical(colnames(native),as.character(spots$barcode))&&identical(rownames(native),as.character(genes$feature_id)),"Actual RDS count dimensions/barcodes/gene identity mismatch")
    check(ncol(reducedDim(native,"PCA"))==15,"Actual stored PCA dimension differs; PCA was not rerun")
    labels<-as.numeric(as.character(colData(native)$spatial.cluster))
    check(length(labels)==status$n_spots&&all(is.finite(labels))&&all(labels==floor(labels))&&all(labels>=1&labels<=status$q),"Actual stored cluster labels invalid")
    hdf<-register(file.path(directory,paste0(stage,"_50000_chain.h5")))
    check(normalizePath(metadata(native)$chain.h5,mustWork=TRUE)==hdf,"Actual native RDS chain pointer differs")
    listing<-rhdf5::h5ls(hdf,recursive=TRUE);listing<-listing[listing$otype=="H5I_DATASET",]
    dimensions<-lapply(listing$dim,function(x)as.integer(strsplit(gsub("[[:space:]]","",x),"x",fixed=TRUE)[[1]]));names(dimensions)<-listing$name
    expected.dimensions<-list(z=c(50000L,status$n_spots),weights=c(50000L,status$n_spots),mu=c(50000L,status$q*15L),lambda=c(50000L,225L),plogLik=c(50000L,1L))
    check(all(vapply(names(expected.dimensions),function(n)identical(as.integer(dimensions[[n]]),as.integer(expected.dimensions[[n]])),TRUE)),"Independent actual five native HDF dimensions differ")
    # Read the complete small LL vector only; do not load z/weights/mu/lambda.
    ll<-as.numeric(rhdf5::h5read(hdf,"plogLik"))
    check(length(ll)==50000L&&is.na(ll[1L])&&identical(as.integer(which(!is.finite(ll))),1L)&&all(is.finite(ll[2:50000]))&&all(is.finite(ll[1001:50000])),"Actual native likelihood has an unexpected nonfinite or incomplete state")
    certificate<-js(file.path(directory,paste0(stage,"_CHAIN_ARRAY_AUDIT.json")))
    check(certificate$status=="CHAIN_ARRAY_AUDIT_PASS"&&certificate$initial_undefined_plogLik_expected_from_official_source&&certificate$plogLik_R2_to_nrep_all_finite&&certificate$posterior_plogLik_all_finite&&certificate$posterior_diagnostic_states==49000&&identical(as.integer(certificate$posterior_diagnostic_R_indices),c(1001L,50000L))&&isFALSE(certificate$native_LL_array_modified)&&isFALSE(certificate$native_labels_recalculated),"Existing certificate disagrees with independent native inspection")
    check(identical(as.integer(certificate$full_plogLik_nonfinite_R_indices),1L)&&identical(as.integer(certificate$initial_plogLik_R_index),1L),"Certificate nonfinite window differs")
    diagnostics.path<-file.path(directory,paste0(stage,"_chain_diagnostics.tsv"))
    diagnostics<-read.delim(register(diagnostics.path),check.names=FALSE,stringsAsFactors=FALSE,na.strings=character())
    check(nrow(diagnostics)==1L&&diagnostics$donor==donor&&diagnostics$stage==stage&&diagnostics$nrep==50000&&diagnostics$burn.in==1000&&diagnostics$seed==seed,"Native existing diagnostic identity mismatch")
    prior<-list(applicable=FALSE)
    if(stage=="formal"&&donor%in%c("P1","P10")) {
      prior.hdf<-original[[relative(hdf)]]
      unvalidated<-file.path(directory,"formal_50000_UNVALIDATED.rds")
      prior.rds<-original[[relative(unvalidated)]]
      check(!is.null(prior.hdf)&&!is.null(prior.rds),"Preserved original formal hashes missing")
      check(sha(hdf)==prior.hdf$sha256&&sha(unvalidated)==prior.rds$sha256,"Original failed-audit native formal bytes changed")
      parse.time<-function(x)as.POSIXct(x,format="%Y-%m-%dT%H:%M:%OS",tz="UTC")
      check(parse.time(success$finish_iso8601)<parse.time(resume$started_at),"P1/P10 formal API call was repeated after explicit recovery began")
      prior<-list(applicable=TRUE,original_formal_HDF_unchanged=TRUE,original_UNVALIDATED_RDS_unchanged=TRUE,
                  canonical_RDS_SHA_equals_preserved_native_RDS_SHA=(sha(rds.path)==prior.rds$sha256),
                  original_HDF_sha256=prior.hdf$sha256,original_UNVALIDATED_RDS_sha256=prior.rds$sha256,
                  original_formal_API_finished_before_resume_started=TRUE)
    }
    results[[stage]]<-list(status="INDEPENDENT_NATIVE_AUDIT_PASS",read_started_at=stage.start,read_finished_at=iso(),
       actual_native_RDS_metadata=metadata.actual,native_RDS_sha256=sha(rds.path),native_HDF_sha256=sha(hdf),
       n_spots=ncol(native),n_genes=nrow(native),raw_counts_nonnegative_integer=TRUE,raw_nonzero_entries=sum(values!=0),
       cluster_NA_or_nonfinite_count=sum(!is.finite(labels)),q=status$q,d=15L,stored_PCA_dimensions=ncol(reducedDim(native,"PCA")),
       actual_HDF_dimensions_R_rhdf5=dimensions,LL_total_states=length(ll),LL_initial_R_index=1L,
       LL_initial_kind=if(is.nan(ll[1L]))"NaN" else "NA_R",LL_nonfinite_R_indices=as.list(which(!is.finite(ll))),
       LL_R2_to_50000_all_finite=TRUE,posterior_R_indices=c(1001L,50000L),posterior_finite_states=49000L,
       successful_native_API_call_count=nrow(success),all_native_stage_calls=calls,native_successful_call=success,
       original_formal_reuse=prior,native_diagnostics_not_recomputed=diagnostics)
    rows[[stage]]<-data.frame(donor=donor,stage=stage,status="INDEPENDENT_NATIVE_AUDIT_PASS",actual_iterations=50000L,
       posterior_finite_states=49000L,n_spots=status$n_spots,q=status$q,d=15L,seed=seed,
       native_successful_API_calls=nrow(success),original_formal_reused=(stage=="formal"&&donor%in%c("P1","P10")),
       independent_seed_ARI=status$independent_seed_ARI,stability_warning=status$stability_warning,
       effective_size_pseudologlik=diagnostics$effective_size_pseudologlik,
       first_half_mean=diagnostics$first_half_mean,second_half_mean=diagnostics$second_half_mean,
       posterior_mean=diagnostics$posterior_mean,posterior_sd=diagnostics$posterior_sd,stringsAsFactors=FALSE)
    rm(native,raw,values,ll,labels);gc(verbose=FALSE)
  }
  check(digest::digest(file=native.status.path,algo="sha256")==sha(native.status.path),"Native donor status changed during audit")
  check(digest::digest(file=within(protocol.path),algo="sha256")==protocol.sha,"Frozen protocol changed during audit")
  table<-do.call(rbind,rows)
  tsv<-file.path(out,paste0(donor,"_NATIVE_RUNTIME_QA.tsv"))
  if(file.exists(tsv))stop("Preserve existing per-donor TSV")
  write.table(table,tsv,sep="\t",quote=TRUE,row.names=FALSE,na="NA")
  record$status<-"INDEPENDENT_NATIVE_AUDIT_PASS";record$stages<-results
  record$independent_seed_ARI<-status$independent_seed_ARI;record$stability_warning<-status$stability_warning
  record$stability_scope<-"Native independent-seed diagnostic, not a convergence proof; all values retained without changing q"
  record$source_files<-sources;record$summary_table<-relative(tsv);record$summary_table_sha256<-digest::digest(file=tsv,algo="sha256")
  record$audit_finished_at<-iso();record$audit_elapsed_seconds<-unname(proc.time()["elapsed"]-tick)
  atomic.json(record,output)
  cat(toJSON(list(status=record$status,donor=donor,output=output,finished_at=record$audit_finished_at),auto_unbox=TRUE),"\n")
},error=function(e) {
  record$status<-"INDEPENDENT_NATIVE_AUDIT_FAILED_NEEDS_REVIEW";record$error<-conditionMessage(e)
  record$source_files<-sources;record$audit_finished_at<-iso();record$audit_elapsed_seconds<-unname(proc.time()["elapsed"]-tick)
  if(!file.exists(output))atomic.json(record,output)
  message(conditionMessage(e));quit(save="no",status=1L)
})
