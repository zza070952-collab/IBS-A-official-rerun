#!/usr/bin/env Rscript
# Fresh official BayesSpace models, restricted to this clean project's inputs.
# Usage: env/r_exec.sh scripts/spatial_bayesspace.R --root ROOT [--donors P1,P5]
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(Matrix); library(BayesSpace)
  library(scran); library(BiocParallel); library(jsonlite)
})
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
register(SerialParam())
argv <- commandArgs(trailingOnly=TRUE)
arg <- function(key, default=NULL) { i <- match(key, argv); if (is.na(i)) default else argv[i+1L] }
root <- normalizePath(arg("--root"), mustWork=TRUE)
all.donors <- c("P1","P5","P8","P9","P10","P11","P12","P13")
donors <- strsplit(arg("--donors", paste(all.donors, collapse=",")), ",", fixed=TRUE)[[1]]
stopifnot(all(donors %in% all.donors))
version <- as.character(packageVersion("BayesSpace"))
if (version != "1.10.1") stop("Locked stable BayesSpace 1.10.1 required; actual ", version)
outroot <- file.path(root,"results/spatial/BayesSpace")
dir.create(outroot, recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(sessionInfo()), file.path(outroot,"sessionInfo.txt"))
sha <- function(f) digest::digest(file=f, algo="sha256")
iso <- function(x=Sys.time()) format(x, "%Y-%m-%dT%H:%M:%OS6Z", tz="UTC")
json <- function(x,p) { tmp <- paste0(p,".",Sys.getpid(),".tmp"); write_json(x,tmp,auto_unbox=TRUE,pretty=TRUE,na="null"); if(!file.rename(tmp,p)) stop("Atomic JSON rename failed") }
protocol <- list(version=version, input_counts="unaltered_integer_UMI_at_author_primary_spots",
  input_logcounts="log2(1+CP10000(background_corrected_counts)); background on original rows before clip",
  spatialPreprocess=list(log.normalize=FALSE,assay.type="logcounts",n.PCs=15,n.HVGs=2000),
  qTune=list(qs=4:8,nrep=2000,burn.in=500),
  q_selection="largest signed normalized improvement above endpoint line, interior q5:7; ties lower q; flat/no improvement choose q4 with warning",
  formal=list(nrep=50000,burn.in=1000,save.chain=TRUE),
  stability=list(nrep=50000,burn.in=1000,save.chain=TRUE),
  model="t",precision="equal",gamma=3,init.method="mclust",threads=1,
  seeds="20261105+100*zero_based_fixed_donor_index; stability+1; qTune+2; PCA+3",
  source_manual_exclusions="all_applied_before_this_model",candidate_blind=TRUE)
protocol.path <- file.path(outroot,"FROZEN_MODEL_PROTOCOL.json")
if(file.exists(protocol.path)) {
  old <- read_json(protocol.path,simplifyVector=FALSE)
  if(!identical(toJSON(old,auto_unbox=TRUE),toJSON(protocol,auto_unbox=TRUE))) stop("Frozen model protocol differs")
} else json(protocol,protocol.path)
protocol.sha <- sha(protocol.path)
trace.path <- file.path(outroot,"MODEL_EXECUTION_TRACE.tsv")
timed <- function(donor,stage,input.sha,params,expr) {
  start <- Sys.time(); elapsed <- proc.time()["elapsed"]; code <- 0L; err <- NULL
  result <- tryCatch(force(expr), error=function(e) {code <<- 1L; err <<- e; NULL})
  finish <- Sys.time()
  api <- if(stage %in% c("formal","stability")) "BayesSpace::spatialCluster" else if(stage=="findMarkers") "scran::findMarkers" else paste0("BayesSpace::",stage)
  row <- data.frame(module=paste("BayesSpace",donor,stage,sep=":"),
    source_command=paste0(api," ",toJSON(params,auto_unbox=TRUE)),
    input_hashes=input.sha,start_iso8601=iso(start),finish_iso8601=iso(finish),
    duration_seconds=unname(proc.time()["elapsed"]-elapsed),exit_code=code,
    hostname=Sys.info()[["nodename"]],software_version=version,
    run_id=paste0(donor,"_",stage,"_",format(start,"%Y%m%dT%H%M%OS6",tz="UTC"),"_",Sys.getpid()),
    is_new_fit=as.integer(stage %in% c("qTune","formal","stability")),
    reused_only_immutable_raw_inputs=TRUE,notes=if(is.null(err)) "actual_API_call_complete" else conditionMessage(err))
  # Each donor worker owns its trace, so independent workers never race an append.
  p <- file.path(outroot,donor,"MODEL_EXECUTION_TRACE.tsv")
  write.table(row,p,sep="\t",quote=TRUE,row.names=FALSE,col.names=!file.exists(p),append=file.exists(p))
  if(!is.null(err)) stop(err)
  result
}
atomic.rds <- function(x,p) {tmp <- paste0(p,".tmp"); saveRDS(x,tmp,compress=FALSE); if(!file.rename(tmp,p)) stop("Atomic RDS rename failed")}
inspect.chain <- function(sce,nrep) {
  h <- metadata(sce)$chain.h5
  if(is.null(h)||!file.exists(h)||file.info(h)$size<=0) stop("Native chain file is missing or empty")
  params <- metadata(sce)$CLEAN_RERUN
  if(is.null(params$q)||is.null(params$d)||params$nrep!=nrep) stop("Model-call metadata missing for actual chain inspection")
  listing <- rhdf5::h5ls(h,recursive=TRUE)
  listing <- listing[listing$otype=="H5I_DATASET",]
  dimensions <- lapply(listing$dim,function(x)as.integer(strsplit(gsub("[[:space:]]","",x),"x",fixed=TRUE)[[1]]))
  names(dimensions) <- listing$name
  # BayesSpace 1.10.1 clustering chains are unthinned.  rhdf5 reports R dimension
  # order; h5py sees the reversed HDF storage axes, which must not be confused.
  expected <- list(z=c(nrep,ncol(sce)),weights=c(nrep,ncol(sce)),
    mu=c(nrep,params$q*params$d),lambda=c(nrep,params$d*params$d),plogLik=c(nrep,1L))
  for(n in names(expected)) {
    if(is.null(dimensions[[n]])||!identical(as.integer(dimensions[[n]]),as.integer(expected[[n]])))
      stop("Actual native HDF dataset shape mismatch: ",n)
  }
  ll <- rhdf5::h5read(h,"plogLik")
  if(length(ll)!=nrep) stop("Native pseudo-log-likelihood chain is incomplete")
  ll <- as.numeric(ll)  # read-only vector view for inspection; never written back
  burn <- as.integer(params$burn.in)
  if(is.na(burn)||burn<1L||burn>=nrep) stop("Invalid recorded burn-in for chain inspection")
  # Official 1.10.1 src/cluster.cpp initializes plogLik[0]=NA_REAL and only
  # assigns i=1..nrep-1. Preserve that R-index1 NA, not an imputed value.
  if(!is.na(ll[1L])) stop("Native initial plogLik differs from official expected NA initialization")
  if(!all(is.finite(ll[seq.int(2L,nrep)]))) stop("Native plogLik has unexpected nonfinite values after initialization")
  posterior.indices <- seq.int(burn+1L,nrep)
  if(!all(is.finite(ll[posterior.indices]))) stop("Native posterior plogLik is nonfinite")
  list(status="CHAIN_ARRAY_AUDIT_PASS",file=h,file_bytes=unname(file.info(h)$size),
    version=version,formal_model_call_metadata=params,n_spots=ncol(sce),
    actual_saved_iterations_R_dimension_order=nrep,saved_chain_thinning=1L,
    dimensions_R_rhdf5=dimensions,expected_dimensions_R_rhdf5=expected,
    pseudo_loglikelihood_all_finite=all(is.finite(ll)),
    initial_plogLik_R_index=1L,initial_plogLik_CPP_index=0L,
    initial_plogLik_value_kind=if(is.nan(ll[1L]))"NaN" else "NA_R",
    initial_undefined_plogLik_expected_from_official_source=TRUE,
    full_plogLik_nonfinite_R_indices=which(!is.finite(ll)),
    plogLik_R2_to_nrep_all_finite=TRUE,posterior_plogLik_all_finite=TRUE,
    posterior_diagnostic_R_indices=c(burn+1L,nrep),posterior_diagnostic_states=length(posterior.indices),
    official_label_aggregation_R_indices=c(burn+1L,nrep),official_label_aggregation_states=nrep-burn,
    official_label_window_note="spatialCluster message prints burn.in but actual source uses seq(burn.in+1,nrep); LL diagnostics use the same actual window",
    native_LL_array_modified=FALSE,native_labels_recalculated=FALSE,
    audit_scope="all_five_native_dataset_shapes; expected_initial_NA_retained; R2:nrep_and_posterior_plogLik_finite; other_large_parameters_not_fully_scanned")
}
chain.valid <- function(sce,nrep) {
  tryCatch({inspect.chain(sce,nrep);TRUE},error=function(e){message("Chain audit: ",conditionMessage(e));FALSE})
}
load.checkpoint <- function(p,expected,chain=FALSE) {
  if(!file.exists(p)) return(NULL)
  x <- readRDS(p); actual <- metadata(x)$CLEAN_RERUN
  if(!identical(actual,expected)) stop("Existing NEW-project checkpoint metadata mismatch: ",p)
  if(chain && !chain.valid(x,50000L)) stop("Existing model checkpoint chain missing/incomplete: ",p)
  cat(iso(),"REUSE_VERIFIED_NEW_PROJECT_CHECKPOINT",p,"\n"); x
}
for(donor in donors) {
  od <- file.path(outroot,donor); dir.create(od,recursive=TRUE,showWarnings=FALSE)
  if(!dir.create(file.path(od,"worker.lock"),showWarnings=FALSE)) stop("Donor lock exists; do not duplicate: ",donor)
  tryCatch({
    ip <- file.path(root,"data/processed/spatial",donor)
    st <- read_json(file.path(ip,"PREPROCESS_STATUS.json"),simplifyVector=TRUE)
    if(st$status!="PREPROCESS_COMPLETE" || st$manual_retained!=0) stop("Author preprocessing gate failed")
    infiles <- file.path(ip,c("raw_integer_UMI_genes_by_spots.mtx","author_background_corrected_genes_by_spots.mtx","spots.tsv","genes.tsv"))
    hashes <- vapply(infiles,sha,""); input.sha <- paste(names(hashes),hashes,collapse=";")
    seed <- as.integer(20261105+100*(match(donor,all.donors)-1L))
    base.expected <- list(donor=donor,input_hashes=unname(hashes),protocol_sha256=protocol.sha,version=version)
    prep.expected <- c(base.expected,list(stage="spatialPreprocess",seed=seed+3L))
    pp <- file.path(od,"preprocessed.rds")
    sce <- load.checkpoint(pp,prep.expected)
    if(is.null(sce)) {
      spots <- read.delim(infiles[3],check.names=FALSE,stringsAsFactors=FALSE,na.strings=character())
      genes <- read.delim(infiles[4],check.names=FALSE,stringsAsFactors=FALSE,na.strings=character())
      raw <- as(readMM(infiles[1]),"dgCMatrix"); corrected <- as(readMM(infiles[2]),"dgCMatrix")
      stopifnot(identical(dim(raw),dim(corrected)),ncol(raw)==nrow(spots),nrow(raw)==nrow(genes),
                all(raw@x>=0),all(raw@x==floor(raw@x)),all(corrected@x>=0),all(is.finite(corrected@x)))
      rownames(raw) <- rownames(corrected) <- genes$feature_id
      colnames(raw) <- colnames(corrected) <- spots$barcode
      totals <- colSums(corrected); stopifnot(all(totals>0))
      lc <- corrected %*% Diagonal(x=10000/totals); lc@x <- log2(1+lc@x)
      rownames(spots) <- spots$barcode; rownames(genes) <- genes$feature_id
      # Include both documented array_* and the row/col fields required by 1.10.1.
      spots$row <- as.numeric(spots$array_row); spots$col <- as.numeric(spots$array_col)
      stopifnot(all(is.finite(spots$row)),all(is.finite(spots$col)),!anyDuplicated(spots[,c("row","col")]))
      sce <- SingleCellExperiment(assays=list(counts=raw,logcounts=lc,author_background_corrected=corrected),
                                  colData=S4Vectors::DataFrame(spots),rowData=S4Vectors::DataFrame(genes))
      set.seed(seed+3L)
      sce <- timed(donor,"spatialPreprocess",input.sha,protocol$spatialPreprocess,
        spatialPreprocess(sce,platform="Visium",n.PCs=15,n.HVGs=2000,log.normalize=FALSE,assay.type="logcounts"))
      metadata(sce)$CLEAN_RERUN <- prep.expected
      atomic.rds(sce,pp)
      write.table(data.frame(gene_symbol=genes$gene_symbol,feature_id=genes$feature_id,is_HVG=rowData(sce)$is.HVG),file.path(od,"HVGs.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
      pca <- reducedDim(sce,"PCA")
      write.table(data.frame(PC=seq_len(ncol(pca)),empirical_score_variance=apply(pca,2,var)),file.path(od,"PCA_score_variance.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    }
    tune.expected <- c(base.expected,list(stage="qTune",seed=seed+2L,nrep=2000L,burn.in=500L,qs=4:8))
    tp <- file.path(od,"qtune_2000.rds"); tuned <- load.checkpoint(tp,tune.expected)
    if(is.null(tuned)) {
      set.seed(seed+2L)
      tuned <- timed(donor,"qTune",input.sha,tune.expected,
        qTune(sce,qs=4:8,nrep=2000,burn.in=500,platform="Visium",d=15,model="t",precision="equal",gamma=3,init.method="mclust"))
      metadata(tuned)$CLEAN_RERUN <- tune.expected; atomic.rds(tuned,tp)
    }
    curve <- attr(tuned,"q.logliks"); stopifnot(all(is.finite(curve$loglik)))
    curve <- curve[order(curve$q),]; nll <- -curve$loglik
    if(nll[1]>nll[length(nll)] && diff(range(nll))>1e-8) {
      gain <- (nll[1]-nll)/(nll[1]-nll[length(nll)])
      curve$elbow_score <- gain-(curve$q-min(curve$q))/diff(range(curve$q))
      interior <- which(curve$q>min(curve$q)&curve$q<max(curve$q))
      q <- curve$q[interior[which.max(curve$elbow_score[interior])]]
      selection.warning <- if(max(curve$elbow_score[interior])<=0) "No convex elbow; fixed interior argmax is weakly supported" else "none"
    } else {q <- 4L; curve$elbow_score <- NA_real_; selection.warning <- "Flat/nonimproving tuning curve; predetermined smallest q4"}
    write.table(curve,file.path(od,"qTune_curve.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    pdf(file.path(od,"qPlot_official.pdf"),width=6,height=4); print(qPlot(tuned,qs=4:8)); dev.off()
    json(list(q=q,rule=protocol$q_selection,warning=selection.warning,candidate_information_used=FALSE),file.path(od,"Q_SELECTION.json"))
    fits <- list()
    for(stage in c("formal","stability")) {
      fitseed <- seed+if(stage=="stability")1L else 0L
      expected <- c(base.expected,list(stage=stage,seed=fitseed,q=as.integer(q),d=15L,nrep=50000L,burn.in=1000L))
      fp <- file.path(od,paste0(stage,"_50000.rds")); fit <- load.checkpoint(fp,expected,chain=TRUE)
      if(is.null(fit)) {
        native.path <- file.path(od,paste0(stage,"_50000_UNVALIDATED.rds"))
        validation.path <- file.path(od,paste0(stage,"_CHECKPOINT_VALIDATION.json"))
        if(file.exists(native.path)) {
          fit <- readRDS(native.path)
          if(!identical(metadata(fit)$CLEAN_RERUN,expected))
            stop("UNVALIDATED native return metadata differs; preserve without rerunning MCMC: ",native.path)
          cat(iso(),"RECOVER_UNVALIDATED_NATIVE_RETURN_METADATA_MATCH_ONLY_NOT_PASS",native.path,"\n")
        } else {
          hp <- file.path(od,paste0(stage,"_50000_chain.h5"))
          if(file.exists(hp)) stop("H5 exists without a saved native return; preserve and investigate before explicit recovery: ",hp)
          set.seed(fitseed)
          fit <- timed(donor,stage,input.sha,expected,
            spatialCluster(sce,q=q,platform="Visium",d=15,nrep=50000,burn.in=1000,
              model="t",precision="equal",gamma=3,init.method="mclust",save.chain=TRUE,chain.fname=hp))
          metadata(fit)$CLEAN_RERUN <- expected
          # Persist the completed native return BEFORE any audit/diagnostic can
          # fail. This file alone never satisfies a chain or dataset PASS gate.
          atomic.rds(fit,native.path)
          json(list(status="UNVALIDATED_NATIVE_RETURN_SAVED_NOT_PASS",native_return=native.path,
            canonical_checkpoint=fp,model_call_metadata=expected,saved_at=iso()),validation.path)
        }
        tryCatch(inspect.chain(fit,50000L),error=function(e) {
          json(list(status="CHAIN_VALIDATION_FAILED_NATIVE_RETURN_RETAINED_NOT_PASS",native_return=native.path,
            canonical_checkpoint=fp,error=conditionMessage(e),checked_at=iso()),validation.path)
          stop(e)
        })
        atomic.rds(fit,fp)
        json(list(status="CHAIN_VALIDATED_CANONICAL_WRITTEN",native_return=native.path,
          canonical_checkpoint=fp,model_call_metadata=expected,validated_at=iso(),
          downstream_BAYESSPACE_COMPLETE_still_required=TRUE),validation.path)
      }
      fits[[stage]] <- fit
      json(inspect.chain(fit,50000L),file.path(od,paste0(stage,"_CHAIN_ARRAY_AUDIT.json")))
      ll <- as.numeric(mcmcChain(fit,params="plogLik"))
      post <- ll[1001:50000]
      dg <- data.frame(donor=donor,stage=stage,seed=fitseed,nrep=length(ll),burn.in=1000,
                       effective_size_pseudologlik=as.numeric(coda::effectiveSize(coda::mcmc(post))),
                       first_half_mean=mean(post[1:24500]),second_half_mean=mean(post[24501:49000]),
                       posterior_mean=mean(post),posterior_sd=sd(post))
      write.table(dg,file.path(od,paste0(stage,"_chain_diagnostics.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
      write.table(data.frame(iteration=seq_along(ll),pseudologlik=ll),file.path(od,paste0(stage,"_pseudologlik_trace.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
      pdf(file.path(od,paste0(stage,"_chain_trace.pdf")),width=7,height=4); plot(seq_along(ll),ll,type="l",xlab="Actual iteration",ylab="Pseudo-log-likelihood",main=paste(donor,stage)); abline(v=1000,col="red"); dev.off()
    }
    formal <- fits$formal; stability <- fits$stability
    ari <- mclust::adjustedRandIndex(formal$spatial.cluster,stability$spatial.cluster)
    domains <- as.data.frame(colData(formal)); domains$donor <- donor
    write.table(domains,file.path(od,"domains.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    write.table(as.data.frame(table(domain=formal$spatial.cluster)),file.path(od,"domain_sizes.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    marker.rds <- file.path(od,"markers_scran_findMarkers.rds")
    if(!file.exists(marker.rds)) {
      marker <- timed(donor,"findMarkers",input.sha,list(test.type="t",direction="up",pval.type="any",assay.type="logcounts"),
        findMarkers(formal,groups=factor(formal$spatial.cluster),assay.type="logcounts",direction="up",pval.type="any",BPPARAM=SerialParam()))
      atomic.rds(marker,marker.rds)
    } else marker <- readRDS(marker.rds)
    # Preserve unequal per-comparison logFC columns rather than losing any tests.
    marker.frames <- lapply(names(marker),function(k) {z<-as.data.frame(marker[[k]]);z$feature_id<-rownames(z);z$domain<-k;z$gene_symbol<-rowData(formal)$gene_symbol[match(z$feature_id,rownames(formal))];z})
    fields <- unique(unlist(lapply(marker.frames,names)))
    marker.frames <- lapply(marker.frames,function(z){for(n in setdiff(fields,names(z)))z[[n]]<-NA;z[,fields,drop=FALSE]})
    write.table(do.call(rbind,marker.frames),file.path(od,"markers_all.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    normalized <- assay(formal,"logcounts")
    means <- lapply(sort(unique(formal$spatial.cluster)),function(k) {
      ii<-which(formal$spatial.cluster==k)
      data.frame(donor=donor,domain=k,feature_id=rownames(formal),gene_symbol=rowData(formal)$gene_symbol,
        n_spots=length(ii),n_positive_spots=Matrix::rowSums(assay(formal,"author_background_corrected")[,ii,drop=FALSE]>0),
        mean_log2_CP10000=Matrix::rowMeans(normalized[,ii,drop=FALSE]))})
    write.table(do.call(rbind,means),file.path(od,"domain_gene_expression.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    json(list(status="BAYESSPACE_COMPLETE",donor=donor,n_spots=ncol(formal),q=q,n_actual_domains=length(unique(formal$spatial.cluster)),d=15,nrep_formal=50000,
      burn_in_formal=1000,nrep_stability=50000,burn_in_stability=1000,seed_formal=seed,seed_stability=seed+1L,
      seed_qTune=seed+2L,qTune_nrep=2000,qTune_burn_in=500,version=version,
      independent_seed_ARI=ari,stability_warning=if(ari<0.6)"LOW_STABILITY_EXPLORATORY_DOMAINS" else "none",
      independent_seed_stability_cohort="same_primary_author_QC_spots_as_formal",
      strict_z3_5_QC_mask_BayesSpace_sensitivity="NOT_PERFORMED; descriptive_expression_mask_tables_are_not_model_robustness",
      q_selection_warning=selection.warning,input_sha256=unname(hashes),protocol_sha256=protocol.sha,
      formal_rds_sha256=sha(file.path(od,"formal_50000.rds")),stability_rds_sha256=sha(file.path(od,"stability_50000.rds")),
      actual_chain_files_verified=TRUE,model_metadata_records_parameters=TRUE,is_new_fit=TRUE,finished_at=iso()),file.path(od,"BAYESSPACE_STATUS.json"))
    cat(iso(),"BAYESSPACE_COMPLETE",donor,"spots",ncol(formal),"q",q,"ARI",ari,"\n")
  },finally={unlink(file.path(od,"worker.lock"),recursive=TRUE)})
}
# Canonical final status only after every donor's actual model/chain status exists.
status.paths <- file.path(outroot,all.donors,"BAYESSPACE_STATUS.json")
if(all(file.exists(status.paths))) {
  ss <- lapply(status.paths,read_json,simplifyVector=TRUE)
  if(all(vapply(ss,function(x)identical(x$status,"BAYESSPACE_COMPLETE")&&x$nrep_formal==50000&&x$nrep_stability==50000,TRUE)))
    json(list(status="BAYESSPACE_COMPLETE",n_donors=8,donors=ss,finished_at=iso()),file.path(outroot,"BAYESSPACE_ALL_DONORS_STATUS.json"))
}
