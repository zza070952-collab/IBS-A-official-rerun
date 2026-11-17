# R-only quantitative figure builders. No export device is opened here.
# Source tables retain author/original statistics; no biological model is fit.
local({
  f <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
  assign(".jpr6_source_home", if (!is.null(f)) dirname(normalizePath(f, mustWork=FALSE)) else getwd(), envir=.GlobalEnv)
})

jpr6_source_dir <- function(ctx) {
  p <- NULL
  if (!is.null(ctx$figure6_source_dir)) p <- ctx$figure6_source_dir
  if (is.null(p) && !is.null(ctx$data$figure6_source_dir)) p <- ctx$data$figure6_source_dir
  if (is.null(p)) p <- file.path(.jpr6_source_home, "source_data")
  if (!dir.exists(p)) stop("Figure 6 source_data directory is missing: ", p)
  p
}
jpr6_read <- function(ctx, name) {
  p <- file.path(jpr6_source_dir(ctx), name)
  if (!file.exists(p)) stop("Required real source table is absent: ", p)
  read.delim(p, stringsAsFactors=FALSE, check.names=FALSE,
             na.strings=character(), quote="", comment.char="")
}
jpr6_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
jpr6_family <- function(ctx) if (!is.null(ctx$font_family)) ctx$font_family else "Arial"
jpr6_text <- function(label, x, y, size=8.2, color="#26333F", just="centre", face="plain", rot=0, family="Arial") {
  grid::textGrob(label, x=x, y=y, default.units="npc", just=just, rot=rot,
                 gp=grid::gpar(fontfamily=family, fontsize=size, col=color,
                               fontface=face, lineheight=1.08))
}
jpr6_seg <- function(x0,y0,x1,y1,col="#50606D",lwd=.55,lty=1) {
  grid::segmentsGrob(x0,y0,x1,y1,default.units="npc",gp=grid::gpar(col=col,lwd=lwd,lty=lty))
}
jpr6_rect <- function(x,y,w,h,fill="#FFFFFF",col=NA,lwd=.4) {
  grid::rectGrob(x=x,y=y,width=w,height=h,default.units="npc",gp=grid::gpar(fill=fill,col=col,lwd=lwd))
}
jpr6_points <- function(x,y,pch=21,col="#8064A2",fill=col,size=1.8,lwd=.5) {
  nx <- length(x); ny <- length(y)
  if (nx==0L || ny==0L || !(nx==ny || nx==1L || ny==1L)) stop("Point-coordinate lengths differ without a scalar coordinate")
  n <- max(nx,ny)
  x <- rep_len(x,n); y <- rep_len(y,n)
  grid::pointsGrob(x=x,y=y,default.units="npc",pch=pch,size=grid::unit(size,"mm"),gp=grid::gpar(col=col,fill=fill,lwd=lwd))
}
jpr6_panel <- function(children,x,y,w,h,name) {
  do.call(grid::grobTree,c(children,list(vp=grid::viewport(x=x,y=y,width=w,height=h,clip="off"),name=name)))
}
jpr6_tree <- function(children,name,width=180,height=180,sources=character()) {
  ans <- do.call(grid::grobTree,c(children,list(name=name)))
  attr(ans,"width_mm") <- width; attr(ans,"height_mm") <- height
  attr(ans,"source_files") <- sources
  attr(ans,"backend") <- "R/grid"
  ans
}
jpr6_q <- function(x,digits=2) {
  n <- jpr6_num(x)
  if (!is.finite(n)) return("NA")
  if (n==0) return("0")
  if (n < .001) return(formatC(n,format="e",digits=digits-1))
  formatC(n,format="fg",digits=digits,flag="#")
}
jpr6_fc_color <- function(x) {
  n <- jpr6_num(x)
  if (!is.finite(n)) return("#F4F5F6")
  v <- max(-1.5,min(1.5,n))
  edge <- if (v<0) "#24578A" else "#B76D9A"
  # Muted fills keep all 8-point native-value labels at readable contrast.
  grDevices::colorRampPalette(c("#FFFFFF",edge))(101)[1+round(60*abs(v)/1.5)]
}
jpr6_axis <- function(xleft,xright,ybottom,ytop,xlim,ylim,xticks,yticks,
                      xlabels=as.character(xticks),family="Arial",xlabel=NULL,ylabel=NULL) {
  sx <- function(v) xleft+(v-xlim[1])/diff(xlim)*(xright-xleft)
  sy <- function(v) ybottom+(v-ylim[1])/diff(ylim)*(ytop-ybottom)
  L <- list(jpr6_seg(xleft,ybottom,xright,ybottom,col="#26333F",lwd=.65),
            jpr6_seg(xleft,ybottom,xleft,ytop,col="#26333F",lwd=.65))
  for (i in seq_along(xticks)) L <- c(L,list(jpr6_seg(sx(xticks[i]),ybottom,sx(xticks[i]),ybottom-.012),jpr6_text(xlabels[i],sx(xticks[i]),ybottom-.041,size=8,family=family)))
  for (v in yticks) L <- c(L,list(jpr6_seg(xleft-.012,sy(v),xleft,sy(v)),jpr6_text(format(v,trim=TRUE),xleft-.026,sy(v),size=8,just="right",family=family)))
  if (!is.null(xlabel)) L <- c(L,list(jpr6_text(xlabel,(xleft+xright)/2,ybottom-.105,size=8.2,family=family)))
  if (!is.null(ylabel)) L <- c(L,list(jpr6_text(ylabel,xleft-.12,(ybottom+ytop)/2,size=8.2,rot=90,family=family)))
  list(grobs=L,sx=sx,sy=sy)
}

build_figure6 <- function(ctx) {
  fam <- jpr6_family(ctx)
  m <- jpr6_read(ctx,"FIGURE6A_EVIDENCE_MATRIX.tsv")
  f <- jpr6_read(ctx,"FIGURE6B_NCAM1_PATIENT_EFFECTS.tsv")
  s <- jpr6_read(ctx,"QPCR_SAMPLE_DERIVED_VALUES.tsv")
  st <- jpr6_read(ctx,"QPCR_PAIRED_STATISTICS.tsv")
  gs <- jpr6_read(ctx,"QPCR_GROUP_SUMMARY.tsv")
  genes <- c("NCAM1","PCLO","SKIC2","C4B","VWA7","RBM33")
  comparisons <- c("Duo_IBSD_vs_HV","Duo_IBSC_vs_HV","Jej_IBSD_vs_HV","Jej_IBSC_vs_HV")
  cohorts <- c("main","qc_low_750k","qc_high_1250k")
  stopifnot(identical(as.character(m$gene_symbol),genes),nrow(m)==6,nrow(f)==12,nrow(s)==12,nrow(st)==2)
  # A. Each column retains a distinct evidence type; no sum or new rank.
  L <- list(jpr6_text("A",.006,.973,size=11.5,face="bold",just="left",family=fam),
            jpr6_text("Six candidates retain distinct evidence profiles",.055,.973,size=9.3,face="bold",just="left",family=fam))
  left <- .145; right <- .995; nc <- 8; cw <- (right-left)/nc
  xc <- left+cw*(seq_len(nc)-.5)
  headers <- c("MAGMA\nBH-FDR","Duodenum\nIBS-D","Duodenum\nIBS-C","Jejunum\nIBS-D","Jejunum\nIBS-C","QC\ndirection","Spatial\nspots / donors","LR\napplicability")
  for (j in seq_len(nc)) L <- c(L,list(jpr6_text(headers[j],xc[j],.864,size=8.2,face="bold",family=fam)))
  L <- c(L,list(jpr6_text("Original\nrank / gene",.005,.864,size=8.2,just="left",family=fam)))
  yh <- .102; yc <- .765-(seq_len(6)-1)*yh
  for (i in seq_len(nrow(m))) {
    r <- m[i,,drop=FALSE]
    L <- c(L,list(jpr6_text(paste0(r$original_computational_rank,"  ",r$gene_symbol),.006,yc[i],size=8.5,just="left",family=fam)))
    for (j in seq_len(nc)) L <- c(L,list(jpr6_rect(xc[j],yc[i],cw*.955,yh*.93,fill=if(j==1)"#EDF3F8" else if(j==7)"#FAF4E8" else "#F5F6F7",col="#FFFFFF")))
    L <- c(L,list(jpr6_text(jpr6_q(r$MAGMA_BH_FDR,3),xc[1],yc[i],size=8.2,family=fam)))
    for (j in seq_along(comparisons)) {
      c <- comparisons[j];lf <- r[[paste0(c,"__main__log2FoldChange")]];q <- r[[paste0(c,"__main__padj")]]
      L <- c(L,list(jpr6_rect(xc[j+1],yc[i],cw*.955,yh*.93,fill=jpr6_fc_color(lf),col="#FFFFFF"),
                   jpr6_text(sprintf("%+.2f",jpr6_num(lf)),xc[j+1],yc[i]+.020,size=8.2,family=fam),
                   jpr6_text(paste0("q=",jpr6_q(q,2)),xc[j+1],yc[i]-.012,size=8,family=fam)))
      for (k in seq_along(cohorts)) {
        qq <- jpr6_num(r[[paste0(c,"__",cohorts[k],"__padj")]])
        flag <- as.character(r[[paste0(c,"__",cohorts[k],"__DEG_FDR05")]])
        px <- xc[j+1]+(k-2)*cw*.17
        if (!is.finite(qq)) L <- c(L,list(jpr6_text("x",px,yc[i]-.040,size=8,color="#67727C",family=fam)))
        else L <- c(L,list(jpr6_points(px,yc[i]-.040,pch=21,size=.75,col="#26333F",fill=if(tolower(flag)=="true")"#26333F" else "#FFFFFF",lwd=.4)))
      }
    }
    all_dirs <- vapply(comparisons,function(c) {
      ds <- vapply(cohorts,function(co) as.character(r[[paste0(c,"__",co,"__direction")]]),character(1))
      length(unique(ds))==1
    },logical(1))
    # Direction agreement is a descriptive equality of existing signs, not a test.
    L <- c(L,list(jpr6_text(paste0(sum(all_dirs),"/4 contrasts\nsame sign"),xc[6],yc[i],size=8,family=fam),
                 jpr6_text(paste0(r$spatial_positive_spots," spots\n",r$spatial_donors_detected,"/",r$spatial_donors," donors"),xc[7],yc[i],size=8,family=fam)))
    lrlabel <- if (r$LR_database_status=="NOT_IN_DATABASE") "Not in LR\ncatalogue" else paste0(r$LR_global_significant_donor_pairs,"/",r$LR_global_tested_donor_pairs," global\n",r$LR_local_supported_donor_pairs,"/",r$LR_local_tested_donor_pairs," local")
    L <- c(L,list(jpr6_text(lrlabel,xc[8],yc[i],size=8,family=fam)))
  }
  L <- c(L,list(jpr6_text("Patient cells: primary log2FC and q; QC marks below are main / 750k / 1250k.",.006,.097,size=8,just="left",family=fam),
               jpr6_text("Filled: FDR<0.05; open: tested, nonsignificant; x: FDR NA. Spatial denominator: 20,971 spots.",.006,.052,size=8,just="left",family=fam),
               jpr6_text("Blue < 0; mauve > 0 (log2FC). Direction agreement is not independent replication or an evidence score.",.006,.010,size=8,just="left",family=fam)))
  pa <- jpr6_panel(L,.5,.78,.98,.42,"figure6_panel_A")
  # B. Existing coefficients and lfcSE; the interval is a descriptive Wald display.
  co_col <- c(main="#168B8A",qc_low_750k="#B78B48",qc_high_1250k="#24578A")
  co_pch <- c(main=21,qc_low_750k=24,qc_high_1250k=22)
  f <- f[order(match(f$contrast,comparisons),match(f$QC_cohort,cohorts)),,drop=FALSE]
  lo <- jpr6_num(f$Wald_95CI_lower_descriptive);hi <- jpr6_num(f$Wald_95CI_upper_descriptive)
  xt <- pretty(range(c(lo,hi,0)),n=4);xr <- range(xt)
  sx <- function(z) .285+(z-xr[1])/diff(xr)*.49
  yy <- seq(.83,.235,length.out=nrow(f))
  L <- list(jpr6_text("B",.006,.975,size=11.5,face="bold",just="left",family=fam),
            jpr6_text("NCAM1 across patient comparisons",.072,.975,size=9.1,face="bold",just="left",family=fam),
            jpr6_text("Existing log2FC +/- 1.96 SE",.285,.913,size=8.2,just="left",family=fam),
            jpr6_text("BH-FDR",.895,.913,size=8.2,family=fam),
            jpr6_seg(sx(0),.19,sx(0),.865,col="#9AA3AB",lwd=.6,lty=2),
            jpr6_seg(.285,.19,.775,.19,col="#26333F",lwd=.65))
  group_labels <- c("Duo IBS-D","Duo IBS-C","Jej IBS-D","Jej IBS-C")
  for (i in seq_len(nrow(f))) {
    co <- as.character(f$QC_cohort[i]);q <- f$padj[i]
    L <- c(L,list(jpr6_seg(sx(lo[i]),yy[i],sx(hi[i]),yy[i],col=co_col[co],lwd=.9),
                 jpr6_points(sx(jpr6_num(f$log2FoldChange[i])),yy[i],pch=co_pch[co],size=1.8,col=co_col[co],fill=co_col[co]),
                 jpr6_text(jpr6_q(q,3),.895,yy[i],size=8.1,face=if(is.finite(jpr6_num(q))&&jpr6_num(q)<.05)"bold" else "plain",family=fam)))
  }
  for (g in seq_len(4)) L <- c(L,list(jpr6_text(group_labels[g],.008,mean(yy[(g-1)*3+seq_len(3)]),size=8.2,just="left",family=fam)))
  for (v in xt) L <- c(L,list(jpr6_seg(sx(v),.19,sx(v),.178),jpr6_text(format(v,trim=TRUE),sx(v),.145,size=8,family=fam)))
  for (k in seq_along(cohorts)) L <- c(L,list(jpr6_points(.10+(k-1)*.30,.082,pch=co_pch[cohorts[k]],size=1.6,col=co_col[cohorts[k]],fill=co_col[cohorts[k]]),jpr6_text(c("Main","750k","1250k")[k],.14+(k-1)*.30,.082,size=8,just="left",family=fam)))
  L <- c(L,list(jpr6_text("Overlapping QC cohorts; NA FDR is retained.",.008,.023,size=8,just="left",family=fam)))
  pb <- jpr6_panel(L,.305,.285,.59,.50,"figure6_panel_B")
  # C. One wet-lab panel only; horizontal offsets never modify observed values.
  gnames <- c("G1:siNC","G2:siNCAM1-1","G3:siNCAM1-2")
  xpos <- c(.25,.56,.86);offs <- c(-.040,-.013,.013,.040)
  sy <- function(y) .25+y/1.5*.59
  L <- list(jpr6_text("C",.006,.975,size=11.5,face="bold",just="left",family=fam),
            jpr6_text("NCAM1 RNA follow-up",.092,.975,size=9.1,face="bold",just="left",family=fam),
            jpr6_seg(.15,.25,.96,.25,col="#26333F",lwd=.65),
            jpr6_seg(.15,.25,.15,.84,col="#26333F",lwd=.65))
  for (v in c(0,.5,1,1.5)) L <- c(L,list(jpr6_seg(.137,sy(v),.15,sy(v)),jpr6_text(format(v,trim=TRUE),.12,sy(v),size=8,just="right",family=fam)))
  for (k in seq_along(gnames)) {
    rr <- s[s$Group==gnames[k],,drop=FALSE];rr <- rr[order(rr$Batch_ID),,drop=FALSE]
    stopifnot(nrow(rr)==4)
    summary <- gs[gs$Group==gnames[k],,drop=FALSE];stopifnot(nrow(summary)==1)
    avg <- jpr6_num(summary$mean_fold);sd <- jpr6_num(summary$SD_fold)
    col <- if(k==1) "#747B83" else "#8064A2"
    L <- c(L,list(jpr6_seg(xpos[k],sy(avg-sd),xpos[k],sy(avg+sd),col=col,lwd=.8),
                 jpr6_seg(xpos[k]-.055,sy(avg),xpos[k]+.055,sy(avg),col=col,lwd=1.4),
                 jpr6_points(xpos[k]+offs,sy(jpr6_num(rr$Relative_expression)),pch=c(21,22,24,23),size=1.8,col=col,fill=col),
                 jpr6_text(c("siNC","siNCAM1-1","siNCAM1-2")[k],xpos[k],.185,size=8,family=fam)))
  }
  # Values and P labels are read from the author-statistics table, never typed in.
  for (k in seq_len(2)) {
    p <- jpr6_num(st$P_raw[st$experimental_group==gnames[k+1]])
    stopifnot(length(p)==1,is.finite(p))
    y <- sy(c(1.16,1.38)[k]);xb <- xpos[k+1]
    L <- c(L,list(jpr6_seg(xpos[1],y-.018,xpos[1],y,col="#4B3C5E",lwd=.6),jpr6_seg(xpos[1],y,xb,y,col="#4B3C5E",lwd=.6),jpr6_seg(xb,y,xb,y-.018,col="#4B3C5E",lwd=.6),jpr6_text(paste0("P=",jpr6_q(p,3)),mean(c(xpos[1],xb)),y+.034,size=8,family=fam)))
  }
  L <- c(L,list(jpr6_text("Relative NCAM1 mRNA",.012,.555,size=8.2,rot=90,family=fam),
               jpr6_text("Mean +/- SD; biological n=4/group",.06,.103,size=8,just="left",family=fam),
               jpr6_text("Paired DeltaCt; siNC calibrated to 1",.06,.045,size=8,just="left",family=fam)))
  pc <- jpr6_panel(L,.81,.3685,.36,.333,"figure6_panel_C")
  note <- jpr6_text("Neuroblastoma surrogate\nTarget mRNA only\nNo pain-function endpoint",.635,.112,size=8,just="left",family=fam)
  ans <- jpr6_tree(list(jpr6_rect(.5,.5,1,1,"#FFFFFF"),pa,pb,pc,note),"Figure_6",180,180,
                   c("FIGURE6A_EVIDENCE_MATRIX.tsv","FIGURE6B_NCAM1_PATIENT_EFFECTS.tsv","QPCR_SAMPLE_DERIVED_VALUES.tsv","QPCR_GROUP_SUMMARY.tsv","QPCR_PAIRED_STATISTICS.tsv"))
  attr(ans,"panel_count") <- 3L
  attr(ans,"no_workflow_or_mechanism_panels") <- TRUE
  ans
}

build_s8 <- function(ctx) {
  fam <- jpr6_family(ctx)
  samples <- jpr6_read(ctx,"QPCR_SAMPLE_DERIVED_VALUES.tsv")
  tech <- jpr6_read(ctx,"QPCR_TECHNICAL_REPLICATE_SUMMARY.tsv")
  ctl <- jpr6_read(ctx,"QPCR_NEGATIVE_CONTROL_SUMMARY.tsv")
  groups <- c("G1:siNC","G2:siNCAM1-1","G3:siNCAM1-2")
  batches <- c("B01","B02","B03","B04")
  stopifnot(nrow(samples)==12,nrow(tech)==36,nrow(ctl)==24)
  bcols <- c("#595F66","#8064A2","#A994BF","#9E7FAD")
  bpch <- c(21,22,24,23)
  # A. Real paired DeltaCt, including original biological variation in controls.
  yt <- pretty(range(jpr6_num(samples$DeltaCt)),n=4);yr <- range(yt)
  ax <- jpr6_axis(.18,.96,.21,.85,c(1,3),yr,1:3,yt,c("siNC","siNCAM1-1","siNCAM1-2"),fam,ylabel="NCAM1 DeltaCt")
  L <- c(list(jpr6_text("A",.006,.975,size=11.5,face="bold",just="left",family=fam),jpr6_text("Batch-paired DeltaCt",.075,.975,size=9.3,face="bold",just="left",family=fam),jpr6_text("Four matched biological batches",.18,.913,size=8.2,just="left",family=fam)),ax$grobs)
  for (b in seq_along(batches)) {
    rr <- samples[samples$Batch_ID==batches[b],,drop=FALSE];rr<-rr[order(match(rr$Group,groups)),,drop=FALSE]
    xx<-ax$sx(1:3)+(b-2.5)*.010; yy<-ax$sy(jpr6_num(rr$DeltaCt))
    L<-c(L,list(grid::polylineGrob(x=xx,y=yy,id=rep(1,3),default.units="npc",gp=grid::gpar(col=bcols[b],lwd=.65)),jpr6_points(xx,yy,pch=bpch[b],size=1.8,col=bcols[b],fill=bcols[b])))
    L<-c(L,list(jpr6_points(.13+(b-1)*.23,.063,pch=bpch[b],size=1.6,col=bcols[b],fill=bcols[b]),jpr6_text(batches[b],.16+(b-1)*.23,.063,size=8,just="left",family=fam)))
  }
  pa<-jpr6_panel(L,.25,.755,.48,.46,"S8_A")
  # B. Both reference transcripts, original sample-level technical means.
  rt<-tech[tech$Target %in% c("RPLP0","HPRT1"),,drop=FALSE]
  yt<-pretty(range(jpr6_num(rt$mean_Ct)),n=4);yr<-range(yt)
  L<-list(jpr6_text("B",.006,.975,size=11.5,face="bold",just="left",family=fam),jpr6_text("Reference-gene Ct values",.075,.975,size=9.3,face="bold",just="left",family=fam),jpr6_text("Each point is one biological sample",.18,.913,size=8.2,just="left",family=fam))
  sy<-function(z).21+(z-yr[1])/diff(yr)*.59
  L<-c(L,list(jpr6_seg(.18,.21,.18,.80),jpr6_seg(.18,.21,.96,.21),jpr6_text("Mean technical Ct",.035,.53,size=8.2,rot=90,family=fam)))
  for(v in yt)L<-c(L,list(jpr6_seg(.167,sy(v),.18,sy(v)),jpr6_text(format(v,trim=TRUE),.145,sy(v),size=8,just="right",family=fam)))
  for(ti in seq_along(c("RPLP0","HPRT1"))){
    gene<-c("RPLP0","HPRT1")[ti];x0<-c(.20,.62)[ti];xs<-x0+c(.035,.145,.255)
    L<-c(L,list(jpr6_text(gene,mean(xs),.846,size=8.5,face="bold",family=fam)))
    for(g in seq_along(groups)){
      rr<-rt[rt$Target==gene&rt$Group==groups[g],,drop=FALSE];rr<-rr[order(rr$Batch_ID),,drop=FALSE]
      L<-c(L,list(jpr6_points(xs[g]+c(-.018,-.006,.006,.018),sy(jpr6_num(rr$mean_Ct)),pch=bpch,size=1.7,col=bcols,fill=bcols),jpr6_text(c("NC","si1","si2")[g],xs[g],.157,size=8,family=fam)))
    }
  }
  L<-c(L,list(jpr6_text("No claim of independently verified efficiency or stability.",.045,.057,size=8,just="left",family=fam)))
  pb<-jpr6_panel(L,.75,.755,.48,.46,"S8_B")
  # C. All 36 real sample-by-target technical spreads; not pseudo-replicate n.
  yt<-pretty(c(0,max(jpr6_num(tech$technical_range_Ct))),n=5);yr<-range(yt)
  ax<-jpr6_axis(.18,.96,.24,.84,c(1,3),yr,1:3,yt,c("NCAM1","RPLP0","HPRT1"),fam,ylabel="Technical Ct range")
  L<-c(list(jpr6_text("C",.006,.975,size=11.5,face="bold",just="left",family=fam),jpr6_text("All technical-replicate spreads",.075,.975,size=9.3,face="bold",just="left",family=fam),jpr6_text("One point = sample x target; three wells",.18,.913,size=8.2,just="left",family=fam)),ax$grobs)
  gcols<-c("#747B83","#8064A2","#B59AC7")
  for(ti in seq_along(c("NCAM1","RPLP0","HPRT1"))){
    gene<-c("NCAM1","RPLP0","HPRT1")[ti]
    for(g in seq_along(groups)){
      rr<-tech[tech$Target==gene&tech$Group==groups[g],,drop=FALSE];rr<-rr[order(rr$Batch_ID),,drop=FALSE]
      xx<-ax$sx(ti)+(g-2)*.050+c(-.012,-.004,.004,.012)
      L<-c(L,list(jpr6_points(xx,ax$sy(jpr6_num(rr$technical_range_Ct)),pch=bpch,size=1.7,col=gcols[g],fill=gcols[g])))
    }
  }
  for(g in seq_along(groups))L<-c(L,list(jpr6_points(.18+(g-1)*.25,.075,pch=21,size=1.6,col=gcols[g],fill=gcols[g]),jpr6_text(c("siNC","si1","si2")[g],.21+(g-1)*.25,.075,size=8,just="left",family=fam)))
  L<-c(L,list(jpr6_text("All wells retained; no post hoc outlier removal.",.045,.018,size=8,just="left",family=fam)))
  pc<-jpr6_panel(L,.25,.265,.48,.46,"S8_C")
  # D. Recorded status counts; never draw Undetermined as Ct=0 or Ct=40.
  L<-list(jpr6_text("D",.006,.975,size=11.5,face="bold",just="left",family=fam),jpr6_text("Recorded negative controls",.075,.975,size=9.3,face="bold",just="left",family=fam),jpr6_text("Undetermined / reaction count (not Ct)",.045,.913,size=8.2,just="left",family=fam))
  genes<-c("NCAM1","RPLP0","HPRT1");xc<-c(.45,.66,.87)
  for(j in seq_along(genes))L<-c(L,list(jpr6_text(genes[j],xc[j],.842,size=8.2,face="bold",family=fam)))
  ri<-0
  for(b in batches)for(control in c("NTC","NoRT")){
    ri<-ri+1;y<-.77-(ri-1)*.071
    L<-c(L,list(jpr6_text(paste(b,control),.045,y,size=8.1,just="left",family=fam)))
    for(j in seq_along(genes)){
      rr<-ctl[ctl$Batch_ID==b&ctl$Control_Type==control&ctl$Target==genes[j],,drop=FALSE];stopifnot(nrow(rr)==1)
      L<-c(L,list(jpr6_rect(xc[j],y,.195,.064,fill="#EDF3F0",col="#FFFFFF"),jpr6_text(paste0(rr$n_Undetermined,"/",rr$n_reactions),xc[j],y,size=8.2,family=fam)))
    }
  }
  L<-c(L,list(jpr6_text("48 recorded control reactions; no Ct imputation.\nRaw amplification/melting curves were not supplied.",.045,.095,size=8,just="left",family=fam)))
  pd<-jpr6_panel(L,.75,.265,.48,.46,"S8_D")
  ans<-jpr6_tree(list(jpr6_rect(.5,.5,1,1,"#FFFFFF"),pa,pb,pc,pd),"Supplementary_Figure_S8",180,180,
                 c("QPCR_SAMPLE_DERIVED_VALUES.tsv","QPCR_TECHNICAL_REPLICATE_SUMMARY.tsv","QPCR_NEGATIVE_CONTROL_SUMMARY.tsv"))
  attr(ans,"panel_count")<-4L
  attr(ans,"fold_data_repeated")<-FALSE
  ans
}
