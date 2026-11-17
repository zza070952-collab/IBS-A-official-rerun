"""Read frozen spatial results and prepare manuscript source data; NO model fitting.

No spatial or gene-test package is called. Python is used for non-visual file
inspection/data conversion only; all drawing/export remains in the R backend.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,shutil,time,datetime,collections,io,struct
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

DONORS=['P1','P5','P8','P9','P10','P11','P12','P13']
GENES=['NCAM1','PCLO','SKIC2','C4B','VWA7','RBM33']
PROGRAMS={
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19','VIL1','ALPI','APOA4','SLC5A1','SI','ANPEP'],
 'Immune':['PTPRC','CD74','HLA-DRA','HLA-DPA1','LST1','CD3D','CD3E','CD79A','MS4A1'],
 'Stromal':['COL1A1','COL1A2','COL3A1','DCN','LUM','COL6A1','COL6A2','PDGFRA'],
 'Endothelial':['PECAM1','VWF','EMCN','KDR','ENG'],
 'Contractile':['ACTA2','TAGLN','MYL9','CNN1','DES'],
 'Goblet':['MUC2','SPINK4','CLCA1','AGR2','FCGBP'],
 'Paneth':['DEFA5','DEFA6','LYZ','REG3A']}
HEATMAP_MARKERS=['EPCAM','VIL1','ALPI','APOA4','DEFA5','DEFA6','MUC2','PTPRC','CD74','COL1A1','DCN','PECAM1','ACTA2','TAGLN']
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def truth(v):return str(v).lower() in ('true','1')
def fnum(v):
 try:return float(v)
 except (ValueError,TypeError):return np.nan
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--out',required=True);args=ap.parse_args()
 root=Path(args.root).resolve();out=Path(args.out).resolve()
 assert root.name=='IBS_A_OFFICIAL_RERUN_20261105'
 assert out.is_relative_to(root/'deliverables') and out.name.startswith('MANUSCRIPT_V1_')
 sources=out/'02_Main_Figures/source_data/spatial';snap=out/'05_Code_and_Inputs/source_snapshot';review=out/'04_Review'
 for p in [sources,snap,review,sources/'images']:p.mkdir(parents=True,exist_ok=True)
 assert not (sources/'SPATIAL_EXTRACTION_STATUS.json').exists(),'Existing extraction status; do not silently repeat or replace'
 start=time.monotonic();manifest={};derived=[]
 def track(rel,role,copy=False,compress=False):
  p=root/rel;assert p.is_file(),str(p)
  if rel in manifest:return p
  h=sha(p);dest='';stored_sha='';encoding='ORIGINAL_BYTES'
  if copy:
   target=snap/rel
   if compress and p.suffix not in ('.gz','.png','.jpg','.jpeg'):
    target=target.with_name(target.name+'.gz');encoding='GZIP_WRAPPED_ORIGINAL_BYTES'
   target.parent.mkdir(parents=True,exist_ok=True)
   assert not target.exists(),str(target)
   if encoding=='GZIP_WRAPPED_ORIGINAL_BYTES':
    with p.open('rb') as src,target.open('xb') as raw:
     with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=9) as dst:shutil.copyfileobj(src,dst,4*1024*1024)
   else:
    with p.open('rb') as src,target.open('xb') as dst:shutil.copyfileobj(src,dst,4*1024*1024)
   dest=str(target.relative_to(out));stored_sha=sha(target)
  manifest[rel]=dict(source_path=str(p),role=role,bytes=p.stat().st_size,sha256_before=h,sha256_after='',copied_path=dest,stored_sha256=stored_sha,restore_encoding=encoding,provenance_note='Native official rerun; no P/FDR/statistic changes or model fit')
  return p
 def tab(rel,role='native table',copy=True):
  p=track(rel,role,copy,compress=True)
  return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
 def js(rel,role='native status',copy=True):return json.loads(track(rel,role,copy).read_text())
 def write(name,df,rule,source_paths):
  p=sources/name
  assert not p.exists(),str(p)
  df.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression='gzip' if name.endswith('.gz') else None)
  derived.append(dict(output=str(p.relative_to(out)),rows=len(df),columns=len(df.columns),sha256=sha(p),filter_or_derivation=rule,sources=';'.join(source_paths)))

 aliases={g:{g} for g in GENES}
 alias_rel='deliverables/A6_FINAL_COMPARISON/Source_Evidence/A6_IDENTITY_MAPPED_ALIASES.tsv'
 at=tab(alias_rel,'existing exact identity-to-feature alias mapping')
 for _,r in at.iterrows():
  if r.canonical_symbol in aliases and r.mapping_status=='MAPPED':aliases[r.canonical_symbol].add(r.exact_mapped_value)
 assert 'SKIV2L' in aliases['SKIC2'],'Must use verified SKIV2L alias, never SKIV2L2'
 expr_rows=[];all_domains=[];detection=[];diagnostics=[];annotation=[];markers_top=[];marker_heat=[];summaries=[];lr_all=[];lr_elig=[];lr_selection=[];local_summaries=[];pair_detail=[];qc_summaries=[]
 for donor in DONORS:
  drel=f'results/spatial/BayesSpace/{donor}/domains.tsv';dom=tab(drel,'all actual domain spots with full-resolution image coordinates')
  assert not dom.barcode.duplicated().any()
  qrel=f'results/spatial/BayesSpace/{donor}/BAYESSPACE_STATUS.json';status=js(qrel)
  prep_rel=f'data/processed/spatial/{donor}/PREPROCESS_STATUS.json';prep=js(prep_rel)
  for stage in ['formal','stability']:
   rel=f'results/spatial/BayesSpace/{donor}/{stage}_chain_diagnostics.tsv';c=tab(rel,'native actual 50000-iteration chain diagnostics');diagnostics.append(c)
  image_rel=f'data/raw/spatial/extracted/{donor}_tissue_hires_image.png';scale_rel=f'data/raw/spatial/extracted/{donor}_scalefactors_json.json'
  image=track(image_rel,'original author tissue image; no cropping/recolouring',copy=False);scale=js(scale_rel,'actual full-resolution to hires scale factors')
  image_out=sources/'images'/f'{donor}_tissue_hires_image.png'
  assert not image_out.exists();shutil.copyfile(image,image_out);assert sha(image_out)==sha(image)
  dom['donor']=donor;dom['hires_x']=pd.to_numeric(dom.pixel_col)*float(scale['tissue_hires_scalef']);dom['hires_y']=pd.to_numeric(dom.pixel_row)*float(scale['tissue_hires_scalef'])
  dom['domain']=dom['spatial.cluster'];dom['image_filename']=str(image_out.relative_to(sources));all_domains.append(dom)
  with image.open('rb') as f:header=f.read(24)
  assert header[:8]==b'\x89PNG\r\n\x1a\n' and header[12:16]==b'IHDR'
  image_width,image_height=struct.unpack('>II',header[16:24])
  outside=((dom.hires_x<0)|(dom.hires_x>=image_width)|(dom.hires_y<0)|(dom.hires_y>=image_height))
  summaries.append(dict(donor=donor,n_spots=len(dom),q=status['q'],independent_seed_ARI=status['independent_seed_ARI'],tissue_hires_scalef=scale['tissue_hires_scalef'],spot_diameter_fullres=scale.get('spot_diameter_fullres','NA'),image_width_pixels=image_width,image_height_pixels=image_height,n_spots_outside_hires_image=int(outside.sum()),scale_bar_units='500 full-resolution image pixels; pixel-to-micron conversion not independently verified',n_input_UTT=prep['n_input_UTT'],n_retained=prep['n_retained'],image_filename=str(image_out.relative_to(sources)),image_sha256=sha(image_out),representative_selection='P8/P9 provisional from candidate-blind domain stability; root must inspect histology before final use'))
  qc_rel=f'data/processed/spatial/{donor}/UTT_QC_all_spots.tsv';qc=tab(qc_rel,'all raw/retained/QC-excluded spot records');qc['donor']=donor
  kept=qc[qc.entered_primary_analysis.map(truth)]
  qc_summaries.append(dict(donor=donor,UTT_spots=len(qc),retained_spots=len(kept),median_raw_UMI=pd.to_numeric(kept.raw_UMI).median(),median_raw_genes=pd.to_numeric(kept.raw_nGenes).median(),median_corrected_total=pd.to_numeric(kept.background_corrected_total).median(),median_corrected_genes=pd.to_numeric(kept.corrected_nGenes).median(),median_corrected_mito_fraction=pd.to_numeric(kept.corrected_mito_fraction).median(),exclusion_counts_json=json.dumps(qc.loc[~qc.entered_primary_analysis.map(truth),'exclusion_reason'].value_counts().to_dict()),source_path=qc_rel))
  detrel=f'data/processed/spatial/{donor}/gene_detection.tsv';dt=tab(detrel,'all native corrected and raw donor gene detection statistics',copy=False)
  dt=dt[dt.gene_symbol.isin(set().union(*aliases.values()))].copy()
  for gene in GENES:
   hit=dt[dt.gene_symbol.isin(aliases[gene])].copy();assert len(hit)==1,(donor,gene,len(hit));hit['canonical_symbol']=gene;detection.append(hit)

  hrel=f'data/processed/spatial/{donor}/author_reimplementation_primary.h5ad';hp=track(hrel,'existing expression object read solely for six-gene plotting values; object not copied')
  obj=ad.read_h5ad(hp)
  assert len(obj)==len(dom) and obj.obs_names.is_unique
  assert 'counts_author_background_corrected' in obj.layers
  layer=obj.layers['counts_author_background_corrected'];total=np.asarray(layer.sum(axis=1)).ravel()
  assert np.all(total>0)
  symbols=np.asarray(obj.var['gene_symbol'].astype(str))
  obs=pd.DataFrame({'barcode':obj.obs['barcode'].astype(str).to_numpy()})
  obs=obs.merge(dom[['barcode','domain','pixel_col','pixel_row','hires_x','hires_y']],on='barcode',how='left',validate='one_to_one');assert obs.domain.notna().all()
  for gene in GENES:
   ids=np.flatnonzero(np.isin(symbols,list(aliases[gene])));assert len(ids)>0,(donor,gene)
   vals=np.asarray(layer[:,ids].sum(axis=1)).ravel();log=np.log1p(10000*vals/total)
   frame=obs.copy();frame['donor']=donor;frame['gene_symbol']=gene;frame['raw_symbols']=';'.join(sorted(set(symbols[ids])));frame['n_original_feature_rows']=len(ids);frame['background_corrected_counts']=vals;frame['all_gene_corrected_total']=total;frame['log1p_CP10000']=log;frame['corrected_positive']=vals>0
   expr_rows.append(frame)
   expected=dt[dt.gene_symbol.isin(aliases[gene])].iloc[0]
   assert int((vals>0).sum())==int(expected.n_positive_spots),(donor,gene,'detection discrepancy')
   assert np.isclose(float(log.mean()),float(expected.mean_log_expression),rtol=1e-10,atol=1e-12),(donor,gene,'mean-expression discrepancy')
  del obj,layer
  mrel=f'results/spatial/BayesSpace/{donor}/markers_all.tsv';mk=tab(mrel,'all domain marker statistics including non-significant and NA')
  erel=f'results/spatial/BayesSpace/{donor}/domain_gene_expression.tsv';exp=tab(erel,'all existing domain gene means, log2(1+CP10000) distinct from natural log detection',copy=True)
  exp=exp[exp.gene_symbol.isin(set().union(*aliases.values()))].copy();exp['donor']=donor
  write(f'{donor}_A6_DOMAIN_EXPRESSION_NATIVE.tsv.gz',exp,'Exact mapped six-gene rows from every existing domain; no test or changed values',[erel])
  for domain in sorted(dom.domain.unique(),key=int):
   sub=mk[mk.domain.astype(str)==str(domain)].copy();sub['_fdr']=pd.to_numeric(sub.FDR,errors='coerce');sub['_fc']=pd.to_numeric(sub['summary.logFC'],errors='coerce')
   top=sub.sort_values('_fdr',kind='stable').head(15).drop(columns=['_fdr','_fc']).copy();top['donor']=donor;markers_top.append(top)
   for gene in HEATMAP_MARKERS:
    hits=sub[sub.gene_symbol==gene]
    if len(hits)==1:
     r=hits.iloc[0];marker_heat.append(dict(donor=donor,domain=domain,marker=gene,native_summary_logFC=r['summary.logFC'],native_FDR=r['FDR'],source_path=mrel,source_feature_id=r.get('feature_id','')))
    else:marker_heat.append(dict(donor=donor,domain=domain,marker=gene,native_summary_logFC='NA',native_FDR='NA',source_path=mrel,source_feature_id='MISSING_OR_MULTIPLE_NATIVE_ROWS'))
   positive=sub[(sub._fdr<0.05)&(sub._fc>0)]
   support={program:sorted(set(positive.gene_symbol)&set(genes)) for program,genes in PROGRAMS.items()}
   enriched={p:g for p,g in support.items() if len(g)>=2}
   # Descriptive post hoc screening only, with histology explicitly pending.
   if len(enriched)==1:label=next(iter(enriched))+'-marker-enriched candidate'
   elif len(enriched)>1:label='Mixed marker programs: '+' / '.join(enriched)
   else:label='Domain '+str(domain)+' / unresolved'
   annotation.append(dict(donor=donor,domain=domain,supporting_markers=json.dumps(enriched,ensure_ascii=False),countermarkers=json.dumps({k:v for k,v in support.items() if v and k not in enriched},ensure_ascii=False),histology_support='NOT_YET_REVIEWED',source_reference='Harnik 2024 atlas DOI 10.1038/s41586-024-07793-3; native marker statistics',proposed_label=label,confidence='TENTATIVE_LOW_UNTIL_HISTOLOGY_REVIEW',limitations='Spots contain mixtures; >=2 native BH<0.05 positive min-logFC markers is a descriptive screening convention, not new enrichment testing or cell-type assignment; histology still pending',marker_source_path=mrel,marker_source_sha256=manifest[mrel]['sha256_before'],selection_rule='Fixed manuscript-stage post hoc marker list covering epithelial/immune/stromal/endothelial/contractile/goblet/Paneth programs; use native FDR<0.05 and summary.logFC>0; not an originally preregistered annotation rule; no NCAM1-driven labels'))
  base=f'results/spatial/SpatialDM/{donor}'
  glob=tab(base+'/global_LR.tsv','all tested global ligand-receptor statistics, native permutation P/BH unchanged');glob['donor']=donor;lr_all.append(glob)
  elig=tab(base+'/ALL_LR_ELIGIBILITY.tsv','all database hypotheses and eligibility, including untested');elig['donor']=donor;lr_elig.append(elig)
  sel=tab(base+'/LOCAL_SELECTION_FROZEN.tsv','all frozen global-to-local selection statuses, cap exclusions remain untested');sel['donor']=donor;lr_selection.append(sel)
  ls=js(base+'/SPATIALDM_STATUS.json');summaries[-1].update(global_tested=ls['n_global_pairs_tested'],global_significant=ls['n_global_permutation_BH_significant'],local_tested=ls['n_local_pairs'],local_significant_pair_spot_cells=ls['n_local_pair_spot_cells_significant'])
  for suffix in ['/GLOBAL_CHECKPOINT.json','/LOCAL_CHECKPOINT.json','/primary_global_MODEL_ARRAY_AUDIT.json','/local_MODEL_ARRAY_AUDIT.json','/MODEL_EXECUTION_TRACE.tsv']:
   track(base+suffix,'actual model runtime/checkpoint/array evidence, read-only',copy=True,compress=suffix.endswith('.tsv'))
  lrel=base+'/local_spot_statistics.tsv.gz';lp=track(lrel,'complete actual local pair-by-spot statistics, all tested negative and positive cells retained',copy=True)
  per_pair=collections.defaultdict(lambda:{'n_spots_tested':0,'n_BH_significant_spots':0,'min_native_padj':np.inf})
  with gzip.open(lp,'rt',newline='') as fh:
   for r in csv.DictReader(fh,delimiter='\t'):
    s=per_pair[r['pair_id']];s['n_spots_tested']+=1;s['n_BH_significant_spots']+=truth(r['significant']);p=fnum(r['permutation_padj']);s['min_native_padj']=min(s['min_native_padj'],p) if np.isfinite(p) else s['min_native_padj']
  assert len(per_pair)==ls['n_local_pairs']
  for pair,s in per_pair.items():local_summaries.append(dict(donor=donor,pair_id=pair,**s,source_path=lrel))

 domains=pd.concat(all_domains,ignore_index=True);expressions=pd.concat(expr_rows,ignore_index=True)
 write('ALL_DOMAIN_SPOTS.tsv.gz',domains,'All actual frozen primary spots; image x=pixel_col*actual_scale, y=pixel_row*actual_scale; no domain fitting',[f'results/spatial/BayesSpace/{d}/domains.tsv' for d in DONORS])
 write('A6_SPOT_EXPRESSION.tsv.gz',expressions,'Six mapped genes; sum original corrected counts after background clipping then log1p(10000*gene/sum_all_original_gene_corrected_counts); zeros retained; detection/mean equal native gene_detection',[f'data/processed/spatial/{d}/author_reimplementation_primary.h5ad' for d in DONORS])
 write('A6_DONOR_DETECTION_NATIVE.tsv',pd.concat(detection,ignore_index=True),'Exact native gene-detection records; full precision original values',[f'data/processed/spatial/{d}/gene_detection.tsv' for d in DONORS])
 write('DONOR_SUMMARY.tsv',pd.DataFrame(summaries),'Native counts/status plus actual image scale and source hash',[f'results/spatial/SpatialDM/{d}/SPATIALDM_STATUS.json' for d in DONORS])
 write('DONOR_QC_SUMMARY.tsv',pd.DataFrame(qc_summaries),'Descriptive medians of retained existing UTT_QC rows; not tests or new filtering',[f'data/processed/spatial/{d}/UTT_QC_all_spots.tsv' for d in DONORS])
 write('DOMAIN_TOP15_MARKERS_NATIVE.tsv.gz',pd.concat(markers_top,ignore_index=True),'Each donor/domain existing FDR ascending, stable original row-order ties, top15 for readability; full original marker tables copied',[f'results/spatial/BayesSpace/{d}/markers_all.tsv' for d in DONORS])
 write('FIXED_MARKER_NATIVE_EFFECTS.tsv',pd.DataFrame(marker_heat),'Fixed basic-marker panel across all domains; preserve existing minimum pairwise log2 FC and FDR; no test recalculation',[f'results/spatial/BayesSpace/{d}/markers_all.tsv' for d in DONORS])
 write('ALL_CHAIN_DIAGNOSTICS_NATIVE.tsv',pd.concat(diagnostics,ignore_index=True),'All 16 native chain rows; no convergence threshold invented',[f'results/spatial/BayesSpace/{d}/formal_chain_diagnostics.tsv' for d in DONORS]+[f'results/spatial/BayesSpace/{d}/stability_chain_diagnostics.tsv' for d in DONORS])
 write('ALL_GLOBAL_LR_NATIVE.tsv.gz',pd.concat(lr_all,ignore_index=True),'All tested global LR records; no new BH',[f'results/spatial/SpatialDM/{d}/global_LR.tsv' for d in DONORS])
 write('ALL_LR_ELIGIBILITY_NATIVE.tsv.gz',pd.concat(lr_elig,ignore_index=True),'All database eligibility rows, not only tested/positive',[f'results/spatial/SpatialDM/{d}/ALL_LR_ELIGIBILITY.tsv' for d in DONORS])
 write('ALL_LOCAL_SELECTION_NATIVE.tsv.gz',pd.concat(lr_selection,ignore_index=True),'All native selection statuses; untested cap-excluded records retained',[f'results/spatial/SpatialDM/{d}/LOCAL_SELECTION_FROZEN.tsv' for d in DONORS])
 local=pd.DataFrame(local_summaries);write('ALL_TESTED_LOCAL_PAIR_SUMMARY.tsv',local,'Count all actually tested local spot rows and native significant flags, minimum original BH; not a new test',[f'results/spatial/SpatialDM/{d}/local_spot_statistics.tsv.gz' for d in DONORS])
 g=pd.concat(lr_all,ignore_index=True);e=pd.concat(lr_elig,ignore_index=True)
 # Rank only existing LR summaries for a target-blind visual selection.
 rank=[]
 for pair,rr in g.groupby('pair_id',sort=False):
  sig=rr[rr.significant.map(truth)]
  effects=pd.to_numeric(rr.global_I,errors='coerce')
  rank.append(dict(pair_id=pair,n_donors_tested=len(rr),n_donors_significant=len(sig),median_existing_global_I=effects.median(),median_existing_padj=pd.to_numeric(rr.padj,errors='coerce').median()))
 ranking=pd.DataFrame(rank).sort_values(['n_donors_significant','median_existing_global_I','pair_id'],ascending=[False,False,True],kind='stable')
 write('GLOBAL_PAIR_PLOT_RANKING.tsv',ranking,'All native tested pair identities ranked by significant-donor count then median existing global_I then lexical pairID; descriptive plot order only',['ALL_GLOBAL_LR_NATIVE.tsv.gz'])
 selected=ranking.head(10).pair_id.tolist();matrix=[]
 for pair in selected:
  for donor in DONORS:
   rr=g[(g.donor==donor)&(g.pair_id==pair)];ee=e[(e.donor==donor)&(e.pair_id==pair)]
   matrix.append(dict(donor=donor,pair_id=pair,test_status='TESTED' if len(rr) else 'UNTESTED',pvalue=rr.iloc[0].pvalue if len(rr) else 'NA',padj=rr.iloc[0].padj if len(rr) else 'NA',global_I=rr.iloc[0].global_I if len(rr) else 'NA',significant=rr.iloc[0].significant if len(rr) else 'NA',eligibility_status=ee.iloc[0].test_status if len(ee) else 'NOT_IN_DATABASE'))
 write('GLOBAL_TOP10_MATRIX.tsv',pd.DataFrame(matrix),'Top10 from target-blind existing summary ranking, all eight donor test statuses retained',['GLOBAL_PAIR_PLOT_RANKING.tsv','ALL_GLOBAL_LR_NATIVE.tsv.gz','ALL_LR_ELIGIBILITY_NATIVE.tsv.gz'])
 # Local display chosen without candidate identity from strongest actually local-supported pair.
 ranked_local=local[local.n_BH_significant_spots>0].sort_values(['n_BH_significant_spots','donor','pair_id'],ascending=[False,True,True],kind='stable')
 if len(ranked_local)==0:ranked_local=local.sort_values(['donor','pair_id'],kind='stable')
 chosen=ranked_local.iloc[0];lrel=f'results/spatial/SpatialDM/{chosen.donor}/local_spot_statistics.tsv.gz';rows=[]
 with gzip.open(root/lrel,'rt',newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['pair_id']==chosen.pair_id:rows.append(r)
 write('REPRESENTATIVE_ACTUAL_LOCAL.tsv',pd.DataFrame(rows),'Actually tested local pair with most native BH-significant spots; ties lexical donor/pair; candidate-blind, no statistical refit',[lrel])
 a6rel='deliverables/A6_FINAL_COMPARISON/Source_Evidence/A6_SPATIALDM_PAIR_DETAIL.tsv.gz';a6=tab(a6rel,'complete existing exact-component six-candidate pair joins and all NA statuses')
 write('A6_EXACT_PAIR_DETAIL_NATIVE.tsv',a6,'All native existing theoretical candidate pairs; no substring matching or different-pair attribution',[a6rel])
 annotation=pd.DataFrame(annotation);assert len(annotation)==44
 annotation.to_csv(review/'DOMAIN_ANNOTATION_EVIDENCE.tsv',sep='\t',index=False)
 annotation.to_csv(sources/'DOMAIN_ANNOTATION_EVIDENCE.tsv',sep='\t',index=False)
 for rel,r in manifest.items():
  r['sha256_after']=sha(root/rel);assert r['sha256_before']==r['sha256_after'],rel
 pd.DataFrame(manifest.values()).to_csv(review/'SPATIAL_INPUT_MANIFEST.tsv',sep='\t',index=False)
 pd.DataFrame(derived).to_csv(sources/'SPATIAL_DERIVED_SOURCE_MAP.tsv',sep='\t',index=False)
 status=dict(status='SPATIAL_READONLY_EXTRACTION_PASS_HISTOLOGY_ADJUDICATION_PENDING',started_or_finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start,n_spatial_spots=len(domains),n_candidate_spot_rows=len(expressions),n_domains=len(annotation),n_donors=8,n_actual_local_pairs=len(local),all_sources_unchanged=True,models_run=False,P_FDR_recomputed=False,images_modified=False,python_visuals_generated=False,representative_donors=['P8','P9'],representative_status='PROVISIONAL_PENDING_ROOT_HISTOLOGY_REVIEW',expression_units='natural log1p(CP10000(sum same-HGNC original rows after background clipping)); native domain means use log2 and are not mixed',source_files=len(manifest),source_data_tables=len(derived))
 (sources/'SPATIAL_EXTRACTION_STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
 print(json.dumps(status))
if __name__=='__main__':main()
