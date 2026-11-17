"""Postprocessing of existing native results only. No model fitting or BH reruns."""
from pathlib import Path
import argparse,json,csv,gzip,hashlib,shutil,datetime
import pandas as pd
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
root=Path(a.root);out=Path(a.out);dest=out/'02_Main_Figures/source_data/core';dest.mkdir(parents=True,exist_ok=True)
tables=out/'03_Supplementary/native_tables';tables.mkdir(parents=True,exist_ok=True)
assert not (dest/'CORE_SOURCE_STATUS.json').exists()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources={};outputs=[]
def track(rel):
 p=root/rel;assert p.is_file(),str(p)
 if rel not in sources:sources[rel]=dict(source_path=str(p),bytes=p.stat().st_size,sha256_before=sha(p))
 return p
def tab(rel):return pd.read_csv(track(rel),sep='\t',dtype=str,keep_default_na=False)
def js(rel):return json.loads(track(rel).read_text())
def save(name,df,rule,inputs):
 p=dest/name;assert not p.exists();df.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression='gzip' if name.endswith('.gz') else None)
 outputs.append(dict(output=str(p.relative_to(out)),rows=len(df),columns=len(df.columns),sha256=sha(p),filter_rule=rule,source_paths=';'.join(inputs)))
def copy_native(rel):
 p=track(rel);target=tables/rel
 if p.suffix!='.gz':target=target.with_name(target.name+'.gz')
 target.parent.mkdir(parents=True,exist_ok=True)
 assert not target.exists(),str(target)
 if p.suffix=='.gz':shutil.copyfile(p,target)
 else:
  with p.open('rb') as f,target.open('xb') as z:
   with gzip.GzipFile(fileobj=z,mode='wb',filename='',mtime=0,compresslevel=9) as g:shutil.copyfileobj(f,g)
 # The byte stream restored from gzip is verified against the original source.
 if p.suffix!='.gz':
  h=hashlib.sha256()
  with gzip.open(target,'rb') as g:
   for b in iter(lambda:g.read(4*1024*1024),b''):h.update(b)
  assert h.hexdigest()==sources[rel]['sha256_before']
 return target

qc=js('results/gwas/qc/summary.json');magma=tab('results/gwas/magma_all_genes.tsv')
assert len(magma)==18333 and sum(pd.to_numeric(magma.BH_FDR)<.05)==77
freeze=js('results/integrated/CANDIDATE_FREEZE.json')
ranking_rel=freeze['statistical_ranking_file'];ranking=tab(ranking_rel)
aa=ranking[ranking.evidence_grade.eq('A')].sort_values('statistical_rank',key=lambda x:pd.to_numeric(x))
assert aa.gene_symbol.tolist()==['NCAM1','PCLO','SKIC2','C4B','VWA7','RBM33']
genes=aa.gene_symbol.tolist();ids=dict(zip(aa.gene_symbol,aa.identity_key))
copy_native(ranking_rel)
wide='results/integrated/continuous/ALL_GENE_EVIDENCE_WIDE.tsv.gz'
if (root/wide).exists():copy_native(wide)
for rel in ['results/gwas/magma_all_genes.tsv','results/gwas/magma_BH_FDR05_genes.tsv']:
 if (root/rel).exists():copy_native(rel)
save('MAGMA_ALL.tsv',magma,'All original tested genes; unchanged P/BH/rank',['results/gwas/magma_all_genes.tsv'])
save('MAGMA_FDR05.tsv',magma[pd.to_numeric(magma.BH_FDR)<.05],'Native BH_FDR < 0.05; no recalculation',['results/gwas/magma_all_genes.tsv'])
save('A6_ORIGINAL_RANK.tsv',aa,'Original six A-grade rows, all original native/statistical columns',[ranking_rel])
projection_cols=[k for k in ranking.columns if k in ['statistical_rank','identity_key','gene_symbol','hgnc_entrez_id','hgnc_ensembl_gene_id','evidence_grade','selection_eligible','selection_eligibility_unknown','MAGMA_P','MAGMA_BH_FDR','MAGMA_rank','observed_modules'] or (k.endswith('_log2FoldChange') or k.endswith('_padj') or k.endswith('_measurement_status'))]
save('ALL_CANDIDATE_RANKING_PROJECTION.tsv.gz',ranking[projection_cols],'All frozen ranked identities; audit projection only; full original wide and full ranking retained separately',[ranking_rel])

snp_rel='results/gwas/qc/manhattan_points.tsv.gz';snp=tab(snp_rel);qq_rel='results/gwas/qc/qq_quantiles.tsv';qq=tab(qq_rel)
chrommax={}
for c in range(1,23):
 a1=pd.to_numeric(snp.loc[pd.to_numeric(snp.CHR)==c,'BP'])
 a2=pd.to_numeric(magma.loc[pd.to_numeric(magma.CHR)==c,'STOP'])
 chrommax[c]=max(a1.max(),a2.max())
offset={};total=0;centers=[]
for c in range(1,23):offset[c]=total;centers.append(dict(CHR=c,center_Mb=(total+chrommax[c]/2)/1e6));total+=chrommax[c]+5e6
for df,bp in [(snp,'BP'),(magma,'START')]:
 xx=pd.to_numeric(df[bp]) if bp=='BP' else (pd.to_numeric(df.START)+pd.to_numeric(df.STOP))/2
 df['plot_x_Mb']=(xx+pd.to_numeric(df.CHR).map(offset))/1e6
 df['plot_y_minus_log10_P']=-np.log10(pd.to_numeric(df.P));assert np.isfinite(df.plot_y_minus_log10_P).all()
save('FIG2A_GWAS_POINTS.tsv.gz',snp,'Existing fixed CRC32 modulo-1000 background plus every P<5e-8; no SNP selected by gene or manuscript narrative',[snp_rel,'results/gwas/qc/summary.json'])
save('FIG2B_GWAS_QQ.tsv',qq,'Existing full-universe order statistics; no refit or independent-SNP confidence band',[qq_rel])
save('FIG2C_MAGMA_POINTS.tsv',magma,'All 18333 native genes; original P and BH, observed-coordinate cumulative plotting only',['results/gwas/magma_all_genes.tsv'])
save('GENOME_AXIS_CENTERS.tsv',pd.DataFrame(centers),'Observed chromosome maxima from original SNP cache and gene table; 5Mb visual separators',[snp_rel,'results/gwas/magma_all_genes.tsv'])
selected=magma[pd.to_numeric(magma['rank'])<=12].copy()
selected=pd.concat([selected,magma[magma.hgnc_id.isin(ids.values())]]).drop_duplicates('gene_id').sort_values('P',key=lambda x:pd.to_numeric(x))
selected['is_A6']=selected.hgnc_id.isin(ids.values());save('FIG2D_GENE_DISPLAY.tsv',selected,'Top12 original MAGMA P ranks union all six frozen A candidates; original rank never changed',['results/gwas/magma_all_genes.tsv',ranking_rel])

comparisons=['Duo_IBSD_vs_HV','Duo_IBSC_vs_HV','Jej_IBSD_vs_HV','Jej_IBSC_vs_HV']
cohorts=['main','qc_low_750k','qc_high_1250k'];deg_counts=[];all_gsea=[]
for cohort in cohorts:
 folder='results/bulk/main' if cohort=='main' else 'results/bulk/sensitivity/'+cohort
 for contrast in comparisons:
  rel=folder+'/'+contrast+'.tsv';df=tab(rel);assert len(df)==64253
  copy_native(rel)
  p=pd.to_numeric(df.padj,errors='coerce');lfc=pd.to_numeric(df.log2FoldChange,errors='coerce');sig=p<.05
  deg_counts.append(dict(cohort=cohort,contrast=contrast,total_feature_rows=len(df),native_finite_P=pd.to_numeric(df.pvalue,errors='coerce').notna().sum(),native_finite_FDR=p.notna().sum(),up=int((sig&(lfc>0)).sum()),down=int((sig&(lfc<0)).sum()),FDR05=int(sig.sum()),native_FDR_NA=int(p.isna().sum()),rule='native BH<0.05; no |LFC| threshold; NA not counted as tested-negative',source_path=rel,source_sha256=sources[rel]['sha256_before']))
  # Full figure source rows retained, with no new fitting, NA imputation or filtering.
  full=df.copy();full['cohort']=cohort;full['contrast']=contrast
  save('VOLCANO_'+cohort+'_'+contrast+'.tsv.gz',full,'All 64253 original rows; rendering only finite native LFC and FDR; NA rows remain in source',[rel])
  gr=folder+'/'+contrast+'_GSEA_ALL.tsv'
  if (root/gr).exists():
   gd=tab(gr);copy_native(gr);gd['cohort']=cohort;gd['contrast']=contrast;all_gsea.append(gd)
for tissue in ['Duo','Jej']:
 for suffix in ['PCA_SAMPLE_QC.tsv','PCA_VARIANCE.tsv','FITTED_SAMPLE_METADATA.tsv','MODEL_AUDIT.json']:
  rel='results/bulk/main/'+tissue+'_'+suffix
  if suffix.endswith('.json'):copy_native(rel)
  else:save(tissue+'_'+suffix,tab(rel),'All original primary sample/PCA/model metadata; no PCA rerun or outlier deletion',[rel])
save('DEG_COUNTS_ALL_QC.tsv',pd.DataFrame(deg_counts),'Count original native BH<.05 by LFC sign; all primary and two fixed cohorts',[r['source_path'] for r in deg_counts])
gsea=pd.concat(all_gsea,ignore_index=True);save('GSEA_ALL_COHORTS.tsv.gz',gsea,'Every declared/tested/untested pathway row for every existing cohort; native correction columns unchanged',['results/bulk/'+c for c in cohorts])
primary=gsea[gsea.cohort.eq('main')].copy();chosen=set()
for (contrast,collection),sub in primary.groupby(['contrast','collection']):
 n=sub.copy();n['_q']=pd.to_numeric(n.padj,errors='coerce');n=n[np.isfinite(n._q)].sort_values(['_q','pathway'],kind='stable')
 if len(n):chosen.add((collection,n.iloc[0].pathway))
display=primary[[ (c,p) in chosen for c,p in zip(primary.collection,primary.pathway)]].copy()
save('FIG3D_GSEA_DISPLAY.tsv',display,'Manuscript-stage display: lowest finite native pooled BH within each contrast and collection, union <=8 pathways; all four contrasts retained; no pain-term selection',['GSEA_ALL_COHORTS.tsv.gz'])
sens='results/bulk/sensitivity/QC_SENSITIVITY_COMPARISON_SUMMARY.tsv';save('QC_SENSITIVITY_SUMMARY_NATIVE.tsv',tab(sens),'Complete existing native sensitivity summaries; overlapping cohorts not replication',[sens])
sample_rel='data/metadata/bulk/SAMPLE_MANIFEST.tsv';sm=tab(sample_rel);save('BULK_SAMPLE_MANIFEST_NATIVE.tsv',sm,'All71 count columns plus GEO mapping and inclusion/exclusion; no inferred subject IDs',[sample_rel])
groups=tab('data/metadata/bulk/GROUP_COUNTS_AND_SENSITIVITY.tsv');save('BULK_GROUP_COUNTS_NATIVE.tsv',groups,'All six deposited/published/main/QC group count records',['data/metadata/bulk/GROUP_COUNTS_AND_SENSITIVITY.tsv'])
flow_rel='data/metadata/SPATIAL_SPOT_FLOW.tsv';flow=tab(flow_rel);save('SPATIAL_SPOT_FLOW_NATIVE.tsv',flow,'All original author-logic preprocessing/QC stages; no spot filtering rerun',[flow_rel])
for rel in ['data/metadata/SPATIAL_CLINICAL_METADATA_ORIGINAL.tsv','data/metadata/bulk/GEO_samples_all.tsv','data/metadata/bulk/QC_PROTOCOL_FROZEN.json','results/integrated/candidate_review/A_CANDIDATES_ALL_NATIVE_COLUMNS.tsv.gz']:
 if (root/rel).is_file():copy_native(rel)
f1=[dict(resource='genetics',accession='GCST90016564',quantity='cases_max',value=qc['numeric_ranges']['N_CASE'][1],source='results/gwas/qc/summary.json:numeric_ranges.N_CASE.maximum; source-study sample count'),dict(resource='genetics',accession='GCST90016564',quantity='controls_max',value=qc['numeric_ranges']['N_CONTROL'][1],source='results/gwas/qc/summary.json:numeric_ranges.N_CONTROL.maximum; source-study sample count'),dict(resource='patient',accession='GSE166869',quantity='main_tissue_libraries',value=sum(sm.main_analysis_pass.str.lower().eq('true')),source=sample_rel),dict(resource='spatial',accession='Zenodo10715015',quantity='donors',value=8,source=flow_rel),dict(resource='candidate',accession='HGNC harmonized frozen rank',quantity='A_grade_genes',value=len(aa),source=ranking_rel)]
save('FIGURE1_RESOURCES.tsv',pd.DataFrame(f1),'Actual recorded resource counts; genetic max N agrees with referenced source-study maximum, not per-SNP constant N',['results/gwas/qc/summary.json',sample_rel,flow_rel,ranking_rel])
facts=dict(GWAS_counts=qc['counts'],GWAS_lambda_GC=qc['plotting']['lambda_gc'],MAGMA_n=18333,MAGMA_FDR05=77,MAGMA_native_source_sha256=sources['results/gwas/magma_all_genes.tsv']['sha256_before'],DEG_counts=deg_counts,GSEA_counts=[dict(cohort=co,contrast=ct,declared=len(rr),tested=int(rr.tested.str.lower().eq('true').sum()),FDR05=int((pd.to_numeric(rr.padj,errors='coerce')<.05).sum())) for (co,ct),rr in gsea.groupby(['cohort','contrast'])],A6_order=genes,MHC_reference=dict(assembly='GRCh37',chromosome=6,start=28477797,stop=33448354,url='https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37',interpretation='GRC assembly region annotation, not an LD block or fine-mapping result'))
(dest/'CORE_FACTS.json').write_text(json.dumps(facts,indent=2,default=int))
for rel,row in sources.items():row['sha256_after']=sha(root/rel);assert row['sha256_before']==row['sha256_after']
pd.DataFrame(sources.values()).to_csv(out/'04_Review/CORE_INPUT_MANIFEST.tsv',sep='\t',index=False)
pd.DataFrame(outputs).to_csv(dest/'CORE_FIGURE_SOURCE_MAP.tsv',sep='\t',index=False)
status=dict(status='CORE_READONLY_POSTPROCESSING_PASS',source_files=len(sources),outputs=len(outputs),original_sources_unchanged=True,models_run=False,original_P_FDR_unchanged=True,finished_at=datetime.datetime.now().astimezone().isoformat())
(dest/'CORE_SOURCE_STATUS.json').write_text(json.dumps(status,indent=2));print(json.dumps(status))
