from pathlib import Path
import csv,json,hashlib,collections
work=Path(__file__).resolve().parent
source=work/'draft_review_2/04_Review/DOMAIN_ANNOTATION_EVIDENCE.tsv'
rows=list(csv.DictReader(source.open(encoding='utf-8'),delimiter='\t'));assert len(rows)==44
out=work/'final_domain_annotation';out.mkdir(exist_ok=True)
observed={
'P1':'Open folded C/S-shaped tissue outline; mucosal folds and intervening spaces visible. The domain map does not establish a single lineage.',
'P5':'Coiled section with an evident mucosal rim; gross structure is interpretable. Low formal-chain ESS remains a separate limitation.',
'P8':'Near-closed coiled section with clear mucosal folds and an inner opening; retained as a representative together with P9, independently of candidate expression.',
'P9':'Coiled folded section with a clear dense mucosal rim and inner spaces; retained with P8 for architectural display and high archived seed concordance.',
'P10':'Coiled tissue with uneven inner folded structures; no gene-based region cropping. The strongest NCAM1 donor was not used to choose Figure4 representatives.',
'P11':'Large central opening and a folded mucosal rim; lower archived independent-seed ARI limits domain-specific interpretation.',
'P12':'Open curved section with a partial rim and fewer retained spots; image geometry is retained rather than cropped to a favorable region.',
'P13':'Curved open mucosal arch with uneven retained coverage; no precise cellular identity inferred from a domain color.'}
for r in rows:
 support=json.loads(r['supporting_markers']);other=json.loads(r['countermarkers'])
 r['original_screening_label']=r['proposed_label']
 if len(support)>1:r['proposed_label']='Mixed marker-enriched region'
 elif len(support)==0:r['proposed_label']='Domain '+r['domain']+' / unresolved'
 else:
  key=next(iter(support))
  if key=='Paneth' and {'DEFA5','DEFA6'}.issubset(support[key]):r['proposed_label']='DEFA5/DEFA6-associated transcript program; cell identity unresolved'
  else:r['proposed_label']=key+'-associated marker program; cell identity unresolved'
 r['countermarkers']=json.dumps(dict(competing_supported_programs=support if len(support)>1 else {},other_singleton_support=other,interpretation='Competing marker support is not a validated negative lineage-exclusion test'),ensure_ascii=False)
 r['histology_support']='Gross H&E and registered domain map actually reviewed: '+observed[r['donor']]
 r['confidence']='DESCRIPTIVE_MARKER_CONTEXT_ONLY; no confirmed domain-specific cell identity'
 r['source_reference']='Native scran 1.28.2 findMarkers(direction=up,pval.type=any); Harnik 2024 DOI10.1038/s41586-024-07793-3 provides tissue/reference context, not verified author labels for these 44 domains'
 r['limitations']='summary.logFC is the effect from the lowest-P pairwise comparison under pval.type=any, not a minimum or an exclusive-domain effect. Cluster-derived FDR is used descriptively. Visium spots mix cells; gross image review is not specialist histopathology or cell-type validation.'
 r['selection_rule']='Fixed manuscript-stage post hoc marker lists; >=2 original native BH<0.05 positive summary.logFC markers per program; no new enrichment test or model fit'
 r['annotation_stage']='POST_HOC_MANUSCRIPT_20261108'
 r['visual_review_source']='R-rendered S4_1.png / S4_2.png; full original image orientation and frozen domain coordinates; donor-level architectural assessment only'
keys=list(rows[0])
with (out/'DOMAIN_ANNOTATION_EVIDENCE.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
counts=collections.Counter(r['proposed_label'] for r in rows)
result=dict(status='POSTHOC_MARKER_CONTEXT_AND_GROSS_IMAGE_REVIEW_COMPLETE',n_domains=44,label_counts=dict(counts),representative_donors=['P8','P9'],selection_basis='Candidate-blind choice supported by coherent folded tissue architecture and archived independent-seed ARI, not NCAM1 expression',screening_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),final_annotation_sha256=hashlib.sha256((out/'DOMAIN_ANNOTATION_EVIDENCE.tsv').read_bytes()).hexdigest(),all_original_marker_statistics_unchanged=True,models_run=False,cell_type_identities_confirmed=False,histology_scope='Eight complete donor H&E/domain panels visually inspected; no microscopic lineage or pathology certification')
(out/'DOMAIN_ANNOTATION_ADJUDICATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=True))
