#!/usr/bin/env python3
"""Append-only command evidence. This runner never imports analysis results."""
import argparse,datetime,fcntl,hashlib,json,os,shlex,socket,subprocess,sys,time,uuid
from pathlib import Path
FIELDS=['event','module','source_command','input_hashes','start_iso8601','finish_iso8601','duration_seconds','exit_code','hostname','software_version','run_id','is_new_fit','reused_only_immutable_raw_inputs','output_hashes','notes','stdout_log','stderr_log']
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 a=argparse.ArgumentParser();a.add_argument('--root',type=Path,required=True);a.add_argument('--module',required=True)
 a.add_argument('--input',action='append',default=[]);a.add_argument('--output',action='append',default=[])
 a.add_argument('--software-version',default='see env locks');a.add_argument('--is-new-fit',choices=['0','1'],required=True)
 a.add_argument('--notes',default='');a.add_argument('--env',action='append',default=[]);a.add_argument('command',nargs=argparse.REMAINDER)
 x=a.parse_args();root=x.root.resolve();cmd=x.command[1:] if x.command[:1]==['--'] else x.command
 assert root.name=='IBS_A_OFFICIAL_RERUN_20261105' and cmd
 forbidden=['/IBS_A_computation/results/','/IBS_A_computation/figures/','/IBS_A_computation/reports/']
 assert not any(k in s for s in cmd for k in forbidden),'Old analysis input forbidden'
 def hashes(names,missing=False):
  out={}
  for name in names:
   p=Path(name);p=p if p.is_absolute() else root/p;p=p.resolve()
   assert p==root or root in p.parents,('Input/output escaped new root',p)
   if not p.exists():
    if missing:out[str(p.relative_to(root))]='MISSING';continue
    raise FileNotFoundError(p)
   files=sorted(q for q in p.rglob('*') if q.is_file()) if p.is_dir() else [p]
   for q in files:out[q.relative_to(root).as_posix()]=digest(q)
  return out
 runid=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8]
 logs=root/'logs'/x.module;logs.mkdir(parents=True,exist_ok=True)
 stdout=logs/(runid+'.stdout.log');stderr=logs/(runid+'.stderr.log')
 input_digests=hashes(x.input)
 row=dict.fromkeys(FIELDS,'NA');start=datetime.datetime.now(datetime.timezone.utc).isoformat();mono=time.monotonic()
 row.update(event='START',module=x.module,source_command=shlex.join(cmd),input_hashes=json.dumps(input_digests,sort_keys=True),
            start_iso8601=start,hostname=socket.gethostname(),software_version=x.software_version,run_id=runid,is_new_fit=x.is_new_fit,
            reused_only_immutable_raw_inputs='TRUE',notes=x.notes,stdout_log=str(stdout.relative_to(root)),stderr_log=str(stderr.relative_to(root)))
 trace=root/'EXECUTION_TRACE.tsv'
 def append(record):
  import csv
  with trace.open('a+',newline='') as f:
   fcntl.flock(f,fcntl.LOCK_EX);f.seek(0,2);new=f.tell()==0
   w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\n')
   if new:w.writeheader()
   w.writerow(record);f.flush();os.fsync(f.fileno());fcntl.flock(f,fcntl.LOCK_UN)
 append(row);(logs/(runid+'.launch.json')).write_text(json.dumps(row,indent=2)+'\n')
 env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']})
 for item in x.env:
  k,v=item.split('=',1);env[k]=v
 try:
  with stdout.open('wb') as so,stderr.open('wb') as se:
   p=subprocess.run(cmd,cwd=root,env=env,stdout=so,stderr=se)
  rc=p.returncode
 except BaseException as e:
  rc=255
  with stderr.open('ab') as f:f.write((repr(e)+'\n').encode())
 row.update(event='FINISH',finish_iso8601=datetime.datetime.now(datetime.timezone.utc).isoformat(),duration_seconds=f'{time.monotonic()-mono:.6f}',exit_code=rc,output_hashes=json.dumps(hashes(x.output,missing=True),sort_keys=True))
 append(row);(logs/(runid+'.finish.json')).write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True)
 raise SystemExit(rc)
if __name__=='__main__':main()
