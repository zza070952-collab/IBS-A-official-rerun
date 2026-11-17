#!/usr/bin/env python3
"""Read-only check of archived source identity; does not import analysis code."""
from pathlib import Path
import csv,hashlib,sys

root=Path(__file__).resolve().parents[1]
with (root/'provenance/SOURCE_CODE_MANIFEST.tsv').open(encoding='utf-8',newline='') as stream:
    rows=list(csv.DictReader(stream,delimiter='\t'))
failures=[]
for row in rows:
    relative=Path(row['repository_path'])
    if relative.is_absolute() or '..' in relative.parts:
        failures.append((str(relative),'unsafe path'));continue
    path=root/relative
    if not path.is_file():
        failures.append((str(relative),'missing'));continue
    with path.open('rb') as stream:
        digest=hashlib.sha256()
        for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
    if path.stat().st_size!=int(row['bytes']) or digest.hexdigest()!=row['sha256']:
        failures.append((str(relative),'size or SHA256 mismatch'))
if failures:
    for path,reason in failures:print('FAIL',path,reason)
    sys.exit(1)
print('PASS:',len(rows),'archived source/parameter/environment files; no models run')
