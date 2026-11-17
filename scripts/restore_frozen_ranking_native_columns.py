#!/usr/bin/env python3
"""One explicit serialization repair; no ranking, statistics, models or selection.

Only --apply mutates the two canonical ranking files, after complete staged
string equality checks. Failure artifacts and the one-use guard are retained.
Run through run_logged.py with is-new-fit=0. Never rerun the ranker for this fix.
"""
import argparse
import collections
import copy
import csv
import datetime as dt
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import socket
import sqlite3
import sys
import time
import traceback

ROOT = '/home/asset/projects/IBS_A_OFFICIAL_RERUN_20261105'
WIDE = 'results/integrated/continuous/CONTINUOUS_EVIDENCE_WIDE.tsv.gz'
RANK = 'results/integrated/statistical_ranking/STATISTICAL_RANKING_ALL_GENES.tsv.gz'
FREEZE = 'results/integrated/statistical_ranking/STATISTICAL_RANKING_FREEZE.json'
RULE = 'results/integrated/CANDIDATE_RANKING_RULE.json'
SOURCE = 'results/integrated/continuous/INTEGRATION_SOURCE_FREEZE.json'
RANKER = 'scripts/rank_and_freeze_candidates.py'
PROOF = 'logs/PKG_RESOURCES_COMPATIBILITY_FIX_20261105/ENVIRONMENT_FIX_STATUS.json'
AUDIT = 'logs/RANKING_NATIVE_COLUMN_RESTORATION_20261105'
GUARD = 'results/integrated/statistical_ranking/NATIVE_COLUMN_RESTORATION_ONCE.json'
GATES = ('PHASE0_GATE', 'PHASE1_GATE', 'SOFTWARE_GATE', 'GWAS_MAGMA_GATE', 'BULK_GATE',
         'BAYESSPACE_GATE', 'SPATIALDM_GATE', 'INTEGRATION_CONTINUOUS_GATE')
EXPECTED = {
 WIDE: 'af3a66731a0a5050bd92c3e7dfc7a6bf47c63da3a174e5d2d70a6cf1c16e2e76',
 RANK: 'a1d1c5a8c65c8257e79ad04caf0999896f9676b26c90cb327ada185b1ef8f86a',
 FREEZE: 'd67070f50ce81dfc41c017d6fb7c4387de45e2e9a8e511262363706745bacad9',
 RULE: 'a6cb1f8371d69e92c8ec3328b8d7be5f49226873763cfbdbd68b79982ec7d204',
 SOURCE: '07ae51d34b0b926d1f915972d914d17df563dc9f6af1ea24161ff37d834e5102',
 RANKER: 'a5b496d2be52c6de2e4f132aa7e01ad1ee411d9174c782498c7d6fd5e0774118',
 PROOF: 'b92d8fcaf983295a43d81d2cab6a23991d757f9d6b5ac93c2cd543463c88c199',
 'scripts/spatial_bayesspace.R': 'a6f0d8f2cde01b1edb7dd34b92b66b753a89425727488e1d76246d5c2eee8d13',
 'scripts/spatial_spatialdm.py': 'cb6d291612062b339e5e40c0162c4b76923780eb9002e455862bc71c107206ab',
 'scripts/integrate_continuous_evidence.py': '01459c92879f311170f3b7a2fbb42825be7f267bda6bffff4b531f538113dc5a',
 'scripts/run_logged.py': '4458bf4600baf0420b7ebdb10f2e4a274fce9d9fb4c3c293b7d11e0c2231bc56',
 'env/r_exec.sh': '0e217a53f327e20280a912fc59f722035f19ccab053326fb31cf4a0dffdbc68d',
}


def utc(): return dt.datetime.now(dt.timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b''): h.update(b)
    return h.hexdigest()


def need(test, message):
    if not test: raise ValueError(message)


def path(root, relative):
    p = root / relative
    need(p.resolve() == p and root in p.parents and not p.is_symlink(), 'Noncanonical/symlink path: ' + relative)
    return p


def js(p): return json.loads(p.read_text(encoding='utf-8-sig'))


def dump(p, data):
    temporary = p.with_name(p.name + '.tmp')
    with temporary.open('x', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(temporary, p)


def rows(p):
    with gzip.open(p, 'rt', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter='\t', strict=True)
        header = next(reader)
        need(len(header) == len(set(header)), 'Duplicate table columns')
        yield header
        for row in reader:
            need(len(row) == len(header), 'Ragged table row')
            yield row


def load_wide(db, wide, expected_rows, expected_columns):
    source = rows(wide); header = next(source)
    need(len(header) == expected_columns and 'identity_key' in header, 'Wide schema/count differs')
    key_at = header.index('identity_key')
    db.execute('CREATE TABLE native (identity TEXT PRIMARY KEY, cells TEXT NOT NULL)')
    count = 0
    for row in source:
        need(bool(row[key_at]), 'Empty canonical identity')
        db.execute('INSERT INTO native VALUES (?,?)', (row[key_at], json.dumps(row, ensure_ascii=False, separators=(',', ':'))))
        count += 1
    db.commit()
    need(count == expected_rows, 'Wide row count differs')
    return header


def native_row(db, identity):
    hit = db.execute('SELECT cells FROM native WHERE identity=?', (identity,)).fetchone()
    need(hit is not None, 'Rank identity absent from native wide: ' + identity)
    return json.loads(hit[0])


def restore_to_stage(db, wide_header, old_rank, stage, expected_rows, expected_columns):
    old = rows(old_rank); header = next(old)
    need(len(header) == expected_columns and set(wide_header).issubset(header), 'Rank schema/count differs')
    original_at = [header.index(name) for name in wide_header]
    need([name for name in header if name in set(wide_header)] == wide_header, 'Original column names/order differ')
    key_at = header.index('identity_key')
    differences = collections.Counter(); seen = set(); count = 0
    with stage.open('xb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0, compresslevel=6) as gz:
            import io
            with io.TextIOWrapper(gz, encoding='utf-8', newline='') as text:
                writer = csv.writer(text, delimiter='\t', lineterminator='\n'); writer.writerow(header)
                for row in old:
                    identity = row[key_at]
                    need(identity not in seen, 'Duplicate ranked identity'); seen.add(identity)
                    native = native_row(db, identity); replacement = list(row)
                    for i, rank_at in enumerate(original_at):
                        if row[rank_at] != native[i]: differences[wide_header[i]] += 1
                        replacement[rank_at] = native[i]
                    writer.writerow(replacement); count += 1
    need(count == expected_rows, 'Ranking row count differs')
    need(len(seen) == db.execute('SELECT COUNT(*) FROM native').fetchone()[0], 'Identity universes differ')
    return header, dict(differences)


def validate_stage(db, wide_header, old_rank, stage, expected_rows, expected_columns):
    old, new = rows(old_rank), rows(stage)
    old_header, new_header = next(old), next(new)
    need(old_header == new_header and len(new_header) == expected_columns, 'Staged column/order differs')
    original_at = [new_header.index(name) for name in wide_header]
    original_set = set(original_at); other_at = [i for i in range(len(new_header)) if i not in original_set]
    key_at = new_header.index('identity_key'); rank_at = new_header.index('statistical_rank'); grade_at = new_header.index('evidence_grade')
    old_order, new_order = hashlib.sha256(), hashlib.sha256()
    count = 0; seen = set()
    for before, after in itertools.zip_longest(old, new):
        need(before is not None and after is not None, 'Staged row count differs')
        need(before[key_at] == after[key_at] and before[key_at] not in seen, 'Staged identity/order changed')
        seen.add(before[key_at]); native = native_row(db, before[key_at])
        need(all(after[j] == native[i] for i, j in enumerate(original_at)), 'Staged native source cell mismatch')
        need(all(after[i] == before[i] for i in other_at), 'Nonoriginal frozen rank/feature/grade cell changed')
        for row, h in ((before, old_order), (after, new_order)):
            h.update((json.dumps([row[i] for i in other_at], ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf-8'))
        need(after[rank_at] == str(count + 1), 'Original ranks are not the frozen contiguous order')
        need(before[grade_at] == after[grade_at], 'Evidence grade changed')
        count += 1
    need(count == expected_rows and len(seen) == db.execute('SELECT COUNT(*) FROM native').fetchone()[0], 'Full identity/row coverage differs')
    need(old_order.hexdigest() == new_order.hexdigest(), 'Frozen nonoriginal ordered-cell digest changed')
    return {'rows': count, 'columns': len(new_header), 'native_columns': len(wide_header), 'nonoriginal_columns': len(other_at),
            'all_native_cells_equal_wide_strings': True, 'all_nonoriginal_cells_equal_old_rank_strings': True,
            'ordered_identity_rank_grade_unchanged': True, 'ordered_nonoriginal_cells_sha256': new_order.hexdigest()}


def preflight(root):
    need(str(root) == ROOT and os.name == 'posix', 'Actual repair is restricted to the exact authorized Linux new root')
    need(not path(root, GUARD).exists() and not path(root, AUDIT).exists(), 'Restoration was already attempted; no automatic reentry')
    for rel in ('results/integrated/CANDIDATE_FREEZE.json', 'manifest/CANDIDATE_FREEZE_GATE.json'):
        need(not path(root, rel).exists(), 'Candidate decision already exists; stop for explicit review')
    protected = dict(EXPECTED)
    for rel, expected in EXPECTED.items(): need(digest(path(root, rel)) == expected, 'Expected SHA differs: ' + rel)
    for gate in GATES:
        rel = 'manifest/' + gate + '.json'; p = path(root, rel)
        need(js(p).get('status') == 'PASS', 'Required gate not PASS: ' + gate); protected[rel] = digest(p)
    coordinator_rel = 'manifest/COORDINATOR_STATUS.json'; coordinator = js(path(root, coordinator_rel))
    need(coordinator.get('status') == 'ALL_MODELS_AND_CONTINUOUS_COMPLETE_WAITING_POSTFREEZE_HUMAN_REVIEW', 'Coordinator not at actual completed model/continuous boundary')
    protected[coordinator_rel] = digest(path(root, coordinator_rel))
    protected['EXECUTION_TRACE.tsv'] = digest(path(root, 'EXECUTION_TRACE.tsv'))
    old = js(path(root, FREEZE)); gate = js(path(root, 'manifest/INTEGRATION_CONTINUOUS_GATE.json'))
    need(old['status'] == 'FROZEN' and old['all_canonical_rows_retained'] == 45289 and old['candidate_selection_performed'] is False, 'Unexpected ranking freeze')
    need('native_column_restoration' not in old, 'Restoration metadata already exists')
    for field, rel in (('ranking_rule', RULE), ('source_freeze', SOURCE), ('wide_table', WIDE), ('ranking_table', RANK)):
        file_field = field + '_file' if field in ('ranking_rule', 'source_freeze') else field
        need(old[file_field] == rel and old[field + '_sha256'] == EXPECTED[rel], 'Freeze binding differs: ' + field)
    need(gate['wide_table'] == WIDE and gate['wide_table_sha256'] == EXPECTED[WIDE]
         and gate['source_freeze'] == SOURCE and gate['source_freeze_sha256'] == EXPECTED[SOURCE], 'Continuous gate binding differs')
    source = js(path(root, SOURCE)); need(source['status'] == 'FROZEN' and source['ranking_performed'] is False, 'Continuous source not frozen unranked')
    # Check the complete previously frozen statistical input catalogue without
    # importing any model or changing a statistic. Never rewrite these sources.
    for rel, expected in old['input_sha256'].items():
        need(digest(path(root, rel)) == expected, 'Frozen ranking input changed: ' + rel)
        protected[rel] = expected
    return old, protected


def unchanged(root, protected, exclude=()):
    for rel, expected in protected.items():
        if rel not in exclude: need(digest(path(root, rel)) == expected, 'Protected source changed during repair: ' + rel)


def run(root, apply):
    old_freeze, protected = preflight(root)
    if not apply:
        print(json.dumps({'status': 'PREFLIGHT_VALIDATED_NO_RESTORATION_OR_OUTPUT', 'ranking_evaluated': False})); return
    started, tick = utc(), time.monotonic()
    audit = path(root, AUDIT); audit.mkdir(exist_ok=False)
    state = {'status': 'STARTED', 'started_at': started, 'source_command': sys.argv, 'hostname': socket.gethostname(),
             'helper_sha256': digest(Path(__file__)), 'ranking_evaluated': False, 'new_inferential_tests': False,
             'models_run': False, 'candidate_selection_performed': False, 'protected_inputs_sha256': protected}
    guard = path(root, GUARD)
    with guard.open('x', encoding='utf-8') as f: json.dump({'status': 'ONE_EXPLICIT_ATTEMPT_STARTED', 'started_at': started, 'audit': AUDIT}, f)
    state_path = audit / 'RESTORATION_AUDIT.json'
    try:
        archive = audit / 'original'; archive.mkdir()
        preserved = {}
        for rel in (RANK, FREEZE, RANKER):
            dest = archive / Path(rel).name
            with path(root, rel).open('rb') as src, dest.open('xb') as dst: shutil.copyfileobj(src, dst, 4 * 1024 * 1024)
            need(digest(dest) == EXPECTED[rel], 'Original archive SHA differs')
            preserved[rel] = {'path': dest.relative_to(root).as_posix(), 'sha256': digest(dest), 'bytes': dest.stat().st_size}
        dump(audit / 'ORIGINAL_ARCHIVE_SHA256.json', preserved)
        staged = audit / 'staged'; staged.mkdir(); stage_rank = staged / 'RESTORED_RANK.tsv.gz'
        with sqlite3.connect(str(audit / 'native_strings.sqlite')) as db:
            header = load_wide(db, path(root, WIDE), 45289, 487)
            rank_header, counts = restore_to_stage(db, header, path(root, RANK), stage_rank, 45289, 543)
            validation = validate_stage(db, header, path(root, RANK), stage_rank, 45289, 543)
        with (audit / 'RESTORED_COLUMN_DIFF_COUNTS.tsv').open('x', encoding='utf-8', newline='') as f:
            writer = csv.writer(f, delimiter='\t', lineterminator='\n'); writer.writerow(['column', 'cells_restored_from_wide'])
            for column in header: writer.writerow([column, counts.get(column, 0)])
        restored_at = utc(); new_hash = digest(stage_rank)
        new_freeze = copy.deepcopy(old_freeze)
        new_freeze['ranking_table_sha256'] = new_hash
        new_freeze['native_column_restoration'] = {'restoration_at': restored_at, 'audit_path': AUDIT,
            'old_ranking_table_sha256': EXPECTED[RANK], 'old_ranking_freeze_sha256': EXPECTED[FREEZE],
            'original_archive': preserved, 'original_ranking_frozen_at_unchanged': old_freeze['frozen_at'],
            'restored_native_columns': 487, 'all_other_56_columns_unchanged': True,
            'row_order_identity_rank_grade_unchanged': True, 'no_reranking_or_feature_recalculation': True,
            'restoration_validation': validation}
        need(all(new_freeze[k] == v for k, v in old_freeze.items() if k != 'ranking_table_sha256'), 'Unrelated freeze metadata changed')
        stage_freeze = staged / 'RESTORED_FREEZE.json'; dump(stage_freeze, new_freeze)
        need(js(stage_freeze) == new_freeze and new_freeze['frozen_at'] == old_freeze['frozen_at'], 'Staged freeze failed roundtrip')
        state.update(status='ALL_STAGED_STRING_AUDITS_PASS_BEFORE_PROMOTION', archive=preserved, validation=validation,
                     changed_native_cells=sum(counts.values()), restored_at=restored_at, new_rank_sha256=new_hash,
                     new_freeze_sha256=digest(stage_freeze), original_frozen_at_unchanged=old_freeze['frozen_at'])
        dump(state_path, state)
        unchanged(root, protected)
        for rel in ('results/integrated/CANDIDATE_FREEZE.json', 'manifest/CANDIDATE_FREEZE_GATE.json'):
            need(not path(root, rel).exists(), 'Candidate decision appeared during audit; refuse promotion')
        # Two files cannot be atomically replaced together. Each replace is
        # atomic; a crash between them is fail-closed, guarded and documented.
        os.replace(stage_rank, path(root, RANK))
        state['status'] = 'RANK_PROMOTED_FREEZE_PENDING'; dump(state_path, state)
        os.replace(stage_freeze, path(root, FREEZE))
        need(digest(path(root, RANK)) == new_hash and digest(path(root, FREEZE)) == state['new_freeze_sha256'], 'Canonical post-promotion SHA mismatch')
        unchanged(root, protected, exclude=(RANK, FREEZE))
        state.update(status='PASS_NATIVE_STRINGS_RESTORED_FROZEN_ORDER_UNCHANGED', finished_at=utc(),
                     duration_seconds=time.monotonic() - tick, helper_did_not_write_execution_trace=True)
        dump(state_path, state)
        (audit / 'RESTORATION_AUDIT_CN.md').write_text(
            '# 冻结排名原列字符串恢复审计\n\n仅恢复宽表487个原列；56个非原列、45289行身份/行序/rank/grade均逐格不变。没有调用排名、统计、模型或候选选择。\n\n'
            + '原冻结时间保留：' + old_freeze['frozen_at'] + '\n\n恢复真实时间：' + restored_at
            + '\n\n新rank SHA256：' + new_hash + '\n\n完整真实命令、起止/时长、原归档SHA和全表验证见 RESTORATION_AUDIT.json。'
            + '\n\nhelper未写EXECUTION_TRACE；外层run_logged正常追加此次文件修复真实START/FINISH，既有历史不改。'
            + '\n\n两文件分别原子替换，不能声称跨文件原子事务。单次guard与全部失败/partial保留；禁止自动重入。\n', encoding='utf-8')
        print(json.dumps({k: state[k] for k in ('status', 'changed_native_cells', 'new_rank_sha256', 'new_freeze_sha256', 'validation')}, ensure_ascii=False))
    except BaseException as exc:
        failure_stage = state['status']
        state.update(status='FAILED_OR_PARTIALLY_PROMOTED_NEEDS_EXPLICIT_REVIEW_NO_RETRY', failure_stage=failure_stage, failed_at=utc(),
                     elapsed_seconds=time.monotonic() - tick, error=repr(exc), traceback=traceback.format_exc())
        dump(state_path, state)
        (audit / 'RESTORATION_AUDIT_CN.md').write_text(
            '# 原列字符串恢复未完成：需要明确审查\n\n失败阶段：' + failure_stage + '\n\n真实开始：' + started
            + '\n\n错误：' + repr(exc) + '\n\n全部原归档/partial/单次guard保留。不要自动重新执行，不要据此声称排名恢复完成。详见 RESTORATION_AUDIT.json。\n', encoding='utf-8')
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--root', required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True); mode.add_argument('--apply', action='store_true'); mode.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    try: csv.field_size_limit(sys.maxsize)
    except OverflowError: csv.field_size_limit(2 ** 31 - 1)
    run(args.root.resolve(), args.apply)


if __name__ == '__main__': main()
