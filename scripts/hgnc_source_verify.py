#!/usr/bin/env python3
"""Authenticate an unchanged HGNC reference against a new official snapshot.

The old reference is opened read-only. This operation never changes a model,
identity mapping, statistical table, SNP/gene inclusion rule or project gate.
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import socket
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
ORIGINAL = "data/reference/HGNC/hgnc_complete_set.txt"
MANIFEST_COLUMNS = ("kind", "path", "source_url", "immutable_source", "bytes", "sha256", "validation", "copied_at_utc")


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for data in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def atomic_json(path, value):
    temp = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def validate_hgnc(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8-sig")
    required = ("hgnc_id", "symbol", "status")
    missing_fields = [c for c in required if c not in df.columns]
    result = {"columns": df.columns.tolist(), "data_rows": len(df), "required_fields": list(required),
              "missing_required_fields": missing_fields}
    if missing_fields:
        result["status"] = "SCHEMA_FAIL"
        return result
    approved = df[df.status.eq("Approved")]
    result.update(status_counts={str(k): int(v) for k, v in df.status.value_counts(dropna=False).items()},
                  approved_rows=len(approved), all_hgnc_id_nonempty=bool(df.hgnc_id.ne("").all()),
                  all_hgnc_id_format_valid=bool(df.hgnc_id.str.fullmatch(r"HGNC:[0-9]+").all()),
                  all_hgnc_ids_unique=bool(df.hgnc_id.is_unique),
                  approved_hgnc_ids_unique=bool(approved.hgnc_id.is_unique),
                  approved_symbols_unique_case_sensitive=bool(approved.symbol.is_unique),
                  approved_symbols_nonempty=bool(approved.symbol.ne("").all()),
                  approved_duplicate_hgnc_id_rows=int(approved.hgnc_id.duplicated(keep=False).sum()),
                  approved_duplicate_symbol_rows=int(approved.symbol.duplicated(keep=False).sum()))
    result["status"] = ("APPROVED_IDENTITIES_PASS" if len(approved) > 0 and
                        all(result[c] for c in ("all_hgnc_id_nonempty", "all_hgnc_id_format_valid", "all_hgnc_ids_unique",
                            "approved_hgnc_ids_unique", "approved_symbols_unique_case_sensitive", "approved_symbols_nonempty"))
                        else "IDENTITY_VALIDATION_FAIL")
    if result["status"] != "APPROVED_IDENTITIES_PASS":
        duplicate = approved.hgnc_id.duplicated(keep=False) | approved.symbol.duplicated(keep=False)
        result["duplicate_examples_first20"] = approved.loc[duplicate, list(required)].head(20).to_dict("records")
    return result


def append_authentication(root, original, snapshot, original_sha, original_bytes):
    import fcntl
    path = root / "manifest/DATA_PROVENANCE.tsv"
    if not path.exists():
        raise FileNotFoundError("Existing provenance manifest required; refusing to create or replace it")
    row = {"kind": "OFFICIAL_REFERENCE_AUTHENTICATED_BY_REFETCH", "path": original.relative_to(root).as_posix(),
           "source_url": URL, "immutable_source": str(snapshot), "bytes": str(original_bytes), "sha256": original_sha,
           "validation": "authenticated_by_same_SHA_new_official_refetch;original_retrieval_time_not_recorded_NA;see_data/reference/HGNC/SOURCE_VERIFICATION.json",
           "copied_at_utc": "NA"}
    with path.open("a+", encoding="utf-8", newline="") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.seek(0)
        rows = list(csv.reader(f, delimiter="\t"))
        if not rows:
            raise ValueError("Existing provenance manifest has no header")
        header = rows[0]
        if header:
            header[0] = header[0].lstrip("\ufeff")
        if tuple(header) != MANIFEST_COLUMNS:
            raise ValueError("Actual provenance columns differ from the supplied existing schema")
        values = [row[c] for c in header]
        if values in rows[1:]:
            result = "ALREADY_PRESENT_EQUIVALENT_AUTH_ROW"
        else:
            f.seek(0, 2)
            csv.writer(f, delimiter="\t", lineterminator="\n").writerow(values)
            f.flush()
            os.fsync(f.fileno())
            result = "ONE_AUTHENTICATION_ROW_APPENDED_NO_EXISTING_ROWS_REWRITTEN"
        fcntl.flock(f, fcntl.LOCK_UN)
    return {"path": path.relative_to(root).as_posix(), "append_status": result, "appended_values": row,
            "columns_unchanged": True, "original_copied_at_representation": "copied_at_utc=NA"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--curl", default="/home/asset/miniconda3/bin/curl")
    args = ap.parse_args()
    root = args.root.resolve()
    if root.name != "IBS_A_OFFICIAL_RERUN_20261105":
        raise ValueError("New clean-rerun root required")
    original = root / ORIGINAL
    directory = original.parent
    snapshot = directory / "official_refetch_20261105.txt"
    part = snapshot.with_suffix(snapshot.suffix + ".part")
    headers = directory / "official_refetch_20261105.http_headers.txt"
    stderr_path = directory / "official_refetch_20261105.curl_stderr.txt"
    status_path = directory / "SOURCE_VERIFICATION.json"
    if snapshot.exists() or part.exists():
        raise FileExistsError("Preserve existing refetch/partial snapshot; audit must not silently replace it")
    original_sha, original_bytes = sha(original), original.stat().st_size
    curl = args.curl if Path(args.curl).is_file() else shutil.which("curl")
    if not curl:
        raise FileNotFoundError("Official HTTPS curl client is unavailable")
    curl_version = subprocess.check_output([curl, "--version"], text=True).splitlines()[0]
    command = [curl, "--location", "--fail", "--show-error", "--silent", "--noproxy", "*",
               "--retry", "1", "--retry-delay", "2", "--connect-timeout", "15", "--max-time", "180",
               "--dump-header", str(headers), "--output", str(part),
               "--write-out", "%{http_code}\t%{url_effective}\t%{size_download}\t%{content_type}\t%{time_total}\n", URL]
    start, mono = now(), time.monotonic()
    process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    finish, duration = now(), time.monotonic() - mono
    stderr_path.write_text(process.stderr, encoding="utf-8")
    record = {"status": "FETCH_FAILED", "scope": "official_raw_reference_source_authentication_only_not_a_model",
              "is_new_fit": False, "source_url": URL,
              "original_reference": {"path": ORIGINAL, "bytes": original_bytes, "sha256": original_sha,
                  "initial_retrieved_at": "NA", "original_copied_at": "NA",
                  "initial_time_missing_reason": "The first genuine curl download start/finish was not independently logged; no timestamp is inferred from file mtime or later operations."},
              "refetch": {"requested_url": URL, "start_iso8601": start, "finish_iso8601": finish,
                  "duration_seconds": duration, "curl_exit_code": process.returncode, "curl_version": curl_version,
                  "command_argv": command, "curl_writeout": process.stdout.strip(),
                  "HTTP_headers_path": headers.relative_to(root).as_posix(),
                  "HTTP_headers_sha256": sha(headers) if headers.exists() else None,
                  "stderr_path": stderr_path.relative_to(root).as_posix(), "stderr_sha256": sha(stderr_path)},
              "hostname": socket.gethostname(), "software": {"Python": platform.python_version(), "pandas": pd.__version__},
              "script_sha256": sha(Path(__file__)), "any_model_mapping_statistic_or_gate_changed": False,
              "manifest_authentication_row_appended": False}
    if process.returncode:
        record["original_reference_unchanged"] = sha(original) == original_sha
        atomic_json(status_path, record)
        raise SystemExit(process.returncode)
    fields = process.stdout.strip().splitlines()[-1].split("\t")
    if len(fields) != 5 or fields[0] != "200":
        record["status"] = "UNEXPECTED_HTTP_RESULT"
        record["original_reference_unchanged"] = sha(original) == original_sha
        atomic_json(status_path, record)
        raise RuntimeError("Expected a full HTTP200 official snapshot")
    part.replace(snapshot)
    snapshot_sha, snapshot_bytes = sha(snapshot), snapshot.stat().st_size
    record["refetch"].update(HTTP_status=int(fields[0]), effective_url=fields[1], size_download_reported=fields[2],
            content_type=fields[3], curl_time_total_seconds=float(fields[4]),
            snapshot_path=snapshot.relative_to(root).as_posix(), bytes=snapshot_bytes, sha256=snapshot_sha)
    record["original_reference_unchanged"] = sha(original) == original_sha
    record["same_SHA256"] = snapshot_sha == original_sha
    record["original_reference_schema"] = validate_hgnc(original)
    record["current_snapshot_schema"] = validate_hgnc(snapshot)
    if not record["original_reference_unchanged"]:
        record["status"] = "ORIGINAL_REFERENCE_CHANGED_EXTERNALLY_DURING_VERIFICATION"
        atomic_json(status_path, record)
        raise RuntimeError("Original changed during verification; this script never writes it")
    schemas_pass = all(record[k]["status"] == "APPROVED_IDENTITIES_PASS" for k in ("original_reference_schema", "current_snapshot_schema"))
    if not record["same_SHA256"]:
        record["status"] = "DIFFERENT_CURRENT_SNAPSHOT"
        record["interpretation"] = "Both snapshots retained separately. This temporal source difference does not change the original reference, mappings, statistics, gates or model-completion status."
    elif not schemas_pass:
        record["status"] = "SOURCE_MATCHED_BUT_REFERENCE_SCHEMA_REQUIRES_REVIEW"
    else:
        record["status"] = "OFFICIAL_REFERENCE_AUTHENTICATED_BY_REFETCH"
        # Commit source evidence before linking it from the append-only manifest.
        atomic_json(status_path, record)
        record["manifest"] = append_authentication(root, original, snapshot, original_sha, original_bytes)
        record["manifest_authentication_row_appended"] = record["manifest"]["append_status"].startswith("ONE_")
    record["finished_at"] = now()
    record["output_hashes_sha256"] = {p.relative_to(root).as_posix(): sha(p) for p in (snapshot, headers, stderr_path)}
    record["JSON_self_hash"] = "Recorded by run_logged FINISH output_hashes and post-write stdout; not recursively embedded"
    atomic_json(status_path, record)
    print(json.dumps({"status": record["status"], "original_SHA256": original_sha,
        "snapshot_SHA256": snapshot_sha, "original_unchanged": record["original_reference_unchanged"],
        "original_bytes": original_bytes, "snapshot_bytes": snapshot_bytes,
        "approved_rows": record["original_reference_schema"].get("approved_rows"),
        "refetch_start": start, "refetch_finish": finish, "refetch_seconds": duration,
        "SOURCE_VERIFICATION_JSON_SHA256": sha(status_path),
        "manifest_authentication_row_appended": record["manifest_authentication_row_appended"]}), flush=True)


if __name__ == "__main__":
    main()
