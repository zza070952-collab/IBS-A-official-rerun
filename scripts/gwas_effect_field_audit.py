#!/usr/bin/env python3
"""Independent descriptive audit of raw GWAS effect fields; never filters SNPs.

Reads every input row once in chunks. P/N, duplicate keys, LD matching, effect
calibration and association tests are deliberately outside this audit's scope.
"""
import argparse
import csv
import hashlib
import io
import json
import os
import platform
import socket
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

NUMERIC = ("beta", "standard_error", "effect_allele_frequency")
ALLELES = ("effect_allele", "other_allele")
REQUIRED = ("variant_id", "p_value", "chromosome", "base_pair_location", *ALLELES,
            "effect_allele_frequency", "beta", "standard_error", "N_CASE", "N_CONTROL")
READ = ("variant_id", "chromosome", "base_pair_location", *ALLELES, *NUMERIC)
MISSING = {"", ".", "na", "n/a", "null", "none", "nan", "+nan", "-nan"}
MAX_EXAMPLES = 1000


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for data in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


class DigestingReader(io.RawIOBase):
    def __init__(self, raw):
        self.raw, self.hasher, self.bytes_read = raw, hashlib.sha256(), 0

    def readable(self):
        return True

    def readinto(self, buffer):
        data = self.raw.read(len(buffer))
        n = len(data)
        buffer[:n] = data
        self.hasher.update(data)
        self.bytes_read += n
        return n


def numeric_chunk(series, field):
    text = series.astype(str).str.strip()
    missing = text.str.lower().isin(MISSING).to_numpy()
    values = pd.to_numeric(text, errors="coerce").to_numpy(dtype=np.float64, copy=True)
    nonnumeric = np.isnan(values) & ~missing
    # pandas may coerce numeric overflow (e.g. 1e999) to NaN. Distinguish it
    # from genuinely nonnumeric text with Python's float parser on exceptions.
    if nonnumeric.any():
        for token in pd.unique(text[nonnumeric]):
            mask = (text == token).to_numpy() & nonnumeric
            try:
                value = float(token)
            except (ValueError, OverflowError):
                continue
            values[mask] = value
            nonnumeric[mask] = False
    finite = np.isfinite(values) & ~missing & ~nonnumeric
    nonfinite = ~np.isfinite(values) & ~missing & ~nonnumeric
    counts = {"rows": len(series), "missing_token": int(missing.sum()),
              "nonnumeric_nonmissing": int(nonnumeric.sum()), "nonfinite_nonmissing": int(nonfinite.sum()),
              "positive_infinity": int((nonfinite & np.isposinf(values)).sum()),
              "negative_infinity": int((nonfinite & np.isneginf(values)).sum()),
              "other_nonfinite_nan": int((nonfinite & np.isnan(values)).sum()),
              "finite": int(finite.sum()), "finite_negative": int((finite & (values < 0)).sum()),
              "finite_zero": int((finite & (values == 0)).sum()), "finite_positive": int((finite & (values > 0)).sum())}
    assert counts["missing_token"] + counts["nonnumeric_nonmissing"] + counts["nonfinite_nonmissing"] + counts["finite"] == len(series)
    flags = {field + ":missing_token": missing, field + ":nonnumeric": nonnumeric,
             field + ":nonfinite_nonmissing": nonfinite}
    if field == "standard_error":
        nonpositive = finite & (values <= 0)
        counts["finite_nonpositive"] = int(nonpositive.sum())
        flags[field + ":finite_nonpositive"] = nonpositive
    if field == "effect_allele_frequency":
        for name, mask in (("finite_below0", finite & (values < 0)), ("finite_equal0", finite & (values == 0)),
                           ("finite_between0and1", finite & (values > 0) & (values < 1)),
                           ("finite_equal1", finite & (values == 1)), ("finite_above1", finite & (values > 1))):
            counts[name] = int(mask.sum())
        flags[field + ":finite_outside_closed_0_1"] = finite & ((values < 0) | (values > 1))
    limits = (float(values[finite].min()), float(values[finite].max())) if finite.any() else (None, None)
    return counts, limits, flags


def allele_chunk(series, field):
    raw = series.astype(str)
    text = raw.str.strip().str.upper()
    missing = text.str.lower().isin(MISSING).to_numpy()
    single = text.str.fullmatch("[ACGT]", na=False).to_numpy() & ~missing
    dna = text.str.fullmatch("[ACGT]+", na=False).to_numpy() & ~missing
    multi = dna & ~single
    symbolic = (text.isin(["I", "D", "+", "-", "INS", "DEL", "*"]) |
                text.str.fullmatch("<[^<>]+>", na=False)).to_numpy() & ~missing
    ambiguous = text.str.fullmatch("[ACGTRYSWKMBDHVN]+", na=False).to_numpy() & ~missing & ~dna & ~symbolic
    other = ~(missing | single | multi | symbolic | ambiguous)
    masks = {"missing_token": missing, "single_ACGT": single, "multi_base_ACGT_sequence": multi,
             "symbolic_or_I_D_token": symbolic, "ambiguous_IUPAC_sequence": ambiguous,
             "other_string": other}
    counts = {k: int(v.sum()) for k, v in masks.items()}
    counts.update({"rows": len(series), "surrounding_whitespace": int(raw.ne(raw.str.strip()).sum()),
                   "case_changed_by_classification_only": int(raw.str.strip().ne(text).sum())})
    for base in "ACGT":
        counts["single_" + base] = int(text.eq(base).sum())
    assert sum(counts[k] for k in masks) == len(series)
    flags = {field + ":missing_token": missing, field + ":ambiguous_IUPAC_sequence": ambiguous,
             field + ":other_string": other}
    # Symbolic tokens and multibase DNA are catalogued, not labelled invalid.
    return text, masks, counts, flags


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--input", default="data/raw/gwas/GCST90016564_buildGRCh37.tsv")
    ap.add_argument("--chunksize", type=int, default=100000)
    args = ap.parse_args()
    root = args.root.resolve()
    if root.name != "IBS_A_OFFICIAL_RERUN_20261105" or args.chunksize < 1:
        raise ValueError("Only the new clean project and a positive chunk size are permitted")
    source = (root / args.input).resolve()
    if root not in source.parents:
        raise ValueError("Input escaped new root")
    output = root / "results/gwas/qc"
    output.mkdir(parents=True, exist_ok=True)
    started, mono = now(), time.monotonic()
    before = source.stat()
    with source.open(encoding="utf-8-sig", newline="") as f:
        columns = next(csv.reader(f, delimiter="\t"))
    if len(columns) != len(set(columns)) or not set(REQUIRED).issubset(columns):
        raise ValueError("Original required header absent or duplicated")
    numeric = {f: {"counts": Counter(), "finite_min": None, "finite_max": None} for f in NUMERIC}
    alleles = {f: Counter() for f in ALLELES}
    pairs, examples, total, flagged_rows = Counter(), [], 0, 0
    with source.open("rb", buffering=0) as raw:
        reader = DigestingReader(raw)
        with io.BufferedReader(reader, buffer_size=1024 * 1024) as stream:
            chunks = pd.read_csv(stream, sep="\t", encoding="utf-8-sig", usecols=list(READ),
                                 dtype=str, na_filter=False, keep_default_na=False,
                                 skip_blank_lines=False, chunksize=args.chunksize)
            for chunk in chunks:
                chunk = chunk.fillna("")
                flags = {}
                for field in NUMERIC:
                    counts, limits, field_flags = numeric_chunk(chunk[field], field)
                    numeric[field]["counts"].update(counts)
                    if limits[0] is not None:
                        old = numeric[field]
                        old["finite_min"] = limits[0] if old["finite_min"] is None else min(old["finite_min"], limits[0])
                        old["finite_max"] = limits[1] if old["finite_max"] is None else max(old["finite_max"], limits[1])
                    flags.update(field_flags)
                classified = {}
                for field in ALLELES:
                    text, masks, counts, field_flags = allele_chunk(chunk[field], field)
                    alleles[field].update(counts)
                    classified[field] = (text, masks)
                    flags.update(field_flags)
                e, em = classified[ALLELES[0]]
                o, om = classified[ALLELES[1]]
                both_single = em["single_ACGT"] & om["single_ACGT"]
                both_dna = (em["single_ACGT"] | em["multi_base_ACGT_sequence"]) & (om["single_ACGT"] | om["multi_base_ACGT_sequence"])
                identical = e.eq(o).to_numpy()
                equal_length = e.str.len().eq(o.str.len()).to_numpy()
                missing_pair = em["missing_token"] | om["missing_token"]
                pm = {"both_single_distinct": both_single & ~identical,
                      "both_single_identical": both_single & identical,
                      "other_DNA_different_lengths": both_dna & ~both_single & ~equal_length,
                      "other_DNA_equal_lengths_distinct": both_dna & ~both_single & equal_length & ~identical,
                      "other_DNA_identical": both_dna & ~both_single & identical,
                      "at_least_one_missing_token": missing_pair,
                      "nonmissing_symbolic_ambiguous_or_other_pair": ~both_dna & ~missing_pair}
                assert sum(int(m.sum()) for m in pm.values()) == len(chunk)
                pairs.update({k: int(v.sum()) for k, v in pm.items()})
                flags["allele_pair:identical_DNA_strings"] = both_dna & identical
                flagged = np.logical_or.reduce(list(flags.values()))
                flagged_rows += int(flagged.sum())
                if len(examples) < MAX_EXAMPLES:
                    for i in np.flatnonzero(flagged)[:MAX_EXAMPLES - len(examples)]:
                        row = {"data_record_index_1based": total + int(i) + 1,
                               **{c: str(chunk.iloc[i][c]) for c in READ},
                               "descriptive_flags": ";".join(k for k, m in flags.items() if m[i])}
                        examples.append(row)
                total += len(chunk)
                if total % 1000000 == 0:
                    print(now(), "DESCRIPTIVE_EFFECT_ROWS_READ", total, flush=True)
            # Force EOF to cover any parser buffering in the streaming input hash.
            while stream.read(1024 * 1024):
                pass
        input_sha, input_bytes = reader.hasher.hexdigest(), reader.bytes_read
    after = source.stat()
    if input_bytes != before.st_size or before.st_size != after.st_size or before.st_mtime_ns != after.st_mtime_ns:
        raise RuntimeError("Input file changed or was not fully read")
    if not total:
        raise ValueError("No original data records")
    for field in NUMERIC:
        if numeric[field]["counts"]["rows"] != total:
            raise RuntimeError("Numeric field did not cover every row")
    rows = []
    def add(section, field, metric, value):
        rows.append({"section": section, "field": field, "metric": metric, "value": value})
    add("input", "all", "data_rows", total)
    add("input", "all", "input_bytes_fully_streamed", input_bytes)
    add("input", "all", "input_sha256", input_sha)
    for field, obj in numeric.items():
        for metric, value in obj["counts"].items():
            add("numeric", field, metric, value)
        for metric in ("finite_min", "finite_max"):
            add("numeric", field, metric, obj[metric])
    for field, counts in alleles.items():
        for metric, value in counts.items():
            add("allele_string", field, metric, value)
    for metric, value in pairs.items():
        add("allele_pair", "effect_allele_vs_other_allele", metric, value)
    add("sampling", "all", "flagged_data_rows", flagged_rows)
    add("sampling", "all", "examples_saved", len(examples))
    add("sampling", "all", "examples_truncated", flagged_rows > len(examples))
    tsv = output / "RAW_EFFECT_FIELD_AUDIT.tsv"
    pd.DataFrame(rows).to_csv(tsv, sep="\t", index=False, na_rep="NA")
    sample_path = output / "RAW_EFFECT_FIELD_ANOMALIES_FIRST1000.tsv"
    pd.DataFrame(examples, columns=["data_record_index_1based", *READ, "descriptive_flags"]).to_csv(sample_path, sep="\t", index=False)
    result = {"status": "DESCRIPTIVE_AUDIT_COMPLETE", "is_new_fit": False,
              "scope": "beta_SE_EAF_and_allele_strings_only; no_P_N_duplicate_LD_tests_or_effect_calibration",
              "main_Phase0_gate_modified": False, "any_row_or_SNP_inclusion_changed": False,
              "input_path": str(source), "input_relative_path": source.relative_to(root).as_posix(),
              "input_sha256": input_sha, "input_bytes_fully_streamed": input_bytes, "data_rows": total,
              "original_header": columns, "numeric_fields": numeric, "allele_string_fields": alleles,
              "allele_pair_categories": pairs,
              "classification_definitions": {"missing_tokens_case_insensitive_after_whitespace_strip": sorted(MISSING),
                  "numeric_count_partition": "missing_token+nonnumeric_nonmissing+nonfinite_nonmissing+finite=all_rows",
                  "finite_ranges": "float64-parsed observed finite values only; infinities/missing/nonnumeric counted separately",
                  "SE_flag": "finite_SE<=0; descriptive_only_not_a_filter", "EAF_flag": "finite_EAF<0_or>1; endpoints0and1_permitted",
                  "allele_classification": "uppercase/strip copies only; original input never changed; symbolic/multibase tokens not automatically invalid",
                  "sampling_flags": "missing/nonfinite/nonnumeric, finite SE/EAF range flags, ambiguous/other allele strings, identical DNA allele strings",
                  "no_biological_inference": "string categories do not infer variant normalization, strand correctness or effect validity"},
              "flagged_data_rows": flagged_rows, "anomaly_examples_saved": len(examples),
              "anomaly_example_limit": MAX_EXAMPLES, "anomaly_examples_truncated": flagged_rows > len(examples),
              "anomaly_selection": "first_flagged_data_records_in_original_order; at_most1000; unselected_rows_still_counted",
              "software": {"Python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__},
              "hostname": socket.gethostname(), "script_sha256": sha(Path(__file__)),
              "started_at": started, "finished_at": now(), "duration_seconds": time.monotonic() - mono,
              "output_hashes_sha256": {p.relative_to(root).as_posix(): sha(p) for p in (tsv, sample_path)},
              "JSON_self_hash": "recorded_after_write_in_stdout_and_new_root_EXECUTION_TRACE.tsv_FINISH_output_hashes"}
    jp = output / "RAW_EFFECT_FIELD_AUDIT.json"
    temp = jp.with_suffix(".json.tmp")
    temp.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(jp)
    print(json.dumps({"status": result["status"], "data_rows": total, "exit_code": 0,
                      "input_sha256": input_sha, "output_sha256": {p.relative_to(root).as_posix(): sha(p) for p in (jp, tsv, sample_path)},
                      "flagged_data_rows": flagged_rows, "anomaly_examples_saved": len(examples),
                      "anomaly_examples_truncated": result["anomaly_examples_truncated"]}), flush=True)


if __name__ == "__main__":
    main()
