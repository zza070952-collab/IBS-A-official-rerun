#!/usr/bin/env python3
"""Independent raw GWAS / MAGMA-reference audit for the official IBS rerun.

No association model is fitted here. The on-disk SQLite database bounds RAM and
allows duplicate SNP identifiers to be removed in their entirety, without
choosing the most significant row. All paths are confined to this new project.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import sys
import time
import zlib


EXPECTED_PROJECT = "IBS_A_OFFICIAL_RERUN_20261105"
GWAS_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016564/GCST90016564_buildGRCh37.tsv"
MAGMA_PAGE = "https://cncr.nl/research/magma/"
GENE_ARCHIVE_URL = "https://vu.data.surf.nl/index.php/s/Pj2orwuF2JYyKxq/download"
OFFICIAL_GENE_SHA256 = "bfcb5273801cc3e37dcb98f5465f4498daa8c0cfcda2f4f3fbc51abe02000e7e"


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def confined(root, relative, must_exist=True):
    path = root / relative
    resolved = path.resolve()
    if not resolved.is_relative_to(root) or path.is_symlink():
        raise ValueError(f"Input/output must be a real file inside new root: {path}")
    if must_exist and (not path.is_file() or path.stat().st_size == 0):
        raise ValueError(f"Missing or empty file: {path}")
    return path


def dump_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--gwas", default="data/raw/gwas/GCST90016564_buildGRCh37.tsv")
    parser.add_argument("--bfile", default="data/reference/gwas/g1000_eur")
    parser.add_argument("--gene-loc", default="data/reference/gwas/NCBI37.3.gene.loc")
    parser.add_argument("--magma", default="env/tools/magma/magma")
    parser.add_argument("--chunksize", type=int, default=100000)
    parser.add_argument("--manhattan-modulus", type=int, default=1000,
                        help="CRC32(SNP) modulo sampling; retain every P<5e-8 row regardless")
    return parser.parse_args()


def main():
    args = get_args()
    root = args.root.resolve()
    if root.name != EXPECTED_PROJECT:
        raise ValueError(f"Expected new project root named {EXPECTED_PROJECT}, received {root}")
    if args.chunksize < 1000 or args.manhattan_modulus < 1:
        raise ValueError("Invalid chunk/sampling size")
    import numpy as np
    import pandas as pd
    from scipy.stats import chi2

    started, monotonic = now(), time.monotonic()
    outdir = root / "results/gwas/qc"
    outdir.mkdir(parents=True, exist_ok=True)
    (root / "manifest").mkdir(exist_ok=True)
    gate = root / "manifest/GWAS_PHASE0_GATE.json"
    # Never accept a stale success marker during a failed replacement audit.
    dump_json(gate, {"status": "RUNNING", "start_iso8601": started})
    database = outdir / "phase0.sqlite"
    if database.exists():
        raise FileExistsError("Phase0 database already exists; preserve it and use a new explicitly audited attempt directory, not silent overwrite")
    paths = {"raw_gwas": confined(root, args.gwas),
             "gene_locations": confined(root, args.gene_loc),
             "magma_executable": confined(root, args.magma)}
    paths.update({"ld_" + ext: confined(root, args.bfile + "." + ext) for ext in ("bed", "bim", "fam")})
    inputs = {k: {"path": str(v.relative_to(root)), "sha256": sha256(v), "bytes": v.stat().st_size} for k, v in paths.items()}
    # This exact archive was downloaded anew from the migrated official page;
    # six columns alone are insufficient to authenticate a reference file.
    if inputs["gene_locations"]["sha256"] != OFFICIAL_GENE_SHA256:
        raise ValueError("GRCh37 gene-location SHA256 differs from independently verified official archive")
    genes = pd.read_csv(paths["gene_locations"], sep="\t", header=None, dtype=str)
    if genes.shape[1] != 6 or genes[0].duplicated().any():
        raise ValueError("Invalid or duplicate gene locations")
    gstart = pd.to_numeric(genes[2], errors="coerce")
    gend = pd.to_numeric(genes[3], errors="coerce")
    if (gstart.isna() | gend.isna() | (gstart < 1) | (gend < gstart)).any():
        raise ValueError("Invalid gene intervals")
    if not genes[1].isin([str(i) for i in range(1, 23)] + ["X", "Y"]).all():
        raise ValueError("Unexpected gene chromosome")
    with paths["ld_bed"].open("rb") as handle:
        if handle.read(3) != b"\x6c\x1b\x01":
            raise ValueError("LD BED is not SNP-major PLINK bed")
    n_fam = sum(1 for _ in paths["ld_fam"].open())
    connection = sqlite3.connect(database)
    connection.executescript("""
      PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL; PRAGMA temp_store=FILE;
      PRAGMA cache_size=-131072;
      CREATE TABLE ld (SNP TEXT PRIMARY KEY, CHR INTEGER, BP INTEGER, seen INTEGER NOT NULL);
      CREATE TABLE raw (SNP TEXT PRIMARY KEY, CHR INTEGER, BP INTEGER, P REAL, N INTEGER,
                        eligible INTEGER NOT NULL, seen INTEGER NOT NULL);
    """)
    reference_rows = 0
    for chunk in pd.read_csv(paths["ld_bim"], sep=r"\s+", header=None,
                             names=["CHR", "SNP", "CM", "BP", "A1", "A2"], dtype=str, chunksize=args.chunksize):
        chrn = pd.to_numeric(chunk.CHR.replace({"X": "23", "Y": "24", "MT": "26"}), errors="coerce")
        bp = pd.to_numeric(chunk.BP, errors="coerce")
        if (chrn.isna() | bp.isna() | (bp < 1) | (bp % 1 != 0)).any() or chunk.SNP.isna().any():
            raise ValueError("Malformed PLINK BIM")
        rows = zip(chunk.SNP, chrn.astype(int), bp.astype(int))
        connection.executemany("INSERT INTO ld VALUES(?,?,?,1) ON CONFLICT(SNP) DO UPDATE SET seen=seen+1", rows)
        connection.commit()
        reference_rows += len(chunk)
    expected_bed_bytes = 3 + reference_rows * ((n_fam + 3) // 4)
    if paths["ld_bed"].stat().st_size != expected_bed_bytes:
        raise ValueError("BED/BIM/FAM dimensions disagree")
    raw_header = list(pd.read_csv(paths["raw_gwas"], sep="\t", nrows=0).columns)
    aliases = {
        "SNP": ["variant_id", "rsid", "rs_id", "SNP"],
        "CHR": ["chromosome", "CHR"], "BP": ["base_pair_location", "BP"],
        "P": ["p_value", "P"], "N_CASE": ["N_CASE", "n_cases", "n_case"],
        "N_CONTROL": ["N_CONTROL", "n_controls", "n_control"],
    }
    columns = {}
    for canonical, options in aliases.items():
        found = [name for name in options if name in raw_header]
        if len(found) != 1:
            raise ValueError(f"Expected exactly one {canonical} column, found {found}; raw header={raw_header}")
        columns[canonical] = found[0]
    count = Counter()
    missing = Counter()
    limits = {}
    chr_counts = Counter()
    example_path = outdir / "invalid_row_examples.tsv"
    n_examples = 0
    with example_path.open("w", newline="", encoding="utf-8") as ex:
        writer = csv.writer(ex, delimiter="\t")
        writer.writerow(["raw_row_number", "SNP", "reason"])
        for chunk in pd.read_csv(paths["raw_gwas"], sep="\t", dtype=str, chunksize=args.chunksize):
            offset = count["raw_rows"]
            count["raw_rows"] += len(chunk)
            for col in raw_header:
                missing[col] += int(chunk[col].isna().sum() + chunk[col].eq("").sum())
            values = {}
            for field in ("CHR", "BP", "P", "N_CASE", "N_CONTROL"):
                text_col = chunk[columns[field]]
                if field == "CHR":
                    text_col = text_col.str.replace(r"^chr", "", regex=True).replace({"X": "23", "Y": "24", "MT": "26", "M": "26"})
                values[field] = pd.to_numeric(text_col, errors="coerce")
                good = values[field].dropna()
                if len(good):
                    low, high = float(good.min()), float(good.max())
                    previous = limits.get(field, [low, high])
                    limits[field] = [min(previous[0], low), max(previous[1], high)]
            sid = chunk[columns["SNP"]]
            bad = {
                "missing_snp": sid.isna() | sid.eq("") | sid.eq("."),
                "nonautosomal_or_invalid_chr": ~values["CHR"].isin(range(1, 23)),
                "invalid_bp": values["BP"].isna() | (values["BP"] < 1) | (values["BP"] % 1 != 0),
                "invalid_p": values["P"].isna() | ~np.isfinite(values["P"]) | (values["P"] <= 0) | (values["P"] > 1),
                "invalid_n_case": values["N_CASE"].isna() | ~np.isfinite(values["N_CASE"]) | (values["N_CASE"] <= 0) | (values["N_CASE"] % 1 != 0),
                "invalid_n_control": values["N_CONTROL"].isna() | ~np.isfinite(values["N_CONTROL"]) | (values["N_CONTROL"] <= 0) | (values["N_CONTROL"] % 1 != 0),
            }
            eligible = ~pd.concat(list(bad.values()), axis=1).any(axis=1)
            count["p_equal_zero"] += int(values["P"].eq(0).sum())
            count["eligible_raw_rows_before_dedup"] += int(eligible.sum())
            for reason, mask in bad.items():
                count[reason] += int(mask.sum())
            chr_counts.update(chunk[columns["CHR"]].fillna("MISSING").value_counts().to_dict())
            for i in chunk.index[~eligible][:max(0, 1000 - n_examples)]:
                writer.writerow([offset + int(chunk.index.get_loc(i)) + 2, sid.loc[i], ";".join(k for k, v in bad.items() if bool(v.loc[i]))])
                n_examples += 1
            values["N"] = values["N_CASE"] + values["N_CONTROL"]
            if eligible.any():
                nv = values["N"][eligible]
                previous = limits.get("N_TOTAL", [float(nv.min()), float(nv.max())])
                limits["N_TOTAL"] = [min(previous[0], float(nv.min())), max(previous[1], float(nv.max()))]
            # Store identifiers of invalid rows too: a repeated identifier is
            # excluded altogether even if only one occurrence has valid P/N.
            records = ((str(s), int(c) if ok else None, int(b) if ok else None,
                        float(p) if ok else None, int(n) if ok else None, int(ok))
                       for s, c, b, p, n, ok in zip(sid.fillna(""), values["CHR"], values["BP"], values["P"], values["N"], eligible)
                       if str(s) not in ("", "."))
            connection.executemany("INSERT INTO raw VALUES(?,?,?,?,?,?,1) ON CONFLICT(SNP) DO UPDATE SET seen=seen+1", records)
            connection.commit()
            print(json.dumps({"stage": "raw_qc", "rows_read": count["raw_rows"], "time": now()}), flush=True)
    scalar = lambda sql: int(connection.execute(sql).fetchone()[0])
    count["duplicate_raw_identifiers"] = scalar("SELECT COUNT(*) FROM raw WHERE seen>1")
    count["rows_with_duplicate_raw_identifiers"] = scalar("SELECT COALESCE(SUM(seen),0) FROM raw WHERE seen>1")
    count["duplicate_ld_identifiers"] = scalar("SELECT COUNT(*) FROM ld WHERE seen>1")
    count["valid_unique_autosomal_snps"] = scalar("SELECT COUNT(*) FROM raw WHERE eligible=1 AND seen=1")
    base = "FROM raw r JOIN ld l ON r.SNP=l.SNP WHERE r.eligible=1 AND r.seen=1 AND l.seen=1"
    count["coordinate_conflicts_with_ld"] = scalar("SELECT COUNT(*) " + base + " AND (r.CHR!=l.CHR OR r.BP!=l.BP)")
    count["eligible_exact_ld_matches"] = scalar("SELECT COUNT(*) " + base + " AND r.CHR=l.CHR AND r.BP=l.BP")
    count["valid_absent_from_ld"] = scalar("SELECT COUNT(*) FROM raw r LEFT JOIN ld l ON r.SNP=l.SNP WHERE r.eligible=1 AND r.seen=1 AND l.SNP IS NULL")
    if count["eligible_exact_ld_matches"] == 0:
        raise ValueError("No eligible exact SNP/coordinate overlap with EUR LD")
    score_file = outdir / "magma_snps_p_n.tsv"
    snploc_file = outdir / "magma_snps.loc"
    with score_file.open("w", newline="") as pf, snploc_file.open("w", newline="") as lf:
        pw, lw = csv.writer(pf, delimiter="\t"), csv.writer(lf, delimiter="\t")
        pw.writerow(["SNP", "P", "N"])
        for snp, chrom, bp, p, n in connection.execute("SELECT r.SNP,r.CHR,r.BP,r.P,r.N " + base + " AND r.CHR=l.CHR AND r.BP=l.BP ORDER BY r.CHR,r.BP,r.SNP"):
            pw.writerow([snp, format(p, ".17g"), n])
            lw.writerow([snp, chrom, bp])
    # The Manhattan/QQ universe is the full validated unique autosomal GWAS,
    # not merely the subset available in the LD panel or gene annotation.
    plot_file = outdir / "manhattan_points.tsv.gz"
    retained_plot = 0
    with gzip.open(plot_file, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["SNP", "CHR", "BP", "P", "sampling_reason"])
        for snp, chrom, bp, p in connection.execute("SELECT SNP,CHR,BP,P FROM raw WHERE eligible=1 AND seen=1"):
            significant = p < 5e-8
            if significant or zlib.crc32(snp.encode()) % args.manhattan_modulus == 0:
                writer.writerow([snp, chrom, bp, format(p, ".17g"), "all_P_lt_5e-8" if significant else "deterministic_CRC32_modulo"])
                retained_plot += 1
    n_qq = count["valid_unique_autosomal_snps"]
    positions = set(np.linspace(0, n_qq - 1, min(10000, n_qq), dtype=np.int64).tolist())
    positions.update(range(min(1000, n_qq)))
    midpoint = {n_qq // 2} if n_qq % 2 else {n_qq // 2 - 1, n_qq // 2}
    median_values = []
    with (outdir / "qq_quantiles.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "total_tests", "expected_P", "observed_P"])
        for i, (p,) in enumerate(connection.execute("SELECT P FROM raw WHERE eligible=1 AND seen=1 ORDER BY P")):
            if i in midpoint:
                median_values.append(p)
            if i in positions:
                writer.writerow([i + 1, n_qq, (i + 0.5) / n_qq, format(p, ".17g")])
    median_p = float(np.mean(median_values))
    lambda_gc = float(chi2.isf(median_p, 1) / chi2.ppf(0.5, 1))
    connection.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    connection.close()
    summary = {
        "status": "PASS", "start_iso8601": started, "finish_iso8601": now(),
        "duration_seconds": time.monotonic() - monotonic, "source_url": GWAS_URL,
        "genome_build": "GRCh37", "raw_columns": raw_header, "resolved_columns": columns,
        "missing_by_column": dict(missing), "numeric_ranges": limits,
        "raw_chromosome_counts": dict(chr_counts), "counts": dict(count),
        "reference": {"magma_official_page": MAGMA_PAGE, "gene_archive_url": GENE_ARCHIVE_URL,
                      "gene_rows": len(genes), "gene_symbol_version": "official historical NCBI gene_info, not current HGNC",
                      "ld_rows": reference_rows, "ld_samples": n_fam, "inputs": inputs},
        "policies": {"primary_chromosomes": "autosomes 1-22", "duplicate_ids": "exclude every occurrence, including identical duplicates",
                     "invalid_p": "exclude; never replace P=0", "N": "per-variant N_CASE + N_CONTROL, not effective N",
                     "ld_coordinate_disagreement": "exclude, never alter coordinates", "alleles": "P-only MAGMA does not require signed effect alignment; raw allele columns preserved in source",
                     "gene_window": "0 kb upstream and downstream; gene body only (analyst pre-specified, not a claim of official recommendation)"},
        "plotting": {"universe": "all valid unique autosomal raw GWAS SNPs", "manhattan_modulus": args.manhattan_modulus,
                     "background_sampling": "CRC32(SNP) mod modulus == 0; every P<5e-8 retained", "points_written": retained_plot,
                     "qq_positions": "first 1000 ranks plus up to 10000 evenly spaced order statistics", "lambda_gc": lambda_gc,
                     "lambda_gc_interpretation": "descriptive; case-control meta-analysis polygenicity can affect inflation; not a calibration correction"},
    }
    dump_json(outdir / "summary.json", summary)
    with (root / "manifest/GWAS_INPUT_MANIFEST.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["input", "relative_path", "bytes", "sha256", "genome_build"])
        for name, obj in inputs.items():
            writer.writerow([name, obj["path"], obj["bytes"], obj["sha256"], "GRCh37" if name != "magma_executable" else "NA"])
    outputs = {str(p.relative_to(root)): sha256(p) for p in [score_file, snploc_file, outdir / "summary.json", plot_file, outdir / "qq_quantiles.tsv"]}
    dump_json(gate, {"status": "PASS", "phase": "raw/reference audit only; no association fit",
                     "finish_iso8601": now(), "input_sha256": inputs, "output_sha256": outputs,
                     "eligible_snps_for_magma": count["eligible_exact_ld_matches"]})
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"GWAS_PHASE0_FAIL: {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        raise
