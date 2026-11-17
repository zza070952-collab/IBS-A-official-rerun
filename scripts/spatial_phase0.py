#!/usr/bin/env python3
"""Fresh, candidate-blind audit of the official Harnik raw Visium archive.

Only the named new project and official source snapshot are accepted.  No model,
QC-selected cohort, old project, or candidate table is an input to this program.
The raw-derived sparse cache is created here so Phase 1 need not parse CSV twice.
"""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

DONORS = ("P1", "P5", "P8", "P9", "P10", "P11", "P12", "P13")
COMMIT = "cbccfecca66b7896e0595905ad3c01f4e6484ecf"
MD5 = "a855987d31f5d0f1c082a55dc4763c61"
SOURCE_FILES = ("README.md", "import_and_process_visium_data.m",
                "fun/import_visium_data_funcion.m", "fun/compute_background_counts.m",
                "fun/read_spots2exclude.m")


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def barcode(value):
    """Only the reversible MATLAB hyphen/underscore conversion, case insensitive."""
    value = str(value).strip().strip('"').replace("_", "-").upper()
    if not re.fullmatch(r"[ACGTN]+-[0-9]+", value):
        raise ValueError("Unexpected barcode, refusing heuristic repair: " + repr(value))
    return value


def table_header(path):
    line = Path(path).open(encoding="utf-8-sig").readline()
    sep = ";" if line.count(";") > line.count(",") else ","
    values = next(csv.reader([line], delimiter=sep))
    names = [barcode(x) for x in values[1:]]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate matrix barcode: " + str(path))
    return sep, names


def read_positions(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    if rows and rows[0][1].strip() not in ("0", "1"):
        rows = rows[1:]
    columns = ["barcode", "in_tissue", "array_row", "array_col", "pixel_row", "pixel_col"]
    if any(len(row) != 6 for row in rows):
        raise ValueError("Positions must have six columns")
    df = pd.DataFrame(rows, columns=columns)
    df["barcode"] = df.barcode.map(barcode)
    if df.barcode.duplicated().any():
        raise ValueError("Duplicate positions barcode")
    for col in columns[1:]:
        df[col] = pd.to_numeric(df[col], errors="raise")
    if not np.isfinite(df[columns[1:]].to_numpy()).all():
        raise ValueError("Nonfinite positions")
    if not df.in_tissue.isin([0, 1]).all():
        raise ValueError("Invalid in_tissue annotation")
    if not np.equal(df[["array_row", "array_col"]], np.floor(df[["array_row", "array_col"]])).all().all():
        raise ValueError("Array coordinates are not integers")
    return df


def read_manual(path):
    # Author readtable takes every data row's first column, regardless of reason.
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    if df.shape[1] < 1:
        raise ValueError("Empty manual annotation schema")
    out = pd.DataFrame({"barcode": df.iloc[:, 0].map(barcode)})
    out["manual_reason"] = (df.iloc[:, 1:].astype(str).agg(" | ".join, axis=1)
                             if df.shape[1] > 1 else "author_manual_exclusion")
    return out.groupby("barcode", sort=False).manual_reason.agg(" | ".join).reset_index()


def read_count_cache(utt_path, all_path, out):
    """Validate ALL/UTT integer values and identical gene order, chunk by chunk."""
    us, ub = table_header(utt_path)
    als, ab = table_header(all_path)
    if not set(ub).issubset(ab):
        raise ValueError("UTT barcodes missing from ALL")
    outside = np.array([x not in set(ub) for x in ab])
    if not outside.any():
        raise ValueError("No ALL-minus-UTT spots: author mean background undefined")
    uit = pd.read_csv(utt_path, sep=us, chunksize=256, dtype={0: str})
    ait = pd.read_csv(all_path, sep=als, chunksize=256, dtype={0: str})
    blocks, genes, background = [], [], []
    sentinel = object()
    from itertools import zip_longest
    for u, a in zip_longest(uit, ait, fillvalue=sentinel):
        if u is sentinel or a is sentinel or len(u) != len(a):
            raise ValueError("ALL and UTT gene row count mismatch")
        ug, ag = u.iloc[:, 0].astype(str).tolist(), a.iloc[:, 0].astype(str).tolist()
        if ug != ag:
            raise ValueError("ALL and UTT gene order/names differ; refusing silent alignment")
        uv = u.iloc[:, 1:].to_numpy(dtype=np.float64)
        av = a.iloc[:, 1:].to_numpy(dtype=np.float64)
        for values in (uv, av):
            if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
                raise ValueError("Raw matrix is not finite nonnegative integer UMI")
            if values.max(initial=0) > np.iinfo(np.int32).max:
                raise ValueError("Raw value exceeds int32")
        # Author uses ALL \ UTT barcode set, not a reinterpretation of positions.
        background.append(av[:, outside].mean(axis=1))
        blocks.append(sparse.csr_matrix(uv.astype(np.int32)))
        genes.extend(ug)
    if not genes:
        raise ValueError("No genes in raw matrix")
    raw = sparse.vstack(blocks, format="csr")
    out.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(out / "raw_UTT_genes_by_spots.npz", raw)
    np.save(out / "background_mean_ALL_minus_UTT.npy", np.concatenate(background))
    pd.DataFrame({"original_gene_row": np.arange(len(genes)), "gene_symbol": genes,
                  "feature_id": [f"original_gene_row_{i:05d}" for i in range(len(genes))]}).to_csv(out / "genes.tsv", sep="\t", index=False)
    return ub, ab, genes, int(outside.sum()), int(raw.nnz)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--archive", type=Path)
    ap.add_argument("--author-source", type=Path)
    args = ap.parse_args()
    root = args.root.resolve()
    archive = args.archive or root / "data/raw/spatial/Data.zip"
    source = args.author_source or root / "data/reference/ST_Human_SI"
    meta = root / "data/metadata"
    meta.mkdir(parents=True, exist_ok=True)
    started = now()
    sha, md5 = digest(archive), digest(archive, "md5")
    if md5 != MD5:
        raise ValueError(f"Official ZIP MD5 mismatch: {md5}")
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if commit != COMMIT:
        raise ValueError("Author source commit differs from pre-frozen commit")
    source_rows = [{"file": x, "sha256": digest(source / x), "commit": commit,
                    "url": f"https://github.com/yotamharnik/ST_Human_SI/blob/{commit}/{x}"} for x in SOURCE_FILES]
    pd.DataFrame(source_rows).to_csv(meta / "SPATIAL_AUTHOR_SOURCE.tsv", sep="\t", index=False)
    extracted = root / "data/raw/spatial/extracted"
    with zipfile.ZipFile(archive) as z:
        bad = z.testzip()
        if bad:
            raise ValueError("ZIP CRC failure: " + bad)
        names = [n for n in z.namelist() if not n.endswith("/")]
        # No extraction of supplied processed models: only the explicitly named RAW inputs.
        required = ["metadata.csv"]
        for donor in DONORS:
            required += [f"{donor}_{suffix}" for suffix in ("counts_UTT.csv", "counts_ALL.csv", "tissue_positions_list.csv", "tissue_hires_image.png", "scalefactors_json.json")]
            required += [f"spots2exclude_{donor.lower()}.csv"]
        files = {}
        for filename in required:
            hits = [n for n in names if Path(n).name.lower() == filename.lower() and "/raw/" in ("/" + n.lower())]
            if len(hits) != 1:
                raise ValueError(f"Need exactly one official RAW {filename}; found {hits}")
            member = hits[0]
            path = extracted / filename
            path.parent.mkdir(parents=True, exist_ok=True)
            with z.open(member) as src, open(path, "wb") as dst:
                import shutil
                shutil.copyfileobj(src, dst)
            files[filename] = {"path": str(path), "zip_member": member, "sha256": digest(path)}
        inventory = pd.DataFrame({"zip_member": names})
        inventory["is_raw"] = inventory.zip_member.str.lower().str.contains("/raw/")
        inventory.to_csv(meta / "SPATIAL_ZIP_INVENTORY.tsv", sep="\t", index=False)
    clinical = pd.read_csv(extracted / "metadata.csv")
    clinical.to_csv(meta / "SPATIAL_CLINICAL_METADATA_ORIGINAL.tsv", sep="\t", index=False)
    provenance, flow, statuses = [], [], []
    for donor in DONORS:
        print(now(), "PHASE0_START", donor, flush=True)
        cache = root / "data/processed/spatial/raw_validated" / donor
        ub, ab, genes, nbg, nnz = read_count_cache(extracted / f"{donor}_counts_UTT.csv", extracted / f"{donor}_counts_ALL.csv", cache)
        positions = read_positions(extracted / f"{donor}_tissue_positions_list.csv")
        manual = read_manual(extracted / f"spots2exclude_{donor.lower()}.csv")
        if not set(ab).issubset(set(positions.barcode)):
            raise ValueError(donor + " ALL barcodes lack coordinates")
        if not set(manual.barcode).issubset(set(positions.barcode)):
            raise ValueError(donor + " manual barcodes cannot be matched to full array")
        p = positions.merge(manual, on="barcode", how="left", validate="one_to_one")
        p["donor"] = donor
        p["source_id"] = "zenodo:10715015:Data.zip:" + donor
        p["in_counts_ALL"] = p.barcode.isin(ab)
        p["in_counts_UTT"] = p.barcode.isin(ub)
        p["in_author_manual_exclusion"] = p.manual_reason.notna()
        p["manual_reason"] = p.manual_reason.fillna("")
        p["entered_primary_analysis"] = "PENDING_PHASE1"
        p.to_csv(cache / "all_array_barcodes.tsv", sep="\t", index=False)
        pd.DataFrame({"barcode": ub, "source_UTT_column": np.arange(len(ub))}).to_csv(cache / "UTT_order.tsv", sep="\t", index=False)
        scales = json.loads((extracted / f"{donor}_scalefactors_json.json").read_text())
        if not np.isfinite(float(scales["tissue_hires_scalef"])) or float(scales["tissue_hires_scalef"]) <= 0:
            raise ValueError("Invalid image scale")
        duplicate_genes = int(pd.Series(genes).duplicated(keep=False).sum())
        status = {"status": "RAW_SCHEMA_PASS", "donor": donor, "archive_sha256": sha,
                  "author_commit": commit, "n_raw_gene_rows": len(genes), "n_UTT_spots": len(ub),
                  "n_ALL_spots": len(ab), "n_background_spots": nbg, "raw_UTT_nnz": nnz,
                  "duplicate_symbol_rows": duplicate_genes, "duplicates_policy": "preserve_original_gene_rows; aggregate_after_background_clip_only_for_symbol_based_downstream",
                  "UTT_vs_position_in_tissue_disagreement": int(((p.in_counts_UTT) != p.in_tissue.eq(1)).sum()),
                  "files": {k: v for k, v in files.items() if k.startswith(donor + "_") or k == f"spots2exclude_{donor.lower()}.csv"},
                  "completed_at": now(), "is_new_raw_validation": True, "is_statistical_model": False}
        write_json(cache / "RAW_VALIDATION.json", status)
        statuses.append(status)
        provenance.append(p)
        for stage, n in (("positions_full_array", len(p)), ("counts_ALL", len(ab)), ("counts_UTT", len(ub)), ("background_ALL_minus_UTT", nbg), ("manual_full_array", len(manual)), ("manual_intersect_UTT", int((p.in_counts_UTT & p.in_author_manual_exclusion).sum()))):
            flow.append({"donor": donor, "stage": stage, "n_spots": n, "analysis_stage": "phase0_raw_schema"})
        print(now(), "PHASE0_RAW_SCHEMA_PASS", donor, len(ub), flush=True)
    pd.concat(provenance).to_csv(meta / "SPATIAL_BARCODE_PROVENANCE.tsv", sep="\t", index=False)
    pd.DataFrame(flow).to_csv(meta / "SPATIAL_SPOT_FLOW.tsv", sep="\t", index=False)
    write_json(meta / "SPATIAL_PHASE0_STATUS.json", {"status": "RAW_SCHEMA_PASS", "started_at": started,
               "finished_at": now(), "archive_sha256": sha, "archive_md5": md5, "zip_full_crc": "PASS",
               "author_commit": commit, "donors": statuses, "not_a_QC_or_model_completion_claim": True,
               "author_main_pipeline_scope": "Importer and background/QC will be reimplemented; pixel classification, zonation, apicome and cell signatures are separate author downstream steps, not claimed reproduced."})


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(now(), "PHASE0_FAILED", repr(exc), file=sys.stderr, flush=True)
        raise
