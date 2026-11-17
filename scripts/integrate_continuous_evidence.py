#!/usr/bin/env python3
"""Freeze continuous new-run evidence; deliberately perform no candidate ranking.

One current approved HGNC identity per master row. All original evidence rows,
including unmapped or ambiguous rows, remain in linked long tables. There is no
candidate name, candidate tier, composite score or gene-selection threshold.
Run this postprocessing module through the public run_logged.py, is-new-fit=0.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import time

PROJECT = "IBS_A_OFFICIAL_RERUN_20261105"
DONORS = ("P1", "P5", "P8", "P9", "P10", "P11", "P12", "P13")
CONTRASTS = ("Duo_IBSD_vs_HV", "Duo_IBSC_vs_HV", "Jej_IBSD_vs_HV", "Jej_IBSC_vs_HV")
MAP_FIELDS = ["identity_key", "identity_mapping_status", "identity_mapping_method", "identity_mapping_candidates", "identity_mapping_note"]


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--hgnc", default="data/reference/HGNC/hgnc_complete_set.txt")
    p.add_argument("--bulk-sensitivity-status", default="results/bulk/sensitivity/BULK_SENSITIVITY_STATUS.json")
    return p.parse_args()


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def json_write(path, value):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def clean(value):
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in ("", "nan", "None", "NA", "<NA>") else text


def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def is_true(value):
    return str(value).lower() in ("true", "1", "t")


class Identity:
    """Exact primary identifiers; explicit ambiguity, controlled symbol fallback."""
    def __init__(self, table):
        self.ensembl, self.entrez, self.symbol, self.legacy = (defaultdict(set) for _ in range(4))
        self.cache = {}
        self.master = {}
        for row in table.to_dict("records"):
            if row["status"] != "Approved":
                continue
            key = row["hgnc_id"]
            if key in self.master or not re.fullmatch(r"HGNC:\d+", key):
                raise ValueError("Duplicate/invalid approved HGNC identity")
            self.master[key] = {"identity_key": key, "gene_symbol": row["symbol"],
                                "hgnc_entrez_id": clean(row.get("entrez_id")),
                                "hgnc_ensembl_gene_id": clean(row.get("ensembl_gene_id")),
                                "hgnc_locus_group": clean(row.get("locus_group")),
                                "hgnc_locus_type": clean(row.get("locus_type")),
                                "observed_modules": set()}
            self.symbol[row["symbol"]].add(key)
            for field, lookup in (("entrez_id", self.entrez), ("ensembl_gene_id", self.ensembl)):
                for token in clean(row.get(field)).split("|"):
                    if token:
                        lookup[token].add(key)
            for field in ("prev_symbol", "alias_symbol"):
                for token in clean(row.get(field)).split("|"):
                    if token:
                        self.legacy[token].add(key)

    def map(self, ensembl="", entrez="", symbol=""):
        ensembl = re.sub(r"\.\d+$", "", clean(ensembl))
        entrez, symbol = clean(entrez), clean(symbol)
        cache_key = (ensembl, entrez, symbol)
        if cache_key in self.cache:
            return self.cache[cache_key]
        strong = [("Ensembl", self.ensembl.get(ensembl, set()))] if ensembl else []
        if entrez:
            strong.append(("Entrez", self.entrez.get(entrez, set())))
        union = set().union(*(v for _, v in strong)) if strong else set()
        exact = self.symbol.get(symbol, set()) if symbol else set()
        old = self.legacy.get(symbol, set()) if symbol and not exact else set()
        symbol_set = exact or old
        note = ""
        if any(len(v) > 1 for _, v in strong) or len(union) > 1:
            result = ("", "AMBIGUOUS_PRIMARY_ID", "NO_MAPPING", "|".join(sorted(union)), "No symbol-based rescue of ambiguous/conflicting primary IDs")
        elif len(union) == 1:
            chosen = next(iter(union))
            if symbol_set and chosen not in symbol_set:
                result = ("", "PRIMARY_ID_SYMBOL_CONFLICT", "NO_MAPPING", "|".join(sorted(union | symbol_set)), "Provided symbol disagrees with mapped stable ID; retained as unmapped evidence")
            else:
                if len(symbol_set) > 1:
                    note = "Ambiguous legacy symbol does not choose identity; unique primary ID is authoritative"
                method = "+".join(name for name, values in strong if values)
                result = (chosen, "MAPPED", "approved_" + method, chosen, note)
        elif len(symbol_set) == 1:
            chosen = next(iter(symbol_set))
            note = "Provided primary ID not found in approved HGNC; exact controlled symbol fallback" if strong else ""
            result = (chosen, "MAPPED", "approved_symbol" if exact else "unique_prev_or_alias_symbol", chosen, note)
        elif len(symbol_set) > 1:
            result = ("", "AMBIGUOUS_SYMBOL", "NO_MAPPING", "|".join(sorted(symbol_set)), "Do not pick an alias by expression or statistical significance")
        else:
            result = ("", "UNMAPPED", "NO_MAPPING", "", "No approved HGNC exact primary-ID/current-symbol/unique legacy-symbol match")
        self.cache[cache_key] = result
        return result


def main():
    args = arguments()
    root = args.root.resolve()
    if root.name != PROJECT:
        raise ValueError("Only the newly authorized clean-rerun project is allowed")
    import numpy as np
    import pandas as pd

    started, tick = now(), time.monotonic()
    out = root / "results/integrated/continuous"
    if out.exists() and any(out.iterdir()):
        raise FileExistsError("Continuous integration attempt exists: preserve/audit it explicitly, never silently overwrite frozen evidence")
    out.mkdir(parents=True, exist_ok=True)
    sources, schemas, declared_hashes = {}, {}, {}

    def path(rel):
        p = root / rel
        if not p.is_file() or p.is_symlink() or p.resolve() != p or not p.resolve().is_relative_to(root):
            raise ValueError(f"Source must be a real file inside this new root: {p}")
        return p

    def register(rel):
        p = path(rel)
        key = p.relative_to(root).as_posix()
        if key not in sources:
            sources[key] = {"sha256": digest(p), "bytes": p.stat().st_size}
        if key in declared_hashes and sources[key]["sha256"] != declared_hashes[key]:
            raise RuntimeError(f"Source differs from its completed model/gate hash: {key}")
        return p

    def load_json(rel):
        return json.loads(register(rel).read_text(encoding="utf-8-sig"))

    def load_table(rel, **kwargs):
        p = register(rel)
        frame = pd.read_csv(p, sep="\t", **kwargs)
        schemas[p.relative_to(root).as_posix()] = list(frame.columns)
        return frame

    for gate in ("PHASE1_GATE", "GWAS_MAGMA_GATE", "BULK_GATE", "BAYESSPACE_GATE", "SPATIALDM_GATE"):
        gate_record = load_json(f"manifest/{gate}.json")
        if gate_record.get("status") != "PASS":
            raise RuntimeError(f"Required new-run gate not PASS: {gate}")
        declared = dict(gate_record.get("output_sha256", {}))
        if isinstance(gate_record.get("evidence"), str) and isinstance(gate_record.get("sha256"), str):
            declared[gate_record["evidence"]] = gate_record["sha256"]
        for name, expected_hash in declared.items():
            if not isinstance(expected_hash, str) or not re.fullmatch("[0-9a-fA-F]{64}", expected_hash):
                continue
            p = path(name)
            declared_hashes[p.relative_to(root).as_posix()] = expected_hash.lower()
    if load_json("results/gwas/MAGMA_SUMMARY.json").get("is_new_fit") is not True:
        raise RuntimeError("MAGMA must certify a new fit")
    raw_effect_audit = load_json("results/gwas/qc/RAW_EFFECT_FIELD_AUDIT.json")
    raw_qc = load_json("results/gwas/qc/summary.json")
    if raw_effect_audit.get("status") != "DESCRIPTIVE_AUDIT_COMPLETE" or raw_effect_audit["data_rows"] != raw_qc["counts"]["raw_rows"]:
        raise RuntimeError("Full descriptive beta/SE/EAF/allele audit and actual raw QC row counts disagree")
    if raw_effect_audit["input_sha256"] != raw_qc["reference"]["inputs"]["raw_gwas"]["sha256"]:
        raise RuntimeError("Descriptive and core GWAS audits used different immutable raw files")
    for rel, expected_hash in raw_effect_audit["output_hashes_sha256"].items():
        if digest(register(rel)) != expected_hash:
            raise RuntimeError("Descriptive raw-field audit table hash changed")
    bulk_status = load_json("results/bulk/BULK_STATUS.json")
    if bulk_status.get("status") != "ANALYSIS_COMPLETE" or bulk_status.get("models_all_new_fit") is not True:
        raise RuntimeError("Bulk models are not all new/complete")
    bulk_hash_manifest = "results/bulk/main/OUTPUT_SHA256.tsv"
    declared_hashes[bulk_hash_manifest] = bulk_status["output_hashes_manifest_sha256"]
    for row in load_table(bulk_hash_manifest).to_dict("records"):
        p = path(row["path"])
        declared_hashes[p.relative_to(root).as_posix()] = row["sha256"]
    bsp_all = load_json("results/spatial/BayesSpace/BAYESSPACE_ALL_DONORS_STATUS.json")
    sdm_all = load_json("results/spatial/SpatialDM/SPATIALDM_ALL_DONORS_STATUS.json")
    if bsp_all.get("status") != "BAYESSPACE_COMPLETE" or bsp_all.get("n_donors") != 8:
        raise RuntimeError("BayesSpace all-donor status is incomplete")
    if sdm_all.get("status") != "SPATIALDM_COMPLETE" or sdm_all.get("n_donors") != 8:
        raise RuntimeError("SpatialDM all-donor status is incomplete")
    sensitivity_status = load_json(args.bulk_sensitivity_status)
    if sensitivity_status.get("status") != "ANALYSIS_COMPLETE":
        raise RuntimeError("Fixed bulk QC-sensitivity audit must finish before continuous integration")
    spatial_preprocess_protocol = load_json("data/metadata/SPATIAL_PREPROCESS_PROTOCOL.json")
    load_json("results/spatial/BayesSpace/FROZEN_MODEL_PROTOCOL.json")
    load_json("results/spatial/SpatialDM/FROZEN_MODEL_PROTOCOL.json")
    hgnc = load_table(args.hgnc, dtype=str, keep_default_na=False)
    hgnc_source_verification = load_json("data/reference/HGNC/SOURCE_VERIFICATION.json")
    if (hgnc_source_verification.get("status") != "OFFICIAL_REFERENCE_AUTHENTICATED_BY_REFETCH"
            or hgnc_source_verification.get("same_SHA256") is not True
            or hgnc_source_verification.get("original_reference_unchanged") is not True
            or hgnc_source_verification["original_reference"]["path"] != args.hgnc
            or digest(register(args.hgnc)) != hgnc_source_verification["original_reference"]["sha256"]):
        raise RuntimeError("Current HGNC identity reference is not the independently authenticated original snapshot")
    hgnc_headers = hgnc_source_verification["refetch"]["HTTP_headers_path"]
    if digest(register(hgnc_headers)) != hgnc_source_verification["refetch"]["HTTP_headers_sha256"]:
        raise RuntimeError("Official HGNC reference-authentication HTTP record changed")
    load_table("manifest/DATA_PROVENANCE.tsv", dtype=str, keep_default_na=False)
    # These are source-cohort limitations, not gene evidence or selection rules.
    # Bind them to the same provenance freeze without changing any spot mask.
    for rel in ("data/metadata/SPATIAL_AUTHOR_COHORT_SOURCE_AUDIT_20261105.md",
                "data/metadata/Harnik_published_cohort_audit/ARTICLE_COHORT_AUDIT.md",
                "data/metadata/Harnik_published_cohort_audit/ARTICLE_COHORT_AUDIT.json",
                "data/metadata/Harnik_published_cohort_audit/PUBLISHED_DONOR_COUNTS_AUDIT.tsv",
                "data/metadata/Harnik_published_cohort_audit/SOURCE_MANIFEST.tsv"):
        register(rel)
    for required in ("hgnc_id", "symbol", "status", "ensembl_gene_id", "entrez_id", "prev_symbol", "alias_symbol"):
        if required not in hgnc:
            raise ValueError(f"Official HGNC complete-set column missing: {required}")
    identity = Identity(hgnc)
    master = identity.master
    mapping_counts, evidence_counts = Counter(), Counter()
    representatives = {}
    written_headers = {}
    spatial_qc_thresholds, spatial_domain_ids = {}, {}

    def append_table(frame, filename, columns=None):
        if columns is not None:
            if set(frame.columns) - set(columns):
                raise ValueError(f"Unregistered fields in {filename}: {sorted(set(frame.columns) - set(columns))}")
            frame = frame.reindex(columns=columns)
        if filename in written_headers and list(frame.columns) != written_headers[filename]:
            raise ValueError(f"Source schema changed across rows of {filename}; no columns may be silently lost")
        first = filename not in written_headers
        frame.to_csv(out / filename, sep="\t", index=False, mode="w" if first else "a",
                     compression="gzip" if filename.endswith(".gz") else None,
                     header=first, na_rep="NA", float_format="%.17g")
        written_headers[filename] = list(frame.columns)

    def add_mapping(frame, module, source, feature_col, symbol_col=None, ensembl_col=None, entrez_col=None):
        mapped = [identity.map(r.get(ensembl_col, "") if ensembl_col else "",
                               r.get(entrez_col, "") if entrez_col else "",
                               r.get(symbol_col, "") if symbol_col else "") for r in frame.to_dict("records")]
        frame = frame.copy()
        frame[MAP_FIELDS] = pd.DataFrame(mapped, index=frame.index)
        unique = frame.drop_duplicates([feature_col] + ([symbol_col] if symbol_col else []))
        m = pd.DataFrame({"module": module, "source_path": source,
                          "original_feature_id": unique[feature_col].astype(str),
                          "provided_symbol": unique[symbol_col].astype(str) if symbol_col else "",
                          "provided_ensembl": unique[ensembl_col].astype(str) if ensembl_col else "",
                          "provided_entrez": unique[entrez_col].astype(str) if entrez_col else ""})
        for field in MAP_FIELDS:
            m[field] = unique[field]
        append_table(m.reset_index(drop=True), "FEATURE_IDENTITY_MAP.tsv.gz")
        mapping_counts.update(m.identity_mapping_status.value_counts().to_dict())
        unmapped = frame.loc[frame.identity_key.eq("")].copy()
        if len(unmapped):
            u = pd.DataFrame({"module": module, "source_path": source,
                              "original_feature_id": unmapped[feature_col].astype(str),
                              "provided_symbol": unmapped[symbol_col].astype(str) if symbol_col else ""})
            for field in MAP_FIELDS:
                u[field] = unmapped[field]
            append_table(u.reset_index(drop=True), "UNMAPPED_OR_AMBIGUOUS_SOURCE_RECORDS.tsv.gz")
        return frame

    fig_path = out / "figure_evidence_long.tsv.gz"
    fig_handle = gzip.open(fig_path, "wt", newline="", encoding="utf-8")
    figure = csv.writer(fig_handle, delimiter="\t")
    figure.writerow(["identity_key", "gene_symbol", "evidence_axis", "value", "display", "measurement_status", "source_path", "source_record_id", "evidence_table_path"])
    figure_rows = 0
    figure_axes = set()

    def fig(key, axis, value, status, source, record, display=None):
        nonlocal figure_rows
        if not key:
            return
        value_string = format(float(value), ".17g") if finite(value) else "NA"
        shown = display if display is not None else (format(float(value), ".4g") if finite(value) else "NA: " + status)
        long_name = next((name for prefix, name in (("MAGMA_", "MAGMA_GENE_EVIDENCE_LONG.tsv.gz"),
                         ("BULK_", "BULK_FEATURE_EVIDENCE_LONG.tsv.gz"), ("DETECT_", "SPATIAL_CORRECTED_DETECTION_LONG.tsv.gz"),
                         ("DOMAIN_", "SPATIAL_DOMAIN_AND_MARKER_EVIDENCE_LONG.tsv.gz"), ("LR_", "LR_COMPONENT_EVIDENCE_LONG.tsv.gz"))
                         if str(record).startswith(prefix)), "")
        evidence_table = "results/integrated/continuous/" + long_name if long_name else ""
        figure.writerow([key, master[key]["gene_symbol"], axis, value_string, shown, status, source, record, evidence_table])
        figure_rows += 1
        figure_axes.add(axis)

    # MAGMA remains all tested genes. Multiple Entrez rows mapped to one HGNC
    # never receive the minimum P as a fabricated single-gene statistic.
    gwas_source = "results/gwas/magma_all_genes.tsv"
    genes = load_table(gwas_source, dtype={"gene_id": str})
    genes["source_row_number"] = np.arange(2, len(genes) + 2)
    genes["source_path"] = gwas_source
    genes["source_record_id"] = [f"MAGMA_{i:06d}" for i in range(len(genes))]
    genes = add_mapping(genes, "MAGMA", gwas_source, "gene_id", "gene_symbol" if "gene_symbol" in genes else None, entrez_col="gene_id")
    append_table(genes, "MAGMA_GENE_EVIDENCE_LONG.tsv.gz")
    evidence_counts["MAGMA_source_rows"] = len(genes)
    for key, group in genes.loc[genes.identity_key.ne("")].groupby("identity_key", sort=False):
        row = master[key]
        row["observed_modules"].add("MAGMA")
        records = "|".join(group.source_record_id)
        row.update(MAGMA_source_record_ids=records, MAGMA_source_rows_n=len(group), MAGMA_source_path=gwas_source)
        if len(group) == 1:
            r = group.iloc[0]
            for field in ("P", "BH_FDR", "rank", "NSNPS", "NPARAM", "N", "ZSTAT"):
                if field in r:
                    row["MAGMA_" + field] = r[field]
            row["MAGMA_measurement_status"] = "SIGNIFICANT" if float(r.BH_FDR) < .05 else "NONSIGNIFICANT"
            fig(key, "MAGMA_P", r.P, row["MAGMA_measurement_status"], gwas_source, records)
            fig(key, "MAGMA_BH_FDR", r.BH_FDR, row["MAGMA_measurement_status"], gwas_source, records)
        else:
            row["MAGMA_measurement_status"] = "MULTIPLE_SOURCE_RECORDS"
            for field in ("P", "BH_FDR"):
                fig(key, "MAGMA_" + field, None, "MULTIPLE_SOURCE_RECORDS", gwas_source, records,
                    "all source values=" + "|".join(group[field].map(str)))

    # Four fixed main contrasts, representative chosen by expression rather
    # than P, and every original feature retained in the long table.
    bulk_long_columns = None
    for contrast in CONTRASTS:
        source = f"results/bulk/main/{contrast}.tsv"
        frame = load_table(source, dtype={"gene_id": str, "gene_id_stable": str})
        for field in ("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "P_NA_reason", "padj_NA_reason"):
            if field not in frame:
                raise ValueError(f"Bulk contrast lacks {field}: {source}")
        frame["source_path"] = source
        frame["source_row_number"] = np.arange(2, len(frame) + 2)
        frame["source_record_id"] = [f"BULK_{contrast}_{i:06d}" for i in range(len(frame))]
        frame = add_mapping(frame, "BULK_" + contrast, source, "gene_id", "gene_symbol", "gene_id")
        mapped = frame.loc[frame.identity_key.ne("")].copy()
        mapped["representative_baseMean_sort"] = pd.to_numeric(mapped.baseMean, errors="coerce")
        selected = mapped.sort_values(["identity_key", "representative_baseMean_sort", "gene_id"],
                                      ascending=[True, False, True], na_position="last", kind="mergesort").drop_duplicates("identity_key")
        selected_ids = set(selected.source_record_id)
        frame["canonical_representative"] = frame.source_record_id.isin(selected_ids)
        frame["canonical_representative_rule"] = "maximum finite baseMean; then lexical original gene_id; if all baseMean NA use lexical ID with explicit NA"
        if bulk_long_columns is None:
            bulk_long_columns = list(frame.columns)
        append_table(frame, "BULK_FEATURE_EVIDENCE_LONG.tsv.gz", bulk_long_columns)
        evidence_counts["BULK_" + contrast + "_source_rows"] = len(frame)
        all_records = mapped.groupby("identity_key", sort=False).source_record_id.agg("|".join).to_dict()
        for r in selected.to_dict("records"):
            key = r["identity_key"]
            row = master[key]
            row["observed_modules"].add("BULK")
            representatives[(contrast, key)] = r["gene_id"]
            for field in ("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "P_NA_reason", "padj_NA_reason", "direction"):
                if field in r:
                    row[contrast + "_" + field] = r[field]
            row[contrast + "_all_feature_record_ids"] = all_records[key]
            row[contrast + "_representative_source_record_id"] = r["source_record_id"]
            row[contrast + "_source_path"] = source
            measurement = (("SIGNIFICANT" if float(r["padj"]) < .05 else "NONSIGNIFICANT") if finite(r["padj"])
                           else "FDR_UNAVAILABLE" if finite(r["pvalue"]) else "UNTESTED")
            row[contrast + "_measurement_status"] = measurement
            row[contrast + "_native_pvalue_tested"] = finite(r["pvalue"])
            row[contrast + "_native_FDR_available"] = finite(r["padj"])
            for field in ("log2FoldChange", "padj"):
                status = (measurement if measurement in ("SIGNIFICANT", "NONSIGNIFICANT") else "OBSERVED_DESCRIPTIVE") if finite(r[field]) else "NA_NATIVE_STATISTIC"
                reason = r.get("padj_NA_reason", "NATIVE_NA") if field == "padj" else r.get("P_NA_reason", "NATIVE_NA")
                fig(key, "BULK_" + contrast + "_" + field, r[field], status, source, r["source_record_id"],
                    display=None if finite(r[field]) else "NA: " + clean(reason))

    # Fixed bulk sensitivity uses the SAME main representative feature when
    # examining a changed cohort. Identical current-run cohorts are linked,
    # rather than performing or claiming a redundant statistical refit.
    sens_rows = []
    for label, record in sensitivity_status.get("analyses", {}).items():
        native = record.get("status")
        identical = native == "IDENTICAL_TO_NEW_MAIN_INPUT_NO_REDUNDANT_REFIT"
        if not identical and native != "ANALYSIS_COMPLETE":
            raise RuntimeError(f"Bulk sensitivity incomplete: {label}")
        for contrast in CONTRASTS:
            source = f"results/bulk/main/{contrast}.tsv" if identical else f"results/bulk/sensitivity/{label}/{contrast}.tsv"
            sensitivity = load_table(source, dtype={"gene_id": str}).set_index("gene_id", drop=False)
            if sensitivity.index.duplicated().any():
                raise ValueError("Duplicate feature ID in bulk sensitivity")
            for (main_contrast, key), feature in representatives.items():
                if main_contrast != contrast:
                    continue
                r = sensitivity.loc[feature] if feature in sensitivity.index else {}
                z = {"identity_key": key, "gene_symbol": master[key]["gene_symbol"], "contrast": contrast,
                     "sensitivity": label, "original_feature_id": feature, "source_path": source,
                     "sensitivity_native_status": native, "representative_rule": "retain primary max-baseMean representative; no sensitivity P-based replacement"}
                for field in ("baseMean", "log2FoldChange", "pvalue", "padj", "P_NA_reason", "padj_NA_reason"):
                    z[field] = r.get(field, "FEATURE_NOT_PRESENT")
                    master[key][contrast + "_" + label + "_" + field] = z[field]
                sens_rows.append(z)
    if not sens_rows:
        raise RuntimeError("Bulk sensitivity status has no declared analyses")
    append_table(pd.DataFrame(sens_rows), "BULK_QC_SENSITIVITY_EVIDENCE_LONG.tsv.gz")

    # Spatial domains/markers retain all native feature-by-domain statistics.
    # Native scran comparison columns vary with donor q, so create a union
    # schema before appending any donor, never discard unequal comparisons.
    domain_columns = set()
    detection_columns = None
    for donor in DONORS:
        expression_path = f"results/spatial/BayesSpace/{donor}/domain_gene_expression.tsv"
        marker_path = f"results/spatial/BayesSpace/{donor}/markers_all.tsv"
        expression_header = load_table(expression_path, nrows=0).columns
        marker_header = load_table(marker_path, nrows=0).columns
        domain_columns.update(expression_header)
        domain_columns.update("marker_" + x for x in marker_header if x not in ("feature_id", "domain", "gene_symbol"))
    domain_columns.update(MAP_FIELDS + ["source_path", "source_row_number", "marker_source_path", "marker_source_row_number", "source_record_id"])
    domain_columns = ["donor", "domain", "feature_id", "gene_symbol"] + sorted(domain_columns - {"donor", "domain", "feature_id", "gene_symbol"})
    for donor in DONORS:
        bsp_source = f"results/spatial/BayesSpace/{donor}/BAYESSPACE_STATUS.json"
        bsp = load_json(bsp_source)
        if bsp.get("status") != "BAYESSPACE_COMPLETE" or bsp.get("is_new_fit") is not True:
            raise RuntimeError(f"BayesSpace native donor gate failed: {donor}")
        if bsp.get("nrep_formal") != 50000 or bsp.get("nrep_stability") != 50000:
            raise RuntimeError("Actual BayesSpace chains differ from fixed protocol")
        pre = f"data/processed/spatial/{donor}"
        qc_scope = load_json(f"{pre}/QC_SENSITIVITY_SCOPE.json")
        if qc_scope.get("status") != "DESCRIPTIVE_QC_SENSITIVITY_COMPLETE":
            raise RuntimeError(f"Descriptive QC sensitivity audit incomplete: {donor}")
        qc = load_table(f"{pre}/QC_SENSITIVITY.tsv")
        thresholds = sorted(set(float(v) for v in qc.absolute_z_threshold))
        spatial_qc_thresholds[donor] = thresholds
        for label, rel in [("primary", f"{pre}/gene_detection.tsv")] + [(f"QC_z{v:.1f}", f"{pre}/gene_detection_sensitivity_z{v:.1f}.tsv") for v in thresholds]:
            frame = load_table(rel, dtype={"gene_symbol": str})
            if frame.gene_symbol.duplicated().any():
                raise ValueError("Symbol-aggregated spatial detection has duplicate source symbols")
            frame["source_path"] = rel
            frame["source_row_number"] = np.arange(2, len(frame) + 2)
            frame["source_record_id"] = [f"DETECT_{donor}_{label}_{i:06d}" for i in range(len(frame))]
            frame["QC_label"] = label
            frame["absolute_z_threshold"] = np.nan if label == "primary" else float(label.split("z")[-1])
            frame = add_mapping(frame, "SPATIAL_DETECTION_" + donor + "_" + label, rel, "gene_symbol", "gene_symbol")
            if detection_columns is None:
                detection_columns = list(frame.columns)
            append_table(frame, "SPATIAL_CORRECTED_DETECTION_LONG.tsv.gz", detection_columns)
            for key, group in frame.loc[frame.identity_key.ne("")].groupby("identity_key", sort=False):
                row = master[key]
                row["observed_modules"].add("SPATIAL")
                prefix = donor + "_" + label
                refs = "|".join(group.source_record_id)
                row[prefix + "_source_record_ids"] = refs
                row[prefix + "_source_path"] = rel
                row[prefix + "_source_symbols_n"] = len(group)
                row[prefix + "_measurement_status"] = ("OBSERVED_DESCRIPTIVE" if int(group.iloc[0].n_spots) > 0 else "UNTESTED_empty_QC_cohort") if len(group) == 1 else "MULTIPLE_SOURCE_RECORDS"
                if label == "primary":
                    row[donor + "_BayesSpace_independent_seed_ARI"] = bsp.get("independent_seed_ARI")
                    row[donor + "_BayesSpace_stability_warning"] = bsp.get("stability_warning")
                    row[donor + "_BayesSpace_q"] = bsp.get("q")
                    row[donor + "_BayesSpace_status_source"] = bsp_source
                    row[donor + "_QC_sensitivity_scope_source"] = f"{pre}/QC_SENSITIVITY_SCOPE.json"
                    row[donor + "_QC_sensitivity_primary_rule"] = qc_scope.get("primary_rule")
                    row[donor + "_QC_sensitivity_strict_rule"] = qc_scope.get("strict_rule")
                    row[donor + "_BayesSpace_strict_QC_mask_sensitivity_status"] = qc_scope.get("BayesSpace_on_strict_mask", "NOT_REPORTED")
                for source_field, suffix in (("detection_fraction", "corrected_detection_fraction"), ("mean_log_expression", "mean_log1p_CP10000")):
                    axis = "SPATIAL_" + donor + "_" + ("" if label == "primary" else label + "_") + suffix
                    if len(group) == 1:
                        value = group.iloc[0][source_field]
                        row[prefix + "_" + suffix] = value
                        if label == "primary":
                            row[donor + "_primary_n_positive_spots"] = group.iloc[0].n_positive_spots
                            row[donor + "_primary_n_spots"] = group.iloc[0].n_spots
                        fig(key, axis, value, "OBSERVED_DESCRIPTIVE" if finite(value) else "UNTESTED", rel, refs,
                            display=None if finite(value) else "NA: EMPTY_QC_COHORT_OR_NATIVE_NA")
                    else:
                        # Summing detection counts or log means across aliases
                        # would be incorrect without returning to corrected counts.
                        fig(key, axis, None, "MULTIPLE_SOURCE_RECORDS", rel, refs,
                            "not collapsed; source values=" + "|".join(group[source_field].map(str)))
        qc_effects_source = f"{pre}/QC_SENSITIVITY_GENE_EFFECTS.tsv"
        qc_effects = load_table(qc_effects_source, dtype={"gene_symbol": str})
        qc_effects["source_path"] = qc_effects_source
        qc_effects["source_row_number"] = np.arange(2, len(qc_effects) + 2)
        qc_effects = add_mapping(qc_effects, "SPATIAL_QC_EFFECT_" + donor, qc_effects_source, "gene_symbol", "gene_symbol")
        append_table(qc_effects, "SPATIAL_QC_SENSITIVITY_GENE_EFFECTS_LONG.tsv.gz")
        expr_source = f"results/spatial/BayesSpace/{donor}/domain_gene_expression.tsv"
        marker_source = f"results/spatial/BayesSpace/{donor}/markers_all.tsv"
        expression = load_table(expr_source, dtype={"feature_id": str, "domain": str})
        spatial_domain_ids[donor] = sorted(set(expression.domain), key=lambda x: int(x))
        markers = load_table(marker_source, dtype={"feature_id": str, "domain": str})
        expression["source_row_number"] = np.arange(2, len(expression) + 2)
        markers["source_row_number"] = np.arange(2, len(markers) + 2)
        markers = markers.rename(columns={c: "marker_" + c for c in markers if c not in ("feature_id", "domain")}).drop(columns="marker_gene_symbol", errors="ignore")
        domain = expression.merge(markers, on=["feature_id", "domain"], how="outer", validate="one_to_one", indicator=True)
        if not domain._merge.eq("both").all():
            raise ValueError("Domain expression / native marker feature-domain universes disagree")
        domain = domain.drop(columns="_merge")
        domain["source_path"], domain["marker_source_path"] = expr_source, marker_source
        domain["source_record_id"] = [f"DOMAIN_{donor}_{i:07d}" for i in range(len(domain))]
        domain = add_mapping(domain, "SPATIAL_DOMAIN_" + donor, expr_source, "feature_id", "gene_symbol")
        append_table(domain, "SPATIAL_DOMAIN_AND_MARKER_EVIDENCE_LONG.tsv.gz", domain_columns)
        for key, group in domain.loc[domain.identity_key.ne("")].groupby("identity_key", sort=False):
            master[key][donor + "_domain_marker_source_record_ids"] = "|".join(group.source_record_id)
            master[key][donor + "_domain_marker_source_path"] = "results/integrated/continuous/SPATIAL_DOMAIN_AND_MARKER_EVIDENCE_LONG.tsv.gz"
            master[key][donor + "_domain_marker_records_n"] = len(group)
            for domain_id, domain_group in group.groupby("domain", sort=False):
                if "marker_FDR" not in domain_group:
                    raise ValueError("scran native marker table lacks FDR")
                refs = "|".join(domain_group.source_record_id)
                if len(domain_group) == 1:
                    q = domain_group.iloc[0].marker_FDR
                    status = ("SIGNIFICANT" if float(q) < .05 else "NONSIGNIFICANT") if finite(q) else "UNTESTED"
                    fig(key, f"SPATIAL_{donor}_domain_{domain_id}_marker_FDR", q, status, marker_source, refs)
                else:
                    fig(key, f"SPATIAL_{donor}_domain_{domain_id}_marker_FDR", None, "MULTIPLE_SOURCE_RECORDS", marker_source, refs,
                        "not collapsed; source FDR=" + "|".join(domain_group.marker_FDR.map(str)))
        evidence_counts["SPATIAL_" + donor + "_feature_domain_rows"] = len(domain)

    # All catalogue pair eligibility rows, not only significant/selected pairs.
    # Components of complete complexes are attributed separately, while each
    # gene's displayed pair counts are deduplicated by donor+pair_id.
    lr_columns = None
    for donor in DONORS:
        prefix = f"results/spatial/SpatialDM/{donor}"
        native = load_json(f"{prefix}/SPATIALDM_STATUS.json")
        if native.get("status") != "SPATIALDM_COMPLETE" or native.get("is_new_fit") is not True or native.get("n_perm_global") != 1000:
            raise RuntimeError(f"Native SpatialDM donor status incomplete: {donor}")
        elig_source, global_source = f"{prefix}/ALL_LR_ELIGIBILITY.tsv", f"{prefix}/global_LR.tsv"
        geometry_source = f"{prefix}/global_LR_geometry_sensitivity_radius1.5.tsv"
        local_source, plan_source = f"{prefix}/local_spot_statistics.tsv.gz", f"{prefix}/LOCAL_SELECTION_FROZEN.tsv"
        elig = load_table(elig_source, dtype={"pair_id": str}, keep_default_na=False)
        glob = load_table(global_source, dtype={"pair_id": str})
        if sources[global_source]["sha256"] != native["global_table_sha256"]:
            raise RuntimeError("Global LR table differs from completed native model status")
        geometry = load_table(geometry_source, dtype={"pair_id": str})
        plan = load_table(plan_source, dtype={"pair_id": str})
        for frame, name in ((elig, elig_source), (glob, global_source), (geometry, geometry_source), (plan, plan_source)):
            if frame.pair_id.duplicated().any():
                raise ValueError(f"Duplicate pair in native donor table: {name}")
        if set(glob.pair_id) != set(geometry.pair_id) or set(glob.pair_id) != set(plan.pair_id):
            raise ValueError("Global, geometry-sensitivity and local-plan pair universes disagree")
        glob["source_row_number"] = np.arange(2, len(glob) + 2)
        geometry["source_row_number"] = np.arange(2, len(geometry) + 2)
        elig["eligibility_source_row_number"] = np.arange(2, len(elig) + 2)
        glob = glob.rename(columns={c: "global_" + c for c in glob if c != "pair_id"})
        geometry = geometry.rename(columns={c: "geometry_" + c for c in geometry if c != "pair_id"})
        pairs = elig.merge(glob, on="pair_id", how="left", validate="one_to_one").merge(geometry, on="pair_id", how="left", validate="one_to_one")
        pairs = pairs.merge(plan[["pair_id", "selected_for_local", "local_status"]], on="pair_id", how="left", validate="one_to_one")
        eligible = pairs.eligible.map(is_true)
        if not pd.to_numeric(pairs.loc[eligible, "global_padj"], errors="coerce").notna().all():
            raise ValueError("Eligible LR pair lacks an actual global test")
        if pairs.loc[~eligible, "global_padj"].notna().any():
            raise ValueError("Ineligible LR pair unexpectedly has a global test")
        local_file = register(local_source)
        local_agg = {}
        for chunk in pd.read_csv(local_file, sep="\t", dtype={"pair_id": str}, chunksize=200000):
            schemas[local_source] = list(chunk.columns)
            for pair_id, group in chunk.groupby("pair_id", sort=False):
                v = local_agg.setdefault(pair_id, [0, 0])
                v[0] += len(group)
                v[1] += int(group.significant.map(is_true).sum())
        if set(local_agg) != set(plan.loc[plan.selected_for_local.map(is_true), "pair_id"]):
            raise ValueError("Local actual result pair set differs from frozen local selection")
        records = []
        for pair in pairs.to_dict("records"):
            tested = is_true(pair["eligible"])
            supported = tested and float(pair["global_padj"]) < .05
            lp = local_agg.get(pair["pair_id"], [0, 0])
            for role in ("ligand", "receptor"):
                for symbol in dict.fromkeys(clean(pair[role + "_components"]).split(";")):
                    if not symbol:
                        raise ValueError("Empty required LR component")
                    r = dict(pair)
                    r.update(donor=donor, component_role=role, component_symbol=symbol,
                             global_tested=tested, global_supported=supported,
                             global_evidence_status="SUPPORT" if supported else "NO_SUPPORT" if tested else "UNTESTED_" + clean(pair["test_status"]),
                             local_tested=bool(lp[0]), local_pair_spot_tests_n=lp[0], local_pair_spot_BH_supported_n=lp[1],
                             eligibility_source_path=elig_source, global_source_path=global_source,
                             geometry_source_path=geometry_source, local_source_path=local_source, local_selection_source_path=plan_source)
                    records.append(r)
        frame = pd.DataFrame(records)
        frame["source_record_id"] = [f"LR_{donor}_{i:07d}" for i in range(len(frame))]
        frame = add_mapping(frame, "LR_COMPONENT_" + donor, elig_source, "component_symbol", "component_symbol")
        if lr_columns is None:
            lr_columns = list(frame.columns)
        append_table(frame, "LR_COMPONENT_EVIDENCE_LONG.tsv.gz", lr_columns)
        for key, group in frame.loc[frame.identity_key.ne("")].groupby("identity_key", sort=False):
            distinct = group.drop_duplicates("pair_id")
            vals = {"global_pairs_tested_n": int(distinct.global_tested.sum()),
                    "global_pairs_BH_supported_n": int(distinct.global_supported.sum()),
                    "global_pairs_tested_no_support_n": int((distinct.global_tested & ~distinct.global_supported).sum()),
                    "global_pairs_untested_n": int((~distinct.global_tested).sum()),
                    "local_pairs_tested_n": int(distinct.local_tested.sum()),
                    "local_pairs_with_supported_spots_n": int((distinct.local_pair_spot_BH_supported_n > 0).sum()),
                    "geometry_sensitivity_pairs_BH_supported_n": int((pd.to_numeric(distinct.geometry_padj, errors="coerce") < .05).sum())}
            row = master[key]
            row["observed_modules"].add("LR_DATABASE_COMPONENT")
            row["LR_" + donor + "_strict_QC_mask_global_sensitivity_status"] = native.get("strict_z3_5_QC_mask_global_LR_sensitivity", "NOT_REPORTED")
            row["LR_" + donor + "_strict_QC_mask_local_sensitivity_status"] = native.get("strict_z3_5_QC_mask_local_LR_sensitivity", "NOT_REPORTED")
            row["LR_" + donor + "_geometry_sensitivity_completed"] = native.get("geometry_sensitivity_completed")
            row["LR_" + donor + "_source_record_ids"] = "|".join(group.source_record_id)
            row["LR_" + donor + "_eligibility_reasons_json"] = json.dumps(distinct.test_status.value_counts().to_dict(), sort_keys=True)
            for field, value in vals.items():
                row["LR_" + donor + "_" + field] = value
                fig(key, "LR_" + donor + "_" + field, value, "OBSERVED_DESCRIPTIVE", global_source if "global" in field else geometry_source if "geometry" in field else local_source,
                    row["LR_" + donor + "_source_record_ids"])
            row["LR_" + donor + "_global_measurement_status"] = "TESTED_WITH_SUPPORT" if vals["global_pairs_BH_supported_n"] else "TESTED_NO_SUPPORT" if vals["global_pairs_tested_n"] else "UNTESTED"
        evidence_counts["LR_" + donor + "_component_evidence_rows"] = len(frame)

    # Explicit absence semantics for every registry gene. Registry-only rows
    # remain visible but never inflate an asserted candidate/gene discovery count.
    for key, row in master.items():
        modules = row["observed_modules"]
        row["new_source_observed"] = bool(modules)
        row["observed_modules"] = "|".join(sorted(modules)) if modules else "NO_NEW_EVIDENCE"
        row.setdefault("MAGMA_measurement_status", "NOT_OBSERVED_IN_TESTED_MAGMA_UNIVERSE")
        if modules and row["MAGMA_measurement_status"].startswith("NOT_OBSERVED"):
            for field in ("P", "BH_FDR"):
                fig(key, "MAGMA_" + field, None, "NOT_OBSERVED", gwas_source, "", "NA: NOT_OBSERVED_IN_TESTED_MAGMA_UNIVERSE")
        for contrast in CONTRASTS:
            row.setdefault(contrast + "_measurement_status", "NOT_OBSERVED")
            for field in ("P_NA_reason", "padj_NA_reason"):
                row.setdefault(contrast + "_" + field, "NO_UNIQUE_HGNC_MAPPED_FEATURE_IN_BULK_SOURCE")
            if modules and row[contrast + "_measurement_status"] == "NOT_OBSERVED":
                for field in ("log2FoldChange", "padj"):
                    fig(key, "BULK_" + contrast + "_" + field, None, "NOT_OBSERVED", f"results/bulk/main/{contrast}.tsv", "",
                        "NA: NO_UNIQUE_HGNC_MAPPED_FEATURE_IN_BULK_SOURCE")
        for donor in DONORS:
            row.setdefault(donor + "_primary_measurement_status", "NOT_OBSERVED")
            row.setdefault("LR_" + donor + "_global_measurement_status", "NOT_APPLICABLE_not_a_uniquely_mapped_database_component")
            if modules and row[donor + "_primary_measurement_status"] == "NOT_OBSERVED":
                for field in ("corrected_detection_fraction", "mean_log1p_CP10000"):
                    fig(key, "SPATIAL_" + donor + "_" + field, None, "NOT_OBSERVED", f"data/processed/spatial/{donor}/gene_detection.tsv", "",
                        "NA: NO_UNIQUE_HGNC_MAPPED_SYMBOL_IN_SPATIAL_SOURCE")
                for threshold in spatial_qc_thresholds[donor]:
                    for field in ("corrected_detection_fraction", "mean_log1p_CP10000"):
                        fig(key, f"SPATIAL_{donor}_QC_z{threshold:.1f}_" + field, None, "NOT_OBSERVED",
                            f"data/processed/spatial/{donor}/gene_detection_sensitivity_z{threshold:.1f}.tsv", "",
                            "NA: NO_UNIQUE_HGNC_MAPPED_SYMBOL_IN_SPATIAL_QC_SOURCE")
                for domain_id in spatial_domain_ids[donor]:
                    fig(key, f"SPATIAL_{donor}_domain_{domain_id}_marker_FDR", None, "NOT_OBSERVED",
                        f"results/spatial/BayesSpace/{donor}/markers_all.tsv", "",
                        "NA: NO_UNIQUE_HGNC_MAPPED_FEATURE_IN_NATIVE_DOMAIN_MARKER_SOURCE")
            if modules and row["LR_" + donor + "_global_measurement_status"].startswith("NOT_APPLICABLE"):
                for field in ("global_pairs_tested_n", "global_pairs_BH_supported_n", "global_pairs_tested_no_support_n", "global_pairs_untested_n", "local_pairs_tested_n", "local_pairs_with_supported_spots_n", "geometry_sensitivity_pairs_BH_supported_n"):
                    fig(key, "LR_" + donor + "_" + field, None, "NOT_APPLICABLE", f"results/spatial/SpatialDM/{donor}/ALL_LR_ELIGIBILITY.tsv", "",
                        "NA: NOT_A_UNIQUELY_MAPPED_DATABASE_COMPONENT; NOT_A_NEGATIVE_TEST")
        row["LR_donors_tested_n"] = sum(row.get("LR_" + donor + "_global_pairs_tested_n", 0) > 0 for donor in DONORS)
        row["LR_donors_with_supported_pairs_n"] = sum(row.get("LR_" + donor + "_global_pairs_BH_supported_n", 0) > 0 for donor in DONORS)
        row["LR_total_tested_donor_pair_records_n"] = sum(row.get("LR_" + donor + "_global_pairs_tested_n", 0) for donor in DONORS)
        row["LR_total_supported_donor_pair_records_n"] = sum(row.get("LR_" + donor + "_global_pairs_BH_supported_n", 0) for donor in DONORS)
    fig_handle.close()
    wide = pd.DataFrame(master.values()).sort_values("identity_key", key=lambda s: s.str.replace("HGNC:", "", regex=False).astype(int))
    if wide.identity_key.duplicated().any() or wide.gene_symbol.duplicated().any():
        raise ValueError("Master is not one approved canonical gene per row")
    append_table(wide, "CONTINUOUS_EVIDENCE_WIDE.tsv.gz")
    # Re-read hashes of every registered statistical source to exclude an
    # integration racing a later writer, even when a gate was already present.
    for rel, record in sources.items():
        if digest(path(rel)) != record["sha256"]:
            raise RuntimeError(f"A source changed during integration: {rel}")
    freeze = {"status": "FROZEN", "phase6_continuous_complete": True, "frozen_at": now(),
              "new_root": str(root), "source_files": sources, "source_schemas": schemas,
              "identity_reference": args.hgnc, "identity_key": "approved HGNC ID",
              "identity_reference_source_verification": "data/reference/HGNC/SOURCE_VERIFICATION.json",
              "identity_reference_source_url": hgnc_source_verification["source_url"],
              "identity_reference_original_download_time": hgnc_source_verification["original_reference"]["initial_retrieved_at"],
              "spatial_preprocessing_protocol": spatial_preprocess_protocol,
              "raw_effect_field_audit_scope": raw_effect_audit["scope"],
              "figure_evidence_axes": sorted(figure_axes),
              "ranking_performed": False, "candidate_tiers_defined": False, "target_selected": False,
              "no_old_results_or_rankings_read": True,
              "bulk_representative_rule": "max baseMean per canonical HGNC and contrast; lexical original gene ID breaks ties; full features linked",
              "ambiguous_identity_policy": "unmapped evidence retained outside canonical master; no candidate-count inflation",
              "uncollapsible_statistics_policy": "do not combine multiple gene P, log means or detection fractions; preserve every source row and explicit MULTIPLE_SOURCE_RECORDS status",
              "LR_count_unit": "within donor, unique pair per canonical gene regardless of duplicate roles/components; across donors, donor-pair records, not independent loci/genes",
              "interpretation": "continuous evidence catalogue only; primary biological replicate is donor; marker/local spot evidence is not donor-level disease association; spatial coexpression is not direct binding or pain causality"}
    freeze_file = out / "INTEGRATION_SOURCE_FREEZE.json"
    json_write(freeze_file, freeze)
    summary = {"status": "PASS", "start_iso8601": started, "finish_iso8601": now(), "duration_seconds": time.monotonic() - tick,
               "approved_HGNC_registry_rows": len(wide), "registry_genes_with_any_new_source_record": int(wide.new_source_observed.sum()),
               "registry_only_rows_not_candidate_discoveries": int((~wide.new_source_observed).sum()),
               "feature_mapping_record_counts_not_unique_gene_counts": dict(mapping_counts),
               "source_evidence_record_counts": dict(evidence_counts), "figure_evidence_records": figure_rows,
               "figure_evidence_axes": sorted(figure_axes),
               "ranking_performed": False, "candidate_freeze_performed": False,
               "source_freeze_sha256": digest(freeze_file), "figure_evidence_file": fig_path.relative_to(root).as_posix(),
               "figure_evidence_sha256": digest(fig_path)}
    json_write(out / "CONTINUOUS_INTEGRATION_SUMMARY.json", summary)
    with (out / "INTEGRATION_OUTPUT_MANIFEST.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["relative_path", "bytes", "sha256"])
        for p in sorted(out.iterdir()):
            if p.is_file() and p.name != "INTEGRATION_OUTPUT_MANIFEST.tsv":
                writer.writerow([p.relative_to(root).as_posix(), p.stat().st_size, digest(p)])
    json_write(root / "manifest/INTEGRATION_CONTINUOUS_GATE.json",
               {"status": "PASS", "phase6_continuous_complete": True, "ranking_performed": False,
                "finish_iso8601": now(), "source_freeze": freeze_file.relative_to(root).as_posix(),
                "source_freeze_sha256": digest(freeze_file), "evidence_file": fig_path.relative_to(root).as_posix(),
                "evidence_sha256": digest(fig_path), "wide_table": "results/integrated/continuous/CONTINUOUS_EVIDENCE_WIDE.tsv.gz",
                "figure_evidence_axes": sorted(figure_axes),
                "wide_table_sha256": digest(out / "CONTINUOUS_EVIDENCE_WIDE.tsv.gz")})
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"CONTINUOUS_INTEGRATION_FAIL: {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        raise
