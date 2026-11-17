#!/usr/bin/env python3
"""Independent, explicit Python reimplementation of the frozen author importer.

Not an execution of MATLAB.  Background and QC operate on original gene rows.
Raw integer UMI and noninteger background-corrected expression remain distinct.
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from spatial_phase0 import COMMIT, DONORS, digest, now, write_json

PROTOCOL = {"implementation": "independent_Python_reimplementation_of_author_importer_not_MATLAB",
            "author_commit": COMMIT, "Z_THRESH_SUM_UMI": -4.0, "Z_THRESH_MITO": 4.0,
            "HEMO_FRAC": "Inf", "BOUNDARY_PRCTILE": 100, "DIST_THRESH": "Inf_passed_but_unused_in_author_function",
            "standard_deviation": "sample_ddof_1; propagate_nonfinite_as_MATLAB_default",
            "mitochondrial_pattern": "case_insensitive_contains_mt-", "hemoglobin_patterns": "case_sensitive_Hba-_or_Hbb-",
            "manual_exclusions": "all_rows_all_reasons", "sensitivity_absolute_z": [3.5],
            "boundary": "author_axiswise_distance_to_midpoint; percentile_100_removes_no_finite_spot",
            "candidate_blind": True, "publication_spot_target": "none_not_targeted"}


def zscore_matlab(values):
    with np.errstate(invalid="ignore", divide="ignore"):
        return (values - np.mean(values)) / np.std(values, ddof=1)


def normalized_log(matrix):
    sums = np.asarray(matrix.sum(axis=1)).ravel()
    if (sums <= 0).any() or not np.isfinite(sums).all():
        raise ValueError("Retained spot has undefined corrected normalization")
    norm = sparse.diags(10000.0 / sums) @ matrix
    norm = norm.tocsr()
    norm.data = np.log1p(norm.data)
    return norm


def aggregate_symbol(matrix, symbols):
    labels = pd.Index(pd.unique(symbols))
    codes = labels.get_indexer(symbols)
    mapping = sparse.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(labels)))
    return (matrix @ mapping).tocsr(), labels


def detection_table(donor, raw_spots_genes, corrected_spots_genes, symbols, retained):
    rr, labels = aggregate_symbol(raw_spots_genes[retained], symbols)
    cc, _ = aggregate_symbol(corrected_spots_genes[retained], symbols)
    nr = np.asarray((rr > 0).sum(axis=0)).ravel()
    nc = np.asarray((cc > 0).sum(axis=0)).ravel()
    n = len(retained)
    if n:
        ll = normalized_log(cc)
        mean_log = np.asarray(ll.mean(axis=0)).ravel()
        mean_raw = np.asarray(rr.mean(axis=0)).ravel()
        mean_corrected = np.asarray(cc.mean(axis=0)).ravel()
    else:
        # Empty cohorts have no defined mean/detection fraction; never call this zero expression.
        mean_log = mean_raw = mean_corrected = np.full(len(labels), np.nan)
    return pd.DataFrame({"donor": donor, "gene_symbol": labels, "n_spots": len(retained),
                         "n_positive_spots": nc, "detection_fraction": nc / n if n else np.nan,
                         "n_positive_spots_raw_UMI": nr,
                         "detection_fraction_raw_UMI": nr / n if n else np.nan,
                         "mean_log_expression": mean_log,
                         "mean_raw_UMI": mean_raw, "mean_background_corrected_counts": mean_corrected,
                         "expression_layer": "background_corrected_CP10000_log1p; detection_corrected_gt0"})


def sensitivity_gene_comparison(primary, strict):
    """Unranked all-symbol paired descriptive changes, not independent donor tests."""
    columns = ["n_spots", "n_positive_spots", "detection_fraction", "n_positive_spots_raw_UMI",
               "detection_fraction_raw_UMI", "mean_log_expression", "mean_raw_UMI", "mean_background_corrected_counts"]
    keys = ["donor", "gene_symbol"]
    a = primary[keys + columns].rename(columns={c: c + "_primary_z4" for c in columns})
    b = strict[keys + columns].rename(columns={c: c + "_strict_z3_5" for c in columns})
    result = a.merge(b, on=keys, how="outer", validate="one_to_one", indicator=True)
    if not result._merge.eq("both").all():
        raise ValueError("Primary/strict sensitivity tables do not contain the same full gene universe")
    result = result.drop(columns="_merge")
    for c in columns:
        result[c + "_delta_strict_minus_primary"] = result[c + "_strict_z3_5"] - result[c + "_primary_z4"]
    result["comparison_scope"] = "descriptive_expression_and_detection_only; no_BayesSpace_or_LR_model_sensitivity_claim"
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--donors", nargs="+", choices=DONORS, default=list(DONORS))
    args = ap.parse_args()
    root = args.root.resolve()
    gate = json.loads((root / "manifest/PHASE0_GATE.json").read_text())
    if gate.get("status") != "PASS":
        raise ValueError("Global Phase 0 gate is not PASS")
    phase0 = json.loads((root / "data/metadata/SPATIAL_PHASE0_STATUS.json").read_text())
    if phase0.get("status") != "RAW_SCHEMA_PASS" or phase0["author_commit"] != COMMIT:
        raise ValueError("New spatial raw audit/source gate failed")
    meta, results = root / "data/metadata", root / "results/spatial"
    results.mkdir(parents=True, exist_ok=True)
    # Freeze before loading gene expression or applying QC. Refuse protocol changes.
    protocol_path = meta / "SPATIAL_PREPROCESS_PROTOCOL.json"
    if protocol_path.exists() and json.loads(protocol_path.read_text()) != PROTOCOL:
        raise ValueError("Preprocessing protocol already frozen differently")
    write_json(protocol_path, PROTOCOL)
    all_flow, all_provenance, all_detection, all_sensitivity, all_sensitivity_masks, summaries = [], [], [], [], [], []
    for donor in args.donors:
        start = now()
        cache = root / "data/processed/spatial/raw_validated" / donor
        out = root / "data/processed/spatial" / donor
        out.mkdir(parents=True, exist_ok=True)
        validation = json.loads((cache / "RAW_VALIDATION.json").read_text())
        if (validation.get("status") != "RAW_SCHEMA_PASS" or validation.get("author_commit") != COMMIT
                or validation["archive_sha256"] != phase0["archive_sha256"]):
            raise ValueError("Raw cache belongs to different input archive")
        raw = sparse.load_npz(cache / "raw_UTT_genes_by_spots.npz").tocsr()
        genes = pd.read_csv(cache / "genes.tsv", sep="\t", keep_default_na=False)
        bg = np.load(cache / "background_mean_ALL_minus_UTT.npy")
        utt = pd.read_csv(cache / "UTT_order.tsv", sep="\t")
        full = pd.read_csv(cache / "all_array_barcodes.tsv", sep="\t", keep_default_na=False)
        if (raw.shape != (validation["n_raw_gene_rows"], validation["n_UTT_spots"])
                or raw.nnz != validation["raw_UTT_nnz"] or not np.issubdtype(raw.dtype, np.integer)
                or (raw.data < 0).any() or not np.isfinite(raw.data).all()):
            raise ValueError("Raw cache dimensions/nnz/integer nonnegative UMI contract failed")
        if bg.shape != (raw.shape[0],) or not np.isfinite(bg).all() or (bg < 0).any():
            raise ValueError("Author mean-background vector shape/nonnegative finite contract failed")
        if (not np.array_equal(genes.original_gene_row.to_numpy(), np.arange(raw.shape[0]))
                or not genes.feature_id.is_unique
                or genes.feature_id.tolist() != [f"original_gene_row_{i:05d}" for i in range(raw.shape[0])]):
            raise ValueError("Original gene-row bookkeeping was reordered or is not unique")
        if (not utt.barcode.is_unique or not full.barcode.is_unique
                or not np.array_equal(utt.source_UTT_column.to_numpy(), np.arange(raw.shape[1]))):
            raise ValueError("Raw UTT column-order/barcode provenance contract failed")
        obs = utt.merge(full, on="barcode", how="left", validate="one_to_one", sort=False)
        if obs.pixel_row.isna().any() or len(obs) != raw.shape[1] or raw.shape[0] != len(genes):
            raise ValueError("Raw/cache barcode or feature mismatch")
        if (obs.barcode.tolist() != utt.barcode.tolist()
                or not obs.in_counts_UTT.astype(str).str.lower().eq("true").all()
                or not np.isfinite(obs[["pixel_col", "pixel_row", "array_row", "array_col"]].to_numpy(float)).all()):
            raise ValueError("Raw UTT barcode/coordinate mapping is not exact and finite")
        # Subtract on each ORIGINAL row, clip then preserve sparsity; zero raw entries stay zero.
        corrected = raw.astype(np.float64)
        corrected.data = np.maximum(corrected.data - np.repeat(bg, np.diff(corrected.indptr)), 0.0)
        corrected.eliminate_zeros()
        corrected_sums = np.asarray(corrected.sum(axis=0)).ravel()
        raw_sums = np.asarray(raw.sum(axis=0)).ravel()
        symbols = genes.gene_symbol.astype(str)
        mito = symbols.str.lower().str.contains("mt-", regex=False).to_numpy()
        hemo = (symbols.str.contains("Hba-", regex=False) | symbols.str.contains("Hbb-", regex=False)).to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            log_sum = np.log10(corrected_sums)
            fraction_mito = np.asarray(corrected[mito].sum(axis=0)).ravel() / corrected_sums
            fraction_hemo = np.asarray(corrected[hemo].sum(axis=0)).ravel() / corrected_sums
        zu, zm = zscore_matlab(log_sum), zscore_matlab(fraction_mito)
        if not np.isfinite(zu).all() or not np.isfinite(zm).all():
            write_json(out / "PREPROCESS_FAILED.json", {"status": "STOP_UNDEFINED_AUTHOR_ZSCORES", "donor": donor,
                        "zero_corrected_total": int((corrected_sums == 0).sum()), "nonfinite_UMI_z": int((~np.isfinite(zu)).sum()),
                        "nonfinite_mito_z": int((~np.isfinite(zm)).sum()), "note": "No nanmean/nanstd fallback: original author MATLAB defaults would propagate nonfinite statistics."})
            raise ValueError(donor + " author default zscores undefined; audit required")
        # All finite coordinates pass the author 100th-percentile boundary rule.
        xy = obs[["pixel_col", "pixel_row"]].to_numpy(float)
        midpoint = xy.min(axis=0) + (xy.max(axis=0) - xy.min(axis=0)) / 2
        distances = np.abs(xy - midpoint)
        boundary = (distances > np.percentile(distances, 100, axis=0)).any(axis=1)
        manual = obs.in_author_manual_exclusion.astype(str).str.lower().eq("true").to_numpy()
        step_umi, step_mito = zu > -4, zm < 4
        step_hemo = fraction_hemo < np.inf
        retained_mask = step_umi & step_mito & step_hemo & ~boundary & ~manual
        strict_umi, strict_mito = zu > -3.5, zm < 3.5
        strict_mask = strict_umi & strict_mito & step_hemo & ~boundary & ~manual
        if (strict_mask & ~retained_mask).any():
            raise ValueError("Strict predeclared QC mask is not a subset of primary")
        keep = np.flatnonzero(retained_mask)
        if len(keep) < 50:
            raise ValueError(donor + " fewer than 50 primary spots after exact author logic")
        obs["raw_UMI"] = raw_sums
        obs["raw_nGenes"] = np.asarray((raw > 0).sum(axis=0)).ravel()
        obs["background_corrected_total"] = corrected_sums
        obs["corrected_nGenes"] = np.asarray((corrected > 0).sum(axis=0)).ravel()
        obs["corrected_mito_fraction"] = fraction_mito
        obs["raw_mito_fraction"] = np.asarray(raw[mito].sum(axis=0)).ravel() / raw_sums
        obs["author_Z_log10_corrected_total"] = zu
        obs["author_Z_corrected_mito_fraction"] = zm
        obs["author_corrected_hemoglobin_fraction_exact_pattern"] = fraction_hemo
        obs["passed_UMI_Z"] = step_umi
        obs["passed_mito_Z"] = step_mito
        obs["passed_Hb_Inf"] = step_hemo
        obs["boundary_excluded"] = boundary
        obs["entered_primary_analysis"] = retained_mask
        obs["passed_strict_z3_5_UMI"] = strict_umi
        obs["passed_strict_z3_5_mito"] = strict_mito
        obs["entered_strict_z3_5_sensitivity"] = strict_mask
        reasons = []
        strict_reasons = []
        for i in range(len(obs)):
            r = [name for flag, name in ((not step_umi[i], "author_low_corrected_UMI_Z"), (not step_mito[i], "author_high_mito_Z"), (not step_hemo[i], "author_Hb"), (boundary[i], "author_boundary"), (manual[i], "author_manual:" + str(obs.manual_reason.iloc[i]))) if flag]
            reasons.append(";".join(r))
            sr = [name for flag, name in ((not strict_umi[i], "strict_low_corrected_UMI_Z_le_minus3.5"),
                  (not strict_mito[i], "strict_high_mito_Z_ge3.5"), (not step_hemo[i], "author_Hb"),
                  (boundary[i], "author_boundary"), (manual[i], "author_manual:" + str(obs.manual_reason.iloc[i]))) if flag]
            strict_reasons.append(";".join(sr))
        obs["exclusion_reason"] = reasons
        obs["strict_z3_5_exclusion_reason"] = strict_reasons
        obs.to_csv(out / "UTT_QC_all_spots.tsv", sep="\t", index=False)
        raw_sg, corrected_sg = raw.T.tocsr(), corrected.T.tocsr()
        kept_obs = obs.iloc[keep].copy().set_index("barcode", drop=False)
        kept_obs.index.name = "spot_id"
        var = genes.set_index("feature_id", drop=False)
        obj = ad.AnnData(normalized_log(corrected_sg[keep]), obs=kept_obs, var=var)
        obj.layers["counts_raw_integer_UMI"] = raw_sg[keep]
        obj.layers["counts_author_background_corrected"] = corrected_sg[keep]
        fraction = sparse.diags(1.0 / corrected_sums[keep]) @ corrected_sg[keep]
        obj.layers["author_count_fraction"] = fraction.tocsr()
        obj.obsm["spatial"] = xy[keep]
        scales = json.loads((root / "data/raw/spatial/extracted" / f"{donor}_scalefactors_json.json").read_text())
        obj.obsm["spatial_hires"] = xy[keep] * float(scales["tissue_hires_scalef"])
        obj.uns["preprocess"] = PROTOCOL
        obj.uns["donor"] = donor
        obj.uns["archive_sha256"] = phase0["archive_sha256"]
        obj.uns["layer_contract"] = {"X": "log1p(CP10000(author_background_corrected_counts))",
                    "counts_raw_integer_UMI": "unaltered integer UTT counts at retained spots",
                    "counts_author_background_corrected": "max(original_gene_row_raw - ALL_minus_UTT_mean, 0); noninteger",
                    "gene_duplicates": "original gene rows retained; aggregate only after background clip for symbol based tests"}
        obj.write_h5ad(out / "author_reimplementation_primary.h5ad", compression="gzip")
        from scipy.io import mmwrite
        mmwrite(out / "raw_integer_UMI_genes_by_spots.mtx", raw[:, keep])
        mmwrite(out / "author_background_corrected_genes_by_spots.mtx", corrected[:, keep])
        genes.to_csv(out / "genes.tsv", sep="\t", index=False)
        kept_obs.to_csv(out / "spots.tsv", sep="\t", index=False)
        pd.DataFrame({"original_gene_row": genes.original_gene_row, "gene_symbol": symbols, "background_mean": bg}).to_csv(out / "background_by_original_gene_row.tsv", sep="\t", index=False)
        primary_detection = detection_table(donor, raw_sg, corrected_sg, symbols, keep)
        primary_detection.to_csv(out / "gene_detection.tsv", sep="\t", index=False)
        sensitivity = []
        strict_detection = None
        for threshold in (3.5, 4.0):
            ids = np.flatnonzero((zu > -threshold) & (zm < threshold) & step_hemo & ~boundary & ~manual)
            detected = detection_table(donor, raw_sg, corrected_sg, symbols, ids)
            detected.assign(absolute_z_threshold=threshold).to_csv(out / f"gene_detection_sensitivity_z{threshold:.1f}.tsv", sep="\t", index=False)
            if threshold == 3.5:
                strict_detection = detected
            sensitivity.append({"donor": donor, "absolute_z_threshold": threshold, "n_spots": len(ids), "same_all_manual_exclusions": True})
        pd.DataFrame(sensitivity).to_csv(out / "QC_SENSITIVITY.tsv", sep="\t", index=False)
        gene_sensitivity = sensitivity_gene_comparison(primary_detection, strict_detection)
        gene_sensitivity.to_csv(out / "QC_SENSITIVITY_GENE_EFFECTS.tsv", sep="\t", index=False)
        write_json(out / "QC_SENSITIVITY_SCOPE.json", {"status": "DESCRIPTIVE_QC_SENSITIVITY_COMPLETE", "donor": donor,
                   "primary_rule": "log10_corrected_UMI_Z_gt_minus4_and_mito_Z_lt4", "strict_rule": "log10_corrected_UMI_Z_gt_minus3.5_and_mito_Z_lt3.5",
                   "same_author_background_and_all_manual_exclusions": True, "n_primary_spots": int(retained_mask.sum()),
                   "n_strict_spots": int(strict_mask.sum()), "n_symbols_all_reported": len(gene_sensitivity),
                   "gene_subset_or_candidate_selection": "none", "BayesSpace_on_strict_mask": "NOT_PERFORMED",
                   "global_LR_on_strict_mask": "NOT_PERFORMED", "local_LR_on_strict_mask": "NOT_PERFORMED",
                   "do_not_infer": "Expression/detection and mask sensitivity does not establish domain, LR-P, or candidate-ranking model robustness.",
                   "SpatialDM_geometry_sensitivity": "separate_planned_primary_cohort_secreted_radius1.5; completion_requires_its_own_model_status"})
        counts = [("counts_UTT", len(obs)), ("after_corrected_UMI_Z_gt_minus4", int(step_umi.sum())),
                  ("after_mito_Z_lt4", int((step_umi & step_mito).sum())),
                  ("after_Hb_Inf_and_boundary100", int((step_umi & step_mito & step_hemo & ~boundary).sum())),
                  ("after_all_author_manual_exclusions_PRIMARY", len(keep))]
        flow = pd.DataFrame([{"donor": donor, "stage": name, "n_spots": n, "analysis_stage": "phase1_author_independent_reimplementation"} for name, n in counts])
        flow.to_csv(out / "SPOT_FLOW.tsv", sep="\t", index=False)
        primary = obs[["barcode", "entered_primary_analysis", "exclusion_reason", "entered_strict_z3_5_sensitivity", "strict_z3_5_exclusion_reason"]]
        full = full.drop(columns=["entered_primary_analysis"]).merge(primary, on="barcode", how="left", validate="one_to_one")
        full["entered_primary_analysis"] = full.entered_primary_analysis.fillna(False)
        full["exclusion_reason"] = full.exclusion_reason.fillna("not_in_author_UTT_input")
        full["entered_strict_z3_5_sensitivity"] = full.entered_strict_z3_5_sensitivity.fillna(False)
        full["strict_z3_5_exclusion_reason"] = full.strict_z3_5_exclusion_reason.fillna("not_in_author_UTT_input")
        full.to_csv(out / "BARCODE_PROVENANCE.tsv", sep="\t", index=False)
        sensitivity_mask = full[["donor", "barcode", "in_counts_ALL", "in_counts_UTT", "in_author_manual_exclusion", "manual_reason",
                                 "entered_primary_analysis", "exclusion_reason", "entered_strict_z3_5_sensitivity", "strict_z3_5_exclusion_reason"]]
        sensitivity_mask.to_csv(out / "QC_SENSITIVITY_BARCODE_MASK.tsv", sep="\t", index=False)
        status = {"status": "PREPROCESS_COMPLETE", "donor": donor, "n_input_UTT": len(obs), "n_retained": len(keep),
                  "manual_intersect_UTT": int(manual.sum()), "manual_retained": int((manual & retained_mask).sum()),
                  "n_original_gene_rows": raw.shape[0], "background_nonzero_genes": int((bg > 0).sum()),
                  "started_at": start, "finished_at": now(), "input_archive_sha256": phase0["archive_sha256"],
                  "output_h5ad_sha256": digest(out / "author_reimplementation_primary.h5ad"),
                  "mode": PROTOCOL["implementation"], "is_new_fit": False, "is_new_preprocessing": True,
                  "publication_15145_reproduced": None, "no_forced_spot_count": True,
                  "n_strict_z3_5_sensitivity_spots": int(strict_mask.sum()),
                  "QC_sensitivity_scope": "all_gene_expression_detection_and_full_array_barcode_masks_only; strict_mask_BayesSpace_and_LR_NOT_PERFORMED",
                  "clinical_interpretation": "surgical proximal small intestine reference; not IBS patients or unscreened healthy volunteers"}
        write_json(out / "PREPROCESS_STATUS.json", status)
        summaries.append(status)
        all_flow.append(flow)
        all_provenance.append(full)
        all_detection.append(pd.read_csv(out / "gene_detection.tsv", sep="\t"))
        all_sensitivity.append(gene_sensitivity)
        all_sensitivity_masks.append(sensitivity_mask)
        print(now(), "PREPROCESS_COMPLETE", donor, len(keep), flush=True)
    # Only canonical eight-donor aggregates get a completed all-donor status.
    if tuple(args.donors) == DONORS:
        raw_flow = pd.read_csv(meta / "SPATIAL_SPOT_FLOW.tsv", sep="\t")
        raw_flow = raw_flow[raw_flow.analysis_stage.eq("phase0_raw_schema")]
        pd.concat([raw_flow] + all_flow).to_csv(meta / "SPATIAL_SPOT_FLOW.tsv", sep="\t", index=False)
        pd.concat(all_provenance).to_csv(meta / "SPATIAL_BARCODE_PROVENANCE.tsv", sep="\t", index=False)
        pd.concat(all_detection).to_csv(results / "gene_detection_by_donor.tsv", sep="\t", index=False)
        pd.concat(all_sensitivity).to_csv(results / "QC_SENSITIVITY_GENE_EFFECTS_all_donors.tsv", sep="\t", index=False)
        pd.concat(all_sensitivity_masks).to_csv(meta / "SPATIAL_QC_SENSITIVITY_BARCODE_MASK_all_donors.tsv", sep="\t", index=False)
        write_json(meta / "SPATIAL_PREPROCESS_ALL_DONORS_STATUS.json", {"status": "PREPROCESS_COMPLETE", "n_donors": 8, "donors": summaries,
                    "total_retained": sum(x["n_retained"] for x in summaries), "author_pipeline_scope": "independent_reimplementation_of_import_background_QC_not_full_zonation_or_MATLAB_pipeline"})


if __name__ == "__main__":
    main()
