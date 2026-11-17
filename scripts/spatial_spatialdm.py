#!/usr/bin/env python3
"""Fresh official SpatialDM 0.3.1, all complete/detectable human LR pairs.

Primary inference is 1000-permutation donor-wise testing, analytic results are
separate.  Local selection is the predeclared global-BH/P/pair-ID top-100 cap.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import argparse
import csv
import gzip
import importlib.metadata
import json
import pickle
import random
import socket
import time
import traceback
from contextlib import contextmanager
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
import spatialdm as sdm

from spatial_phase0 import DONORS, digest, now, write_json
from spatial_preprocess import aggregate_symbol, normalized_log

VERSION = importlib.metadata.version("spatialdm")
PROTOCOL = {"version": "0.3.1", "global_method": "both_primary_permutation_secondary_zscore",
            "global_n_perm": 1000, "global_BH_threshold": .05, "local_n_perm": 1000,
            "local_selection": "all_global_permutation_BH_lt0.05_capped100_sorted_BH_P_pairID",
            "local_cap": 100, "local_FDR_family": "all_tested_local_pair_by_spot_cells_per_donor",
            "complete_complexes": True, "minimum_positive_spots_per_required_subunit": 3,
            "coordinates": "full_resolution_xy_divided_by_median_donor_nearest_spacing",
            "secreted_radius_spacing": 2.0, "secreted_cutoff": .2,
            "contact_radius_spacing": 1.2, "single_cell": False,
            "geometry_sensitivity_secreted_radius_spacing": 1.5,
            "X": "log1p(CP10000(background_corrected_counts)); sum_duplicate_symbols_after_original_row_clip",
            "local_raw_interface": "nonlog_CP10000_of_author_background_corrected_expression_using_ALL_gene_total_then_subset_required_genes; NOT_raw_integer_UMI",
            "local_library_normalization": "project_specific_predeclared_author_mat_norm_alignment; permitted_by_API_but_NOT_claimed_as_official_rawcount_tutorial_default",
            "seed": "20261105+100*fixed_zero_based_donor_index+50; local+1; geometry_sensitivity+2",
            "nproc": 1, "candidate_or_pathway_shortlist": "none"}


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if not np.isfinite(p).all() or (p < 0).any() or (p > 1).any():
        raise ValueError("Nonfinite or invalid p-values; not silently replacing")
    return multipletests(p, method="fdr_bh")[1]


def local_cp10000(corrected, gene_indices):
    """Normalize by ALL corrected genes per spot BEFORE selecting LR subunits.

    The official local API applies per-gene max/SD scaling, not library-size
    normalization.  Its .raw slot here is deliberately normalized nonlog
    expression, not integer UMI and not the raw-count tutorial default.
    """
    totals = np.asarray(corrected.sum(axis=1)).ravel()
    if not np.isfinite(totals).all() or (totals <= 0).any():
        raise ValueError("Undefined whole-transcriptome corrected library normalization")
    selected = (sparse.diags(10000.0 / totals) @ corrected[:, gene_indices]).tocsr()
    if not np.isfinite(selected.data).all() or (selected.data < 0).any():
        raise ValueError("Invalid local nonlog CP10000 expression")
    return selected, totals


def dump_model(path, obj):
    temp = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(temp, "wb", compresslevel=1) as f:
        pickle.dump(obj, f, protocol=4)
    temp.replace(path)


def read_model(path):
    # The only pickles accepted are this script's independently fitted new-project checkpoints.
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def audit_global_model(obj, out, label, seed):
    """Verify actual official null-array dimensions and the empirical-P identity."""
    null = np.asarray(obj.uns["global_stat"]["perm"]["global_perm"])
    observed = np.asarray(obj.uns["global_I"])
    pairs = obj.uns["global_res"].index.astype(str).tolist()
    if null.shape != (len(pairs), 1000) or observed.shape != (len(pairs),):
        raise ValueError("Global actual permutation-array shape differs from 1000 per tested pair")
    if not np.isfinite(null).all() or not np.isfinite(observed).all():
        raise ValueError("Nonfinite official global permutation statistic")
    # Official global source computes 1 - count(I > null)/1000, including ties in the upper tail.
    reconstructed = np.mean(null >= observed[:, None], axis=1)
    reported = obj.uns["global_res"]["perm_pval"].to_numpy(float)
    if not np.allclose(reconstructed, reported, rtol=0, atol=1e-12):
        raise ValueError("Official global empirical P cannot be reconstructed from saved null statistics")
    record = {"stage": "global", "version": VERSION, "seed": int(seed), "method_requested": "both",
              "primary_P_method": "permutation", "nproc": 1,
              "nproc_source_semantics": "BLAS_thread_limit_only; no_fork_or_multiprocessing_in_0.3.1_test_path",
              "n_perm_requested": 1000, "actual_null_statistics_per_pair": int(null.shape[1]),
              "shape_global_perm_pairs_by_permutations": list(null.shape), "global_perm_dtype": str(null.dtype),
              "shape_global_I": list(observed.shape), "n_spots": obj.n_obs,
              "saved_empirical_P_reconstructed": True,
              "rng_source": "NumPy_global_RandomState_seeded_immediately_before_official_call; advancing_permutation_loop",
              "global_null_semantics": "joint_same_index_permutation_of_L_and_R; one_1000_iteration_loop_shared_across_LR_family",
              "unique_draw_guarantee": "not_claimed; random_sampling_may_naturally_repeat_a_permutation; no_fork_cloning_path",
              "package_main_sha256": digest(Path(sdm.__file__).parent / "main.py"),
              "package_utils_sha256": digest(Path(sdm.__file__).parent / "utils.py")}
    obj.uns.setdefault("CLEAN_RERUN_EXECUTION", {})[label] = record
    write_json(out / f"{label}_MODEL_ARRAY_AUDIT.json", record)
    return record


def audit_local_model(obj, out, seed):
    """Check 1000 null statistics per pair/spot without allocating a full 3D sum."""
    left = np.asarray(obj.uns["local_stat"]["local_permI"])
    right = np.asarray(obj.uns["local_stat"]["local_permI_R"])
    p = obj.uns["unadjusted_local_perm_p_saved_before_sig_spots"]
    n_pairs, n_spots = p.shape
    if left.shape != (n_pairs, 1000, n_spots) or right.shape != left.shape:
        raise ValueError("Local actual permutation arrays differ from 1000 per tested pair/spot")
    direct = np.asarray(obj.uns["local_stat"]["local_I"])
    reverse = np.asarray(obj.uns["local_stat"]["local_I_R"])
    if direct.shape != (n_spots, n_pairs) or reverse.shape != direct.shape:
        raise ValueError("Local observed statistics have incorrect orientation")
    # The official source sets P=1 where neither standardized component is positive.
    # These masked cells cannot be reconstructed from null arrays alone; verify all
    # other cells exactly and retain the full masked P table and component model.
    verified = 0
    for j in range(n_pairs):
        if not np.isfinite(left[j]).all() or not np.isfinite(right[j]).all():
            raise ValueError("Nonfinite local null statistics")
        reconstructed = np.mean((left[j] + right[j]) >= (direct[:, j] + reverse[:, j])[None, :], axis=0)
        reported = p.iloc[j].to_numpy(float)
        # Official positivity masking only changes cells to one; no other discrepancy is permitted.
        matched = np.isclose(reconstructed, reported, rtol=0, atol=1e-12)
        if not np.all(matched | (reported == 1)):
            raise ValueError("Official local empirical P differs from saved null statistics beyond positivity masking")
        verified += int(matched.sum())
    annotations = obj.uns["geneInter"].loc[p.index, "annotation"]
    groups = int((annotations != "Secreted Signaling").any()) + int((annotations == "Secreted Signaling").any())
    record = {"stage": "local", "version": VERSION, "seed": int(seed), "method_requested": "both",
              "primary_P_method": "permutation", "nproc": 1, "n_perm_requested": 1000,
              "actual_null_statistics_per_tested_pair_spot": int(left.shape[1]),
              "shape_local_permI_pairs_by_permutations_by_spots": list(left.shape),
              "shape_local_permI_R": list(right.shape), "shape_local_I_spots_by_pairs": list(direct.shape),
              "shape_saved_unadjusted_P_pairs_by_spots": list(p.shape),
              "shape_local_FDR_pairs_by_spots": list(obj.uns["local_stat"]["local_fdr"].shape),
              "n_empirical_P_cells_reconstructed_exactly": verified,
              "unit_P_cells_differing_from_unmasked_null_tail": int(n_pairs * n_spots - verified),
              "unit_P_difference_interpretation": "compatible_with_official_positivity_mask; positivity_indicator_not_independently_reconstructed",
              "expected_NumPy_shuffle_calls_from_official_source": 1000 * groups,
              "shuffle_call_count_is_source_inference_not_instrumented_timer": True,
              "null_semantics": "each_nonempty_contact_or_secreted_group_has_1000_advancing_NumPy_shuffle_iterations; both_directional_statistics_share_each_group_iteration_index",
              "local_FDR_family": "all_selected_pairs_times_spots_within_donor",
              "raw_interface_layer": "noninteger_nonlog_CP10000_of_author_background_corrected_expression; ALL_gene_library_denominator_before_gene_subset; not_raw_UMI",
              "nproc_source_semantics": "BLAS_thread_limit_only; no_fork_or_multiprocessing_in_0.3.1_test_path",
              "package_main_sha256": digest(Path(sdm.__file__).parent / "main.py"),
              "package_utils_sha256": digest(Path(sdm.__file__).parent / "utils.py")}
    obj.uns.setdefault("CLEAN_RERUN_EXECUTION", {})["local"] = record
    write_json(out / "local_MODEL_ARRAY_AUDIT.json", record)
    return record


@contextmanager
def timing(out, donor, stage, input_sha, parameters):
    start_iso, start = now(), time.monotonic()
    code, note = 0, "actual_official_API_call"
    try:
        yield
    except BaseException as e:
        code, note = 1, repr(e)
        raise
    finally:
        p = out / "MODEL_EXECUTION_TRACE.tsv"
        api = "spatialdm_global" if stage == "global_geometry_sensitivity" else stage
        record = {"module": f"SpatialDM:{donor}:{stage}", "source_command": "spatialdm." + api + " " + json.dumps(parameters, sort_keys=True),
                  "input_hashes": input_sha, "start_iso8601": start_iso, "finish_iso8601": now(),
                  "duration_seconds": time.monotonic() - start, "exit_code": code,
                  "hostname": socket.gethostname(), "software_version": VERSION,
                  "run_id": f"{donor}_{stage}_{start_iso}_{os.getpid()}", "is_new_fit": int(stage != "extract_lr"),
                  "reused_only_immutable_raw_inputs": True, "notes": note}
        exists = p.exists()
        with p.open("a", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(record), delimiter="\t")
            if not exists:
                w.writeheader()
            w.writerow(record)


def database(out):
    import importlib.resources as resources
    base = resources.files("spatialdm").joinpath("datasets", "LR_data")
    ip = base.joinpath("human-interaction_input_CellChatDB.csv.gz")
    cp = base.joinpath("human-complex_input_CellChatDB.csv")
    with ip.open("rb") as f:
        inter = pd.read_csv(f, index_col=0, compression="gzip")
    with cp.open("rb") as f:
        complexes = pd.read_csv(f, index_col=0)
    inter.index = inter.index.astype(str)
    if not inter.index.is_unique:
        raise ValueError("CellChatDB pair indices are not unique")
    inter.to_csv(out / "CellChatDB_human_interactions_full.tsv", sep="\t")
    complexes.to_csv(out / "CellChatDB_human_complexes_full.tsv", sep="\t")
    write_json(out / "DATABASE_PROVENANCE.json", {"source": "official_SpatialDM_0.3.1_bundled_CellChatDB_human",
               "interaction_file": str(ip), "interaction_sha256": digest(str(ip)),
               "complex_file": str(cp), "complex_sha256": digest(str(cp)), "package_version": VERSION})
    return inter, complexes


def eligibility(obj, inter, complexes, corrected):
    def components(value):
        return ([str(x) for x in complexes.loc[value].dropna().tolist() if str(x)]
                if value in complexes.index else [str(value)])
    positive = np.asarray((corrected > 0).sum(axis=0)).ravel()
    expression = dict(zip(obj.var_names, positive))
    rows = []
    for pair_id, row in inter.iterrows():
        lig, rec = components(row.ligand), components(row.receptor)
        missing = [g for g in lig + rec if g not in expression]
        low = [g for g in lig + rec if expression.get(g, 0) < 3]
        variance_ok = True
        if not missing and not low:
            for group in (lig, rec):
                values = np.asarray(obj[:, group].X.mean(axis=1)).ravel()
                if np.std(values) == 0 or not np.isfinite(values).all():
                    variance_ok = False
        reason = ("missing_required_complex_subunit" if missing else
                  "required_subunit_detected_in_fewer_than3spots" if low else
                  "constant_or_nonfinite_complex_expression" if not variance_ok else "eligible")
        rows.append({"pair_id": pair_id, "ligand": row.ligand, "receptor": row.receptor,
                     "ligand_components": ";".join(lig), "receptor_components": ";".join(rec),
                     "complete_complex": not missing, "minimum_subunit_positive_spots": min([expression.get(g, 0) for g in lig + rec]),
                     "missing_components": ";".join(missing), "low_detection_components": ";".join(low),
                     "eligible": reason == "eligible", "test_status": reason,
                     "annotation": row.annotation, "pathway_name": row.get("pathway_name", "")})
    return pd.DataFrame(rows).set_index("pair_id", drop=False)


def weights(obj, radius, out, label, fixed_primary_contact=None):
    cutoff = .2
    lengthscale = radius / np.sqrt(-2 * np.log(cutoff))
    # Generate enough official kNN candidates for the registered distance caps.
    # The release API applies cutoff only to secreted weights, not contact weights.
    sdm.weight_matrix(obj, l=lengthscale, cutoff=cutoff,
                      n_neighbors=min(186, obj.n_obs), n_nearest_neighbors=min(186, obj.n_obs), single_cell=False)
    contact = obj.obsp["nearest_neighbors"].tocoo()
    xy = np.asarray(obj.obsm["spatial"])
    distance = np.linalg.norm(xy[contact.row] - xy[contact.col], axis=1)
    keep = distance <= 1.2 + 1e-8
    contact = sparse.csr_matrix((contact.data[keep], (contact.row[keep], contact.col[keep])), shape=contact.shape)
    contact.eliminate_zeros()
    contact *= obj.n_obs / contact.sum()
    if fixed_primary_contact is not None:
        # This is a SECRETED-distance sensitivity; preserve primary contact weights exactly.
        contact = fixed_primary_contact.copy()
    obj.obsp["nearest_neighbors"] = contact
    audit = []
    for name in ("weight", "nearest_neighbors"):
        w = obj.obsp[name].tocsr()
        w.eliminate_zeros()
        coo = w.tocoo()
        ds = np.linalg.norm(xy[coo.row] - xy[coo.col], axis=1)
        if not np.isfinite(w.data).all() or (w.data < 0).any():
            raise ValueError("Invalid weight matrix")
        sparse.save_npz(out / f"{label}_{name}.npz", w)
        audit.append({"matrix": name, "n_spots": obj.n_obs, "nnz": w.nnz,
                      "sum_weights": float(w.sum()), "diagonal_nonzero": int((w.diagonal() != 0).sum()),
                      "maximum_edge_distance_spacing": float(ds.max()),
                      "median_degree_including_diagonal": float(np.median(np.diff(w.indptr))),
                      "asymmetry_max": float(abs(w - w.T).max()),
                      "lengthscale": float(2.0 / np.sqrt(-2 * np.log(cutoff)) if name == "nearest_neighbors" and fixed_primary_contact is not None else lengthscale),
                      "secreted_cutoff": cutoff,
                      "distance_cap_spacing": radius if name == "weight" else 1.2,
                      "weight_definition": "official_RBF_and_global_sum_rescaling; contact_further_geometric_cap_then_same_rescaling"})
    pd.DataFrame(audit).to_csv(out / f"{label}_WEIGHT_AUDIT.tsv", sep="\t", index=False)


def results_table(obj, elig, donor, label):
    x = obj.uns["global_res"].copy()
    x["pair_id"] = x.index.astype(str)
    x["pvalue"] = x.perm_pval
    x["padj"] = bh(x.perm_pval)
    x["analytic_pvalue"] = x.z_pval
    x["analytic_padj"] = bh(x.z_pval)
    x["significant"] = x.padj < .05
    x["global_I"] = obj.uns["global_I"]
    x["donor"] = donor
    x["n_spots"] = obj.n_obs
    x["primary_method"] = "official_permutation_1000"
    x["geometry"] = label
    x["pvalue_resolution"] = .001
    x["zero_pvalue_note"] = "official_empirical_zero_means_no_null_exceedance_not_probability_zero"
    return x.reset_index(drop=True).merge(elig.drop(columns=["pair_id"]), left_on="pair_id", right_index=True, validate="one_to_one")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--donors", nargs="+", choices=DONORS, default=list(DONORS))
    args = ap.parse_args()
    root = args.root.resolve()
    if VERSION != "0.3.1":
        raise ValueError("Stable locked SpatialDM0.3.1 required")
    pre = json.loads((root / "data/metadata/SPATIAL_PREPROCESS_ALL_DONORS_STATUS.json").read_text())
    if pre.get("status") != "PREPROCESS_COMPLETE" or pre.get("n_donors") != 8:
        raise ValueError("All eight new author preprocessing inputs must complete first")
    outroot = root / "results/spatial/SpatialDM"
    outroot.mkdir(parents=True, exist_ok=True)
    protocol = outroot / "FROZEN_MODEL_PROTOCOL.json"
    if protocol.exists() and json.loads(protocol.read_text()) != PROTOCOL:
        raise ValueError("SpatialDM frozen protocol changed")
    write_json(protocol, PROTOCOL)
    protocol_sha = digest(protocol)
    inter, complexes = database(outroot)
    for donor in args.donors:
        out = outroot / donor
        out.mkdir(exist_ok=True)
        lock = out / "worker.lock"
        lock.mkdir()  # fail closed on concurrent or incompletely audited prior worker
        try:
            ip = root / "data/processed/spatial" / donor / "author_reimplementation_primary.h5ad"
            input_sha = digest(ip)
            expected = {"input_sha256": input_sha, "protocol_sha256": protocol_sha, "version": VERSION, "donor": donor}
            final = out / "SPATIALDM_STATUS.json"
            if final.exists():
                previous = json.loads(final.read_text())
                if previous.get("status") != "SPATIALDM_COMPLETE" or any(previous.get(k) != v for k, v in expected.items()):
                    raise ValueError("Existing NEW-project final status differs from current inputs")
                if digest(out / "global_LR.tsv") != previous["global_table_sha256"]:
                    raise ValueError("Existing final table checksum differs")
                print(now(), "REUSE_VERIFIED_NEW_PROJECT_COMPLETE", donor, flush=True)
                continue
            start = now()
            source = ad.read_h5ad(ip)
            raw = source.layers["counts_raw_integer_UMI"]
            if not np.equal(raw.data, np.floor(raw.data)).all():
                raise ValueError("Raw integer layer violated")
            corrected, symbols = aggregate_symbol(source.layers["counts_author_background_corrected"], source.var.gene_symbol.astype(str))
            obj = ad.AnnData(normalized_log(corrected), obs=source.obs.copy(), var=pd.DataFrame(index=symbols))
            xy = np.asarray(source.obsm["spatial"], dtype=float)
            ds, _ = cKDTree(xy).query(xy, k=2)
            spacing = float(np.median(ds[:, 1]))
            if spacing <= 0 or not np.isfinite(xy).all():
                raise ValueError("Invalid coordinates/nearest spacing")
            obj.obsm["spatial"] = xy / spacing
            pd.DataFrame({"barcode": obj.obs_names, "pixel_x": xy[:, 0], "pixel_y": xy[:, 1],
                          "spacing_x": xy[:, 0] / spacing, "spacing_y": xy[:, 1] / spacing}).to_csv(out / "coordinates.tsv", sep="\t", index=False)
            elig = eligibility(obj, inter, complexes, corrected)
            elig.assign(donor=donor).to_csv(out / "ALL_LR_ELIGIBILITY.tsv", sep="\t", index=False)
            with timing(out, donor, "extract_lr", input_sha, {"official_API": "extract_lr", "species": "human", "min_cell": 3, "datahost": "package"}):
                sdm.extract_lr(obj, "human", mean="algebra", min_cell=3, datahost="package")
            ids = [x for x in obj.uns["geneInter"].index if bool(elig.loc[x, "eligible"])]
            if set(ids) != set(elig.index[elig.eligible]):
                raise ValueError("Official extracted pairs differ from independent eligible audit")
            for key in ("geneInter", "ligand", "receptor"):
                obj.uns[key] = obj.uns[key].loc[ids]
            obj.uns["num_pairs"] = len(ids)
            seed = 20261105 + 100 * DONORS.index(donor) + 50
            gp, gs = out / "global_model.pkl.gz", out / "GLOBAL_CHECKPOINT.json"
            if gp.exists() and gs.exists():
                checkpoint = json.loads(gs.read_text())
                if checkpoint["expected"] != expected or digest(gp) != checkpoint["model_sha256"]:
                    raise ValueError("Global checkpoint input/model checksum mismatch")
                obj = read_model(gp)
                audit_global_model(obj, out, "primary_global", seed)
                print(now(), "REUSE_VERIFIED_NEW_GLOBAL_CHECKPOINT", donor, flush=True)
            else:
                weights(obj, 2.0, out, "primary")
                np.random.seed(seed); random.seed(seed)
                with timing(out, donor, "spatialdm_global", input_sha, {"n_perm": 1000, "method": "both", "nproc": 1, "seed": seed, "n_pairs": len(ids)}):
                    sdm.spatialdm_global(obj, n_perm=1000, specified_ind=None, method="both", nproc=1)
                    sdm.sig_pairs(obj, method="permutation", fdr=True, threshold=.05)
                audit_global_model(obj, out, "primary_global", seed)
                dump_model(gp, obj)
                write_json(gs, {"expected": expected, "model_sha256": digest(gp), "actual_n_perm": 1000, "completed_at": now()})
            global_table = results_table(obj, elig, donor, "primary_secreted2_contact1.2")
            global_table.to_csv(out / "global_LR.tsv", sep="\t", index=False)
            global_table[["pair_id", "analytic_pvalue", "analytic_padj", "z", "donor"]].to_csv(out / "global_LR_analytic_supplement.tsv", sep="\t", index=False)
            selected = global_table[global_table.padj.lt(.05)].sort_values(["padj", "pvalue", "pair_id"], kind="stable").head(100).pair_id.tolist()
            local_plan = global_table[["pair_id", "pvalue", "padj"]].copy()
            local_plan["selected_for_local"] = local_plan.pair_id.isin(selected)
            local_plan["local_status"] = np.where(local_plan.selected_for_local, "selected_by_prefixed_global_rule", np.where(local_plan.padj.lt(.05), "global_significant_beyond_fixed100cap", "global_not_significant"))
            local_plan.to_csv(out / "LOCAL_SELECTION_FROZEN.tsv", sep="\t", index=False)
            lp, ls = out / "local_model.pkl.gz", out / "LOCAL_CHECKPOINT.json"
            nlocal = 0
            if selected:
                if lp.exists() and ls.exists():
                    checkpoint = json.loads(ls.read_text())
                    if checkpoint["expected"] != expected or checkpoint["selected_pairs"] != selected or digest(lp) != checkpoint["model_sha256"]:
                        raise ValueError("Local checkpoint differs")
                    local = read_model(lp)
                    audit_local_model(local, out, seed + 1)
                else:
                    local = obj.copy()
                    # Official per-gene max and SD scaling is independent of other genes.
                    # Restrict .raw to actually used complete subunits to avoid a whole-gene dense allocation.
                    required = sorted(set(g for pair in selected for side in ("ligand_components", "receptor_components") for g in str(elig.loc[pair, side]).split(";")))
                    jj = obj.var_names.get_indexer(required)
                    if (jj < 0).any():
                        raise ValueError("Selected local complex missing gene")
                    local_input, corrected_all_gene_totals = local_cp10000(corrected, jj)
                    local.raw = ad.AnnData(local_input, obs=local.obs.copy(), var=pd.DataFrame(index=required))
                    local.uns["LOCAL_INPUT_LAYER_CONTRACT"] = {
                        "AnnData_raw_slot": "library_normalized_background_corrected_nonlog_CP10000; NOT_raw_UMI",
                        "normalization_denominator": "sum_author_background_corrected_expression_over_ALL_genes_per_spot_before_subsetting",
                        "required_gene_subset_only_for_memory": True,
                        "official_local_further_scaling": "per_gene_max_then_scanpy_scale_zero_center_False",
                        "adaptation_scope": "project_specific_author_mat_norm_alignment; not_official_tutorial_rawcount_default"}
                    pd.DataFrame({"barcode": local.obs_names,
                                  "background_corrected_total_ALL_genes": corrected_all_gene_totals,
                                  "CP10000_scaling_factor": 10000.0 / corrected_all_gene_totals,
                                  "selected_required_component_genes_n": len(required),
                                  "selected_component_sum_CP10000": np.asarray(local_input.sum(axis=1)).ravel()}).to_csv(out / "LOCAL_INPUT_NORMALIZATION.tsv", sep="\t", index=False)
                    # Preserve official category ordering: local's implementation splits contact/secreted by contiguous count.
                    selected_official_order = [x for x in obj.uns["geneInter"].index if x in set(selected)]
                    np.random.seed(seed + 1); random.seed(seed + 1)
                    with timing(out, donor, "spatialdm_local", input_sha, {"n_perm": 1000, "method": "both", "nproc": 1, "seed": seed + 1, "selected_pairs": selected_official_order, "scale_X": True}):
                        sdm.spatialdm_local(local, n_perm=1000, method="both", specified_ind=selected_official_order, nproc=1, scale_X=True)
                        local.uns["unadjusted_local_perm_p_saved_before_sig_spots"] = local.uns["local_perm_p"].copy()
                        sdm.sig_spots(local, method="permutation", fdr=True, threshold=.05)
                    audit_local_model(local, out, seed + 1)
                    dump_model(lp, local)
                    write_json(ls, {"expected": expected, "model_sha256": digest(lp), "selected_pairs": selected, "actual_n_perm": 1000, "completed_at": now()})
                p = local.uns["unadjusted_local_perm_p_saved_before_sig_spots"]
                q = local.uns["local_stat"]["local_fdr"]
                z = local.uns["local_z_p"]
                local_rows = []
                for pair in p.index:
                    t = pd.DataFrame({"donor": donor, "pair_id": pair, "barcode": p.columns,
                                      "permutation_pvalue": p.loc[pair].to_numpy(float),
                                      "permutation_padj": q.loc[pair].to_numpy(float),
                                      "analytic_pvalue": z.loc[pair].to_numpy(float),
                                      "significant": q.loc[pair].to_numpy(float) < .05,
                                      "local_I": local.uns["local_stat"]["local_I"][:, list(p.index).index(pair)],
                                      "local_I_R": local.uns["local_stat"]["local_I_R"][:, list(p.index).index(pair)]})
                    t["pixel_x"], t["pixel_y"] = xy[:, 0], xy[:, 1]
                    local_rows.append(t)
                local_table = pd.concat(local_rows, ignore_index=True)
                local_table.to_csv(out / "local_spot_statistics.tsv.gz", sep="\t", index=False, compression="gzip")
                nlocal = int(local_table.significant.sum())
                del local
            else:
                pd.DataFrame(columns=["donor", "pair_id", "barcode", "permutation_pvalue", "permutation_padj", "analytic_pvalue", "significant", "pixel_x", "pixel_y"]).to_csv(out / "local_spot_statistics.tsv.gz", sep="\t", index=False, compression="gzip")
            # Fixed distance sensitivity: same complete global family, own actual 1000 permutations.
            sp, ss = out / "sensitivity_radius1.5_model.pkl.gz", out / "SENSITIVITY_CHECKPOINT.json"
            if sp.exists() and ss.exists():
                ck = json.loads(ss.read_text())
                if ck["expected"] != expected or digest(sp) != ck["model_sha256"]:
                    raise ValueError("Geometry sensitivity checkpoint differs")
                sens = read_model(sp)
                audit_global_model(sens, out, "geometry_radius1.5_global", seed + 2)
            else:
                sens = obj.copy()
                weights(sens, 1.5, out, "sensitivity_radius1.5", fixed_primary_contact=obj.obsp["nearest_neighbors"])
                np.random.seed(seed + 2); random.seed(seed + 2)
                with timing(out, donor, "global_geometry_sensitivity", input_sha, {"n_perm": 1000, "method": "both", "seed": seed + 2, "secreted_radius": 1.5, "contact_radius": 1.2}):
                    sdm.spatialdm_global(sens, n_perm=1000, method="both", nproc=1)
                    sdm.sig_pairs(sens, method="permutation", fdr=True, threshold=.05)
                audit_global_model(sens, out, "geometry_radius1.5_global", seed + 2)
                dump_model(sp, sens)
                write_json(ss, {"expected": expected, "model_sha256": digest(sp), "completed_at": now()})
            results_table(sens, elig, donor, "sensitivity_secreted1.5_contact1.2").to_csv(out / "global_LR_geometry_sensitivity_radius1.5.tsv", sep="\t", index=False)
            write_json(final, {**expected, "status": "SPATIALDM_COMPLETE", "n_spots": obj.n_obs,
                       "n_eligible_complete_pairs": len(ids), "n_global_pairs_tested": len(global_table),
                       "n_global_permutation_BH_significant": int(global_table.significant.sum()), "n_local_pairs": len(selected),
                       "n_local_pair_spot_cells_significant": nlocal, "n_perm_global": 1000,
                       "n_perm_local": 1000 if selected else None, "n_perm_local_requested": 1000,
                       "local_execution_status": "OFFICIAL_LOCAL_COMPLETE" if selected else "NOT_RUN_NO_GLOBAL_SIGNIFICANT_PAIRS",
                       "primary_P_type": "permutation", "analytic_P_separate": True, "local_cap": 100,
                       "median_nearest_spacing_pixels": spacing, "seed_global": seed, "seed_local": seed + 1,
                       "geometry_sensitivity_completed": True, "all_complexes_complete": True,
                       "geometry_sensitivity_cohort": "same_primary_author_QC_spots; secreted_radius_only",
                       "strict_z3_5_QC_mask_global_LR_sensitivity": "NOT_PERFORMED",
                       "strict_z3_5_QC_mask_local_LR_sensitivity": "NOT_PERFORMED",
                       "global_table_sha256": digest(out / "global_LR.tsv"), "is_new_fit": True,
                       "started_at": start, "finished_at": now(), "no_candidate_or_pathway_preselection": True})
            print(now(), "SPATIALDM_COMPLETE", donor, "global", len(global_table), "local", len(selected), flush=True)
        except Exception as e:
            write_json(out / "LAST_FAILURE.json", {"status": "FAILED", "at": now(), "error": repr(e), "traceback": traceback.format_exc()})
            raise
        finally:
            lock.rmdir()
    paths = [outroot / d / "SPATIALDM_STATUS.json" for d in DONORS]
    if all(p.exists() for p in paths):
        statuses = [json.loads(p.read_text()) for p in paths]
        if all(s.get("status") == "SPATIALDM_COMPLETE" and s.get("n_perm_global") == 1000 for s in statuses):
            pd.concat([pd.read_csv(outroot / d / "global_LR.tsv", sep="\t") for d in DONORS], ignore_index=True).to_csv(root / "results/spatial/global_LR_all_donors.tsv", sep="\t", index=False)
            pd.concat([pd.read_csv(outroot / d / "ALL_LR_ELIGIBILITY.tsv", sep="\t") for d in DONORS], ignore_index=True).to_csv(root / "results/spatial/ALL_LR_ELIGIBILITY_all_donors.tsv", sep="\t", index=False)
            write_json(outroot / "SPATIALDM_ALL_DONORS_STATUS.json", {"status": "SPATIALDM_COMPLETE", "n_donors": 8, "donors": statuses,
                       "finished_at": now(), "biological_replicate_unit": "donor_not_spot", "global_and_local_are_spatial_coexpression_not_direct_binding_proof": True})


if __name__ == "__main__":
    main()
