#!/usr/bin/env python3
"""Build evidence-bound final reports from the completed NEW project only.

This generator fits nothing and nominates no target. Its default never opens an
old-project path. Old comparison is a separately requested post-freeze audit.
Run normally through run_logged.py, then run --refresh-manifest-only once after
that logger has finished, without another project logger wrapping the refresh.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import time

PROJECT = "IBS_A_OFFICIAL_RERUN_20261105"
DONORS = ("P1", "P5", "P8", "P9", "P10", "P11", "P12", "P13")
CONTRASTS = ("Duo_IBSD_vs_HV", "Duo_IBSC_vs_HV", "Jej_IBSD_vs_HV", "Jej_IBSC_vs_HV")
REPORTS = ("Computational_Methods_EN.md", "Computational_Results_EN.md", "RESULTS_REVIEW_CN.md",
           "DATA_SOURCE_AND_QC_CN.md", "FIGURE_LEGENDS_EN.md", "LIMITATIONS_AND_FAILURES.md",
           "EXECUTION_TRACE.tsv", "FILES_MANIFEST.tsv", "PROJECT_REBUILD_DIFF_VS_OLD.md")


def args_parse():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--context", default="data/metadata/FINAL_REPORT_CONTEXT.json")
    p.add_argument("--freeze", default="results/integrated/CANDIDATE_FREEZE.json")
    p.add_argument("--ranking-rule", default="results/integrated/CANDIDATE_RANKING_RULE.json")
    p.add_argument("--figure-dir", default="figures/main")
    p.add_argument("--out", default="reports/final")
    p.add_argument("--old-comparison", action="store_true", help="Explicit post-new-freeze read-only audit; never an analysis input")
    p.add_argument("--old-context", help="New-root JSON containing exact permitted old files, SHA256 and metric selectors")
    p.add_argument("--refresh-manifest-only", action="store_true", help="After report logger FINISH: refresh trace snapshot and recursive manifest only; no model/report re-fit")
    return p.parse_args()


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for b in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def write_json(path, obj):
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def require(obj, key, where):
    if key not in obj or obj[key] is None:
        raise ValueError(f"Required field absent: {where}:{key}")
    return obj[key]


def tf(value):
    return str(value).strip().lower() in ("true", "1", "t")


def number(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def md(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return str(value).replace("|", "\\|").replace("\n", " ")


def md_table(records, fields):
    rows = ["| " + " | ".join(md(label) for _, label in fields) + " |", "| " + " | ".join("---" for _ in fields) + " |"]
    rows.extend("| " + " | ".join(md(row.get(key, "NA")) for key, _ in fields) + " |" for row in records)
    return "\n".join(rows) + "\n"


def bulk_qc_comparison(d, labels, bsummary):
    """Read the whole native summary; validate membership, never select favorable rows."""
    rel = "results/bulk/sensitivity/QC_SENSITIVITY_COMPARISON_SUMMARY.tsv"
    table = d.tab(rel, dtype=str, keep_default_na=False, na_filter=False)
    needed = {"sensitivity", "contrast", "n_features_union", "n_finite_LFC_both", "LFC_spearman",
              "LFC_sign_agreement_fraction", "n_main_FDR05", "n_sensitivity_FDR05", "n_both_FDR05",
              "n_main_significant_sensitivity_NA"}
    expected = {(label, contrast) for label in labels for contrast in CONTRASTS}
    if not needed.issubset(table.columns) or table.duplicated(["sensitivity", "contrast"]).any():
        raise ValueError("Native bulk QC comparison schema/unique cohort-contrast rows are incomplete")
    if len(table) != len(expected) or set(zip(table.sensitivity, table.contrast)) != expected:
        raise ValueError("All fixed bulk-QC cohort-by-contrast rows must be reported; no favorable subset")
    main_counts = bsummary.set_index("contrast")["n_DEG_FDR05"].to_dict()
    changed = False
    for row in table.to_dict("records"):
        counts = {}
        for field in ("n_features_union", "n_finite_LFC_both", "n_main_FDR05", "n_sensitivity_FDR05",
                      "n_both_FDR05", "n_main_significant_sensitivity_NA"):
            value = row[field]
            if not number(value) or float(value) < 0 or float(value) != int(float(value)):
                raise ValueError("Invalid native bulk sensitivity count: " + field)
            counts[field] = int(float(value))
        if counts["n_main_FDR05"] != int(main_counts[row["contrast"]]):
            raise ValueError("QC comparison main FDR count differs from the actual main contrast summary")
        if counts["n_both_FDR05"] > min(counts["n_main_FDR05"], counts["n_sensitivity_FDR05"]) or counts["n_main_significant_sensitivity_NA"] > counts["n_main_FDR05"]:
            raise ValueError("Native sensitivity intersection/NA count is inconsistent")
        changed |= (counts["n_main_FDR05"] != counts["n_sensitivity_FDR05"] or
                    counts["n_both_FDR05"] != counts["n_main_FDR05"] or counts["n_main_significant_sensitivity_NA"] > 0)
    fields = [(column, column) for column in table.columns]
    table_md = md_table(table.to_dict("records"), fields)
    en = ("Corrected-significance calls changed across the prespecified QC cohorts." if changed else
          "The native summary did not show a change in corrected-significance membership across these QC cohorts.")
    en += (" Continuous log2-fold-change correlation/sign agreement and corrected-significance overlap are different properties. "
           "A gain or loss of FDR significance is not evidence of whole-transcriptome direction reversal, and all primary DEGs are not thereby labelled robust. "
           "Native adjusted-P NA is unavailable, not a nonsignificant zero. The complete cohort-by-contrast table below is read verbatim as strings from `" + rel +
           "`; no thresholds, samples, P values, models or rankings were changed by reporting.\n\n")
    cn = ("本次校正显著性集合随预先固定的QC cohort改变。" if changed else
          "该原生汇总表中未见这些QC cohort之间校正显著性集合改变。")
    cn += ("连续log2FC相关/方向一致性与FDR显著性交集是不同性质；显著性增减不代表全转录组方向整体翻转，" 
           "也不能据此把主分析全部DEG称为稳健。原生q值NA仍为不可用，不能填成不显著或0。" 
           "以下逐项保留完整固定cohort×四比较的原始字符串值，来源为`" + rel +
           "`；报告不改变阈值、样本、P值、模型或排名。\n\n")
    return table, en + table_md, cn + table_md


def frozen_target_wide(d, frozen, continuous_gate, target_keys):
    """Pure row selection of the SHA-bound wide catalogue; all values stay strings."""
    import pandas as pd
    rel = require(continuous_gate, "wide_table", "continuous integration gate")
    expected = require(continuous_gate, "wide_table_sha256", "continuous integration gate")
    source = d.register(rel)
    if sha(source) != expected:
        raise ValueError("Complete continuous wide table differs from its frozen integration gate")
    ranking = d.js(require(frozen, "statistical_ranking_freeze_file", "candidate freeze"))
    if (ranking.get("status") != "FROZEN" or ranking.get("wide_table") != rel or ranking.get("wide_table_sha256") != expected or
            ranking.get("ranking_rule_sha256") != frozen["ranking_rule_sha256"] or
            ranking.get("source_freeze_sha256") != frozen["continuous_source_freeze_sha256"] or
            ranking.get("ranking_table_sha256") != frozen["statistical_ranking_sha256"]):
        raise ValueError("Candidate statistical-ranking freeze does not bind this complete continuous wide table")
    parts, columns, source_rows = [], None, 0
    for chunk in pd.read_csv(source, sep="\t", chunksize=10000, dtype=str, keep_default_na=False, na_filter=False):
        if columns is None:
            columns = list(chunk.columns)
        elif list(chunk.columns) != columns:
            raise ValueError("Complete wide source schema changed while streaming")
        if "identity_key" not in chunk:
            raise ValueError("Continuous wide table has no canonical identity key")
        source_rows += len(chunk)
        parts.append(chunk.loc[chunk.identity_key.isin(target_keys)].copy())
    if columns is None:
        columns = list(pd.read_csv(source, sep="\t", nrows=0, dtype=str, keep_default_na=False, na_filter=False).columns)
    selected = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)
    if selected.identity_key.duplicated().any() or set(selected.identity_key) != set(target_keys):
        raise ValueError("Frozen target identities must map to exactly one full wide source row each")
    return selected, {"source_path": rel, "source_sha256": expected, "source_rows": source_rows,
        "source_columns": len(columns), "all_original_columns": columns, "selected_HGNC_identities": list(selected.identity_key),
        "selected_rows": len(selected), "row_order": "original_source_order_not_candidate_rank_order",
        "all_source_values_read_as_strings": True, "numeric_recalculation_or_new_score": False,
        "no_target_reason": frozen.get("no_target_reason") if not target_keys else None}


def target_bulk_qc_view(wide, labels):
    """Readable projection; the separate target-wide export keeps every source column."""
    import pandas as pd
    fields = ["identity_key", "gene_symbol", "contrast", "sensitivity", "original_feature_id",
              "main_log2FC", "main_P", "main_padj", "main_P_NA_reason", "main_padj_NA_reason",
              "QC_log2FC", "QC_P", "QC_padj", "QC_P_NA_reason", "QC_padj_NA_reason"]
    records = []
    required = {"identity_key", "gene_symbol"}
    for contrast in CONTRASTS:
        required.add(contrast + "_gene_id")
        for field in ("log2FoldChange", "pvalue", "padj", "P_NA_reason", "padj_NA_reason"):
            required.add(contrast + "_" + field)
            required.update(contrast + "_" + label + "_" + field for label in labels)
    if not required.issubset(wide.columns):
        raise ValueError("Full target-wide catalogue lacks the registered primary/fixed-QC bulk columns")
    for source in wide.to_dict("records"):
        for contrast in CONTRASTS:
            for label in labels:
                row = {"identity_key": source["identity_key"], "gene_symbol": source["gene_symbol"],
                       "contrast": contrast, "sensitivity": label,
                       "original_feature_id": source.get(contrast + "_gene_id", "NA")}
                for output, native in (("log2FC", "log2FoldChange"), ("P", "pvalue"), ("padj", "padj"),
                                       ("P_NA_reason", "P_NA_reason"), ("padj_NA_reason", "padj_NA_reason")):
                    row["main_" + output] = source.get(contrast + "_" + native, "NA")
                    row["QC_" + output] = source.get(contrast + "_" + label + "_" + native, "NA")
                records.append(row)
    return pd.DataFrame(records, columns=fields)


def frozen_target_lr_scope(d, targets, wide, wide_provenance):
    """Report frozen targets' native self-pairs beside the distinct gene-component scope.

    All native cells stay strings. Only exact identity row selection and descriptive
    counts of existing native calls occur; no model, threshold, P adjustment or score.
    """
    import pandas as pd
    meta = ["report_identity_key", "report_gene_symbol", "report_verified_spatial_symbol",
            "report_source_path", "report_source_sha256", "report_source_native_columns_json",
            "report_scope_definition"]
    required = {"donor", "pair_id", "ligand", "receptor", "ligand_components", "receptor_components",
                "pvalue", "padj", "significant"}
    native_tables, native_columns, native_records = {}, [], []
    for donor in DONORS:
        rel = f"results/spatial/SpatialDM/{donor}/global_LR.tsv"
        table = d.tab(rel, dtype=str, keep_default_na=False, na_filter=False)
        if not required.issubset(table.columns) or any(c in table.columns for c in meta):
            raise ValueError("Native global LR schema is incomplete or collides with report metadata: " + rel)
        if not table.empty and set(table.donor) != {donor}:
            raise ValueError("Native global LR donor identity differs from its source path: " + rel)
        for column in table.columns:
            if column not in native_columns:
                native_columns.append(column)
        native_tables[donor] = (table, rel, d.files[rel]["sha256"])
    summary_columns = ["row_type", "identity_key", "gene_symbol", "verified_spatial_symbol", "donor",
        "any_component_global_pairs_tested_n", "any_component_global_pairs_supported_n",
        "any_component_global_measurement_status", "any_component_global_donors_tested_n",
        "any_component_global_donors_supported_n", "local_any_component_pairs_tested_n",
        "local_any_component_pairs_with_supported_spots_n", "local_any_component_measurement_status",
        "local_any_component_donors_tested_n", "local_any_component_donors_supported_n",
        "local_any_component_untested_donors_json", "homotypic_pair_rows_n",
        "homotypic_native_significant_true_rows_n", "homotypic_native_significant_false_rows_n",
        "homotypic_native_significance_unavailable_rows_n", "homotypic_measurement_status",
        "homotypic_donors_tested_n", "homotypic_donors_supported_n", "homotypic_donors_nonsignificant_n",
        "homotypic_donors_native_significance_unavailable_n", "homotypic_donors_no_matching_test_n",
        "homotypic_supported_donors_json", "homotypic_nonsignificant_donors_json",
        "homotypic_native_significance_unavailable_donors_json", "homotypic_no_matching_test_donors_json",
        "homotypic_native_P_q_significant_json", "homotypic_source_paths_json", "homotypic_source_SHA256_json",
        "gene_component_wide_source_path", "gene_component_wide_source_sha256", "scope_definition"]
    records, totals = [], []
    definition = ("ANY_COMPONENT uses the original HGNC-wide gene-involving-pair catalogue; "
        "SELF_PAIR requires exact ligand=receptor=ligand_components=receptor_components=one verified raw spatial symbol. "
        "Native significant calls only; no P/FDR recomputation, cross-pair P combination or direct-binding inference. "
        "LOCAL_ANY_COMPONENT is gene-involving catalogue evidence, not a local self-pair test.")

    def count(value):
        text = str(value)
        return int(text) if re.fullmatch(r"[0-9]+", text) else None

    def call(value):
        text = str(value).strip().lower()
        if text in ("true", "1", "t"):
            return "TRUE"
        if text in ("false", "0", "f"):
            return "FALSE"
        return "UNAVAILABLE"

    def total_counts(rows, field):
        values = [count(row[field]) for row in rows]
        return str(sum(values)) if all(v is not None for v in values) else "NA"

    def donors_positive(rows, field):
        values = [count(row[field]) for row in rows]
        return str(sum(v > 0 for v in values)) if all(v is not None for v in values) else "NA"

    for target in targets:
        identity, gene = target["identity_key"], target["gene_symbol"]
        symbols = target.get("spatial_symbols")
        if not isinstance(symbols, list) or any(not isinstance(s, str) or not s for s in symbols) or len(set(symbols)) != len(symbols):
            raise ValueError("Frozen spatial symbols must be an explicit unique verified list")
        hit = wide.loc[wide.identity_key.eq(identity)]
        if len(hit) != 1:
            raise ValueError("Exactly one frozen HGNC-wide source row is required for LR scope reporting")
        source = hit.iloc[0]
        wide_fields = ["LR_donors_tested_n", "LR_donors_with_supported_pairs_n",
                       "LR_total_tested_donor_pair_records_n", "LR_total_supported_donor_pair_records_n"]
        wide_fields += [f"LR_{donor}_{suffix}" for donor in DONORS for suffix in (
            "global_pairs_tested_n", "global_pairs_BH_supported_n", "global_measurement_status",
            "local_pairs_tested_n", "local_pairs_with_supported_spots_n")]
        if any(field not in source.index for field in wide_fields):
            raise ValueError("Original wide LR scope fields are incomplete; do not infer missing counts")
        for symbol in symbols or [""]:
            donor_rows = []
            for donor in DONORS:
                table, rel, checksum = native_tables[donor]
                matched = table.loc[table.ligand.eq(symbol) & table.receptor.eq(symbol) &
                                    table.ligand_components.eq(symbol) & table.receptor_components.eq(symbol)].copy() if symbol else table.iloc[:0].copy()
                states = [call(v) for v in matched.significant]
                n_true, n_false, n_unknown = states.count("TRUE"), states.count("FALSE"), states.count("UNAVAILABLE")
                if not symbol:
                    status = "NOT_APPLICABLE_NO_VERIFIED_SPATIAL_SYMBOL"
                elif matched.empty:
                    status = "NOT_APPLICABLE_NO_MATCHING_TESTED_SELF_PAIR"
                elif n_true:
                    status = "TESTED_WITH_NATIVE_SUPPORT"
                elif n_unknown:
                    status = "TESTED_NATIVE_SIGNIFICANCE_UNAVAILABLE"
                else:
                    status = "TESTED_NATIVE_NONSIGNIFICANT"
                for original in matched.to_dict("records"):
                    native_records.append({"report_identity_key": identity, "report_gene_symbol": gene,
                        "report_verified_spatial_symbol": symbol, "report_source_path": rel,
                        "report_source_sha256": checksum, "report_source_native_columns_json": json.dumps(list(table.columns)),
                        "report_scope_definition": "STRICT_SINGLE_COMPONENT_SELF_PAIR_NATIVE_ROW", **original})
                lt = source[f"LR_{donor}_local_pairs_tested_n"]
                ls = source[f"LR_{donor}_local_pairs_with_supported_spots_n"]
                local_tested, local_support = count(lt), count(ls)
                local_status = ("NA_NATIVE_STATISTIC" if local_tested is None or local_support is None else
                    "UNTESTED_NOT_A_NEGATIVE" if local_tested == 0 else
                    "TESTED_WITH_SUPPORTED_SPOTS" if local_support > 0 else "TESTED_NO_SUPPORTED_SPOTS")
                row = {column: "" for column in summary_columns}
                row.update(row_type="DONOR", identity_key=identity, gene_symbol=gene, verified_spatial_symbol=symbol, donor=donor,
                    any_component_global_pairs_tested_n=source[f"LR_{donor}_global_pairs_tested_n"],
                    any_component_global_pairs_supported_n=source[f"LR_{donor}_global_pairs_BH_supported_n"],
                    any_component_global_measurement_status=source[f"LR_{donor}_global_measurement_status"],
                    local_any_component_pairs_tested_n=lt, local_any_component_pairs_with_supported_spots_n=ls,
                    local_any_component_measurement_status=local_status,
                    homotypic_pair_rows_n=str(len(matched)), homotypic_native_significant_true_rows_n=str(n_true),
                    homotypic_native_significant_false_rows_n=str(n_false), homotypic_native_significance_unavailable_rows_n=str(n_unknown),
                    homotypic_measurement_status=status, homotypic_donors_tested_n=str(int(not matched.empty)),
                    homotypic_donors_supported_n=str(int(n_true > 0)),
                    homotypic_donors_nonsignificant_n=str(int(bool(states) and n_false == len(states))),
                    homotypic_donors_native_significance_unavailable_n=str(int(n_unknown > 0 and n_true == 0)),
                    homotypic_donors_no_matching_test_n=str(int(matched.empty)),
                    homotypic_native_P_q_significant_json=json.dumps(matched[["pair_id", "pvalue", "padj", "significant"]].to_dict("records"), ensure_ascii=False),
                    homotypic_source_paths_json=json.dumps([rel]), homotypic_source_SHA256_json=json.dumps({rel: checksum}),
                    gene_component_wide_source_path=wide_provenance["source_path"],
                    gene_component_wide_source_sha256=wide_provenance["source_sha256"], scope_definition=definition)
                records.append(row); donor_rows.append(row)
            total = {column: "" for column in summary_columns}
            total.update(row_type="EIGHT_DONOR_SUMMARY", identity_key=identity, gene_symbol=gene, verified_spatial_symbol=symbol,
                donor="ALL8", any_component_global_pairs_tested_n=source["LR_total_tested_donor_pair_records_n"],
                any_component_global_pairs_supported_n=source["LR_total_supported_donor_pair_records_n"],
                any_component_global_measurement_status="ORIGINAL_HGNC_CATALOGUE_COUNTS",
                any_component_global_donors_tested_n=source["LR_donors_tested_n"],
                any_component_global_donors_supported_n=source["LR_donors_with_supported_pairs_n"],
                local_any_component_pairs_tested_n=total_counts(donor_rows,"local_any_component_pairs_tested_n"),
                local_any_component_pairs_with_supported_spots_n=total_counts(donor_rows,"local_any_component_pairs_with_supported_spots_n"),
                local_any_component_measurement_status="DONOR_SCOPE_DETAIL_RETAINED_NOT_A_SELF_PAIR_LOCAL_TEST",
                local_any_component_donors_tested_n=donors_positive(donor_rows,"local_any_component_pairs_tested_n"),
                local_any_component_donors_supported_n=donors_positive(donor_rows,"local_any_component_pairs_with_supported_spots_n"),
                local_any_component_untested_donors_json=json.dumps([r["donor"] for r in donor_rows if r["local_any_component_measurement_status"] == "UNTESTED_NOT_A_NEGATIVE"]),
                homotypic_measurement_status="DESCRIPTIVE_NATIVE_CALL_COUNTS" if symbol else "NOT_APPLICABLE_NO_VERIFIED_SPATIAL_SYMBOL",
                homotypic_native_P_q_significant_json=json.dumps({r["donor"]: json.loads(r["homotypic_native_P_q_significant_json"]) for r in donor_rows}, ensure_ascii=False),
                homotypic_source_paths_json=json.dumps([native_tables[v][1] for v in DONORS]),
                homotypic_source_SHA256_json=json.dumps({native_tables[v][1]: native_tables[v][2] for v in DONORS}),
                gene_component_wide_source_path=wide_provenance["source_path"], gene_component_wide_source_sha256=wide_provenance["source_sha256"],
                scope_definition=definition)
            for field in ("homotypic_pair_rows_n", "homotypic_native_significant_true_rows_n", "homotypic_native_significant_false_rows_n",
                          "homotypic_native_significance_unavailable_rows_n", "homotypic_donors_tested_n", "homotypic_donors_supported_n",
                          "homotypic_donors_nonsignificant_n", "homotypic_donors_native_significance_unavailable_n", "homotypic_donors_no_matching_test_n"):
                total[field] = total_counts(donor_rows, field)
            for status_name, field in (("TESTED_WITH_NATIVE_SUPPORT", "homotypic_supported_donors_json"),
                ("TESTED_NATIVE_NONSIGNIFICANT", "homotypic_nonsignificant_donors_json"),
                ("TESTED_NATIVE_SIGNIFICANCE_UNAVAILABLE", "homotypic_native_significance_unavailable_donors_json")):
                total[field] = json.dumps([r["donor"] for r in donor_rows if r["homotypic_measurement_status"] == status_name])
            total["homotypic_no_matching_test_donors_json"] = json.dumps([r["donor"] for r in donor_rows if r["homotypic_measurement_status"].startswith("NOT_APPLICABLE")])
            records.append(total); totals.append(total)
    native = pd.DataFrame(native_records, columns=meta + native_columns).fillna("")
    summary = pd.DataFrame(records, columns=summary_columns)
    en = ("## Frozen-target LR scope: gene-involving pairs versus single-component self-pairs\n\n"
          "A gene-level donor supported by any pair containing that gene as a ligand/receptor component is a different scope from its same-symbol single-component self-pair. "
          "The former is copied from the SHA-bound HGNC-wide catalogue, and may include heterotypic pairs and verified aliases. "
          "The latter requires exact ligand, receptor, ligand_components and receptor_components equality to each frozen verified raw spatial symbol. "
          "Counts below use the native `significant` field, never a new P/FDR cutoff; native permutation and analytic statistics remain separate. "
          "No matching tested self-pair is NOT_APPLICABLE/untested, not a negative test; unavailable native calls are not counted as nonsignificant. "
          "The complete matching native rows, including every negative and NA field, are in `FROZEN_TARGET_HOMOTYPIC_LR_NATIVE.tsv.gz`; "
          "all eight donor scopes and source SHA256 are in `FROZEN_TARGET_LR_SCOPE_SUMMARY.tsv`. No P values were combined and no target or rank was changed.\n\n")
    cn = ("## 冻结靶点LR范围：任意涉及基因的配对与单组分同名自配对\n\n"
          "任意含该基因组分的配对有支持，不等于该基因同名自配对有支持。前者原样读取SHA绑定HGNC宽表，可能包含异名配对及已核实别名；"
          "后者逐个冻结verified spatial symbol，严格要求ligand、receptor、ligand_components、receptor_components四列均精确等于该单一原始符号。"
          "以下计数仅按原生significant字段，不另设P/FDR阈值，置换与解析统计保持分离；无匹配已检验自配对为NOT_APPLICABLE/未检验，不能当阴性；"
          "原生显著性不可用也不能计作不显著。完整原生自配对行（含全部阴性及NA字段）另存FROZEN_TARGET_HOMOTYPIC_LR_NATIVE.tsv.gz，"
          "八供者范围及来源SHA256另存FROZEN_TARGET_LR_SCOPE_SUMMARY.tsv；没有合并P值、重新选靶或改排名。\n\n")
    for row in totals:
        name = row["gene_symbol"] + " (" + row["identity_key"] + ")"
        raw = row["verified_spatial_symbol"] or "NO_VERIFIED_SPATIAL_SYMBOL"
        en += (f"**{name}; raw symbol {raw}.** Any-component global support was {row['any_component_global_donors_supported_n']}/{row['any_component_global_donors_tested_n']} tested donors in the original catalogue. "
               f"For the exact single-component {raw} | {raw} self-pair, {row['homotypic_donors_supported_n']}/{row['homotypic_donors_tested_n']} tested donors had native support; "
               f"supported donors={row['homotypic_supported_donors_json']}; {row['homotypic_donors_nonsignificant_n']} were tested-nonsignificant ({row['homotypic_nonsignificant_donors_json']}); "
               f"{row['homotypic_donors_native_significance_unavailable_n']} had unavailable native significance and {row['homotypic_donors_no_matching_test_n']} had no matching tested self-pair. "
               f"The separate gene-involving local scope contained {row['local_any_component_pairs_with_supported_spots_n']}/{row['local_any_component_pairs_tested_n']} supported/tested donor-pair records, "
               f"across {row['local_any_component_donors_tested_n']} tested donor(s); {row['local_any_component_donors_supported_n']} donor(s) had supported spots. "
               f"Locally untested donors={row['local_any_component_untested_donors_json']}. This local gene-component count does not validate the gene's self-pair.\n\n")
        cn += (f"**{name}；原始符号{raw}。** 原宽表任意涉及基因组分的global配对支持为{row['any_component_global_donors_supported_n']}/{row['any_component_global_donors_tested_n']}个已检验供者。"
               f"严格单组分{raw} | {raw}自配对仅{row['homotypic_donors_supported_n']}/{row['homotypic_donors_tested_n']}个已检验供者原生显著，"
               f"支持供者={row['homotypic_supported_donors_json']}；另有{row['homotypic_donors_nonsignificant_n']}个已检验不显著（{row['homotypic_nonsignificant_donors_json']}），"
               f"{row['homotypic_donors_native_significance_unavailable_n']}个原生显著性不可用，{row['homotypic_donors_no_matching_test_n']}个没有匹配已检验自配对。"
               f"另行定义的gene-involving local范围为{row['local_any_component_pairs_with_supported_spots_n']}/{row['local_any_component_pairs_tested_n']}个有支持/已检验donor-pair记录，"
               f"涉及{row['local_any_component_donors_tested_n']}个已检验供者、其中{row['local_any_component_donors_supported_n']}个有显著spots；"
               f"local未检验供者={row['local_any_component_untested_donors_json']}。这些基因组分local计数不是同名自配对的local验证。\n\n")
    detail_fields = [(c, c) for c in ("gene_symbol", "verified_spatial_symbol", "donor", "any_component_global_pairs_supported_n",
        "any_component_global_pairs_tested_n", "homotypic_measurement_status", "homotypic_native_P_q_significant_json",
        "local_any_component_pairs_tested_n", "local_any_component_pairs_with_supported_spots_n", "local_any_component_measurement_status")]
    if totals:
        details = md_table([r for r in records if r["row_type"] == "DONOR"], detail_fields)
        en += details; cn += details
    else:
        en += "The fresh freeze nominated no target. Both auxiliary exports have genuine headers and zero selected target rows; no target-specific positive or negative record was invented.\n\n"
        cn += "本次freeze未提名靶点；两个辅助导出均有真实表头且选中靶点行数为0，不制造靶点正负记录。\n\n"
    en += "Gene-component support, self-pair coexpression, local selection and direct molecular binding are distinct claims. None establishes directed signalling, IBS-specific disease change, a pain mechanism or treatment efficacy.\n"
    cn += "基因组分配对支持、同名自配对空间共表达、local选择与直接分子结合是不同主张；均不能据此证明定向信号、IBS特异变化、疼痛机制或疗效。\n"
    return native, summary, en, cn


def frozen_feasibility_text(d, frozen, targets):
    """Reproduce the separate frozen review; do not select or infer target biology."""
    if not targets:
        reason = require(frozen, "no_target_reason", "no-target freeze")
        return ("## Frozen candidate rationale and proposed follow-up\n\nNo main target or backup was nominated: " + md(reason) +
                ". No human-model or downstream-qPCR choice is fabricated for this no-target decision.\n",
                "## 冻结候选的选择理由及拟议后续\n\n本次没有选定主靶点或备选：" + md(reason) +
                "。不为无靶点结论制造人源模型或下游qPCR选择。\n")
    sources = require(frozen, "review_source_evidence", "candidate freeze")
    if not isinstance(sources, dict) or not sources:
        raise ValueError("Selected targets require the actual frozen source-cited feasibility review")
    for identifier, source in sources.items():
        path = require(source, "path", "frozen feasibility source")
        if sha(d.register(path)) != require(source, "sha256", "frozen feasibility source"):
            raise ValueError("Frozen feasibility source changed: " + identifier)
        urls = require(source, "urls", "frozen feasibility source")
        if not isinstance(urls, list) or not urls or not all(isinstance(u, str) and u.startswith("https://") for u in urls):
            raise ValueError("Frozen feasibility source needs its recorded HTTPS citations")
    def cited(claim):
        statement = require(claim, "statement", "frozen feasibility claim")
        ids = require(claim, "source_ids", "frozen feasibility claim")
        if not isinstance(ids, list) or not ids or not set(ids).issubset(sources):
            raise ValueError("Frozen feasibility claim lacks its declared source IDs")
        links = "; ".join(" / ".join(f"[{identifier}]({url})" for url in sources[identifier]["urls"]) for identifier in ids)
        return statement, links
    en = ("## Frozen candidate rationale and proposed follow-up\n\nThe following reproduces the source-cited post-ranking review already bound to the candidate freeze. "
          "It does not alter statistical ranks, promote a backup or infer a new mechanism. Statements and proposed transcript directions are retained as recorded; "
          "model/reagent access and expression caveats are part of those statements. These are proposed siRNA/qPCR checks, not performed experiments or evidence of pain causality or therapeutic efficacy.\n\n")
    cn = ("## 冻结候选的选择理由及拟议后续\n\n以下逐项展示候选freeze已绑定的、统计排名完成后的来源支持可行性审核。"
          "不改变统计排名、不提升备选、不新增机制推断；陈述及拟议下游方向按冻结记录保留，模型/试剂获取与表达限制也原样保留。"
          "这些仅为拟议siRNA/qPCR验证，不是已完成实验，也不是疼痛因果或治疗有效证据。\n\n")
    components = (("human_model", "Human model", "人源模型"), ("target_expression", "Target expression", "靶点表达"),
                  ("sirna_availability", "siRNA feasibility", "siRNA可行性"),
                  ("IBS_or_visceral_sensitivity_rationale", "IBS/visceral-sensitivity rationale", "IBS或内脏敏感性依据"))
    for index, target in enumerate(targets):
        role_en, role_cn = ("Primary", "主靶点") if index == 0 else ("Backup", "备选")
        name = md(require(target, "gene_symbol", "frozen target"))
        identity = md(require(target, "identity_key", "frozen target"))
        rank = md(require(target, "statistical_rank", "frozen target"))
        grade = md(require(target, "evidence_grade", "frozen target"))
        reason = md(require(target, "selection_reason", "frozen target"))
        en += f"### {role_en}: {name} ({identity})\n\nFrozen statistical rank: {rank}; evidence grade: {grade}. Selection reason: {reason}\n\n"
        cn += f"### {role_cn}：{name}（{identity}）\n\n冻结统计排名：{rank}；证据等级：{grade}。选择理由：{reason}\n\n"
        feasibility = require(target, "feasibility", "frozen target")
        rows = []
        for field, label, label_cn in components:
            statement, links = cited(require(feasibility, field, "frozen target feasibility"))
            rows.append({"component": label, "component_cn": label_cn, "statement": statement, "sources": links})
        en += md_table(rows, [("component", "Component"), ("statement", "Frozen statement"), ("sources", "Sources")]) + "\n"
        cn += md_table(rows, [("component_cn", "项目"), ("statement", "冻结陈述"), ("sources", "来源")]) + "\n"
        endpoints = require(feasibility, "qpcr_downstream", "frozen target feasibility")
        if not isinstance(endpoints, list) or not 2 <= len(endpoints) <= 4:
            raise ValueError("Selected target must retain two to four frozen downstream qPCR transcripts")
        symbols = [require(item, "gene_symbol", "frozen qPCR endpoint") for item in endpoints]
        if len(set(symbols)) != len(symbols) or target["gene_symbol"] in symbols:
            raise ValueError("Frozen downstream qPCR transcripts must be distinct and not the knockdown target")
        rows = []
        for endpoint in endpoints:
            statement, links = cited(endpoint)
            direction = require(endpoint, "direction", "frozen qPCR endpoint")
            if direction not in {"increase_expected", "decrease_expected", "uncertain", "context_dependent"}:
                raise ValueError("Frozen proposed qPCR direction must retain its explicit uncertainty category")
            rows.append({"gene_symbol": endpoint["gene_symbol"], "direction": direction, "statement": statement, "sources": links})
        en += "Proposed downstream qPCR transcripts (not measured in this project):\n\n" + md_table(rows, [
            ("gene_symbol", "Proposed downstream transcript"), ("direction", "Recorded direction"),
            ("statement", "Frozen rationale"), ("sources", "Sources")]) + "\n"
        cn += "拟议下游qPCR转录本（本项目尚未实验测量）：\n\n" + md_table(rows, [
            ("gene_symbol", "拟议下游转录本"), ("direction", "记录方向"), ("statement", "冻结依据"), ("sources", "来源")]) + "\n"
    return en, cn


class Inputs:
    def __init__(self, root):
        self.root, self.files = root, {}

    def path(self, rel, exists=True):
        p = self.root / rel
        if not p.resolve().is_relative_to(self.root) or p.is_symlink() or p.resolve() != p:
            raise ValueError(f"Not a real new-root path: {p}")
        if exists and (not p.is_file() or not p.stat().st_size):
            raise FileNotFoundError(p)
        return p

    def register(self, rel):
        p = self.path(rel)
        key = p.relative_to(self.root).as_posix()
        if key not in self.files:
            self.files[key] = {"sha256": sha(p), "bytes": p.stat().st_size}
        return p

    def js(self, rel):
        return json.loads(self.register(rel).read_text(encoding="utf-8-sig"))

    def tab(self, rel, **kwargs):
        import pandas as pd
        return pd.read_csv(self.register(rel), sep="\t", **kwargs)

    def text(self, rel):
        return self.register(rel).read_text(encoding="utf-8", errors="replace")


def canonical_trace(d):
    trace = d.tab("EXECUTION_TRACE.tsv", dtype=str, keep_default_na=False)
    needed = {"event", "module", "source_command", "input_hashes", "start_iso8601", "finish_iso8601",
              "duration_seconds", "exit_code", "hostname", "software_version", "run_id", "is_new_fit", "output_hashes", "notes"}
    if not needed.issubset(trace.columns):
        raise ValueError("Canonical trace lacks required actual-execution fields")
    finished = trace.loc[trace.event.eq("FINISH")].copy()
    if finished.empty:
        raise ValueError("No completed execution records")
    for row in finished.to_dict("records"):
        for field in ("start_iso8601", "finish_iso8601"):
            datetime.fromisoformat(row[field].replace("Z", "+00:00"))
        if not number(row["duration_seconds"]) or float(row["duration_seconds"]) < 0:
            raise ValueError("Invalid actual monotonic duration in canonical trace")
        int(row["exit_code"])
    return trace, finished


def native_trace(d, rel):
    import pandas as pd
    parser_dialect = "standard_tsv"
    # The frozen BayesSpace producer uses R write.table(qmethod="escape").
    # Its JSON source_command quotes are backslash-escaped, unlike Python CSV
    # and bulk fwrite traces. Decode only this exact native producer dialect;
    # retain/hash the original bytes and never rewrite any source trace.
    if re.fullmatch(r"results/spatial/BayesSpace/(?:P1|P5|P8|P9|P10|P11|P12|P13)/MODEL_EXECUTION_TRACE\.tsv", rel):
        parser_dialect = "R_write_table_qmethod_escape_known_BayesSpace_native_trace"
        with d.register(rel).open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", escapechar="\\", doublequote=False, strict=True)
            header = next(reader)
            rows = list(reader)
        expected = ["module", "source_command", "input_hashes", "start_iso8601", "finish_iso8601", "duration_seconds",
                    "exit_code", "hostname", "software_version", "run_id", "is_new_fit", "reused_only_immutable_raw_inputs", "notes"]
        if header != expected or any(len(row) != len(header) for row in rows):
            raise ValueError("Native R BayesSpace trace dialect/schema mismatch: " + rel)
        frame = pd.DataFrame(rows, columns=header)
        for command in frame.source_command:
            api, separator, parameters = command.partition(" ")
            if not separator or not api.startswith(("BayesSpace::", "scran::")) or not isinstance(json.loads(parameters), dict):
                raise ValueError("Native R source_command did not decode to its actual API/JSON parameters")
    else:
        frame = d.tab(rel, dtype=str, keep_default_na=False)
    for field in ("module", "source_command", "start_iso8601", "finish_iso8601", "duration_seconds", "exit_code", "software_version"):
        if field not in frame:
            raise ValueError(f"Native trace lacks {field}: {rel}")
    for row in frame.to_dict("records"):
        datetime.fromisoformat(row["start_iso8601"].replace("Z", "+00:00"))
        datetime.fromisoformat(row["finish_iso8601"].replace("Z", "+00:00"))
        if not number(row["duration_seconds"]) or float(row["duration_seconds"]) < 0:
            raise ValueError(f"Native model duration is unavailable: {rel}")
    frame["source_trace_path"] = rel
    frame["source_trace_row"] = range(2, len(frame) + 2)
    frame["source_trace_sha256"] = d.files[rel]["sha256"]
    frame["source_trace_parser_dialect"] = parser_dialect
    return frame


def versions(d, figure_status):
    # The exact native installed.packages lock is unquoted R write.table
    # output with field-internal literal newlines, not standard quoted TSV.
    def expect(condition, message):
        if not condition:
            raise ValueError(message)
    R_PACKAGE_LOCK_HEADER = ["Package", "LibPath", "Version", "Priority", "Depends", "Imports",
        "LinkingTo", "Suggests", "Enhances", "License", "License_is_FOSS", "License_restricts_use",
        "OS_type", "MD5sum", "NeedsCompilation", "Built"]

    def native_r_unquoted_package_lock_rows(text):
        """Read only the two known 16-column R locks without changing source bytes.

        R wrote field-internal newlines literally and tab-separated 16 fields.
        Accumulate physical lines until all 15 tab boundaries have been consumed;
        retain every internal newline/quote/backslash in its original cell. A
        surplus tab or an incomplete final record is an error, never truncated.
        """
        pieces = []
        for line in text:
            pieces.append(line)
            record = "".join(pieces)
            n_tabs = record.count("\t")
            expect(n_tabs <= 15, "Unexpected embedded/excess tab in native R package lock")
            if n_tabs < 15:
                continue
            # Remove exactly the logical record terminator, not internal newlines.
            if record.endswith("\r\n"):
                record = record[:-2]
            elif record.endswith(("\n", "\r")):
                record = record[:-1]
            yield record.split("\t")
            pieces = []
        expect(not pieces, "Incomplete final native R package-lock record")
    r_lock_rel = "env/R_FINAL_PACKAGE_LOCK.tsv"
    with d.register(r_lock_rel).open(encoding="utf-8-sig", newline="") as handle:
        records = list(native_r_unquoted_package_lock_rows(handle))
    expect(bool(records) and records[0] == R_PACKAGE_LOCK_HEADER,
           "Unexpected native installed.packages lock schema")
    by_name, seen_library_records = {}, set()
    for row in records[1:]:
        expect(len(row) == 16 and re.fullmatch(r"[A-Za-z][A-Za-z0-9.]*", row[0]) is not None
               and row[1].startswith("/") and not any(c in row[1] for c in "\r\n\t")
               and row[2] not in ("", "NA") and re.fullmatch(r"[^\s]+", row[2]) is not None
               and re.fullmatch(r"[0-9]+\.[0-9]+(?:\.[0-9]+)?", row[-1]) is not None,
               "Unexpected native R package-lock identity/version/Built fields")
        key = (row[0], row[1])
        expect(key not in seen_library_records, "Duplicate native R package/library record: " + str(key))
        seen_library_records.add(key)
        by_name.setdefault(row[0], []).append(row)
    project_library = str(d.root / "env/R_library")
    r, selected_libraries = {}, {}
    for name, candidates in by_name.items():
        preferred = [row for row in candidates if row[1] == project_library]
        expect(len(preferred) <= 1, "Ambiguous project-library R package version: " + name)
        chosen = preferred[0] if preferred else candidates[0]
        r[name], selected_libraries[name] = chosen[2], chosen[1]
    software_gate = d.js("manifest/SOFTWARE_GATE.json")
    for item in require(software_gate, "R_packages", "actual R software gate"):
        name = require(item, "package", "actual R software package")
        expect(name in r and selected_libraries[name] == project_library
               and require(item, "path", "actual R software package") == str(d.root / "env/R_library" / name)
               and r[name].replace("-", ".") == require(item, "version", "actual R software package").replace("-", "."),
               "Native R lock differs from actual isolated package-version evidence: " + name)
    for name in ("DESeq2", "fgsea", "msigdbr", "BayesSpace", "scran", "DirichletReg", "RcppDist", "Matrix"):
        require(r, name, "final R lock")
    # pip freeze retains PEP508 direct URLs for conda-built distributions.
    # Resolve their versions only from the already recorded actual distribution
    # ledger; do not infer versions from build paths or alter the original lock.
    lock_rel = "env/PYTHON_SOFTWARE_LOCK_FINAL.txt"
    ledger_rel = "env/PYTHON_DISTRIBUTIONS_AFTER_PKG_RESOURCES_FIX_20261105.json"
    ledger = d.js(ledger_rel)
    prefix = d.root / "env/conda"
    if ledger.get("prefix") != str(prefix) or ledger.get("python_executable") != str(prefix / "bin/python"):
        raise ValueError("Actual Python distribution ledger identifies a different environment")
    if str(ledger.get("python_runtime", "")).split()[0:1] != [figure_status.get("python")]:
        raise ValueError("Actual Python ledger runtime differs from the figure runtime")
    required_packages = {"numpy", "pandas", "scipy", "anndata", "matplotlib", "spatialdm"}
    python = require(ledger, "versions", "actual Python distribution ledger")
    installed = {}
    for item in require(ledger, "distributions", "actual Python distribution ledger"):
        key = re.sub(r"[-_.]+", "-", require(item, "name", "installed distribution")).lower()
        if key not in required_packages:
            continue
        version = require(item, "version", "installed distribution")
        location = Path(require(item, "location", "installed distribution"))
        if (key != item.get("normalized_name") or key in installed or
                (location != prefix and prefix not in location.parents) or ".." in location.parts):
            raise ValueError("Ambiguous or non-isolated actual Python distribution: " + key)
        installed[key] = version
    if not isinstance(python, dict) or installed != {k: python.get(k) for k in required_packages}:
        raise ValueError("Actual Python distribution ledger version maps disagree")
    lock_names, direct_url_names = set(), []
    for raw_line in d.text(lock_rel).splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        match = re.fullmatch(r"([A-Za-z0-9][A-Za-z0-9._-]*)(?:==([^\s]+)|\s+@\s+(\S+))", line)
        if match is None:
            raise ValueError("Unsupported final Python lock format: " + line)
        key = re.sub(r"[-_.]+", "-", match[1]).lower()
        if key not in required_packages:
            continue
        if key in lock_names or key not in installed:
            raise ValueError("Duplicate or absent actual Python lock distribution: " + key)
        lock_names.add(key)
        if match[2] is not None and match[2] != installed[key]:
            raise ValueError("Pinned and actual Python package versions disagree: " + key)
        if match[3] is not None:
            direct_url_names.append(key)
    if lock_names != required_packages:
        raise ValueError("Required actual Python packages absent from original pip lock")
    for name in ("numpy", "pandas", "scipy", "anndata", "matplotlib", "spatialdm"):
        require(python, name, "final Python lock")
    session = d.text("env/R_FINAL_SESSIONINFO.txt")
    rline = next((line for line in session.splitlines() if line.startswith("R version ")), None)
    if rline is None:
        raise ValueError("Final sessionInfo does not identify the actual R runtime")
    for name in ("numpy", "pandas", "scipy", "matplotlib"):
        if figure_status.get("package_versions", {}).get(name) != python[name]:
            raise ValueError("Actual figure package version differs from final Python lock: " + name)
    return {"R_runtime": rline, "R_packages": r,
            "R_version_source": {"path": r_lock_rel, "sha256": d.files[r_lock_rel]["sha256"],
                "parser": "R_INSTALLED_PACKAGES_UNQUOTED_MULTILINE", "logical_package_rows": len(records) - 1,
                "unique_package_names": len(by_name), "selected_libraries": selected_libraries,
                "version_selection": "Explicit project R_library takes precedence; otherwise first native lock record. This is installed-package inventory, not a claim that every package was loaded.",
                "all_library_records_retained_in_native_source": True,
                "columns": R_PACKAGE_LOCK_HEADER, "native_16_fields_and_internal_newlines_preserved": True,
                "native_lock_bytes_unchanged": True},
            "Python_packages": python, "Python_runtime": require(figure_status, "python", "figure status"),
            "Python_version_sources": {"original_pip_freeze": lock_rel, "actual_installed_distribution_ledger": ledger_rel,
                "ledger_sha256": d.files[ledger_rel]["sha256"], "direct_url_versions_resolved_from_actual_ledger": direct_url_names,
                "original_lock_unchanged": True, "version_inference_from_build_paths": False}}


def clinical_table(d, context):
    frame = d.tab(require(context, "clinical_metadata_path", "report context"), dtype=str, keep_default_na=False)
    if "clinical_donor_column" in context:
        donor_column = context["clinical_donor_column"]
        if donor_column not in frame:
            raise ValueError("Explicit clinical donor column absent")
    else:
        candidates = [c for c in frame if set(DONORS).issubset(set(frame[c].astype(str)))]
        if len(candidates) != 1:
            raise ValueError("Cannot uniquely identify actual clinical donor column; specify clinical_donor_column from the official metadata")
        donor_column = candidates[0]
    selected = frame.loc[frame[donor_column].isin(DONORS)].copy()
    if len(selected) != 8 or selected[donor_column].duplicated().any():
        raise ValueError("Official clinical metadata does not identify one row for each of the eight actual donors")
    selected = selected.set_index(donor_column).loc[list(DONORS)].reset_index()
    return selected, donor_column


def old_comparison(d, args, new_frozen):
    """The only function allowed to open old files, after explicit new freeze."""
    if not args.old_comparison:
        return ("# New/old descriptive comparison\n\nNOT_RUN: `--old-comparison` was not requested. No old-project statistical file, model, figure or ranking was opened by this report generator. The new model and candidate freezes are recorded independently.\n", {})
    if not new_frozen or not args.old_context:
        raise ValueError("Old comparison requires all new model/candidate freezes and explicit --old-context")
    context = d.js(args.old_context)
    if context.get("authorised_read_only_comparison") is not True:
        raise ValueError("Old comparison context does not explicitly authorise read-only comparison")
    old_root = Path(require(context, "old_root", "old context")).resolve()
    if old_root == d.root or old_root.is_relative_to(d.root) or not old_root.is_dir():
        raise ValueError("Invalid explicitly supplied old archive root")
    imports = {}
    for item in require(context, "files", "old context"):
        p = old_root / require(item, "path", "old file")
        if p.is_symlink() or p.resolve() != p or not p.resolve().is_relative_to(old_root) or not p.is_file():
            raise ValueError("Old comparison file escapes its explicitly authorised archive")
        if p.suffix.lower() not in (".json", ".tsv", ".csv", ".txt"):
            raise ValueError("Old comparison accepts statistical text only, never executable model objects")
        observed = sha(p)
        if observed != require(item, "sha256", "old file"):
            raise ValueError("Old comparison checksum differs from explicitly permitted file")
        imports[item["path"]] = {"path": p, "sha256": observed}
    def value(selector, old):
        import pandas as pd
        rel, kind = require(selector, "file", "comparison selector"), require(selector, "kind", "comparison selector")
        p = imports[rel]["path"] if old else d.register(rel)
        if kind == "json_pointer":
            obj = json.loads(p.read_text(encoding="utf-8-sig"))
            for part in require(selector, "pointer", "comparison selector").strip("/").split("/"):
                obj = obj[int(part)] if isinstance(obj, list) else obj[part.replace("~1", "/").replace("~0", "~")]
            return json.dumps(obj, ensure_ascii=False) if isinstance(obj, (dict, list)) else obj
        if kind in ("table_rows", "table_lookup"):
            table = pd.read_csv(p, sep="," if p.suffix == ".csv" else "\t", dtype=str, keep_default_na=False)
            if kind == "table_rows":
                return len(table)
            for key, expected in require(selector, "where", "comparison selector").items():
                table = table.loc[table[key].eq(str(expected))]
            column = require(selector, "column", "comparison selector")
            return "NO_MATCH" if table.empty else table.iloc[0][column] if len(table) == 1 else "MULTIPLE_MATCHES:" + json.dumps(table[column].tolist(), ensure_ascii=False)
        raise ValueError("Unsupported explicit comparison selector")
    records = [{"metric": require(c, "metric", "comparison"), "old": value(c["old"], True), "new": value(c["new"], False)} for c in require(context, "comparisons", "old context")]
    text = "# New/old descriptive comparison\n\nThe new models, continuous source statistics, ranking rule and candidate decision had already been frozen before the following old statistical text files were read. These comparisons did not correct, select or re-fit any new result. A changed count or rank is descriptive, not proof that one analysis is more valid.\n\n"
    text += md_table(records, [("metric", "Metric"), ("old", "Old archive"), ("new", "New frozen run")])
    text += "\nExplicit old input checksums:\n\n" + md_table([{"path": k, "sha256": v["sha256"]} for k, v in imports.items()], [("path", "Old relative path"), ("sha256", "SHA256")])
    for item in imports.values():
        if sha(item["path"]) != item["sha256"]:
            raise ValueError("Read-only old comparison source changed during comparison")
    return text, {k: {"absolute_path": str(v["path"]), "sha256": v["sha256"]} for k, v in imports.items()}


def manifest(root, out, trace, known_hashes, final):
    """Inventory every file; exclusion affects review packaging, not coverage."""
    manifest_path = out / "FILES_MANIFEST.tsv"
    self_metadata = {manifest_path.relative_to(root).as_posix(), (out / "FINAL_REPORTS_STATUS.json").relative_to(root).as_posix(), "manifest/FINAL_REPORTS_GATE.json"}
    active = set(trace.loc[trace.event.eq("START"), "run_id"]) - set(trace.loc[trace.event.eq("FINISH"), "run_id"])
    mutable_logs = set()
    for r in trace.loc[trace.run_id.isin(active)].to_dict("records"):
        mutable_logs.update(r.get(k, "") for k in ("stdout_log", "stderr_log"))
    records = []
    for directory, dirs, files in os.walk(root, followlinks=False):
        for name in [n for n in dirs if (Path(directory) / n).is_symlink()]:
            p = Path(directory) / name
            records.append({"relative_path": p.relative_to(root).as_posix() + "/", "bytes": p.lstat().st_size,
                            "sha256": "NA", "sha256_basis": "SYMLINK_DIRECTORY_TARGET_NOT_READ",
                            "review_package_policy": "EXCLUDED_SYMLINK_DIRECTORY_NO_TARGET_READ", "retained_on_server": True})
        dirs[:] = [n for n in dirs if not (Path(directory) / n).is_symlink()]
        dirs[:] = sorted(dirs)
        for name in sorted(files):
            p = Path(directory) / name
            rel = p.relative_to(root).as_posix()
            size = p.lstat().st_size
            reason = ""
            if p.is_symlink():
                reason = "EXCLUDED_SYMLINK_NO_TARGET_READ"
            elif rel.startswith("data/raw/"):
                reason = "EXCLUDED_RAW_DATA_RETAINED_ON_SERVER"
            elif rel.startswith(("env/conda/", "env/R_library/", "env/R_empty_site/")):
                reason = "EXCLUDED_SOFTWARE_INSTALLATION"
            elif re.search(r"\.(rds|rda|rdata|h5|h5ad|hdf5|mtx|npy|npz|pkl|pickle)(\.gz)?$", rel, re.I) or rel.endswith("phase0.sqlite") or ".sqlite-" in rel:
                reason = "EXCLUDED_MATRIX_MODEL_OR_DATABASE_RETAINED_ON_SERVER"
            elif re.search(r"\.(zip|7z|whl|tar|tgz|tar\.gz|tar\.bz2|tar\.xz|deb|rpm)$", rel, re.I):
                reason = "EXCLUDED_INSTALLER_OR_ORIGINAL_ARCHIVE"
            elif rel.startswith("data/reference/gwas/g1000_eur"):
                reason = "EXCLUDED_LD_REFERENCE_RETAINED_ON_SERVER"
            elif size > 64 * 1024 * 1024 and not (rel.startswith("results/") and re.search(r"\.(tsv|csv)(\.gz)?$", rel)):
                reason = "EXCLUDED_LARGE_NONCORE_ARTIFACT_RETAINED_ON_SERVER"
            if rel in self_metadata:
                checksum, basis = "NA", "SELF_REFERENTIAL_MANIFEST_METADATA_NOT_HASHED"
                if rel == manifest_path.relative_to(root).as_posix():
                    size = "NA"
            elif not final and (rel == "EXECUTION_TRACE.tsv" or rel in mutable_logs):
                checksum, basis = "NA", "PENDING_CURRENT_LOGGER_FINISH_REQUIRES_REFRESH"
            elif reason:
                checksum = known_hashes.get(rel, "NA")
                basis = "RECORDED_CURRENT_RUN_HASH_NOT_REHASHED" if checksum != "NA" else "NOT_HASHED_PACKAGING_EXCLUSION_EXPLICIT"
            else:
                checksum, basis = sha(p), "REHASHED_REPORT_MANIFEST"
            records.append({"relative_path": rel, "bytes": size, "sha256": checksum, "sha256_basis": basis,
                            "review_package_policy": reason or "INCLUDE", "retained_on_server": True})
    # Include the manifest itself even before it exists. Its self hash/size is
    # deliberately not fabricated; the final ZIP obtains an external checksum.
    self_rel = manifest_path.relative_to(root).as_posix()
    if not any(r["relative_path"] == self_rel for r in records):
        records.append({"relative_path": self_rel, "bytes": "NA", "sha256": "NA", "sha256_basis": "SELF_MANIFEST", "review_package_policy": "INCLUDE", "retained_on_server": True})
    with manifest_path.open("w", newline="", encoding="utf-8") as h:
        writer = csv.DictWriter(h, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(sorted(records, key=lambda r: r["relative_path"]))
    return {"files_listed": len(records), "bytes_listed_excluding_self_unknown": sum(r["bytes"] for r in records if isinstance(r["bytes"], int)),
            "review_included_files": sum(r["review_package_policy"] == "INCLUDE" for r in records),
            "excluded_files_individually_listed": sum(r["review_package_policy"] != "INCLUDE" for r in records)}


def main():
    args = args_parse()
    root = args.root.resolve()
    if root.name != PROJECT:
        raise ValueError("Only the authorised new project can be reported")
    d = Inputs(root)
    out = d.path(args.out, exists=False)
    if not out.is_relative_to(root / "reports"):
        raise ValueError("Reports must remain in this new project's reports directory")
    out.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    started, tick = now(), time.monotonic()
    status_path, gate_path = out / "FINAL_REPORTS_STATUS.json", root / "manifest/FINAL_REPORTS_GATE.json"
    trace, finished = canonical_trace(d)
    known_hashes = {}
    for row in trace.to_dict("records"):
        for field in ("input_hashes", "output_hashes"):
            try:
                known_hashes.update({k: v for k, v in json.loads(row[field]).items() if re.fullmatch(r"[0-9a-f]{64}", str(v))})
            except (ValueError, TypeError):
                pass
    if args.refresh_manifest_only:
        state = d.js(str(status_path.relative_to(root)))
        if state.get("status") not in ("REPORTS_GENERATED_PENDING_TRACE_REFRESH", "REPORTS_COMPLETE"):
            raise ValueError("No completed report-generation attempt is available to refresh")
        for name in REPORTS:
            if name != "FILES_MANIFEST.tsv" and not (out / name).is_file():
                raise FileNotFoundError(out / name)
        for rel, expected in require(state, "target_evidence_export_sha256", "completed report-generation state").items():
            if sha(d.path(rel)) != expected:
                raise ValueError("A required target-evidence export changed before final report closure: " + rel)
        chain_export = require(state, "native_chain_diagnostics_export", "completed report-generation state")
        if chain_export.get("donor_stage_rows") != 16 or sha(d.path(chain_export["path"])) != chain_export["sha256"]:
            raise ValueError("The complete native sixteen-chain diagnostic export changed before final report closure")
        for rel, obj in state["source_hashes"].items():
            if rel != "EXECUTION_TRACE.tsv" and sha(d.path(rel)) != obj["sha256"]:
                raise ValueError(f"Source changed after report generation: {rel}; audit rather than refresh away the change")
        # Report's own enclosing logger must now have a FINISH record.
        starts = set(state["report_logger_run_ids_at_generation"])
        if not starts.issubset(set(finished.run_id)):
            raise RuntimeError("Previously active logger(s) still lack FINISH evidence; do not claim final inventory")
        shutil.copyfile(d.path("EXECUTION_TRACE.tsv"), out / "EXECUTION_TRACE.tsv")
        state.update(status="REPORTS_COMPLETE", trace_snapshot_refreshed_at=now(), canonical_trace_sha256=sha(d.path("EXECUTION_TRACE.tsv")))
        write_json(status_path, state)
        write_json(gate_path, {"status": "PASS", "reports_complete": True, "no_models_run_by_reporter": True,
                              "candidate_freeze_sha256": state["candidate_freeze_sha256"], "created_at": now(),
                              "manifest_self_hash_policy": "self-referential metadata explicit NA; final archive requires external SHA256"})
        counts = manifest(root, out, trace, known_hashes, True)
        print(json.dumps({"status": "REPORTS_COMPLETE", **counts}, indent=2)); return
    if status_path.exists():
        raise FileExistsError("Report attempt exists; preserve it. Use --refresh-manifest-only only for its final logger/manifest closure")
    d.register("scripts/write_final_reports.py")
    protocols = [rel for rel in ("PRE_RUN_PROTOCOL.md", "specification/PRE_RUN_PROTOCOL.md") if (root / rel).is_file()]
    if not protocols:
        raise FileNotFoundError("A real pre-run protocol must exist before final reporting")
    for rel in protocols:
        d.register(rel)
    for name in ("PHASE0_GATE", "PHASE1_GATE", "SOFTWARE_GATE", "GWAS_MAGMA_GATE", "BULK_GATE", "BAYESSPACE_GATE", "SPATIALDM_GATE", "INTEGRATION_CONTINUOUS_GATE"):
        if d.js(f"manifest/{name}.json").get("status") != "PASS":
            raise RuntimeError(f"Required gate not PASS: {name}")
    context = d.js(args.context)
    required_context = ("matlab_detected", "matlab_used", "implementation", "hardware_source", "software_gate",
                        "clinical_metadata_path", "author_processed_cohort_available", "author_cohort_evidence", "published_spots", "published_spots_cohort_description",
                        "no_forced_count", "experiments_performed", "experimental_data_generated", "wetlab_only_future_siRNA_qPCR", "figure_visual_QA_path")
    for key in required_context:
        require(context, key, "FINAL_REPORT_CONTEXT")
    if not isinstance(context["published_spots"], int) or context["published_spots"] < 1:
        raise ValueError("published_spots must be an authenticated positive integer, not a guessed placeholder")
    d.register(context["hardware_source"]); d.register(context["author_cohort_evidence"])
    d.register("local_tools/SpatialDM_0.3.1_INPUT_LAYER_DECISION.md")
    provenance = d.tab("manifest/DATA_PROVENANCE.tsv", dtype=str, keep_default_na=False)
    for record in provenance.to_dict("records"):
        if not re.fullmatch(r"[0-9a-f]{64}", record.get("sha256", "")):
            raise ValueError("Immutable-input provenance has an invalid SHA256")
        known_hashes[record["path"]] = record["sha256"]
    hgnc_verification = d.js("data/reference/HGNC/SOURCE_VERIFICATION.json")
    if (hgnc_verification.get("status") != "OFFICIAL_REFERENCE_AUTHENTICATED_BY_REFETCH" or
            hgnc_verification.get("same_SHA256") is not True or hgnc_verification.get("original_reference_unchanged") is not True):
        raise ValueError("Current HGNC original-reference authentication is incomplete or changed the reference")
    hgnc_original = require(hgnc_verification, "original_reference", "HGNC source verification")
    if sha(d.register(hgnc_original["path"])) != hgnc_original["sha256"]:
        raise ValueError("The authenticated original HGNC reference SHA differs")
    hgnc_approved_rows = require(require(hgnc_verification, "original_reference_schema", "HGNC verification"),
                                 "approved_rows", "original HGNC reference schema")
    software = d.js(context["software_gate"])
    if software.get("status") != "PASS" or context["matlab_used"] is not False or context["experiments_performed"] is not False or context["experimental_data_generated"] is not False or context["wetlab_only_future_siRNA_qPCR"] is not True:
        raise ValueError("This report template requires verified no-MATLAB execution and no generated wet-lab data; changed execution requires a documented template revision")
    if context["no_forced_count"] is not True:
        raise ValueError("An unaudited forced spatial cohort count cannot be reported as valid")
    frozen = d.js(args.freeze)
    if frozen.get("status") != "FROZEN" or frozen.get("phase6_complete") is not True:
        raise RuntimeError("Report generation requires the completed new candidate decision, including an explicit no-target freeze if applicable")
    for field in ("evidence_file", "evidence_sha256", "ranking_rule_sha256", "main_target", "backup_targets"):
        if field not in frozen:
            raise ValueError(f"Fresh candidate freeze lacks {field}")
    if sha(d.register(args.ranking_rule)) != frozen["ranking_rule_sha256"] or sha(d.register(frozen["evidence_file"])) != frozen["evidence_sha256"]:
        raise ValueError("Ranking rule / continuous evidence differ from the frozen candidate decision")
    targets = ([frozen["main_target"]] if frozen["main_target"] is not None else []) + frozen["backup_targets"]
    if len(frozen["backup_targets"]) > 2 or len({t["identity_key"] for t in targets}) != len(targets):
        raise ValueError("Candidate identity uniqueness / at-most-two backups violated")
    if frozen["main_target"] is None and (targets or not frozen.get("no_target_reason")):
        raise ValueError("No target requires an explicit reason and no silent backup promotion")
    source_freeze = d.js("results/integrated/continuous/INTEGRATION_SOURCE_FREEZE.json")
    if source_freeze.get("status") != "FROZEN" or source_freeze.get("ranking_performed") is not False:
        raise ValueError("Continuous source catalogue was not separately frozen before ranking")
    source_freeze_time = datetime.fromisoformat(require(source_freeze, "frozen_at", "continuous source freeze").replace("Z", "+00:00"))
    ranking_time = datetime.fromisoformat(require(frozen, "ranking_rule_frozen_at", "candidate freeze").replace("Z", "+00:00"))
    candidate_time = datetime.fromisoformat(require(frozen, "frozen_at", "candidate freeze").replace("Z", "+00:00"))
    if not source_freeze_time <= ranking_time <= candidate_time:
        raise ValueError("Recorded source-freeze / ranking-rule / candidate timestamps are not in the required order")
    if require(frozen, "continuous_source_freeze_sha256", "candidate freeze") != sha(d.register("results/integrated/continuous/INTEGRATION_SOURCE_FREEZE.json")):
        raise ValueError("Candidate freeze does not bind the newly frozen continuous source catalogue")
    continuous_gate = d.js("manifest/INTEGRATION_CONTINUOUS_GATE.json")
    if continuous_gate.get("source_freeze_sha256") != frozen["continuous_source_freeze_sha256"]:
        raise ValueError("Continuous integration gate and candidate decision bind different source freezes")
    for rel, item in source_freeze["source_files"].items():
        if sha(d.register(rel)) != item["sha256"]:
            raise RuntimeError(f"Frozen continuous evidence source changed: {rel}")
    render = d.js(args.figure_dir + "/FIGURE_RENDER_STATUS.json")
    if render.get("status") not in ("RENDERED_PENDING_VISUAL_QA", "RENDERED_VISUAL_QA_PASS") or render.get("figures") != ["Fig1", "Fig2", "Fig3", "Fig4", "Fig5"]:
        raise ValueError("Five completed fresh figure renders are required")
    if render["freeze_sha256"] != sha(d.register(args.freeze)):
        raise ValueError("Figures were not generated from the current fresh target freeze")
    qa = d.js(context["figure_visual_QA_path"])
    if qa.get("status") != "PASS":
        raise ValueError("Final reports require actual completed visual QA, not rendering alone")
    for field in ("reviewed_at", "reviewer", "method", "checked_png_sha256"):
        require(qa, field, "visual QA evidence")
    for fig in ("Fig1", "Fig2", "Fig3", "Fig4", "Fig5"):
        for ext in ("pdf", "png"):
            p = d.register(f"{args.figure_dir}/{fig}.{ext}")
            if render["output_sha256"].get(p.relative_to(root).as_posix()) != sha(p):
                raise ValueError("A required figure differs from its fresh render manifest")
            if ext == "png" and qa["checked_png_sha256"].get(p.relative_to(root).as_posix()) != sha(p):
                raise ValueError("Visual QA did not inspect this exact final PNG")
    figure_notes = d.js(args.figure_dir + "/PANEL_METHODS_AND_CAUTIONS.json")
    panels = d.tab(args.figure_dir + "/source_data/PANEL_PROVENANCE.tsv")
    ver = versions(d, render)
    gqc, gwas = d.js("results/gwas/qc/summary.json"), d.js("results/gwas/MAGMA_SUMMARY.json")
    effect_audit = d.js("results/gwas/qc/RAW_EFFECT_FIELD_AUDIT.json")
    if effect_audit.get("status") != "DESCRIPTIVE_AUDIT_COMPLETE" or effect_audit["data_rows"] != gqc["counts"]["raw_rows"] or effect_audit["input_sha256"] != gqc["reference"]["inputs"]["raw_gwas"]["sha256"]:
        raise ValueError("Full raw effect-field and core GWAS audit provenance/counts disagree")
    if gwas.get("is_new_fit") is not True:
        raise ValueError("Missing new MAGMA fit evidence")
    magma = d.tab("results/gwas/magma_all_genes.tsv", dtype={"gene_id": str})
    if len(magma) != gwas["genes_tested"] or int((magma.BH_FDR < .05).sum()) != gwas["genes_BH_FDR_lt_0.05"]:
        raise ValueError("MAGMA summary and actual complete gene table disagree")
    raw_bulk = d.js("data/metadata/bulk/RAW_COUNT_VALIDATION.json")
    bulk_status = d.js("results/bulk/BULK_STATUS.json")
    recovery = d.js("data/metadata/bulk/OFFICIAL60_RECOVERY_STATUS.json")
    bprotocol = d.js("data/metadata/bulk/QC_PROTOCOL_FROZEN.json")
    samples = d.tab("data/metadata/bulk/SAMPLE_MANIFEST.tsv", dtype=str, keep_default_na=False)
    counts = d.tab("data/metadata/bulk/GROUP_COUNTS_AND_SENSITIVITY.tsv")
    bsummary = d.tab("results/bulk/main/BULK_CONTRAST_SUMMARY.tsv")
    if set(bsummary.contrast) != set(CONTRASTS):
        raise ValueError("All four actual bulk contrasts are required")
    main_samples = samples.loc[samples.main_analysis_pass.map(tf)]
    if len(main_samples) != raw_bulk["n_main_libraries"]:
        raise ValueError("Sample identities and reported retained-library count disagree")
    audits = {t: d.js(f"results/bulk/main/{t}_MODEL_AUDIT.json") for t in ("Duo", "Jej")}
    for audit in audits.values():
        p = d.path(audit["fitted_object"])
        known_hashes[p.relative_to(root).as_posix()] = audit["fitted_object_sha256"]
    if sum(a["n_samples"] for a in audits.values()) != len(main_samples):
        raise ValueError("Actual new DESeq2 model cohorts disagree with frozen sample manifest")
    bnative = native_trace(d, "results/bulk/main/BULK_MODEL_EXECUTION.tsv")
    for tissue in ("Duo", "Jej"):
        if not ((bnative.module == "DESEQ2_NEW_FIT") & (bnative.tissue == tissue) & (bnative.exit_code == "0")).any():
            raise ValueError(f"No actual successful new DESeq2 API trace for {tissue}")
    sensitivity = d.js("results/bulk/sensitivity/BULK_SENSITIVITY_STATUS.json")
    if sensitivity.get("status") != "ANALYSIS_COMPLETE":
        raise ValueError("Fixed bulk QC sensitivity has not completed")
    sensitivity_rows, sensitivity_traces = [], []
    for label, item in require(sensitivity, "analyses", "bulk sensitivity status").items():
        sensitivity_rows.append({"label": label, "native_status": item["status"], "libraries": item["n_libraries"], "new_fit": item["is_new_fit"]})
        if item["status"] == "ANALYSIS_COMPLETE":
            current = d.js(f"results/bulk/sensitivity/{label}/BULK_STATUS.json")
            if current.get("models_all_new_fit") is not True:
                raise ValueError("Changed QC cohort lacks verified new DESeq2 fits")
            snative = native_trace(d, f"results/bulk/sensitivity/{label}/BULK_MODEL_EXECUTION.tsv")
            for tissue in ("Duo", "Jej"):
                if not (snative.module.eq("DESEQ2_NEW_FIT") & snative.tissue.eq(tissue) & snative.exit_code.eq("0")).any():
                    raise ValueError("Actual changed-cohort sensitivity DESeq2 trace absent")
            sensitivity_traces.append(snative)
        elif item["status"] != "IDENTICAL_TO_NEW_MAIN_INPUT_NO_REDUNDANT_REFIT":
            raise ValueError("QC sensitivity analysis is neither completed nor verified identical to current main")
    bulk_comparison, bulk_sensitivity_en, bulk_sensitivity_cn = bulk_qc_comparison(d, list(sensitivity["analyses"]), bsummary)
    pprotocol = d.js("data/metadata/SPATIAL_PREPROCESS_PROTOCOL.json")
    pre = d.js("data/metadata/SPATIAL_PREPROCESS_ALL_DONORS_STATUS.json")
    flow = d.tab("data/metadata/SPATIAL_SPOT_FLOW.tsv")
    clinical, clinical_id = clinical_table(d, context)
    bayes_protocol = d.js("results/spatial/BayesSpace/FROZEN_MODEL_PROTOCOL.json")
    sdm_protocol = d.js("results/spatial/SpatialDM/FROZEN_MODEL_PROTOCOL.json")
    bayes_rows, lr_rows, native_frames, native_chain_rows = [], [], [bnative] + sensitivity_traces, []
    pre_by = {p["donor"]: p for p in pre["donors"]}
    for donor in DONORS:
        b = d.js(f"results/spatial/BayesSpace/{donor}/BAYESSPACE_STATUS.json")
        s = d.js(f"results/spatial/SpatialDM/{donor}/SPATIALDM_STATUS.json")
        if b.get("status") != "BAYESSPACE_COMPLETE" or b.get("is_new_fit") is not True or b.get("actual_chain_files_verified") is not True:
            raise ValueError("Native new BayesSpace fit/chain completion is not verified")
        if s.get("status") != "SPATIALDM_COMPLETE" or s.get("is_new_fit") is not True or s.get("geometry_sensitivity_completed") is not True:
            raise ValueError("Native new SpatialDM fit/sensitivity completion is not verified")
        for stage in ("formal", "stability"):
            chain = d.path(f"results/spatial/BayesSpace/{donor}/{stage}_50000_chain.h5")
            certificate = d.js(f"results/spatial/BayesSpace/{donor}/{stage}_CHAIN_ARRAY_AUDIT.json")
            initial_index = certificate.get("initial_plogLik_R_index")
            if not isinstance(initial_index, list):
                initial_index = [initial_index]  # R auto_unbox retains the single expected index as a scalar.
            nonfinite_indices = certificate.get("full_plogLik_nonfinite_R_indices")
            if not isinstance(nonfinite_indices, list):
                nonfinite_indices = [nonfinite_indices]
            dimensions = certificate.get("dimensions_R_rhdf5", {})
            expected_dimensions = {"z": [50000, b["n_spots"]], "weights": [50000, b["n_spots"]],
                                   "mu": [50000, b["q"] * 15], "lambda": [50000, 225], "plogLik": [50000, 1]}
            if not (b.get("version") == "1.10.1" and b.get("d") == 15 and b.get("burn_in_" + stage) == 1000
                    and certificate.get("status") == "CHAIN_ARRAY_AUDIT_PASS"
                    and certificate.get("version") == "1.10.1"
                    and certificate.get("actual_saved_iterations_R_dimension_order") == 50000
                    and certificate.get("saved_chain_thinning") == 1
                    and all(dimensions.get(name) == shape for name, shape in expected_dimensions.items())
                    and certificate.get("initial_undefined_plogLik_expected_from_official_source") is True
                    and initial_index == [1] and nonfinite_indices == [1]
                    and certificate.get("plogLik_R2_to_nrep_all_finite") is True
                    and certificate.get("posterior_plogLik_all_finite") is True
                    and certificate.get("posterior_diagnostic_states") == 49000
                    and certificate.get("posterior_diagnostic_R_indices") == [1001, 50000]
                    and certificate.get("official_label_aggregation_states") == 49000
                    and certificate.get("official_label_aggregation_R_indices") == [1001, 50000]
                    and certificate.get("native_LL_array_modified") is False
                    and certificate.get("native_labels_recalculated") is False
                    and chain.stat().st_size == certificate.get("file_bytes")):
                raise ValueError(f"Corrected native BayesSpace chain certificate/file size differs: {donor} {stage}")
            p = d.path(f"results/spatial/BayesSpace/{donor}/{stage}_50000.rds")
            known_hashes[p.relative_to(root).as_posix()] = b[stage + "_rds_sha256"]
            diagnostic_rel = f"results/spatial/BayesSpace/{donor}/{stage}_chain_diagnostics.tsv"
            diagnostic = d.tab(diagnostic_rel, dtype=str, keep_default_na=False, na_filter=False)
            required_diag = {"donor", "stage", "seed", "nrep", "burn.in", "effective_size_pseudologlik",
                             "first_half_mean", "second_half_mean", "posterior_mean", "posterior_sd"}
            if len(diagnostic) != 1 or not required_diag.issubset(diagnostic.columns):
                raise ValueError("Exactly one complete native chain-diagnostic row is required: " + diagnostic_rel)
            diagnostic_row = diagnostic.iloc[0].to_dict()
            if (diagnostic_row["donor"] != donor or diagnostic_row["stage"] != stage or
                    diagnostic_row["nrep"] != str(b["nrep_" + stage]) or
                    diagnostic_row["burn.in"] != str(b["burn_in_" + stage]) or
                    diagnostic_row["seed"] != str(b["seed_" + stage])):
                raise ValueError("Native diagnostic identity/seed/iteration fields differ from the actual fit: " + diagnostic_rel)
            diagnostic_row.update(source_path=diagnostic_rel, source_sha256=d.files[diagnostic_rel]["sha256"])
            native_chain_rows.append(diagnostic_row)
        for model in ("global_model.pkl.gz", "sensitivity_radius1.5_model.pkl.gz"):
            d.path(f"results/spatial/SpatialDM/{donor}/{model}")
        if s["n_local_pairs"]:
            d.path(f"results/spatial/SpatialDM/{donor}/local_model.pkl.gz")
        if b["n_spots"] != pre_by[donor]["n_retained"] or b["nrep_formal"] != 50000 or b["nrep_stability"] != 50000:
            raise ValueError("Actual BayesSpace new input/chain metadata is inconsistent")
        bt = native_trace(d, f"results/spatial/BayesSpace/{donor}/MODEL_EXECUTION_TRACE.tsv")
        for stage in ("formal", "stability"):
            if not (bt.module.eq(f"BayesSpace:{donor}:{stage}") & bt.exit_code.eq("0")).any():
                raise ValueError("A real successful BayesSpace chain trace is absent")
        st = native_trace(d, f"results/spatial/SpatialDM/{donor}/MODEL_EXECUTION_TRACE.tsv")
        if not (st.module.eq(f"SpatialDM:{donor}:spatialdm_global") & st.exit_code.eq("0")).any():
            raise ValueError("A real successful SpatialDM global API trace is absent")
        global_array = d.js(f"results/spatial/SpatialDM/{donor}/primary_global_MODEL_ARRAY_AUDIT.json")
        if global_array["actual_null_statistics_per_pair"] != 1000 or global_array["saved_empirical_P_reconstructed"] is not True:
            raise ValueError("Actual saved SpatialDM global null-array audit incomplete")
        if s["n_local_pairs"]:
            local_array = d.js(f"results/spatial/SpatialDM/{donor}/local_MODEL_ARRAY_AUDIT.json")
            if local_array["actual_null_statistics_per_tested_pair_spot"] != 1000:
                raise ValueError("Actual local null-array audit incomplete")
        native_frames.extend([bt, st])
        bayes_rows.append({"donor": donor, "spots": b["n_spots"], "q": b["q"], "ARI": b["independent_seed_ARI"],
                           "warning": b["stability_warning"], "formal": b["nrep_formal"], "stability": b["nrep_stability"],
                           "burn_formal": b["burn_in_formal"], "burn_stability": b["burn_in_stability"]})
        lr_rows.append({"donor": donor, "tested": s["n_global_pairs_tested"], "supported": s["n_global_permutation_BH_significant"],
                        "local_pairs": s["n_local_pairs"], "local_supported_cells": s["n_local_pair_spot_cells_significant"],
                        "global_permutations": s["n_perm_global"], "local_permutations": s.get("n_perm_local"),
                        "local_status": require(s, "local_execution_status", donor + " SpatialDM status"),
                        "geometry": s["geometry_sensitivity_completed"],
                        "strict_QC_global": s.get("strict_z3_5_QC_mask_global_LR_sensitivity", "NOT_REPORTED"),
                        "strict_QC_local": s.get("strict_z3_5_QC_mask_local_LR_sensitivity", "NOT_REPORTED")})
    native_chain_diagnostics = pd.DataFrame(native_chain_rows)
    if len(native_chain_diagnostics) != 16 or set(zip(native_chain_diagnostics.donor, native_chain_diagnostics.stage)) != {(donor, stage) for donor in DONORS for stage in ("formal", "stability")}:
        raise ValueError("All sixteen actual donor-by-stage diagnostics must be retained")
    magma_commands = finished.loc[finished.module.isin(("GWAS_MAGMA_NEW_ANNOTATION", "GWAS_MAGMA_NEW_GENE_ASSOCIATION"))].copy()
    magma_commands["source_trace_path"] = "EXECUTION_TRACE.tsv"
    native_all = pd.concat(native_frames + [magma_commands], ignore_index=True, sort=False)
    if any(datetime.fromisoformat(t.replace("Z", "+00:00")) > source_freeze_time for t in native_all.finish_iso8601):
        raise ValueError("A model API call finished after the declared pre-ranking source freeze")
    # Validate MAGMA actual command exits rather than infer success from a table.
    for module in ("GWAS_MAGMA_NEW_ANNOTATION", "GWAS_MAGMA_NEW_GENE_ASSOCIATION"):
        if not (finished.module.eq(module) & finished.exit_code.eq("0")).any():
            raise ValueError("Actual successful logged MAGMA command absent: " + module)
    frozen_keys = {t["identity_key"] for t in targets}
    target_parts = []
    for chunk in pd.read_csv(d.register(frozen["evidence_file"]), sep="\t", chunksize=100000, dtype=str, keep_default_na=False, na_filter=False):
        target_parts.append(chunk.loc[chunk.identity_key.isin(frozen_keys)])
    target_evidence = pd.concat(target_parts, ignore_index=True) if target_parts else pd.DataFrame()
    if targets and set(target_evidence.identity_key) != frozen_keys:
        raise ValueError("Frozen targets are absent from actual continuous evidence")
    target_wide, target_wide_provenance = frozen_target_wide(d, frozen, continuous_gate, frozen_keys)
    target_bulk_qc = target_bulk_qc_view(target_wide, list(sensitivity["analyses"]))
    target_self_lr, target_lr_scope, target_lr_scope_en, target_lr_scope_cn = frozen_target_lr_scope(d, targets, target_wide, target_wide_provenance)
    old_text, old_files = old_comparison(d, args, True)
    total_spots = sum(p["n_retained"] for p in pre_by.values())
    published_bulk_n = sum(recovery["published_counts"].values())
    group_fields = [(c, c) for c in counts.columns]
    clinical_fields = [(c, c) for c in clinical.columns]
    target_names = ", ".join(f"{t['gene_symbol']} ({t['identity_key']})" for t in targets) or "none; " + frozen["no_target_reason"]
    feasibility_en, feasibility_cn = frozen_feasibility_text(d, frozen, targets)
    failures = finished.loc[finished.exit_code.ne("0")].to_dict("records")
    unfinished = sorted(set(trace.loc[trace.event.eq("START"), "run_id"]) - set(finished.run_id))
    report_logger_runs = trace.loc[trace.run_id.isin(unfinished) & trace.source_command.str.contains("write_final_reports.py", regex=False), "run_id"].tolist()
    checkpoint_mentions = []
    for rel in sorted({r for r in finished.get("stdout_log", []) if r and r != "NA"}):
        p = d.path(rel, exists=False)
        if not p.is_file():
            raise FileNotFoundError("Completed-command stdout log is missing: " + rel)
        if p.stat().st_size == 0:
            continue  # A quiet command legitimately has an empty stdout file.
        log = d.text(rel)
        hits = [line for line in log.splitlines() if "REUSE_VERIFIED_NEW_" in line]
        checkpoint_mentions.extend({"source": rel, "message": line} for line in hits)
    # Terminology is fixed before composing prose; no target biology is inferred.
    ledger = "| Term | Meaning used here |\n| --- | --- |\n| Library | RNA-seq count library; not automatically a unique person |\n| Donor | One of the eight independent spatial reference donors |\n| Raw UMI | Unaltered nonnegative integer spatial counts |\n| Corrected expression | Author-background-subtracted, clipped noninteger expression |\n| Marker FDR | Native within-donor scran result, not a disease-association FDR |\n| LR support | Existing permutation BH support for a tested spatial pair |\n| Untested | No corresponding statistical test; never a zero effect |\n"
    source_note = "All numerical statements below are generated from the linked current-run tables/status files. Commands and model timing are preserved in EXECUTION_TRACE.tsv and NATIVE_MODEL_EXECUTION_TRACE.tsv; nested/concurrent durations must not be summed as wall time.\n"
    chain_audit_en = ("The formal and stability chain certificates for all eight donors were checked against their native HDF files. "
                      "Official BayesSpace 1.10.1 leaves the initial pseudo-log-likelihood at R index 1 undefined (NA); that value was retained. "
                      "R indices 2–50,000 were finite, and the 49,000 posterior diagnostic states at R indices 1,001–50,000 were finite and matched the official label-aggregation window. "
                      "No native likelihood array or labels were rewritten. The documented initial-NA audit failure arose from an overstrict validation guard after the native API returned, not a native API failure. "
                      "Passing this array audit does not establish convergence or imply high independent-seed ARI.\n")
    chain_audit_cn = ("已逐一核验八供者formal/stability链证书与原生HDF文件。官方BayesSpace 1.10.1将R索引1的初始伪对数似然留为未定义NA，该值原样保留；"
                      "R索引2–50,000均有限，R索引1,001–50,000的49,000个后验诊断状态均有限，且与官方标签聚合窗口一致。"
                      "没有改写原生似然数组或标签。已记录的初始NA审计失败来自原生API返回后过严的验证guard，不是原生API调用失败；"
                      "数组审计通过不等于证明收敛，也不意味着独立seed的ARI较高。\n")
    diagnostic_table = md_table(native_chain_diagnostics.to_dict("records"), [(c, c) for c in native_chain_diagnostics.columns])
    chain_diagnostics_en = ("## Native chain diagnostics: all sixteen formal/stability chains\n\n"
        "Every saved effective sample size (ESS), half-chain mean, posterior mean and standard deviation below is read as its original string from the native diagnostic table, without recalculation or a new cutoff. "
        "Low ESS and moderate independent-seed ARI remain visible. The native warning value `none` means that the implementation's ARI<0.6 warning was not triggered; "
        "it is not a certificate of convergence, mixing or biological reproducibility. Complete finite chains do not rule out autocorrelation. "
        "All sixteen rows and source SHA256 are also retained in `BAYESSPACE_NATIVE_CHAIN_DIAGNOSTICS.tsv`.\n\n" + diagnostic_table)
    chain_diagnostics_cn = ("## 全部十六条正式/稳定性链的原生诊断\n\n"
        "下表逐项从原生诊断表读取ESS、前后半链均值、后验均值和标准差的原始字符串，没有重新计算或新增阈值。"
        "低ESS与中等独立seed ARI均不删除。原生warning为`none`仅表示代码的ARI<0.6提示未触发，不表示已证明收敛、混合充分或生物学可复现；"
        "完整且有限的链仍可能有自相关。全部十六行及来源SHA256另存`BAYESSPACE_NATIVE_CHAIN_DIAGNOSTICS.tsv`。\n\n" + diagnostic_table)
    methods = f"""# Computational methods

This independent computational reanalysis combined fixed IBS GWAS GCST90016564, GSE166869 small-intestinal patient RNA-seq and the eight-donor Harnik adult proximal small-intestinal spatial reference. Data streams were analysed separately before continuous evidence integration; the spatial reference was not an IBS patient cohort. {source_note}

## Inputs, identities and reproducibility

Official raw/reference inputs, new-source hashes and source-code snapshots are archived in `manifest/DATA_PROVENANCE.tsv`, the module manifests and `results/integrated/continuous/INTEGRATION_SOURCE_FREEZE.json`. A copied immutable official download/reference or software archive was not treated as reuse of an earlier statistical model. Approved HGNC IDs defined the canonical identity. Ambiguous identifiers remained unresolved; all original feature rows were linked. Bulk representatives used maximum baseMean, then lexical original feature ID, rather than minimum P. The continuous source catalogue was frozen before the separate recorded ranking rule `{args.ranking_rule}` was applied. The resulting decision file was `{args.freeze}`; this report did not nominate targets.

The current approved-human HGNC complete-set reference contained {hgnc_approved_rows} approved rows, read from the authenticated source schema. An independent official refetch reproduced the same SHA256 (`{hgnc_original['sha256']}`); authentication status was `{hgnc_verification['status']}` and the original reference remained unchanged. This was provenance authentication, not a replacement identity mapping. The original first-download time was unavailable and was not reconstructed from filesystem mtime. The actual refetch interval and source details are retained in `data/reference/HGNC/SOURCE_VERIFICATION.json` and the source-verification appendix below.

## GWAS and gene-level association

The GRCh37 GWAS was audited across {gqc['counts']['raw_rows']} raw rows. Primary autosomal variants required finite P in (0,1], valid coordinates and positive integer case/control counts. Every repeated raw SNP identifier was removed entirely, and only unique exact SNP-ID/chr/BP matches to the EUR 1000 Genomes Phase3 reference were admitted. Per-SNP total N was N_CASE+N_CONTROL, with observed valid-raw range {gqc['numeric_ranges']['N_TOTAL']}; no article-wide or effective-N constant replaced the supplied values. These sample-size values are reported without an unsupported explanation for any difference from the source publication. New official MAGMA v1.10 annotation used the verified historical NCBI37.3 protein-coding gene coordinates and {gwas['annotation_window_kb']} kb flanks. This gene-body window was an analyst choice made before fitting. The SNP-wise mean model used `use=SNP,P ncol=N duplicate=error`; `synonyms=0` preserved the frozen exact-ID/coordinate universe. All {gwas['genes_tested']} tested gene P values entered BH correction before expression or HGNC filters. Gene association was not interpreted as causal-gene identification or an independent-locus count. Manhattan data retained all P<5e-8 variants plus the declared deterministic CRC32 background sample; QQ quantiles used the full valid unique autosomal GWAS. Lambda-GC was descriptive and was not applied as a correction.

## Patient RNA-seq

Only verified nonnegative integer raw counts were used. The actual new cohort comprised {len(main_samples)} libraries; the source publication's group totals sum to {published_bulk_n}. Official individual-ID recovery status was `{recovery['status']}`. When the original list was not authenticated, no samples were discarded to satisfy group quotas; the fixed independent thresholds were {bprotocol['main_min_non_rrna_counts']} non-rRNA counts and {bprotocol['main_min_detected_non_rrna_genes']} detected non-rRNA genes. Tissues were fitted separately, without inventing subject pairing or batch IDs. Low-expression filtering was `{bprotocol['low_expression_filter']}`. Actual designs were Duo `{audits['Duo']['design']}` and Jej `{audits['Jej']['design']}`; matrix ranks and covariate decisions are preserved in their MODEL_AUDIT.json files. DESeq2 used HV as reference, Wald tests, fitType={bprotocol['deseq2']['fitType']}, sfType={bprotocol['deseq2']['sfType']}, betaPrior={bprotocol['deseq2']['betaPrior']}, and the recorded default Cook's and independent-filtering behaviour. Full P, adjusted P, effect estimates and native NA reasons were retained. This was an independent DESeq2 reanalysis, not an exact reproduction of the source-paper glmLRT workflow. The fixed QC sensitivity compared the declared {bprotocol['sensitivity_min_non_rrna_counts']} thresholds without significance-driven cohort selection. GSEA used {bprotocol['gsea']['method']} with finite signed Wald statistics, maximum-baseMean duplicate-symbol representatives, {bprotocol['gsea']['collections']}, minSize={bprotocol['gsea']['minSize']}, maxSize={bprotocol['gsea']['maxSize']}, eps={bprotocol['gsea']['eps']} and nPermSimple={bprotocol['gsea']['nPermSimple']}; primary BH pooled both tested collections within each contrast.

## Full fixed-cohort patient QC comparison

The native full-feature summary reports Spearman correlation and sign agreement over finite common LFCs, both cohorts' original corrected-significance counts, their intersection and main-significant entries whose QC adjusted P is unavailable. Reporting reads all fixed cohort-by-contrast rows without recomputing a test or choosing a preferred QC result.

{bulk_sensitivity_en}

## Spatial preprocessing and domain analysis

The author source commit was `{pprotocol['author_commit']}`. MATLAB detected={context['matlab_detected']}; MATLAB used={context['matlab_used']}. The implementation was `{context['implementation']}`, explicitly an independent reimplementation rather than an official MATLAB execution. Background means were calculated from counts_ALL barcodes outside counts_UTT on original gene rows, subtracted and clipped at zero before duplicate-symbol aggregation. Primary thresholds were corrected log10-UMI Z>{pprotocol['Z_THRESH_SUM_UMI']} and mitochondrial-fraction Z<{pprotocol['Z_THRESH_MITO']}; source-defined Hb/boundary rules and every matched author manual exclusion were applied. Unmatched barcode problems were audited instead of retained by default. Raw integer UMI and noninteger corrected expression remained distinct. The primary author-QC reimplementation cohort contained {total_spots} spots. The published comparator of {context['published_spots']} refers to: {context['published_spots_cohort_description']}. These different analysis cohorts have not been shown equivalent; the published number is not asserted to be every spot after the author's importer/QC. The provenance, source-access scope and unresolved cohort-transformation details are recorded in `{context['author_cohort_evidence']}`. No count was forced or cropped to that comparator. The spatial reference was surgical tissue, not a set of unscreened healthy volunteers; its actual archived clinical fields are reproduced below in the data/QC report.

Official BayesSpace {ver['R_packages']['BayesSpace']} was called separately for each donor. SpatialPreprocess settings were `{json.dumps(bayes_protocol['spatialPreprocess'], ensure_ascii=False)}`; qTune settings were `{json.dumps(bayes_protocol['qTune'], ensure_ascii=False)}`. Donor q followed `{bayes_protocol['q_selection']}`, independent of target statistics. The formal and independent-seed stability chains each used 50,000 iterations; their actual burn-ins and seed ARIs are reported per donor. Model/precision/gamma/init settings and all recorded parameters are reproduced in the exact parameter appendix. Clusters were donor-specific domains, not automatically equivalent cell types. Native scran markers were within-donor domain tests; they were not IBS-versus-control tests.

{chain_audit_en}

{chain_diagnostics_en}

## Spatial ligand-receptor coexpression

Official SpatialDM {ver['Python_packages']['spatialdm']} with bundled human CellChatDB tested every complete/detectable eligible pair within each donor. Global input X was `{sdm_protocol['X']}`; the local raw interface was `{sdm_protocol['local_raw_interface']}`. In this project's adaptation, the local nonlog layer used author-background-corrected CP10000 with the full-gene per-spot denominator applied before restricting to required LR genes. This aligned the source author's normalization and is not a claim that SpatialDM officially recommends replacing raw UMI with this adaptation. Unaltered integer UMI was preserved separately. Coordinates were scaled using actual median spot-neighbour spacing; primary secreted/contact geometry and the distance sensitivity were fixed in FROZEN_MODEL_PROTOCOL.json. Primary global P values came from 1,000 actual permutations per tested pair; saved null-array audits reconstructed the reported empirical P. Analytic z-score P/FDR were kept separate. Local tests used the frozen global-based selection and cap and, when performed, 1,000 permutations. A donor with no qualifying global pair had no local test, rather than an asserted completed local permutation run. Local BH covered all selected pair-by-spot cells within that donor. Empirical zero P was not interpreted as probability zero; figure-only floors did not change canonical statistics. A geometry sensitivity changed the secreted radius on the same primary cohort. Strict-QC expression/detection sensitivity did not re-fit BayesSpace or global/local LR models and must not be described as such. Same-gene pairs denote spatial coexpression, not direct molecular binding or demonstrated sender-to-receiver flux. Across-donor summaries used donors as biological replicates; spots were not pooled as independent donors.

## Evidence integration and experimental boundary

All tested/available genes were retained as continuous evidence before ranking. Native FDR-significant, nonsignificant, unavailable and untested states were kept distinct, including all four patient contrasts, donor expression/domain markers, fixed expression-QC sensitivity and LR eligibility/global/local evidence. No mixed-unit score was created by the integration module. The separate new decision froze: {target_names}. Experimental data generated={context['experimental_data_generated']}; any siRNA+qPCR work is future feasibility work and cannot establish visceral pain causality or therapeutic benefit in this computation-only run.

## Exact parameter and version appendix

```json
{json.dumps({'bulk_protocol': bprotocol, 'spatial_preprocessing': pprotocol, 'BayesSpace': bayes_protocol, 'SpatialDM': sdm_protocol, 'versions': ver, 'software_gate': software, 'HGNC_reference_verification': hgnc_verification}, indent=2, ensure_ascii=False)}
```

## Terminology

{ledger}
"""
    results = f"# Computational results\n\n{source_note}\n## Data and cohort audit\n\nThe new patient analysis retained {len(main_samples)} verified count libraries, compared with the published group total {published_bulk_n}; original individual-cohort recovery was `{recovery['status']}`. These library counts are not a claim of {len(main_samples)} independent unique subjects. Unmapped count libraries quarantined: {raw_bulk['n_quarantined_unmapped_libraries']}. Explicit subject IDs present in source metadata: {raw_bulk['n_explicit_subject_ids']}.\n\n" + md_table(counts.to_dict("records"), group_fields)
    results += f"\n## GWAS and patient expression\n\nThe new MAGMA model tested {len(magma)} genes; {gwas['genes_BH_FDR_lt_0.05']} passed BH FDR<0.05 and {gwas['genes_bonferroni_0.05']} passed its recorded Bonferroni threshold. The complete table, including all nonsignificant genes, is `results/gwas/magma_all_genes.tsv`. No independent-locus count is inferred.\n\n" + md_table(bsummary.to_dict("records"), [(c, c) for c in bsummary.columns])
    results += f"\nThe primary spatial cohort retained {total_spots} spots. The published {context['published_spots']} comparator is described as: {context['published_spots_cohort_description']}. This is a different analysis-cohort comparison, not an established equality with all post-importer/QC spots; the new cohort was not cropped to match it. All eight donor-specific domain fits and independent-seed chains completed with the following actual values; warnings are retained.\n\n" + md_table(bayes_rows, [(k, k) for k in bayes_rows[0]])
    results += "\n" + chain_audit_en
    results += "\n" + chain_diagnostics_en
    results += "\n## Spatial coexpression and untested analyses\n\n" + md_table(lr_rows, [(k, k) for k in lr_rows[0]])
    results += f"\nThese pair/spot counts are not gene-level effect sizes or independent donor counts. Analytic and permutation statistics are separate. Strict-QC-mask BayesSpace/global-LR/local-LR model sensitivity was not performed; only descriptive expression/detection sensitivity and the separately recorded geometry model sensitivity were available.\n\n## Fresh target decision\n\nFrozen identities: {target_names}. The decision and its separate ranking rule are `{args.freeze}` and `{args.ranking_rule}`. `FROZEN_TARGET_ALL_CONTINUOUS_EVIDENCE.tsv.gz` is the interface-axis long export, not an exhaustive record of every fixed-bulk-QC column. `FROZEN_TARGET_WIDE_ALL_MODULE_COLUMNS.tsv.gz` preserves all original wide-catalogue columns and string values for the selected HGNC identities, including fixed-QC LFC/P/padj/NA fields; its SHA-bound source and pure row-selection provenance are recorded separately. No backup was promoted, no new score was calculated and no target was ranked by this report. This candidate decision is exploratory and does not establish causality, direct binding, pain mechanism or a treatment effect.\n"
    results += "\n" + feasibility_en
    results += "\n## Fixed patient QC sensitivity\n\n" + md_table(sensitivity_rows, [(k, k) for k in ("label", "native_status", "libraries", "new_fit")]) + "\n" + bulk_sensitivity_en
    timing_fields = [(k, k) for k in ("module", "start_iso8601", "finish_iso8601", "duration_seconds", "exit_code", "source_trace_path")]
    results += "\n## Actual API execution times\n\nThese are actual recorded starts, finishes and monotonic elapsed seconds. Concurrent and nested intervals overlap; their sum is not the end-to-end wall-clock duration. Interrupted attempts with no finish are not assigned an inferred duration.\n\n" + md_table(native_all.to_dict("records"), timing_fields)
    failed_table = md_table(failures, [(k, k) for k in ("module", "run_id", "start_iso8601", "finish_iso8601", "duration_seconds", "exit_code", "stderr_log", "notes")]) if failures else "No nonzero canonical FINISH exits were recorded in this trace snapshot.\n"
    checkpoint_text = md_table(checkpoint_mentions, [("source", "Current-run stdout path"), ("message", "Exact verified-new reuse message")]) if checkpoint_mentions else "No verified-new-project reuse message was found in completed stdout logs.\n"
    limitations = f"""# Limitations, negative evidence, failures and recovery

The source-paper original patient-cohort identity status is `{recovery['status']}`. The new cohort contains {len(main_samples)} libraries, not an assumed exact reproduction of {published_bulk_n} individual samples. Group totals do not identify which libraries to remove. Unknown pairing and technical batches were not fabricated; per-tissue designs and their limitations are preserved in MODEL_AUDIT.json.

The author processed spatial cohort was available={context['author_processed_cohort_available']}; new retained author-QC-reimplementation spots={total_spots}. The published comparator={context['published_spots']} refers to {context['published_spots_cohort_description']}, not an automatically equivalent total post-importer/QC cohort. The difference between these analysis cohorts remains unresolved within the source-access scope documented in `{context['author_cohort_evidence']}`; an inaccessible full Methods section is not evidence that the paper contains no explanation. The preprocessing is `{context['implementation']}` and was not official MATLAB execution. The surgical reference donors are not an IBS cohort and are not labelled healthy volunteers; original clinical metadata are reproduced in DATA_SOURCE_AND_QC_CN.md. Background-corrected values can be noninteger and are never called raw integer UMI. Donor-specific domains/markers do not directly establish disease-state cell types. Low independent-seed stability warnings, if present in the complete donor table, are not suppressed.

{chain_audit_en}

Every native ESS and half-chain/posterior diagnostic is shown unchanged in the all-sixteen-chain table in Computational_Results_EN.md and BAYESSPACE_NATIVE_CHAIN_DIAGNOSTICS.tsv. A native ARI-warning value of `none` only means the source code's ARI<0.6 flag was not triggered; it does not certify convergence, adequate mixing or biological reproducibility. Low ESS and moderate ARI are not hidden by the successful-execution gate, and no new threshold or extra model fit was introduced by reporting.

Strict-QC-mask model sensitivity was not performed for BayesSpace or global/local SpatialDM. The available strict-QC comparison concerns expression/detection and inclusion masks, and the separate LR geometry sensitivity uses the same primary cohort. Local input was the project-specific author-normalized background-corrected nonlog CP10000 layer (`{sdm_protocol['local_raw_interface']}`), not unaltered raw UMI or an asserted SpatialDM default recommendation. Local pairs not selected by the pre-recorded cap were not tested locally. No-global-significance donors have an explicit no-local-run state. Permutation and analytic P are separate, and an empirical zero does not represent probability zero. Same-gene LR coexpression is not direct molecular binding or proven directed signalling. Overlapping LR components and spots do not add independent donors.

Nonsignificant FDR, native NA reasons, identity ambiguity, absent/unavailable catalogue entries and untested analyses remain in the full source tables and target evidence export. A missing adjusted P is never converted to nonsignificance. No siRNA, qPCR, functional assay, pain behaviour or treatment experiment was performed. Future siRNA+qPCR alone would not establish pain causality or treatment efficacy. Whole-gene MAGMA association does not establish a causal target, and symbol/reference versions have their stated limits.

## Complete fixed-cohort patient QC sensitivity

{bulk_sensitivity_en}

Corrected-significance sensitivity and continuous-effect agreement must be retained together in target interpretation. A particular target's primary and fixed-QC P/adjusted-P/NA fields remain available in `FROZEN_TARGET_WIDE_ALL_MODULE_COLUMNS.tsv.gz`; the figure-interface long export alone is not asserted to contain every QC field. Source-column values are selected as strings without statistical recalculation or candidate reordering.

## Actual failed attempts retained

{failed_table}
## Trace and checkpoint scope

At prose-generation time, unmatched START records were {json.dumps(unfinished)}. They were incomplete logger records, not evidence of a completed or lost model. The enclosing report logger may have been active then; final trace/manifest closure appends its later FINISH before the final trace snapshot. Verified checkpoint-reuse messages refer only to this new project's own checked fits; copied official raw/reference archives and software are a separate provenance category. No file modification timestamp was used as a runtime.

{checkpoint_text}
"""
    qc_cn = f"# 数据来源与QC\n\n本报告只读取本次新工程。GWAS实际原始行数 {gqc['counts']['raw_rows']}，可用于MAGMA的精确LD匹配SNP数 {gwas['input_eligible_snps']}；N实际范围 {gqc['numeric_ranges']['N_TOTAL']}，定义为每SNP病例加对照。P=0行数 {gqc['counts']['p_equal_zero']}；完整排除原因见 results/gwas/qc/summary.json。\n\n患者主分析 {len(main_samples)} 个已核验文库；发表组别合计 {published_bulk_n}。原始60份逐ID状态：{recovery['status']}。未映射隔离列 {raw_bulk['n_quarantined_unmapped_libraries']}，不按配额任意删样本，不把文库数写成独立患者数。\n\n" + md_table(counts.to_dict("records"), group_fields)
    qc_cn += f"\nbeta/SE/EAF与等位基因字符串另做全量描述性审计，实际覆盖 {effect_audit['data_rows']} 行，flagged_data_rows={effect_audit['flagged_data_rows']}。其范围是 `{effect_audit['scope']}`，没有更改P/N、SNP纳入、LD匹配或效应校准。全部类型/缺失/范围与分类见 results/gwas/qc/RAW_EFFECT_FIELD_AUDIT.json / .tsv。\n"
    qc_cn += f"\n空间新主分析保留作者QC独立重实现的 {total_spots} spots；原文对照数 {context['published_spots']} 的实际cohort说明为：{context['published_spots_cohort_description']}。不同分析cohort尚未证明等价，不能无条件把该发表数字称为全部QC后spots，也未裁剪新cohort去凑该数。来源访问范围和未解决的cohort转换规则见 `{context['author_cohort_evidence']}`；完整Methods不可访问不等于论文没有解释。全部作者人工排除、背景校正及条码来源见 SPATIAL_SPOT_FLOW.tsv / SPATIAL_BARCODE_PROVENANCE.tsv。MATLAB detected={context['matlab_detected']}，used={context['matlab_used']}；为作者MATLAB流程独立Python重实现。严格QC仅表达/检出率描述性敏感性；没有严格mask下的域或LR模型重跑。\n\n以下为官方Data.zip metadata.csv经本次Phase0导出的八位供者原始临床字段，逐列原样展示；供者列为 `{clinical_id}`。这些是手术获得近端小肠参考，不应写为健康志愿者或IBS患者。\n\n" + md_table(clinical.to_dict("records"), clinical_fields)
    review_cn = f"# 结果复核\n\n本次新模型及候选已冻结后才生成报告；不从旧结果修正新结果。主/备冻结身份：{target_names}。MAGMA {gwas['genes_tested']} 个实际检验基因，BH<0.05 为 {gwas['genes_BH_FDR_lt_0.05']}；四对比bulk统计见下表，包括全部NA原因的完整源表。\n\n" + md_table(bsummary.to_dict("records"), [(c, c) for c in bsummary.columns])
    review_cn += "\n## 完整预设患者QC敏感性\n\n" + bulk_sensitivity_cn
    review_cn += "\n八供者域模型与独立seed稳定性：\n\n" + md_table(bayes_rows, [(k, k) for k in bayes_rows[0]])
    review_cn += "\n" + chain_audit_cn
    review_cn += "\n" + chain_diagnostics_cn
    review_cn += "\nLR真实置换及local未运行状态：\n\n" + md_table(lr_rows, [(k, k) for k in lr_rows[0]])
    if not target_evidence.empty:
        target_summary = target_evidence.loc[target_evidence.evidence_axis.str.match(r"^(MAGMA_(P|BH_FDR)|BULK_.*_(log2FoldChange|padj)|SPATIAL_.*_corrected_detection_fraction|LR_.*_global_pairs_(tested_n|BH_supported_n|untested_n))$")]
        review_cn += "\n冻结靶点的固定证据轴摘要（按轴类型展示，未依P值筛行；全部域/局部细节仍保留源表）：\n\n" + md_table(target_summary.to_dict("records"), [(k, k) for k in ("identity_key", "gene_symbol", "evidence_axis", "display", "measurement_status", "source_path")])
    target_wide_note_en = ("The target-wide export retains every original source column/string and the original source row order, not candidate rank order. "
                          "It is a SHA-bound HGNC row filter, not a new score or a replacement for all-gene sources. "
                          "The readable bulk-QC projection below retains native P/adjusted-P/NA values; NA is not a negative test.\n\n")
    target_wide_note_cn = ("`FROZEN_TARGET_ALL_CONTINUOUS_EVIDENCE.tsv.gz`仅为图接口证据轴的长表，不能据此声称覆盖全部bulk QC列。" 
                          "新增`FROZEN_TARGET_WIDE_ALL_MODULE_COLUMNS.tsv.gz`按冻结HGNC身份筛行，保留原宽表全部列、全部字符串值及源行顺序；" 
                          "原源SHA和筛行出处另存provenance。没有新评分、统计重算或候选重排。下列可读bulk-QC表仍保留原生P/q/NA，不把NA填成阴性。\n\n")
    if not target_bulk_qc.empty:
        target_qc_table = md_table(target_bulk_qc.to_dict("records"), [(c, c) for c in target_bulk_qc.columns])
        results += "\n## Frozen targets: complete patient-QC projection\n\n" + target_wide_note_en + target_qc_table
        review_cn += "\n## 冻结靶点全部患者QC细节\n\n" + target_wide_note_cn + target_qc_table
    else:
        results += "\n## No selected target rows\n\n" + target_wide_note_en + "The fresh freeze nominated no target: " + md(frozen["no_target_reason"]) + ". Both target-wide and readable target-QC outputs have real headers and zero selected rows.\n"
        review_cn += "\n## 无冻结靶点行\n\n" + target_wide_note_cn + "本次freeze没有选定靶点，原因：" + md(frozen["no_target_reason"]) + "。两个导出均保留真实表头、选中行数为0，不制造target记录。\n"
    results += "\n" + target_lr_scope_en
    review_cn += "\n" + target_lr_scope_cn
    limitations += "\n" + target_lr_scope_en + "\n" + target_lr_scope_cn
    review_cn += "\n" + feasibility_cn
    review_cn += "\n必须保留：所有未显著统计、bulk原生NA/独立过滤、未检验LR、无数据库对应、低ESS/中等ARI、严格QC仅表达而非模型的范围限制。图仅采用本次新freeze靶点。以上已展示冻结时的来源支持可行性判断及其限制；后续实验前仍须核实实际模型/试剂获取和目标模型中的表达。未来siRNA+qPCR为验证起点，当前没有实验数据，不能写成疼痛因果或治疗有效。\n"
    legends = "# Figure legends\n\nAll panels derive from the new run and the current candidate freeze. Displayed targets: " + target_names + ". Biological replicate units, statistical families, empirical-P display censoring and missing-state semantics are preserved below.\n\n"
    fixed_legends = {"Fig1": "Data provenance and audited analysis units. Published cohort totals are descriptive comparators and do not identify a sample list. Library and spot counts are distinct from independent persons/donors.",
                     "Fig2": f"GWAS and gene association. The Manhattan background uses CRC32 modulo {gqc['plotting']['manhattan_modulus']} sampling while retaining all P<5e-8 points; the underlying universe contains {gqc['counts']['valid_unique_autosomal_snps']} valid unique autosomal variants. QQ uses full-universe order statistics. MAGMA gene results use all tested-gene BH correction, not a causal-gene or locus count.",
                     "Fig3": "Patient RNA-seq and pathway evidence. The four fixed contrasts are relative to HV within tissue. Log2 fold change, Wald statistics, native adjusted P and pathway enrichment retain their units. Source-paper original-cohort identity recovery and actual new cohort sizes are shown explicitly. Panel g retains every fixed cohort/contrast's native FDR-significant count, primary-set overlap, continuous LFC correlation/sign agreement and primary-hit sensitivity-q NA count. Agreement of continuous effects does not establish a robust corrected-significance set; no cohort or target was reselected. Panel h retains the original pooled-FDR GSEA results without a new display-subset adjustment.",
                     "Fig4": "Eight-donor domain analysis and fresh target expression. Domains are donor-specific BayesSpace clusters, not matched cell types; formal and independent-seed chains each used 50,000 iterations. Expression is from the stated corrected layer, not raw integer UMI. Spatial reference donors are surgical and not IBS patients or healthy volunteers.",
                     "Fig5": "SpatialDM and continuous target evidence. Global permutation BH and analytic supplementary statistics are distinguished. Representative local panels follow the fixed documented selection; unselected pairs are not local negatives. Count summaries are donor-pair records. Empirical zero P is display-censored only, and native statistics are unchanged. Original evidence units and explicit missing/untested/descriptive states prevent a mixed-unit composite score. Same-gene pairs do not establish direct binding or pain causality."}
    for fig in ("Fig1", "Fig2", "Fig3", "Fig4", "Fig5"):
        if fig not in figure_notes:
            raise ValueError("Actual figure panel caution metadata missing: " + fig)
        legends += f"## {fig}\n\n{fixed_legends[fig]}\n\nActual rendering/panel rules:\n\n```json\n{json.dumps(figure_notes[fig], indent=2, ensure_ascii=False)}\n```\n\n"
    # Freeze metadata only after every required value has been obtained; no
    # placeholder number or fabricated successful stage can reach a report.
    write_json(status_path, {"status": "WRITING", "started_at": started})
    texts = {"Computational_Methods_EN.md": methods, "Computational_Results_EN.md": results,
             "RESULTS_REVIEW_CN.md": review_cn, "DATA_SOURCE_AND_QC_CN.md": qc_cn,
             "FIGURE_LEGENDS_EN.md": legends, "LIMITATIONS_AND_FAILURES.md": limitations,
             "PROJECT_REBUILD_DIFF_VS_OLD.md": old_text}
    for filename, text in texts.items():
        (out / filename).write_text(text, encoding="utf-8")
    shutil.copyfile(d.path("EXECUTION_TRACE.tsv"), out / "EXECUTION_TRACE.tsv")
    native_all.to_csv(out / "NATIVE_MODEL_EXECUTION_TRACE.tsv", sep="\t", index=False, na_rep="NA")
    chain_diagnostic_path = out / "BAYESSPACE_NATIVE_CHAIN_DIAGNOSTICS.tsv"
    native_chain_diagnostics.to_csv(chain_diagnostic_path, sep="\t", index=False)
    diagnostic_roundtrip = pd.read_csv(chain_diagnostic_path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)
    if list(diagnostic_roundtrip.columns) != list(native_chain_diagnostics.columns) or diagnostic_roundtrip.to_dict("records") != native_chain_diagnostics.to_dict("records"):
        raise ValueError("Native sixteen-chain diagnostic strings changed during report export")
    target_evidence.to_csv(out / "FROZEN_TARGET_ALL_CONTINUOUS_EVIDENCE.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    target_wide_path = out / "FROZEN_TARGET_WIDE_ALL_MODULE_COLUMNS.tsv.gz"
    target_wide.to_csv(target_wide_path, sep="\t", index=False, compression="gzip")
    target_wide_roundtrip = pd.read_csv(target_wide_path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)
    if list(target_wide_roundtrip.columns) != list(target_wide.columns) or target_wide_roundtrip.to_dict("records") != target_wide.to_dict("records"):
        raise ValueError("Target-wide all-column string values changed during export")
    target_qc_path = out / "FROZEN_TARGET_BULK_QC_EVIDENCE.tsv"
    target_bulk_qc.to_csv(target_qc_path, sep="\t", index=False)
    target_self_lr_path = out / "FROZEN_TARGET_HOMOTYPIC_LR_NATIVE.tsv.gz"
    target_lr_scope_path = out / "FROZEN_TARGET_LR_SCOPE_SUMMARY.tsv"
    for table, path, compression in ((target_self_lr, target_self_lr_path, "gzip"), (target_lr_scope, target_lr_scope_path, None)):
        table.to_csv(path, sep="\t", index=False, compression=compression)
        restored = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)
        if list(restored.columns) != list(table.columns) or restored.to_dict("records") != table.to_dict("records"):
            raise ValueError("Frozen-target native LR/scope string values changed during report export: " + str(path))
    target_wide_provenance.update(output_path=target_wide_path.relative_to(root).as_posix(),
        output_sha256=sha(target_wide_path), readable_QC_output=target_qc_path.relative_to(root).as_posix(),
        readable_QC_output_sha256=sha(target_qc_path), readable_QC_rows=len(target_bulk_qc),
        exported_all_column_string_roundtrip_verified=True)
    write_json(out / "FROZEN_TARGET_WIDE_SELECTION_PROVENANCE.json", target_wide_provenance)
    panels.to_csv(out / "FIGURE_PANEL_SOURCE_INDEX.tsv", sep="\t", index=False)
    # Exclude the moving canonical logger trace from source-immutability checks;
    # it is deliberately snapshotted and then refreshed after this command exits.
    for rel, item in d.files.items():
        if rel == "EXECUTION_TRACE.tsv" or rel in set(trace.loc[trace.run_id.isin(unfinished), "stdout_log"]):
            continue
        if sha(d.path(rel)) != item["sha256"]:
            raise ValueError("Report source changed during generation: " + rel)
    state = {"status": "REPORTS_GENERATED_PENDING_TRACE_REFRESH", "start_iso8601": started, "finish_iso8601": now(),
             "duration_seconds": time.monotonic() - tick, "no_models_run_by_reporter": True,
             "source_hashes": d.files, "candidate_freeze_sha256": sha(d.path(args.freeze)),
             "trace_active_run_ids_at_generation": unfinished, "old_comparison_requested": args.old_comparison,
             "target_evidence_export_sha256": {p.relative_to(root).as_posix(): sha(p) for p in (
                 out / "FROZEN_TARGET_ALL_CONTINUOUS_EVIDENCE.tsv.gz", target_wide_path, target_qc_path,
                 out / "FROZEN_TARGET_WIDE_SELECTION_PROVENANCE.json", target_self_lr_path, target_lr_scope_path)},
             "native_chain_diagnostics_export": {"path": chain_diagnostic_path.relative_to(root).as_posix(),
                 "sha256": sha(chain_diagnostic_path), "donor_stage_rows": len(native_chain_diagnostics),
                 "all_native_values_retained_as_strings": True, "diagnostics_recalculated": False,
                 "new_warning_threshold_applied": False},
             "report_logger_run_ids_at_generation": report_logger_runs,
             "old_comparison_files_read_only": old_files, "reports": list(REPORTS),
             "manifest_finalisation_command": f"env/conda/bin/python scripts/write_final_reports.py --root {root} --out {args.out} --refresh-manifest-only"}
    write_json(status_path, state)
    counts_manifest = manifest(root, out, trace, known_hashes, False)
    print(json.dumps({"status": state["status"], "out": str(out), "reports": len(REPORTS), **counts_manifest}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FINAL_REPORTS_FAIL: {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        raise
