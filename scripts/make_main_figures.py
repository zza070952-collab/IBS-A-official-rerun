#!/usr/bin/env python3
"""Five auditable Python-only main figures, exclusively from the clean rerun.

No models, scientific P adjustments, target selection, or sample filters run here.
Formal rendering requires the fresh Phase6 freeze and completed fresh modules.
Use --validate-only to audit interfaces without opening a graphics device.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = "IBS_A_OFFICIAL_RERUN_20261105"
DONORS = ("P1", "P5", "P8", "P9", "P10", "P11", "P12", "P13")
CONTRASTS = ("Duo_IBSD_vs_HV", "Duo_IBSC_vs_HV", "Jej_IBSD_vs_HV", "Jej_IBSC_vs_HV")
COLORS = {"gray": "#767676", "light": "#D8D8D8", "dark": "#272727",
          "blue": "#0F4D92", "blue2": "#3775BA", "orange": "#C77932",
          "violet": "#9A4D8E", "na": "#F1F1F1", "border": "#C8C8C8"}
GROUP_COLORS = {"HV": COLORS["gray"], "IBSD": COLORS["blue"], "IBSC": COLORS["violet"]}
DOMAIN_COLORS = ("#4477AA", "#EE9977", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB", "#882255")
EVIDENCE_STATUS = {
    "SIGNIFICANT": ("#DCE8F4", "Canonical corrected significance"),
    "NONSIGNIFICANT": ("#F3E6D9", "Tested; not corrected-significant"),
    "OBSERVED_DESCRIPTIVE": ("#E9E9E9", "Observed descriptive measurement"),
    "UNTESTED": ("#FFFFFF", "Untested"),
    "NOT_OBSERVED": ("#F3F3F3", "Not observed / not measurable"),
    "NOT_APPLICABLE": ("#FFFFFF", "Not applicable"),
    "MULTIPLE_SOURCE_RECORDS": ("#F0EAF2", "Multiple source records; no arbitrary collapse"),
    "NA_NATIVE_STATISTIC": ("#F3F3F3", "Native statistic NA; original reason retained"),
}


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def need_columns(tab, columns, label):
    missing = set(columns) - set(tab.columns)
    if missing:
        raise ValueError(f"{label}: missing explicit columns {sorted(missing)}")


def bools(values):
    return pd.Series(values).astype(str).str.lower().isin(["true", "1"])


def numeric(values):
    return pd.to_numeric(values, errors="coerce")


def p_display(values, cap=16.0, floor=1e-300):
    """Display transform only; zeros/NA and every clipping flag remain explicit."""
    p = np.asarray(numeric(pd.Series(values)), dtype=float)
    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    out = np.full(len(p), np.nan)
    out[valid] = -np.log10(np.maximum(p[valid], floor))
    clipped = valid & (out > cap)
    out[valid] = np.minimum(out[valid], cap)
    return out, clipped, valid & (p == 0)


class Inputs:
    def __init__(self, root, out):
        self.root = root.resolve()
        if self.root.name != PROJECT:
            raise ValueError("Only the authorized new project can be plotted")
        self.out = self.confine(out, must_exist=False)
        try:
            self.out.relative_to(self.root / "figures")
        except ValueError:
            raise ValueError("Figure outputs must stay in this NEW project's figures directory")
        self.source = self.out / "source_data"
        self.files = {}
        self.panels = []
        self.notes = {}

    def confine(self, value, must_exist=True):
        p = Path(value)
        p = (self.root / p).resolve() if not p.is_absolute() else p.resolve()
        try:
            p.relative_to(self.root)
        except ValueError:
            raise ValueError(f"Input/output escapes NEW project: {value}")
        if must_exist and not p.is_file():
            raise FileNotFoundError(p)
        return p

    def register(self, value):
        p = self.confine(value)
        rel = p.relative_to(self.root).as_posix()
        if rel not in self.files:
            self.files[rel] = {"path": rel, "bytes": p.stat().st_size, "sha256": sha(p)}
        return p

    def js(self, value):
        return json.loads(self.register(value).read_text(encoding="utf-8-sig"))

    def tab(self, value, columns=(), **kwargs):
        p = self.register(value)
        tab = pd.read_csv(p, sep="\t", **kwargs)
        need_columns(tab, columns, str(value))
        return tab

    def export(self, panel, tab, sources, notes="", compressed=False):
        self.source.mkdir(parents=True, exist_ok=True)
        suffix = ".tsv.gz" if compressed else ".tsv"
        path = self.source / (panel + suffix)
        if not isinstance(tab, pd.DataFrame):
            tab = pd.DataFrame(tab)
        tab.to_csv(path, sep="\t", index=False, na_rep="NA", compression="gzip" if compressed else None,
                   float_format="%.17g")
        refs = []
        for value in sources:
            p = self.register(value)
            refs.append(p.relative_to(self.root).as_posix())
        self.panels.append({"panel": panel, "source_data": path.relative_to(self.root).as_posix(),
                            "source_data_sha256": sha(path), "n_source_rows": len(tab),
                            "input_paths": json.dumps(refs), "notes": notes})
        return path

    def finish(self):
        for rel, metadata in self.files.items():
            if sha(self.root / rel) != metadata["sha256"]:
                raise ValueError(f"Source changed while figures were being rendered: {rel}")
        self.source.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(self.panels).to_csv(self.source / "PANEL_PROVENANCE.tsv", sep="\t", index=False)
        write_json(self.out / "INPUT_SHA256.json", {"inputs": self.files})
        write_json(self.out / "PANEL_METHODS_AND_CAUTIONS.json", self.notes)


def validate_gates(d, freeze_file):
    global_gate = d.js("manifest/PHASE0_GATE.json")
    if global_gate.get("status") != "PASS":
        raise ValueError("Global Phase0 has not PASSed")
    gwas = d.js("manifest/GWAS_MAGMA_GATE.json")
    if gwas.get("status") != "PASS" or gwas.get("is_new_fit") is not True:
        raise ValueError("Fresh MAGMA gate is not PASS/new-fit")
    if not all(code == 0 for code in gwas.get("model_process_exit_codes", [])) or not gwas.get("model_process_exit_codes"):
        raise ValueError("MAGMA does not provide actual successful process exit codes")
    bulk = d.js("results/bulk/BULK_STATUS.json")
    if bulk.get("status") != "ANALYSIS_COMPLETE" or bulk.get("models_all_new_fit") is not True or bulk.get("n_contrasts") != 4:
        raise ValueError("Fresh four-contrast bulk branch is incomplete")
    pre = d.js("data/metadata/SPATIAL_PREPROCESS_ALL_DONORS_STATUS.json")
    if pre.get("status") != "PREPROCESS_COMPLETE" or pre.get("n_donors") != 8:
        raise ValueError("Fresh eight-donor author preprocessing is incomplete")
    if {x["donor"] for x in pre["donors"]} != set(DONORS) or any(x.get("manual_retained") != 0 for x in pre["donors"]):
        raise ValueError("Preprocessing must include the fixed8 with no author-manual excluded barcode retained")
    bayes = d.js("results/spatial/BayesSpace/BAYESSPACE_ALL_DONORS_STATUS.json")
    if bayes.get("status") != "BAYESSPACE_COMPLETE" or bayes.get("n_donors") != 8:
        raise ValueError("Fresh eight-donor BayesSpace is incomplete")
    if {x["donor"] for x in bayes["donors"]} != set(DONORS):
        raise ValueError("BayesSpace donor set differs")
    for x in bayes["donors"]:
        if x.get("nrep_formal") != 50000 or x.get("nrep_stability") != 50000 or x.get("is_new_fit") is not True:
            raise ValueError("BayesSpace formal/stability actual50k new fits are required")
    spatialdm = d.js("results/spatial/SpatialDM/SPATIALDM_ALL_DONORS_STATUS.json")
    if spatialdm.get("status") != "SPATIALDM_COMPLETE" or spatialdm.get("n_donors") != 8:
        raise ValueError("Fresh eight-donor SpatialDM is incomplete")
    if {x["donor"] for x in spatialdm["donors"]} != set(DONORS):
        raise ValueError("SpatialDM donor set differs")
    if any(x.get("n_perm_global") != 1000 or x.get("primary_P_type") != "permutation" for x in spatialdm["donors"]):
        raise ValueError("Figure requires actual primary1000-permutation global results")
    frozen = d.js(freeze_file)
    if frozen.get("status") != "FROZEN" or frozen.get("phase6_complete") is not True:
        raise ValueError("No formal figure before the fresh Phase6 candidate freeze")
    for key in ("evidence_file", "evidence_sha256", "ranking_rule_file", "ranking_rule_sha256",
                "figure_evidence_axes", "measurement_status_legend"):
        if key not in frozen or frozen[key] is None:
            raise ValueError(f"Candidate freeze is missing {key}")
    if not isinstance(frozen["figure_evidence_axes"], list) or not isinstance(frozen["measurement_status_legend"], dict):
        raise ValueError("Fresh freeze must explicitly declare its evidence-axis list and measurement-status legend")
    for path_key, hash_key in (("evidence_file", "evidence_sha256"), ("ranking_rule_file", "ranking_rule_sha256")):
        if sha(d.register(frozen[path_key])) != frozen[hash_key]:
            raise ValueError(f"Phase6 frozen {path_key} SHA256 mismatch")
    main = frozen.get("main_target")
    backup = frozen.get("backup_targets", [])
    if not isinstance(backup, list) or len(backup) > 2:
        raise ValueError("At most two freshly frozen backup targets are allowed")
    if main is None and not frozen.get("no_target_reason"):
        raise ValueError("No main target requires an explicit fresh-freeze no_target_reason")
    targets = ([main] if main is not None else []) + backup
    if main is None and backup:
        raise ValueError("Do not silently promote a backup after a no-target main freeze")
    identities = []
    for target in targets:
        for key in ("identity_key", "gene_symbol", "spatial_symbols"):
            if key not in target:
                raise ValueError(f"Fresh target identity must explicitly include {key}")
        if not isinstance(target["spatial_symbols"], list):
            raise ValueError("spatial_symbols must explicitly list verified raw symbols, even if empty")
        identities.append(target["identity_key"])
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate target identity in freeze")
    return frozen, targets, {"bulk": bulk, "pre": pre, "bayes": bayes, "spatialdm": spatialdm}


def plotting_setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        font = "Arial"
    except ValueError:
        font = "DejaVu Sans"
    matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [font],
        "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
        "ytick.labelsize": 6, "legend.fontsize": 6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6,
        "savefig.facecolor": "white", "figure.facecolor": "white", "text.color": COLORS["dark"],
        "axes.labelcolor": COLORS["dark"], "xtick.color": COLORS["dark"], "ytick.color": COLORS["dark"]})
    return plt, font


def label(ax, text, title=None):
    ax.text(-.09, 1.08, text, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
    if title:
        ax.set_title(title, loc="left", pad=7)


def no_measurement(ax, text):
    ax.set_axis_off()
    ax.text(.5, .5, text, transform=ax.transAxes, ha="center", va="center", color=COLORS["gray"], wrap=True)


def save_figure(fig, d, name, formats):
    d.out.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        path = d.out / f"{name}.{fmt}"
        # dpi controls only rasterized point collections in vector exports;
        # text/axes remain editable vector objects.
        kwargs = {"facecolor": "white", "dpi": 300}
        if fmt == "png":
            kwargs["dpi"] = 300
        elif fmt in ("tif", "tiff"):
            kwargs.update(dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        fig.savefig(path, **kwargs)


def figure1(d, frozen, targets, modules, plt):
    import textwrap
    from matplotlib.patches import FancyBboxPatch
    fig = plt.figure(figsize=(183/25.4, 170/25.4))
    grid = fig.add_gridspec(3, 2, left=.105, right=.975, bottom=.065, top=.94,
                          height_ratios=[1.0, 1.0, 1.18], hspace=.62, wspace=.30)
    ax = fig.add_subplot(grid[0, :]); label(ax, "a", "Public inputs → audited analysis units → new evidence freeze")
    ax.set_axis_off()
    gwas_source = "results/gwas/qc/summary.json"
    raw_source = "data/metadata/bulk/RAW_COUNT_VALIDATION.json"
    gwas, raw = d.js(gwas_source), d.js(raw_source)
    main_n = int(raw["n_main_libraries"])
    pre = modules["pre"]
    flow_rows = [
        {"stream": "GWAS", "stage": "raw_rows", "n": gwas["counts"]["raw_rows"], "unit": "SNP rows"},
        {"stream": "GWAS", "stage": "valid_unique_autosomal_snps", "n": gwas["counts"]["valid_unique_autosomal_snps"], "unit": "SNPs"},
        {"stream": "GWAS", "stage": "eligible_exact_ld_matches", "n": gwas["counts"]["eligible_exact_ld_matches"], "unit": "SNPs"},
        {"stream": "bulk", "stage": "count_columns", "n": raw["n_count_columns"], "unit": "tissue libraries"},
        {"stream": "bulk", "stage": "GSM_mapped", "n": raw["n_GSM"], "unit": "tissue libraries"},
        {"stream": "bulk", "stage": "primary_libraries", "n": main_n, "unit": "tissue libraries"},
        {"stream": "spatial", "stage": "UTT_input", "n": sum(x["n_input_UTT"] for x in pre["donors"]), "unit": "spots within8donors"},
        {"stream": "spatial", "stage": "primary_spots", "n": pre["total_retained"], "unit": "spots within8donors"},
        {"stream": "freeze", "stage": "main_target", "n": int(frozen.get("main_target") is not None), "unit": "target identities"},
        {"stream": "freeze", "stage": "backup_targets", "n": len(frozen.get("backup_targets", [])), "unit": "target identities"},
    ]
    texts = [
        ("IBS GWAS", f"{gwas['counts']['raw_rows']:,} raw SNP rows\n{gwas['counts']['valid_unique_autosomal_snps']:,} valid unique autosomal SNPs\n{gwas['counts']['eligible_exact_ld_matches']:,} exact EUR-LD matches", COLORS["blue"]),
        ("Patient bulk RNA-seq", f"{raw['n_count_columns']} count libraries → {raw['n_GSM']} GSM mapped\n{main_n} primary tissue libraries\nDuo/Jej modeled separately", COLORS["violet"]),
        ("Surgical small-intestine reference", f"8 donors; {sum(x['n_input_UTT'] for x in pre['donors']):,} UTT spots\n{pre['total_retained']:,} primary spots\nauthor exclusions/background applied", COLORS["orange"]),
    ]
    for i, (heading, content, color) in enumerate(texts):
        x = .005 + i*.337
        centre = x+.155
        ax.text(centre, .97, textwrap.fill(heading, width=28, break_long_words=False),
                transform=ax.transAxes, ha="center", va="top", fontsize=6.5,
                fontweight="bold", color=color, linespacing=1.1)
        ax.add_patch(FancyBboxPatch((x, .23), .31, .51, boxstyle="round,pad=.006",
                                   transform=ax.transAxes, facecolor="#FAFAFA", edgecolor=color, linewidth=.7))
        wrapped = "\n".join(textwrap.fill(line, width=29, break_long_words=False,
                                         break_on_hyphens=False) for line in content.splitlines())
        ax.text(centre, .485, wrapped, transform=ax.transAxes, ha="center", va="center",
                fontsize=6.0, linespacing=1.15)
        ax.annotate("", xy=(.5,.105), xytext=(centre,.22), xycoords="axes fraction", arrowprops={"arrowstyle": "-", "color": COLORS["border"]})
    target_text = ", ".join(x["gene_symbol"] for x in targets) if targets else "No target nominated: " + frozen["no_target_reason"]
    nomination = "Frozen computational nomination: " + target_text + ". Human-model and assay feasibility remain conditional."
    ax.text(.5, .065, textwrap.fill(nomination, width=108, break_long_words=False),
            transform=ax.transAxes, ha="center", va="top", fontsize=6.2, linespacing=1.15)
    d.export("Fig1a_flow", flow_rows, [gwas_source, raw_source, "data/metadata/SPATIAL_PREPROCESS_ALL_DONORS_STATUS.json"],
             "Arrowless stream joins are data flow, not biological causality; tissue-library counts are not subject counts.")
    group_source = "data/metadata/bulk/GROUP_COUNTS_AND_SENSITIVITY.tsv"
    groups = d.tab(group_source, ["tissue", "group", "main_analysis_pass", "published_group_count_not_identity"])
    recovery = d.js("data/metadata/bulk/OFFICIAL60_RECOVERY_STATUS.json")
    for j, tissue in enumerate(("Duo", "Jej")):
        sub = groups[groups.tissue.eq(tissue)].set_index("group").loc[["HV", "IBSD", "IBSC"]].reset_index()
        bx = fig.add_subplot(grid[1, j]); label(bx, chr(ord("b")+j), f"{tissue}: actual retained tissue libraries")
        x = np.arange(3)
        bx.bar(x, sub.main_analysis_pass, width=.55, color=[GROUP_COLORS[g] for g in sub.group])
        bx.bar(x, sub.published_group_count_not_identity, width=.75, color="none", edgecolor=COLORS["gray"], linewidth=.8, linestyle="--")
        for xx, n in zip(x, sub.main_analysis_pass):
            bx.text(xx, n+.2, str(int(n)), ha="center", fontsize=7)
        bx.set_xticks(x, sub.group.str.replace("IBS", "IBS-", regex=False)); bx.set_ylabel("Tissue libraries")
        bx.set_ylim(0, max(sub.main_analysis_pass.max(), sub.published_group_count_not_identity.max())+5)
        bx.text(.02, .96, "Dashed outline: published group total\nIndividual original60 IDs: " + ("authenticated" if recovery["status"] == "AUTHENTICATED_60_GSM_LIST_SUPPLIED" else "unrecovered"),
                transform=bx.transAxes, fontsize=5.5, va="top", linespacing=1.15)
        d.export(f"Fig1{chr(ord('b')+j)}_bulk_counts", sub, [group_source, "data/metadata/bulk/OFFICIAL60_RECOVERY_STATUS.json"], "Reference outlines do not imply the same individual library membership.")
    spatial_grid = grid[2, :].subgridspec(2, 1, height_ratios=[.23, 1], hspace=.12)
    legend_ax = fig.add_subplot(spatial_grid[0, 0]); legend_ax.set_axis_off()
    label(legend_ax, "d", "Every donor: actual author-processing spot flow")
    sx = fig.add_subplot(spatial_grid[1, 0])
    source = "data/metadata/SPATIAL_SPOT_FLOW.tsv"
    flow = d.tab(source, ["donor", "stage", "n_spots", "analysis_stage"])
    flow = flow[flow.analysis_stage.eq("phase1_author_independent_reimplementation")].copy()
    order = ["counts_UTT", "after_corrected_UMI_Z_gt_minus4", "after_mito_Z_lt4", "after_Hb_Inf_and_boundary100", "after_all_author_manual_exclusions_PRIMARY"]
    stages = [s for s in order if s in set(flow.stage)]
    if len(stages) != 5 or flow.duplicated(["donor", "stage"]).any():
        raise ValueError("Explicit Phase1 spot-flow schema changed; update the documented figure adapter")
    pivot = flow.pivot(index="donor", columns="stage", values="n_spots").reindex(DONORS)
    if pivot[stages].isna().any().any():
        raise ValueError("Missing donor/stage spot count")
    for k, stage in enumerate(stages):
        sx.plot(np.arange(8)+(k-2)*.09, pivot[stage], marker="o", ms=3, lw=.8,
                color=plt.get_cmap("Blues")(.25+.14*k), label=["UTT input", "Corrected UMI QC", "Mito QC", "Hb/boundary", "All manual exclusions"][k])
    sx.set_xticks(np.arange(8), DONORS); sx.set_ylabel("Spots per donor")
    handles, legend_names = sx.get_legend_handles_labels()
    legend_ax.legend(handles, legend_names, ncol=3, loc="center left", fontsize=5.6,
                     frameon=False, borderaxespad=0, columnspacing=1.1, handletextpad=.5)
    d.export("Fig1d_spatial_flow", flow, [source], "Spots are subsamples; biological reference n=8 donors; all actual counts retained.")
    d.notes["Fig1"] = {"unit": "SNPs; tissue libraries; spots nested within8 surgical donors", "statistics": "Descriptive counts only", "clinical_caution": "Spatial reference is surgical tissue, not IBS patient replication or unscreened healthy volunteers"}
    return fig


def manhattan_coordinates(tab, chromosome, bp):
    out = tab.copy()
    out["plot_chr"] = numeric(out[chromosome].astype(str).str.replace("chr", "", regex=False))
    out["plot_bp"] = numeric(out[bp])
    if out[["plot_chr", "plot_bp"]].isna().any().any():
        raise ValueError("Unrecognized chromosome/position; no guessed genome coordinates")
    if not out.plot_chr.between(1, 22).all():
        raise ValueError("Expected the prespecified autosomal plotting universe")
    offsets, centers, offset = {}, {}, 0.0
    for chrom in sorted(out.plot_chr.unique()):
        length = out.loc[out.plot_chr.eq(chrom), "plot_bp"].max()
        offsets[chrom], centers[chrom] = offset, offset + length/2
        offset += length + 5e6
    out["plot_cumulative_bp"] = out.plot_bp + out.plot_chr.map(offsets)
    return out, centers


def figure2(d, plt):
    fig = plt.figure(figsize=(183/25.4, 160/25.4), layout="constrained")
    grid = fig.add_gridspec(2, 3, width_ratios=[1.25, 1.25, 1])
    ax = fig.add_subplot(grid[0, :2]); qqax = fig.add_subplot(grid[0, 2])
    gx = fig.add_subplot(grid[1, :2]); gqq = fig.add_subplot(grid[1, 2])
    source = "results/gwas/qc/manhattan_points.tsv.gz"
    raw = d.tab(source, ["SNP", "CHR", "BP", "P", "sampling_reason"])
    summary = d.js("results/gwas/qc/summary.json")
    points, centers = manhattan_coordinates(raw, "CHR", "BP")
    points["plot_minus_log10_P"], points["plot_y_clipped"], points["P_is_zero"] = p_display(points.P, cap=350)
    if points.P_is_zero.any() or not np.isfinite(points.plot_minus_log10_P).all():
        raise ValueError("Raw GWAS plot input contains invalid/zero P; inspect the Phase0 cache")
    for chrom, part in points.groupby("plot_chr"):
        ax.scatter(part.plot_cumulative_bp/1e6, part.plot_minus_log10_P, s=2,
                   color=COLORS["blue"] if int(chrom)%2 else COLORS["gray"], linewidths=0, rasterized=True)
    ax.axhline(-math.log10(5e-8), color=COLORS["orange"], ls="--", lw=.7)
    ax.set_xticks([v/1e6 for v in centers.values()], [str(int(k)) for k in centers], rotation=90)
    ax.tick_params(axis="x", labelsize=5.5, pad=2); ax.set_ylabel("−log10 SNP P")
    ax.set_xlabel("Autosomal position (GRCh37; chromosome)")
    label(ax,"a", "IBS GWAS: all genome-wide significant SNPs\n+ deterministic background")
    ax.text(.01,.98,f"{summary['counts']['valid_unique_autosomal_snps']:,} valid unique SNPs; {len(raw):,} drawn\nBackground CRC32 modulo {summary['plotting']['manhattan_modulus']}; every P<5×10⁻⁸ retained", transform=ax.transAxes,va="top",fontsize=5.5)
    d.export("Fig2a_GWAS_Manhattan", points, [source,"results/gwas/qc/summary.json"], "Display cache selection fixed in new raw QC; not all valid SNPs are drawn.", compressed=True)
    qsource = "results/gwas/qc/qq_quantiles.tsv"
    qq = d.tab(qsource,["rank","total_tests","expected_P","observed_P"])
    qq["plot_expected"] = -np.log10(numeric(qq.expected_P)); qq["plot_observed"] = -np.log10(numeric(qq.observed_P))
    qqax.scatter(qq.plot_expected,qq.plot_observed,s=3,color=COLORS["blue"],linewidths=0,rasterized=True)
    limit=max(qq.plot_expected.max(),qq.plot_observed.max()); qqax.plot([0,limit],[0,limit],color=COLORS["border"],lw=.8)
    qqax.set_xlabel("Expected −log10 P"); qqax.set_ylabel("Observed −log10 P"); label(qqax,"b","GWAS QQ")
    qqax.text(.03,.98,f"n={int(qq.total_tests.iloc[0]):,} tests\nSelected order statistics",transform=qqax.transAxes,va="top",fontsize=5.5)
    d.export("Fig2b_GWAS_QQ",qq,[qsource],"Precomputed full-universe order statistics; no confidence band assuming independent LD-correlated tests.")
    msource="results/gwas/magma_all_genes.tsv"
    genes=d.tab(msource,["gene_id","CHR","START","STOP","P","BH_FDR","rank","gene_symbol"])
    genes["gene_midpoint_bp"]=(numeric(genes.START)+numeric(genes.STOP))/2
    genes,centers=manhattan_coordinates(genes,"CHR","gene_midpoint_bp")
    genes["plot_minus_log10_P"],genes["plot_y_clipped"],genes["P_is_zero"] = p_display(genes.P,cap=350)
    if genes.P_is_zero.any() or not np.isfinite(genes.plot_minus_log10_P).all():
        raise ValueError("MAGMA has invalid P in canonical new results")
    significant=numeric(genes.BH_FDR)<.05
    gx.scatter(genes.plot_cumulative_bp/1e6,genes.plot_minus_log10_P,s=3,c=np.where(significant,COLORS["orange"],COLORS["gray"]),linewidths=0,rasterized=True)
    gx.axhline(-math.log10(.05/len(genes)),lw=.7,ls="--",color=COLORS["blue"])
    gx.set_xticks([v/1e6 for v in centers.values()],[str(int(k)) for k in centers], rotation=90)
    gx.tick_params(axis="x", labelsize=5.5, pad=2); gx.set_ylabel("−log10 gene P"); gx.set_xlabel("Gene-body midpoint (GRCh37; chromosome)")
    label(gx,"c",f"New MAGMA: {len(genes):,} tested genes; {int(significant.sum())} BH-FDR<0.05")
    gx.text(.01,.98,"Orange: BH-FDR<0.05\nDashed: Bonferroni0.05 / tested genes",transform=gx.transAxes,va="top",fontsize=5.5)
    d.export("Fig2c_MAGMA_Manhattan",genes,[msource],"Full tested-gene universe; gene association does not establish causality or independent loci.",compressed=True)
    qqg=genes.sort_values(["P","gene_id"],kind="mergesort").copy()
    qqg["plot_expected"]=-np.log10((np.arange(len(qqg))+.5)/len(qqg)); qqg["plot_observed"]=-np.log10(numeric(qqg.P))
    gqq.scatter(qqg.plot_expected,qqg.plot_observed,s=3,color=COLORS["gray"],linewidths=0,rasterized=True)
    limit=max(qqg.plot_expected.max(),qqg.plot_observed.max()); gqq.plot([0,limit],[0,limit],color=COLORS["border"],lw=.8)
    gqq.set_xlabel("Expected −log10 P"); gqq.set_ylabel("Observed −log10 P"); label(gqq,"d","MAGMA gene QQ")
    d.export("Fig2d_MAGMA_QQ",qqg,[msource],"All tested genes; no independent-test interval.",compressed=True)
    d.notes["Fig2"]={"statistics":"Canonical raw SNP P and new MAGMA gene P/BH FDR; plotting does not re-adjust", "raw_sampling":summary["plotting"],"genome_build":"GRCh37","inferential_unit":"SNP/gene association tests; no causal-gene or independent-locus claim"}
    return fig


def bulk_qc_panel_inputs(d):
    """Read native descriptive comparisons only; no gene test/join or refit."""
    main_source = "results/bulk/main/BULK_CONTRAST_SUMMARY.tsv"
    qc_source = "results/bulk/sensitivity/QC_SENSITIVITY_COMPARISON_SUMMARY.tsv"
    sample_source = "data/metadata/bulk/SAMPLE_MANIFEST.tsv"
    protocol_source = "data/metadata/bulk/QC_PROTOCOL_FROZEN.json"
    main_status_source = "results/bulk/BULK_STATUS.json"
    sensitivity_status_source = "results/bulk/sensitivity/BULK_SENSITIVITY_STATUS.json"
    sources = [main_source, qc_source, sample_source, protocol_source, main_status_source, sensitivity_status_source]
    main = d.tab(main_source, ["contrast", "n_samples", "n_DEG_FDR05"], dtype=str, keep_default_na=False)
    required = ["sensitivity", "contrast", "n_features_union", "n_finite_LFC_both", "LFC_spearman",
                "LFC_sign_agreement_fraction", "n_main_FDR05", "n_sensitivity_FDR05", "n_both_FDR05",
                "n_main_significant_sensitivity_NA"]
    comparison = d.tab(qc_source, required, dtype=str, keep_default_na=False)
    protocol = d.js(protocol_source)
    primary = d.js(main_status_source)
    sensitivity = d.js(sensitivity_status_source)
    labels = sorted(protocol["sensitivity_min_non_rrna_counts"], key=protocol["sensitivity_min_non_rrna_counts"].get)
    if len(labels) != 2 or set(sensitivity.get("analyses", {})) != set(labels):
        raise ValueError("All and only the two frozen QC-sensitivity cohorts are required")
    if primary.get("status") != "ANALYSIS_COMPLETE" or primary.get("models_all_new_fit") is not True or sensitivity.get("status") != "ANALYSIS_COMPLETE":
        raise ValueError("Main and fixed bulk-QC native analyses must be complete before Figure3")
    selection = {"main": primary["selection_column"], **{label: label for label in labels}}
    samples = d.tab(sample_source, ["library_id", "tissue"] + list(selection.values()), dtype=str, keep_default_na=False)
    if samples.library_id.duplicated().any() or main.contrast.duplicated().any() or set(main.contrast) != set(CONTRASTS):
        raise ValueError("Cohort membership or native main contrast identity is not unique/complete")
    expected = {(label, contrast) for label in labels for contrast in CONTRASTS}
    if comparison.duplicated(["sensitivity", "contrast"]).any() or len(comparison) != len(expected) or set(zip(comparison.sensitivity, comparison.contrast)) != expected:
        raise ValueError("All eight native QC comparison rows are mandatory; no favorable subset")

    def count(value, field):
        number = float(value)
        if not math.isfinite(number) or number < 0 or number != int(number):
            raise ValueError("Invalid native QC count: " + field)
        return int(number)

    sizes = {}
    for cohort, column in selection.items():
        if not samples[column].str.lower().isin(["true", "false", "1", "0"]).all():
            raise ValueError("Unknown frozen cohort-membership flag: " + column)
        mask = bools(samples[column])
        sizes[cohort] = {"total": int(mask.sum()), **{t: int((mask & samples.tissue.eq(t)).sum()) for t in ("Duo", "Jej")}}
        if cohort != "main" and count(sensitivity["analyses"][cohort]["n_libraries"], "n_libraries") != sizes[cohort]["total"]:
            raise ValueError("Native sensitivity and frozen membership library totals disagree")
    main_by = main.set_index("contrast")
    rows = []
    for contrast in CONTRASTS:
        raw = main_by.loc[contrast].to_dict()
        tissue = contrast.split("_", 1)[0]
        if count(raw["n_samples"], "n_samples") != sizes["main"][tissue]:
            raise ValueError("Actual main contrast tissue-library count disagrees with frozen membership")
        rows.append({"cohort": "main", "contrast": contrast, "n_cohort_libraries_total": sizes["main"]["total"],
                     "n_contrast_tissue_libraries": sizes["main"][tissue], "n_FDR05": count(raw["n_DEG_FDR05"], "n_DEG_FDR05"),
                     "n_intersection_with_main": "NOT_APPLICABLE_REFERENCE", "LFC_spearman": "NOT_APPLICABLE_REFERENCE",
                     "LFC_sign_agreement_fraction": "NOT_APPLICABLE_REFERENCE", "n_main_significant_sensitivity_NA": "NOT_APPLICABLE_REFERENCE",
                     "source_path": main_source, "source_record_key": contrast, "native_comparison": "PRIMARY_REFERENCE_NOT_A_RETEST"})
    for _, raw in comparison.iterrows():
        cohort, contrast = raw["sensitivity"], raw["contrast"]
        main_n = count(raw["n_main_FDR05"], "n_main_FDR05")
        sensitivity_n = count(raw["n_sensitivity_FDR05"], "n_sensitivity_FDR05")
        overlap = count(raw["n_both_FDR05"], "n_both_FDR05")
        q_na = count(raw["n_main_significant_sensitivity_NA"], "n_main_significant_sensitivity_NA")
        if main_n != count(main_by.loc[contrast, "n_DEG_FDR05"], "n_DEG_FDR05") or overlap > min(main_n, sensitivity_n) or q_na > main_n:
            raise ValueError("Native main/sensitivity FDR counts, intersection or q-NA counts disagree")
        for field, lower, upper in (("LFC_spearman", -1, 1), ("LFC_sign_agreement_fraction", 0, 1)):
            value = pd.to_numeric(raw[field], errors="coerce")
            if pd.notna(value) and (not np.isfinite(value) or not lower <= value <= upper):
                raise ValueError("Out-of-range native continuous agreement: " + field)
        rows.append({"cohort": cohort, "contrast": contrast, "n_cohort_libraries_total": sizes[cohort]["total"],
                     "n_contrast_tissue_libraries": sizes[cohort][contrast.split("_", 1)[0]], "n_FDR05": sensitivity_n,
                     "n_intersection_with_main": overlap, "LFC_spearman": raw["LFC_spearman"],
                     "LFC_sign_agreement_fraction": raw["LFC_sign_agreement_fraction"],
                     "n_main_significant_sensitivity_NA": raw["n_main_significant_sensitivity_NA"],
                     "source_path": qc_source, "source_record_key": cohort + ";" + contrast,
                     "native_comparison": "ORIGINAL_NATIVE_DESCRIPTIVE_COMPARISON_NOT_RECALCULATED"})
    return pd.DataFrame(rows), main, comparison, samples, sizes, labels, protocol, sources


def figure3(d,targets,plt):
    import textwrap
    from matplotlib.colors import TwoSlopeNorm,LinearSegmentedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    qc,main_summary,qc_native,cohort_samples,cohort_sizes,qc_labels,bulk_protocol,qc_sources=bulk_qc_panel_inputs(d)
    fig=plt.figure(figsize=(183/25.4,260/25.4))
    # Separate blank rows reserve real space for axis labels and panel titles.
    grid=fig.add_gridspec(7,4,left=.105,right=.955,bottom=.045,top=.94,
                         height_ratios=[1.05,.60,1.18,.65,1.18,.20,2.80],
                         hspace=0,wspace=.62)
    sex_handles=[Line2D([],[],marker=marker,ls="",color=COLORS["gray"],
                        label="sex "+sex,ms=3.5)
                 for sex,marker in (("F","o"),("M","^"),("NA","s"))]
    fig.legend(handles=sex_handles,loc="upper center",bbox_to_anchor=(.55,.991),
               ncol=3,frameon=False,fontsize=5.5,handletextpad=.35,columnspacing=1.1)
    for i,tissue in enumerate(("Duo","Jej")):
        src=f"results/bulk/main/{tissue}_PCA_SAMPLE_QC.tsv"; varsrc=f"results/bulk/main/{tissue}_PCA_VARIANCE.tsv"
        pc=d.tab(src,["library_id","PC1","PC2","group","sex"]); var=d.tab(varsrc,["PC","variance_fraction"])
        ax=fig.add_subplot(grid[0,i*2:(i+1)*2]);label(ax,chr(97+i),f"{tissue}: blind VST PCA, n={len(pc)} tissue libraries")
        for group in ("HV","IBSD","IBSC"):
            for sex,marker in (("F","o"),("M","^"),("NA","s")):
                part=pc[pc.group.eq(group) & pc.sex.fillna("NA").eq(sex)]
                if len(part): ax.scatter(part.PC1,part.PC2,s=15,c=GROUP_COLORS[group],marker=marker,edgecolors="white",linewidths=.2)
        explained=var.set_index("PC").variance_fraction
        ax.set_xlabel(f"PC1 ({float(explained['PC1'])*100:.1f}% variance)");ax.set_ylabel(f"PC2 ({float(explained['PC2'])*100:.1f}% variance)")
        handles=[Line2D([],[],marker="o",ls="",color=GROUP_COLORS[g],label=g.replace("IBS","IBS-"),ms=4) for g in ("HV","IBSD","IBSC")]
        ax.legend(handles=handles,loc="best",ncol=3,fontsize=5.6,
                  handletextpad=.35,columnspacing=.8,borderpad=.35)
        d.export(f"Fig3{chr(97+i)}_{tissue}_PCA",pc,[src,varsrc],"Blind VST/top500 variable genes; no automatic deletion of flagged PCA samples; marker circle=F triangle=M square=unknown.")
    gsea_tables=[]
    for i,contrast in enumerate(CONTRASTS):
        src=f"results/bulk/main/{contrast}.tsv";tab=d.tab(src,["gene_symbol","gene_id","log2FoldChange","padj","P_NA_reason","padj_NA_reason"])
        tab["plot_x"]=numeric(tab.log2FoldChange).clip(-6,6)
        tab["plot_x_clipped"]=numeric(tab.log2FoldChange).abs()>6
        tab["plot_y"],tab["plot_y_clipped"],tab["padj_is_zero"] = p_display(tab.padj,cap=16)
        tab["plotted"]=np.isfinite(tab.plot_x)&np.isfinite(tab.plot_y)
        significant=numeric(tab.padj)<.05
        ax=fig.add_subplot(grid[2,i]);label(ax,chr(99+i),contrast.replace("_vs_HV"," vs HV").replace("_"," ").replace("IBS","IBS-"))
        part=tab[tab.plotted]
        ax.scatter(part.plot_x,part.plot_y,c=np.where(significant[tab.plotted],COLORS["orange"],COLORS["light"]),s=2,linewidths=0,rasterized=True)
        cap=part.plot_x_clipped|part.plot_y_clipped
        if cap.any(): ax.scatter(part.loc[cap,"plot_x"],part.loc[cap,"plot_y"],s=6,marker="^",c=COLORS["gray"],linewidths=0,rasterized=True)
        ax.axhline(-math.log10(.05),color=COLORS["gray"],lw=.5,ls="--");ax.axvline(0,color=COLORS["border"],lw=.5)
        ax.set_xlim(-6.5,6.5);ax.set_ylim(-.15,16.8);ax.set_xlabel("log2 fold change");ax.set_ylabel("−log10 BH-FDR" if i==0 else "")
        ax.text(.02,.98,f"q<0.05: {int(significant.sum())}\nMissing q: {int(tab.padj.isna().sum())}\nDisplay caps: ±6 / 16",transform=ax.transAxes,va="top",fontsize=5)
        target_rows=[]
        for target in targets:
            symbols=target.get("bulk_symbols",[target["gene_symbol"]])
            hits=tab[tab.gene_symbol.isin(symbols) & tab.plotted]
            for _,r in hits.iterrows():
                ax.scatter([r.plot_x],[r.plot_y],s=14,facecolors="none",edgecolors=COLORS["blue"],lw=.6)
                ax.annotate(target["gene_symbol"],(r.plot_x,r.plot_y),xytext=(2,3),textcoords="offset points",fontsize=4.7)
                target_rows.append(r.gene_id)
        d.export(f"Fig3{chr(99+i)}_{contrast}_volcano",tab,[src],"All features retained; only finite LFC/q plotted; q NA not imputed; arrows mark visual caps; no candidate-dependent filter.",compressed=True)
        gs=f"results/bulk/main/{contrast}_GSEA_ALL.tsv";gsea=d.tab(gs,["collection","pathway","NES","padj","tested","NA_reason"])
        gsea["contrast"]=contrast;gsea_tables.append(gsea)
    qc_grid=grid[4,:].subgridspec(2,1,height_ratios=[1,.27],hspace=.12)
    qax=fig.add_subplot(qc_grid[0,0]);label(qax,"g","Fixed QC cohorts: corrected-significance sets and continuous LFC agreement")
    qax.set_axis_off();qax.set_xlim(-1.65,3.54);qax.set_ylim(2.54,-1.08)
    qc_note=fig.add_subplot(qc_grid[1,0]);qc_note.set_axis_off()
    cohort_order=["main"]+qc_labels
    qc_lookup=qc.set_index(["cohort","contrast"])
    row_colors=["#EDF3FA","#F5F5F5","#FBF3EA"]
    for x,contrast in enumerate(CONTRASTS):
        column_name=contrast.replace("_vs_HV","").replace("_"," ").replace("IBS","IBS-")+"\nvs HV"
        qax.text(x,-.78,column_name,ha="center",va="center",fontsize=5.8,linespacing=1.05)
    for y,cohort in enumerate(cohort_order):
        name="Primary" if cohort=="main" else f"QC ≥{bulk_protocol['sensitivity_min_non_rrna_counts'][cohort]/1e6:g}M"
        count_label=f"n={cohort_sizes[cohort]['total']} libraries"
        qax.text(-.58,y,name+"\n"+count_label,ha="right",va="center",fontsize=5.7,linespacing=1.1)
        for x,contrast in enumerate(CONTRASTS):
            row=qc_lookup.loc[(cohort,contrast)]
            qax.add_patch(Rectangle((x-.48,y-.44),.96,.88,facecolor=row_colors[y],edgecolor="white",lw=1))
            if cohort=="main":
                display=f"{int(row.n_FDR05)}\nprimary q<0.05 set"
            else:
                rho=pd.to_numeric(row.LFC_spearman,errors="coerce")
                agreement=pd.to_numeric(row.LFC_sign_agreement_fraction,errors="coerce")
                rho_label=f"{rho:.3f}" if pd.notna(rho) else "NA"
                sign_label=f"{agreement*100:.1f}%" if pd.notna(agreement) else "NA"
                display=(f"{int(row.n_FDR05)} (∩ main {int(row.n_intersection_with_main)})\n"
                         f"ρ {rho_label} | sign {sign_label}\n"
                         f"main hits q NA: {row.n_main_significant_sensitivity_NA}")
            qax.text(x,y,display,ha="center",va="center",fontsize=5.7,linespacing=1.12)
    qc_note.text(0,1,"Counts: native q<0.05; ∩: overlap with primary. Direction agreement does not establish robust FDR membership.\n"
             "ρ/sign: common finite LFCs; q NA remains unavailable. No refit or primary-cohort/target reselection.",
             transform=qc_note.transAxes,ha="left",va="top",fontsize=5.3,linespacing=1.15)
    d.export("Fig3g_QC_cohort_comparison",qc,qc_sources,
             "All three frozen cohorts/four contrasts; original native FDR counts, overlap, LFC correlation/sign agreement and primary-hit q-NA. No refit/retest, new sample choice, significance adjustment or target selection.")
    d.export("Fig3g_QC_native_main_summary_ALL_COLUMNS",main_summary,[qc_sources[0]],"Every original native main summary column and all four rows; no statistical recalculation.")
    d.export("Fig3g_QC_native_comparison_ALL_ROWS_COLUMNS",qc_native,[qc_sources[1]],"All eight original native sensitivity rows/columns, including negative/NA evidence; values read as strings.")
    d.export("Fig3g_QC_frozen_cohort_membership",cohort_samples,[qc_sources[2]],"Full frozen sample manifest; cohort counts are tissue libraries, not verified unique participants.")
    allg=pd.concat(gsea_tables,ignore_index=True)
    chosen=set()
    for contrast in CONTRASTS:
        sub=allg[allg.contrast.eq(contrast)&bools(allg.tested).to_numpy()&np.isfinite(numeric(allg.NES))&np.isfinite(numeric(allg.padj))]
        for sign in (1,-1):
            picked=sub[numeric(sub.NES)*sign>0].sort_values(["padj","collection","pathway"],kind="mergesort").head(2)
            chosen.update(zip(picked.collection,picked.pathway))
    gsea_host=fig.add_subplot(grid[6,:]);gsea_host.set_axis_off()
    label(gsea_host,"h","Full-rank GSEA: fixed top-two per NES sign and contrast\nAll contrast outcomes retained")
    if not chosen:
        no_measurement(gsea_host,"No pathway has a finite eligible NES and primary pooled FDR; inspect full GSEA NA reasons.")
        d.export("Fig3h_GSEA_selected",allg.iloc[0:0], [f"results/bulk/main/{c}_GSEA_ALL.tsv" for c in CONTRASTS],"Explicit no-finite-pathway outcome; no replacement with fabricated pathways.")
    else:
        subset=allg[[tuple(x) in chosen for x in allg[["collection","pathway"]].itertuples(index=False,name=None)]].copy()
        keys=sorted(chosen,key=lambda k:(numeric(subset[(subset.collection.eq(k[0]))&subset.pathway.eq(k[1])].padj).min(),k))
        if subset.duplicated(["contrast","collection","pathway"]).any() or len(subset)!=len(keys)*4:
            raise ValueError("GSEA display union must contain all four explicit contrast/pathway records")
        ylookup={key:i for i,key in enumerate(keys)};xlookup={c:i for i,c in enumerate(CONTRASTS)}
        subset["plot_x"]=subset.contrast.map(xlookup);subset["plot_y"]=[ylookup[x] for x in zip(subset.collection,subset.pathway)]
        subset["plot_size_logq"],subset["plot_size_clipped"],subset["padj_is_zero"] = p_display(subset.padj,cap=8)
        finite=np.isfinite(numeric(subset.NES))&np.isfinite(subset.plot_size_logq)
        norm=TwoSlopeNorm(vcenter=0,vmin=-max(1,float(numeric(subset.NES).abs().max())),vmax=max(1,float(numeric(subset.NES).abs().max())))
        cmap=LinearSegmentedColormap.from_list("signed_effect",[COLORS["blue"],"#FFFFFF",COLORS["orange"]])
        gsea_grid=grid[6,:].subgridspec(2,3,height_ratios=[1,.16],
                                     width_ratios=[2.45,3.60,.14],wspace=.08,hspace=.17)
        name_ax=fig.add_subplot(gsea_grid[0,0]);name_ax.set_axis_off()
        ax=fig.add_subplot(gsea_grid[0,1]);cax=fig.add_subplot(gsea_grid[0,2])
        legend_ax=fig.add_subplot(gsea_grid[1,:]);legend_ax.set_axis_off()
        names=[k[1].removeprefix("HALLMARK_").removeprefix("REACTOME_").replace("_"," ") for k in keys]
        wrapped=[textwrap.fill(n,width=42,break_long_words=False,break_on_hyphens=False) for n in names]
        # Expand categorical label spacing only. The exported plot_y indices,
        # pathway order, every record, statistic, dot area and colour stay fixed.
        row_heights=np.array([max(2,name.count("\n")+1)+.50 for name in wrapped])
        row_centres=np.cumsum(row_heights)-row_heights/2
        display_y={i:value for i,value in enumerate(row_centres)}
        name_ax.set_xlim(0,1);name_ax.set_ylim(row_heights.sum(),0)
        for y,name in zip(row_centres,wrapped):
            name_ax.text(1,y,name,ha="right",va="center",fontsize=5.6,linespacing=1.08)
        points=ax.scatter(subset.loc[finite,"plot_x"],subset.loc[finite,"plot_y"].map(display_y),c=numeric(subset.loc[finite,"NES"]),
            s=8+subset.loc[finite,"plot_size_logq"]*12,cmap=cmap,norm=norm,edgecolors=np.where(numeric(subset.loc[finite,"padj"])<.05,COLORS["dark"],"none"),linewidths=.5)
        if (~finite).any():ax.scatter(subset.loc[~finite,"plot_x"],subset.loc[~finite,"plot_y"].map(display_y),marker="x",s=18,c=COLORS["gray"],lw=.6)
        ax.set_yticks([]);ax.set_ylim(row_heights.sum(),0)
        ax.set_xticks(range(4),[c.replace("_vs_HV","").replace("_","\n").replace("IBS","IBS-") for c in CONTRASTS]);ax.set_xlim(-.6,3.6)
        ax.tick_params(axis="x",labelsize=5.7)
        fig.colorbar(points,cax=cax)
        cax.set_title("NES",fontsize=6,pad=5)
        cax.tick_params(labelsize=5.5)
        legend_ax.text(0,1,"Bubble area scales with −log10 pooled FDR (display cap 8). Dark outline: q<0.05.\n"
                       "×: untested/NA. Colour: NES for IBS subgroup vs HV.",
                       transform=legend_ax.transAxes,ha="left",va="top",fontsize=5.6,linespacing=1.15)
        d.export("Fig3h_GSEA_selected",subset,[f"results/bulk/main/{c}_GSEA_ALL.tsv" for c in CONTRASTS],"Top-two per sign per contrast by original pooled FDR; no subset FDR recalculation; complete outcomes shown.")
    d.notes["Fig3"]={"unit":"tissue libraries, not independent participants across tissues","tests":"New DESeq2 Wald with canonical BH; new fgseaMultilevel with pooled-across-two-collections BH per contrast","volcano_caps":{"abs_log2FC":6,"minus_log10_FDR":16},"GSEA_size_cap":8,"missing":"Excluded/NA reasons preserved in full source tables; not assigned zero",
                     "figure_size_mm":[183,260],"fixed_QC_panel":"g; complete native counts/primary overlap/LFC agreement/q-NA; not re-estimated",
                     "layout_only":"Explicit row gaps; independent QC/GSEA annotation axes; full pathway labels wrap without truncation; GSEA categorical row spacing expands by label lines while exported plot_y indices and all data stay unchanged; shared sex-shape key F=circle/M=triangle/NA=square matches unchanged PCA markers",
                     "fixed_QC_cohort_sizes":cohort_sizes,"fixed_QC_thresholds":bulk_protocol["sensitivity_min_non_rrna_counts"],
                     "fixed_QC_boundary":"LFC rho/sign agreement is among common finite LFCs. Continuous agreement and corrected-significance membership are different; primary DEGs are not all called robust; q NA is unavailable, not zero; no primary cohort or target changed"}
    return fig


def target_expression(d,target,donor):
    from scipy.io import mmread
    from scipy import sparse
    folder=f"data/processed/spatial/{donor}"
    genepath=folder+"/genes.tsv";spotpath=folder+"/spots.tsv";matrixpath=folder+"/author_background_corrected_genes_by_spots.mtx"
    genes=d.tab(genepath,["gene_symbol","feature_id"]);spots=d.tab(spotpath,["barcode","pixel_col","pixel_row"])
    symbols=target["spatial_symbols"]
    idx=np.flatnonzero(genes.gene_symbol.isin(symbols).to_numpy())
    if not len(idx):
        return None,spots,[genepath,spotpath],"NO_VERIFIED_RAW_SYMBOL_ROW_MEASURABLE"
    matrix=mmread(d.register(matrixpath)).tocsr()
    if matrix.shape!=(len(genes),len(spots)):
        raise ValueError("Candidate expression matrix/barcode/feature schema differs")
    totals=np.asarray(matrix.sum(axis=0)).ravel()
    if (totals<=0).any() or not np.isfinite(totals).all():
        raise ValueError("Invalid corrected spot totals; no fallback to raw UMI")
    value=np.asarray(matrix[idx,:].sum(axis=0)).ravel()
    out=spots.copy();out["identity_key"]=target["identity_key"];out["gene_symbol"]=target["gene_symbol"]
    out["raw_symbols_from_freeze"]=";".join(symbols);out["n_original_feature_rows_aggregated"]=len(idx)
    out["corrected_counts_aggregated_after_background_clip"]=value
    out["plot_log1p_CP10000"]=np.log1p(10000*value/totals)
    return out,spots,[genepath,spotpath,matrixpath],"MEASURED_BACKGROUND_CORRECTED"


def pixel_map(ax,tab,values,cmap=None,norm=None,colors=None,size=2):
    need_columns(tab,["pixel_col","pixel_row"],"pixel map")
    pts=ax.scatter(tab.pixel_col,tab.pixel_row,c=colors if colors is not None else values,s=size,cmap=cmap,norm=norm,linewidths=0,rasterized=True)
    ax.invert_yaxis();ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_visible(False)
    return pts


def figure4(d,targets,modules,plt):
    from matplotlib.colors import ListedColormap,BoundaryNorm
    from matplotlib.patches import Patch
    fig=plt.figure(figsize=(183/25.4,220/25.4),layout="constrained")
    grid=fig.add_gridspec(4,4,height_ratios=[.85,1.15,1.15,1.0])
    status={x["donor"]:x for x in modules["bayes"]["donors"]};pre={x["donor"]:x for x in modules["pre"]["donors"]}
    qc=pd.DataFrame([{"donor":donor,"UTT_input":pre[donor]["n_input_UTT"],"primary_spots":pre[donor]["n_retained"],
        "q":status[donor]["q"],"ARI":status[donor]["independent_seed_ARI"],"stability_warning":status[donor]["stability_warning"],
        "nrep_formal":status[donor]["nrep_formal"],"nrep_stability":status[donor]["nrep_stability"]} for donor in DONORS])
    ax=fig.add_subplot(grid[0,:2]);label(ax,"a","Audited spot cohorts: every donor")
    x=np.arange(8);ax.bar(x,qc.UTT_input,color=COLORS["light"],width=.65,label="UTT input");ax.bar(x,qc.primary_spots,color=COLORS["blue2"],width=.38,label="Primary")
    ax.set_xticks(x,DONORS,rotation=45);ax.set_ylabel("Spots");ax.legend(ncol=2,fontsize=5.5)
    d.export("Fig4a_spot_QC",qc,["data/metadata/SPATIAL_PREPROCESS_ALL_DONORS_STATUS.json","results/spatial/BayesSpace/BAYESSPACE_ALL_DONORS_STATUS.json"],"Spots nested in8donors; no spot-level biological replication claim.")
    ar=fig.add_subplot(grid[0,2:]);label(ar,"b","Independent 50k-chain stability\nDonor-local q")
    ar.scatter(x,qc.ARI,c=np.where(qc.ARI<.6,COLORS["orange"],COLORS["blue"]),s=16)
    ar.axhline(.6,color=COLORS["gray"],ls="--",lw=.6)
    for xx,r in qc.iterrows(): ar.text(xx,min(1.04,float(r.ARI)+.04),f"q={int(r.q)}",ha="center",fontsize=5.3)
    ar.set_xticks(x,DONORS,rotation=45);ar.set_ylim(-.05,1.12);ar.set_ylabel("Independent-seed ARI")
    d.export("Fig4b_stability",qc,["results/spatial/BayesSpace/BAYESSPACE_ALL_DONORS_STATUS.json"],"ARI<0.6 is the saved protocol warning, not proof all other maps are biologically validated.")
    cmap=ListedColormap(DOMAIN_COLORS);norm=BoundaryNorm(np.arange(.5,9.5),8)
    for i,donor in enumerate(DONORS):
        src=f"results/spatial/BayesSpace/{donor}/domains.tsv"
        tab=d.tab(src,["barcode","pixel_col","pixel_row","spatial.cluster"])
        if tab.barcode.duplicated().any() or len(tab)!=int(status[donor]["n_spots"]):raise ValueError("Domain map row count/uniqueness differs from actual model")
        donor_grid=grid[1+i//4,i%4].subgridspec(3,1,height_ratios=[.20,1,.34],hspace=.02)
        head=fig.add_subplot(donor_grid[0]);head.set_axis_off()
        head.text(0,1,chr(99+i),transform=head.transAxes,fontsize=8,fontweight="bold",va="top")
        head.text(.13,1,f"{donor} | q={status[donor]['q']}\nARI={status[donor]['independent_seed_ARI']:.2f}",transform=head.transAxes,fontsize=6.1,va="top",linespacing=1.1)
        ax=fig.add_subplot(donor_grid[1])
        pixel_map(ax,tab,numeric(tab["spatial.cluster"]),cmap,norm,size=2)
        foot=fig.add_subplot(donor_grid[2]);foot.set_axis_off()
        foot.text(.5,1,f"n={len(tab):,} spots | donor-local domain",transform=foot.transAxes,fontsize=4.9,ha="center",va="top")
        domain_handles=[Patch(facecolor=DOMAIN_COLORS[j-1],edgecolor="none",label=str(j)) for j in range(1,int(status[donor]["q"])+1)]
        foot.legend(handles=domain_handles,ncol=4,loc="lower center",bbox_to_anchor=(.5,0),frameon=False,fontsize=5.1,handlelength=.8,handletextpad=.3,columnspacing=.7,borderaxespad=0,labelspacing=.3)
        d.export(f"Fig4{chr(99+i)}_{donor}_domains",tab,[src],"Formal50k fit; full-resolution pixels, image-y inverted; numeric domains/colors are local to donor, not aligned cell types.",compressed=True)
    target_grid=grid[3,:2].subgridspec(2,1,height_ratios=[1,.14],hspace=.02)
    cx=fig.add_subplot(target_grid[0]);target_caption=fig.add_subplot(target_grid[1]);target_caption.set_axis_off()
    dx=fig.add_subplot(grid[3,2:]);label(cx,"k","Primary target expression\nFirst fixed-order reference donor");label(dx,"l","Corrected detection\nAll 8 reference donors")
    if not targets:
        no_measurement(cx,"No primary target nominated by the new freeze");no_measurement(dx,"No target identities to display")
        d.export("Fig4k_no_target",pd.DataFrame(columns=["reason"]),[],"Fresh freeze contains no nominated target.")
        d.export("Fig4l_no_target",pd.DataFrame(columns=["reason"]),[],"Fresh freeze contains no nominated target.")
    else:
        expr,spots,sources,why=target_expression(d,targets[0],DONORS[0])
        if expr is None:
            no_measurement(cx,f"{targets[0]['gene_symbol']} | {DONORS[0]}\nNo verified measurable raw-symbol row")
            record=spots.copy();record["measurement_status"]=why
        else:
            point=pixel_map(cx,expr,expr.plot_log1p_CP10000,cmap="magma",size=3)
            fig.colorbar(point,ax=cx,label="log1p(CP10000\ncorrected counts)",fraction=.035,pad=.015)
            target_caption.text(0,1,f"{targets[0]['gene_symbol']} | {DONORS[0]} | Full-resolution pixel geometry",transform=target_caption.transAxes,fontsize=5.3,va="top")
            record=expr
        d.export("Fig4k_primary_target_expression",record,sources,why+"; no candidate-driven donor selection; no invented micrometre calibration.",compressed=True)
        source="results/spatial/gene_detection_by_donor.tsv"
        detection=d.tab(source,["donor","gene_symbol","n_spots","n_positive_spots","detection_fraction","mean_log_expression"])
        values=np.full((len(targets),8),np.nan);rows=[]
        for iy,target in enumerate(targets):
            for ix,donor in enumerate(DONORS):
                hit=detection[detection.donor.eq(donor)&detection.gene_symbol.isin(target["spatial_symbols"])].copy()
                hit["identity_key"]=target["identity_key"];hit["review_canonical_symbol"]=target["gene_symbol"]
                if hit.empty:
                    hit=pd.DataFrame([{"donor":donor,"gene_symbol":"NA","identity_key":target["identity_key"],
                        "review_canonical_symbol":target["gene_symbol"],"plot_measurement_status":"NOT_OBSERVED",
                        "n_spots":np.nan,"n_positive_spots":np.nan,"detection_fraction":np.nan,
                        "mean_log_expression":np.nan,"reason":"No verified raw-symbol detection row"}])
                else:
                    hit["plot_measurement_status"]="OBSERVED_DESCRIPTIVE" if len(hit)==1 else "MULTIPLE_SOURCE_RECORDS"
                rows.append(hit)
                # Distinct raw aliases cannot be collapsed by averaging their
                # detection fractions, as they may be positive in the same spots.
                if len(hit)==1 and hit.plot_measurement_status.iloc[0]=="OBSERVED_DESCRIPTIVE":
                    values[iy,ix]=float(hit.detection_fraction.iloc[0])
                if len(hit)>1: dx.text(ix,iy,"multi",ha="center",va="center",fontsize=5)
        image=dx.imshow(np.ma.masked_invalid(values),vmin=0,vmax=1,cmap="Blues",aspect="auto")
        image.cmap.set_bad(COLORS["na"])
        dx.set_xticks(range(8),DONORS,rotation=45);dx.set_yticks(range(len(targets)),[t["gene_symbol"] for t in targets])
        for iy in range(len(targets)):
            for ix in range(8):
                if np.isfinite(values[iy,ix]):dx.text(ix,iy,f"{values[iy,ix]:.2f}",ha="center",va="center",fontsize=5.4,color="white" if values[iy,ix]>.55 else COLORS["dark"])
        fig.colorbar(image,ax=dx,label="Corrected expression>0\nfraction",fraction=.035,pad=.015)
        source_rows=pd.concat(rows,ignore_index=True) if rows else detection.iloc[:0]
        d.export("Fig4l_target_donor_detection",source_rows,[source],"Single verified raw-symbol detection is plotted; multiple raw alias records remain explicit multi/NA rather than an invalid average; missing is not0.")
    d.notes["Fig4"]={"biological_n":8,"subsampling_unit":"spots within donor","clinical_context":"Surgical proximal-small-intestine reference, not IBS cases or unscreened healthy volunteers", "domain_identity":"Numeric domains are donor-specific and not matched cell-type labels","display_coordinates":"Full-resolution pixel_col/pixel_row; equal aspect; inverted y; no histology or micrometre scale invented","expression_layer":"Author background-corrected original rows, summed after clipping, CP10000 log1p; not raw integer UMI"}
    return fig


def read_target_evidence(d,frozen,targets):
    path=d.register(frozen["evidence_file"])
    required=["identity_key","gene_symbol","evidence_axis","value","display","measurement_status","source_path","source_record_id"]
    ids={x["identity_key"] for x in targets};parts=[]
    for chunk in pd.read_csv(path,sep="\t",dtype=str,keep_default_na=False,chunksize=100000):
        need_columns(chunk,required,str(path))
        parts.append(chunk[chunk.identity_key.isin(ids)].copy())
    tab=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame(columns=required)
    if tab.duplicated(["identity_key","evidence_axis"]).any():raise ValueError("Evidence has duplicate identity/axis rows; never choose the most favorable record")
    unknown=set(tab.measurement_status)-set(EVIDENCE_STATUS)
    if unknown:raise ValueError(f"Unknown explicit measurement status: {unknown}")
    axes=frozen.get("figure_evidence_axes",[])
    if axes:
        expected={(identity,axis) for identity in ids for axis in axes}
        observed=set(zip(tab.identity_key,tab.evidence_axis))
        missing=expected-observed
        if missing:raise ValueError(f"Missing target/evidence-axis records, including explicit NA outcomes: {sorted(missing)[:10]}")
    for value in tab.source_path.unique():
        if value and value!="NA":
            # A source locator must be an explicit single path or a JSON list;
            # do not invent separators or silently resolve outside the project.
            refs=json.loads(value) if value.startswith("[") else [value]
            for ref in refs:d.register(ref)
    return tab


def evidence_summaries(evidence,targets):
    """Prospective descriptive10-axis presentation; no test, score or ranking."""
    rows=[]
    for target in targets:
        sub=evidence[evidence.identity_key.eq(target["identity_key"])].set_index("evidence_axis")
        def get(axis):
            if axis not in sub.index:raise ValueError(f"Missing explicit evidence axis {target['identity_key']}/{axis}")
            row=sub.loc[axis]
            value=pd.to_numeric(pd.Series([row["value"]]),errors="coerce").iloc[0]
            return value,row["measurement_status"],row["display"]
        def add(axis,display,status,source_axes):
            rows.append({"identity_key":target["identity_key"],"gene_symbol":target["gene_symbol"],"figure_axis":axis,
                         "display":display,"measurement_status":status,"source_axes":";".join(source_axes),
                         "summary_method":"descriptive_original_units_no_new_test_or_score"})
        def absence(source_axes):
            statuses={get(a)[1] for a in source_axes}
            if statuses=={"NOT_APPLICABLE"}:return "NOT_APPLICABLE"
            for state in ("MULTIPLE_SOURCE_RECORDS","NA_NATIVE_STATISTIC","UNTESTED","NOT_OBSERVED"):
                if state in statuses:return state
            return "NOT_OBSERVED"
        q,status,display=get("MAGMA_BH_FDR")
        add("MAGMA\nBH-FDR",f"{q:.2g}" if np.isfinite(q) else display,status,["MAGMA_BH_FDR"])
        for contrast in CONTRASTS:
            prefix="BULK_"+contrast
            l,ls,ld=get(prefix+"_log2FoldChange");q,qs,qd=get(prefix+"_padj")
            text=(f"{l:+.2f}" if np.isfinite(l) else "LFC=NA")+"\n"+(f"q={q:.2g}" if np.isfinite(q) else "q=NA")
            add(contrast.replace("_vs_HV","").replace("_","\n"),text,qs,[prefix+"_log2FoldChange",prefix+"_padj"])
        primary_axes=[f"SPATIAL_{x}_corrected_detection_fraction" for x in DONORS]
        primary=[get(a)[0] for a in primary_axes];finite=np.asarray(primary,dtype=float);finite=finite[np.isfinite(finite)]
        add("Detection\n8 donors",f"{finite.min():.2f}–{finite.max():.2f}\n{len(finite)}/8 measured" if len(finite) else absence(primary_axes).replace("_"," ").title(),
            "OBSERVED_DESCRIPTIVE" if len(finite) else absence(primary_axes),primary_axes)
        sensitivity_axes=[f"SPATIAL_{x}_QC_z3.5_corrected_detection_fraction" for x in DONORS]
        sensitivity=[get(a)[0] for a in sensitivity_axes];finite_s=np.asarray(sensitivity,dtype=float);finite_s=finite_s[np.isfinite(finite_s)]
        add("QC sensitivity\nUMI Z>−3.5\nmito Z<3.5",f"{finite_s.min():.2f}–{finite_s.max():.2f}\n{len(finite_s)}/8 measured" if len(finite_s) else absence(sensitivity_axes).replace("_"," ").title(),
            "OBSERVED_DESCRIPTIVE" if len(finite_s) else absence(sensitivity_axes),sensitivity_axes)
        tested_axes=[f"LR_{x}_global_pairs_tested_n" for x in DONORS]
        support_axes=[f"LR_{x}_global_pairs_BH_supported_n" for x in DONORS]
        tested=[get(a)[0] for a in tested_axes];support=[get(a)[0] for a in support_axes]
        measured=np.isfinite(tested)&np.isfinite(support)
        if np.any(measured):
            tt=np.asarray(tested)[measured];ss=np.asarray(support)[measured]
            donors_tested=int((tt>0).sum());donors_supported=int((ss>0).sum())
            add("LR donors\nsupported/tested",f"{donors_supported}/{donors_tested}\n{int(measured.sum())}/8 measured","OBSERVED_DESCRIPTIVE",tested_axes+support_axes)
            add("Global LR\ndonor-pairs",f"{int(ss.sum())}/{int(tt.sum())}\nsupported/tested","OBSERVED_DESCRIPTIVE",tested_axes+support_axes)
        else:
            missing_status=absence(tested_axes+support_axes)
            add("LR donors\nsupported/tested",missing_status.replace("_"," ").title(),missing_status,tested_axes+support_axes)
            add("Global LR\ndonor-pairs",missing_status.replace("_"," ").title(),missing_status,tested_axes+support_axes)
        local_t=[f"LR_{x}_local_pairs_tested_n" for x in DONORS];local_s=[f"LR_{x}_local_pairs_with_supported_spots_n" for x in DONORS]
        lt=np.asarray([get(a)[0] for a in local_t]);ls=np.asarray([get(a)[0] for a in local_s]);valid=np.isfinite(lt)&np.isfinite(ls)
        missing_status=absence(local_t+local_s)
        add("Local LR\ndonor-pairs",f"{int(ls[valid].sum())}/{int(lt[valid].sum())}\nsupported/tested" if valid.any() else missing_status.replace("_"," ").title(),
            "OBSERVED_DESCRIPTIVE" if valid.any() else missing_status,local_t+local_s)
    return pd.DataFrame(rows)


def representative_local(d):
    for donor in DONORS:
        selection=f"results/spatial/SpatialDM/{donor}/LOCAL_SELECTION_FROZEN.tsv"
        plan=d.tab(selection,["pair_id","pvalue","padj","selected_for_local"])
        plan=plan[bools(plan.selected_for_local).to_numpy()].sort_values(["padj","pvalue","pair_id"],kind="mergesort")
        if plan.empty:continue
        pair=str(plan.pair_id.iloc[0]);path=f"results/spatial/SpatialDM/{donor}/local_spot_statistics.tsv.gz"
        pieces=[]
        for chunk in pd.read_csv(d.register(path),sep="\t",chunksize=100000):
            need_columns(chunk,["pair_id","barcode","permutation_pvalue","permutation_padj","pixel_x","pixel_y"],path)
            pieces.append(chunk[chunk.pair_id.astype(str).eq(pair)])
        local=pd.concat(pieces,ignore_index=True)
        if local.empty or local.barcode.duplicated().any():raise ValueError("Saved selected local pair lacks a unique actual spot table")
        status=d.js(f"results/spatial/SpatialDM/{donor}/SPATIALDM_STATUS.json")
        if len(local)!=status["n_spots"]:raise ValueError("Representative local pair does not retain all actual tested spots")
        return donor,pair,local,[selection,path]
    return None,None,None,[]


def figure5(d,frozen,targets,modules,plt):
    from matplotlib.colors import ListedColormap,BoundaryNorm
    from matplotlib.patches import Rectangle,Patch
    import textwrap
    n_display_targets=max(1,len(targets))
    fig=plt.figure(figsize=(183/25.4,(170+12*(n_display_targets-1))/25.4),layout="constrained")
    grid=fig.add_gridspec(3,1,height_ratios=[3.3,.36*n_display_targets,.50],hspace=.10)
    top=grid[0].subgridspec(1,2,width_ratios=[3.3,1.6],wspace=.12)
    ax=fig.add_subplot(top[0])
    local_grid=top[1].subgridspec(3,1,height_ratios=[.30,1,.55],hspace=.04)
    local_head=fig.add_subplot(local_grid[0]);local_head.set_axis_off()
    local_head.text(0,1,"b",transform=local_head.transAxes,fontsize=8,fontweight="bold",va="top")
    local_head.text(.15,1,"Representative\nlocal test",transform=local_head.transAxes,fontsize=7,va="top")
    lx=fig.add_subplot(local_grid[1]);local_caption=fig.add_subplot(local_grid[2]);local_caption.set_axis_off()
    ex=fig.add_subplot(grid[1]);legend_ax=fig.add_subplot(grid[2]);legend_ax.set_axis_off()
    source="results/spatial/global_LR_all_donors.tsv"
    global_lr=d.tab(source,["donor","pair_id","pvalue","padj","ligand","receptor","complete_complex"])
    if set(global_lr.donor)!=set(DONORS) or global_lr.duplicated(["donor","pair_id"]).any():raise ValueError("Global LR donor/pair universe is incomplete or duplicated")
    if not bools(global_lr.complete_complex).all():raise ValueError("A tested incomplete complex is present")
    global_lr["BH_supported"]=numeric(global_lr.padj)<.05
    ranked=global_lr.groupby("pair_id",as_index=False).agg(n_tested_donors=("donor","nunique"),
        n_supported_donors=("BH_supported","sum"),median_fdr=("padj","median"))
    ranked=ranked.sort_values(["n_supported_donors","median_fdr","pair_id"],ascending=[False,True,True],kind="mergesort")
    top=ranked.head(12).pair_id.tolist();selected=global_lr[global_lr.pair_id.isin(top)].copy()
    pivot=selected.pivot(index="pair_id",columns="donor",values="padj").reindex(index=top,columns=DONORS)
    q=pivot.to_numpy(float);zeros=np.isfinite(q)&(q==0);valid=np.isfinite(q)&(q>=0)&(q<=1)
    if np.any(np.isfinite(q)&~valid):raise ValueError("Invalid canonical LR FDR")
    floor=1/1001.;display=np.full_like(q,np.nan);display[valid]=-np.log10(np.maximum(q[valid],floor))
    image=ax.imshow(np.ma.masked_invalid(display),vmin=0,vmax=-math.log10(floor),cmap="Blues",aspect="auto")
    image.cmap.set_bad(COLORS["na"])
    names=[]
    for pair in top:
        roles=global_lr[global_lr.pair_id.eq(pair)][["ligand","receptor"]].drop_duplicates()
        if len(roles)!=1:raise ValueError("The same pair ID has conflicting ligand/receptor identities")
        names.append(str(roles.ligand.iloc[0])+" | "+str(roles.receptor.iloc[0]))
    ax.set_xticks(range(8),DONORS,rotation=45);ax.set_yticks(range(len(top)),names)
    for iy in range(len(top)):
        for ix in range(8):
            if zeros[iy,ix]:ax.plot(ix,iy,marker="v",ms=3,color=COLORS["dark"],ls="")
            elif np.isfinite(q[iy,ix]) and q[iy,ix]<.05:ax.plot(ix,iy,marker="o",ms=1.7,color=COLORS["dark"],ls="")
            elif not np.isfinite(q[iy,ix]):ax.text(ix,iy,"×",ha="center",va="center",fontsize=6,color=COLORS["gray"])
    fig.colorbar(image,ax=ax,label="−log10 permutation BH-FDR\nZero censored at 1/1001",fraction=.035,pad=.02)
    label(ax,"a","Global LR: top 12 by donor support\nThen median FDR and pair ID")
    expanded=[]
    lookup=selected.set_index(["pair_id","donor"])
    for iy,pair in enumerate(top):
        for ix,donor in enumerate(DONORS):
            if (pair,donor) in lookup.index:r=lookup.loc[(pair,donor)].to_dict()
            else:r={"pair_id":pair,"donor":donor,"measurement_status":"UNTESTED_NOT_A_ZERO"}
            r.update(pair_id=pair,donor=donor,plot_minus_log10_FDR=display[iy,ix],empirical_zero_censored=bool(zeros[iy,ix]))
            expanded.append(r)
    d.export("Fig5a_global_LR_top12",pd.DataFrame(expanded),[source],"Exact top12 rule; missing donor/pair=untested; ▼ empirical zero null exceedance, not exact probability0; • positive nonzero FDR<0.05. No role-direction or binding inference.")
    d.export("Fig5a_all_global_pair_ranking",ranked,[source],"Display-only pair order, not gene/target ranking; full tested global universe retained.",compressed=True)
    donor,pair,local,refs=representative_local(d)
    if local is None:
        no_measurement(lx,"No local pair selected by the frozen global-significance rule")
        d.export("Fig5b_no_local_selection",pd.DataFrame(columns=["reason"]),[],"Explicit negative/no-local selection, no fabricated map.")
    else:
        local=local.copy();local["plot_BH_supported"]=numeric(local.permutation_padj)<.05
        lx.scatter(local.pixel_x,local.pixel_y,c=np.where(local.plot_BH_supported,COLORS["orange"],COLORS["light"]),s=3,linewidths=0,rasterized=True)
        lx.invert_yaxis();lx.set_aspect("equal");lx.set_xticks([]);lx.set_yticks([])
        roles=global_lr[global_lr.pair_id.astype(str).eq(pair)][["ligand","receptor"]].drop_duplicates()
        rolelabel=(str(roles.ligand.iloc[0])+" | "+str(roles.receptor.iloc[0])) if len(roles)==1 else "Pair "+pair
        local_caption_lines=textwrap.wrap(f"{donor}: {rolelabel}",width=23,break_long_words=False,break_on_hyphens=False)
        local_caption_lines += [f"{int(local.plot_BH_supported.sum())}/{len(local)} spots q<0.05","Permutation test","All tested spots shown"]
        local_caption.text(0,1,"\n".join(local_caption_lines),transform=local_caption.transAxes,fontsize=5.5,va="top",linespacing=1.2)
        d.export("Fig5b_representative_local",local,refs,"First fixed-order donor with local-selected pairs; lowest global FDR/P/pair ID; existing local BH family includes all selected pairs times spots within donor; spatial coexpression only.",compressed=True)
    label(ex,"c","Frozen-target evidence (original units)\nDescriptive counts are not gene-level P values")
    evidence=read_target_evidence(d,frozen,targets)
    d.export("Fig5c_targets_ALL_EVIDENCE",evidence,[frozen["evidence_file"]],"Every target record, including all donor/domain/role-specific negative, NA and untested evidence; no best-P collapse.",compressed=True)
    if not targets:
        no_measurement(ex,"No target nominated by the fresh Phase6 freeze:\n"+frozen["no_target_reason"])
        summary=pd.DataFrame(columns=["identity_key","gene_symbol","figure_axis","display","measurement_status"])
    else:
        summary=evidence_summaries(evidence,targets)
        axes=summary.figure_axis.drop_duplicates().tolist()
        if len(axes)!=10:raise ValueError("Prospective evidence summary must contain exactly10 axes")
        ex.set_xlim(-.5,len(axes)-.5);ex.set_ylim(len(targets)-.5,-.5)
        for iy,target in enumerate(targets):
            sub=summary[summary.identity_key.eq(target["identity_key"])].set_index("figure_axis")
            for ix,axis in enumerate(axes):
                r=sub.loc[axis];status=str(r.measurement_status);color=EVIDENCE_STATUS[status][0]
                hatch="//" if status in ("UNTESTED","NOT_OBSERVED","NOT_APPLICABLE","MULTIPLE_SOURCE_RECORDS","NA_NATIVE_STATISTIC") else None
                ex.add_patch(Rectangle((ix-.5,iy-.5),1,1,facecolor=color,edgecolor=COLORS["border"],lw=.45,hatch=hatch))
                display_text=str(r.display).replace("supported/tested","supported/\ntested")
                display_lines=[]
                for line in display_text.splitlines():
                    display_lines.extend(textwrap.wrap(line,width=13,break_long_words=False,break_on_hyphens=False) or [""])
                ex.text(ix,iy,"\n".join(display_lines),ha="center",va="center",fontsize=5.2,linespacing=1.05)
        axis_labels=[a.replace("QC sensitivity","QC\nsensitivity").replace("supported/tested","supported/\ntested") for a in axes]
        ex.set_xticks(range(len(axes)),axis_labels);ex.tick_params(axis="x",length=0,pad=5,labelsize=5.2)
        ex.set_yticks(range(len(targets)),[t["gene_symbol"]+("\n(primary)" if i==0 else "\n(backup)") for i,t in enumerate(targets)]);ex.tick_params(axis="y",length=0,labelsize=5.6,pad=5)
        used=summary.measurement_status.drop_duplicates().tolist()
        handles=[Patch(facecolor=EVIDENCE_STATUS[s][0],edgecolor=COLORS["border"],label=EVIDENCE_STATUS[s][1],hatch="//" if s in ("UNTESTED","NOT_OBSERVED","NOT_APPLICABLE","MULTIPLE_SOURCE_RECORDS","NA_NATIVE_STATISTIC") else None) for s in used]
        legend_ax.legend(handles=handles,loc="upper center",ncol=2,fontsize=5.5,frameon=False,handlelength=1.2,handletextpad=.5,columnspacing=1.2,borderaxespad=0)
    d.export("Fig5c_descriptive_10axis_summary",summary,[frozen["evidence_file"]],"Four LFC/q combinations + genetic q + all8donor detection/range + fixed QC sensitivity + donor and donor-pair LR counts; no new test, ranking or composite score.")
    d.notes["Fig5"]={"unit":"8 donor-level global tests; spots only within a donor's local test","global_selection":"Descending supported donors, ascending median canonical FDR, ascending pair ID; top12","local_selection":"First fixed-order donor with a local-selected pair, then ascending global FDR/P/pair ID","zero_policy":"Only visualization floors empirical0 at1/1001 and marks censoring; canonical P/FDR unchanged","biological_boundary":"Spatial coexpression does not prove directed signaling, direct molecular binding, pain mechanism or therapeutic effect","candidate_summary":"Original units displayed; statuses distinguish significance, descriptive observations, negatives and untested; no mixed-unit score"}
    return fig


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",required=True,type=Path)
    ap.add_argument("--freeze",default="results/integrated/CANDIDATE_FREEZE.json")
    ap.add_argument("--out",default="figures/main")
    ap.add_argument("--formats",default="pdf,png",help="pdf,png,svg,tiff; all drawing/export uses Python")
    ap.add_argument("--validate-only",action="store_true")
    ap.add_argument("--overwrite-render",action="store_true",help="Explicitly replace only this program's prior figure outputs; never analysis inputs")
    args=ap.parse_args()
    d=Inputs(args.root,args.out)
    d.register("data/metadata/FIGURE_CONTRACT.md")
    d.register(Path(__file__).resolve())
    frozen,targets,modules=validate_gates(d,args.freeze)
    # Long-table and all required fixed axes are checked before the first figure.
    evidence=read_target_evidence(d,frozen,targets)
    if targets:evidence_summaries(evidence,targets)
    formats=[x.strip().lower() for x in args.formats.split(",")]
    if not formats or set(formats)-{"pdf","png","svg","tif","tiff"}:raise ValueError("Unsupported explicit export format")
    if args.validate_only:
        print(json.dumps({"status":"INTERFACES_VALIDATED_NO_RENDER","freeze_sha256":sha(d.register(args.freeze)),
                          "fresh_target_identities":[x["identity_key"] for x in targets],"evidence_records":len(evidence)},ensure_ascii=False))
        return
    status_path=d.out/"FIGURE_RENDER_STATUS.json"
    if status_path.exists() and not args.overwrite_render:raise RuntimeError("Existing figure attempt/status retained; use a new --out or explicit --overwrite-render")
    plt,font=plotting_setup()
    begin=utc();d.out.mkdir(parents=True,exist_ok=True)
    write_json(status_path,{"status":"RENDERING","start_utc":begin,"backend":"python","model_fits_run":False})
    try:
        builders=[("Fig1",lambda:figure1(d,frozen,targets,modules,plt)),("Fig2",lambda:figure2(d,plt)),
                  ("Fig3",lambda:figure3(d,targets,plt)),("Fig4",lambda:figure4(d,targets,modules,plt)),
                  ("Fig5",lambda:figure5(d,frozen,targets,modules,plt))]
        for name,builder in builders:
            fig=builder();save_figure(fig,d,name,formats);plt.close(fig)
        d.finish()
        versions={p:importlib.metadata.version(p) for p in ("matplotlib","numpy","pandas","scipy","Pillow")}
        outputs=[p for p in sorted(d.out.rglob("*")) if p.is_file() and p!=status_path]
        write_json(status_path,{"status":"RENDERED_PENDING_VISUAL_QA","start_utc":begin,"finish_utc":utc(),
            "backend":"python","python":platform.python_version(),"package_versions":versions,"font":font,
            "formats":formats,"figures":["Fig1","Fig2","Fig3","Fig4","Fig5"],"model_fits_run":False,
            "candidate_selection_run":False,"freeze_path":str(d.register(args.freeze).relative_to(d.root)),
            "freeze_sha256":sha(d.register(args.freeze)),"output_sha256":{str(p.relative_to(d.root)):sha(p) for p in outputs},
            "required_human_or_visual_QA":"Inspect all PNGs at final size, overlap/crop/colorbars, and selectable PDF text before publication release"})
        print(json.dumps({"status":"RENDERED_PENDING_VISUAL_QA","out":str(d.out),"figures":5,"formats":formats}))
    except Exception as exc:
        plt.close("all")
        write_json(status_path,{"status":"FAILED","start_utc":begin,"finish_utc":utc(),"error":repr(exc),"backend":"python","model_fits_run":False})
        raise


if __name__=="__main__":
    main()
