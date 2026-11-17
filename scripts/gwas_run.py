#!/usr/bin/env python3
"""New MAGMA annotation and gene association; no reuse of previous model files.

Run this driver through run_logged.py too. Each external MAGMA command is in
turn run by the shared logger, which owns stdout/stderr and the execution trace.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys
import time


PROJECT_NAME = "IBS_A_OFFICIAL_RERUN_20261105"


def get_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--phase1-gate", required=True,
                   help="New-project spatial-preprocessing JSON gate; requires status PASS")
    p.add_argument("--bfile", default="data/reference/gwas/g1000_eur")
    p.add_argument("--gene-loc", default="data/reference/gwas/NCBI37.3.gene.loc")
    p.add_argument("--magma", default="env/tools/magma/magma")
    p.add_argument("--logger", default="scripts/run_logged.py")
    p.add_argument("--hgnc", help="Optional official HGNC complete-set TSV, within new root")
    p.add_argument("--threads", type=int, default=2, help="Recorded thread environment cap; MAGMA is not invoked with an unsupported thread flag")
    return p.parse_args()


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def confined(root, rel, exists=True):
    p = root / rel
    if p.is_symlink() or not p.resolve().is_relative_to(root):
        raise ValueError(f"Path is not a real new-root file: {p}")
    if exists and (not p.is_file() or not p.stat().st_size):
        raise FileNotFoundError(p)
    return p


def write_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def check_gate(root, rel):
    path = confined(root, rel)
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if obj.get("status") != "PASS":
        raise RuntimeError(f"Required gate is not PASS: {rel}")
    return path, obj


def bh_fdr(values):
    import numpy as np
    p = np.asarray(values, dtype=float)
    order = np.argsort(p, kind="mergesort")
    ordered_q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty(len(p), dtype=float)
    q[order] = np.minimum(ordered_q, 1.0)
    return q


def current_hgnc_by_entrez(path):
    """Preserve all many-to-one ambiguity; never infer a match from rank/P."""
    import pandas as pd
    h = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for name in ("hgnc_id", "symbol", "status", "entrez_id"):
        if name not in h:
            raise ValueError(f"HGNC complete set lacks {name}")
    result = {}
    for row in h.to_dict("records"):
        if row["status"] != "Approved" or not row["entrez_id"]:
            continue
        # The full HGNC table normally has one Entrez value. Explicitly retain
        # all values if the downloaded source contains a pipe-separated entry.
        for identifier in row["entrez_id"].split("|"):
            result.setdefault(identifier, []).append(row)
    return result


def main():
    args = get_args()
    root = args.root.resolve()
    if root.name != PROJECT_NAME or not root.is_dir():
        raise ValueError(f"Only the authorized new root named {PROJECT_NAME} is accepted")
    if not 1 <= args.threads <= 4:
        raise ValueError("Keep concurrent GWAS thread cap between one and four")
    import numpy as np
    import pandas as pd
    import os

    started, tick = utcnow(), time.monotonic()
    global_gate, _ = check_gate(root, "manifest/PHASE0_GATE.json")
    spatial_gate, _ = check_gate(root, args.phase1_gate)
    phase0_gate, phase0 = check_gate(root, "manifest/GWAS_PHASE0_GATE.json")
    qc_summary = json.loads(confined(root, "results/gwas/qc/summary.json").read_text())
    paths = {
        "magma": confined(root, args.magma), "logger": confined(root, args.logger),
        "gene_locations": confined(root, args.gene_loc),
        "p_values": confined(root, "results/gwas/qc/magma_snps_p_n.tsv"),
        "snp_locations": confined(root, "results/gwas/qc/magma_snps.loc"),
        **{"ld_" + ext: confined(root, args.bfile + "." + ext) for ext in ("bed", "bim", "fam")},
    }
    inputs = {name: {"path": str(path.relative_to(root)), "sha256": digest(path)} for name, path in paths.items()}
    reference_key = {"magma": "magma_executable", "gene_locations": "gene_locations", "ld_bed": "ld_bed", "ld_bim": "ld_bim", "ld_fam": "ld_fam"}
    for key, qckey in reference_key.items():
        if inputs[key]["sha256"] != qc_summary["reference"]["inputs"][qckey]["sha256"]:
            raise ValueError(f"Reference changed after Phase0: {key}")
    for key in ("p_values", "snp_locations"):
        rel = str(paths[key].relative_to(root)).replace("\\", "/")
        if inputs[key]["sha256"] != phase0["output_sha256"].get(rel):
            raise ValueError(f"Validated MAGMA input changed after Phase0: {key}")
    hgnc_map = None
    if args.hgnc:
        hgnc_path = confined(root, args.hgnc)
        inputs["hgnc"] = {"path": str(hgnc_path.relative_to(root)), "sha256": digest(hgnc_path)}
        hgnc_map = current_hgnc_by_entrez(hgnc_path)
    outdir = root / "results/gwas/magma"
    outdir.mkdir(parents=True, exist_ok=True)
    # A failed invocation leaves evidence and is not mistaken for a checkpoint.
    # A rerun requires an explicitly preserved/new attempt directory, not an
    # automatic fallback to an existing .genes.out file.
    attempt = outdir / "NEW_FIT_STARTED.json"
    if attempt.exists() or any(outdir.glob("*.genes.*")):
        raise FileExistsError("A MAGMA attempt already exists in this new project. Preserve/audit it explicitly before a further run; this driver never reuses a fitted model")
    write_json(attempt, {"start_iso8601": started, "inputs": inputs, "status": "STARTED", "is_new_fit": True})
    annotation = outdir / "IBS_GRCh37_gene_body"
    analysis = outdir / "IBS_GCST90016564"
    commands = []
    env = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[name] = str(args.threads)

    def logged_magma(module, command, input_files, output_files, is_fit, notes):
        run = [sys.executable, str(paths["logger"]), "--root", str(root), "--module", module]
        for path in input_files:
            run += ["--input", str(path.relative_to(root))]
        for path in output_files:
            run += ["--output", str(path.relative_to(root))]
        run += ["--software-version", "MAGMA v1.10; executable SHA256=" + inputs["magma"]["sha256"],
                "--is-new-fit", str(int(is_fit)), "--notes", notes, "--"] + command
        record = {"module": module, "source_command": shlex.join(command),
                  "logger_command": shlex.join(run), "start_iso8601": utcnow()}
        commands.append(record)
        write_json(outdir / "MAGMA_COMMANDS.json", commands)
        result = subprocess.run(run, cwd=root, env=env, check=False)
        record.update(finish_iso8601=utcnow(), logger_exit_code=result.returncode)
        write_json(outdir / "MAGMA_COMMANDS.json", commands)
        if result.returncode:
            raise RuntimeError(f"Shared logger / MAGMA failed: {module}, exit={result.returncode}; retain its stdout/stderr and trace")
        for output in output_files:
            if not output.is_file() or not output.stat().st_size:
                raise RuntimeError(f"MAGMA returned zero but expected output is empty/missing: {output}")

    annotate_command = [str(paths["magma"]), "--annotate", "window=0", "--snp-loc", str(paths["snp_locations"]),
                        "--gene-loc", str(paths["gene_locations"]), "--out", str(annotation)]
    logged_magma("GWAS_MAGMA_NEW_ANNOTATION", annotate_command,
                 [paths["magma"], paths["snp_locations"], paths["gene_locations"], global_gate, spatial_gate, phase0_gate],
                 [Path(str(annotation) + ".genes.annot"), Path(str(annotation) + ".log")], False,
                 "New SNP-to-gene annotation; GRCh37 gene body; zero kb flanks chosen before fit; exact validated LD matches only")
    association_command = [str(paths["magma"]), "--bfile", str(root / args.bfile), "synonyms=0",
                           "--gene-annot", str(annotation) + ".genes.annot", "--pval", str(paths["p_values"]),
                           "use=SNP,P", "ncol=N", "duplicate=error", "--gene-model", "snp-wise=mean", "--out", str(analysis)]
    logged_magma("GWAS_MAGMA_NEW_GENE_ASSOCIATION", association_command,
                 [paths["magma"], paths["p_values"], paths["ld_bed"], paths["ld_bim"], paths["ld_fam"], Path(str(annotation) + ".genes.annot")],
                 [Path(str(analysis) + ".genes.out"), Path(str(analysis) + ".genes.raw"), Path(str(analysis) + ".log")], True,
                 "New full gene association; per-SNP total N=N_CASE+N_CONTROL; no old results; snp-wise mean; synonyms suppressed after exact-ID/chr/bp QC")
    gene_raw = Path(str(analysis) + ".genes.out")
    gene_log = Path(str(analysis) + ".log").read_text(errors="replace")
    if "v1.10" not in gene_log:
        raise ValueError("Actual MAGMA log does not confirm expected version v1.10; inspect binary provenance before interpreting output")
    genes = pd.read_csv(gene_raw, sep=r"\s+", comment="#", dtype={"GENE": str})
    required = {"GENE", "CHR", "START", "STOP", "NSNPS", "P"}
    if not required.issubset(genes.columns) or genes.empty or genes.GENE.duplicated().any():
        raise ValueError(f"Invalid MAGMA gene output schema/uniqueness: {list(genes.columns)}")
    p = pd.to_numeric(genes.P, errors="coerce")
    if (~np.isfinite(p) | (p <= 0) | (p > 1)).any():
        raise ValueError("MAGMA gene output contains missing/zero/out-of-range P; no invented replacement")
    n_snps = pd.to_numeric(genes.NSNPS, errors="coerce")
    if (n_snps.isna() | (n_snps < 1)).any():
        raise ValueError("Tested gene without a positive SNP count")
    loc = pd.read_csv(paths["gene_locations"], sep="\t", header=None, dtype=str,
                      names=["gene_id", "reference_chr", "reference_start", "reference_stop", "reference_strand", "reference_symbol"])
    loc = loc.set_index("gene_id")
    output = genes.rename(columns={"GENE": "gene_id"}).copy()
    output["P"] = p
    output["BH_FDR"] = bh_fdr(p)
    output["reference_symbol"] = output.gene_id.map(loc.reference_symbol)
    if output.reference_symbol.isna().any():
        raise ValueError("MAGMA emitted a gene outside the validated official gene-location reference")
    mapped_rows = []
    for identifier, old_symbol in zip(output.gene_id, output.reference_symbol):
        choices = hgnc_map.get(identifier, []) if hgnc_map is not None else []
        unique = {r["hgnc_id"]: r for r in choices}
        if len(unique) == 1:
            row = next(iter(unique.values()))
            mapped_rows.append((row["symbol"], row["hgnc_id"], row.get("ensembl_gene_id", ""), "approved_entrez_unique"))
        elif len(unique) > 1:
            mapped_rows.append(("", "|".join(sorted(unique)), "", "ambiguous_entrez_to_hgnc"))
        else:
            mapped_rows.append(("", "", "", "no_approved_entrez_match" if hgnc_map is not None else "HGNC_reference_not_supplied"))
    output[["gene_symbol", "hgnc_id", "ensembl_gene_id", "HGNC_mapping_status"]] = mapped_rows
    output["significant_BH_0.05"] = output.BH_FDR < 0.05
    output["bonferroni_0.05"] = output.P < (0.05 / len(output))
    output["gene_id_numeric_order"] = pd.to_numeric(output.gene_id)
    output = output.sort_values(["P", "gene_id_numeric_order"], kind="mergesort").drop(columns="gene_id_numeric_order")
    output.insert(0, "rank", np.arange(1, len(output) + 1))
    all_file = root / "results/gwas/magma_all_genes.tsv"
    fdr_file = root / "results/gwas/magma_BH_FDR_lt_0.05.tsv"
    output.to_csv(all_file, sep="\t", index=False, float_format="%.17g")
    output.loc[output["significant_BH_0.05"]].to_csv(fdr_file, sep="\t", index=False, float_format="%.17g")
    summary = {"status": "PASS", "start_iso8601": started, "finish_iso8601": utcnow(),
               "driver_duration_seconds": time.monotonic() - tick, "is_new_fit": True,
               "genome_build": "GRCh37", "gene_model": "snp-wise=mean", "annotation_window_kb": [0, 0],
               "snp_synonyms": "disabled; Phase0 exact SNP ID plus chr/bp intersection", "N_definition": "per-SNP N_CASE+N_CONTROL",
               "observed_N_range_in_valid_raw": qc_summary["numeric_ranges"]["N_TOTAL"],
               "input_eligible_snps": phase0["eligible_snps_for_magma"],
               "genes_tested": len(output), "genes_BH_FDR_lt_0.05": int(output["significant_BH_0.05"].sum()),
               "genes_bonferroni_0.05": int(output["bonferroni_0.05"].sum()),
               "BH_universe": "every valid tested gene, independent of HGNC mapping or downstream expression",
               "HGNC_mapping_status_counts": output.HGNC_mapping_status.value_counts().to_dict(),
               "actual_commands": commands, "input_hashes": inputs,
               "top_genes_by_statistical_P_only": output.head(20).replace({np.nan: None}).to_dict("records"),
               "interpretation": "Gene-level association, not causal-gene identification or independent-locus count; no final candidate selected by this module"}
    summary_file = root / "results/gwas/MAGMA_SUMMARY.json"
    write_json(summary_file, summary)
    model_outputs = [all_file, fdr_file, summary_file, gene_raw, Path(str(analysis) + ".genes.raw"),
                     Path(str(analysis) + ".log"), Path(str(annotation) + ".genes.annot"), Path(str(annotation) + ".log")]
    with (root / "manifest/GWAS_MAGMA_OUTPUT_MANIFEST.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["relative_path", "bytes", "sha256"])
        for path in model_outputs:
            writer.writerow([str(path.relative_to(root)), path.stat().st_size, digest(path)])
    write_json(root / "manifest/GWAS_MAGMA_GATE.json",
               {"status": "PASS", "is_new_fit": True, "finish_iso8601": utcnow(),
                "genes_tested": len(output), "BH_significant": int(output["significant_BH_0.05"].sum()),
                "model_process_exit_codes": [c["logger_exit_code"] for c in commands],
                "output_sha256": {str(p.relative_to(root)): digest(p) for p in model_outputs}})
    print(json.dumps({k: v for k, v in summary.items() if k not in ("input_hashes", "top_genes_by_statistical_P_only", "actual_commands")}, indent=2), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"GWAS_MAGMA_FAIL: {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        raise
