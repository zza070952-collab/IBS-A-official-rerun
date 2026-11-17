#!/usr/bin/env python3
"""Fresh, candidate-blind GSE166869 raw-input audit. No fitted data are read.

This script deliberately does NOT infer subject IDs from library numbering,
GSM ordering, sex, BMI, or SRA sample aliases. It does not select a 60-sample
cohort from group quotas. It can be executed before the global Phase-0 gate.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import decimal
import gzip
import hashlib
import io
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

PROJECT_NAME = "IBS_A_OFFICIAL_RERUN_20261105"
SOURCE_URLS = {
    "merged": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166869/suppl/GSE166869_merged_gene_count_rpkm_by_gene.tsv.gz",
    "soft": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166869/soft/GSE166869_family.soft.gz",
    "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166869",
    "paper": "https://eprints.whiterose.ac.uk/id/eprint/194614/7/Manuscript%20R2.pdf",
    "publisher": "https://onlinelibrary.wiley.com/doi/10.1111/apt.16591",
}
PUBLISHED_GROUP_COUNTS = {"Duo:HV": 11, "Duo:IBSC": 10, "Duo:IBSD": 12,
                          "Jej:HV": 12, "Jej:IBSC": 7, "Jej:IBSD": 8}
RRNA_TYPES = {"rrna", "mt_rrna", "rrna_pseudogene"}
SUBJECT_KEYS = {"subject_id", "subject id", "patient_id", "patient id", "donor_id",
                "donor id", "participant_id", "participant id", "individual_id", "individual id"}
BATCH_KEYS = {"batch", "sequencing batch", "sequencing_batch", "library batch",
              "library_batch", "sequencing lane", "sequencing_lane"}


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def jwrite(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def twrite(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        w.writeheader()
        for row in rows:
            w.writerow({k: ("NA" if row.get(k) is None else row.get(k, "")) for k in fields})


def text_open(path):
    return gzip.open(path, "rt", encoding="utf-8-sig", newline="") if str(path).endswith(".gz") else open(path, encoding="utf-8-sig", newline="")


def safe_project_input(root, path):
    # Check the logical input path; immutable raw symlinks may legitimately
    # point to an independently verified archive managed by the parent runner.
    p = Path(os.path.abspath(path))
    if os.path.commonpath([str(root), str(p)]) != str(root):
        raise ValueError(f"Input must be explicitly staged beneath the NEW project: {p}")
    if not p.is_file():
        raise FileNotFoundError(p)
    return p


def canonical_key(value):
    return re.sub(r"[^a-z0-9]", "", value.lower())


def first_column(header, aliases, required=True):
    normalized = {canonical_key(c): c for c in header}
    hits = [normalized[a] for a in aliases if a in normalized]
    hits = list(dict.fromkeys(hits))
    if len(hits) != 1:
        if not required and not hits:
            return None
        raise ValueError(f"Need one explicit annotation column for {aliases}; found {hits}; header={header}")
    return hits[0]


def normalize_group(s):
    s = re.sub(r"[^A-Z]", "", str(s).upper())
    return {"HV": "HV", "HEALTHYVOLUNTEER": "HV", "HEALTHYCONTROL": "HV",
            "CONTROL": "HV", "IBSC": "IBSC", "IBSD": "IBSD"}.get(s)


def normalize_tissue(s):
    s = re.sub(r"[^A-Z]", "", str(s).upper())
    return {"DUO": "Duo", "DUODENUM": "Duo", "DUODENAL": "Duo",
            "JEJ": "Jej", "JEJUNUM": "Jej", "JEJUNAL": "Jej"}.get(s)


def normalize_sex(s):
    s = str(s).strip().lower()
    return {"female": "F", "f": "F", "male": "M", "m": "M"}.get(s)


def parse_soft(path):
    samples, series, current = [], collections.defaultdict(list), None
    with text_open(path) as f:
        for line_number, raw in enumerate(f, 1):
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                current = {"gsm": line.split(" = ", 1)[1], "soft_line": line_number,
                           "fields": collections.defaultdict(list)}
                samples.append(current)
            elif line.startswith("^"):
                current = None
            elif " = " in line and line.startswith("!Sample_") and current is not None:
                key, value = line.split(" = ", 1)
                current["fields"][key.removeprefix("!Sample_")].append(value)
            elif " = " in line and line.startswith("!Series_"):
                key, value = line.split(" = ", 1)
                series[key.removeprefix("!Series_")].append(value)
    if len(samples) != 69 or len({x["gsm"] for x in samples}) != 69:
        raise ValueError(f"Expected 69 unique GSM records in the independently supplied SOFT; found {len(samples)}")
    result = []
    for sample in samples:
        fields = sample["fields"]
        title = ";".join(fields.get("title", []))
        libraries = re.findall(r"(?<![A-Za-z0-9])(\d+_S\d+)(?![A-Za-z0-9])", title)
        if len(libraries) != 1:
            raise ValueError(f"Title must explicitly identify one library: {sample['gsm']} {title}")
        chars = collections.defaultdict(list)
        for value in fields.get("characteristics_ch1", []):
            if ":" in value:
                k, v = value.split(":", 1)
                chars[k.strip().lower()].append(v.strip())
        tg = re.search(r"\b(Duo|Jej)[.\s_-]+(HV|IBS[-_ ]?[CD])\b", title, re.I)
        tissue_candidates = {normalize_tissue(tg.group(1))} if tg else set()
        group_candidates = {normalize_group(tg.group(2))} if tg else set()
        for k, values in chars.items():
            if k in {"tissue", "tissue type", "organ", "source tissue"}:
                tissue_candidates.update(normalize_tissue(v) for v in values)
            if k in {"group", "disease", "disease status", "diagnosis", "condition", "phenotype"}:
                group_candidates.update(normalize_group(v) for v in values)
        tissue_candidates.discard(None)
        group_candidates.discard(None)
        if len(tissue_candidates) != 1 or len(group_candidates) != 1:
            raise ValueError(f"Unknown/conflicting tissue/group for {sample['gsm']}: {title} {dict(chars)}")
        sex_values = {normalize_sex(v) for k, vs in chars.items() if k in {"sex", "gender"} for v in vs}
        sex_values.discard(None)
        if len(sex_values) > 1:
            raise ValueError(f"Conflicting explicit sex fields: {sample['gsm']}")
        subject_values = {v for k, vs in chars.items() if k in SUBJECT_KEYS for v in vs if v.strip()}
        if len(subject_values) > 1:
            raise ValueError(f"Conflicting explicit subject IDs: {sample['gsm']}")
        batch = {k: vs for k, vs in chars.items() if k in BATCH_KEYS}
        result.append({"gsm": sample["gsm"], "library_id": libraries[0], "title": title,
                       "tissue": next(iter(tissue_candidates)), "group": next(iter(group_candidates)),
                       "sex": next(iter(sex_values), "NA"), "subject_id": next(iter(subject_values), "NA"),
                       "subject_id_basis": "EXPLICIT_GEO_CHARACTERISTIC" if subject_values else "UNAVAILABLE_NOT_INFERRED",
                       "batch": json.dumps(batch, ensure_ascii=False, sort_keys=True) if batch else "NA",
                       "characteristics_json": json.dumps(dict(chars), ensure_ascii=False, sort_keys=True),
                       "relations_json": json.dumps(fields.get("relation", []), ensure_ascii=False),
                       "soft_line": sample["soft_line"]})
    if len({s["library_id"] for s in result}) != len(result):
        raise ValueError("GEO titles map multiple GSM records to the same count-library ID")
    return result, dict(series)


def download(url, dest):
    dest = Path(dest)
    event = {"url": url, "path": str(dest), "start_utc": now()}
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "IBS-clean-reanalysis-source-audit/1.0"})
        with urllib.request.urlopen(req, timeout=25) as response:
            event["http_status"] = response.status
            data = response.read(12 * 1024 * 1024 + 1)
            event["final_url"] = response.url
            event["content_type"] = response.headers.get("Content-Type", "")
        if len(data) > 12 * 1024 * 1024:
            raise ValueError("Source document exceeds the audit download size limit")
        if dest.suffix.lower() == ".pdf" and not data.startswith(b"%PDF"):
            raise ValueError("Response is not a PDF")
        if dest.suffix.lower() in {"docx", ".docx", ".xlsx"} and not data.startswith(b"PK"):
            raise ValueError("Response is not an Office ZIP container")
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
        event.update(status="DOWNLOADED", bytes=len(data), sha256=sha256(dest))
    except Exception as exc:
        event.update(status="ACCESS_FAILED", error=f"{type(exc).__name__}: {exc}")
    event["finish_utc"] = now()
    return event


def extract_document_identifiers(path):
    """Record possible source identifiers for human verification, not cohort selection."""
    result = {"path": str(path), "text_extraction": "UNAVAILABLE", "gsm_hits": [], "library_hits": []}
    text = ""
    try:
        if path.suffix == ".docx":
            with zipfile.ZipFile(path) as z:
                text = " ".join(ET.fromstring(z.read("word/document.xml")).itertext())
        elif path.suffix == ".pdf":
            try:
                from pypdf import PdfReader
                text = "\n".join(p.extract_text() or "" for p in PdfReader(path).pages)
            except ImportError:
                result["text_extraction"] = "PYPDF_NOT_INSTALLED_MANUAL_REVIEW_REQUIRED"
        elif path.suffix == ".xlsx":
            # XLSX shared strings and inline XML retain IDs without needing an
            # Excel engine. This is only an identifier search, never a table parse.
            with zipfile.ZipFile(path) as z:
                text = " ".join(z.read(n).decode("utf-8") for n in z.namelist() if n.endswith(".xml"))
        elif path.suffix in {".xml", ".csv", ".txt", ".json"}:
            text = path.read_text(encoding="utf-8", errors="replace")
        if text:
            result.update(text_extraction="SEARCHED", gsm_hits=sorted(set(re.findall(r"GSM\d+", text))),
                          library_hits=sorted(set(re.findall(r"(?<![A-Za-z0-9])\d+_S\d+(?![A-Za-z0-9])", text))))
    except Exception as exc:
        result["text_extraction"] = f"FAILED: {exc}"
    return result


def recover_sources(meta, soft_samples, recover):
    folder = meta / "original_sources"
    folder.mkdir(parents=True, exist_ok=True)
    events = []
    if recover:
        events.append(download(SOURCE_URLS["paper"], folder / "Grover_2021_accepted.pdf"))
        names = ["apt16591-sup-0002-MethodS1.docx", "apt16591-sup-0003-TableS1.pdf",
                 "apt16591-sup-0004-TableS2.pdf", "apt16591-sup-0005-TableS3.xlsx",
                 "apt16591-sup-0006-TableS4.xlsx", "apt16591-sup-0007-TableS5.xlsx",
                 "apt16591-sup-0008-TableS6.xlsx"]
        for name in names:
            url = "https://onlinelibrary.wiley.com/action/downloadSupplement?" + urllib.parse.urlencode({"doi": "10.1111/apt.16591", "file": name})
            event = download(url, folder / name)
            events.append(event)
            if event["status"] != "DOWNLOADED":
                direct = f"https://onlinelibrary.wiley.com/doi/suppl/10.1111/apt.16591/supinfo/{name}"
                events.append(download(direct, folder / name))
        search = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=GSE166869&retmax=200&retmode=json"
        event = download(search, folder / "SRA_esearch.json")
        events.append(event)
        if event["status"] == "DOWNLOADED":
            try:
                ids = json.loads((folder / "SRA_esearch.json").read_text())["esearchresult"]["idlist"]
                if ids:
                    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=" + ",".join(ids)
                    events.append(download(base + "&rettype=runinfo&retmode=text", folder / "SRA_RunInfo.csv"))
                    events.append(download(base + "&retmode=xml", folder / "SRA_full_metadata.xml"))
            except Exception as exc:
                events.append({"status": "SRA_PARSE_FAILED", "error": str(exc), "recorded_utc": now()})
    jwrite(meta / "SOURCE_ACCESS_ATTEMPTS.json", {"created_utc": now(), "recovery_requested": recover, "attempts": events})
    scans = [extract_document_identifiers(p) for p in sorted(folder.iterdir()) if p.is_file()]
    jwrite(meta / "SOURCE_IDENTIFIER_SEARCH.json", {"sources": scans,
           "interpretation": "Identifier hits are NOT an official retained60 list and never determine sample inclusion."})
    # Only exact explicit subject attributes in SRA SAMPLE records are eligible.
    sra_subjects = {}
    xml_path = folder / "SRA_full_metadata.xml"
    if xml_path.exists():
        try:
            tree = ET.parse(xml_path)
            for sample in tree.findall(".//SAMPLE"):
                explicit = set()
                for a in sample.findall(".//SAMPLE_ATTRIBUTE"):
                    key, value = (a.findtext("TAG") or "").strip().lower(), (a.findtext("VALUE") or "").strip()
                    if key in SUBJECT_KEYS and value:
                        explicit.add(value)
                gsms = set(re.findall(r"GSM\d+", ET.tostring(sample, encoding="unicode")))
                if len(explicit) == 1 and len(gsms) == 1:
                    gsm, subject = next(iter(gsms)), next(iter(explicit))
                    if gsm in sra_subjects and sra_subjects[gsm] != subject:
                        raise ValueError(f"Conflicting SRA subject IDs for {gsm}")
                    sra_subjects[gsm] = subject
        except Exception as exc:
            raise ValueError(f"Downloaded SRA subject metadata need review: {exc}") from exc
    for sample in soft_samples:
        subject = sra_subjects.get(sample["gsm"])
        if subject:
            if sample["subject_id"] != "NA" and sample["subject_id"] != subject:
                raise ValueError(f"GEO/SRA explicit subject conflict: {sample['gsm']}")
            sample["subject_id"] = subject
            sample["subject_id_basis"] = "EXPLICIT_SRA_SAMPLE_ATTRIBUTE"
    return events, scans


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--soft", required=True)
    ap.add_argument("--recover-sources", action="store_true")
    ap.add_argument("--official60", help="Optional manually authenticated TSV with exactly 60 GSM rows; never infer by quota")
    ap.add_argument("--official60-source", help="Authoritative provenance URI/document locator for supplied official list")
    args = ap.parse_args()
    root = Path(os.path.abspath(args.root))
    if root.name != PROJECT_NAME:
        raise ValueError(f"This program only writes the authorized new project {PROJECT_NAME}")
    meta, processed = root / "data/metadata/bulk", root / "data/processed/bulk"
    meta.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)
    status_path = meta / "BULK_PHASE0_STATUS.json"
    if status_path.exists() and json.loads(status_path.read_text())["status"] == "PASS":
        raise RuntimeError("Phase0 already PASS in this NEW run. Audit existing checkpoint rather than overwriting it.")
    start = now()
    merged, soft = safe_project_input(root, args.merged), safe_project_input(root, args.soft)
    inputs = {str(merged): sha256(merged), str(soft): sha256(soft)}
    # Freeze thresholds BEFORE opening expression values, extracting sample
    # library totals, fitting any model, or reading any downstream target.
    protocol = {
        "schema_version": 1, "frozen_utc": now(), "input_sha256": inputs,
        "selection_basis": "OFFICIAL60_IF_AUTHENTICATED_ELSE_INDEPENDENT_CANDIDATE_BLIND_QC",
        "main_min_non_rrna_counts": 1000000, "main_min_detected_non_rrna_genes": 5000,
        "sensitivity_min_non_rrna_counts": {"qc_low_750k": 750000, "qc_high_1250k": 1250000},
        "sensitivity_other_rules_unchanged": True,
        "rrna_gene_types_casefold_exact": sorted(RRNA_TYPES),
        "valid_mapping": "one count library exact title ID to exactly one GSM, explicit consistent tissue/group",
        "unmapped_count_columns": "QUARANTINE_NO_INFERRED_DISEASE_TISSUE_OR_PATIENT",
        "low_expression_filter": "per tissue, count>=10 in at least the smallest retained group size",
        "model": "~sex+group if every sex known, both sexes represented and full rank; else ~group; HV reference",
        "batch": "record available technical batch; no fabricated batch; unadjusted batch is a disclosed limitation",
        "subject_policy": "only explicit GEO/SRA subject identifiers; no S-number/GSM/phenotype inference; tissues fitted separately",
        "univariate_outlier_policy": "PCA distance flag at median+4MAD; report only; no automatic sample removal",
        "deseq2": {"test": "Wald", "fitType": "parametric", "sfType": "ratio", "betaPrior": False,
                    "minReplicatesForReplace": 7, "results_alpha": 0.05, "BH": True,
                    "cooksCutoff": "DESeq2 default", "independentFiltering": True, "lfcThreshold": 0},
        "gsea": {"method": "fgseaMultilevel", "rank": "finite Wald statistics with finite DESeq2 pvalue, no FDR filter",
                 "duplicate_symbol": "highest baseMean then lexical original gene_id; no P-based selection",
                 "collections": ["Hallmark", "Reactome"], "minSize": 15, "maxSize": 500, "eps": 0,
                 "nPermSimple": 1000, "nproc": 1, "seed_base": 20261105,
                 "FDR": "Primary BH pooled over both tested collections within contrast; within-collection BH is supplementary"},
        "prohibited_selection_inputs": "No candidate names, prior project models/results/rankings, group quotas, or DE significance",
        "method_boundary": "DESeq2 independent reanalysis is not exact reproduction of the source paper glmLRT workflow",
    }
    protocol_path = meta / "QC_PROTOCOL_FROZEN.json"
    if protocol_path.exists():
        old = json.loads(protocol_path.read_text())
        for key in ["frozen_utc", "input_sha256"]:
            old.pop(key, None)
        compare = dict(protocol)
        for key in ["frozen_utc", "input_sha256"]:
            compare.pop(key, None)
        if old != compare:
            raise RuntimeError("A different QC protocol already exists; append an explicit amendment instead of overwriting")
        protocol = json.loads(protocol_path.read_text())
        if protocol["input_sha256"] != inputs:
            raise RuntimeError("Frozen Phase0 inputs changed")
    else:
        jwrite(protocol_path, protocol)
    try:
        samples, series = parse_soft(soft)
        events, scans = recover_sources(meta, samples, args.recover_sources)
        jwrite(meta / "GEO_SERIES_RAW_FIELDS.json", series)
        twrite(meta / "GEO_samples_all.tsv", samples)
        by_library = {s["library_id"]: s for s in samples}
        official = None
        if args.official60:
            if not args.official60_source:
                raise ValueError("An official60 list requires explicit authoritative provenance")
            official_path = safe_project_input(root, args.official60)
            with official_path.open(encoding="utf-8-sig", newline="") as f:
                official_rows = list(csv.DictReader(f, delimiter="\t"))
            official = {r["gsm"] for r in official_rows}
            if len(official_rows) != 60 or len(official) != 60 or not official.issubset({s["gsm"] for s in samples}):
                raise ValueError("Authenticated official60 file must list 60 unique deposited GSM records")
            actual = collections.Counter(f"{s['tissue']}:{s['group']}" for s in samples if s["gsm"] in official)
            if dict(actual) != PUBLISHED_GROUP_COUNTS:
                raise ValueError("Authenticated official60 list conflicts with published tissue/group counts")
        recovery = {"recorded_utc": now(), "status": "AUTHENTICATED_60_GSM_LIST_SUPPLIED" if official else "NOT_RECOVERED_FROM_CURRENT_VERIFIABLE_SOURCES",
                    "official60_source": args.official60_source, "official60_file": args.official60,
                    "published_counts": PUBLISHED_GROUP_COUNTS, "source_attempts_file": str(meta / "SOURCE_ACCESS_ATTEMPTS.json"),
                    "identifier_search_file": str(meta / "SOURCE_IDENTIFIER_SEARCH.json"),
                    "limits": "No list is inferred from counts, numbering, group quotas or candidate expression. Failed downloads remain failed; unparsed files require manual review.",
                    "exact_extra3_identity": "UNVERIFIABLE_WITHOUT_OFFICIAL60_INDIVIDUAL_IDS" if official is None else "COMPUTABLE_BY_EXACT_ID_SET_COMPARISON"}
        jwrite(meta / "OFFICIAL60_RECOVERY_STATUS.json", recovery)
        with text_open(merged) as f:
            reader = csv.DictReader(f, delimiter="\t")
            header = reader.fieldnames
            if not header or len(set(header)) != len(header):
                raise ValueError("Missing or duplicate raw merged header fields")
            gene_column = first_column(header, ["geneid", "ensemblgeneid", "ensemblid", "gene"])
            symbol_column = first_column(header, ["genesymbol", "genename", "symbol", "name"])
            type_column = first_column(header, ["genetype", "genebiotype", "featuretype", "biotype"])
            count_columns = [c for c in header if re.fullmatch(r"\d+_S\d+_counts?", c, re.I)]
            rpkm_columns = [c for c in header if re.fullmatch(r"\d+_S\d+_rpkm", c, re.I)]
            if not count_columns or len(count_columns) != len(rpkm_columns):
                raise ValueError("Raw schema needs matching explicit library_count and library_rpkm columns")
            library_ids = [re.sub(r"_counts?$", "", c, flags=re.I) for c in count_columns]
            if len(set(library_ids)) != len(library_ids):
                raise ValueError("Duplicate raw count library IDs")
            if {re.sub(r"_rpkm$", "", c, flags=re.I) for c in rpkm_columns} != set(library_ids):
                raise ValueError("Count/RPKM sample columns do not match exactly")
            missing = set(by_library) - set(library_ids)
            if missing:
                raise ValueError(f"Deposited GSM libraries missing from count table: {sorted(missing)}")
            column_rows = []
            for i, c in enumerate(header, 1):
                kind = "RAW_COUNT" if c in count_columns else "RPKM_NOT_FOR_DESEQ2" if c in rpkm_columns else "FEATURE_ANNOTATION"
                lib = re.sub(r"_(?:counts?|rpkm)$", "", c, flags=re.I) if kind != "FEATURE_ANNOTATION" else "NA"
                column_rows.append({"raw_column_number_1based": i, "raw_column_name": c, "kind": kind,
                                    "library_id": lib, "gsm": by_library.get(lib, {}).get("gsm", "NA")})
            twrite(meta / "COUNT_COLUMNS.tsv", column_rows)
            jwrite(meta / "RAW_SCHEMA.json", {"header": header, "gene_id_column": gene_column,
                   "gene_symbol_column": symbol_column, "gene_type_column": type_column, "n_counts": len(count_columns),
                   "n_rpkm": len(rpkm_columns), "only_count_columns_used_for_DESeq2": True})
            sums = [0] * len(library_ids)
            detected = [0] * len(library_ids)
            full_sums = [0] * len(library_ids)
            full_detected = [0] * len(library_ids)
            seen, n_rows, n_rrna, min_value, max_value = set(), 0, 0, None, 0
            type_counts = collections.Counter()
            genes = []
            out_path = processed / "raw_non_rrna_counts.tsv.gz"
            out_tmp = out_path.with_name(out_path.name + ".tmp")
            with gzip.open(out_tmp, "wt", encoding="utf-8", newline="") as of:
                writer = csv.writer(of, delimiter="\t", lineterminator="\n")
                writer.writerow(["gene_id", "gene_symbol", "gene_type"] + library_ids)
                for row in reader:
                    if None in row:
                        raise ValueError(f"Inconsistent column count at gene row {n_rows + 1}")
                    gene = row[gene_column]
                    if not re.fullmatch(r"ENSG\d+(?:\.\d+)?", gene) or gene in seen:
                        raise ValueError(f"Non-Ensembl or duplicate feature ID: {gene}")
                    seen.add(gene)
                    gene_type = row[type_column]
                    rrna = gene_type.lower() in RRNA_TYPES
                    type_counts[gene_type] += 1
                    n_rows += 1
                    n_rrna += int(rrna)
                    values = []
                    for i, c in enumerate(count_columns):
                        text = row[c]
                        try:
                            if re.fullmatch(r"\d+", text):
                                value = int(text)
                            else:
                                exact = decimal.Decimal(text)
                                if not exact.is_finite() or exact < 0 or exact != exact.to_integral_value():
                                    raise ValueError("not a nonnegative integer")
                                value = int(exact)
                        except Exception as exc:
                            raise ValueError(f"Invalid raw count {gene}/{c}: {text}") from exc
                        if value < 0 or value > 2147483647:
                            raise ValueError(f"Raw count exceeds DESeq2 integer storage: {gene}/{c}={value}")
                        values.append(value)
                        full_sums[i] += value
                        full_detected[i] += int(value > 0)
                        min_value = value if min_value is None else min(min_value, value)
                        max_value = max(max_value, value)
                        if not rrna:
                            sums[i] += value
                            detected[i] += int(value > 0)
                    genes.append({"gene_id": gene, "gene_symbol": row[symbol_column], "gene_type": gene_type, "rrna_excluded": rrna})
                    if not rrna:
                        writer.writerow([gene, row[symbol_column], gene_type] + values)
            os.replace(out_tmp, out_path)
        twrite(processed / "gene_annotation.tsv", genes)
        manifest = []
        for i, lib in enumerate(library_ids):
            s = by_library.get(lib)
            row = dict(s) if s else {"gsm": "NA", "library_id": lib, "title": "NA", "tissue": "NA", "group": "NA", "sex": "NA",
                    "subject_id": "NA", "subject_id_basis": "UNMAPPED_LIBRARY_NOT_INFERRED", "batch": "NA",
                    "characteristics_json": "{}", "relations_json": "[]", "soft_line": "NA"}
            eligible = s is not None
            independent_pass = eligible and sums[i] >= 1000000 and detected[i] >= 5000
            row.update(mapping_status="UNIQUE_EXACT_GEO_TITLE" if eligible else "QUARANTINED_NO_GSM",
                       raw_count_column=count_columns[i], raw_total_counts=full_sums[i], raw_detected_genes=full_detected[i],
                       non_rrna_counts=sums[i], detected_non_rrna_genes=detected[i], independent_QC_pass=independent_pass,
                       main_analysis_pass=(row["gsm"] in official) if official is not None else independent_pass,
                       qc_low_750k=eligible and sums[i] >= 750000 and detected[i] >= 5000,
                       qc_high_1250k=eligible and sums[i] >= 1250000 and detected[i] >= 5000,
                       main_selection_basis="AUTHENTICATED_OFFICIAL60" if official is not None else "FROZEN_INDEPENDENT_QC",
                       original60_membership=(row["gsm"] in official) if official is not None else "UNVERIFIABLE",
                       exclusion_reason="UNMAPPED_LIBRARY" if not eligible else "PASS" if independent_pass else
                         ";".join((["NON_RRNA_COUNTS_BELOW_1000000"] if sums[i] < 1000000 else []) +
                                  (["DETECTED_NON_RRNA_GENES_BELOW_5000"] if detected[i] < 5000 else [])))
            manifest.append(row)
        twrite(meta / "SAMPLE_MANIFEST.tsv", manifest)
        group_rows = []
        for tissue in ["Duo", "Jej"]:
            for group in ["HV", "IBSC", "IBSD"]:
                subset = [r for r in manifest if r["tissue"] == tissue and r["group"] == group]
                record = {"tissue": tissue, "group": group, "n_GEO_deposited": len(subset),
                          "published_group_count_not_identity": PUBLISHED_GROUP_COUNTS[f"{tissue}:{group}"]}
                for selection in ["main_analysis_pass", "independent_QC_pass", "qc_low_750k", "qc_high_1250k"]:
                    record[selection] = sum(r[selection] for r in subset)
                group_rows.append(record)
        twrite(meta / "GROUP_COUNTS_AND_SENSITIVITY.tsv", group_rows)
        # Duplicate verified subject IDs within the SAME tissue invalidate the
        # simple independent-sample model and require an explicit design decision.
        verified = collections.Counter((r["tissue"], r["subject_id"]) for r in manifest
                     if r["main_analysis_pass"] and r["subject_id"] != "NA")
        duplicates = [list(k) for k, v in verified.items() if v > 1]
        if duplicates:
            raise ValueError(f"Verified repeated subject within tissue: {duplicates}; stop for appropriate model, do not silently treat as independent")
        if any(r["main_analysis_pass"] < 2 for r in group_rows):
            raise ValueError("A main tissue/group has fewer than two libraries")
        for label in ["qc_low_750k", "qc_high_1250k"]:
            if any(r[label] < 2 for r in group_rows):
                raise ValueError(f"Frozen sensitivity {label} has <2 libraries in a tissue/group")
        validation = {"source": str(merged), "sha256": inputs[str(merged)], "n_gene_rows": n_rows,
                      "n_count_columns": len(count_columns), "n_rpkm_columns": len(rpkm_columns),
                      "count_entries_checked": n_rows * len(count_columns), "all_raw_count_values_nonnegative_integer": True,
                      "min_raw_count": min_value, "max_raw_count": max_value, "excluded_rrna_gene_rows": n_rrna,
                      "retained_non_rrna_gene_rows": n_rows - n_rrna, "n_GSM": len(samples),
                      "n_exact_mapped_libraries": len(samples), "n_quarantined_unmapped_libraries": len(manifest) - len(samples),
                      "n_main_libraries": sum(r["main_analysis_pass"] for r in manifest),
                      "n_explicit_subject_ids": sum(r["subject_id"] != "NA" for r in samples),
                      "n_missing_sex_mapped_libraries": sum(r["sex"] == "NA" for r in samples),
                      "gene_types": dict(type_counts), "sample_independence": "Tissues modeled separately; no unverified subject pairing inferred",
                      "gzip_integrity": "PASS: complete streaming read reached EOF and checked CRC"}
        jwrite(meta / "RAW_COUNT_VALIDATION.json", validation)
        twrite(meta / "DATA_PROVENANCE.tsv", [
            {"resource": "GSE166869_merged", "staged_path": str(merged), "official_source_url": SOURCE_URLS["merged"], "sha256": inputs[str(merged)]},
            {"resource": "GSE166869_SOFT", "staged_path": str(soft), "official_source_url": SOURCE_URLS["soft"], "sha256": inputs[str(soft)]}])
        twrite(meta / "DATA_VALIDATION.tsv", [
            {"check": "raw_count_integer_audit", "status": "PASS", "value": n_rows * len(count_columns)},
            {"check": "unique_GSM_mapping", "status": "PASS", "value": len(samples)},
            {"check": "GSM_unmapped_columns_quarantined", "status": "PASS", "value": len(manifest) - len(samples)},
            {"check": "official60_individual_ID_list", "status": "AVAILABLE" if official else "NOT_RECOVERED_LIMITATION", "value": len(official) if official else "NA"},
            {"check": "subject_pairing", "status": "EXPLICIT_ONLY" if any(r["subject_id"] != "NA" for r in samples) else "UNAVAILABLE_NOT_INFERRED", "value": sum(r["subject_id"] != "NA" for r in samples)}])
        outputs = [p for p in sorted(meta.rglob("*")) if p.is_file() and p != status_path]
        outputs += [p for p in sorted(processed.iterdir()) if p.is_file()]
        jwrite(status_path, {"status": "PASS", "module": "bulk_phase0", "start_utc": start, "finish_utc": now(),
               "input_sha256": inputs, "protocol_sha256": sha256(protocol_path),
               "manifest_sha256": sha256(meta / "SAMPLE_MANIFEST.tsv"),
               "matrix_sha256": sha256(out_path), "annotation_sha256": sha256(processed / "gene_annotation.tsv"),
               "n_main_libraries": validation["n_main_libraries"], "official60_status": recovery["status"],
               "interpretation": "Independent candidate-blind reanalysis if original60 IDs unavailable; no claim of exact source-paper cohort/method reproduction",
               "new_downstream_models_fitted": False,
               "output_sha256": {str(p): sha256(p) for p in outputs}})
        print(json.dumps({"status": "PASS", "validation": validation, "status_file": str(status_path)}, ensure_ascii=False))
    except Exception as exc:
        jwrite(status_path, {"status": "FAILED", "module": "bulk_phase0", "start_utc": start, "finish_utc": now(),
                            "input_sha256": inputs, "error": f"{type(exc).__name__}: {exc}"})
        raise


if __name__ == "__main__":
    main()
