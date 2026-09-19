# -*- coding: utf-8 -*-
"""GSE309037 / GSE309038 / GSE309039 raw-count input chain (added in C4.4G).

The three series are author-processed count matrices published as GEO series
supplementary files. No decompression or format conversion step exists: each
file is a gzip-compressed, SEMICOLON-separated CSV (header row present, first
column = gene identifier, index_col=0) that the DESeq2 scripts read directly
with pandas:

    pd.read_csv(<file>, compression="gzip", sep=";", index_col=0)

This tool makes that chain auditable and repeatable:
  --list                 print the frozen accession -> file -> URL -> sha256 table
  --check-local DIR      verify sha256 of files found under DIR (offline)
  --download DIR         download each file from the GEO FTP/HTTPS mirror, then verify
  --parse DIR            parse each file and print shape / first columns (offline)
  --project-root DIR     root that holds temp/GSE30903x/raw (default: package parent)

Nothing is downloaded unless --download is given; nothing outside DIR is written.

FROZEN TABLE (do not edit; sha256 computed from the original-run downloads)
  accession  geo file name                    bytes      sha256
  GSE309037  GSE309037_raw_counts.csv.gz       508232     598c95060c49de1e0e5667cdcf2ecb9561fbc9fe08d7101bffb84595e8173bb0
  GSE309038  GSE309038_raw_counts.csv.gz       451756     b943820976ab1954fa58235be94b7ee3505d89b9da0966320d3e1ac0cc9e4103
  GSE309039  GSE309039_data1_raw_counts.csv.gz 1005871    b94cffe93670ff1bafd742139005763aae705f27d1c7336d3f938c92d782929e
  GSE309039  GSE309039_data2_raw_counts.csv.gz 801308     08aae847c8e98e14fc6dacb9a0c51f40d0830c6b4580f7f2740589c60ff6c5ef
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.request

PROJ_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/{acc}/suppl/{name}"
CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"

# accession, geo file name, local relative path used by the DESeq2 scripts, byte size, sha256
FROZEN_INPUTS = [
    ("GSE309037", "GSE309037_raw_counts.csv.gz",
     os.path.join("temp", "GSE309037", "raw", "GSE309037_raw_counts.csv.gz"),
     508232, "598c95060c49de1e0e5667cdcf2ecb9561fbc9fe08d7101bffb84595e8173bb0"),
    ("GSE309038", "GSE309038_raw_counts.csv.gz",
     os.path.join("temp", "GSE309038", "raw", "GSE309038_raw_counts.csv.gz"),
     451756, "b943820976ab1954fa58235be94b7ee3505d89b9da0966320d3e1ac0cc9e4103"),
    ("GSE309039", "GSE309039_data1_raw_counts.csv.gz",
     os.path.join("temp", "GSE309039", "raw", "GSE309039_data1_raw_counts.csv.gz"),
     1005871, "b94cffe93670ff1bafd742139005763aae705f27d1c7336d3f938c92d782929e"),
    ("GSE309039", "GSE309039_data2_raw_counts.csv.gz",
     os.path.join("temp", "GSE309039", "raw", "GSE309039_data2_raw_counts.csv.gz"),
     801308, "08aae847c8e98e14fc6dacb9a0c51f40d0830c6b4580f7f2740589c60ff6c5ef"),
]

CONSUMERS = {
    "GSE309039": "scripts/02_b1_deseq2.py (stage B1, both data1 and data2 matrices)",
    "GSE309037": "scripts/08_b2_deseq2.py (stage B2)",
    "GSE309038": "scripts/12_b3_deseq2.py (stage B3)",
}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def table() -> list:
    rows = []
    for acc, name, rel, size, digest in FROZEN_INPUTS:
        rows.append({
            "accession": acc,
            "geo_supplementary_file": name,
            "download_url": FTP_BASE.format(acc=acc, name=name),
            "series_page": CGI.format(acc=acc),
            "delimiter": ";",
            "compression": "gzip",
            "local_relative_path": rel,
            "bytes": size,
            "sha256": digest,
            "required_by": CONSUMERS[acc],
            "conversion_step": "none (author-supplied count matrix; pandas reads gzip+; directly)",
        })
    return rows


def cmd_list() -> int:
    print(json.dumps(table(), indent=2, ensure_ascii=False))
    return 0


def check_dir(directory: str, parse: bool) -> int:
    ok = True
    results = []
    for row in table():
        cand = os.path.join(directory, os.path.basename(row["local_relative_path"]))
        if not os.path.isfile(cand):
            results.append({"file": row["geo_supplementary_file"], "status": "MISSING", "path_searched": cand})
            ok = False
            continue
        got = sha256_file(cand)
        status = "OK" if got == row["sha256"] else "SHA256_MISMATCH"
        entry = {"file": row["geo_supplementary_file"], "status": status, "sha256": got}
        if status != "OK":
            ok = False
        if parse and status == "OK":
            import pandas as pd
            df = pd.read_csv(cand, compression="gzip", sep=";", index_col=0, nrows=5)
            entry["parse_ok"] = True
            entry["n_columns_first5rows"] = int(df.shape[1])
            entry["first_columns"] = [str(c) for c in list(df.columns)[:3]]
        results.append(entry)
    print(json.dumps(results, indent=2, ensure_ascii=False))
    return 0 if ok else 2


def check_project(root: str, parse: bool) -> int:
    """Verify the four files at their canonical temp/GSE30903x/raw/ locations."""
    ok = True
    results = []
    for row in table():
        cand = os.path.join(root, row["local_relative_path"])
        if not os.path.isfile(cand):
            results.append({"file": row["geo_supplementary_file"], "status": "MISSING", "path_searched": cand})
            ok = False
            continue
        got = sha256_file(cand)
        status = "OK" if got == row["sha256"] else "SHA256_MISMATCH"
        entry = {"file": row["geo_supplementary_file"], "status": status, "sha256": got, "path": cand}
        if status != "OK":
            ok = False
        if parse and status == "OK":
            import pandas as pd
            df = pd.read_csv(cand, compression="gzip", sep=";", index_col=0, nrows=5)
            entry["parse_ok"] = True
            entry["n_columns_first5rows"] = int(df.shape[1])
            entry["first_columns"] = [str(c) for c in list(df.columns)[:3]]
        results.append(entry)
    print(json.dumps(results, indent=2, ensure_ascii=False))
    return 0 if ok else 2


def cmd_download(directory: str) -> int:
    os.makedirs(directory, exist_ok=True)
    for row in table():
        dest = os.path.join(directory, row["geo_supplementary_file"])
        print("[GET] %s -> %s" % (row["download_url"], dest))
        urllib.request.urlretrieve(row["download_url"], dest)
        got = sha256_file(dest)
        if got != row["sha256"]:
            print("[FAIL] sha256 mismatch for %s: %s" % (dest, got), file=sys.stderr)
            return 2
        print("[OK] sha256 verified: %s" % got)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--list", action="store_true", help="print the frozen input table")
    ap.add_argument("--check-local", metavar="DIR", help="verify sha256 of the four files inside DIR")
    ap.add_argument("--check-project", action="store_true",
                    help="verify the four files at their canonical temp/GSE30903x/raw/ paths under --project-root")
    ap.add_argument("--parse", action="store_true", help="with --check-local: also parse the first rows")
    ap.add_argument("--download", metavar="DIR", help="download from GEO and verify (network access required)")
    ap.add_argument("--project-root", default=PROJ_ROOT, help="repository root holding temp/GSE30903x/raw")
    args = ap.parse_args()

    if args.download:
        return cmd_download(args.download)
    if args.check_project:
        return check_project(os.path.abspath(args.project_root), args.parse)
    if args.check_local:
        return check_dir(args.check_local, args.parse)
    if not args.list:
        # default: report the expected local locations derived from --project-root
        rows = table()
        for r in rows:
            r["expected_local_path"] = os.path.join(os.path.abspath(args.project_root), r["local_relative_path"])
        print(json.dumps(rows, indent=2, ensure_ascii=False))
        return 0
    return cmd_list()


if __name__ == "__main__":
    raise SystemExit(main())
