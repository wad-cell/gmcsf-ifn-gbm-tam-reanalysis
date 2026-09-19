# -*- coding: utf-8 -*-
"""B7 derived object: GSE182109 myeloid barcode-aligned expression matrix.

Replaces the former ad-hoc workspace script (C4.4F marked this object
"NOT BUNDLED / recipe only"). It deterministically rebuilds

    temp/B7/gse182109_myeloid_expr_aligned.pkl

from two documented inputs:

    in1 = temp/B7/gse182109_myeloid_expr.pkl      (cells x genes, produced by
                                                    scripts/b7_01_prep_data.py)
    in2 = output/B6/GSE182109_myeloid_annot.csv   (barcode/sample/patient/lineage
                                                    table, produced by
                                                    scripts/b6_m8_analysis_v2.py)

Consumers: scripts/b7_03b_classical_ifn.py, scripts/b7_06_rescue.py,
scripts/prep_patient_paired.py (all read the aligned object from temp/B7/).

FROZEN ALIGNMENT RULES (do not change; they reproduce the original run object)
  R1 Row order baseline: the row order of the annotation table column `barcode`
     is the single baseline. No re-sorting, no alphabetical ordering, no
     set-based ordering is applied.
  R2 Barcode intersection: only barcodes present in BOTH the expression matrix
     and the annotation table survive.
  R3 Duplicate handling: when a barcode occurs k times in the expression matrix,
     its rows are assigned -- in ascending original row order -- to the 1st, 2nd,
     ... k-th occurrence of that barcode in the annotation table
     (first-come-first-served consumption). Extra expression rows that are never
     requested are dropped; they do not raise an error (R2/R4 only require every
     requested barcode to be satisfiable).
  R4 Result shape/order: the aligned object must have exactly the same number of
     rows as the annotation table, and its index must equal the annotation
     `barcode` sequence element-wise.
  R5 Failure conditions -> hard error (SystemExit, non-zero, no file written):
     (a) either input file is missing;
     (b) in1 is not a DataFrame, or has no columns;
     (c) in2 has no `barcode` column, or contains null / empty barcode values;
     (d) some requested barcode has fewer available rows in in1 than requested
         (i.e. the intersection under R2/R3 cannot fulfil the annotation order);
     (e) R4 shape or order verification fails.
     Partial or "best effort" output is never produced.

Usage
-----
  python scripts/generate_gse182109_aligned_object.py                 # write to
        <PROJ_ROOT>/temp/B7/gse182109_myeloid_expr_aligned.pkl
  python scripts/generate_gse182109_aligned_object.py --check-only    # rebuild
        in memory, verify against --verify-sha256, write nothing
  python scripts/generate_gse182109_aligned_object.py --out <path> \
        --verify-sha256 <hex> --expect-rows 59118

The script prints a JSON block with input sizes, sha256 of the written file,
matrix checksum and shape, so runs can be audited or compared with the frozen
object's sha256.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

PROJ_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

EXPR_PKL = os.path.join("temp", "B7", "gse182109_myeloid_expr.pkl")
ANNOT_CSV = os.path.join("output", "B6", "GSE182109_myeloid_annot.csv")
DEFAULT_OUT = os.path.join("temp", "B7", "gse182109_myeloid_expr_aligned.pkl")


def fail(msg: str) -> "NoReturn":  # noqa: F821
    print("[FAIL] " + msg, file=sys.stderr)
    sys.exit(2)


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def matrix_checksum(df: pd.DataFrame) -> str:
    """Deterministic content checksum: index + columns + values (row-major)."""
    h = hashlib.sha256()
    h.update("\n".join(map(str, df.index)).encode("utf-8"))
    h.update(b"\x00")
    h.update("\n".join(map(str, df.columns)).encode("utf-8"))
    h.update(b"\x00")
    vals = np.ascontiguousarray(df.to_numpy())
    h.update(np.asarray(vals.shape, dtype="<i8").tobytes())
    h.update(vals.tobytes())
    return h.hexdigest()


def build(expr_path: str, annot_path: str, expect_rows=None):
    if not os.path.isfile(expr_path):
        fail("expression matrix not found: %s" % expr_path)
    if not os.path.isfile(annot_path):
        fail("annotation table not found: %s" % annot_path)

    expr = pd.read_pickle(expr_path)
    annot = pd.read_csv(annot_path)

    if not isinstance(expr, pd.DataFrame):
        fail("expression object is %s, expected pandas.DataFrame" % type(expr).__name__)
    if expr.shape[1] == 0:
        fail("expression matrix has 0 gene columns")
    if "barcode" not in annot.columns:
        fail("annotation table has no 'barcode' column (columns=%s)" % list(annot.columns))

    bc_annot = annot["barcode"]
    if bc_annot.isna().any():
        fail("annotation table contains %d null barcode value(s)" % int(bc_annot.isna().sum()))
    if (bc_annot.astype(str).str.strip() == "").any():
        fail("annotation table contains empty-string barcode value(s)")
    if expect_rows is not None and len(annot) != int(expect_rows):
        fail("annotation row count %d != expected %s" % (len(annot), expect_rows))

    # R1/R2/R3: positional consumption of expression rows in annotation order.
    positions = {}
    for i, bc in enumerate(expr.index):
        positions.setdefault(bc, []).append(i)
    used = {}
    idx = []
    for bc in bc_annot:
        bucket = positions.get(bc)
        k = used.get(bc, 0)
        if bucket is None or k >= len(bucket):
            have = 0 if bucket is None else len(bucket)
            fail(
                "barcode '%s' requested %d time(s) but the expression matrix "
                "provides %d row(s) (R2/R3 unsatisfiable)" % (bc, k + 1, have)
            )
        idx.append(bucket[k])
        used[bc] = k + 1

    aligned = expr.iloc[np.asarray(idx, dtype="int64")]

    # R4: shape and element-wise order verification.
    if len(aligned) != len(annot):
        fail("aligned rows %d != annotation rows %d" % (len(aligned), len(annot)))
    if not np.array_equal(aligned.index.to_numpy(), bc_annot.to_numpy()):
        fail("aligned index does not match the annotation barcode order (R4)")
    if aligned.index.has_duplicates:
        print("[WARN] aligned index contains duplicate barcodes (annotation table does)")

    return expr, annot, aligned, sorted(used.keys())


def main() -> int:
    ap = argparse.ArgumentParser(description="Rebuild the GSE182109 barcode-aligned matrix.")
    ap.add_argument("--project-root", default=PROJ_ROOT, help="repository root (default: package parent)")
    ap.add_argument("--out", default=None, help="output .pkl (default: <root>/temp/B7/gse182109_myeloid_expr_aligned.pkl)")
    ap.add_argument("--check-only", action="store_true", help="rebuild and verify, write nothing")
    ap.add_argument("--verify-sha256", default=None, help="hex sha256 the written file must match")
    ap.add_argument("--expect-rows", type=int, default=None, help="expected annotation/aligned row count")
    args = ap.parse_args()

    root = os.path.abspath(args.project_root)
    expr_path = os.path.join(root, EXPR_PKL)
    annot_path = os.path.join(root, ANNOT_CSV)
    out_path = args.out or os.path.join(root, DEFAULT_OUT)

    expr, annot, aligned, matched = build(expr_path, annot_path, args.expect_rows)

    report = {
        "expr_input": os.path.relpath(expr_path, root),
        "expr_input_sha256": sha256_file(expr_path),
        "expr_shape": list(expr.shape),
        "annot_input": os.path.relpath(annot_path, root),
        "annot_input_sha256": sha256_file(annot_path),
        "annot_shape": list(annot.shape),
        "distinct_barcodes_matched": len(matched),
        "aligned_shape": list(aligned.shape),
        "aligned_index_matches_annot_order": True,
        "aligned_matrix_sha256": matrix_checksum(aligned),
        "output": None if args.check_only else os.path.abspath(out_path),
        "output_sha256": None,
        "check_only": bool(args.check_only),
    }

    if not args.check_only:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        aligned.to_pickle(out_path)
        report["output_sha256"] = sha256_file(out_path)

    if args.verify_sha256:
        got = report["output_sha256"]
        if got is None:
            # check-only mode: compare via the matrix checksum of the referenced object
            if os.path.isfile(args.verify_sha256):
                got = sha256_file(args.verify_sha256)
            else:
                fail("--verify-sha256 without --out requires an existing file path")
        if got.lower() != args.verify_sha256.lower():
            fail("sha256 mismatch: got %s expected %s" % (got, args.verify_sha256))
        report["sha256_verified"] = True

    print("[OK] aligned object rebuilt")
    print("```json")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print("```")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
