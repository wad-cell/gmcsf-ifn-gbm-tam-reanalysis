# -*- coding: utf-8 -*-
"""
B4d: 外部 GBM scRNA (GSM4972211, GSE163120) 验证 ISG 抑制签名
- 用 B4b 的 68 基因 ISG 核心集做签名投影
- 计算每细胞 ISG 签名得分，按细胞类型/样本比较
- 输入: aria2 下载的完整 matrix.csv.gz + annot.csv.gz
- 输出到 output/B4/
"""
import gzip, os, sys
import numpy as np
import pandas as pd

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

TEMP = os.path.join(PROJ_ROOT, "temp", "GSM4972211")
OUT = os.path.join(PROJ_ROOT, "output", "B4")
CORE = os.path.join(PROJ_ROOT, "output", "B4", "B4b_ISG_core_set.csv")
MATRIX = os.path.join(TEMP, "matrix.csv.gz")
EXPECT = 41389236  # 完整 gz 文件字节数

def load_matrix():
    print("Loading matrix...", flush=True)
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as f:
        df = pd.read_csv(f, index_col=0)
    print(f"Matrix shape: {df.shape}", flush=True)
    return df

def main():
    if not os.path.exists(MATRIX):
        print("matrix.csv.gz not found", flush=True)
        sys.exit(2)
    size = os.path.getsize(MATRIX)
    print(f"Matrix size: {size} (expect {EXPECT})", flush=True)
    if size < EXPECT:
        print("INCOMPLETE - still downloading", flush=True)
        sys.exit(2)

    # 加载 annot
    annot = pd.read_csv(os.path.join(TEMP, "annot.csv.gz"), compression="gzip")
    print(f"Annot: {annot.shape}", flush=True)

    # 加载矩阵
    mat = load_matrix()

    # 加载 ISG 核心集
    core = pd.read_csv(CORE)
    print(f"ISG core genes: {len(core)}", flush=True)

    # 检查矩阵行名格式
    row0 = mat.index[:5].tolist()
    print(f"Row names sample: {row0}", flush=True)
    is_ensembl = any(str(r).startswith("ENSG") for r in row0)

    # 匹配基因
    if is_ensembl:
        # 行名是 ensembl，去掉版本号匹配
        mat_idx = {str(r).split(".")[0]: r for r in mat.index}
        core["ens"] = core["gene"].str.split(".").str[0]
        found_ens = [e for e in core["ens"] if e in mat_idx]
        found = [mat_idx[e] for e in found_ens]
        print(f"ISG matched (ensembl): {len(found)}/{len(core)}", flush=True)
    else:
        # 行名是 symbol
        found = [g for g in core["symbol"] if g in mat.index]
        print(f"ISG matched (symbol): {len(found)}/{len(core)}", flush=True)

    if len(found) < 10:
        print("Too few matched - check gene naming", flush=True)
        sys.exit(3)

    sub = mat.loc[found]
    maxv = sub.values.max()
    print(f"Max value in ISG submatrix: {maxv}", flush=True)
    if maxv > 50:
        sub = np.log1p(sub)
        print("Applied log1p (raw counts)", flush=True)
    else:
        print("Data appears already log-transformed", flush=True)

    score = sub.mean(axis=0)
    score.name = "ISG_score"

    # 对齐 annot
    cells = score.index.tolist()
    annot_cells = set(annot["cell"].tolist())
    common = [c for c in cells if c in annot_cells]
    print(f"Cells matched: {len(common)}/{len(cells)}", flush=True)
    if len(common) < 1000:
        print("Cell barcode mismatch - check naming", flush=True)
        sys.exit(4)

    s = score.loc[common]
    a = annot.set_index("cell").loc[common]

    res = pd.DataFrame({"ISG_score": s.values}, index=common)
    res["cluster"] = a["cluster"].values
    res["sample"] = a["sample"].values
    res["ident"] = a["ident"].values

    # 按细胞类型汇总
    summ = res.groupby("cluster")["ISG_score"].agg(["mean", "median", "std", "count"])
    summ = summ.sort_values("mean", ascending=False)
    print("\nISG score by cell type:")
    print(summ.round(4), flush=True)

    # 按样本汇总
    summ_s = res.groupby("sample")["ISG_score"].agg(["mean", "median", "std", "count"])
    summ_s = summ_s.sort_values("mean", ascending=False)
    print("\nISG score by sample:")
    print(summ_s.round(4), flush=True)

    # 保存
    os.makedirs(OUT, exist_ok=True)
    res.to_csv(os.path.join(OUT, "B4d_cell_ISG_scores.csv"))
    summ.to_csv(os.path.join(OUT, "B4d_ISG_score_by_celltype.csv"))
    summ_s.to_csv(os.path.join(OUT, "B4d_ISG_score_by_sample.csv"))
    pd.DataFrame({"gene": found}).to_csv(os.path.join(OUT, "B4d_ISG_genes_matched.csv"), index=False)
    print("\nSaved to output/B4/", flush=True)

if __name__ == "__main__":
    main()
