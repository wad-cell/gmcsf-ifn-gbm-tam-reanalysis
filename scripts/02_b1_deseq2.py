# -*- coding: utf-8 -*-
"""
B1: GSE309039 主分析 —— 配对 DESeq2（data1 + data2），C1-C12 对比
模型: ~ donor + condition (4 水平)，数值 contrast 向量提取主效应/交互/简单效应
可续跑: 每步结果落盘，已存在则跳过
输出: output/B1/
"""
import os, sys, gzip, json, time
import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

RAW = os.path.join(PROJ_ROOT, "temp", "GSE309039", "raw")
OUT = os.path.join(PROJ_ROOT, "output", "B1")
os.makedirs(OUT, exist_ok=True)

N_CPUS = 4

def load_counts(fname):
    df = pd.read_csv(os.path.join(RAW, fname), sep=";", index_col=0)
    df = df[~df.index.duplicated(keep="first")]
    df = df.astype(int)
    return df

def build_meta(cols):
    meta = pd.DataFrame(index=cols)
    meta["condition"] = [c.rsplit("_n", 1)[0] for c in cols]
    meta["donor"] = ["donor_" + c.split("_n")[1] for c in cols]
    return meta

def run_block(block, counts, meta, contrasts, outdir):
    """拟合 DESeq2 并提取所有对比"""
    dds_path = os.path.join(outdir, f"{block}_dds.pkl")
    if os.path.exists(dds_path):
        print(f"[{block}] 加载已拟合 dds")
        import pickle
        with open(dds_path, "rb") as f:
            dds = pickle.load(f)
    else:
        print(f"[{block}] 拟合 DESeq2: {counts.shape[0]} genes x {counts.shape[1]} samples")
        t0 = time.time()
        dds = DeseqDataSet(counts=counts.T, metadata=meta,
                           design_factors=["donor", "condition"],
                           refit_cooks=True, n_cpus=N_CPUS)
        dds.deseq2()
        print(f"[{block}] 拟合完成, 耗时 {time.time()-t0:.1f}s")
        import pickle
        with open(dds_path, "wb") as f:
            pickle.dump(dds, f)

    # 设计矩阵列名（formulaic treatment coding）
    dm_cols = list(dds.obsm["design_matrix"].columns)
    print(f"[{block}] 设计矩阵列: {dm_cols}")
    colmap = {c: i for i, c in enumerate(dm_cols)}

    # 保存设计矩阵列名供 contrast 构建
    with open(os.path.join(outdir, f"{block}_design_cols.json"), "w") as f:
        json.dump(dm_cols, f)

    # 提取每个对比
    for cid, level_coefs in contrasts.items():
        out_csv = os.path.join(outdir, f"{block}_{cid}.csv")
        if os.path.exists(out_csv):
            print(f"[{block}] {cid} 已存在, 跳过")
            continue
        # 数值 contrast 向量（长度=设计矩阵列数）；参考水平无列，系数自动为 0
        v = np.zeros(len(dm_cols))
        for level, coef in level_coefs.items():
            col = f"condition[T.{level}]"
            if col in colmap:
                v[colmap[col]] = coef
            else:
                print(f"[{block}] {cid}: 水平 {level} 为参考水平(系数0), 忽略")
        print(f"[{block}] 计算 {cid}: {v}")
        t0 = time.time()
        stats = DeseqStats(dds, contrast=v, n_cpus=N_CPUS, quiet=True)
        stats.summary()
        res = stats.results_df.copy()
        res.insert(0, "gene", res.index)
        res.to_csv(out_csv, index=False)
        print(f"[{block}] {cid} 完成: {res.shape[0]} genes, {time.time()-t0:.1f}s")

# ============ data1 ============
print("===== data1 =====")
d1 = load_counts("GSE309039_data1_raw_counts.csv.gz")
meta1 = build_meta(d1.columns)
meta1.to_csv(os.path.join(OUT, "data1_metadata.csv"), index=False)

# 数值 contrast（{水平: 系数}，参考水平自动忽略）
contrasts1 = {
    "C1_culture_main": {"CoCul_IgG": 0.5, "CoCul_aGMCSF": 0.5, "MoCul_IgG": -0.5, "MoCul_aGMCSF": -0.5},
    "C2_antibody_main": {"MoCul_aGMCSF": 0.5, "CoCul_aGMCSF": 0.5, "MoCul_IgG": -0.5, "CoCul_IgG": -0.5},
    "C3_interaction":   {"CoCul_aGMCSF": 1.0, "MoCul_aGMCSF": -1.0, "CoCul_IgG": -1.0, "MoCul_IgG": 1.0},
    "C4_co_vs_mono_IgG": {"CoCul_IgG": 1.0, "MoCul_IgG": -1.0},
    "C5_co_vs_mono_aGMCSF": {"CoCul_aGMCSF": 1.0, "MoCul_aGMCSF": -1.0},
    "C6_aGMCSF_vs_IgG_mono": {"MoCul_aGMCSF": 1.0, "MoCul_IgG": -1.0},
    "C7_aGMCSF_vs_IgG_co":   {"CoCul_aGMCSF": 1.0, "CoCul_IgG": -1.0},
}
run_block("data1", d1, meta1, contrasts1, OUT)

# ============ data2 ============
print("===== data2 =====")
d2 = load_counts("GSE309039_data2_raw_counts.csv.gz")
meta2 = build_meta(d2.columns)
meta2.to_csv(os.path.join(OUT, "data2_metadata.csv"), index=False)

contrasts2 = {
    "C8_GMCSF_vs_ctrl":   {"MoCul_GMCSF": 1.0, "MoCul": -1.0},
    "C9_LN229_vs_mono":   {"CoCul_LN229": 1.0, "MoCul": -1.0},
    "C10_U251_vs_mono":   {"CoCul_U251": 1.0, "MoCul": -1.0},
    "C11_LN229_vs_U251":  {"CoCul_LN229": 1.0, "CoCul_U251": -1.0},
    "C12_pooled_co_vs_mono": {"CoCul_LN229": 0.5, "CoCul_U251": 0.5, "MoCul": -1.0},
}
run_block("data2", d2, meta2, contrasts2, OUT)

print("===== B1 全部完成 =====")
