# -*- coding: utf-8 -*-
"""B2: GSE309037 分析（IFN-β 金标准 ISG + 直接共培养验证）
模型: ~ donor + condition (3 供者配对, 3 条件)
对比:
  B2_IFN_vs_MoCul:  IFN-β 1ng/mL 6h vs 单核对照 (金标准 ISG 诱导)
  B2_CoCul_vs_MoCul: LN229 直接共培养 48h vs 单核对照
  B2_CoCul_vs_IFN:  共培养 vs IFN-β 处理
"""
import os, time, pickle, json
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

RAW = os.path.join(PROJ_ROOT, "temp", "GSE309037", "raw")
OUT = os.path.join(PROJ_ROOT, "output", "B2")
ANNOT = os.path.join(PROJ_ROOT, "output", "B1", "gene_annotation.csv")
os.makedirs(OUT, exist_ok=True)
N_CPUS = 4

def load_counts(fname):
    return pd.read_csv(os.path.join(RAW, fname), sep=";", index_col=0)

def build_meta(cols):
    meta = pd.DataFrame(index=cols)
    meta["condition"] = [c.rsplit("_n", 1)[0] for c in cols]
    meta["donor"] = ["donor_" + c.split("_n")[1] for c in cols]
    return meta

def run_block(block, counts, meta, contrasts, outdir):
    dds_path = os.path.join(outdir, f"{block}_dds.pkl")
    if os.path.exists(dds_path):
        with open(dds_path, "rb") as f:
            dds = pickle.load(f)
        print(f"[{block}] 加载已拟合 dds")
    else:
        print(f"[{block}] 拟合 DESeq2: {counts.shape[0]} genes x {counts.shape[1]} samples")
        dds = DeseqDataSet(counts=counts.T, metadata=meta,
                           design_factors=["donor", "condition"],
                           refit_cooks=True, n_cpus=N_CPUS)
        dds.deseq2()
        with open(dds_path, "wb") as f:
            pickle.dump(dds, f)
        print(f"[{block}] 拟合完成")

    dm_cols = list(dds.obsm["design_matrix"].columns)
    print(f"[{block}] 设计矩阵列: {dm_cols}")
    colmap = {c: i for i, c in enumerate(dm_cols)}
    with open(os.path.join(outdir, f"{block}_design_cols.json"), "w") as f:
        json.dump(dm_cols, f)

    for cid, level_coefs in contrasts.items():
        out_csv = os.path.join(outdir, f"{block}_{cid}.csv")
        if os.path.exists(out_csv):
            print(f"[{block}] {cid} 已存在, 跳过")
            continue
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

# ============ 主流程 ============
print("===== GSE309037 (B2) =====")
d = load_counts("GSE309037_raw_counts.csv.gz")
meta = build_meta(d.columns)
meta.to_csv(os.path.join(OUT, "GSE309037_metadata.csv"), index=False)
print("metadata:")
print(meta.to_string(index=False))

contrasts = {
    "B2_IFN_vs_MoCul":  {"MoCul_IFN": 1.0, "MoCul": -1.0},
    "B2_CoCul_vs_MoCul": {"CoCul": 1.0, "MoCul": -1.0},
    "B2_CoCul_vs_IFN":  {"CoCul": 1.0, "MoCul_IFN": -1.0},
}
run_block("GSE309037", d, meta, contrasts, OUT)

# 注释基因符号
annot = pd.read_csv(ANNOT)
for cid in contrasts:
    f = os.path.join(OUT, f"GSE309037_{cid}.csv")
    df = pd.read_csv(f)
    if "symbol" in df.columns:
        continue
    df["ensembl_id"] = df["gene"].str.split(".").str[0]
    df = df.merge(annot, on="ensembl_id", how="left")
    cols = ["gene", "symbol"] + [c for c in df.columns if c not in ("gene", "symbol", "ensembl_id")]
    df = df[cols]
    df.to_csv(f, index=False)
    print(f"注释 {cid}: {len(df)} 行")

print("===== B2 完成 =====")
