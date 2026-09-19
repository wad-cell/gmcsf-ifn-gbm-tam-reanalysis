# -*- coding: utf-8 -*-
"""B3: GSE309038 细胞系表达谱分析
模型: ~ condition (4 细胞系, 各 n=3 独立培养, 非供者配对)
对比 (参考水平 LN229):
  B3_T98G_vs_LN229
  B3_U87_vs_LN229
  B3_U251_vs_LN229
重点: CSF2(GM-CSF)/TGFB1/2/3/IFNAR1/2/TGFBR1/2/SMAD 等分泌因子与通路基因表达,
     解释 LN229 vs U251 对单核细胞 ISG 抑制强度差异。
"""
import os, time, pickle, json
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

RAW = os.path.join(PROJ_ROOT, "temp", "GSE309038", "raw")
OUT = os.path.join(PROJ_ROOT, "output", "B3")
ANNOT = os.path.join(PROJ_ROOT, "output", "B1", "gene_annotation.csv")
os.makedirs(OUT, exist_ok=True)
N_CPUS = 4

# 关键基因（分泌因子 / 受体 / 通路）
KEY_GENES = ["CSF2","CSF1","IL6","IL1B","CCL2","CXCL10","TGFB1","TGFB2","TGFB3",
             "IFNAR1","IFNAR2","TGFBR1","TGFBR2","TGFBR3","SMAD2","SMAD3","SMAD4",
             "STAT1","STAT2","IRF1","IRF7","IFNB1","IFNA1","IFNG","TNF","VEGFA","MIF"]

def load_counts(fname):
    return pd.read_csv(os.path.join(RAW, fname), sep=";", index_col=0)

def build_meta(cols):
    meta = pd.DataFrame(index=cols)
    meta["condition"] = [c.rsplit("_n", 1)[0] for c in cols]
    meta["rep"] = [c.rsplit("_n", 1)[1] for c in cols]
    return meta

# ============ 主流程 ============
print("===== GSE309038 (B3) =====")
d = load_counts("GSE309038_raw_counts.csv.gz")
meta = build_meta(d.columns)
meta.to_csv(os.path.join(OUT, "GSE309038_metadata.csv"), index=False)
print("metadata:")
print(meta.to_string(index=False))

# 参考水平 LN229
dds_path = os.path.join(OUT, "GSE309038_dds.pkl")
if os.path.exists(dds_path):
    with open(dds_path, "rb") as f:
        dds = pickle.load(f)
    print("加载已拟合 dds")
else:
    print(f"拟合 DESeq2: {d.shape[0]} genes x {d.shape[1]} samples")
    dds = DeseqDataSet(counts=d.T, metadata=meta,
                       design_factors=["condition"],
                       refit_cooks=True, n_cpus=N_CPUS)
    dds.deseq2()
    with open(dds_path, "wb") as f:
        pickle.dump(dds, f)
    print("拟合完成")

dm_cols = list(dds.obsm["design_matrix"].columns)
print("设计矩阵列:", dm_cols)
colmap = {c: i for i, c in enumerate(dm_cols)}
with open(os.path.join(OUT, "GSE309038_design_cols.json"), "w") as f:
    json.dump(dm_cols, f)

contrasts = {
    "B3_T98G_vs_LN229": {"T98G": 1.0, "LN229": -1.0},
    "B3_U87_vs_LN229":  {"U87": 1.0, "LN229": -1.0},
    "B3_U251_vs_LN229": {"U251": 1.0, "LN229": -1.0},
}
for cid, level_coefs in contrasts.items():
    out_csv = os.path.join(OUT, f"GSE309038_{cid}.csv")
    if os.path.exists(out_csv):
        print(f"{cid} 已存在, 跳过")
        continue
    v = np.zeros(len(dm_cols))
    for level, coef in level_coefs.items():
        col = f"condition[T.{level}]"
        if col in colmap:
            v[colmap[col]] = coef
        else:
            print(f"{cid}: 水平 {level} 为参考水平(系数0), 忽略")
    print(f"计算 {cid}: {v}")
    t0 = time.time()
    stats = DeseqStats(dds, contrast=v, n_cpus=N_CPUS, quiet=True)
    stats.summary()
    res = stats.results_df.copy()
    res.insert(0, "gene", res.index)
    res.to_csv(out_csv, index=False)
    print(f"{cid} 完成: {res.shape[0]} genes, {time.time()-t0:.1f}s")

# 注释基因符号
annot = pd.read_csv(ANNOT)
for cid in contrasts:
    f = os.path.join(OUT, f"GSE309038_{cid}.csv")
    df = pd.read_csv(f)
    if "symbol" in df.columns:
        continue
    df["ensembl_id"] = df["gene"].str.split(".").str[0]
    df = df.merge(annot, on="ensembl_id", how="left")
    cols = ["gene", "symbol"] + [c for c in df.columns if c not in ("gene", "symbol", "ensembl_id")]
    df = df[cols]
    df.to_csv(f, index=False)
    print(f"注释 {cid}: {len(df)} 行")

# ============ 分泌因子表达矩阵 ============
# 用 DESeq2 归一化计数 (median-of-ratios) 计算各细胞系平均表达
norm = dds.layers["normed_counts"]  # samples x genes
norm_df = pd.DataFrame(norm, index=dds.obsm["design_matrix"].index, columns=dds.var_names)
norm_df = norm_df.T  # genes x samples
norm_df["ensembl_id"] = norm_df.index.str.split(".").str[0]
norm_df = norm_df.merge(annot, on="ensembl_id", how="left")

# 关键基因表达矩阵（每细胞系均值 + 每样本）
key_rows = norm_df[norm_df["symbol"].isin(KEY_GENES)].copy()
key_rows = key_rows.set_index("symbol")
key_rows = key_rows.drop(columns=["ensembl_id"], errors="ignore")
key_rows = key_rows[list(meta.index)]  # 按样本顺序
key_rows.to_csv(os.path.join(OUT, "GSE309038_key_genes_norm_counts.csv"))

# 每细胞系均值（log2(norm+1)）
log2n = np.log2(key_rows + 1.0)
cell_means = pd.DataFrame(index=log2n.index)
for cl in ["LN229", "T98G", "U87", "U251"]:
    cols = [c for c in log2n.columns if c.startswith(cl + "_n")]
    cell_means[cl] = log2n[cols].mean(axis=1)
cell_means.to_csv(os.path.join(OUT, "GSE309038_key_genes_cell_mean_log2.csv"))

print("\n===== 关键基因各细胞系平均表达 (log2(norm+1)) =====")
print(cell_means.round(2).to_string())

# 全基因细胞系均值（供 B4c 关联）
all_log2 = np.log2(norm_df.set_index("symbol").drop(columns=["ensembl_id"], errors="ignore") + 1.0)
all_means = pd.DataFrame(index=all_log2.index)
for cl in ["LN229", "T98G", "U87", "U251"]:
    cols = [c for c in all_log2.columns if c.startswith(cl + "_n")]
    all_means[cl] = all_log2[cols].mean(axis=1)
all_means.to_csv(os.path.join(OUT, "GSE309038_all_genes_cell_mean_log2.csv"))
print(f"\n全基因细胞系均值矩阵: {all_means.shape}")

print("===== B3 完成 =====")
