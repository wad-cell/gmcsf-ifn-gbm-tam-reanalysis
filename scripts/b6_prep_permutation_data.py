# -*- coding: utf-8 -*-
"""B6 模块一/二: 准备数据 (修正版)"""
import os, json, pickle
import numpy as np, pandas as pd


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
TEMP = os.path.join(WS, "temp")
B2 = os.path.join(WS, "output", "B2")
B1 = os.path.join(WS, "output", "B1")

deseq = pd.read_csv(os.path.join(B2, "GSE309037_B2_CoCul_vs_MoCul.csv"))
deseq2 = deseq.dropna(subset=["baseMean", "log2FoldChange"]).copy()
deseq2["symbol"] = deseq2["symbol"].astype(str)
deseq_dedup = deseq2.sort_values("baseMean", ascending=False).drop_duplicates("symbol", keep="first").copy()
print("dedup bg:", len(deseq_dedup))

# ENSG -> symbol 映射
ga = pd.read_csv(os.path.join(B1, "gene_annotation.csv"))
ens2sym = ga.set_index("ensembl_id")["symbol"].to_dict()

# 检测率
with open(os.path.join(B2, "GSE309037_dds.pkl"), "rb") as f:
    dds = pickle.load(f)
cnt = np.asarray(dds.X)  # (9, 62703)
detect = (cnt > 0).mean(axis=0)
detect_df = pd.DataFrame({"ens": dds.var_names.str.split(".").str[0], "detection_rate": detect})
detect_df["symbol"] = detect_df["ens"].map(ens2sym)
detect_df = detect_df.dropna(subset=["symbol"])
detect_dedup = detect_df.sort_values("detection_rate", ascending=False).drop_duplicates("symbol", keep="first")
detect_map = detect_dedup.set_index("symbol")["detection_rate"]
deseq_dedup["detection_rate"] = deseq_dedup["symbol"].map(detect_map)
print("detection mapped:", deseq_dedup["detection_rate"].notna().mean())

# 基因长度
with open(os.path.join(TEMP, "gene_lengths_cache.json")) as f:
    cache = json.load(f)
sym2ens = ga.set_index("symbol")["ensembl_id"].to_dict()
deseq_dedup["ens"] = deseq_dedup["symbol"].map(sym2ens)
deseq_dedup["gene_len"] = deseq_dedup["ens"].map(cache)
print("gene_len mapped:", deseq_dedup["gene_len"].notna().mean())

# 经验离散度: counts 转置 (62703, 9), 按 symbol 聚合
cntT = pd.DataFrame(cnt.T, index=dds.var_names)  # (62703, 9)
cntT["ens"] = cntT.index.str.split(".").str[0]
cntT["symbol"] = cntT["ens"].map(ens2sym)
cntT = cntT.dropna(subset=["symbol"])
num_cols = [c for c in cntT.columns if c not in ("ens", "symbol")]
cnt_sym = cntT.groupby("symbol")[num_cols].sum()
means = cnt_sym.mean(axis=1)
vars_ = cnt_sym.var(axis=1, ddof=1)
disp = (vars_ - means) / (means ** 2)
disp = disp.replace([np.inf, -np.inf], np.nan)
deseq_dedup["emp_dispersion"] = deseq_dedup["symbol"].map(disp)
print("dispersion mapped:", deseq_dedup["emp_dispersion"].notna().mean())

deseq_dedup.to_csv(os.path.join(TEMP, "b6_bg_genes_annotated.csv"), index=False)
print("saved, shape:", deseq_dedup.shape)
print(deseq_dedup[["symbol","baseMean","log2FoldChange","detection_rate","gene_len","emp_dispersion"]].head(5).to_string())
