# -*- coding: utf-8 -*-
"""B5 模块二：68 基因核心集来源审计
输出 core68_evidence_matrix.csv：每个核心基因在所有 contrast 的
baseMean / log2FC / lfcSE / stat / pvalue / padj / direction
不允许选择性报告，输出完整结果。
"""
import pandas as pd, os

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
OUT = os.path.join(WS, "output", "B5")
B1 = os.path.join(WS, "output", "B1")
B2 = os.path.join(WS, "output", "B2")

# 68 核心集
core = pd.read_csv(os.path.join(B2, "overlap_IFN_ISG_vs_CoCul_down.csv"))
core_syms = sorted(core["symbol"].dropna().unique())
print("68 core genes:", len(core_syms))

# 所有 contrast 文件
contrast_files = [
    # (block, contrast_id, 文件路径, 方向语义)
    ("GSE309037", "B2_IFN_vs_MoCul", os.path.join(B2, "GSE309037_B2_IFN_vs_MoCul.csv")),
    ("GSE309037", "B2_CoCul_vs_MoCul", os.path.join(B2, "GSE309037_B2_CoCul_vs_MoCul.csv")),
    ("GSE309037", "B2_CoCul_vs_IFN", os.path.join(B2, "GSE309037_B2_CoCul_vs_IFN.csv")),
    ("GSE309039", "C1_culture_main", os.path.join(B1, "data1_C1_culture_main.csv")),
    ("GSE309039", "C2_antibody_main", os.path.join(B1, "data1_C2_antibody_main.csv")),
    ("GSE309039", "C3_interaction", os.path.join(B1, "data1_C3_interaction.csv")),
    ("GSE309039", "C4_co_vs_mono_IgG", os.path.join(B1, "data1_C4_co_vs_mono_IgG.csv")),
    ("GSE309039", "C5_co_vs_mono_aGMCSF", os.path.join(B1, "data1_C5_co_vs_mono_aGMCSF.csv")),
    ("GSE309039", "C6_aGMCSF_vs_IgG_mono", os.path.join(B1, "data1_C6_aGMCSF_vs_IgG_mono.csv")),
    ("GSE309039", "C7_aGMCSF_vs_IgG_co", os.path.join(B1, "data1_C7_aGMCSF_vs_IgG_co.csv")),
    ("GSE309039", "C8_GMCSF_vs_ctrl", os.path.join(B1, "data2_C8_GMCSF_vs_ctrl.csv")),
    ("GSE309039", "C9_LN229_vs_mono", os.path.join(B1, "data2_C9_LN229_vs_mono.csv")),
    ("GSE309039", "C10_U251_vs_mono", os.path.join(B1, "data2_C10_U251_vs_mono.csv")),
    ("GSE309039", "C11_LN229_vs_U251", os.path.join(B1, "data2_C11_LN229_vs_U251.csv")),
    ("GSE309039", "C12_pooled_co_vs_mono", os.path.join(B1, "data2_C12_pooled_co_vs_mono.csv")),
]

# 逐 contrast 提取 68 基因统计量
matrices = {}
for block, cid, path in contrast_files:
    df = pd.read_csv(path)
    df = df.drop_duplicates(subset="symbol", keep="first")
    df = df.set_index("symbol")
    sub = df.loc[df.index.isin(core_syms), ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]].copy()
    sub.columns = [f"{cid}_baseMean", f"{cid}_log2FC", f"{cid}_lfcSE", f"{cid}_stat", f"{cid}_pvalue", f"{cid}_padj"]
    matrices[cid] = sub

# 合并
ev = pd.DataFrame(index=core_syms)
for cid in matrices:
    ev = ev.join(matrices[cid])
ev = ev.reset_index().rename(columns={"index": "symbol"})

# 计算每个 contrast 的 direction (padj<0.05 & |log2FC|>1 为显著; 否则按 padj<0.05 方向; 否则 ns; NaN 为 na)
for _, cid, _ in contrast_files:
    lfc = ev[f"{cid}_log2FC"]; padj = ev[f"{cid}_padj"]
    def direction(l, p):
        if pd.isna(p) or pd.isna(l):
            return "na"
        if p < 0.05 and abs(l) > 1:
            return "up" if l > 0 else "down"
        if p < 0.05:
            return "up_sig" if l > 0 else "down_sig"
        return "ns"
    ev[f"{cid}_direction"] = [direction(l, p) for l, p in zip(lfc, padj)]

# 附加 68 核心集来源信息
core_info = core[["symbol", "IFN_log2FC", "IFN_padj", "CoCul_log2FC", "CoCul_padj"]].drop_duplicates(subset="symbol")
ev = ev.merge(core_info, on="symbol", how="left")

# 保存
ev.to_csv(os.path.join(OUT, "core68_evidence_matrix.csv"), index=False)
print("core68_evidence_matrix.csv:", ev.shape)
print("columns:", len(ev.columns))
# 摘要：每个 contrast 的 68 基因方向分布
print("\n=== 68 基因在各 contrast 的方向分布 ===")
for _, cid, _ in contrast_files:
    vc = ev[f"{cid}_direction"].value_counts()
    print(f"{cid}: " + ", ".join(f"{k}={v}" for k, v in vc.items()))
