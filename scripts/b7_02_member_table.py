# -*- coding: utf-8 -*-
"""B7 模块一：core68 成员表 core68_reannotation.csv
为 68 基因标注：Hallmark IFN-alpha/gamma、体外 IFN-beta 处理证据(B2 gold standard)、
IFN-beta log2FC/padj、共培养抑制、GM-CSF 处理、alphaGM-CSF rescue、两队列检测率。
"""
import pandas as pd, numpy as np, os, json

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
TMP = os.path.join(WS, "temp", "B7")
OUT = os.path.join(WS, "output", "B7")
os.makedirs(OUT, exist_ok=True)

# 68 基因证据矩阵
ev = pd.read_csv(os.path.join(WS, "output", "B5", "core68_evidence_matrix.csv"))
# B2 gold standard (体外 IFN-beta 上调)
gs = pd.read_csv(os.path.join(WS, "output", "B2", "IFN_gold_standard_ISG_up.csv"))
gs_sym = set(gs["symbol"].dropna())
# Hallmark
hs = json.load(open(os.path.join(TMP, "hallmark_sets.json")))
ifn_a = set(hs["IFN_ALPHA"]); ifn_g = set(hs["IFN_GAMMA"])

# 检测率
mat163 = pd.read_pickle(os.path.join(TMP, "gse163120_tam_expr.pkl"))
mat182 = pd.read_pickle(os.path.join(TMP, "gse182109_myeloid_expr.pkl"))

rows = []
for _, r in ev.iterrows():
    sym = r["symbol"]
    det163 = float((mat163[sym] > 0).mean()) if sym in mat163.columns else np.nan
    det182 = float((mat182[sym] > 0).mean()) if sym in mat182.columns else np.nan
    rows.append({
        "symbol": sym,
        "module": "conserved_canonical_IFN_15" if sym in ifn_a else "context_dependent_53",
        "in_hallmark_IFNa": sym in ifn_a,
        "in_hallmark_IFNg": sym in ifn_g,
        "in_vitro_IFNb_up_B2_gold_standard": sym in gs_sym,
        "IFNb_log2FC": r.get("IFN_log2FC", np.nan),
        "IFNb_padj": r.get("IFN_padj", np.nan),
        "CoCul_log2FC": r.get("CoCul_log2FC", np.nan),
        "CoCul_padj": r.get("CoCul_padj", np.nan),
        "GMCSF_log2FC": r.get("C8_GMCSF_vs_ctrl_log2FC", np.nan),
        "GMCSF_padj": r.get("C8_GMCSF_vs_ctrl_padj", np.nan),
        "aGMCSF_rescue_co_log2FC": r.get("C7_aGMCSF_vs_IgG_co_log2FC", np.nan),
        "aGMCSF_rescue_co_padj": r.get("C7_aGMCSF_vs_IgG_co_padj", np.nan),
        "aGMCSF_mono_log2FC": r.get("C6_aGMCSF_vs_IgG_mono_log2FC", np.nan),
        "aGMCSF_mono_padj": r.get("C6_aGMCSF_vs_IgG_mono_padj", np.nan),
        "detection_GSE163120": det163,
        "detection_GSE182109": det182,
    })
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "core68_reannotation.csv"), index=False)
print("core68_reannotation.csv 已保存:", df.shape)
print("模块分布:", df["module"].value_counts().to_dict())
print("Hallmark IFN-gamma 命中:", df["in_hallmark_IFNg"].sum())
print("体外 IFN-beta 上调命中:", df["in_vitro_IFNb_up_B2_gold_standard"].sum())
print("\n15 基因模块成员:")
print(df[df["module"]=="conserved_canonical_IFN_15"]["symbol"].tolist())
print("\n53 基因模块成员:")
print(df[df["module"]=="context_dependent_53"]["symbol"].tolist())
