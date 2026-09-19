# -*- coding: utf-8 -*-
"""B7 模块六：rescue 模块分析（探索性）
strict16 / phenotypic31 / mod15∩rescue 在 MO-TAM vs MG-TAM 的患者内配对差异
探索性标注，不得根据患者结果反向选择最佳模块
"""
import os, pickle
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
WS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(WS, "output", "B7")
TMP = os.path.join(WS, "temp", "B7")

ra = pd.read_csv(os.path.join(WS, "output", "B5", "rescue_analysis_full.csv"))
ann68 = pd.read_csv(os.path.join(OUT, "core68_reannotation.csv"))
m = ann68.merge(ra[["symbol","rescue_class","rescue_significant","rescue_index"]], on="symbol", how="left")
strict16 = m[m["rescue_significant"]==True]["symbol"].tolist()
pheno31 = m[m["rescue_class"]=="fully_rescued"]["symbol"].tolist()
mod15 = m[m["module"]=="conserved_canonical_IFN_15"]["symbol"].tolist()
rescue_in_mod15 = [g for g in mod15 if g in set(strict16)]
pheno_in_mod15 = [g for g in mod15 if g in set(pheno31)]
MODULES = {"strict16": strict16, "pheno31": pheno31,
           "mod15_cap_strict": rescue_in_mod15, "mod15_cap_pheno": pheno_in_mod15}

def pseudobulk_paired(expr, annot, lineage_col, sample_col, mod_genes, label):
    """患者内 MO-TAM vs MG-TAM pseudobulk 配对差异"""
    g = [x for x in mod_genes if x in expr.columns]
    rows = []
    for sample in sorted(annot[sample_col].unique()):
        for lin in ["MO-TAM", "MG-TAM"]:
            cells = annot[(annot[sample_col]==sample) & (annot[lineage_col]==lin)].index
            cells = [c for c in cells if c in expr.index]
            if len(cells) < 3:
                continue
            score = expr.loc[cells, g].mean().mean()  # 平均表达
            rows.append({"sample": sample, "lineage": lin, "n_cells": len(cells), "score": score})
    df = pd.DataFrame(rows)
    piv = df.pivot(index="sample", columns="lineage", values="score")
    piv = piv.dropna()
    if len(piv) < 3:
        return None, piv
    diff = piv["MO-TAM"] - piv["MG-TAM"]
    try:
        w, p = wilcoxon(diff)
    except ValueError:
        w, p = np.nan, np.nan
    return {"module": label, "n_patients": len(piv), "mean_diff_MO_minus_MG": diff.mean(),
            "median_diff": diff.median(), "n_MO_gt_MG": int((diff>0).sum()),
            "wilcoxon_p": p}, piv

# ============ GSE163120 ============
with open(os.path.join(TMP, "gse163120_tam_expr.pkl"), "rb") as f:
    e163 = pickle.load(f)
a163 = pd.read_csv(os.path.join(TMP, "gse163120_tam_annot.csv")).set_index("Unnamed: 0")
res163 = []
for name, genes in MODULES.items():
    r, piv = pseudobulk_paired(e163, a163, "lineage", "sample", genes, name)
    if r:
        res163.append(r)
        print(f"[GSE163120] {name}: n_pat={r['n_patients']}, diff={r['mean_diff_MO_minus_MG']:.4f}, "
              f"MO>MG {r['n_MO_gt_MG']}/{r['n_patients']}, p={r['wilcoxon_p']:.4f}")
df163 = pd.DataFrame(res163)
df163.to_csv(os.path.join(OUT, "B7_m6_rescue_GSE163120.csv"), index=False)

# ============ GSE182109 ============
with open(os.path.join(TMP, "gse182109_myeloid_expr_aligned.pkl"), "rb") as f:
    e182 = pickle.load(f)
a182 = pd.read_csv(os.path.join(TMP, "gse182109_myeloid_annot.csv"))
a182 = a182.set_index("barcode")
res182 = []
for name, genes in MODULES.items():
    r, piv = pseudobulk_paired(e182, a182, "lineage", "patient", genes, name)
    if r:
        res182.append(r)
        print(f"[GSE182109] {name}: n_pat={r['n_patients']}, diff={r['mean_diff_MO_minus_MG']:.4f}, "
              f"MO>MG {r['n_MO_gt_MG']}/{r['n_patients']}, p={r['wilcoxon_p']:.4f}")
df182 = pd.DataFrame(res182)
df182.to_csv(os.path.join(OUT, "B7_m6_rescue_GSE182109.csv"), index=False)

# 汇总
df163["cohort"] = "GSE163120"
df182["cohort"] = "GSE182109"
summary = pd.concat([df163, df182], ignore_index=True)
summary.to_csv(os.path.join(OUT, "B7_m6_rescue_summary.csv"), index=False)
print("\n已保存 B7_m6_rescue_*.csv")
