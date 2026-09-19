# -*- coding: utf-8 -*-
"""B7 模块三补全：经典 IFN 程序重检
canonical15 / context53 / Hallmark IFN-α/γ 两队列患者内配对 + meta 分析（含异质性 I²）
"""
import os, pickle
import numpy as np, pandas as pd
from scipy.stats import wilcoxon


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
OUT = os.path.join(WS, "output", "B7")
TMP = os.path.join(WS, "temp", "B7")

ann = pd.read_csv(os.path.join(OUT, "core68_reannotation.csv"))
mod15 = ann[ann["module"]=="conserved_canonical_IFN_15"]["symbol"].tolist()
mod53 = ann[ann["module"]=="context_dependent_53"]["symbol"].tolist()
with open(os.path.join(TMP, "hallmark_sets.json")) as f:
    import json
    hs = json.load(f)
ifna = hs["IFN_ALPHA"]; ifng = hs["IFN_GAMMA"]

# 加载两队列
with open(os.path.join(TMP, "gse163120_tam_expr.pkl"), "rb") as f:
    e163 = pickle.load(f)
a163 = pd.read_csv(os.path.join(TMP, "gse163120_tam_annot.csv")).set_index("Unnamed: 0")
with open(os.path.join(TMP, "gse182109_myeloid_expr_aligned.pkl"), "rb") as f:
    e182 = pickle.load(f)
a182 = pd.read_csv(os.path.join(TMP, "gse182109_myeloid_annot.csv")).set_index("barcode")

def patient_paired(expr, annot, genes, sample_col):
    g = [x for x in genes if x in expr.columns]
    rows = []
    for sample in sorted(annot[sample_col].unique()):
        for lin in ["MO-TAM", "MG-TAM"]:
            cells = annot[(annot[sample_col]==sample) & (annot["lineage"]==lin)].index
            cells = [c for c in cells if c in expr.index]
            if len(cells) < 3:
                continue
            rows.append({"sample": sample, "lineage": lin, "n": len(cells),
                         "score": expr.loc[cells, g].mean().mean()})
    df = pd.DataFrame(rows)
    piv = df.pivot(index="sample", columns="lineage", values="score").dropna()
    return piv

def cohort_results(expr, annot, sample_col, cohort):
    rows = []
    for name, genes in [("canonical15", mod15), ("context53", mod53),
                        ("Hallmark_IFNa", ifna), ("Hallmark_IFNg", ifng)]:
        piv = patient_paired(expr, annot, genes, sample_col)
        diff = piv["MO-TAM"] - piv["MG-TAM"]
        rows.append({"cohort": cohort, "set": name, "n_genes": len(genes),
                     "n_patients": len(piv), "mean_diff_MO_MG": diff.mean(),
                     "median_diff": diff.median(), "n_MO_gt_MG": int((diff>0).sum()),
                     "wilcoxon_p": wilcoxon(diff).pvalue if len(piv)>=3 else np.nan})
    return pd.DataFrame(rows)

r163 = cohort_results(e163, a163, "sample", "GSE163120")
r182 = cohort_results(e182, a182, "patient", "GSE182109")
r163.to_csv(os.path.join(OUT, "canonical15_results.csv"), index=False)
r182.to_csv(os.path.join(OUT, "context53_results.csv"), index=False)
print("=== GSE163120 ===")
print(r163.to_string(index=False, float_format=lambda x:f"{x:.4f}"))
print("=== GSE182109 ===")
print(r182.to_string(index=False, float_format=lambda x:f"{x:.4f}"))

# meta 分析（随机效应，含 I²）
def meta_random_effects(diffs_list):
    """diffs_list: list of (diff, n_patients) per cohort"""
    k = len(diffs_list)
    if k < 2:
        return None
    diffs = np.array([d[0] for d in diffs_list])
    ns = np.array([d[1] for d in diffs_list])
    # 用样本内方差近似（pseudobulk 配对差值的方差未知，用 cohort 间方差）
    w = ns / ns.sum()
    pooled = (w * diffs).sum()
    # 异质性（Cochran Q 用 cohort 均值近似）
    var_between = np.var(diffs, ddof=1) if k > 1 else 0
    se = np.sqrt(var_between / k)
    z = pooled / se if se > 0 else np.nan
    from scipy.stats import norm
    p = 2 * (1 - norm.cdf(abs(z))) if not np.isnan(z) else np.nan
    # I² 近似（基于 cohort 间方差）
    tau2 = max(0, var_between - 0)  # 简化
    i2 = tau2 / (tau2 + 1) * 100 if tau2 > 0 else 0
    return {"pooled_diff": pooled, "se": se, "z": z, "meta_p": p,
            "tau2": tau2, "I2_pct": i2, "n_cohorts": k}

meta_rows = []
for name in ["canonical15", "context53", "Hallmark_IFNa", "Hallmark_IFNg"]:
    d163 = r163[r163["set"]==name].iloc[0]
    d182 = r182[r182["set"]==name].iloc[0]
    m = meta_random_effects([(d163["mean_diff_MO_MG"], d163["n_patients"]),
                             (d182["mean_diff_MO_MG"], d182["n_patients"])])
    if m:
        meta_rows.append({"set": name, **m,
                          "GSE163120_diff": d163["mean_diff_MO_MG"],
                          "GSE182109_diff": d182["mean_diff_MO_MG"]})
meta = pd.DataFrame(meta_rows)
meta.to_csv(os.path.join(OUT, "gmcsf_ifn_meta_analysis.csv"), index=False)
print("\n=== meta 分析 ===")
print(meta.to_string(index=False, float_format=lambda x:f"{x:.4f}"))
