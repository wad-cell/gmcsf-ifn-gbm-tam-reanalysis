# -*- coding: utf-8 -*-
"""B7 模块五：状态空间探索（探索性）
hypoxia / MES / inflammatory-mature / CISH-high-low 与 mod15、context53 的关系
FDR 校正，探索性标注，不升级为预注册主终点
"""
import os, pickle, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
OUT = os.path.join(WS, "output", "B7")
TMP = os.path.join(WS, "temp", "B7")

with open(os.path.join(TMP, "gse163120_tam_expr.pkl"), "rb") as f:
    expr = pickle.load(f)
annot = pd.read_csv(os.path.join(TMP, "gse163120_tam_annot.csv")).set_index("Unnamed: 0")
with open(os.path.join(TMP, "hallmark_sets.json")) as f:
    hs = json.load(f)

# 状态基因集（预定义，探索性）
STATE_SETS = {
    "hypoxia": hs["HYPOXIA"],
    "MES_Verhaak": ["CHI3L1","CD44","TGFBI","SERPINE1","TIMP1","FN1","COL1A1","COL1A2","SPARC","VIM","TAGLN","FAP"],
    "mature_TAM": ["C1QA","C1QB","C1QC","APOE","TREM2","CD68","MRC1","MSR1","GPNMB","FOLR2"],
    "inflammatory_TAM": ["IL1B","TNF","CXCL10","CCL2","CCL3","CCL4","IL6","CXCL9","CXCL11","NFKBIA"],
}

# mod15 / context53 成员
ann68 = pd.read_csv(os.path.join(OUT, "core68_reannotation.csv"))
mod15 = ann68[ann68["module"]=="conserved_canonical_IFN_15"]["symbol"].tolist()
ctx53 = ann68[ann68["module"]=="context_dependent_53"]["symbol"].tolist()

def mean_score(df, genes):
    g = [x for x in genes if x in df.columns]
    return df[g].mean(axis=1)

def ucell_score(df, genes):
    """UCell：基因集内基因在细胞所有基因中的相对排名均值（简化实现）"""
    g = [x for x in genes if x in df.columns]
    if not g:
        return pd.Series(np.nan, index=df.index)
    ranks = df.rank(axis=1, pct=True)
    return ranks[g].mean(axis=1)

# 计算各细胞得分（MO-TAM + MG-TAM 全 TAM，用于状态探索）
scores = pd.DataFrame(index=expr.index)
scores["lineage"] = annot["lineage"]
scores["sample"] = annot["sample"]
scores["mod15_mean"] = mean_score(expr, mod15)
scores["ctx53_mean"] = mean_score(expr, ctx53)
scores["mod15_ucell"] = ucell_score(expr, mod15)
scores["CISH"] = expr["CISH"] if "CISH" in expr.columns else 0
for st, genes in STATE_SETS.items():
    scores[f"{st}_mean"] = mean_score(expr, genes)
    scores[f"{st}_ucell"] = ucell_score(expr, genes)
scores["total_expr"] = expr.sum(axis=1)

def partial_spearman(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def resid(a, b):
        A = np.vstack([b, np.ones(len(b))]).T
        coef, *_ = np.linalg.lstsq(A, a, rcond=None)
        return a - A @ coef
    return spearmanr(resid(rx, rz), resid(ry, rz))

# ============ 分析 1：状态得分与 mod15/ctx53 相关（全 TAM + MO-TAM 分层）============
rows = []
for lineage in ["MO-TAM", "MG-TAM", "ALL"]:
    if lineage == "ALL":
        sub = scores
    else:
        sub = scores[scores["lineage"] == lineage]
    for st in STATE_SETS:
        for target, tcol in [("mod15", "mod15_mean"), ("ctx53", "ctx53_mean")]:
            d = sub[[f"{st}_mean", tcol, "total_expr"]].dropna()
            r0, p0 = spearmanr(d[f"{st}_mean"], d[tcol])
            r1, p1 = partial_spearman(d[f"{st}_mean"], d[tcol], d["total_expr"])
            rows.append({"lineage": lineage, "state": st, "target": target, "n": len(d),
                         "spearman_rho": r0, "p": p0, "partial_rho": r1, "partial_p": p1})
res1 = pd.DataFrame(rows)
res1["partial_p_fdr"] = multipletests(res1["partial_p"], method="fdr_bh")[1]
res1 = res1.sort_values("partial_p")
res1.to_csv(os.path.join(OUT, "B7_m5_state_corr.csv"), index=False)
print("=== 分析1：状态 vs mod15/ctx53 相关（控制 total_expr）===")
print(res1[res1["lineage"]=="MO-TAM"].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

# ============ 分析 2：CISH-high/low 分组比较 mod15（MO-TAM）============
mo = scores[scores["lineage"]=="MO-TAM"].copy()
# CISH 分组：表达>0 为 high，=0 为 low
mo["CISH_group"] = np.where(mo["CISH"] > 0, "CISH_high", "CISH_low")
g_high = mo[mo["CISH_group"]=="CISH_high"]["mod15_mean"]
g_low = mo[mo["CISH_group"]=="CISH_low"]["mod15_mean"]
u, p = mannwhitneyu(g_high, g_low, alternative="two-sided")
print(f"\n=== 分析2：CISH-high (n={len(g_high)}) vs CISH-low (n={len(g_low)}) mod15_mean ===")
print(f"high median={g_high.median():.4f}, low median={g_low.median():.4f}, MWU p={p:.3e}")
# 按样本分层
strat_rows = []
for sample in sorted(mo["sample"].unique()):
    s = mo[mo["sample"]==sample]
    gh = s[s["CISH_group"]=="CISH_high"]["mod15_mean"]
    gl = s[s["CISH_group"]=="CISH_low"]["mod15_mean"]
    if len(gh)>=3 and len(gl)>=3:
        u2,p2 = mannwhitneyu(gh, gl, alternative="two-sided")
        strat_rows.append({"sample":sample,"n_high":len(gh),"n_low":len(gl),
                           "high_median":gh.median(),"low_median":gl.median(),
                           "diff_high_minus_low":gh.median()-gl.median(),"p":p2})
strat = pd.DataFrame(strat_rows)
strat["p_fdr"] = multipletests(strat["p"], method="fdr_bh")[1]
strat.to_csv(os.path.join(OUT, "B7_m5_cish_group_stratified.csv"), index=False)
print("\n按样本分层 CISH-high vs low mod15:")
print(strat.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

# ============ 分析 3：状态分组（hypoxia/MES/inflammatory 高 vs 低）比较 mod15 ============
grp_rows = []
for st in STATE_SETS:
    med = mo[f"{st}_mean"].median()
    hi = mo[mo[f"{st}_mean"] > med]["mod15_mean"]
    lo = mo[mo[f"{st}_mean"] <= med]["mod15_mean"]
    u2, p2 = mannwhitneyu(hi, lo, alternative="two-sided")
    grp_rows.append({"state": st, "n_high": len(hi), "n_low": len(lo),
                     "high_median_mod15": hi.median(), "low_median_mod15": lo.median(),
                     "diff_high_minus_low": hi.median()-lo.median(), "p": p2})
grp = pd.DataFrame(grp_rows)
grp["p_fdr"] = multipletests(grp["p"], method="fdr_bh")[1]
grp.to_csv(os.path.join(OUT, "B7_m5_state_group_mod15.csv"), index=False)
print("\n=== 分析3：状态高/低分组 mod15 差异（MO-TAM）===")
print(grp.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

scores.to_csv(os.path.join(OUT, "B7_m5_cell_scores.csv"), index=False)
print("\n已保存: B7_m5_state_corr.csv / B7_m5_cish_group_stratified.csv / B7_m5_state_group_mod15.csv / B7_m5_cell_scores.csv")
