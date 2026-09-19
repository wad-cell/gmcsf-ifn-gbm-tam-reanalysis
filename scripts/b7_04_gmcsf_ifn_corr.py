# -*- coding: utf-8 -*-
"""B7 模块四：MO-TAM 内部 GM-CSF/STAT5/CISH 与经典 IFN 程序相关性
预定义 GM-CSF 模块（基于体外 C8 GM-CSF 处理显著基因，不用患者数据回测）
分层分析：按样本（ND1-ND7）分层
"""
import os, sys, pickle
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b7_scoring import mean_score, ucell_score


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
OUT = os.path.join(WS, "output", "B7")
os.makedirs(OUT, exist_ok=True)

# ---------- 1. 载入数据 ----------
with open(os.path.join(WS, "temp", "B7", "gse163120_tam_expr.pkl"), "rb") as f:
    expr = pickle.load(f)          # 16878 x 18774, log1p
annot = pd.read_csv(os.path.join(WS, "temp", "B7", "gse163120_tam_annot.csv"))
annot = annot.set_index("Unnamed: 0")
assert (expr.index == annot.index).all()

ann = pd.read_csv(os.path.join(WS, "output", "B7", "core68_reannotation.csv"))
mod15 = ann[ann["module"] == "conserved_canonical_IFN_15"]["symbol"].tolist()

ev = pd.read_csv(os.path.join(WS, "output", "B5", "core68_evidence_matrix.csv"))
c8_sig = ev[(ev["C8_GMCSF_vs_ctrl_padj"] < 0.05)]["symbol"].tolist()  # 22 下调
mod15_set = set(mod15)
gmcsf_nonifn = [g for g in c8_sig if g not in mod15_set]  # 14 非重叠

# STAT5 经典靶基因（GM-CSF 信号正向指标）
stat5_targets = ["CISH", "SOCS3", "PIM1", "MYC", "BCL2", "CCND1"]

# ---------- 2. 仅 MO-TAM 细胞 ----------
mo = annot[annot["lineage"] == "MO-TAM"].index
expr_mo = expr.loc[mo]
annot_mo = annot.loc[mo]
print(f"MO-TAM 细胞数: {len(mo)}")

# ---------- 3. 计算细胞得分 ----------
scores = pd.DataFrame(index=expr_mo.index)
scores["mod15_mean"] = mean_score(expr_mo, mod15)
scores["mod15_ucell"] = ucell_score(expr_mo, mod15)
scores["GMCSF_rep_22_mean"] = mean_score(expr_mo, c8_sig)          # 22 基因（含重叠）
scores["GMCSF_rep_nonIFN_14_mean"] = mean_score(expr_mo, gmcsf_nonifn)  # 14 非重叠
scores["STAT5_targets_mean"] = mean_score(expr_mo, stat5_targets)
scores["CISH"] = expr_mo["CISH"] if "CISH" in expr_mo.columns else np.nan
scores["sample"] = annot_mo["sample"].values

# 各模块实际检测基因数
for name, genes in [("mod15", mod15), ("GMCSF_rep_22", c8_sig),
                    ("GMCSF_rep_nonIFN_14", gmcsf_nonifn), ("STAT5_targets", stat5_targets)]:
    det = [g for g in genes if g in expr_mo.columns]
    print(f"{name}: 定义 {len(genes)} 基因, GSE163120 可检测 {len(det)}")

# ---------- 4. 主分析：全 MO-TAM 相关性 ----------
pairs = [
    ("GMCSF_rep_22_mean", "mod15_mean", "GM-CSF-repressed 22 (含重叠)"),
    ("GMCSF_rep_nonIFN_14_mean", "mod15_mean", "GM-CSF-repressed 14 (非IFN重叠)"),
    ("GMCSF_rep_nonIFN_14_mean", "mod15_ucell", "GM-CSF-repressed 14 vs mod15 UCell"),
    ("STAT5_targets_mean", "mod15_mean", "STAT5 靶基因 (CISH/SOCS3/PIM1/MYC/BCL2/CCND1)"),
    ("STAT5_targets_mean", "mod15_ucell", "STAT5 靶基因 vs mod15 UCell"),
    ("CISH", "mod15_mean", "CISH 单基因"),
]
rows = []
for a, b, label in pairs:
    d = scores[[a, b]].dropna()
    rho, p = spearmanr(d[a], d[b])
    rows.append({"comparison": label, "var_x": a, "var_y": b, "n": len(d),
                 "spearman_rho": rho, "p": p})
res = pd.DataFrame(rows)
# FDR 校正
from statsmodels.stats.multitest import multipletests
res["p_fdr"] = multipletests(res["p"], method="fdr_bh")[1]
res = res.sort_values("p")
print("\n=== 全 MO-TAM 相关性 ===")
print(res.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

# ---------- 5. 分层分析：按样本 ----------
strat_rows = []
for sample in sorted(scores["sample"].unique()):
    sub = scores[scores["sample"] == sample]
    for a, b, label in pairs:
        d = sub[[a, b]].dropna()
        if len(d) < 10:
            continue
        rho, p = spearmanr(d[a], d[b])
        strat_rows.append({"sample": sample, "n": len(d), "comparison": label,
                           "var_x": a, "var_y": b, "spearman_rho": rho, "p": p})
strat = pd.DataFrame(strat_rows)
strat["p_fdr"] = multipletests(strat["p"], method="fdr_bh")[1]
strat = strat.sort_values(["comparison", "sample"])
print("\n=== 按样本分层相关性 ===")
print(strat.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

# ---------- 6. 汇总：分层相关 meta（随机效应） ----------
# 对每个 comparison，用 Fisher z 转换汇总各样本 rho
def fisher_z(r):
    return 0.5 * np.log((1 + r) / (1 - r))
def inv_fisher(z):
    return (np.exp(2 * z) - 1) / (np.exp(2 * z) + 1)

meta_rows = []
for label in strat["comparison"].unique():
    sub = strat[strat["comparison"] == label]
    z = fisher_z(sub["spearman_rho"].values)
    w = sub["n"].values - 3
    zbar = np.sum(w * z) / np.sum(w)
    se = 1 / np.sqrt(np.sum(w))
    z_p = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(zbar / se)))
    # 异质性 Q
    Q = np.sum(w * (z - zbar) ** 2)
    df = len(sub) - 1
    from scipy.stats import chi2
    Q_p = 1 - chi2.cdf(Q, df) if df > 0 else np.nan
    meta_rows.append({"comparison": label, "n_samples": len(sub),
                      "pooled_rho": inv_fisher(zbar), "z_p": z_p,
                      "Q": Q, "Q_df": df, "Q_p": Q_p,
                      "rho_range": f"{sub['spearman_rho'].min():.3f}~{sub['spearman_rho'].max():.3f}"})
meta = pd.DataFrame(meta_rows)
meta["z_p_fdr"] = multipletests(meta["z_p"], method="fdr_bh")[1]
print("\n=== 分层相关随机效应汇总 ===")
print(meta.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

# ---------- 7. 保存 ----------
res.to_csv(os.path.join(OUT, "B7_m4_gmcsf_ifn_corr_overall.csv"), index=False)
strat.to_csv(os.path.join(OUT, "B7_m4_gmcsf_ifn_corr_stratified.csv"), index=False)
meta.to_csv(os.path.join(OUT, "B7_m4_gmcsf_ifn_corr_meta.csv"), index=False)
scores.to_csv(os.path.join(OUT, "B7_m4_cell_scores.csv"), index=False)
print("\n已保存: B7_m4_gmcsf_ifn_corr_overall.csv / _stratified.csv / _meta.csv / B7_m4_cell_scores.csv")
