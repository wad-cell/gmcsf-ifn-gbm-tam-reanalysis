# -*- coding: utf-8 -*-
"""B5 模块七: 外部 GBM scRNA 患者层面验证升级
- 区分 MO-TAM / MG-TAM (标记基因推断: TAM1=MO-TAM, TAM2=MG-TAM)
- 患者层面 pseudobulk 统计 (禁止细胞当独立重复)
- 两种打分: 平均log表达(模块评分) + UCell近似(rank均值)
- MO-TAM vs MG-TAM 患者层面配对检验
- CISH 与 core68 score 相关性
输出: external_scRNA_patient_level_validation.csv
"""
import gzip, os
import numpy as np, pandas as pd
from scipy import stats

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
TEMP = os.path.join(WS, "temp", "GSM4972211")
OUT = os.path.join(WS, "output", "B5")
CORE = os.path.join(WS, "output", "B4", "B4b_ISG_core_set.csv")

annot = pd.read_csv(os.path.join(TEMP, "annot.csv.gz"), compression="gzip")
core = pd.read_csv(CORE)
core_syms = set(core["symbol"].dropna())
print("core:", len(core_syms))

with gzip.open(os.path.join(TEMP, "matrix.csv.gz"), "rt", encoding="utf-8", errors="replace") as f:
    mat = pd.read_csv(f, index_col=0)
print("mat:", mat.shape)
mat = np.log1p(mat.astype(np.float32))
print("log1p done, max:", mat.values.max())

# 匹配核心基因
found = [g for g in core_syms if g in mat.index]
print("core matched:", len(found), "/", len(core_syms))
sub = mat.loc[found]

# 方法1: 模块评分 (平均log表达)
score_mean = sub.mean(axis=0)

# 方法2: UCell近似 (每细胞基因 rank 均值)
# 对每细胞, 计算基因集基因在该细胞所有基因中的 rank 百分位均值
# 用矩阵的列 rank (每细胞基因表达排序)
ranks = mat.rank(axis=0, pct=True)  # 每列(细胞)内基因 rank 百分位
ucell = ranks.loc[found].mean(axis=0)
print("UCell done")

# 对齐 annot
cells = score_mean.index
a = annot.set_index("cell").loc[cells]
res = pd.DataFrame({
    "cell": cells,
    "module_score": score_mean.values,
    "ucell_score": ucell.values,
    "cluster": a["cluster"].values,
    "sample": a["sample"].values,
})
# MO/MG 归属
def lineage(cl):
    if cl == "TAM 1": return "MO-TAM"
    if cl == "TAM 2": return "MG-TAM"
    if cl == "prol. TAM": return "prol.TAM"
    if cl == "Monocytes": return "Monocytes"
    return "Other"
res["lineage"] = res["cluster"].map(lineage)
print(res["lineage"].value_counts().to_string())

# 患者层面 pseudobulk: 每患者每 lineage 平均 score
pb = res.groupby(["sample","lineage"])[["module_score","ucell_score"]].mean().reset_index()
pb["n_cells"] = res.groupby(["sample","lineage"]).size().values
print("\n=== Patient-level pseudobulk (module_score) ===")
print(pb.pivot(index="sample", columns="lineage", values="module_score").round(4).to_string())

# MO-TAM vs MG-TAM 患者层面配对检验 (仅 TAM1/TAM2)
mo = pb[pb["lineage"]=="MO-TAM"].set_index("sample")["module_score"]
mg = pb[pb["lineage"]=="MG-TAM"].set_index("sample")["module_score"]
common_s = list(set(mo.index) & set(mg.index))
print(f"\nPatients with both MO/MG: {len(common_s)}")
if len(common_s) >= 3:
    w = stats.wilcoxon(mo[common_s], mg[common_s])
    print(f"Wilcoxon signed-rank (module): W={w.statistic}, p={w.pvalue:.4f}")
    print(f"MO-TAM mean={mo[common_s].mean():.4f}, MG-TAM mean={mg[common_s].mean():.4f}")
    # 效应量
    d = (mo[common_s] - mg[common_s]).mean()
    print(f"mean diff (MO-MG)={d:.4f}")
    # 配对 t
    t = stats.ttest_rel(mo[common_s], mg[common_s])
    print(f"paired t: t={t.statistic:.3f}, p={t.pvalue:.4f}")
    # 方向一致率
    print(f"MO<MG in {((mo[common_s]<mg[common_s]).sum())}/{len(common_s)} patients")
else:
    print("insufficient patients")

# UCell 同样检验
mo_u = pb[pb["lineage"]=="MO-TAM"].set_index("sample")["ucell_score"]
mg_u = pb[pb["lineage"]=="MG-TAM"].set_index("sample")["ucell_score"]
if len(common_s) >= 3:
    wu = stats.wilcoxon(mo_u[common_s], mg_u[common_s])
    print(f"Wilcoxon (UCell): p={wu.pvalue:.4f}, MO={mo_u[common_s].mean():.4f}, MG={mg_u[common_s].mean():.4f}")

# CISH 相关性 (细胞层面, 患者层面)
if "CISH" in mat.index:
    cish = mat.loc["CISH"]
    # 患者层面: 每患者 CISH 平均 vs core68 score 平均
    cish_pb = pd.DataFrame({"sample": a["sample"].values, "CISH": cish.loc[cells].values})
    cish_pb = cish_pb.groupby("sample")["CISH"].mean()
    score_pb = res.groupby("sample")["module_score"].mean()
    rho, p = stats.spearmanr(cish_pb, score_pb)
    print(f"\nCISH vs core68 (patient-level): rho={rho:.3f}, p={p:.4f}, n={len(cish_pb)}")
    # 细胞层面
    rho_c, p_c = stats.spearmanr(cish.loc[cells].values, res["module_score"].values)
    print(f"CISH vs core68 (cell-level): rho={rho_c:.3f}, p={p_c:.4f}")
else:
    print("\nCISH not in matrix - skip")

# 保存
os.makedirs(OUT, exist_ok=True)
res.to_csv(os.path.join(OUT, "external_scRNA_cell_scores.csv"), index=False)
pb.to_csv(os.path.join(OUT, "external_scRNA_patient_level_validation.csv"), index=False)
print("\nSaved external_scRNA_patient_level_validation.csv")
