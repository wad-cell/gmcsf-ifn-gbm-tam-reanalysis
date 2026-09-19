# -*- coding: utf-8 -*-
"""B6 模块四: 复现原论文 scRNA 分析 (GSM4972211)
原论文 (Meyer et al., EMM 2026) 方法:
- 用 UCell v2.6.2 (Mann-Whitney U 统计) 计算 MSigDB "hallmark interferon-alpha response" 得分
- 比较 MO-TAM vs MG-TAM
本脚本:
1. 严格 UCell 公式复现 hallmark IFN-a 得分
2. 细胞层面 + 患者层面 (pseudobulk) 比较 MO-TAM vs MG-TAM
3. 与 B5 core68 打分结果对比
4. 患者间变异分析 (原论文 Supplementary Fig 1d 提示变异大)
输出: scRNA_reproduction_UCell.csv / scRNA_reproduction_report.md
"""
import gzip, os
import numpy as np, pandas as pd
from scipy import stats

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
TEMP = os.path.join(WS, "temp", "GSM4972211")
OUT = os.path.join(WS, "output", "B6")
os.makedirs(OUT, exist_ok=True)

HALLMARK_IFNA = ("MX1 ISG15 OAS1 IFIT3 IFI44 IFI35 IRF7 RSAD2 IFI44L IFITM1 IFI27 IRF9 OASL "
"EIF2AK2 IFIT2 CXCL10 TAP1 SP110 DDX60 UBE2L6 USP18 PSMB8 IFIH1 BST2 LGALS3BP ADAR ISG20 "
"GBP2 IRF1 PLSCR1 PSMB9 HERC6 SAMD9 CMPK2 IFITM3 RTP4 STAT2 SAMD9L LY6E IFITM2 HELZ2 CXCL11 "
"TRIM21 PARP14 TRIM26 PARP12 NMI RNF31 HLA-C CASP1 TRIM14 TDRD7 DHX58 PARP9 PNPT1 TRIM25 PSME1 "
"WARS1 EPSTI1 UBA7 PSME2 B2M TRIM5 C1S LAP3 LAMP3 GBP4 NCOA7 TMEM140 CD74 GMPR PSMA3 PROCR IL7 "
"IFI30 IRF2 CSF1 IL15 CNP TENT5A IL4R CMTR1 CD47 LPAR6 MOV10 CASP8 TXNIP SLC25A28 SELL TRAFD1 "
"BATF2 RIPK2 CCRL2 NUB1 OGFR MVB12A ELF1").split()
print("hallmark genes:", len(HALLMARK_IFNA))

annot = pd.read_csv(os.path.join(TEMP, "annot.csv.gz"), compression="gzip")
with gzip.open(os.path.join(TEMP, "matrix.csv.gz"), "rt", encoding="utf-8", errors="replace") as f:
    mat = pd.read_csv(f, index_col=0)
print("mat:", mat.shape)
mat = np.log1p(mat.astype(np.float32))

# 严格 UCell: Mann-Whitney U 统计
# 对每个细胞, 基因按表达 rank (1=最低), U = sum(rank_S) - n_S*(n_S+1)/2, score = U/(n_S*n_total)
found = [g for g in HALLMARK_IFNA if g in mat.index]
print("hallmark matched:", len(found), "/", len(HALLMARK_IFNA))
ranks = mat.rank(axis=0, method="average")  # 每列(细胞)内基因 rank
n_total = mat.shape[0]
n_S = len(found)
U = ranks.loc[found].sum(axis=0) - n_S * (n_S + 1) / 2
ucell = U / (n_S * n_total)
print("UCell score range:", ucell.min().round(4), ucell.max().round(4))

cells = ucell.index
a = annot.set_index("cell").loc[cells]
res = pd.DataFrame({
    "cell": cells,
    "ucell_ifna": ucell.values,
    "cluster": a["cluster"].values,
    "sample": a["sample"].values,
})
def lineage(cl):
    if cl == "TAM 1": return "MO-TAM"
    if cl == "TAM 2": return "MG-TAM"
    if cl == "prol. TAM": return "prol.TAM"
    if cl == "Monocytes": return "Monocytes"
    return "Other"
res["lineage"] = res["cluster"].map(lineage)
res.to_csv(os.path.join(OUT, "scRNA_reproduction_UCell.csv"), index=False)

# 细胞层面: MO-TAM vs MG-TAM
mo_c = res[res["lineage"]=="MO-TAM"]["ucell_ifna"]
mg_c = res[res["lineage"]=="MG-TAM"]["ucell_ifna"]
print(f"\nCell-level: MO-TAM n={len(mo_c)} mean={mo_c.mean():.4f} | MG-TAM n={len(mg_c)} mean={mg_c.mean():.4f}")
u_stat, p_cell = stats.mannwhitneyu(mo_c, mg_c, alternative="two-sided")
print(f"Mann-Whitney U={u_stat:.0f}, p={p_cell:.3e}")

# 患者层面 pseudobulk
pb = res.groupby(["sample","lineage"])["ucell_ifna"].mean().reset_index()
pb["n_cells"] = res.groupby(["sample","lineage"]).size().values
piv = pb.pivot(index="sample", columns="lineage", values="ucell_ifna")
print("\nPatient-level UCell (hallmark IFNa):")
print(piv.round(4).to_string())

mo = pb[pb["lineage"]=="MO-TAM"].set_index("sample")["ucell_ifna"]
mg = pb[pb["lineage"]=="MG-TAM"].set_index("sample")["ucell_ifna"]
common_s = list(set(mo.index) & set(mg.index))
print(f"\nPatients with both: {len(common_s)}")
if len(common_s) >= 3:
    w = stats.wilcoxon(mo[common_s], mg[common_s])
    t = stats.ttest_rel(mo[common_s], mg[common_s])
    d = (mo[common_s] - mg[common_s]).mean()
    print(f"Wilcoxon: W={w.statistic}, p={w.pvalue:.4f}")
    print(f"paired t: t={t.statistic:.3f}, p={t.pvalue:.4f}")
    print(f"mean diff (MO-MG)={d:.4f}")
    # 方向一致性
    n_opp = int((mo[common_s] > mg[common_s]).sum())
    print(f"patients where MO>MG (opposite to paper): {n_opp}/{len(common_s)}")

# 患者间变异
print("\nPatient-level MG-TAM SD:", mg.std().round(4), "| MO-TAM SD:", mo.std().round(4))
print("MG-TAM range:", mg.min().round(4), "-", mg.max().round(4))
print("MO-TAM range:", mo.min().round(4), "-", mo.max().round(4))

# 与 B5 core68 结果对比
b5 = pd.read_csv(os.path.join(WS, "output", "B5", "external_scRNA_patient_level_validation.csv"))
print("\nB5 core68 patient-level (module_score):")
print(b5.pivot(index="sample", columns="lineage", values="module_score").round(4).to_string())

# 报告
md = f"""# B6 模块四: 复现原论文 scRNA 分析 (scRNA_reproduction_report.md)

## 1. 复现方法
- 数据: GSM4972211 (GSE163120), 21,303 细胞, ND1-ND7
- 打分: 严格 UCell (Mann-Whitney U 统计), 基因集 = MSigDB HALLMARK_INTERFERON_ALPHA_RESPONSE ({len(HALLMARK_IFNA)} 基因, 矩阵匹配 {len(found)} 个)
- 分组: annot 原始注释 TAM1=MO-TAM, TAM2=MG-TAM (B5 已正交验证)
- 原论文方法: UCell v2.6.2 + hallmark IFN-a (MOESM1 确认)

## 2. 细胞层面结果
- MO-TAM: n={len(mo_c)}, mean={mo_c.mean():.4f}
- MG-TAM: n={len(mg_c)}, mean={mg_c.mean():.4f}
- Mann-Whitney U={u_stat:.0f}, p={p_cell:.3e}
- 方向: {'MO-TAM 低于 MG-TAM (支持原论文)' if mo_c.mean()<mg_c.mean() else 'MO-TAM 高于 MG-TAM (不支持原论文)'}

## 3. 患者层面结果 (pseudobulk, 禁止细胞独立重复)
- 配对患者数: {len(common_s)}
- Wilcoxon signed-rank: W={w.statistic}, p={w.pvalue:.4f}
- paired t: t={t.statistic:.3f}, p={t.pvalue:.4f}
- mean diff (MO-MG)={d:.4f}
- 方向相反患者数: {n_opp}/{len(common_s)} (MO>MG)
- 患者间变异: MG-TAM SD={mg.std():.4f} (range {mg.min():.4f}-{mg.max():.4f}), MO-TAM SD={mo.std():.4f} (range {mo.min():.4f}-{mo.max():.4f})

## 4. 与 B5 core68 结果对比
- B5 用 core68 模块评分, 本模块用 hallmark IFN-a (原论文基因集)
- 两者打分方法不同, 需对比方向一致性

## 5. 结论
- 细胞层面: {'支持' if mo_c.mean()<mg_c.mean() else '不支持'}原论文 (MO-TAM 低 ISG)
- 患者层面: {'支持' if d<0 else '不支持'}原论文
- 患者间变异大, 需模块五核查注释可靠性 + 模块六拆解基因集组成
"""
with open(os.path.join(OUT, "scRNA_reproduction_report.md"), "w", encoding="utf-8") as f:
    f.write(md)
print("\nsaved scRNA_reproduction_report.md")
