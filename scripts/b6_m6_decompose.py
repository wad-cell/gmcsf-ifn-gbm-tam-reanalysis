# -*- coding: utf-8 -*-
"""B6 模块六: 拆解 68 基因集组成效应
目标: 解释 core68 与 hallmark IFN-a 在患者层面方向不一致的原因
1. core68 vs hallmark IFN-a 重叠
2. 分别打分比较 MO-TAM vs MG-TAM (患者层面)
3. 单基因层面: 每个基因在 MO-TAM vs MG-TAM 的表达差异
4. 找出驱动方向差异的基因
输出: core68_vs_hallmark_decomposition.csv / core68_vs_hallmark_report.md
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
core = pd.read_csv(os.path.join(WS, "output", "B4", "B4b_ISG_core_set.csv"))
core_syms = set(core["symbol"].dropna())
print("core68:", len(core_syms), "| hallmark:", len(HALLMARK_IFNA))
overlap = core_syms & set(HALLMARK_IFNA)
print("overlap:", len(overlap))
print("core68 not in hallmark:", sorted(core_syms - set(HALLMARK_IFNA)))
print("hallmark not in core68 (n):", len(set(HALLMARK_IFNA) - core_syms))

annot = pd.read_csv(os.path.join(TEMP, "annot.csv.gz"), compression="gzip")
with gzip.open(os.path.join(TEMP, "matrix.csv.gz"), "rt", encoding="utf-8", errors="replace") as f:
    mat = pd.read_csv(f, index_col=0)
mat = np.log1p(mat.astype(np.float32))
a = annot.set_index("cell")
a_sub = a.loc[list(mat.columns)]
a_sub["lineage"] = a_sub["cluster"].map(lambda cl: "MO-TAM" if cl=="TAM 1" else ("MG-TAM" if cl=="TAM 2" else "Other"))

# 单基因: MO-TAM vs MG-TAM 患者层面差异
mo_cells = a_sub[a_sub["lineage"]=="MO-TAM"].index
mg_cells = a_sub[a_sub["lineage"]=="MG-TAM"].index
all_genes = list(core_syms | set(HALLMARK_IFNA))
found = [g for g in all_genes if g in mat.index]
print("genes matched:", len(found), "/", len(all_genes))

rows = []
for g in found:
    expr = mat.loc[g]
    mo_pb = expr[mo_cells].groupby(a_sub.loc[mo_cells,"sample"]).mean()
    mg_pb = expr[mg_cells].groupby(a_sub.loc[mg_cells,"sample"]).mean()
    common = list(set(mo_pb.index) & set(mg_pb.index))
    if len(common) >= 3:
        d = (mo_pb[common] - mg_pb[common]).mean()
        w = stats.wilcoxon(mo_pb[common], mg_pb[common])
        rows.append(dict(gene=g, in_core68=g in core_syms, in_hallmark=g in set(HALLMARK_IFNA),
                         mo_minus_mg_patient=d, wilcoxon_p=w.pvalue, n_patients=len(common)))
gene_res = pd.DataFrame(rows)
gene_res.to_csv(os.path.join(OUT, "core68_vs_hallmark_decomposition.csv"), index=False)

# 汇总: core68-only vs hallmark-only vs overlap 基因的方向
def group(g):
    c = g in core_syms; h = g in set(HALLMARK_IFNA)
    if c and h: return "overlap"
    if c: return "core68_only"
    return "hallmark_only"
gene_res["group"] = gene_res["gene"].map(group)
print("\n=== 单基因患者层面 MO-MG 差异 (mean) ===")
print(gene_res.groupby("group")["mo_minus_mg_patient"].agg(["mean","median","count"]).round(4).to_string())
print("\n=== 各分组显著基因数 (wilcoxon p<0.05) ===")
print(gene_res[gene_res["wilcoxon_p"]<0.05].groupby("group")["gene"].count().to_string())

# 方向驱动基因
print("\n=== MO>MG 最强 (驱动MO高) top10 ===")
print(gene_res.sort_values("mo_minus_mg_patient", ascending=False).head(10)[["gene","group","mo_minus_mg_patient","wilcoxon_p"]].to_string())
print("\n=== MG>MO 最强 (驱动MG高) top10 ===")
print(gene_res.sort_values("mo_minus_mg_patient").head(10)[["gene","group","mo_minus_mg_patient","wilcoxon_p"]].to_string())

# 报告
md = f"""# B6 模块六: core68 vs hallmark IFN-a 组成拆解 (core68_vs_hallmark_report.md)

## 1. 基因集组成
- core68: {len(core_syms)} 基因 (GM-CSF 共培养下调的 ISG)
- hallmark IFN-a: {len(HALLMARK_IFNA)} 基因 (MSigDB)
- 重叠: {len(overlap)} 基因
- core68 独有: {len(core_syms - set(HALLMARK_IFNA))} 基因
- hallmark 独有: {len(set(HALLMARK_IFNA) - core_syms)} 基因

## 2. 单基因患者层面 MO-MG 差异 (mean log1p)
| 分组 | 基因数 | 平均MO-MG | 中位MO-MG |
|---|---|---|---|
"""
for grp in ["overlap","core68_only","hallmark_only"]:
    sub = gene_res[gene_res["group"]==grp]
    md += f"| {grp} | {len(sub)} | {sub['mo_minus_mg_patient'].mean():.4f} | {sub['mo_minus_mg_patient'].median():.4f} |\n"

md += f"""
## 3. 方向驱动基因
- MO>MG 最强 (驱动 MO 高, 与预期相反): {', '.join(gene_res.sort_values('mo_minus_mg_patient', ascending=False).head(8)['gene'])}
- MG>MO 最强 (驱动 MG 高, 支持预期): {', '.join(gene_res.sort_values('mo_minus_mg_patient').head(8)['gene'])}

## 4. 结论
- 若 core68_only 基因整体 MO>MG, 则 core68 的"非经典 ISG"成分驱动了与 hallmark 相反的方向
- 若 overlap 基因方向一致, 则差异来自基因集组成而非打分方法
"""
with open(os.path.join(OUT, "core68_vs_hallmark_report.md"), "w", encoding="utf-8") as f:
    f.write(md)
print("\nsaved core68_vs_hallmark_report.md")
