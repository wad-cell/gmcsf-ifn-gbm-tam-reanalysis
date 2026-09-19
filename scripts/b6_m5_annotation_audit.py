# -*- coding: utf-8 -*-
"""B6 模块五: 细胞来源注释可靠性检查 (GSM4972211)
1. 标志基因验证: TAM1=MO-TAM (S100A8/VCAN/ITGA4/TGFBI/FPR3), TAM2=MG-TAM (P2RY12/TMEM119/CX3CR1/SLC1A3)
2. Monocytes vs MO-TAM 区分
3. 患者间注释一致性 (每个患者各 lineage 的细胞数)
4. 各 cluster 标志基因表达谱
输出: scRNA_annotation_audit.csv / scRNA_annotation_audit.md
"""
import gzip, os
import numpy as np, pandas as pd

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
TEMP = os.path.join(WS, "temp", "GSM4972211")
OUT = os.path.join(WS, "output", "B6")
os.makedirs(OUT, exist_ok=True)

MO_MARKERS = ["S100A8", "S100A9", "VCAN", "ITGA4", "TGFBI", "FPR3", "CD14", "LYZ"]
MG_MARKERS = ["P2RY12", "TMEM119", "CX3CR1", "SLC1A3", "CSF1R", "TREM2", "C1QA", "C1QB"]
MONO_MARKERS = ["FCGR3A", "FCN1", "S100A8", "S100A9", "CD14"]
LYM_MARKERS = ["CD3G", "GZMA", "KLRB1", "NKG7", "CD79A"]

annot = pd.read_csv(os.path.join(TEMP, "annot.csv.gz"), compression="gzip")
with gzip.open(os.path.join(TEMP, "matrix.csv.gz"), "rt", encoding="utf-8", errors="replace") as f:
    mat = pd.read_csv(f, index_col=0)
mat = np.log1p(mat.astype(np.float32))

def lineage(cl):
    if cl == "TAM 1": return "MO-TAM"
    if cl == "TAM 2": return "MG-TAM"
    if cl == "prol. TAM": return "prol.TAM"
    if cl == "Monocytes": return "Monocytes"
    return "Other"
annot["lineage"] = annot["cluster"].map(lineage)

# 1. 标志基因在各 lineage 的平均表达
all_markers = list(dict.fromkeys(MO_MARKERS + MG_MARKERS + MONO_MARKERS + LYM_MARKERS))
found = [g for g in all_markers if g in mat.index]
print("markers matched:", len(found), "/", len(all_markers))
sub = mat.loc[found]
a = annot.set_index("cell")
# 按 lineage 平均 (确保列对齐)
sub_cols = list(sub.columns)
a_sub = a.loc[sub_cols]
lineage_mean = pd.DataFrame({ln: sub.loc[:, a_sub["lineage"] == ln].mean(axis=1) for ln in a_sub["lineage"].unique()})
print("\n=== Marker expression by lineage (mean log1p) ===")
print(lineage_mean.round(3).to_string())

# 2. 患者 x lineage 细胞数
ct = pd.crosstab(annot["sample"], annot["lineage"])
print("\n=== Patient x lineage cell counts ===")
print(ct.to_string())

# 3. 每个患者 MO-TAM / MG-TAM 比例
tot = ct.sum(axis=1)
print("\n=== MO-TAM fraction per patient ===")
print((ct["MO-TAM"] / tot).round(3).to_string())

# 4. 保存
lineage_mean.to_csv(os.path.join(OUT, "scRNA_annotation_audit.csv"))
ct.to_csv(os.path.join(OUT, "scRNA_patient_lineage_counts.csv"))

# 5. 报告
md = f"""# B6 模块五: 细胞来源注释可靠性审计 (scRNA_annotation_audit.md)

## 1. 方法
- 数据: GSM4972211, 21,303 细胞
- 注释来源: annot.csv.gz (Pombo Antunes 2021 原始注释)
- 标志基因验证: MO-TAM (S100A8/VCAN/ITGA4/TGFBI/FPR3), MG-TAM (P2RY12/TMEM119/CX3CR1/SLC1A3)

## 2. 标志基因表达 (mean log1p)
| 基因 | MO-TAM | MG-TAM | Monocytes | 判定 |
|---|---|---|---|---|
"""
for g in found:
    row = lineage_mean.loc[g]
    mo, mg = row.get("MO-TAM", np.nan), row.get("MG-TAM", np.nan)
    verdict = "MO高" if mo > mg else ("MG高" if mg > mo else "持平")
    md += f"| {g} | {mo:.3f} | {mg:.3f} | {row.get('Monocytes', np.nan):.3f} | {verdict} |\n"

md += f"""
## 3. 患者 x lineage 细胞数
| 患者 | MO-TAM | MG-TAM | Monocytes | 其他 | MO-TAM比例 |
|---|---|---|---|---|---|
"""
for s in ct.index:
    md += f"| {s} | {ct.loc[s,'MO-TAM']} | {ct.loc[s,'MG-TAM']} | {ct.loc[s].get('Monocytes',0)} | {ct.loc[s].sum()-ct.loc[s,'MO-TAM']-ct.loc[s,'MG-TAM']-ct.loc[s].get('Monocytes',0)} | {(ct.loc[s,'MO-TAM']/ct.loc[s].sum()):.3f} |\n"

md += """
## 4. 结论
- 若 MO 标志基因在 TAM1 高、MG 标志基因在 TAM2 高 → 注释可靠
- 若存在混合/异常 → 注释不可靠, 需谨慎解读
- 注意: Monocytes 与 MO-TAM 的区分 (S100A8/9 均高, 需看成熟标志)
"""
with open(os.path.join(OUT, "scRNA_annotation_audit.md"), "w", encoding="utf-8") as f:
    f.write(md)
print("\nsaved scRNA_annotation_audit.md")
