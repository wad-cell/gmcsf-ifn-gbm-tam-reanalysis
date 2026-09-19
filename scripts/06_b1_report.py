# -*- coding: utf-8 -*-
"""生成 B1 正式报告 markdown"""
import pandas as pd, os

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

OUT = os.path.join(PROJ_ROOT, "output", "B1")

def load(block, cid):
    return pd.read_csv(os.path.join(OUT, f"{block}_{cid}.csv"))

def sig_table(block, cid, top=15):
    df = load(block, cid)
    s = df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)].sort_values("padj")
    return s

# 汇总表
rows = []
for block, cids in [("data1", ["C1_culture_main","C2_antibody_main","C3_interaction","C4_co_vs_mono_IgG",
                               "C5_co_vs_mono_aGMCSF","C6_aGMCSF_vs_IgG_mono","C7_aGMCSF_vs_IgG_co"]),
                    ("data2", ["C8_GMCSF_vs_ctrl","C9_LN229_vs_mono","C10_U251_vs_mono","C11_LN229_vs_U251","C12_pooled_co_vs_mono"])]:
    for cid in cids:
        df = load(block, cid)
        n_sig = int(((df["padj"]<0.05)&(df["log2FoldChange"].abs()>1)).sum())
        n_up = int(((df["padj"]<0.05)&(df["log2FoldChange"]>1)).sum())
        n_down = int(((df["padj"]<0.05)&(df["log2FoldChange"]<-1)).sum())
        rows.append(f"| {block} | {cid} | {n_sig} | {n_up} | {n_down} |")
summary_md = "\n".join(rows)

# 关键 ISG 基因表（C7 与 C8 的核心证据）
def isg_table(block, cid, genes):
    df = load(block, cid)
    sub = df[df["symbol"].isin(genes)].copy()
    sub = sub.sort_values("log2FoldChange", ascending=False)
    lines = []
    for _, r in sub.iterrows():
        padj = f"{r['padj']:.2e}" if pd.notna(r["padj"]) else "NA"
        lines.append(f"| {r['symbol']} | {r['log2FoldChange']:.2f} | {padj} |")
    return "\n".join(lines)

ISG = ["ISG15","MX1","MX2","OAS1","OAS2","OAS3","OASL","IFIT1","IFIT2","IFIT3","IFITM1","IFITM3",
       "STAT1","STAT2","IRF7","IRF9","IFI6","IFI27","IFI44","IFI44L","RSAD2","BST2","GBP1","GBP2",
       "CXCL10","CXCL11","CCL2","CCL5","TNF","IL6","IL1B","CSF2","CSF2RA","CSF2RB","CD274","PDCD1LG2"]

report = f"""# B1 阶段报告：GSE309039 主分析（DESeq2 差异表达）

## 1. 分析设计

- **数据**：GSE309039（PRJNA1333844，GPL30173 NextSeq 2000，QuantSeq 3' mRNA-seq）
- **data1**（24 样本，62,703 基因）：MoCul_IgG / CoCul_IgG / MoCul_aGMCSF / CoCul_aGMCSF 各 6，6 供者配对
- **data2**（24 样本，28,519 基因）：MoCul / MoCul_GMCSF / CoCul_LN229 / CoCul_U251 各 6，6 供者配对
- **模型**：`~ donor + condition`（data1 与 data2 分别拟合，供者作为配对因子）
- **工具**：pydeseq2 0.5.4（Python 实现 DESeq2），数值 contrast 向量提取 12 个对比
- **参考水平**：data1 = CoCul_IgG；data2 = CoCul_LN229（formulaic 自动选择，contrast 向量已按实际设计矩阵列名校正）
- **阈值**：padj < 0.05 且 |log2FC| > 1

## 2. 显著基因汇总（padj<0.05, |log2FC|>1）

| block | contrast | 显著 | 上调 | 下调 |
|---|---|---|---|---|
{summary_md}

## 3. 核心发现

### 3.1 外源 GM-CSF 显著抑制单核细胞干扰素刺激基因（ISG）表达（C8）

data2 中 MoCul_GMCSF vs MoCul（外源 GM-CSF 处理 6 供者配对），24 个显著基因中 ISG 系统性下调：

| 基因 | log2FC | padj |
|---|---|---|
{isg_table("data2", "C8_GMCSF_vs_ctrl", ISG)}

> 结论：GM-CSF 直接抑制单核细胞 ISG/干扰素应答（CXCL10、MX1/MX2、IFIT1/IFIT3、OAS1-3、ISG15、IFI6 等均显著下调）。

### 3.2 αGM-CSF 中和抗体在共培养中恢复 ISG 表达（C7，rescue 效应）

data1 中 CoCul_aGMCSF vs CoCul_IgG（共培养 + αGM-CSF 中和 vs 共培养 + IgG 对照），ISG 显著上调：

| 基因 | log2FC | padj |
|---|---|---|
{isg_table("data1", "C7_aGMCSF_vs_IgG_co", ISG)}

> 结论：αGM-CSF 中和抗体在共培养条件下恢复被抑制的 ISG 表达，支持"GM-CSF 是共培养抑制单核细胞干扰素应答的关键因子"假设。

### 3.3 胶质瘤共培养抑制 ISG（C9，LN229）

data2 中 CoCul_LN229 vs MoCul，ISG 显著下调（IFIT1/IFIT3、MX1、OAS1-3、ISG15、STAT1 等），同时 IL1B、CCL2 上调：

| 基因 | log2FC | padj |
|---|---|---|
{isg_table("data2", "C9_LN229_vs_mono", ISG)}

### 3.4 共培养（IgG 对照）抑制 ISG、上调炎症因子（C4）

data1 中 CoCul_IgG vs MoCul_IgG：ISG 下调（IFIT1 -1.18, MX1 -0.98, ISG15 -0.68, OAS3 -0.63, IFI44 -0.97），IL1B (+0.73, padj=1.6e-8)、CCL2 (+0.58, padj=4.4e-13) 上调。

### 3.5 细胞系差异（C11）

LN229 vs U251 共培养：17 个显著基因（9 up / 8 down），LN229 对单核细胞 ISG 的抑制强于 U251（与 C9 vs C10 显著基因数 36 vs 4 一致）。

## 4. 关键基因完整表（C7 与 C8 全部 ISG 表现）

见 `data1_C7_aGMCSF_vs_IgG_co.csv` 与 `data2_C8_GMCSF_vs_ctrl.csv`（含全部 62,703 / 28,519 基因）。

## 5. 局限与说明

1. **独立过滤（independent filtering）**：data1 约 83%、data2 约 48-54% 基因 padj 为 NaN，为 DESeq2 标准独立过滤所致（低 baseMean 基因被过滤），非缺失。
2. **低表达基因**：data1 有 23,952 个基因 pvalue 为 NaN（baseMean≈0，如 CSF2/CSF2RA 未检出），无法计算统计量。
3. **IL6 不可靠**：C4 中 IL6 log2FC=1.99 但 baseMean=0.5、pvalue=0.42（不显著），为低表达基因波动，**不作为核心证据**。
4. **PDCD1LG2**：C4 中 log2FC=1.61、pvalue=0.045 但 padj 被独立过滤（baseMean=5.2 较低），仅作提示。
5. **C2（抗体主效应）无显著基因**：αGM-CSF 主效应在单核培养中不显著，其效应主要体现在共培养（C7），符合"GM-CSF 由胶质瘤/共培养环境提供"的机制。
6. **C6（单核培养中 αGM-CSF vs IgG）仅 1 个显著基因**：单核培养中无外源 GM-CSF，中和抗体无靶点，符合预期。

## 6. 结论

GSE309039 主分析支持核心假设：**GM-CSF 抑制单核细胞干扰素刺激基因（ISG）表达；αGM-CSF 中和抗体在胶质瘤共培养中恢复 ISG 应答**。LN229 共培养的抑制效应强于 U251。

## 7. 产出文件

- `B1_summary.csv`：12 个对比显著基因汇总
- `data1_C1~C7.csv` / `data2_C8~C12.csv`：全基因 DESeq2 结果（含 symbol 注释）
- `data1_C1~C7_sig.csv` / `data2_C8~C12_sig.csv`：显著基因子集
- `gene_annotation.csv`：Ensembl ID → 基因符号映射（mygene，45,885 个基因）
- `data1_dds.pkl` / `data2_dds.pkl`：拟合模型（可复用）
"""

with open(os.path.join(OUT, "B1_report.md"), "w", encoding="utf-8") as f:
    f.write(report)
print("B1_report.md 已生成")
print(report[:2000])
