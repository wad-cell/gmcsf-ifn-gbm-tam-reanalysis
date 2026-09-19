# -*- coding: utf-8 -*-
"""B3 正式报告生成：GSE309038 细胞系表达谱"""
import pandas as pd, os

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

OUT_B3 = os.path.join(PROJ_ROOT, "output", "B3")

def load(cid):
    return pd.read_csv(os.path.join(OUT_B3, f"GSE309038_{cid}.csv"))

t98g = load("B3_T98G_vs_LN229")
u87 = load("B3_U87_vs_LN229")
u251 = load("B3_U251_vs_LN229")

def sig(df, p=0.05, lfc=1.0):
    s = df[(df["padj"]<p)&(df["log2FoldChange"].abs()>lfc)].dropna(subset=["symbol"])
    return s

s_t98g, s_u87, s_u251 = sig(t98g), sig(u87), sig(u251)
print(f"显著: T98G {len(s_t98g)}, U87 {len(s_u87)}, U251 {len(s_u251)}")

# 关键基因表达矩阵
cell_means = pd.read_csv(os.path.join(OUT_B3, "GSE309038_key_genes_cell_mean_log2.csv"), index_col=0)
key_order = ["CSF2","CSF1","IL6","IL1B","CCL2","CXCL10","TGFB1","TGFB2","TGFB3",
             "IFNAR1","IFNAR2","TGFBR1","TGFBR2","TGFBR3","SMAD2","SMAD3","SMAD4",
             "STAT1","STAT2","IRF1","IRF7","MIF"]
cell_means = cell_means.reindex([g for g in key_order if g in cell_means.index])

def md_table(df, cols, n=30):
    lines = []
    for _, r in df.head(n).iterrows():
        vals = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                if abs(v) > 100: vals.append(f"{v:.2e}")
                else: vals.append(f"{v:.2f}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)

def md_table_index(df, cols, n=30):
    lines = []
    for idx, r in df.head(n).iterrows():
        vals = [str(idx)]
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                if abs(v) > 100: vals.append(f"{v:.2e}")
                else: vals.append(f"{v:.2f}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)

# U251 vs LN229 关键基因差异
u251_key = u251[u251["symbol"].isin(key_order)][["symbol","baseMean","log2FoldChange","padj"]].dropna(subset=["symbol"])
u251_key = u251_key.sort_values("padj")

report = f"""# B3 阶段报告：GSE309038 细胞系表达谱（分泌因子归因）

## 1. 分析设计

- **数据**：GSE309038（12 样本，24,795 基因，4 细胞系 × 3 独立培养，非供者配对）
- **分组**：LN229 / T98G / U87 / U251（各 n=3，48h 未处理培养）
- **平台**：QuantSeq 3' mRNA-seq（GPL30173），GRCh38，htseq-count
- **模型**：`~ condition`（pydeseq2 0.5.4），参考水平 LN229
- **对比**：T98G vs LN229、U87 vs LN229、U251 vs LN229
- **目的**：提取细胞系分泌因子（GM-CSF/TGF-β 等）表达，解释 LN229 vs U251 对单核细胞 ISG 抑制强度差异（B1 C9 vs C10）

## 2. 显著基因汇总（padj<0.05, |log2FC|>1）

| 对比 | 显著 | 上调 | 下调 |
|---|---|---|---|
| T98G vs LN229 | {len(s_t98g)} | {(s_t98g['log2FoldChange']>0).sum()} | {(s_t98g['log2FoldChange']<0).sum()} |
| U87 vs LN229 | {len(s_u87)} | {(s_u87['log2FoldChange']>0).sum()} | {(s_u87['log2FoldChange']<0).sum()} |
| U251 vs LN229 | {len(s_u251)} | {(s_u251['log2FoldChange']>0).sum()} | {(s_u251['log2FoldChange']<0).sum()} |

## 3. 关键基因各细胞系平均表达（log2(norm+1)）

| 基因 | LN229 | T98G | U87 | U251 |
|---|---|---|---|---|
{md_table_index(cell_means, ["LN229","T98G","U87","U251"], 30)}

> 注：CSF2（GM-CSF 编码基因）在 GSE309038 转录组中**无表达信号**（不在 24,795 基因索引中），GSE309039 data1/data2 与 GSE309037 中亦为 0 或缺失。**转录组无信号 ≠ 蛋白不分泌**：原论文（Exp Mol Med 2026, PMID 42399654）用 ELISA 检测到 GM-CSF 蛋白在 LN229/T98G/U87 上清中分泌，而 U251 上清未检出（Fig 3c, n=4）。QuantSeq 3' 低深度（5-9M reads）+ 分泌型细胞因子转录本极低，可解释转录组缺失。

## 4. 核心发现

### 4.1 GM-CSF 分泌差异（论文蛋白证据 vs 本转录组）

- 论文 Fig 3c：GM-CSF 蛋白在 **LN229、T98G、U87** 上清检出，**U251 未检出**（ELISA, n=4）。
- 论文 Fig 3b：GM-CSF 靶基因 CISH 在单核细胞共培养中显著升高（LN229/T98G/U87），**U251 共培养不升高**（qPCR, n=7）。
- 本转录组：CSF2 无信号，无法从 mRNA 层面直接验证；但 CISH 在 U251 细胞系自身表达反而更高（log2FC +2.92, padj=0.010），提示 CISH 在细胞系中的表达与单核细胞 GM-CSF 活性标志意义不同，**不能**作为细胞系分泌 GM-CSF 的替代证据。

### 4.2 TGF-β 通路（论文核心机制：GM-CSF → TGF-β → ISG 抑制）

U251 vs LN229 中 TGF-β 通路核心转录因子显著下调：

| 基因 | baseMean | log2FC | padj |
|---|---|---|---|
{md_table(u251_key[u251_key['symbol'].isin(['SMAD2','SMAD3','SMAD4','TGFBR2','TGFB2','TGFB1','TGFB3'])], ["symbol","baseMean","log2FoldChange","padj"], 10)}

- SMAD2/3/4 在 U251 中显著低于 LN229（SMAD4 -2.82, SMAD2 -1.68, SMAD3 -1.56, 均 padj<0.001）。
- TGFB2 在 U251 中反而更高（+2.17, padj<0.001），TGFB1 无差异——与论文 Fig 5c 一致（TGF-β 蛋白在各细胞系上清无显著差异）。
- 提示：U251 的"无 ISG 抑制表型"更可能源于**缺乏 GM-CSF 分泌**（而非 TGF-β 配体缺失），与论文结论一致。

### 4.3 其他分泌因子

- CSF1（M-CSF）：U251 显著高于 LN229（+3.70, padj<0.001）；T98G/U87 亦高。
- MIF：U251 几乎不表达（log2FC -11.0, padj<0.001），LN229/T98G/U87 高表达。
- IL6：U251 低于 LN229（-1.39, 不显著）；IL1B 无显著差异。
- IFNAR1/2：U251 中 IFNAR2 略低（-0.83, padj=0.017），IFNAR1 无差异——与论文"IFNAR 可用性未改变"结论一致。

## 5. 结论

GSE309038 细胞系转录组支持：**LN229/T98G/U87 与 U251 在分泌因子谱上存在系统性差异**，其中 GM-CSF 分泌（蛋白水平，论文 ELISA 证据）是区分"ISG 抑制型"（LN229/T98G/U87）与"非抑制型"（U251）细胞系的关键。本转录组虽无法直接检测 CSF2 mRNA，但 TGF-β 通路（SMAD2/3/4）在 U251 中的下调与论文"GM-CSF 经 TGF-β 通路抑制 ISG"的机制框架一致。**注意**：细胞系转录组为纯生信证据，GM-CSF 分泌结论依赖原论文 ELISA 数据（待核验原文 Fig 3c 数值）。

## 6. 产出文件

- `GSE309038_B3_T98G_vs_LN229.csv` / `_U87_vs_LN229.csv` / `_U251_vs_LN229.csv`：全基因 DESeq2 结果
- `GSE309038_key_genes_norm_counts.csv`：关键基因各样本归一化计数
- `GSE309038_key_genes_cell_mean_log2.csv`：关键基因各细胞系平均表达
- `GSE309038_all_genes_cell_mean_log2.csv`：全基因各细胞系平均表达（供 B4c 关联）
- `GSE309038_dds.pkl`：拟合模型
"""

with open(os.path.join(OUT_B3, "B3_report.md"), "w", encoding="utf-8") as f:
    f.write(report)
print("\nB3_report.md 已生成")
