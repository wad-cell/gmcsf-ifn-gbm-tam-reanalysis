# -*- coding: utf-8 -*-
"""B4 正式报告生成"""
import pandas as pd, os

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

OUT = os.path.join(PROJ_ROOT, "output", "B4")

core = pd.read_csv(f"{OUT}\\B4b_ISG_core_set.csv")
c11 = pd.read_csv(f"{OUT}\\B4c_C11_sig_with_ISG.csv")

def md_table(df, cols, n=25):
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

report = f"""# B4 阶段报告：跨数据集整合（效应量关联 + ISG 核心集 + 分泌因子归因）

## 1. 分析设计

- **B4a 直接 vs 间接共培养效应量关联**：GSE309037 `CoCul vs MoCul`（直接 48h）显著基因，在 GSE309039 data1 `C4 CoCul_IgG vs MoCul_IgG`（间接 24h+IgG）与 data2 `C9 CoCul_LN229 vs MoCul`（间接 24h）中的效应量一致性。
- **B4b ISG 核心集**：IFN-β 金标准 ISG（760）∩ 直接共培养显著下调（338）→ 共培养抑制型 ISG 核心集。
- **B4c 分泌因子-响应关联**：GSE309038 细胞系分泌因子表达 vs GSE309039 `C11 LN229 vs U251 共培养` 单核细胞差异响应。
- **红线**：跨数据集仅做基因集交集/效应量关联/签名投影，不做合并统计（不同批次/时间点）。

## 2. B4a 直接 vs 间接共培养效应量一致性

以 GSE309037 直接共培养显著基因（849 个）为锚，观察其在间接共培养中的效应：

| 锚定基因 | 间接对比 | 中位 log2FC | 同向比例 | Spearman rho (p) |
|---|---|---|---|---|
| CoCul 下调 (340) | C4 间接(IgG) | -0.23 | 77% | — |
| CoCul 上调 (509) | C4 间接(IgG) | +0.18 | 74% | — |
| 全部显著 (849) | C4 间接(IgG) | — | — | 0.504 (5.8e-56) |
| CoCul 下调 (340) | C9 间接(LN229) | -0.23 | 81% | — |
| CoCul 上调 (509) | C9 间接(LN229) | +0.24 | 78% | — |
| 全部显著 (849) | C9 间接(LN229) | — | — | 0.599 (8.0e-84) |

**结论**：直接共培养（48h 接触）与间接共培养（24h transwell）在显著基因层面效应方向高度一致（74-81% 同向，Spearman rho 0.50-0.60），支持**可溶性因子介导**的 ISG 抑制机制（与论文 Fig 2c 结论一致）。注意：全基因关联 rho 仅 0.03-0.10（噪声稀释），显著基因子集关联才有生物学意义。

## 3. B4b 共培养抑制型 ISG 核心集（68 基因）

IFN-β 金标准 ISG（760）∩ 直接共培养显著下调（338）= **68 个核心基因**。这些基因在 IFN-β 诱导下上调、在共培养中被抑制，是"GB 细胞抑制单核细胞 I 型 IFN 应答"的直接分子证据。

Top 25（按 IFN 诱导强度排序）：

| 基因 | IFN_log2FC | CoCul_log2FC | C4_log2FC | C9_log2FC |
|---|---|---|---|---|
{md_table(core, ["symbol","IFN_log2FC","CoCul_log2FC","C4_log2FC","C9_log2FC"], 25)}

> 注：C4/C9 为间接共培养中的 log2FC（负值=同样被抑制）。核心集基因在直接与间接共培养中均呈抑制趋势，交叉验证成立。

## 4. B4c 分泌因子-响应关联（LN229 vs U251）

### 4.1 C11 单核细胞响应（LN229 vs U251 共培养，17 显著基因）

其中 **5 个属 IFN ISG 且全部下调**：SIGLEC1(-1.10)、IFIT3(-1.02)、RSAD2(-1.06)、CXCL10(-1.18)、IFIT1(-1.11)。即：**LN229 共培养相对 U251 显著抑制单核细胞 ISG**，与论文"LN229 分泌 GM-CSF、U251 不分泌"的蛋白证据（ELISA, Fig 3c）一致。

### 4.2 细胞系分泌因子表达（log2 归一化）

| 因子 | LN229 | T98G | U87 | U251 | 备注 |
|---|---|---|---|---|---|
| CSF2 (GM-CSF) | 无信号 | 无信号 | 无信号 | 无信号 | 转录组不可测，蛋白证据见论文 |
| CSF1 (M-CSF) | 5.83 | 9.23 | 8.51 | 9.59 | U251 最高 |
| TGFB1 | 5.84 | 6.77 | 6.46 | 5.91 | 无显著差异 |
| TGFB2 | 7.55 | 10.41 | 6.78 | 9.72 | U251 更高 |
| MIF | 8.77 | 8.52 | 8.55 | 0.00 | U251 缺失 |
| CXCL10 | 0.79 | 4.62 | 0.86 | 0.00 | U251 缺失 |

**解读**：CSF2 转录本在全部数据集无信号（QuantSeq 3' 低深度 + 分泌型细胞因子转录本极低），无法从 mRNA 验证 GM-CSF 分泌；但论文 ELISA 明确 LN229/T98G/U87 分泌 GM-CSF、U251 不分泌。C11 中 LN229 共培养抑制 ISG 而 U251 不抑制，与 GM-CSF 分泌差异方向一致，支持"GM-CSF 分泌 → 单核细胞 ISG 抑制"归因。

## 5. 结论

1. **直接与间接共培养效应高度一致**（B4a）：ISG 抑制主要由可溶性因子介导。
2. **68 个共培养抑制型 ISG 核心集**（B4b）：IFN-β 诱导与共培养抑制的精确交集，构成机制核心证据。
3. **LN229 vs U251 分泌因子归因**（B4c）：LN229 共培养抑制 ISG、U251 不抑制，与 GM-CSF 蛋白分泌差异一致（论文 ELISA 证据），支持"胶质瘤细胞分泌 GM-CSF 抑制单核细胞 I 型 IFN 应答"。

## 6. 产出文件

- `B4a_CoCul_sig_in_C4.csv` / `B4a_CoCul_sig_in_C9.csv`：直接共培养显著基因在间接共培养中的效应
- `B4b_ISG_core_set.csv`：68 个共培养抑制型 ISG 核心集
- `B4c_C11_sig_with_ISG.csv`：C11 显著基因（含 ISG 标注）
"""

with open(os.path.join(OUT, "B4_report.md"), "w", encoding="utf-8") as f:
    f.write(report)
print("B4_report.md 已生成")
