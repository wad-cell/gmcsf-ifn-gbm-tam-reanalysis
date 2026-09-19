# -*- coding: utf-8 -*-
"""B4a 聚焦版：显著基因子集效应量一致性 + B4b ISG 核心集 + B4c 分泌因子-响应关联"""
import pandas as pd, numpy as np, os
from scipy import stats

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

OUT = os.path.join(PROJ_ROOT, "output", "B4")
os.makedirs(OUT, exist_ok=True)
B1 = os.path.join(PROJ_ROOT, "output", "B1")
B2 = os.path.join(PROJ_ROOT, "output", "B2")
B3 = os.path.join(PROJ_ROOT, "output", "B3")

def load(p):
    return pd.read_csv(p)

direct = load(f"{B2}\\GSE309037_B2_CoCul_vs_MoCul.csv")
ind_igG = load(f"{B1}\\data1_C4_co_vs_mono_IgG.csv")
ind_LN = load(f"{B1}\\data2_C9_LN229_vs_mono.csv")
ifn = load(f"{B2}\\GSE309037_B2_IFN_vs_MoCul.csv")
c11 = load(f"{B1}\\data2_C11_LN229_vs_U251.csv")

# ---- B4a：显著基因子集效应量一致性 ----
# 以 B2 CoCul 显著基因为锚，看其在 C4/C9 中的效应
cocul_sig = direct[(direct["padj"]<0.05)&(direct["log2FoldChange"].abs()>1)].dropna(subset=["symbol"])
print(f"B2 CoCul 显著基因: {len(cocul_sig)}")

def anchor_effect(anchor, target, anchor_name, target_name):
    """anchor 显著基因在 target 中的 log2FC 分布"""
    t = target.set_index("gene")[["log2FoldChange","padj"]]
    t.columns = [f"{target_name}_log2FC", f"{target_name}_padj"]
    m = anchor[["gene","symbol","log2FoldChange","padj"]].merge(t, left_on="gene", right_index=True, how="left")
    m.columns = ["gene","symbol",f"{anchor_name}_log2FC",f"{anchor_name}_padj",f"{target_name}_log2FC",f"{target_name}_padj"]
    return m

m_c4 = anchor_effect(cocul_sig, ind_igG, "CoCul_direct", "C4_indIgG")
m_c9 = anchor_effect(cocul_sig, ind_LN, "CoCul_direct", "C9_indLN")

# 方向一致性：CoCul 下调基因在 C4/C9 中是否也下调
for m, lab in [(m_c4,"C4 间接(IgG)"), (m_c9,"C9 间接(LN229)")]:
    down_anchor = m[m["CoCul_direct_log2FC"]<0]
    up_anchor = m[m["CoCul_direct_log2FC"]>0]
    tcol = [c for c in m.columns if c.endswith("_log2FC") and c!="CoCul_direct_log2FC"][0]
    d = down_anchor[tcol].dropna()
    u = up_anchor[tcol].dropna()
    print(f"\n[{lab}] CoCul 下调基因({len(down_anchor)}) 在 {lab} 中: 中位 log2FC={d.median():.2f}, 同向下调比例={(d<0).mean()*100:.0f}%")
    print(f"[{lab}] CoCul 上调基因({len(up_anchor)}) 在 {lab} 中: 中位 log2FC={u.median():.2f}, 同向上调比例={(u>0).mean()*100:.0f}%")
    # 显著子集相关
    sig_both = m.dropna(subset=[tcol])
    sig_both = sig_both[np.isfinite(sig_both[tcol])]
    if len(sig_both)>10:
        r_s, p_s = stats.spearmanr(sig_both["CoCul_direct_log2FC"], sig_both[tcol])
        print(f"[{lab}] 显著基因子集 Spearman rho={r_s:.3f} (p={p_s:.2e}, n={len(sig_both)})")

m_c4.to_csv(f"{OUT}\\B4a_CoCul_sig_in_C4.csv", index=False)
m_c9.to_csv(f"{OUT}\\B4a_CoCul_sig_in_C9.csv", index=False)

# ---- B4b：ISG 核心集 ----
ifn_up = ifn[(ifn["padj"]<0.05)&(ifn["log2FoldChange"]>1)].dropna(subset=["symbol"])
ifn_set = set(ifn_up["symbol"])
cocul_down = direct[(direct["padj"]<0.05)&(direct["log2FoldChange"]<-1)].dropna(subset=["symbol"])
cocul_down_set = set(cocul_down["symbol"])
core = ifn_set & cocul_down_set
print(f"\nIFN ISG(760) ∩ CoCul 下调({len(cocul_down_set)}) = {len(core)}")

core_df = ifn_up[ifn_up["symbol"].isin(core)][["gene","symbol","log2FoldChange","padj"]].copy()
core_df.columns = ["gene","symbol","IFN_log2FC","IFN_padj"]
cd = cocul_down.set_index("symbol")[["log2FoldChange","padj"]]
cd.columns = ["CoCul_log2FC","CoCul_padj"]
core_df = core_df.merge(cd, left_on="symbol", right_index=True, how="left")
# 附加 C4/C9 效应（用 gene 列映射，避免 symbol 重复）
for nm, tgt in [("C4", ind_igG), ("C9", ind_LN)]:
    t = tgt.drop_duplicates("gene").set_index("gene")["log2FoldChange"]
    core_df[f"{nm}_log2FC"] = core_df["gene"].map(t)
core_df = core_df.sort_values("IFN_log2FC", ascending=False)
core_df.to_csv(f"{OUT}\\B4b_ISG_core_set.csv", index=False)
print(f"B4b ISG 核心集已保存: {len(core_df)} 基因")

# ---- B4c：分泌因子-响应关联 ----
# C11 = LN229 vs U251 共培养单核细胞响应（17 显著）
c11_sig = c11[(c11["padj"]<0.05)&(c11["log2FoldChange"].abs()>1)].dropna(subset=["symbol"])
print(f"\nC11 LN229 vs U251 共培养显著: {len(c11_sig)}")
# 其中 ISG（IFN 金标准）下调 → 支持 LN229 分泌因子抑制 ISG
c11_isg = c11_sig[c11_sig["symbol"].isin(ifn_set)]
print(f"C11 显著基因中属 IFN ISG: {len(c11_isg)}")
print(c11_isg[["symbol","log2FoldChange","padj"]].to_string(index=False))

# 细胞系分泌因子表达（LN229 vs U251）
cell = load(f"{B3}\\GSE309038_key_genes_cell_mean_log2.csv").set_index("symbol")
secr = ["CSF2","CSF1","TGFB1","TGFB2","TGFB3","IL6","IL1B","CCL2","MIF","CXCL10"]
print("\n细胞系分泌因子 log2 表达 (LN229 vs U251):")
for g in secr:
    if g in cell.index:
        r = cell.loc[g]
        print(f"  {g}: LN229={r['LN229']:.2f}, T98G={r['T98G']:.2f}, U87={r['U87']:.2f}, U251={r['U251']:.2f}")
    else:
        print(f"  {g}: 无表达信号")

# 保存 C11 显著基因表（含 ISG 标注）
c11_sig["is_IFN_ISG"] = c11_sig["symbol"].isin(ifn_set)
c11_sig.to_csv(f"{OUT}\\B4c_C11_sig_with_ISG.csv", index=False)
print(f"\nB4c C11 显著基因表已保存")
