# -*- coding: utf-8 -*-
"""B2 报告 + 跨数据集 ISG 重叠分析（GSE309037 vs GSE309039 data2 C8）"""
import pandas as pd, os

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

OUT_B2 = os.path.join(PROJ_ROOT, "output", "B2")
OUT_B1 = os.path.join(PROJ_ROOT, "output", "B1")

def load(block, cid, outdir):
    return pd.read_csv(os.path.join(outdir, f"{block}_{cid}.csv"))

# 1. 显著基因汇总
print("===== B2 显著基因汇总 (padj<0.05, |log2FC|>1) =====")
for cid in ["B2_IFN_vs_MoCul", "B2_CoCul_vs_MoCul", "B2_CoCul_vs_IFN"]:
    df = load("GSE309037", cid, OUT_B2)
    n_sig = int(((df["padj"]<0.05)&(df["log2FoldChange"].abs()>1)).sum())
    n_up = int(((df["padj"]<0.05)&(df["log2FoldChange"]>1)).sum())
    n_down = int(((df["padj"]<0.05)&(df["log2FoldChange"]<-1)).sum())
    print(f"{cid}: 显著={n_sig}, 上调={n_up}, 下调={n_down}")

# 2. IFN-β 金标准 ISG 签名（显著上调）
ifn = load("GSE309037", "B2_IFN_vs_MoCul", OUT_B2)
ifn_sig_up = ifn[(ifn["padj"]<0.05)&(ifn["log2FoldChange"]>1)].sort_values("log2FoldChange", ascending=False)
ifn_sig_up.to_csv(os.path.join(OUT_B2, "IFN_gold_standard_ISG_up.csv"), index=False)
print(f"\nIFN-β 金标准 ISG 上调基因数: {len(ifn_sig_up)}")
print("Top 30:")
print(ifn_sig_up[["symbol","log2FoldChange","padj"]].head(30).to_string(index=False))

# 3. 直接共培养 ISG 表现
cocul = load("GSE309037", "B2_CoCul_vs_MoCul", OUT_B2)
cocul_sig = cocul[(cocul["padj"]<0.05)&(cocul["log2FoldChange"].abs()>1)].sort_values("padj")
cocul_sig.to_csv(os.path.join(OUT_B2, "CoCul_vs_MoCul_sig.csv"), index=False)
print(f"\n直接共培养显著基因数: {len(cocul_sig)}")

# 4. 跨数据集重叠：IFN 金标准 ISG vs GM-CSF 抑制 ISG (C8)
c8 = load("data2", "C8_GMCSF_vs_ctrl", OUT_B1)
c8_sig_down = c8[(c8["padj"]<0.05)&(c8["log2FoldChange"]<-1)]
print(f"\nC8 (GM-CSF) 显著下调基因数: {len(c8_sig_down)}")

ifn_set = set(ifn_sig_up["symbol"].dropna())
c8_set = set(c8_sig_down["symbol"].dropna())
overlap = ifn_set & c8_set
print(f"IFN 金标准 ISG 上调 ∩ GM-CSF 下调: {len(overlap)}")
print("重叠基因:", sorted(overlap))

# 5. 关键 ISG 在 B2 三个对比中的表现
ISG = ["ISG15","MX1","MX2","OAS1","OAS2","OAS3","OASL","IFIT1","IFIT2","IFIT3","IFITM1","IFITM3",
       "STAT1","STAT2","IRF7","IRF9","IFI6","IFI27","IFI44","IFI44L","RSAD2","BST2","GBP1","GBP2",
       "CXCL10","CXCL11","CCL2","CCL5","TNF","IL6","IL1B","CSF2","CSF2RA","CSF2RB","CD274","PDCD1LG2"]
print("\n===== 关键基因在 B2 三对比中的表现 =====")
for cid in ["B2_IFN_vs_MoCul", "B2_CoCul_vs_MoCul", "B2_CoCul_vs_IFN"]:
    df = load("GSE309037", cid, OUT_B2)
    sub = df[df["symbol"].isin(ISG)][["symbol","log2FoldChange","padj"]].copy()
    sub = sub.sort_values("log2FoldChange", ascending=False)
    print(f"\n--- {cid} ---")
    print(sub.to_string(index=False))

# 6. 保存重叠分析结果
overlap_df = ifn_sig_up[ifn_sig_up["symbol"].isin(c8_set)][["symbol","log2FoldChange","padj"]].copy()
overlap_df.columns = ["symbol","IFN_log2FC","IFN_padj"]
c8_map = c8.set_index("symbol")[["log2FoldChange","padj"]]
c8_map.columns = ["GMCSF_log2FC","GMCSF_padj"]
overlap_df = overlap_df.merge(c8_map, left_on="symbol", right_index=True, how="left")
overlap_df.to_csv(os.path.join(OUT_B2, "IFN_ISG_vs_GMCSF_overlap.csv"), index=False)
print(f"\n重叠基因表已保存: {len(overlap_df)} 行")
