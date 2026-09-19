# -*- coding: utf-8 -*-
"""B1 结果汇总：显著基因统计 + 关键 ISG 基因检查"""
import os, pandas as pd, numpy as np


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

OUT = os.path.join(PROJ_ROOT, "output", "B1")

contrasts = {
    "data1": ["C1_culture_main","C2_antibody_main","C3_interaction",
              "C4_co_vs_mono_IgG","C5_co_vs_mono_aGMCSF","C6_aGMCSF_vs_IgG_mono","C7_aGMCSF_vs_IgG_co"],
    "data2": ["C8_GMCSF_vs_ctrl","C9_LN229_vs_mono","C10_U251_vs_mono","C11_LN229_vs_U251","C12_pooled_co_vs_mono"],
}

summary_rows = []
for block, cids in contrasts.items():
    for cid in cids:
        f = os.path.join(OUT, f"{block}_{cid}.csv")
        df = pd.read_csv(f)
        # 检查列
        print(f"{block} {cid} cols: {list(df.columns)[:8]}")
        padj_col = [c for c in df.columns if "padj" in c.lower()][0]
        lfc_col = [c for c in df.columns if "log2FoldChange" in c][0]
        n_sig = int(((df[padj_col] < 0.05) & (df[lfc_col].abs() > 1)).sum())
        n_up = int(((df[padj_col] < 0.05) & (df[lfc_col] > 1)).sum())
        n_down = int(((df[padj_col] < 0.05) & (df[lfc_col] < -1)).sum())
        summary_rows.append({"block": block, "contrast": cid, "n_genes": len(df),
                             "sig_padj05_lfc1": n_sig, "up": n_up, "down": n_down})
        # 保存显著基因子集
        sig = df[(df[padj_col] < 0.05) & (df[lfc_col].abs() > 1)].sort_values(padj_col)
        sig.to_csv(os.path.join(OUT, f"{block}_{cid}_sig.csv"), index=False)

summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(OUT, "B1_summary.csv"), index=False)
print("\n===== B1 显著基因汇总 =====")
print(summary.to_string(index=False))

# 关键 ISG 基因检查（C4 co vs mono in IgG, C7 rescue, C8 GMCSF, C9/C10 co-culture）
ISG = ["ISG15","MX1","MX2","OAS1","OAS2","OAS3","OASL","IFIT1","IFIT2","IFIT3","IFITM1","IFITM3",
       "STAT1","STAT2","IRF7","IRF9","IFI6","IFI27","IFI44","IFI44L","RSAD2","BST2","GBP1","GBP2",
       "CXCL10","CXCL11","CCL2","CCL5","TNF","IL6","IL1B","CSF2","CSF2RA","CSF2RB","CD274","PDCD1LG2"]
gene_col = "symbol"
print("\n===== 关键基因在主要对比中的表现 =====")
for cid in ["C4_co_vs_mono_IgG","C7_aGMCSF_vs_IgG_co","C8_GMCSF_vs_ctrl","C9_LN229_vs_mono","C10_U251_vs_mono"]:
    block = "data1" if cid in contrasts["data1"] else "data2"
    f = os.path.join(OUT, f"{block}_{cid}.csv")
    df = pd.read_csv(f)
    padj_col = [c for c in df.columns if "padj" in c.lower()][0]
    lfc_col = [c for c in df.columns if "log2FoldChange" in c][0]
    sub = df[df[gene_col].isin(ISG)][[gene_col, lfc_col, padj_col]].copy()
    sub.columns = ["gene", "log2FC", "padj"]
    sub = sub.sort_values("log2FC", ascending=False)
    print(f"\n--- {cid} ({block}) ---")
    print(sub.to_string(index=False))
