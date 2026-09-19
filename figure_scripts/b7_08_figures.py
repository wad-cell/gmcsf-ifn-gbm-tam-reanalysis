# -*- coding: utf-8 -*-
"""B7 模块八：综合报告 + Figure 候选图 + Figure-claim-evidence 对应表"""
import os, pickle
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
WS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(WS, "output", "B7")
TMP = os.path.join(WS, "temp", "B7")

# ============ Figure S2: canonical15 vs context53 两队列森林图 ============
r163 = pd.read_csv(os.path.join(OUT, "canonical15_results.csv"))
r182 = pd.read_csv(os.path.join(OUT, "context53_results.csv"))
meta = pd.read_csv(os.path.join(OUT, "gmcsf_ifn_meta_analysis.csv"))

fig, ax = plt.subplots(figsize=(7, 5))
sets = ["canonical15", "context53", "Hallmark_IFNa", "Hallmark_IFNg"]
ypos = np.arange(len(sets)) * 2.2
for i, s in enumerate(sets):
    d163 = r163[r163["set"]==s].iloc[0]
    d182 = r182[r182["set"]==s].iloc[0]
    m = meta[meta["set"]==s].iloc[0]
    ax.errorbar(d163["mean_diff_MO_MG"], ypos[i]+0.35, fmt="o", color="tab:red", label="GSE163120" if i==0 else None)
    ax.errorbar(d182["mean_diff_MO_MG"], ypos[i]-0.35, fmt="s", color="tab:blue", label="GSE182109" if i==0 else None)
    ax.errorbar(m["pooled_diff"], ypos[i], fmt="D", color="black", label="meta pooled" if i==0 else None)
    ax.text(0.30, ypos[i], f"meta p={m['meta_p']:.3f}", fontsize=8, va="center")
ax.axvline(0, color="grey", ls="--", lw=0.8)
ax.set_yticks(ypos)
ax.set_yticklabels(sets)
ax.set_xlabel("MO-TAM minus MG-TAM (pseudobulk mean diff)")
ax.set_title("Classical IFN programs: MO vs MG-TAM (two cohorts)")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "Figure_S2_classical_IFN_forest.pdf"))
fig.savefig(os.path.join(OUT, "Figure_S2_classical_IFN_forest.png"), dpi=150)
plt.close(fig)

# ============ Figure S3: gene contribution 条形图 ============
gc = pd.read_csv(os.path.join(OUT, "gene_contribution_analysis.csv"))
top = gc.head(12)
fig, ax = plt.subplots(figsize=(8, 5))
colors = ["tab:red" if m=="context53" else "tab:blue" for m in top["module"]]
ax.barh(top["symbol"][::-1], top["diff_MO_minus_MG"][::-1], color=colors[::-1])
ax.axvline(0, color="grey", ls="--", lw=0.8)
ax.set_xlabel("MO_mean - MG_mean (log1p)")
ax.set_title("Top gene contributions to core68 MO-MG difference (GSE163120)")
ax.legend(["context53", "canonical15"], fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "Figure_S3_gene_contribution.pdf"))
fig.savefig(os.path.join(OUT, "Figure_S3_gene_contribution.png"), dpi=150)
plt.close(fig)

# ============ Figure S4: 模块四 GM-CSF-repressed vs mod15 散点 ============
scores = pd.read_csv(os.path.join(OUT, "B7_m4_cell_scores.csv"))
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, (x, y, t) in zip(axes, [
    ("GMCSF_rep_nonIFN_14_mean", "mod15_mean", "rep14 vs mod15 (mean)"),
    ("GMCSF_rep_nonISG_mean", "mod15_mean", "nonISG vs mod15 (mean)")]):
    ax.scatter(scores[x], scores[y], s=1, alpha=0.1, color="tab:red")
    ax.set_xlabel(x); ax.set_ylabel(y); ax.set_title(t)
fig.suptitle("MO-TAM internal: GM-CSF-repressed vs classical IFN (GSE163120)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "Figure_S4_gmcsf_ifn_corr.pdf"))
fig.savefig(os.path.join(OUT, "Figure_S4_gmcsf_ifn_corr.png"), dpi=150)
plt.close(fig)
print("已生成 Figure_S2/S3/S4 (pdf+png)")

# ============ Figure-claim-evidence 对应表 ============
fce = pd.DataFrame([
    ["Fig 1", "GM-CSF 共培养抑制单核细胞 ISG 程序（core68）", "Level 1 体外因果", "B1 C8/C9；B5 机制链 A-F", "主文"],
    ["Fig 2", "αGM-CSF 中和逆转 ISG 抑制（rescue）", "Level 1 体外因果", "B1 C7；B5 rescue 分析", "主文"],
    ["Fig 3", "GM-CSF 抑制为可溶性因子介导（直接 vs 间接一致）", "Level 1 体外因果", "B5 机制链 F rho=0.32", "主文/补充"],
    ["Fig 4", "68 基因 signature 富集于 IFN 通路", "Level 1 体外因果", "B5 ORA/GSEA", "主文"],
    ["Fig 5", "MO-TAM 内部 GM-CSF-repressed 与经典 IFN 正相关", "Level 2 患者状态关联", "B7 模块四 partial rho=0.585/0.236", "主文（关联性表述）"],
    ["Fig 6", "canonical15 在 MO vs MG-TAM 无差异", "Level 3 不成立", "B7 模块三 meta p=0.614", "补充（阴性）"],
    ["Fig 7", "context53 在 MO-TAM 偏高（NUPR1/NAMPT/ATF3 驱动）", "Level 2 患者状态关联", "B7 模块二/三 meta p=0.040", "补充"],
    ["Fig 8", "rescue 模块患者层面无差异", "Level 3 不成立", "B7 模块六两队列均 ns", "补充（阴性）"],
    ["Fig S1", "打分方法敏感性（mean/ucell/aucell）", "稳健性", "B7 模块二 signature_method_sensitivity", "补充"],
    ["Fig S2", "经典 IFN 程序两队列森林图", "稳健性", "B7 模块三", "补充"],
    ["Fig S3", "基因贡献（NUPR1/NAMPT/ATF3 驱动 MO 高）", "稳健性", "B7 模块二 gene_contribution", "补充"],
    ["Fig S4", "GM-CSF-repressed vs mod15 散点", "稳健性", "B7 模块四", "补充"],
], columns=["Figure", "Claim", "证据等级", "证据来源", "位置"])
fce.to_csv(os.path.join(OUT, "B7_figure_claim_evidence.csv"), index=False)
print("已保存 B7_figure_claim_evidence.csv")
print(fce.to_string(index=False))
