# -*- coding: utf-8 -*-
"""B5 绘图: Figure 1-6 草图 + Supplementary Figure 草图
输出到 output/B5/figures/ 与 output/B5/suppl/
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
WS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(WS, "output", "B5")
FIG = os.path.join(OUT, "figures")
SUP = os.path.join(OUT, "suppl")
os.makedirs(FIG, exist_ok=True); os.makedirs(SUP, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "figure.dpi": 150})

core = pd.read_csv(os.path.join(OUT, "core68_evidence_matrix.csv"))
mech = pd.read_csv(os.path.join(OUT, "mechanism_chain_validation.csv"))
rescue = pd.read_csv(os.path.join(OUT, "rescue_analysis_full.csv"))
loo = pd.read_csv(os.path.join(OUT, "robustness_leave_one_donor_out.csv"))
perm = pd.read_csv(os.path.join(OUT, "permutation_distributions.csv"))
sens_thr = pd.read_csv(os.path.join(OUT, "robustness_threshold_sensitivity.csv"))
sens_rescue = pd.read_csv(os.path.join(OUT, "rescue_sensitivity_analysis.csv"))
sc = pd.read_csv(os.path.join(OUT, "external_scRNA_patient_level_validation.csv"))
ora = pd.read_csv(os.path.join(OUT, "enrichment_ORA_full.csv"))

# ---------- Figure 1: 68 基因核心集证据矩阵热图 ----------
lfc_cols = [c for c in core.columns if c.endswith("_log2FC")]
# 选关键 contrast
key = ["B2_IFN_vs_MoCul_log2FC","B2_CoCul_vs_MoCul_log2FC","C4_co_vs_mono_IgG_log2FC",
       "C7_aGMCSF_vs_IgG_co_log2FC","C8_GMCSF_vs_ctrl_log2FC","C9_LN229_vs_mono_log2FC",
       "C10_U251_vs_mono_log2FC","C11_LN229_vs_U251_log2FC"]
key = [c for c in key if c in core.columns]
sub = core.set_index("symbol")[key].copy()
sub.columns = [c.replace("_log2FC","") for c in sub.columns]
# 按 B2_CoCul 排序
sub = sub.sort_values("B2_CoCul_vs_MoCul")
fig, ax = plt.subplots(figsize=(7, 12))
im = ax.imshow(sub.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
ax.set_yticks(range(len(sub))); ax.set_yticklabels(sub.index, fontsize=6)
ax.set_xticks(range(len(sub.columns))); ax.set_xticklabels(sub.columns, rotation=45, ha="right", fontsize=7)
ax.set_title("Figure 1. Core68 ISG evidence matrix (log2FC)")
plt.colorbar(im, ax=ax, label="log2FC", shrink=0.4)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "Figure1_core68_heatmap.png")); plt.close()

# ---------- Figure 2: 机制链 A-F 验证 ----------
fig, ax = plt.subplots(figsize=(8, 4.5))
rows = mech[mech["ci_lo"].notna()]
labels = rows["label"].tolist()
y = np.arange(len(rows))
ax.errorbar(rows["mean_log2FC"], y, xerr=[rows["mean_log2FC"]-rows["ci_lo"], rows["ci_hi"]-rows["mean_log2FC"]],
            fmt="o", capsize=4)
ax.axvline(0, color="grey", ls="--", lw=0.8)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel("mean log2FC (95% CI)")
ax.set_title("Figure 2. Mechanism chain A-F: gene-set level effects")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "Figure2_mechanism_chain.png")); plt.close()

# ---------- Figure 3: 逆转分析 ----------
fig, ax = plt.subplots(figsize=(6, 5))
d = rescue.dropna(subset=["suppression_effect","rescue_effect"])
colors = {"fully_rescued":"#2ca02c","partially_rescued":"#1f77b4","not_rescued":"#d62728","discordant":"#7f7f7f","na":"#bbbbbb"}
for cls, grp in d.groupby("rescue_class"):
    ax.scatter(grp["suppression_effect"], grp["rescue_effect"], s=25, alpha=0.7,
               label=f"{cls} (n={len(grp)})", color=colors.get(cls,"#999"))
rho, p = stats.spearmanr(d["suppression_effect"], d["rescue_effect"])
ax.set_xlabel("suppression_effect (LN229 CoCul log2FC)")
ax.set_ylabel("rescue_effect (αGM-CSF log2FC)")
ax.axhline(0, color="grey", ls="--", lw=0.8); ax.axvline(0, color="grey", ls="--", lw=0.8)
ax.set_title(f"Figure 3A. Suppression vs rescue (Spearman rho={rho:.3f}, p={p:.2e})")
ax.legend(fontsize=7, loc="best")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "Figure3A_suppression_rescue_scatter.png")); plt.close()

# rescue heatmap: 68 基因 × 关键 contrast
key3 = ["B2_IFN_vs_MoCul_log2FC","B2_CoCul_vs_MoCul_log2FC","C4_co_vs_mono_IgG_log2FC",
        "C7_aGMCSF_vs_IgG_co_log2FC","C8_GMCSF_vs_ctrl_log2FC","C9_LN229_vs_mono_log2FC"]
key3 = [c for c in key3 if c in core.columns]
sub3 = core.set_index("symbol")[key3].copy()
sub3.columns = [c.replace("_log2FC","") for c in sub3.columns]
sub3 = sub3.sort_values("B2_CoCul_vs_MoCul")
fig, ax = plt.subplots(figsize=(6, 12))
im = ax.imshow(sub3.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
ax.set_yticks(range(len(sub3))); ax.set_yticklabels(sub3.index, fontsize=6)
ax.set_xticks(range(len(sub3.columns))); ax.set_xticklabels(sub3.columns, rotation=45, ha="right", fontsize=7)
ax.set_title("Figure 3B. Rescue heatmap (log2FC)")
plt.colorbar(im, ax=ax, label="log2FC", shrink=0.4)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "Figure3B_rescue_heatmap.png")); plt.close()

# ---------- Figure 4: 富集分析 ORA ----------
if "FDR" in ora.columns and "term" in ora.columns:
    top = ora[ora["FDR"]<0.05].sort_values("FDR").head(15)
    if len(top):
        fig, ax = plt.subplots(figsize=(7, 5))
        top = top.iloc[::-1]
        ax.barh(np.arange(len(top)), -np.log10(top["FDR"]), color="#1f77b4")
        ax.set_yticks(np.arange(len(top)))
        ax.set_yticklabels([f"{r['library']}: {r['term'][:40]}" for _, r in top.iterrows()], fontsize=6)
        ax.set_xlabel("-log10(FDR)")
        ax.set_title("Figure 4. Top enriched pathways (ORA, FDR<0.05)")
        plt.tight_layout(); plt.savefig(os.path.join(FIG, "Figure4_enrichment_ORA.png")); plt.close()

# ---------- Figure 5: 稳健性 ----------
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
# 5A: leave-one-donor-out
loo["contrast"] = loo["label"].str.extract(r"_(C\d+|B2)$")[0]
for ct, grp in loo.groupby("contrast"):
    axes[0].plot(grp["mean_log2FC"], "o-", label=ct, markersize=4)
axes[0].axhline(0, color="grey", ls="--", lw=0.8)
axes[0].set_ylabel("mean log2FC (leave-one-donor-out)")
axes[0].set_xlabel("left-out donor")
axes[0].set_title("Figure 5A. Leave-one-donor-out")
axes[0].legend(fontsize=7)
# 5B: permutation 分布
axes[1].hist(perm["perm_mean_log2FC"], bins=50, color="#999", alpha=0.7, label="random gene sets")
axes[1].axvline(-1.878, color="red", lw=2, label="observed (-1.878)")
axes[1].set_xlabel("mean log2FC (permuted)")
axes[1].set_title("Figure 5B. 10,000 permutations")
axes[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "Figure5_robustness.png")); plt.close()

# ---------- Figure 6: 外部 scRNA 患者层面 ----------
fig, ax = plt.subplots(figsize=(6, 4.5))
mo = sc[sc["lineage"]=="MO-TAM"].set_index("sample")["module_score"]
mg = sc[sc["lineage"]=="MG-TAM"].set_index("sample")["module_score"]
common = list(set(mo.index) & set(mg.index))
x = np.arange(len(common))
ax.plot(x, mo[common].values, "o-", color="#d62728", label="MO-TAM")
ax.plot(x, mg[common].values, "s-", color="#1f77b4", label="MG-TAM")
for i in x:
    ax.plot([i,i],[mo[common].values[i], mg[common].values[i]], color="grey", lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(common)
ax.set_ylabel("core68 module score (pseudobulk)")
ax.set_title("Figure 6. Patient-level MO-TAM vs MG-TAM (n=7)")
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "Figure6_scRNA_patient_level.png")); plt.close()

# ---------- Suppl: 阈值敏感性 ----------
fig, ax = plt.subplots(figsize=(6, 4))
piv = sens_thr.pivot_table(index="padj_cutoff", columns="log2FC_cutoff", values="overlap_n")
piv.plot(kind="bar", ax=ax)
ax.set_ylabel("overlap n (core set size)")
ax.set_title("Suppl. Threshold sensitivity of core68")
plt.tight_layout(); plt.savefig(os.path.join(SUP, "Suppl_threshold_sensitivity.png")); plt.close()

# ---------- Suppl: permutation 分布 (frac_negative + rho) ----------
fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
axes[0].hist(perm["perm_frac_negative"], bins=50, color="#999", alpha=0.7)
axes[0].axvline(1.0, color="red", lw=2, label="observed=1.0")
axes[0].set_xlabel("frac negative (permuted)"); axes[0].set_title("Suppl. Permutation: direction concordance")
axes[0].legend(fontsize=7)
axes[1].hist(perm["perm_rho"].dropna(), bins=50, color="#999", alpha=0.7)
axes[1].axvline(0.581, color="red", lw=2, label="observed=0.581")
axes[1].set_xlabel("rho GMCSF vs LN229 (permuted)"); axes[1].set_title("Suppl. Permutation: correlation")
axes[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(SUP, "Suppl_permutation_distributions.png")); plt.close()

# ---------- Suppl: rescue 敏感性 ----------
fig, ax = plt.subplots(figsize=(6, 4))
for ri in sens_rescue["ri_full"].unique():
    s = sens_rescue[sens_rescue["ri_full"]==ri].groupby("suppression_threshold")[["fully","partially","not_rescued"]].mean()
    ax.plot(s.index, s["fully"], "o-", label=f"fully (ri_full={ri})")
ax.set_xlabel("suppression threshold"); ax.set_ylabel("n fully rescued")
ax.set_title("Suppl. Rescue classification sensitivity")
ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(SUP, "Suppl_rescue_sensitivity.png")); plt.close()

print("Figures saved to", FIG)
print("Suppl saved to", SUP)
print(sorted(os.listdir(FIG)))
print(sorted(os.listdir(SUP)))
