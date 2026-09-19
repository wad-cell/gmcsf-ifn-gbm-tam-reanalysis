# -*- coding: utf-8 -*-
"""C2 正式制图：Main Figures 1-5
TIFF 600dpi（投稿）+ PNG 预览（150dpi）+ figure_source_data.csv
所有数值来自 final 结果表。
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.stats import t, spearmanr
WS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(WS, "output", "B7")
B1 = os.path.join(WS, "output", "B1")
B2 = os.path.join(WS, "output", "B2")
B5 = os.path.join(WS, "output", "B5")
B6 = os.path.join(WS, "output", "B6")
C2T = os.path.join(WS, "temp", "C2")
FIG = os.path.join(WS, "output", "C2", "figures")
SRC = os.path.join(WS, "output", "C2", "figure_source_data")
os.makedirs(FIG, exist_ok=True); os.makedirs(SRC, exist_ok=True)

def rd(p, base=None):
    b = base if base else OUT
    return pd.read_csv(os.path.join(b, p))

def save(fig, name, dpi=300):
    fig.savefig(os.path.join(FIG, f"{name}.tif"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, f"{name}.png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, f"{name}.pdf"), bbox_inches="tight")
    plt.close(fig)

def mean_ci(x):
    x = np.asarray(x, dtype=float); x = x[~np.isnan(x)]
    n = len(x); m = x.mean(); s = x.std(ddof=1)
    h = t.ppf(0.975, n-1) * s / np.sqrt(n) if n > 1 else np.nan
    return m, m-h, m+h, n

# =================================================================
# Figure 1. Study design schematic
# =================================================================
fig, ax = plt.subplots(figsize=(11, 6.2))
ax.set_xlim(0, 11); ax.set_ylim(0, 6.2); ax.axis("off")
def box(x, y, w, h, text, fc="#eef4fb", ec="#2b6cb0", fs=8.2, lw=1.2):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec=ec, lw=lw)
    ax.add_patch(p); ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs, wrap=True)

box(0.1, 5.2, 3.0, 0.8, "Perturbation datasets\nGSE309037 (direct co-culture + IFN-β gold standard, 3 donors, paired)\nGSE309039 data1 (indirect + α-GM-CSF, 6 donors, paired)\nGSE309039 data2 (indirect LN229/U251 + exogenous GM-CSF, 6 donors, paired)\nGSE309038 (GBM cell-line baselines)", fs=7.0)
box(3.6, 5.2, 3.0, 0.8, "Contrast definition\nC1–C12 (C13 invalid, excluded)\nDonor-paired designs\nDonor as blocking factor", fs=7.0)
box(7.1, 5.2, 3.4, 0.8, "68-gene core program (core68)\nIFN-β-induced 760 genes ∩\nco-culture-repressed 338 genes\n(hypergeometric P=3.2e-33, 5.9-fold)", fs=7.0)
box(0.1, 3.9, 3.0, 0.8, "Contrast-level evidence matrix\nIFN-β / direct / indirect /\nGM-CSF / neutralization /\nU251 control", fs=7.0)
box(3.6, 3.9, 3.0, 0.8, "Suppression–rescue index\n31 fully / 12 partial /\n8 not / 17 not assessable\n(16 genes padj<0.05)", fs=7.0)
box(7.1, 3.9, 3.4, 0.8, "Decomposition\ncanonical15 (∩ Hallmark IFN-α)\ncontext53 (non-overlap)", fs=7.0)
box(0.1, 2.6, 3.0, 0.8, "Patient-level extrapolation\nGSE163120 (n=7 patients)\nGSE182109 (n=11 patients)\nPseudobulk paired + random-effects meta", fs=7.0)
box(3.6, 2.6, 3.0, 0.8, "Within-MO-TAM exploration\nGM-CSF-repressed signature vs\ncanonical IFN module\n(partial rho=0.585, exploratory)", fs=7.0)
box(7.1, 2.6, 3.4, 0.8, "Evidence levels\nLevel 1: in vitro causal\nLevel 2: patient association\nLevel 3: not supported", fs=7.0)
# arrows
for (x1,y1,x2,y2) in [(1.6,5.2,1.6,4.7),(5.1,5.2,5.1,4.7),(8.8,5.2,8.8,4.7),
                      (1.6,3.9,1.6,3.4),(5.1,3.9,5.1,3.4),(8.8,3.9,8.8,3.4)]:
    ax.annotate("", xy=(x2,y2), xytext=(x1,y1), arrowprops=dict(arrowstyle="->", color="#555", lw=1.4))
# pathway labels
ax.text(0.1, 5.0, "A", fontsize=11, fontweight="bold")
ax.text(0.1, 3.7, "B", fontsize=11, fontweight="bold")
ax.text(0.1, 2.4, "C", fontsize=11, fontweight="bold")
ax.set_title("Figure 1. Study design, datasets, and analysis workflow", fontsize=11)
save(fig, "Figure1_design")

# =================================================================
# Figure 2. IFN-β-defined reference program and co-culture suppression
# =================================================================
ifn = pd.read_csv(os.path.join(B2, "GSE309037_B2_IFN_vs_MoCul.csv"))
coc = pd.read_csv(os.path.join(B2, "GSE309037_B2_CoCul_vs_MoCul.csv"))
core = rd("core68_reannotation.csv")
core68 = core["symbol"].tolist()
isg_up = set(pd.read_csv(os.path.join(B2, "IFN_gold_standard_ISG_up.csv"))["symbol"])
core68set = set(core68)

fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5))
# A: IFN volcano
for df, ax, lab, sigset, col in [
    (ifn, axes[0,0], "IFN-β vs control", isg_up, "#c44e52"),
    (coc, axes[0,1], "Direct co-culture vs control", core68set, "#4c72b0")]:
    d = df.dropna(subset=["log2FoldChange", "padj"]).copy()
    d["sig"] = d["symbol"].isin(sigset)
    d = d[np.isfinite(d["log2FoldChange"])]
    ax.scatter(d.loc[~d["sig"],"log2FoldChange"], -np.log10(d.loc[~d["sig"],"padj"].clip(lower=1e-300)),
               s=1, color="#bbbbbb", alpha=0.4, linewidths=0)
    ax.scatter(d.loc[d["sig"],"log2FoldChange"], -np.log10(d.loc[d["sig"],"padj"].clip(lower=1e-300)),
               s=2, color=col, alpha=0.8, linewidths=0)
    ax.axvline(0, color="grey", lw=0.7, ls="--")
    ax.set_xlabel("log2FC"); ax.set_ylabel("-log10 padj"); ax.set_title(lab, fontsize=9)
# C: heatmap core68 across contrasts
em = rd("core68_evidence_matrix.csv", B5).set_index("symbol")
cols = [("B2_IFN_vs_MoCul_log2FC","IFN-β"),
        ("B2_CoCul_vs_MoCul_log2FC","Direct CoCul"),
        ("C9_LN229_vs_mono_log2FC","Indirect LN229"),
        ("C8_GMCSF_vs_ctrl_log2FC","GM-CSF"),
        ("C7_aGMCSF_vs_IgG_co_log2FC","αGM-CSF (rescue)"),
        ("C10_U251_vs_mono_log2FC","U251 (neg)")]
mat = em[[c for c,_ in cols]].astype(float)
mat = mat.reindex(core68).dropna(how="all")
mat = mat.fillna(0.0)
order = mat["B2_IFN_vs_MoCul_log2FC"].sort_values().index
mat = mat.loc[order]
im = axes[0,2].imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-3.5, vmax=3.5)
axes[0,2].set_xticks(range(len(cols))); axes[0,2].set_xticklabels([l for _,l in cols], rotation=45, ha="right", fontsize=7)
axes[0,2].set_yticks([])
axes[0,2].set_title("core68 (n=68 genes) log2FC across contrasts", fontsize=9)
axes[0,2].text(-2.5, 0, "sorted by IFN-β log2FC", fontsize=7, rotation=90, va="center")
cbar = fig.colorbar(im, ax=axes[0,2], fraction=0.046, pad=0.04); cbar.set_label("log2FC", fontsize=8)
# D: direct vs indirect scatter
d_eff = em["B2_CoCul_vs_MoCul_log2FC"]
i_eff = em["C9_LN229_vs_mono_log2FC"]
ddf = pd.DataFrame({"direct": d_eff, "indirect": i_eff}).dropna()
rho, pv = spearmanr(ddf["direct"], ddf["indirect"])
axes[1,0].scatter(ddf["direct"], ddf["indirect"], s=8, color="#4c72b0", alpha=0.7, edgecolors="white", linewidths=0.3)
axes[1,0].axvline(0, color="grey", lw=0.6, ls="--"); axes[1,0].axhline(0, color="grey", lw=0.6, ls="--")
axes[1,0].set_xlabel("Direct co-culture log2FC"); axes[1,0].set_ylabel("Indirect (LN229) log2FC")
axes[1,0].set_title(f"Direct vs indirect (rho={rho:.3f}, P={pv:.1e}, n={len(ddf)})", fontsize=9)
# E: leave-one-donor
loo = rd("robustness_leave_one_donor_out.csv", B5)
loo = loo[loo["label"].str.contains("C4")].sort_values("label")
axes[1,1].bar(range(len(loo)), loo["mean_log2FC"], color="#4c72b0")
axes[1,1].axhline(-1.878, color="black", ls="--", lw=1)
axes[1,1].text(len(loo)-0.5, -1.88, "full data", fontsize=7, ha="right", va="bottom")
axes[1,1].set_xticks(range(len(loo)))
axes[1,1].set_xticklabels([l.replace("data1_loo_donor_","D").replace("_C4","") for l in loo["label"]], fontsize=7)
axes[1,1].set_ylabel("core68 mean log2FC"); axes[1,1].set_title("Leave-one-donor-out (direct co-culture)", fontsize=9)
axes[1,1].set_ylim(-0.45, -0.2)
# B panel placeholder note (volcano B merged into A panel 2)
axes[1,2].axis("off")
axes[1,2].text(0.02, 0.5, "Panel B: 338 co-culture-repressed genes;\nintersection with IFN-β set = core68 (n=68),\n5.9-fold enriched (hypergeometric P=3.2e-33;\npermutation-confirmed).\n\nObserved core68 suppression: mean log2FC=-1.878,\n100% directionally concordant, P_emp=0.0001\n(three matched-permutation schemes; Fig. S1).",
               fontsize=8.5, va="center", transform=axes[1,2].transAxes)
fig.suptitle("Figure 2. IFN-β-defined reference program and co-culture suppression", fontsize=12)
save(fig, "Figure2_IFN_program")
fig, axes = plt.subplots(); axes.remove(); plt.close(fig)
pd.DataFrame({"symbol": mat.index, **{lbl: mat[c].values for c, lbl in cols}}).to_csv(
    os.path.join(SRC, "fig2_heatmap_source.csv"), index=False)
ddf.assign(direct=ddf["direct"], indirect=ddf["indirect"]).to_csv(
    os.path.join(SRC, "fig2_direct_vs_indirect.csv"), index=False)
loo.to_csv(os.path.join(SRC, "fig2_leave_one_donor.csv"), index=False)
print("Fig2 done")

# =================================================================
# Figure 3. Exogenous GM-CSF partially recapitulates
# =================================================================
fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.5))
# A: C8 vs C9 scatter
c8 = em["C8_GMCSF_vs_ctrl_log2FC"]; c9 = em["C9_LN229_vs_mono_log2FC"]
gdf = pd.DataFrame({"GMCSF": c8, "LN229": c9}).dropna()
rho, pv = spearmanr(gdf["GMCSF"], gdf["LN229"])
conc = ((np.sign(gdf["GMCSF"]) == np.sign(gdf["LN229"])).mean()*100)
axes[0,0].scatter(gdf["GMCSF"], gdf["LN229"], s=9, color="#c44e52", alpha=0.75, edgecolors="white", linewidths=0.3)
axes[0,0].axvline(0, color="grey", lw=0.6, ls="--"); axes[0,0].axhline(0, color="grey", lw=0.6, ls="--")
axes[0,0].set_xlabel("Exogenous GM-CSF (C8) log2FC"); axes[0,0].set_ylabel("LN229 indirect (C9) log2FC")
axes[0,0].set_title(f"GM-CSF vs LN229 co-culture (rho={rho:.3f}, {conc:.0f}% concordant)", fontsize=9)
# B: permutation null
perm = pd.read_csv(os.path.join(B5, "permutation_distributions.csv"))
axes[0,1].hist(perm["perm_rho"], bins=60, color="#8c8c8c", alpha=0.8, edgecolor="white")
axes[0,1].axvline(0.581, color="#c44e52", lw=1.8)
axes[0,1].set_xlabel("Permutation rho (N=10,000)"); axes[0,1].set_ylabel("Count")
axes[0,1].set_title("Permutation null (observed rho=0.581, P_emp=0.0173)", fontsize=9)
# C: GM-CSF effect on core68
m, lo, hi, n = mean_ci(em["C8_GMCSF_vs_ctrl_log2FC"].values)
axes[1,0].bar(["GM-CSF on core68"], [m], yerr=[[m-lo],[hi-m]], color="#c44e52", capsize=6)
axes[1,0].axhline(0, color="grey", lw=0.7)
axes[1,0].set_ylabel("mean log2FC (95% CI)")
axes[1,0].set_title(f"mean={m:.2f} (95% CI {lo:.2f}–{hi:.2f}); {n} genes", fontsize=9)
# D: U251 negative control
m4 = pd.read_csv(os.path.join(B5, "mechanism_chain_validation.csv")).set_index("label")
rows = [("C_U251_indirect_CoCul","U251 indirect"),
        ("B_LN229_indirect_CoCul","LN229 indirect")]
vals = [(m4.loc[l,"mean_log2FC"], m4.loc[l,"ci_lo"], m4.loc[l,"ci_hi"], lab) for l, lab in rows]
for i,(mu,lo_,hi_,lab) in enumerate(vals):
    axes[1,1].errorbar(mu, i, xerr=[[mu-lo_],[hi_-mu]], fmt="o", color="#4c72b0", capsize=5, ms=7)
axes[1,1].axvline(0, color="grey", lw=0.7)
axes[1,1].set_yticks(range(len(vals))); axes[1,1].set_yticklabels([v[3] for v in vals], fontsize=8)
axes[1,1].set_xlabel("core68 mean log2FC (95% CI)")
axes[1,1].set_title("Indirect co-culture: LN229 vs U251", fontsize=9)
axes[1,1].invert_yaxis()
fig.suptitle("Figure 3. Exogenous GM-CSF partially recapitulates the co-culture effect", fontsize=12)
save(fig, "Figure3_GMCSF_recapitulation")
gdf.to_csv(os.path.join(SRC, "fig3_gmcsf_vs_lN229.csv"), index=False)
pd.DataFrame({"perm_rho": perm["perm_rho"]}).to_csv(os.path.join(SRC, "fig3_permutation_null.csv"), index=False)
print("Fig3 done")

# =================================================================
# Figure 4. GM-CSF neutralization partially restores
# =================================================================
res = rd("rescue_analysis_full.csv", B5)
cls_counts = res["rescue_class"].value_counts()
order = ["fully_rescued","partially_rescued","not_rescued","na"]
counts = [int((res["rescue_class"]==c).sum()) for c in order]
n_sig = int(res["rescue_significant"].sum())
# module-level rescue mean (from C7)
r_eff = em["C7_aGMCSF_vs_IgG_co_log2FC"].values
m, lo, hi, n = mean_ci(r_eff)
fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2))
axes[0,0].bar(["αGM-CSF vs IgG\n(within co-culture)"], [m], yerr=[[m-lo],[hi-m]], color="#55a868", capsize=6)
axes[0,0].axhline(0, color="grey", lw=0.7)
axes[0,0].set_ylabel("core68 mean log2FC (95% CI)")
axes[0,0].set_title(f"Module-level restoration: {m:.2f} (95% CI {lo:.2f}–{hi:.2f}); {n} genes", fontsize=9)
# B: classification counts
bcols = ["#55a868","#4c72b0","#c44e52","#cccccc"]
axes[0,1].bar(range(4), counts, color=bcols)
axes[0,1].set_xticks(range(4))
axes[0,1].set_xticklabels([f"{c}\n(fully)" if c=="fully_rescued" else f"{c}\n(partial)" if c=="partially_rescued" else f"{c}\n(not)" if c=="not_rescued" else f"{c}\n(NA)" for c in order], fontsize=7)
for i, v in enumerate(counts):
    axes[0,1].text(i, v+0.3, str(v), ha="center", fontsize=9)
axes[0,1].set_ylabel("Genes (of 68)"); axes[0,1].set_title("Suppression–rescue classification", fontsize=9)
axes[0,1].set_ylim(0, 40)
# C: rescue_index by class (strip)
valid = res.dropna(subset=["rescue_index"])
cmap = {"fully_rescued":"#55a868","partially_rescued":"#4c72b0","not_rescued":"#c44e52","na":"#cccccc"}
for i, c in enumerate(order):
    sub = valid[valid["rescue_class"]==c]
    if len(sub)==0: continue
    xs = np.random.normal(i, 0.08, len(sub))
    axes[1,0].scatter(xs, sub["rescue_index"], s=12, color=cmap[c], alpha=0.8)
axes[1,0].axhline(0.8, color="grey", lw=0.7, ls="--"); axes[1,0].axhline(0.3, color="grey", lw=0.7, ls="--")
axes[1,0].set_xticks(range(4)); axes[1,0].set_xticklabels(["fully","partial","not","NA"], fontsize=8)
axes[1,0].set_ylabel("Rescue index"); axes[1,0].set_title("Rescue index distribution by class", fontsize=9)
# D: sensitivity
sens = rd("rescue_sensitivity_analysis.csv", B5)
sens_g = sens[sens["ri_full"]==0.8]
axes[1,1].plot(sens_g["suppression_threshold"], sens_g["fully"], "o-", color="#55a868", label="fully rescued")
axes[1,1].plot(sens_g["suppression_threshold"], sens_g["partially"], "s-", color="#4c72b0", label="partially")
axes[1,1].plot(sens_g["suppression_threshold"], sens_g["not_rescued"], "^-", color="#c44e52", label="not rescued")
axes[1,1].set_xlabel("Suppression threshold (log2FC)"); axes[1,1].set_ylabel("Genes")
axes[1,1].set_title("Sensitivity (rescue_index thresholds: full≥0.8, partial≥0.3)", fontsize=9)
axes[1,1].legend(fontsize=7)
fig.suptitle("Figure 4. GM-CSF neutralization partially restores the suppressed program", fontsize=12)
save(fig, "Figure4_neutralization")
res[["symbol","suppression_effect","rescue_effect","rescue_index","rescue_class","rescue_significant"]].to_csv(
    os.path.join(SRC, "fig4_rescue_genes.csv"), index=False)
sens.to_csv(os.path.join(SRC, "fig4_sensitivity.csv"), index=False)
print("Fig4 done")

# =================================================================
# Figure 5. Patient-level extrapolation boundaries (balanced)
# =================================================================
pp163 = pd.read_csv(os.path.join(C2T, "patient_paired_GSE163120.csv"))
pp182 = pd.read_csv(os.path.join(C2T, "patient_paired_GSE182109.csv"))
meta = rd("gmcsf_ifn_meta_analysis.csv")
gc = rd("gene_contribution_analysis.csv")
logo = rd("leave_one_gene_out_results.csv")
cs = pd.read_csv(os.path.join(OUT, "B7_m4_cell_scores.csv"))

fig = plt.figure(figsize=(14.5, 9.5))
# A: GSE163120 canonical15 paired
ax = fig.add_subplot(3, 3, 1)
d = pp163[pp163["module"]=="canonical15"].pivot(index="sample", columns="lineage", values="score")
for _, r in d.iterrows():
    ax.plot(["MG-TAM","MO-TAM"], [r["MG-TAM"], r["MO-TAM"]], color="#999999", lw=1)
ax.scatter([0]*len(d), d["MG-TAM"], color="#4c72b0", s=25, zorder=5)
ax.scatter([1]*len(d), d["MO-TAM"], color="#c44e52", s=25, zorder=5)
ax.set_xticks([0,1]); ax.set_xticklabels(["MG-TAM","MO-TAM"], fontsize=8)
ax.set_ylabel("canonical15 score")
ax.set_title("GSE163120 (n=7): diff=-0.065, P=0.156\n(negative result)", fontsize=8.5)
# B: GSE182109 canonical15 paired
ax = fig.add_subplot(3, 3, 2)
d = pp182[pp182["module"]=="canonical15"].pivot(index="sample", columns="lineage", values="score")
for _, r in d.iterrows():
    ax.plot(["MG-TAM","MO-TAM"], [r["MG-TAM"], r["MO-TAM"]], color="#999999", lw=1)
ax.scatter([0]*len(d), d["MG-TAM"], color="#4c72b0", s=25, zorder=5)
ax.scatter([1]*len(d), d["MO-TAM"], color="#c44e52", s=25, zorder=5)
ax.set_xticks([0,1]); ax.set_xticklabels(["MG-TAM","MO-TAM"], fontsize=8)
ax.set_ylabel("canonical15 score")
ax.set_title("GSE182109 (n=11): diff=+0.010, P=0.240\n(negative result)", fontsize=8.5)
# C: meta forest
ax = fig.add_subplot(3, 3, 3)
m2 = meta.sort_values("pooled_diff")
for i, (_, r) in enumerate(m2.iterrows()):
    se = r["se"]; mu = r["pooled_diff"]
    col = "#888888" if r["meta_p"] >= 0.05 else "#c44e52"
    ax.errorbar(mu, i, xerr=1.96*se, fmt="D", color=col, ms=6, capsize=4)
    ax.text(0.03, i, f"P={r['meta_p']:.3f}", fontsize=7.5, va="center")
ax.axvline(0, color="grey", lw=0.8, ls="--")
ax.set_yticks(range(len(m2))); ax.set_yticklabels(m2["set"], fontsize=8)
ax.set_xlabel("Pooled MO-TAM minus MG-TAM diff (95% CI)")
ax.set_title("Random-effects meta-analysis", fontsize=8.5)
# D: per-patient diff distribution (canonical15 vs context53)
ax = fig.add_subplot(3, 3, 4)
allp = pd.concat([pp163.assign(cohort="GSE163120"), pp182.assign(cohort="GSE182109")], ignore_index=True)
diffp = allp.pivot_table(index=["cohort","sample","module"], columns="lineage", values="score")
diffp["diff"] = diffp["MO-TAM"] - diffp["MG-TAM"]
diffp = diffp.reset_index()
xpos = {}
xi = 0
for cohort in ["GSE163120","GSE182109"]:
    for mod in ["canonical15","context53"]:
        xpos[(cohort,mod)] = xi; xi += 1
for (cohort, mod), x in xpos.items():
    sub = diffp[(diffp["cohort"]==cohort)&(diffp["module"]==mod)]
    xs = np.random.normal(x, 0.09, len(sub))
    col = "#4c72b0" if mod=="canonical15" else "#c44e52"
    ax.scatter(xs, sub["diff"], s=14, color=col, alpha=0.8)
    ax.plot([x-0.28, x+0.28], [sub["diff"].mean()]*2, color="black", lw=1.2)
ax.axhline(0, color="grey", lw=0.8, ls="--")
ax.set_xticks(list(xpos.values()))
ax.set_xticklabels([f"{c}\n{m[:7]}" for c,m in xpos], fontsize=7)
ax.set_ylabel("Per-patient MO-MG diff"); ax.set_title("Per-patient module differences", fontsize=8.5)
# E: gene contribution + leave-one-gene
ax = fig.add_subplot(3, 3, 5)
top = gc.sort_values("diff_MO_minus_MG", ascending=False).head(6)
ax.barh(top["symbol"][::-1], top["diff_MO_minus_MG"][::-1], color="#c44e52")
ax.axvline(0, color="grey", lw=0.7, ls="--")
ax.set_xlabel("MO-MG contribution"); ax.set_title("Top contributors (GSE163120)", fontsize=8.5)
ax = fig.add_subplot(3, 3, 6)
logo_t = logo.sort_values("delta_vs_full", ascending=False).head(6)
ax.barh(logo_t["removed_gene"][::-1], logo_t["delta_vs_full"][::-1], color="#4c72b0")
ax.set_xlabel("Δ diff upon gene removal"); ax.set_title("Leave-one-gene-out (top)", fontsize=8.5)
# F: within-MO-TAM exploratory scatter
ax = fig.add_subplot(3, 3, 7)
ax.scatter(cs["GMCSF_rep_nonIFN_14_mean"], cs["mod15_mean"], s=0.6, alpha=0.12, color="#c44e52")
ax.set_xlabel("GM-CSF-repressed sig. (14, non-IFN)"); ax.set_ylabel("canonical15 module")
ax.set_title("Within-MO-TAM (GSE163120, n=4,652 cells)\npartial rho=0.585 — EXPLORATORY", fontsize=8.5, color="#8a6d1a")
# panels G,H: context53 paired (supplement meta context)
ax = fig.add_subplot(3, 3, 8)
d = pp163[pp163["module"]=="context53"].pivot(index="sample", columns="lineage", values="score")
for _, r in d.iterrows():
    ax.plot(["MG-TAM","MO-TAM"], [r["MG-TAM"], r["MO-TAM"]], color="#999999", lw=1)
ax.scatter([0]*len(d), d["MG-TAM"], color="#4c72b0", s=25, zorder=5)
ax.scatter([1]*len(d), d["MO-TAM"], color="#c44e52", s=25, zorder=5)
ax.set_xticks([0,1]); ax.set_xticklabels(["MG-TAM","MO-TAM"], fontsize=8)
ax.set_ylabel("context53 score")
ax.set_title("GSE163120 context53: diff=+0.118, P=0.047\n(exploratory)", fontsize=8.5)
ax = fig.add_subplot(3, 3, 9)
d = pp182[pp182["module"]=="context53"].pivot(index="sample", columns="lineage", values="score")
for _, r in d.iterrows():
    ax.plot(["MG-TAM","MO-TAM"], [r["MG-TAM"], r["MO-TAM"]], color="#999999", lw=1)
ax.scatter([0]*len(d), d["MG-TAM"], color="#4c72b0", s=25, zorder=5)
ax.scatter([1]*len(d), d["MO-TAM"], color="#c44e52", s=25, zorder=5)
ax.set_xticks([0,1]); ax.set_xticklabels(["MG-TAM","MO-TAM"], fontsize=8)
ax.set_ylabel("context53 score")
ax.set_title("GSE182109 context53: diff=+0.046, P=0.005\n(exploratory)", fontsize=8.5)
fig.suptitle("Figure 5. Patient-level extrapolation boundaries (negative and exploratory results shown together)", fontsize=12)
save(fig, "Figure5_patient_boundaries")
pp163.to_csv(os.path.join(SRC, "fig5_paired_GSE163120.csv"), index=False)
pp182.to_csv(os.path.join(SRC, "fig5_paired_GSE182109.csv"), index=False)
meta.to_csv(os.path.join(SRC, "fig5_meta.csv"), index=False)
gc.to_csv(os.path.join(SRC, "fig5_gene_contribution.csv"), index=False)
logo.to_csv(os.path.join(SRC, "fig5_leave_one_gene.csv"), index=False)
cs[["GMCSF_rep_nonIFN_14_mean","mod15_mean","sample"]].to_csv(os.path.join(SRC, "fig5_within_motam.csv"), index=False)
print("Fig5 done. All main figures saved.")
