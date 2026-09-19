# -*- coding: utf-8 -*-
"""C2 正式制图：Supplementary Figures S1-S7
TIFF 600dpi + PNG 预览 + figure_source_data.csv
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
WS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
O = os.path.join(WS, "output")
B5 = os.path.join(O, "B5"); B7 = os.path.join(O, "B7"); B6 = os.path.join(O, "B6")
FIG = os.path.join(O, "C2", "figures"); SRC = os.path.join(O, "C2", "figure_source_data")
os.makedirs(FIG, exist_ok=True); os.makedirs(SRC, exist_ok=True)

def rd(p, b):
    return pd.read_csv(os.path.join(b, p))

def save(fig, name):
    fig.savefig(os.path.join(FIG, f"{name}.tif"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, f"{name}.png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, f"{name}.pdf"), bbox_inches="tight")
    plt.close(fig)

# ============ Figure S1. Permutation schemes ============
ra = rd("robustness_analysis.csv", B5).set_index("test")
perm = rd("permutation_distributions.csv", B5)
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
# left: two statistics observed vs permutation interval
stats = [
    ("frac_negative", "Fraction of\nsuppressed genes", 1.0, 0.450, 0.445, 0.529, "P_emp=0.0001"),
    ("mean_log2FC", "Mean log2FC", -1.878, 0.407, 0.062, 0.422, "P_emp=0.0001"),
]
for i, (_, lab, obs, pm, plo, phi, ptxt) in enumerate(stats):
    axes[0].errorbar(pm, i, xerr=[[pm-plo],[phi-pm]], fmt="o", color="#888888", ms=7, capsize=5)
    axes[0].plot([obs], [i], "rD", ms=9)
    axes[0].text(obs, i+0.06, f"observed {obs:.3f}\n{ptxt}", fontsize=7.5, color="#c44e52", ha="center")
axes[0].axvline(0, color="grey", lw=0.7, ls="--")
axes[0].set_yticks(range(len(stats))); axes[0].set_yticklabels([s[1] for s in stats], fontsize=8)
axes[0].set_xlim(-2.4, 1.6); axes[0].set_xlabel("Statistic value (permutation interval vs observed)")
axes[0].set_title("Matched-permutation schemes (A/B/C):\nobserved far outside null; P_emp=0.0001 all schemes;\nscheme C random mean=-0.820", fontsize=8)
axes[1].hist(perm["perm_rho"], bins=60, color="#8c8c8c", alpha=0.8, edgecolor="white")
axes[1].axvline(0.581, color="#c44e52", lw=1.8)
axes[1].set_xlabel("Permutation rho (GM-CSF vs LN229)"); axes[1].set_ylabel("Count (N=10,000)")
axes[1].set_title("Permutation null (observed rho=0.581; P_emp=0.0173)", fontsize=8)
fig.suptitle("Figure S1. Matched-permutation schemes for co-culture suppression", fontsize=12)
save(fig, "FigureS1_permutation")
pd.DataFrame({"stat": [s[0] for s in stats], "observed":[s[2] for s in stats],
              "perm_mean":[s[3] for s in stats], "perm_lo":[s[4] for s in stats],
              "perm_hi":[s[5] for s in stats]}).to_csv(os.path.join(SRC, "figS1_permutation.csv"), index=False)
print("S1 done")

# ============ Figure S2. Leave-one-donor-out ============
loo = rd("robustness_leave_one_donor_out.csv", B5)
loo["cohort"] = loo["label"].str.split("_").str[0]
loo["contrast"] = loo["label"].str.extract(r"(C\d+|B2)$")[0]
loo["donor"] = loo["label"].str.extract(r"donor_(\d+)")[0]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
g163 = loo[loo["cohort"]=="GSE309037"]
axes[0].bar(range(len(g163)), g163["mean_log2FC"], color="#4c72b0")
axes[0].set_xticks(range(len(g163)))
axes[0].set_xticklabels([f"D{d}" for d in g163["donor"]], fontsize=8)
axes[0].axhline(-1.878, color="black", ls="--", lw=1); axes[0].text(len(g163)-0.3, -1.86, "full data", fontsize=7, ha="right")
axes[0].set_ylabel("core68 mean log2FC (direct, B2)")
axes[0].set_title("GSE309037 direct co-culture (n=3 donors)", fontsize=9)
data1 = loo[loo["cohort"].isin(["data1","data2"])]
for ci, c in enumerate(["C4","C7","C9","C8"]):
    sub = data1[data1["contrast"]==c].sort_values("donor")
    axes[1].bar(np.arange(len(sub))+ci*6.8, sub["mean_log2FC"], width=0.85, label=f"{c}")
axes[1].set_xticks([]); axes[1].set_ylabel("core68 mean log2FC")
axes[1].legend(fontsize=7, ncol=4, loc="lower right")
axes[1].set_title("data1/data2 (n=6 donors per contrast):\nC4 direct, C7 neutralization, C9 indirect LN229, C8 GM-CSF", fontsize=8)
axes[1].axhline(0, color="grey", lw=0.7)
fig.suptitle("Figure S2. Leave-one-donor-out robustness", fontsize=12)
save(fig, "FigureS2_loo_donor")
loo.to_csv(os.path.join(SRC, "figS2_loo_donor.csv"), index=False)
print("S2 done")

# ============ Figure S3. Leave-one-gene-out ============
logo = rd("leave_one_gene_out_results.csv", B7)
fig, ax = plt.subplots(figsize=(9, 4.5))
logo_s = logo.sort_values("delta_vs_full", ascending=False)
cols = ["#4c72b0" if m=="canonical15" else "#c44e52" for m in logo_s["module"]]
ax.bar(range(len(logo_s)), logo_s["delta_vs_full"], color=cols)
for g in ["NUPR1","NAMPT","ATF3"]:
    row = logo_s[logo_s["removed_gene"]==g]
    if len(row):
        i = logo_s.index.get_loc(row.index[0])
        ax.text(i, row["delta_vs_full"].values[0]+0.001, g, fontsize=7.5, ha="center", rotation=90, va="bottom")
ax.set_xlabel("Gene removed (sorted by delta)"); ax.set_ylabel("Δ (MO-MG diff) upon removal")
ax.set_xticks([]); ax.axhline(0, color="grey", lw=0.7)
ax.set_title("Leave-one-gene-out on GSE163120 MO-TAM minus MG-TAM diff\n(blue=canonical15, red=context53; NUPR1/NAMPT/ATF3 largest)", fontsize=9)
fig.suptitle("Figure S3. Leave-one-gene-out analysis", fontsize=12)
save(fig, "FigureS3_loo_gene")
logo.to_csv(os.path.join(SRC, "figS3_loo_gene.csv"), index=False)
print("S3 done")

# ============ Figure S4. Scoring-method sensitivity ============
mc = rd("B7_m2_method_comparison.csv", B7)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
mods = [68, 15, 53]
x = np.arange(len(mods)); w = 0.25
for i, meth in enumerate(["mean","ucell","aucell"]):
    sub = mc[mc["method"]==meth]
    ax.bar(x+i*w, sub.set_index("module").loc[mods,"diff"], width=w, label=meth)
ax.axhline(0, color="grey", lw=0.7)
ax.set_xticks(x+w); ax.set_xticklabels(["core68 (n=68)","canonical15 (n=15)","context53 (n=53)"], fontsize=8)
ax.set_ylabel("MO-TAM minus MG-TAM module diff (GSE163120)")
ax.legend(fontsize=8); ax.set_title("Module-score method sensitivity (mean / UCell / AUCell)\nNote: direction not fully consistent across methods — method-sensitivity flagged; main conclusions rest on mean-based pseudobulk paired meta-analysis", fontsize=8.5)
fig.suptitle("Figure S4. Scoring-method sensitivity", fontsize=12)
save(fig, "FigureS4_method_sensitivity")
mc.to_csv(os.path.join(SRC, "figS4_method.csv"), index=False)
print("S4 done")

# ============ Figure S5. Rescue-module patient-level ============
rsum = rd("B7_m6_rescue_summary.csv", B7)
fig, ax = plt.subplots(figsize=(9, 4.6))
for i, (_, r) in enumerate(rsum.iterrows()):
    col = "#4c72b0" if r["cohort"]=="GSE163120" else "#c44e52"
    ax.plot([r["mean_diff_MO_minus_MG"]], [i], "o", color=col, ms=8)
    ax.text(r["mean_diff_MO_minus_MG"]+0.008, i, f"P={r['wilcoxon_p']:.3f} ({r['n_MO_gt_MG']}/{r['n_patients']})", fontsize=7, va="center")
ax.axvline(0, color="grey", lw=0.8, ls="--")
ax.set_yticks(range(len(rsum)))
ax.set_yticklabels([f"{r['module']} · {r['cohort']}" for _, r in rsum.iterrows()], fontsize=8)
ax.set_xlabel("MO-TAM minus MG-TAM mean diff (pseudobulk paired)")
ax.set_title("Rescue-module patient-level comparisons — all non-significant,\ndirectionally inconsistent across cohorts (Level 3: not supported)", fontsize=9)
fig.suptitle("Figure S5. Rescue-module patient-level analysis", fontsize=12)
save(fig, "FigureS5_rescue_patient")
rsum.to_csv(os.path.join(SRC, "figS5_rescue_patient.csv"), index=False)
print("S5 done")

# ============ Figure S6. Enrichment ============
ora = rd("enrichment_ORA_full.csv", B5)
gsea = rd("enrichment_GSEA_full.csv", B5)
ora_top = ora.sort_values("FDR").head(15)
gsea_top = gsea.sort_values("FDR q-val").head(12)
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
axes[0].barh(range(len(ora_top)), -np.log10(ora_top["FDR"]+1e-300), color="#55a868")
axes[0].set_yticks(range(len(ora_top)))
axes[0].set_yticklabels([t[:42] for t in ora_top["term"]], fontsize=6.5)
axes[0].set_xlabel("-log10 FDR (ORA)"); axes[0].set_title("Over-representation of core68 (GO: BP)", fontsize=9)
axes[1].barh(range(len(gsea_top)), gsea_top["NES"], color=["#c44e52" if v<0 else "#55a868" for v in gsea_top["NES"]])
axes[1].set_yticks(range(len(gsea_top)))
axes[1].set_yticklabels([t[:34] for t in gsea_top["Term"]], fontsize=6.5)
axes[1].axvline(0, color="grey", lw=0.7)
axes[1].set_xlabel("NES (preranked GSEA, FDR q-val<0.01)"); axes[1].set_title("Hallmark pathway perturbation (direct co-culture)", fontsize=9)
fig.suptitle("Figure S6. Full enrichment results (FDR-corrected)", fontsize=12)
save(fig, "FigureS6_enrichment")
ora.sort_values("FDR").to_csv(os.path.join(SRC, "figS6_ORA_full.csv"), index=False)
gsea.sort_values("FDR q-val").to_csv(os.path.join(SRC, "figS6_GSEA_full.csv"), index=False)
print("S6 done")

# ============ Figure S7. Annotation audit ============
aud = rd("scRNA_annotation_audit.csv", B6).set_index("Unnamed: 0")
lineages = ["MG-TAM","MO-TAM","Monocytes","Other","prol.TAM"]
fig, ax = plt.subplots(figsize=(7.6, 5.2))
im = ax.imshow(aud[lineages].values, aspect="auto", cmap="YlOrRd")
ax.set_xticks(range(len(lineages))); ax.set_xticklabels(lineages, fontsize=8, rotation=30)
ax.set_yticks(range(len(aud))); ax.set_yticklabels(aud.index, fontsize=8)
for i in range(len(aud)):
    for j in range(len(lineages)):
        ax.text(j, i, f"{aud[lineages].values[i,j]:.2f}", ha="center", va="center", fontsize=5.5)
ax.set_title("Lineage-marker mean expression by annotated cell type\n(GSE163120; S100A8/CD14 monocyte axis vs P2RY12/TMEM119 microglial axis)", fontsize=8.5)
cbar = fig.colorbar(im, ax=ax, fraction=0.046); cbar.set_label("mean log-normalized expr.", fontsize=8)
fig.suptitle("Figure S7. Annotation marker audit", fontsize=12)
save(fig, "FigureS7_annotation_audit")
aud.reset_index().to_csv(os.path.join(SRC, "figS7_annotation_audit.csv"), index=False)
print("S7 done. All supplementary figures saved.")
