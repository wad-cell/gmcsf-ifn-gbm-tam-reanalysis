# -*- coding: utf-8 -*-
"""Generate Figure5_v3_1 (A-H panels; within-MO removed per C3.2) for manuscript_NOA_v3_1."""
import os
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output", "C3.1")
OUT = os.path.join(BASE, "figures")
os.makedirs(OUT, exist_ok=True)

scores_182 = pd.read_csv(os.path.join(BASE, "A4_module_scores_patientlevel.csv"))
scores_163 = pd.read_csv(os.path.join(BASE, "A5_GSE163120_patient_module_scores.csv"))
corr = pd.read_csv(os.path.join(BASE, "A5_1_one_stage_corrected.csv"))
spec = pd.read_csv(os.path.join(BASE, "A5_1_GMCSF_specificity_audit.csv"))
canon = pd.read_csv(os.path.join(BASE, "A5_1_canonical15_robustness.csv"))

def load_abs(cohort, module):
    out = []
    if cohort == "GSE163120":
        d = scores_163[scores_163["module"] == module][["patient", "mo_mg", "score"]]
    else:
        d = scores_182[scores_182["patient"].str.startswith("ndGBM")][["patient", "mo_mg", module]].rename(columns={module: "score"})
    for p in d["patient"].unique():
        sub = d[d["patient"] == p]
        mo = sub.loc[sub["mo_mg"] == "MO-TAM", "score"]
        mg = sub.loc[sub["mo_mg"] == "MG-TAM", "score"]
        if len(mo) == 1 and len(mg) == 1:
            out.append((str(p), float(mg.iloc[0]), float(mo.iloc[0])))
    return out

def paired_panel(ax, title, module, ylim):
    cohorts = [load_abs("GSE163120", module), load_abs("GSE182109", module)]
    for j, cohort in enumerate(cohorts):
        base = j * 2.6
        for k, (p, mg, mo) in enumerate(cohort):
            y = len(cohort) - k - 1
            ax.plot([base - 0.18, base + 0.18], [mg, mo], color="#666666", lw=0.8, zorder=1)
            ax.scatter(base - 0.18, mg, s=22, color="#66c2a5", edgecolor="k", lw=0.4, zorder=2)
            ax.scatter(base + 0.18, mo, s=22, color="#fc8d62", edgecolor="k", lw=0.4, zorder=2)
        ax.axvline(base, color="grey", lw=0.5, ls=":")
    ax.set_xticks([-0.18, 0.18, 2.6 - 0.18, 2.6 + 0.18])
    ax.set_xticklabels(["MG", "MO", "MG", "MO"], fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.set_ylabel("module score (log2CPM)", fontsize=8)
    if ylim:
        ax.set_ylim(*ylim)
    ym = ax.get_ylim()[1]
    ax.text(0.0, ym * 1.06, "GSE163120 (n=7)", ha="center", fontsize=7.5, style="italic")
    ax.text(2.6, ym * 1.06, "GSE182109 (n=11)", ha="center", fontsize=7.5, style="italic")

def box(ax, x, y, w, h, text, fc="#eaf2f8", fs=7.5):
    ax.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec="k", lw=0.7))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs)
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)

fig = plt.figure(figsize=(19, 10.5))
gs = fig.add_gridspec(2, 4, hspace=0.55, wspace=0.32)

# ============ Panel A: cohort flow ============
axA = fig.add_subplot(gs[0, 0])
axA.axis("off")
axA.set_title("A. Cohort inclusion and patient-level flow", fontsize=10, loc="left", fontweight="bold")
box(axA, 0.2, 8.2, 4.6, 1.5, "GSE163120\nprimary GBM, 7 patients\n(ND1-ND7; 21,303 cells)\n[Pombo Antunes 2021]", "#dbe9f6")
box(axA, 5.2, 8.2, 4.6, 1.5, "GSE182109\nfull-transcriptome reconstruction\n201,893 cells x 37,407 genes\n[Abdelfattah 2022]", "#dbe9f6")
box(axA, 0.2, 6.0, 4.6, 1.6, "MO-TAM / MG-TAM labels\nfrom source annotations\n(patient-level)", "#f2f7fb")
box(axA, 5.2, 6.0, 4.6, 1.6, "11 ndGBM main analysis\n5 rGBM exploratory only\n(MO/MG = operational mapping)", "#f2f7fb")
box(axA, 0.2, 3.9, 4.6, 1.6, "Pseudobulk module scores\nper patient per lineage\nmean log2(CPM+1)", "#fff7ec")
box(axA, 5.2, 3.9, 4.6, 1.6, "Paired difference\nMO-TAM - MG-TAM\n(statistical unit = patient)", "#fff7ec")
box(axA, 0.2, 1.5, 9.6, 1.8, "One-stage patient-level analysis (18 patients, HC3 SE,\nproportion-centered cohort; intercept = overall patient mean;\ncohort-stratified reported in full; two-stage meta as sensitivity;\nBH-FDR over 13 modules)", "#fde8e8")
for (x1,y1,x2,y2) in [(2.5,8.2,2.5,7.6),(7.5,8.2,7.5,7.6),(2.5,6.0,2.5,5.5),(7.5,6.0,7.5,5.5),(2.5,3.9,2.5,3.3),(7.5,3.9,7.5,3.3)]:
    axA.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=9,lw=0.8,color="k"))

# ============ Panel B: GM-CSF-response paired ============
axB = fig.add_subplot(gs[0, 1])
paired_panel(axB, "B. GM-CSF-response\n(MO-TAM vs MG-TAM)", "GM_CSF_response", (3.0, 5.4))

# ============ Panel C: canonical15 paired ============
axC = fig.add_subplot(gs[0, 2])
paired_panel(axC, "C. canonical15 IFN module\n(MO-TAM vs MG-TAM)", "canonical15", (3.6, 6.2))

# ============ Panel D: one-stage forest (C3.2 corrected) ============
axD = fig.add_subplot(gs[0, 3])
axD.set_title("D. One-stage patient-level analysis (18 patients)", fontsize=10, loc="left", fontweight="bold")
m = corr.sort_values("overall_mean_paired_diff_18", ascending=False)
labels = list(m["module"])
ypos = np.arange(len(labels))
ests = m["overall_mean_paired_diff_18"].values
lo = m["overall_HC3_95CI_low"].values
hi = m["overall_HC3_95CI_high"].values
colors = ["#d62728" if x == "GM_CSF_response" else "#1f77b4" for x in labels]
axD.errorbar(ests, ypos, xerr=[ests - lo, hi - ests], fmt="o", ms=4.5, color="k", ecolor="grey", elinewidth=0.8, capsize=2)
axD.scatter(ests, ypos, s=22, c=colors, zorder=3)
axD.axvline(0, color="grey", lw=0.8, ls="--")
axD.set_yticks(ypos); axD.set_yticklabels(labels, fontsize=7)
axD.set_xlabel("overall mean paired difference (MO-MG)", fontsize=8)
axD.text(0.985, 0.02, "red = GM_CSF_response", transform=axD.transAxes, fontsize=7, color="#d62728", ha="right")
axD.text(0.985, 0.055, "HC3 95% CI; exact 2^18 sign-flip P,\nBH-FDR over 13 modules", transform=axD.transAxes, fontsize=6.5, color="grey", ha="right")

# ============ Panel E: GM-CSF overlap-removal sensitivity ============
axE = fig.add_subplot(gs[1, 0])
axE.set_title("E. GM-CSF-response overlap-removal sensitivity", fontsize=10, loc="left", fontweight="bold")
versions = ["A full", "B -IFN", "C -hypoxia", "D -lineage", "E -all overlap"]
vals = [float(spec.loc[spec["item"] == k, "value"].iloc[0]) for k in
        ["A_full_pooled", "B_remove_IFN_genes_pooled", "C_remove_hypoxia_genes_pooled",
         "D_remove_lineage_gene_pooled", "E_remove_all_overlap_pooled"]]
xx = np.arange(len(versions))
axE.bar(xx, vals, width=0.6, color=["#2b8cbe"]*4 + ["#d62728"], alpha=0.85, edgecolor="k", lw=0.5)
for i, v in enumerate(vals):
    axE.text(i, v + 0.006, f"{v:.3f}", ha="center", fontsize=8)
nm = float(spec.loc[spec["item"] == "null_permutation_median", "value"].iloc[0])
axE.axhline(nm, color="grey", ls="--", lw=0.8)
axE.text(len(versions)-0.4, nm + 0.008, f"null median {nm:.4f}\n(1000 matched random sets)", fontsize=6.5, color="grey")
axE.set_xticks(xx); axE.set_xticklabels(versions, fontsize=7.5, rotation=20, ha="right")
axE.set_ylabel("pooled paired difference (MO-MG)", fontsize=8)
axE.set_ylim(0, 0.46)
axE.text(0.02, 0.44, "direction: 7/7 and 10/11 patients MO>MG\npermutation P<0.001", fontsize=7)

# ============ Panel F: canonical15 robustness ============
axF = fig.add_subplot(gs[1, 1])
axF.set_title("F. canonical15 robustness", fontsize=10, loc="left", fontweight="bold")
full = -0.2141
logo = [-0.2635, -0.1787]
loop = [-0.2731, -0.1721]
cats = ["full pooled", "LOGO range (15)", "LOOPO range (18)"]
x2 = np.arange(3)
axF.bar(x2, [full, (logo[0]+logo[1])/2, (loop[0]+loop[1])/2], width=0.5, color="#1f77b4", alpha=0.85, edgecolor="k", lw=0.5)
axF.errorbar([1, 2], [(logo[0]+logo[1])/2, (loop[0]+loop[1])/2],
             yerr=[(logo[1]-logo[0])/2, (loop[1]-loop[0])/2], fmt="none", ecolor="k", capsize=4, lw=1)
axF.text(0, full - 0.015, f"{full:.3f}", ha="center", fontsize=8)
axF.axhline(0, color="grey", lw=0.8)
axF.set_xticks(x2); axF.set_xticklabels(cats, fontsize=7.5)
axF.set_ylabel("pooled paired difference (MO-MG)", fontsize=8)
axF.set_ylim(-0.35, 0.02)
axF.text(0.5, -0.30, "all LOGO<0, all LOOPO<0\nno detection-rate bias (86-100%)\n15/15 genes scored in both cohorts", fontsize=7)

# ============ Panel G: context53 (exploratory; corrected labels) ============
axG = fig.add_subplot(gs[1, 2])
axG.set_title("G. context53 (exploratory decomposition)", fontsize=10, loc="left", fontweight="bold")
def diff_series(cohort, module):
    return [float(mo - mg) for (p, mg, mo) in load_abs(cohort, module)]
g1 = diff_series("GSE163120", "context53"); g2 = diff_series("GSE182109", "context53")
axG.scatter(np.zeros(len(g1)) + 0.15, g1, s=25, color="#2b8cbe", edgecolor="k", lw=0.4, zorder=3)
axG.scatter(np.ones(len(g2)) + 0.15, g2, s=25, color="#fc8d62", edgecolor="k", lw=0.4, zorder=3)
for j, (g, off) in enumerate([(g1, 0.15), (g2, 1.15)]):
    mm = np.mean(g); ss = np.std(g, ddof=1) / math.sqrt(len(g))
    axG.errorbar(off, mm, yerr=1.96*ss, fmt="o", color="k", ms=5, capsize=3, zorder=4)
    axG.text(off, mm + 0.06, f"{mm:+.3f}", ha="center", fontsize=7)
axG.axhline(0, color="grey", ls="--", lw=0.8)
axG.set_xticks([0.15, 1.15]); axG.set_xticklabels(["GSE163120 (n=7)", "GSE182109 (n=11)"], fontsize=8)
axG.set_ylabel("patient-level paired difference", fontsize=8)
axG.text(0.02, 0.97, "overall +0.131; exact raw P=0.023\nFDR=0.059 (not significant after correction)",
         transform=axG.transAxes, fontsize=7, va="top", ha="left", color="#d62728")

# ============ Panel H: evidence boundary diagram (was I) ============
axH = fig.add_subplot(gs[1, 3])
axH.axis("off")
axH.set_title("H. Evidence boundary (in vitro to patient)", fontsize=10, loc="left", fontweight="bold")
box(axH, 0.2, 8.0, 4.6, 1.6, "Level 1 - in vitro perturbation\nco-culture suppression; exogenous\nGM-CSF partial replicate;\nanti-GM-CSF partial rescue", "#d9f0d3")
box(axH, 5.2, 8.0, 4.6, 1.6, "Level 2 - patient association\nGM-CSF-responsive enrichment in\noperationally mapped MO-TAMs;\ncanonical15 lower (C3.2-adjudicated)", "#c6e0f5")
box(axH, 0.2, 5.4, 4.6, 1.9, "Exploratory\ncontext53 (FDR=0.059);\nHallmark_hypoxia (I2=88.6%);\nrGBM (n=5); within-MO-TAM (S8)", "#fee0d2")
box(axH, 5.2, 5.4, 4.6, 1.9, "Caveats (written in text)\nobservational association; operational\ncell-origin mapping; only two cohorts;\nno GM-CSF protein measurement;\nhypoxia & lineage remain alternatives;\nno causality", "#fff7bc")
box(axH, 0.2, 2.6, 9.6, 1.9, "Explicitly NOT claimed\n- globally suppressed type I IFN signaling in patients\n- patient MO-TAMs suppressed by GM-CSF (causality)\n- in vivo validation; direct GM-CSF pathway activation (STAT5/CISH no support)", "#fde8e8")
axH.set_xlim(0, 10); axH.set_ylim(0, 10)

for ax in [axB, axC, axD, axE, axF, axG]:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)

fig.suptitle("Figure 5. Patient-level correlates of the in vitro GM-CSF-associated IFN program (C3.2-corrected, v3_1)",
             fontsize=12, fontweight="bold", y=0.99)
fig.savefig(os.path.join(OUT, "Figure5_v3_1.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure5_v3_1.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure5_v3_1.tiff"), dpi=600, bbox_inches="tight")
print("saved:", [f for f in os.listdir(OUT) if "v3_1" in f])
