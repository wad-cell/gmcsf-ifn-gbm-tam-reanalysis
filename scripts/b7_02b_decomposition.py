# -*- coding: utf-8 -*-
"""B7 模块二补全：core68 在 MO-TAM 偏高来源拆解
1. signature_decomposition.csv: 68/15/53 患者内差值分解
2. gene_contribution_analysis.csv: 基因贡献
3. leave_one_gene_out_results.csv: leave-one-gene-out
4. signature_method_sensitivity.pdf: 4 种打分法比较
"""
import os, pickle
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
OUT = os.path.join(WS, "output", "B7")
TMP = os.path.join(WS, "temp", "B7")

ann = pd.read_csv(os.path.join(OUT, "core68_reannotation.csv"))
mod15 = ann[ann["module"]=="conserved_canonical_IFN_15"]["symbol"].tolist()
mod53 = ann[ann["module"]=="context_dependent_53"]["symbol"].tolist()
mod68 = mod15 + mod53

with open(os.path.join(TMP, "gse163120_tam_expr.pkl"), "rb") as f:
    e163 = pickle.load(f)
a163 = pd.read_csv(os.path.join(TMP, "gse163120_tam_annot.csv")).set_index("Unnamed: 0")

def patient_pseudobulk(expr, annot, genes):
    """患者×lineage pseudobulk 平均表达"""
    g = [x for x in genes if x in expr.columns]
    rows = []
    for sample in sorted(annot["sample"].unique()):
        for lin in ["MO-TAM", "MG-TAM"]:
            cells = annot[(annot["sample"]==sample) & (annot["lineage"]==lin)].index
            cells = [c for c in cells if c in expr.index]
            if len(cells) < 3:
                continue
            rows.append({"sample": sample, "lineage": lin, "n": len(cells),
                         "score": expr.loc[cells, g].mean().mean()})
    df = pd.DataFrame(rows)
    piv = df.pivot(index="sample", columns="lineage", values="score").dropna()
    return piv

# ============ 1. signature_decomposition ============
rows = []
for name, genes in [("core68", mod68), ("canonical15", mod15), ("context53", mod53)]:
    piv = patient_pseudobulk(e163, a163, genes)
    diff = piv["MO-TAM"] - piv["MG-TAM"]
    rows.append({"module": name, "n_genes": len(genes), "n_patients": len(piv),
                 "mean_diff_MO_MG": diff.mean(), "median_diff": diff.median(),
                 "n_MO_gt_MG": int((diff>0).sum()),
                 "wilcoxon_p": wilcoxon(diff).pvalue if len(piv)>=3 else np.nan})
dec = pd.DataFrame(rows)
dec.to_csv(os.path.join(OUT, "signature_decomposition.csv"), index=False)
print("=== signature_decomposition ===")
print(dec.to_string(index=False, float_format=lambda x:f"{x:.4f}"))

# ============ 2. gene_contribution ============
# 每个基因对 MO-MG 差值的贡献（患者内平均表达差）
g_rows = []
for g in mod68:
    if g not in e163.columns:
        continue
    g_rows.append({"symbol": g, "module": "canonical15" if g in mod15 else "context53",
                   "MO_mean": e163.loc[a163[(a163["lineage"]=="MO-TAM")].index, g].mean(),
                   "MG_mean": e163.loc[a163[(a163["lineage"]=="MG-TAM")].index, g].mean()})
gc = pd.DataFrame(g_rows)
gc["diff_MO_minus_MG"] = gc["MO_mean"] - gc["MG_mean"]
gc = gc.sort_values("diff_MO_minus_MG", ascending=False)
gc.to_csv(os.path.join(OUT, "gene_contribution_analysis.csv"), index=False)
print("\n=== gene_contribution top10 (MO-MG 最大) ===")
print(gc.head(10).to_string(index=False, float_format=lambda x:f"{x:.4f}"))
print("=== gene_contribution bottom10 (MG-MO 最大) ===")
print(gc.tail(10).to_string(index=False, float_format=lambda x:f"{x:.4f}"))

# ============ 3. leave_one_gene_out ============
base_piv = patient_pseudobulk(e163, a163, mod68)
base_diff = (base_piv["MO-TAM"] - base_piv["MG-TAM"]).mean()
lo_rows = []
for g in mod68:
    if g not in e163.columns:
        continue
    rest = [x for x in mod68 if x != g]
    piv = patient_pseudobulk(e163, a163, rest)
    d = (piv["MO-TAM"] - piv["MG-TAM"]).mean()
    lo_rows.append({"removed_gene": g, "module": "canonical15" if g in mod15 else "context53",
                    "diff_without_gene": d, "delta_vs_full": d - base_diff})
lo = pd.DataFrame(lo_rows).sort_values("delta_vs_full", ascending=False)
lo.to_csv(os.path.join(OUT, "leave_one_gene_out_results.csv"), index=False)
print(f"\n=== leave_one_gene_out (full68 diff={base_diff:.4f}) ===")
print("去掉后差值下降最多（该基因驱动 MO 高）:")
print(lo.head(8).to_string(index=False, float_format=lambda x:f"{x:.4f}"))
print("去掉后差值上升最多（该基因驱动 MG 高）:")
print(lo.tail(8).to_string(index=False, float_format=lambda x:f"{x:.4f}"))

# ============ 4. signature_method_sensitivity.pdf ============
# 4 种打分法：mean / ucell / aucell / pseudobulk_gsva
import sys
sys.path.insert(0, TMP)
from b7_scoring import mean_score, ucell_score, aucell_score

mo_cells = a163[a163["lineage"]=="MO-TAM"].index
mg_cells = a163[a163["lineage"]=="MG-TAM"].index
methods = {}
for name, genes in [("core68", mod68), ("canonical15", mod15), ("context53", mod53)]:
    g = [x for x in genes if x in e163.columns]
    methods[name] = {
        "mean": mean_score(e163, g),
        "ucell": ucell_score(e163, g),
        "aucell": aucell_score(e163, g),
    }
# pseudobulk_gsva 用 gseapy
try:
    import gseapy as gp
    def gsva_score(expr, genes):
        g = [x for x in genes if x in expr.columns]
        sub = expr[g].T
        gs = gp.gsva(sub, gene_sets={name: g}, method="gsva", verbose=False,
                     outdir=None, min_size=1, max_size=5000)
        return gs.loc["gsva"].values
    for name in ["core68", "canonical15", "context53"]:
        methods[name]["gsva"] = gsva_score(e163, methods[name]["mean"].index if False else e163)
except Exception as ex:
    print("gsva 失败:", ex)

# 绘制 4 种打分法 MO vs MG 分布
fig, axes = plt.subplots(3, 4, figsize=(16, 10))
for i, name in enumerate(["core68", "canonical15", "context53"]):
    for j, m in enumerate(["mean", "ucell", "aucell"]):
        ax = axes[i, j]
        s = methods[name][m]
        if isinstance(s, pd.Series):
            s = s.reindex(e163.index)
        mo = s[mo_cells]; mg = s[mg_cells]
        ax.hist(mo, bins=40, alpha=0.5, label="MO-TAM", color="tab:red")
        ax.hist(mg, bins=40, alpha=0.5, label="MG-TAM", color="tab:blue")
        ax.set_title(f"{name} - {m}")
        ax.legend(fontsize=7)
fig.suptitle("Signature scoring method sensitivity (GSE163120)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "signature_method_sensitivity.pdf"))
fig.savefig(os.path.join(OUT, "signature_method_sensitivity.png"), dpi=150)
plt.close(fig)
print("\n已保存 signature_method_sensitivity.pdf/png")
