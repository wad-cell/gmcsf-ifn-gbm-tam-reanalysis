# -*- coding: utf-8 -*-
"""C2 prep: 计算两队列患者级 MO/MG 配对模块得分（canonical15/context53/Hallmark IFNa/IFNg）
用于 Fig5 配对图与 Figure source data。逻辑与 B7 模块三完全一致。
"""
import os, pickle, json
import numpy as np, pandas as pd


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
OUT = os.path.join(WS, "output", "B7")
TMP = os.path.join(WS, "temp", "B7")
C2T = os.path.join(WS, "temp", "C2")
os.makedirs(C2T, exist_ok=True)

ann = pd.read_csv(os.path.join(OUT, "core68_reannotation.csv"))
mod15 = ann[ann["module"]=="conserved_canonical_IFN_15"]["symbol"].tolist()
mod53 = ann[ann["module"]=="context_dependent_53"]["symbol"].tolist()
with open(os.path.join(TMP, "hallmark_sets.json")) as f:
    hs = json.load(f)
ifna = hs["IFN_ALPHA"]; ifng = hs["IFN_GAMMA"]
modules = {"canonical15": mod15, "context53": mod53,
           "Hallmark_IFNa": ifna, "Hallmark_IFNg": ifng}

def patient_paired(expr, annot, genes, sample_col):
    g = [x for x in genes if x in expr.columns]
    rows = []
    for sample in sorted(annot[sample_col].unique()):
        for lin in ["MO-TAM", "MG-TAM"]:
            cells = annot[(annot[sample_col]==sample) & (annot["lineage"]==lin)].index
            cells = [c for c in cells if c in expr.index]
            if len(cells) < 3:
                continue
            rows.append({"sample": str(sample), "lineage": lin, "n": len(cells),
                         "n_genes_used": len(g),
                         "score": float(expr.loc[cells, g].mean().mean())})
    return pd.DataFrame(rows)

# ---- GSE163120 ----
with open(os.path.join(TMP, "gse163120_tam_expr.pkl"), "rb") as f:
    e163 = pickle.load(f)
a163 = pd.read_csv(os.path.join(TMP, "gse163120_tam_annot.csv")).set_index("Unnamed: 0")
frames163 = []
for name, genes in modules.items():
    d = patient_paired(e163, a163, genes, "sample")
    d = d.assign(module=name, cohort="GSE163120")
    frames163.append(d)
pp163 = pd.concat(frames163, ignore_index=True)
pp163.to_csv(os.path.join(C2T, "patient_paired_GSE163120.csv"), index=False)
print("GSE163120 paired OK:", pp163.shape)
print(pp163[pp163["module"]=="canonical15"].pivot(index="sample",columns="lineage",values="score").round(4).to_string())

# ---- GSE182109 ----
with open(os.path.join(TMP, "gse182109_myeloid_expr_aligned.pkl"), "rb") as f:
    e182 = pickle.load(f)
a182 = pd.read_csv(os.path.join(TMP, "gse182109_myeloid_annot.csv")).set_index("barcode")
frames182 = []
for name, genes in modules.items():
    d = patient_paired(e182, a182, genes, "patient")
    d = d.assign(module=name, cohort="GSE182109")
    frames182.append(d)
pp182 = pd.concat(frames182, ignore_index=True)
pp182.to_csv(os.path.join(C2T, "patient_paired_GSE182109.csv"), index=False)
print("GSE182109 paired OK:", pp182.shape)
print(pp182[pp182["module"]=="canonical15"].pivot(index="sample",columns="lineage",values="score").round(4).to_string())
