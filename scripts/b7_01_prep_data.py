# -*- coding: utf-8 -*-
"""B7 数据准备：加载 GSE163120 与 GSE182109 表达矩阵 + 注释，保存中间产物"""
import pandas as pd, numpy as np, os, gzip, json

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
TMP = os.path.join(WS, "temp", "B7")
os.makedirs(TMP, exist_ok=True)

# ---------- GSE163120 ----------
d = os.path.join(WS, "temp", "GSM4972211")
an = pd.read_csv(os.path.join(d, "annot.csv.gz"))
# 读取矩阵（基因 x 细胞）
mat = pd.read_csv(os.path.join(d, "matrix.csv.gz"), index_col=0)
# mat: 基因(行) x 细胞(列)
mat = mat.T  # 细胞 x 基因
print("GSE163120 mat:", mat.shape)
# 对齐 annot
an = an.set_index("cell")
common = an.index.intersection(mat.index)
mat = mat.loc[common]
an = an.loc[common]
print("GSE163120 aligned:", mat.shape)
# 定义 lineage: TAM 1 = MO-TAM, TAM 2 = MG-TAM
an["lineage"] = an["cluster"].map({"TAM 1": "MO-TAM", "TAM 2": "MG-TAM"})
# 只保留 TAM
tam = an[an["lineage"].isin(["MO-TAM", "MG-TAM"])]
mat_tam = mat.loc[tam.index]
print("GSE163120 TAM cells:", mat_tam.shape)
print("lineage counts:", tam["lineage"].value_counts().to_dict())
print("sample counts:", tam["sample"].value_counts().to_dict())

# 保存
mat_tam.to_pickle(os.path.join(TMP, "gse163120_tam_expr.pkl"))
tam[["cluster", "sample", "lineage"]].to_csv(os.path.join(TMP, "gse163120_tam_annot.csv"))
# 全转录组基因名
all_genes = list(mat.columns)
json.dump(all_genes, open(os.path.join(TMP, "gse163120_all_genes.json"), "w"))

# ---------- GSE182109 ----------
g = os.path.join(WS, "temp", "GSE182109")
go = pd.read_csv(os.path.join(g, "ndGBM_gene_order.txt"), header=None)
genes182 = go[0].astype(str).tolist()
expr = np.load(os.path.join(g, "ndGBM_interest_expr.npy"))  # (164, 132209)
meta = pd.read_csv(os.path.join(g, "ndGBM_meta.csv"))
ma = pd.read_csv(os.path.join(WS, "output", "B6", "GSE182109_myeloid_annot.csv"))
# 对齐：ma 有 barcode/patient/lineage
ma = ma.set_index("barcode")
# 表达矩阵转置为 细胞 x 基因
mat182 = pd.DataFrame(expr.T, index=meta["barcode"], columns=genes182)
# 只保留髓系
common182 = ma.index.intersection(mat182.index)
mat182 = mat182.loc[common182]
ma = ma.loc[common182]
print("\nGSE182109 myeloid mat:", mat182.shape)
print("lineage counts:", ma["lineage"].value_counts().to_dict())
print("patient counts:", ma["patient"].value_counts().to_dict())

mat182.to_pickle(os.path.join(TMP, "gse182109_myeloid_expr.pkl"))
ma[["sample", "patient", "lineage"]].to_csv(os.path.join(TMP, "gse182109_myeloid_annot.csv"))
print("DONE")
