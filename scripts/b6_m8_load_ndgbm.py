# -*- coding: utf-8 -*-
"""B6 模块八: 独立队列 GSE182109 主分析 (初诊 ndGBM)
1. 加载 22 个初诊样本 (10x MTX)
2. 提取髓系细胞 (标志基因)
3. 注释 MO-TAM vs MG-TAM (正交标志基因, 与模块五一致)
4. UCell 计算 core68 与 hallmark IFN-a 得分
5. 患者层面 pseudobulk 配对比较
输出: GSE182109_myeloid_annot.csv / GSE182109_core68_validation.csv / GSE182109_validation_report.md
"""
import os, gzip, re, time
import numpy as np, pandas as pd
from scipy.io import mmread
from scipy import stats


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
EXT = os.path.join(WS, "temp", "GSE182109", "extracted")
OUT = os.path.join(WS, "output", "B6")
os.makedirs(OUT, exist_ok=True)

core = pd.read_csv(os.path.join(WS, "output", "B4", "B4b_ISG_core_set.csv"))
CORE68 = list(core["symbol"].dropna())
HALLMARK_IFNA = ("MX1 ISG15 OAS1 IFIT3 IFI44 IFI35 IRF7 RSAD2 IFI44L IFITM1 IFI27 IRF9 OASL "
"EIF2AK2 IFIT2 CXCL10 TAP1 SP110 DDX60 UBE2L6 USP18 PSMB8 IFIH1 BST2 LGALS3BP ADAR ISG20 "
"GBP2 IRF1 PLSCR1 PSMB9 HERC6 SAMD9 CMPK2 IFITM3 RTP4 STAT2 SAMD9L LY6E IFITM2 HELZ2 CXCL11 "
"TRIM21 PARP14 TRIM26 PARP12 NMI RNF31 HLA-C CASP1 TRIM14 TDRD7 DHX58 PARP9 PNPT1 TRIM25 PSME1 "
"WARS1 EPSTI1 UBA7 PSME2 B2M TRIM5 C1S LAP3 LAMP3 GBP4 NCOA7 TMEM140 CD74 GMPR PSMA3 PROCR IL7 "
"IFI30 IRF2 CSF1 IL15 CNP TENT5A IL4R CMTR1 CD47 LPAR6 MOV10 CASP8 TXNIP SLC25A28 SELL TRAFD1 "
"BATF2 RIPK2 CCRL2 NUB1 OGFR MVB12A ELF1").split()
MYELOID_MARKERS = ["P2RY12","TMEM119","CX3CR1","SLC1A3","CSF1R","C1QA","C1QB","TREM2","CD68",
                   "ITGAM","AIF1","LYZ","CD14","S100A8","S100A9","FCGR3A","FCN1","VCAN","ITGA4","TGFBI","FPR3"]
MO_MARKERS = ["S100A8","S100A9","VCAN","ITGA4","TGFBI","FPR3","CD14","LYZ"]
MG_MARKERS = ["P2RY12","TMEM119","CX3CR1","SLC1A3","CSF1R","TREM2","C1QA","C1QB"]

# 1. 加载初诊样本
nd_samples = sorted([f for f in os.listdir(EXT) if f.endswith("_matrix.mtx.gz") and "ndGBM" in f])
print("ndGBM samples:", len(nd_samples))
# 先读第一个样本的 features
feat0 = os.path.join(EXT, nd_samples[0].replace("_matrix.mtx.gz", "_features.tsv.gz"))
with gzip.open(feat0, "rt") as f:
    feats = [l.rstrip("\n").split("\t") for l in f]
gene_ids = [x[0] for x in feats]
gene_syms = [x[1] for x in feats]
print("genes:", len(gene_ids))

# 感兴趣的基因
interest = list(dict.fromkeys(CORE68 + HALLMARK_IFNA + MYELOID_MARKERS))
sym2row = {}
for i, s in enumerate(gene_syms):
    if s in interest and s not in sym2row:
        sym2row[s] = i
print("interest matched:", len(sym2row), "/", len(interest))
missing = [g for g in interest if g not in sym2row]
print("missing:", missing)

# 2. 逐样本加载, 提取基因子集
t0 = time.time()
cell_meta = []  # (sample, patient, barcode)
matrices = []
for fn in nd_samples:
    sample = fn.replace("_matrix.mtx.gz", "")
    m = sample.split("_")[1]  # ndGBM-01-A
    patient = re.sub(r"-[A-Z0-9]+$", "", m)  # ndGBM-01
    mtx_path = os.path.join(EXT, fn)
    bc_path = os.path.join(EXT, fn.replace("_matrix.mtx.gz", "_barcodes.tsv.gz"))
    with gzip.open(bc_path, "rt") as f:
        bcs = [l.rstrip("\n") for l in f]
    M = mmread(mtx_path).tocsr()
    rows = [sym2row[s] for s in gene_syms if s in sym2row]
    # 需要保持行顺序对应 sym2row 的基因
    sub = M[sorted(sym2row.values()), :]  # 按 sym2row 顺序
    matrices.append(sub)
    for b in bcs:
        cell_meta.append((sample, patient, b))
    print(f"  {sample}: {M.shape[1]} cells", flush=True)
print("loaded in", round(time.time()-t0,1), "s")

# 合并
gene_order = [g for g in interest if g in sym2row]
X = np.hstack([m.toarray() for m in matrices])  # n_genes x n_cells
X = np.log1p(X.astype(np.float32))
meta = pd.DataFrame(cell_meta, columns=["sample","patient","barcode"])
print("X shape:", X.shape, "cells:", meta.shape[0])
np.save(os.path.join(WS,"temp","GSE182109","ndGBM_interest_expr.npy"), X)
meta.to_csv(os.path.join(WS,"temp","GSE182109","ndGBM_meta.csv"), index=False)
with open(os.path.join(WS,"temp","GSE182109","ndGBM_gene_order.txt"),"w") as f:
    f.write("\n".join(gene_order))
print("saved intermediate")
