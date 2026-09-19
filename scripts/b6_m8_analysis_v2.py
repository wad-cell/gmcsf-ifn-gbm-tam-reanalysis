# -*- coding: utf-8 -*-
"""B6 模块八: GSE182109 独立队列主分析 (v2)
修复: 1) 患者提取正则 2) 平均表达打分替代简化UCell 3) 严格标志基因注释
"""
import os, re
import numpy as np, pandas as pd
from scipy import stats


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
TMP = os.path.join(WS, "temp", "GSE182109")
OUT = os.path.join(WS, "output", "B6")
os.makedirs(OUT, exist_ok=True)

X = np.load(os.path.join(TMP, "ndGBM_interest_expr.npy"))
meta = pd.read_csv(os.path.join(TMP, "ndGBM_meta.csv"))
with open(os.path.join(TMP, "ndGBM_gene_order.txt")) as f:
    gene_order = [l.strip() for l in f if l.strip()]
g2i = {g: i for i, g in enumerate(gene_order)}

# 修复患者提取: GSM5518600_ndGBM-01-A -> ndGBM-01 -> 01
def extract_patient(sample):
    m = re.match(r"ndGBM-(\d+)", sample.split("_")[1])
    return m.group(1) if m else sample
meta["patient"] = meta["sample"].map(extract_patient)
print("patients:", sorted(meta["patient"].unique()))

# 严格标志基因 (排除泛髓系 CD14/LYZ/C1QA/C1QB)
MO_MARKERS = ["S100A8","S100A9","VCAN","ITGA4","TGFBI","FPR3"]
MG_MARKERS = ["P2RY12","TMEM119","CX3CR1","SLC1A3","TREM2"]
MYELOID = list(dict.fromkeys(MO_MARKERS + MG_MARKERS + ["CD68","ITGAM","AIF1","CSF1R","C1QA","C1QB","CD14","LYZ","FCGR3A","FCN1"]))
mo_idx = [g2i[g] for g in MO_MARKERS if g in g2i]
mg_idx = [g2i[g] for g in MG_MARKERS if g in g2i]
my_idx = [g2i[g] for g in MYELOID if g in g2i]

# 1. 髓系筛选
my_expr = X[my_idx, :]
n_my = (my_expr > 0).sum(axis=0)
is_myeloid = n_my >= 2
print("myeloid cells:", is_myeloid.sum(), "/", len(is_myeloid))

# 2. 注释: 特异性标志基因平均表达差
mo_score = X[mo_idx, :].mean(axis=0)
mg_score = X[mg_idx, :].mean(axis=0)
diff = mo_score - mg_score
my_cells = np.where(is_myeloid)[0]
lineage = np.where(diff[my_cells] > 0, "MO-TAM", "MG-TAM")
# 排除表达量过低的不确定细胞
my_mean = my_expr.mean(axis=0)[my_cells]
ambig = my_mean < 0.05
lineage[ambig] = "Ambiguous"
print("MO-TAM:", (lineage=="MO-TAM").sum(), "MG-TAM:", (lineage=="MG-TAM").sum(), "Ambiguous:", ambig.sum())

# 3. 程序得分: 平均 log1p 表达 (与 B5 模块评分一致)
core = pd.read_csv(os.path.join(WS,"output","B4","B4b_ISG_core_set.csv"))
core_idx = [g2i[g] for g in core["symbol"].dropna() if g in g2i]
HALLMARK_IFNA = ("MX1 ISG15 OAS1 IFIT3 IFI44 IFI35 IRF7 RSAD2 IFI44L IFITM1 IFI27 IRF9 OASL "
"EIF2AK2 IFIT2 CXCL10 TAP1 SP110 DDX60 UBE2L6 USP18 PSMB8 IFIH1 BST2 LGALS3BP ADAR ISG20 "
"GBP2 IRF1 PLSCR1 PSMB9 HERC6 SAMD9 CMPK2 IFITM3 RTP4 STAT2 SAMD9L LY6E IFITM2 HELZ2 CXCL11 "
"TRIM21 PARP14 TRIM26 PARP12 NMI RNF31 HLA-C CASP1 TRIM14 TDRD7 DHX58 PARP9 PNPT1 TRIM25 PSME1 "
"WARS1 EPSTI1 UBA7 PSME2 B2M TRIM5 C1S LAP3 LAMP3 GBP4 NCOA7 TMEM140 CD74 GMPR PSMA3 PROCR IL7 "
"IFI30 IRF2 CSF1 IL15 CNP TENT5A IL4R CMTR1 CD47 LPAR6 MOV10 CASP8 TXNIP SLC25A28 SELL TRAFD1 "
"BATF2 RIPK2 CCRL2 NUB1 OGFR MVB12A ELF1").split()
hall_idx = [g2i[g] for g in HALLMARK_IFNA if g in g2i]
print("core68 matched:", len(core_idx), "hallmark matched:", len(hall_idx))
core_score = X[core_idx, :].mean(axis=0)
hall_score = X[hall_idx, :].mean(axis=0)

# 4. 髓系细胞表
my_meta = meta.iloc[my_cells].copy()
my_meta["lineage"] = lineage
my_meta["core68_score"] = core_score[my_cells]
my_meta["hallmark_score"] = hall_score[my_cells]
my_meta["mo_marker_mean"] = mo_score[my_cells]
my_meta["mg_marker_mean"] = mg_score[my_cells]
my_meta.to_csv(os.path.join(OUT, "GSE182109_myeloid_annot.csv"), index=False)
print("saved GSE182109_myeloid_annot.csv")

# 5. 患者层面
mo = my_meta[my_meta["lineage"]=="MO-TAM"]
mg = my_meta[my_meta["lineage"]=="MG-TAM"]
ct = pd.crosstab(my_meta["patient"], my_meta["lineage"])
print("\n=== 患者 x lineage 细胞数 ===")
print(ct.to_string())

mo_pb = mo.groupby("patient")["core68_score"].mean()
mg_pb = mg.groupby("patient")["core68_score"].mean()
common = sorted(set(mo_pb.index) & set(mg_pb.index))
print("\npatients with both:", len(common))
d = mo_pb[common] - mg_pb[common]
print("core68: MO-MG per patient:")
for p in common:
    print(f"  {p}: MO={mo_pb[p]:.4f} MG={mg_pb[p]:.4f} diff={d[p]:+.4f}")
if len(common) >= 3:
    w = stats.wilcoxon(mo_pb[common], mg_pb[common])
    t = stats.ttest_rel(mo_pb[common], mg_pb[common])
    print(f"core68: Wilcoxon p={w.pvalue:.4f}, paired t p={t.pvalue:.4f}, mean diff={d.mean():+.4f}, n_support={(d<0).sum()}/{len(common)}")

mo_h = mo.groupby("patient")["hallmark_score"].mean()
mg_h = mg.groupby("patient")["hallmark_score"].mean()
dh = mo_h[common] - mg_h[common]
if len(common) >= 3:
    wh = stats.wilcoxon(mo_h[common], mg_h[common])
    th = stats.ttest_rel(mo_h[common], mg_h[common])
    print(f"hallmark: Wilcoxon p={wh.pvalue:.4f}, paired t p={th.pvalue:.4f}, mean diff={dh.mean():+.4f}, n_support={(dh<0).sum()}/{len(common)}")

res = pd.DataFrame({
    "patient": common,
    "MO_core68": [mo_pb[p] for p in common],
    "MG_core68": [mg_pb[p] for p in common],
    "diff_core68": [d[p] for p in common],
    "MO_hallmark": [mo_h[p] for p in common],
    "MG_hallmark": [mg_h[p] for p in common],
    "diff_hallmark": [dh[p] for p in common],
})
res.to_csv(os.path.join(OUT, "GSE182109_core68_validation.csv"), index=False)
print("\nsaved GSE182109_core68_validation.csv")
