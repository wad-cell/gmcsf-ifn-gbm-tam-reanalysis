# -*- coding: utf-8 -*-
"""B5 模块五：富集与调控分析
对 68 核心基因:
1. GO BP / Reactome / KEGG / Hallmark ORA (自定义超几何, 背景=全部可检测基因)
2. Hallmark + Reactome 预排序 GSEA (B2_CoCul_vs_MoCul 全基因 log2FC)
3. 上游 TF (ChEA_2022 + TRRUST), 重点 STAT1/2, IRF1/7/9
4. STRING PPI (API)
输出: enrichment_ORA_full.csv / enrichment_GSEA_full.csv / upstream_TF_full.csv / string_ppi.csv
"""
import pandas as pd, numpy as np, os, json, urllib.request
from scipy import stats
import gseapy as gp

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
OUT = os.path.join(WS, "output", "B5")
B2 = os.path.join(WS, "output", "B2")

core = pd.read_csv(os.path.join(OUT, "core68_evidence_matrix.csv"))
core_syms = set(core["symbol"].dropna())
print("68 core genes:", len(core_syms))

# 背景: GSE309037 全部可检测基因 (baseMean>0)
deseq = pd.read_csv(os.path.join(B2, "GSE309037_B2_CoCul_vs_MoCul.csv"))
bg = deseq[deseq["baseMean"] > 0]["symbol"].dropna().unique()
print("background genes (baseMean>0):", len(bg))
bg_set = set(bg)

# ---------- 1. ORA: 自定义超几何 ----------
def hypergeom_ora(genes, gene_set, N):
    """超几何检验: 68 基因在通路中的富集"""
    k = len(genes & gene_set)          # 命中数
    K = len(genes)                     # 查询基因数
    M = len(gene_set & bg_set)         # 通路中可检测基因数
    if k == 0 or M == 0: return None
    # 超几何: 从 N 中抽 K, 通路 M 个, 命中 >= k 的概率
    p = stats.hypergeom.sf(k-1, N, M, K)
    return k, K, M, N, p

libs = {}
for name in ["GO_Biological_Process_2023","Reactome_2022","KEGG_2021_Human","MSigDB_Hallmark_2020"]:
    libs[name] = gp.get_library(name=name)
    print(name, "gene sets:", len(libs[name]))

N = len(bg_set)
ora_rows = []
for libname, gs in libs.items():
    for term, genes in gs.items():
        gset = set(genes)
        r = hypergeom_ora(core_syms, gset, N)
        if r is None: continue
        k, K, M, Nn, p = r
        ora_rows.append(dict(library=libname, term=term, hit_genes=k, query_size=K,
                             pathway_size=M, background_size=Nn,
                             gene_ratio=k/K, background_ratio=M/Nn, pvalue=p))
ora = pd.DataFrame(ora_rows)
# FDR (BH) 按库内
for libname in ora["library"].unique():
    sel = ora["library"]==libname
    ora.loc[sel, "FDR"] = stats.false_discovery_control(ora.loc[sel,"pvalue"].values, method="bh")
ora = ora.sort_values(["library","pvalue"])
ora.to_csv(os.path.join(OUT, "enrichment_ORA_full.csv"), index=False)
print("\n=== ORA 显著 (FDR<0.05) 各库 top ===")
for libname in ora["library"].unique():
    sub = ora[(ora["library"]==libname) & (ora["FDR"]<0.05)]
    print(f"\n[{libname}] {len(sub)} significant")
    print(sub.head(8)[["term","hit_genes","gene_ratio","background_ratio","pvalue","FDR"]].to_string(index=False))

# ---------- 2. 预排序 GSEA (Hallmark + Reactome) ----------
rnk = deseq[["symbol","log2FoldChange"]].dropna().drop_duplicates(subset="symbol")
rnk = rnk.set_index("symbol")["log2FoldChange"]
print("\nGSEA ranking genes:", len(rnk))
gsea_res = gp.prerank(rnk=rnk, gene_sets=libs["MSigDB_Hallmark_2020"], min_size=5, max_size=500,
                      permutation_num=1000, outdir=None, seed=42, threads=4)
gsea_hall = gsea_res.res2d
gsea_hall["library"] = "MSigDB_Hallmark_2020"
gsea_res2 = gp.prerank(rnk=rnk, gene_sets=libs["Reactome_2022"], min_size=5, max_size=500,
                       permutation_num=1000, outdir=None, seed=42, threads=4)
gsea_react = gsea_res2.res2d
gsea_react["library"] = "Reactome_2022"
gsea_all = pd.concat([gsea_hall, gsea_react], ignore_index=True)
gsea_all.to_csv(os.path.join(OUT, "enrichment_GSEA_full.csv"), index=False)
print("\n=== GSEA 显著 (FDR<0.05) ===")
sig = gsea_all[gsea_all["FDR q-val"]<0.05].sort_values("FDR q-val")
print(sig[["library","Term","NES","NOM p-val","FDR q-val"]].to_string(index=False))

# ---------- 3. 上游 TF (ChEA + TRRUST) ----------
tf_rows = []
for libname in ["ChEA_2022","TRRUST_Transcription_Factors_2019"]:
    gs = gp.get_library(name=libname)
    for term, genes in gs.items():
        gset = set(genes)
        r = hypergeom_ora(core_syms, gset, N)
        if r is None: continue
        k, K, M, Nn, p = r
        tf_rows.append(dict(library=libname, tf=term, hit_genes=k, query_size=K,
                            pathway_size=M, background_size=Nn, gene_ratio=k/K, pvalue=p))
tf = pd.DataFrame(tf_rows)
for libname in tf["library"].unique():
    sel = tf["library"]==libname
    tf.loc[sel, "FDR"] = stats.false_discovery_control(tf.loc[sel,"pvalue"].values, method="bh")
tf = tf.sort_values("pvalue")
tf.to_csv(os.path.join(OUT, "upstream_TF_full.csv"), index=False)
print("\n=== 上游 TF top20 ===")
print(tf.head(20)[["library","tf","hit_genes","gene_ratio","pvalue","FDR"]].to_string(index=False))
# 重点 TF
focus = ["STAT1","STAT2","IRF1","IRF7","IRF9"]
print("\n=== 重点 TF 检查 ===")
for f in focus:
    hit = tf[tf["tf"].str.contains(f, case=False, na=False)]
    if len(hit):
        print(hit[["library","tf","hit_genes","pvalue","FDR"]].to_string(index=False))
    else:
        print(f, "not found in TF libraries")

# ---------- 4. STRING PPI ----------
print("\n=== STRING PPI ===")
try:
    genes_str = "%0d".join(sorted(core_syms))
    url = ("https://string-db.org/api/tsv/network?identifiers=" + genes_str +
           "&species=9606&required_score=400&caller_identity=gmcsf_ifn_b5")
    req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        txt = resp.read().decode()
    if txt.strip().startswith("Error") or "preferredName" not in txt:
        print("STRING returned no usable network:", txt[:200])
    else:
        lines = [l.split("\t") for l in txt.strip().splitlines()]
        header = lines[0]
        data = [dict(zip(header, r)) for r in lines[1:]]
        ppi = pd.DataFrame(data)
        ppi.to_csv(os.path.join(OUT, "string_ppi.csv"), index=False)
        print("STRING edges:", len(ppi))
        # 节点度
        nodes = pd.concat([ppi["preferredName_A"], ppi["preferredName_B"]]).value_counts()
        nodes.to_csv(os.path.join(OUT, "string_node_degree.csv"))
        print("top degree nodes:", nodes.head(10).to_dict())
except Exception as e:
    print("STRING API error:", e)

print("\nDONE 模块五")
