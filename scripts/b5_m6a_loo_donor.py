# -*- coding: utf-8 -*-
"""B5 模块六 Part A: leave-one-donor-out (可续跑, 参数化)
用法: python b5_m6a_loo_donor.py <block> <contrast>
  block: data1 | data2 | GSE309037
  contrast: C4 | C7 | C9 | C8 | B2
每次拟合后立即落盘, 已存在则跳过。
输出: robustness_leave_one_donor_out.csv
"""
import os, sys, time, json
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
OUT = os.path.join(WS, "output", "B5")
RAW = os.path.join(WS, "temp")
N_CPUS = 4
OUT_CSV = os.path.join(OUT, "robustness_leave_one_donor_out.csv")

annot = pd.read_csv(os.path.join(WS, "output", "B1", "gene_annotation.csv"))
core = pd.read_csv(os.path.join(OUT, "core68_evidence_matrix.csv"))
core_syms = set(core["symbol"].dropna())

def load_counts(path):
    df = pd.read_csv(path, sep=";", index_col=0)
    df = df[~df.index.duplicated(keep="first")].astype(int)
    return df

def build_meta(cols):
    meta = pd.DataFrame(index=cols)
    meta["condition"] = [c.rsplit("_n", 1)[0] for c in cols]
    meta["donor"] = ["donor_" + c.split("_n")[1] for c in cols]
    return meta

def ensg_to_symbol(res):
    res = res.reset_index().rename(columns={"index":"gene"})
    res["gene"] = res["gene"].astype(str).str.split(".").str[0]
    m = annot.set_index("ensembl_id")["symbol"].to_dict()
    res["symbol"] = res["gene"].map(m)
    return res

def fit_and_contrast(counts, meta, contrast_vec):
    dds = DeseqDataSet(counts=counts.T, metadata=meta,
                       design_factors=["donor","condition"],
                       refit_cooks=True, n_cpus=N_CPUS)
    dds.deseq2()
    dm_cols = list(dds.obsm["design_matrix"].columns)
    colmap = {c:i for i,c in enumerate(dm_cols)}
    v = np.zeros(len(dm_cols))
    for level, coef in contrast_vec.items():
        col = f"condition[T.{level}]"
        if col in colmap: v[colmap[col]] = coef
        else: print(f"  warn: {col} not in design (reference?)")
    st = DeseqStats(dds, contrast=v, n_cpus=N_CPUS, quiet=True)
    st.summary()
    return st.results_df

def summarize(res, label):
    res = ensg_to_symbol(res)
    sub = res[res["symbol"].isin(core_syms)].dropna(subset=["log2FoldChange"])
    n = len(sub)
    mean_lfc = sub["log2FoldChange"].mean()
    frac_neg = (sub["log2FoldChange"] < 0).mean()
    frac_sig = (sub["padj"] < 0.05).mean()
    return dict(label=label, n_genes=n, mean_log2FC=mean_lfc,
                frac_negative=frac_neg, frac_sig_padj05=frac_sig)

def append_row(row):
    df = pd.DataFrame([row])
    if os.path.exists(OUT_CSV):
        old = pd.read_csv(OUT_CSV)
        if row["label"] in set(old["label"]):
            print(f"  {row['label']} 已存在, 跳过")
            return
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"  saved {row['label']}")

block, contrast = sys.argv[1], sys.argv[2]

if block == "data1":
    counts = load_counts(os.path.join(RAW, "GSE309039", "raw", "GSE309039_data1_raw_counts.csv.gz"))
    meta = build_meta(counts.columns)
    donors = ["donor_1","donor_2","donor_3","donor_4","donor_5","donor_6"]
    cv = {"C4":{"MoCul_IgG":-1.0}, "C7":{"CoCul_aGMCSF":1.0}}[contrast]
elif block == "data2":
    counts = load_counts(os.path.join(RAW, "GSE309039", "raw", "GSE309039_data2_raw_counts.csv.gz"))
    meta = build_meta(counts.columns)
    donors = ["donor_1","donor_2","donor_3","donor_4","donor_5","donor_6"]
    cv = {"C9":{"MoCul":-1.0}, "C8":{"MoCul_GMCSF":1.0,"MoCul":-1.0}}[contrast]
elif block == "GSE309037":
    counts = load_counts(os.path.join(RAW, "GSE309037", "raw", "GSE309037_raw_counts.csv.gz"))
    meta = build_meta(counts.columns)
    donors = ["donor_1","donor_2","donor_3"]
    cv = {"B2":{"MoCul":-1.0}}[contrast]
else:
    raise SystemExit("bad block")

for donor in donors:
    keep = meta["donor"] != donor
    c = counts.loc[:, keep.values]; m = meta[keep].copy()
    t0 = time.time()
    r = fit_and_contrast(c, m, cv)
    row = summarize(r, f"{block}_loo_{donor}_{contrast}")
    append_row(row)
    print(f"  {block} loo {donor} {contrast} done {time.time()-t0:.0f}s")

print(f"DONE {block} {contrast}")
