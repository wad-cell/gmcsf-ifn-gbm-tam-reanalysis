# -*- coding: utf-8 -*-
"""B1 基因注释：Ensembl ID -> 基因符号（mygene），更新结果表"""
import os, time, pandas as pd, numpy as np, mygene


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

OUT = os.path.join(PROJ_ROOT, "output", "B1")
ANNOT = os.path.join(OUT, "gene_annotation.csv")

# 收集所有基因 ID
all_ids = set()
for block, cid in [("data1", "C1_culture_main"), ("data2", "C8_GMCSF_vs_ctrl")]:
    f = os.path.join(OUT, f"{block}_{cid}.csv")
    df = pd.read_csv(f, usecols=["gene"])
    all_ids.update(df["gene"].str.split(".").str[0])
all_ids = sorted(all_ids)
print(f"总基因数: {len(all_ids)}")

# 检查是否已有注释
if os.path.exists(ANNOT):
    annot = pd.read_csv(ANNOT)
    print(f"已有注释 {len(annot)} 行, 跳过查询")
else:
    mg = mygene.MyGeneInfo()
    results = []
    CHUNK = 1000
    for i in range(0, len(all_ids), CHUNK):
        chunk = all_ids[i:i+CHUNK]
        r = mg.querymany(chunk, scopes="ensembl.gene", fields="symbol", species="human", verbose=False)
        for item in r:
            if "query" in item and "symbol" in item:
                results.append({"ensembl_id": item["query"], "symbol": item["symbol"]})
        print(f"  {min(i+CHUNK, len(all_ids))}/{len(all_ids)}")
        time.sleep(0.3)
    annot = pd.DataFrame(results).drop_duplicates("ensembl_id")
    annot.to_csv(ANNOT, index=False)
    print(f"注释完成: {len(annot)} 个基因有符号")

# 更新所有结果表，添加 symbol 列
for block in ["data1", "data2"]:
    for cid in ["C1_culture_main","C2_antibody_main","C3_interaction","C4_co_vs_mono_IgG",
                "C5_co_vs_mono_aGMCSF","C6_aGMCSF_vs_IgG_mono","C7_aGMCSF_vs_IgG_co",
                "C8_GMCSF_vs_ctrl","C9_LN229_vs_mono","C10_U251_vs_mono","C11_LN229_vs_U251",
                "C12_pooled_co_vs_mono"]:
        f = os.path.join(OUT, f"{block}_{cid}.csv")
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f)
        if "symbol" in df.columns:
            continue
        df["ensembl_id"] = df["gene"].str.split(".").str[0]
        df = df.merge(annot, on="ensembl_id", how="left")
        # 列顺序: gene, symbol, baseMean...
        cols = ["gene", "symbol"] + [c for c in df.columns if c not in ("gene", "symbol", "ensembl_id")]
        df = df[cols]
        df.to_csv(f, index=False)
        print(f"更新 {block}_{cid}: {len(df)} 行")

print("===== 注释完成 =====")
