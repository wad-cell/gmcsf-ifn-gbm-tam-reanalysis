# -*- coding: utf-8 -*-
"""获取背景基因的基因长度 (Ensembl REST 并发批量)"""
import os, json, time, urllib.request
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed


# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
B2 = os.path.join(WS, "output", "B2")
TEMP = os.path.join(WS, "temp")
CACHE = os.path.join(TEMP, "gene_lengths_cache.json")

deseq = pd.read_csv(os.path.join(B2, "GSE309037_B2_CoCul_vs_MoCul.csv"))
bg = deseq.dropna(subset=["baseMean", "log2FoldChange"])
bg_ens = bg["gene"].str.split(".").str[0].tolist()
print("bg genes:", len(bg_ens))

# 已有缓存
cache = {}
if os.path.exists(CACHE):
    with open(CACHE) as f:
        cache = json.load(f)
todo = [g for g in bg_ens if g not in cache]
print("cached:", len(cache), "todo:", len(todo))

def fetch_batch(ids):
    body = json.dumps({"ids": ids}).encode()
    req = urllib.request.Request(
        "https://rest.ensembl.org/lookup/id?content-type=application/json",
        data=body, headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                res = json.loads(r.read().decode())
            out = {}
            for gid, info in res.items():
                if info and "start" in info and "end" in info:
                    out[gid] = info["end"] - info["start"] + 1
            return out
        except Exception as e:
            time.sleep(2 * (attempt + 1))
    return {}

BATCH = 1000
batches = [todo[i:i+BATCH] for i in range(0, len(todo), BATCH)]
print("batches:", len(batches))
t0 = time.time()
with ThreadPoolExecutor(max_workers=8) as ex:
    futs = [ex.submit(fetch_batch, b) for b in batches]
    done = 0
    for f in as_completed(futs):
        cache.update(f.result())
        done += 1
        if done % 5 == 0:
            print(f"  {done}/{len(batches)} batches, {len(cache)} genes, {round(time.time()-t0,1)}s")
with open(CACHE, "w") as f:
    json.dump(cache, f)
print("DONE. total cached:", len(cache), "in", round(time.time()-t0,1), "s")
