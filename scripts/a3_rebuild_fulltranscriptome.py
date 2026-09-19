# -*- coding: utf-8 -*-
"""A3: GSE182109 full-transcriptome object rebuild.
Merge 44 GSM 10x matrices (official QC-filtered cells only) into one AnnData h5ad.
Gene space = union of ENSEMBL IDs across the two feature versions.
obs = official cell-level annotation (Meta_Data_GBMatlas.txt) + QC fields + inclusion flags.
"""
import os, sys, time, gzip, csv, json
import numpy as np
import scipy.io
import scipy.sparse as sp
import pandas as pd
import anndata as ad

t0=time.time()

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

ws = PROJ_ROOT
ext=os.path.join(ws,"temp","GSE182109","extracted")
meta_fp=os.path.join(ws,"temp","GBMatlas","Meta_Data_GBMatlas.txt")
outdir=os.path.join(ws,"output","C3.1")
os.makedirs(outdir,exist_ok=True)
h5ad_out=os.path.join(outdir,"GSE182109_fulltranscriptome.h5ad")

# ---- 1. load official annotation ----
print("[1] load official annotation...",flush=True)
cols_meta=["orig.ident","nCount_RNA","nFeature_RNA","mitoRatio","riboRatio",
           "CopyKatPrediction","Patient","Fragment","Type","Grade","GSMID",
           "barcode","Phase","sex","Cluster","Assignment","SubCluster","SubAssignment"]
official={}  # key -> row dict
qc_head={}
with open(meta_fp,encoding="utf-8",errors="replace") as f:
    f.readline()
    for line in f:
        p=line.rstrip("\n").split("\t")
        if len(p)<19: continue
        # 数据行 = [rownames] + 18 字段（对应表头 18 列）
        key=p[11]+"_"+p[12]
        official[key]=p[1:]
print("  official cells:",len(official),flush=True)

# ---- 2. build union gene space from all features files ----
print("[2] build union gene space...",flush=True)
gene_ids=[]; gene_sym={}
for f in sorted(os.listdir(ext)):
    if not f.endswith("_features.tsv.gz"): continue
    with gzip.open(os.path.join(ext,f),"rt",errors="replace") as fh:
        for line in fh:
            parts=line.rstrip("\n").split("\t")
            if len(parts)<2: continue
            gid=parts[0]; sym=parts[1]
            if gid not in gene_sym:
                gene_sym[gid]=sym
                gene_ids.append(gid)
gene_idx={g:i for i,g in enumerate(gene_ids)}
n_genes=len(gene_ids)
print("  union genes:",n_genes,flush=True)

# ---- 3. per-GSM: read mtx, reindex rows to union, subset cells to official keys ----
print("[3] merge matrices...",flush=True)
blocks=[]; obs_rows=[]
for gsm_idx, f in enumerate(sorted(os.listdir(ext))):
    if not f.endswith("_matrix.mtx.gz"): continue
    gsm=f.split("_")[0]
    bc_f=os.path.join(ext,gsm+"_"+f[len(gsm)+1:].replace("_matrix.mtx.gz","_barcodes.tsv.gz"))
    ft_f=os.path.join(ext,gsm+"_"+f[len(gsm)+1:].replace("_matrix.mtx.gz","_features.tsv.gz"))
    # barcodes
    bcs=[]
    with gzip.open(bc_f,"rt",errors="replace") as fh:
        for line in fh: bcs.append(line.strip())
    # features order
    feat_ids=[]
    with gzip.open(ft_f,"rt",errors="replace") as fh:
        for line in fh:
            parts=line.rstrip("\n").split("\t")
            feat_ids.append(parts[0] if parts else "")
    # matrix
    m=scipy.io.mmread(os.path.join(ext,f)).tocoo()
    assert m.shape[0]==len(feat_ids) and m.shape[1]==len(bcs), (gsm,m.shape,len(feat_ids),len(bcs))
    # map feature rows -> union idx
    row_idx=np.full(len(feat_ids),-1,dtype=np.int64)
    for i,gid in enumerate(feat_ids):
        row_idx[i]=gene_idx.get(gid,-1)
    # remap rows
    remap=row_idx[m.row]
    keep=remap>=0
    row=remap[keep]; col=m.col[keep]; data=m.data[keep]
    csr=sp.coo_matrix((data,(row,col)),shape=(n_genes,m.shape[1])).tocsr()
    # select official cells
    keep_cols=[]
    for j,b in enumerate(bcs):
        if gsm+"_"+b in official:
            keep_cols.append(j)
    if not keep_cols:
        print("  skip",gsm,"no official cells",flush=True); continue
    keep_cols=np.array(keep_cols,dtype=np.int64)
    sub=csr[:,keep_cols].tocsr()
    blocks.append(sub)
    for j in keep_cols:
        b=bcs[j]; key=gsm+"_"+b
        p=official[key]
        obs_rows.append(p)
    if (gsm_idx+1)%8==0 or gsm_idx==43:
        print(f"  [{gsm_idx+1}/44] {gsm}: cells={len(keep_cols)} nnz={sub.nnz} elapsed={time.time()-t0:.0f}s",flush=True)

print("[4] concat...",flush=True)
X=sp.hstack(blocks,format="csr")
del blocks
X=X.astype(np.float32)
print("  final X shape:",X.shape,"nnz:",X.nnz,flush=True)

# ---- 5. build obs ----
print("[5] build obs DataFrame...",flush=True)
cols_full=cols_meta
df_obs=pd.DataFrame(obs_rows,columns=cols_full)
# index = gsm_barcode
idx=df_obs["GSMID"]+"_"+df_obs["barcode"]
df_obs.index=idx.values
# QC 数值列
for c in ["nCount_RNA","nFeature_RNA","mitoRatio","riboRatio"]:
    df_obs[c]=pd.to_numeric(df_obs[c],errors="coerce")
# disease stage / patient / region
df_obs.rename(columns={"Patient":"patient","Fragment":"region","Type":"disease_stage",
                       "Assignment":"cell_type","SubCluster":"subcluster","SubAssignment":"subassignment",
                       "GSMID":"gsm","barcode":"barcode"},inplace=True)
# MO/MG mapping (official labels)
sa_list=df_obs["subassignment"].tolist()
map_list=[]
for sa in sa_list:
    if sa in ("i-microglia","h-microglia","AP-microglia","a-microglia"): map_list.append("MG-TAM")
    elif sa in ("s-mac 1","s-mac 2"): map_list.append("MO-TAM")
    elif sa in ("MDSC","DCs"): map_list.append("other-myeloid")
    elif sa=="NA": map_list.append("NA")
    else: map_list.append("NA")
df_obs["mo_mg"]=map_list
# inclusion flags
def flag(row):
    # primary: Myeloid & MO/MG definite; exclude uncertain(NA), low-quality, doublet-suspect
    if row["cell_type"]!="Myeloid": return "exclude_nonmyeloid"
    if row["mo_mg"] in ("MG-TAM","MO-TAM"): return "include_primary"
    return "exclude_unresolved"
df_obs["primary_flag"]=df_obs.apply(flag,axis=1)
# QC flag (low quality: nFeature<200 or mitoRatio>0.2)
df_obs["qc_flag"]=np.where((df_obs["nFeature_RNA"]<200)|(df_obs["mitoRatio"]>0.2),"low_quality","ok")
print("  obs shape:",df_obs.shape,flush=True)

# ---- 6. var ----
var_df=pd.DataFrame({"ensembl_id":gene_ids,"gene_symbol":[gene_sym[g] for g in gene_ids]})
var_df.index=gene_ids

# ---- 7. save h5ad ----
print("[6] save h5ad...",flush=True)
X=X.T.tocsr().astype(np.float32)  # 转置为 cells × genes
adata=ad.AnnData(X=X,obs=df_obs,var=var_df)
adata.write_h5ad(h5ad_out)
print("  saved:",h5ad_out,flush=True)
print("DONE total elapsed:",time.time()-t0,"s",flush=True)
