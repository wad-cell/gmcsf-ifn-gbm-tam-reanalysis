# -*- coding: utf-8 -*-
"""B7 打分函数库：等权平均 / UCell / AUCell / pseudobulk GSVA"""
import numpy as np, pandas as pd
from scipy.stats import rankdata

def mean_score(mat, genes):
    """等权平均 log1p 表达（B5/B6 一致）"""
    g = [x for x in genes if x in mat.columns]
    if not g:
        return pd.Series(np.nan, index=mat.index)
    return mat[g].mean(axis=1)

def ucell_score(mat, genes, max_rank=1500):
    """UCell (Andreatta & Carmona 2021) Python 实现
    对每个细胞：基因集内基因在细胞所有基因中的 rank（1=最高），
    取前 max_rank 个基因的平均 rank，归一化到 0-1。
    """
    g = [x for x in genes if x in mat.columns]
    if not g:
        return pd.Series(np.nan, index=mat.index)
    sub = mat[g].values  # n_cells x n_genes
    n_cells, n_genes = sub.shape
    # 对每个细胞，计算基因集内基因的 rank（在细胞所有基因中）
    # 用 rankdata 对整行（所有基因）排名，再取基因集内基因
    all_vals = mat.values  # n_cells x n_all
    n_all = all_vals.shape[1]
    # 为节省内存，分批处理
    scores = np.zeros(n_cells)
    batch = 500
    for i in range(0, n_cells, batch):
        idx = slice(i, min(i+batch, n_cells))
        # rankdata 每行：1=最小。我们想要 1=最高，用 n_all - rank + 1
        ranks = n_all - rankdata(all_vals[idx], axis=1) + 1  # 1=最高
        g_ranks = ranks[:, :n_genes]  # 基因集内基因（前 n_genes 列，需保证顺序）
        # 取前 max_rank 个（rank 值最小即最高）
        top = np.sort(g_ranks, axis=1)[:, :min(max_rank, n_genes)]
        scores[idx] = 1 - top.mean(axis=1) / max_rank
    return pd.Series(scores, index=mat.index)

def aucell_score(mat, genes):
    """AUCell (Aibar et al. 2017) Python 实现
    对每个细胞：基因集内基因按 rank 排序，计算 AUC（曲线下面积），归一化。
    """
    g = [x for x in genes if x in mat.columns]
    if not g:
        return pd.Series(np.nan, index=mat.index)
    sub = mat[g].values
    n_cells, n_genes = sub.shape
    all_vals = mat.values
    n_all = all_vals.shape[1]
    scores = np.zeros(n_cells)
    batch = 500
    for i in range(0, n_cells, batch):
        idx = slice(i, min(i+batch, n_cells))
        ranks = n_all - rankdata(all_vals[idx], axis=1) + 1  # 1=最高
        g_ranks = ranks[:, :n_genes]
        # 按 rank 排序（升序=最高优先）
        g_ranks_sorted = np.sort(g_ranks, axis=1)
        # AUC：对排序后的 rank 值，计算累积曲线下面积
        # 归一化：除以 (n_genes * n_all)
        auc = np.cumsum(g_ranks_sorted, axis=1).sum(axis=1) / (n_genes * n_all)
        scores[idx] = auc
    return pd.Series(scores, index=mat.index)

def pseudobulk_gsva(mat, genes, groupby):
    """pseudobulk GSVA：按 groupby 聚合为 pseudobulk，用 gseapy.gsva 打分"""
    import gseapy as gp
    g = [x for x in genes if x in mat.columns]
    if not g:
        return None
    # pseudobulk：按 groupby 求和
    pb = mat.groupby(groupby).sum()
    # 转置为 基因 x 样本
    pb_t = pb.T
    gs = {f"module_{len(g)}": g}
    try:
        res = gp.gsva(pb_t, gene_sets=gs, method="gsva", kcdf="Gaussian", verbose=False)
        # res: 基因集 x 样本
        return res.loc[f"module_{len(g)}"]
    except Exception as e:
        print("GSVA err:", e)
        return None
