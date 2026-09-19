# -*- coding: utf-8 -*-
"""B5 模块六 Part B: 非随机性检验
1. 超几何检验: 68 = IFN ISG(760) ∩ CoCul 下调(338)
2. 阈值敏感性: 不同 padj/log2FC 阈值下交集大小
3. 10,000 次匹配置换: 匹配 baseMean 表达水平的随机 68 基因集,
   检验方向一致率/平均效应/相关系数是否超过随机预期
输出: robustness_analysis.csv / permutation_distributions.csv
"""
import os
import numpy as np, pandas as pd
from scipy import stats
WS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(WS, "output", "B5")
B2 = os.path.join(WS, "output", "B2")
B1 = os.path.join(WS, "output", "B1")

core = pd.read_csv(os.path.join(OUT, "core68_evidence_matrix.csv"))
core_syms = set(core["symbol"].dropna())
print("core:", len(core_syms))

# ---------- 1. 超几何检验 ----------
ifn = pd.read_csv(os.path.join(B2, "IFN_gold_standard_ISG_up.csv"))
cocul_sig = pd.read_csv(os.path.join(B2, "CoCul_vs_MoCul_sig.csv"))
deseq = pd.read_csv(os.path.join(B2, "GSE309037_B2_CoCul_vs_MoCul.csv"))
bg = deseq[deseq["baseMean"] > 0]["symbol"].dropna().unique()
N = len(bg)
A = set(ifn["symbol"].dropna())
B = set(cocul_sig[cocul_sig["log2FoldChange"] < 0]["symbol"].dropna())  # 338 下调
inter = A & B
print(f"背景 N={N}, IFN ISG={len(A)}, CoCul down={len(B)}, 交集={len(inter)}")
# 超几何: 从 N 抽 |A|, 其中 |B| 个标记, 命中 >= |inter|
p_hyper = stats.hypergeom.sf(len(inter)-1, N, len(B), len(A))
print(f"超几何 p = {p_hyper:.3e}")
# 期望交集
exp = len(A)*len(B)/N
print(f"期望交集 = {exp:.1f}, 富集倍数 = {len(inter)/exp:.1f}")

# ---------- 2. 阈值敏感性 ----------
def overlap_size(padj_cut, lfc_cut):
    a = set(ifn[(ifn["padj"]<padj_cut) & (ifn["log2FoldChange"]>lfc_cut)]["symbol"])
    b = set(cocul_sig[(cocul_sig["padj"]<padj_cut) & (cocul_sig["log2FoldChange"]<-lfc_cut)]["symbol"])
    return len(a), len(b), len(a&b)

sens_rows = []
for padj in [0.01, 0.05, 0.1]:
    for lfc in [0.5, 1.0, 1.5]:
        na, nb, ni = overlap_size(padj, lfc)
        sens_rows.append(dict(padj_cutoff=padj, log2FC_cutoff=lfc,
                              IFN_ISG_n=na, CoCul_down_n=nb, overlap_n=ni))
sens = pd.DataFrame(sens_rows)
print("\n=== 阈值敏感性 ===")
print(sens.to_string(index=False))

# ---------- 3. 10,000 次匹配置换 ----------
# 匹配 baseMean: 将背景基因按 baseMean 分箱, 68 核心基因所在箱内随机抽取
deseq2 = deseq.dropna(subset=["baseMean","log2FoldChange"]).copy()
deseq2["symbol"] = deseq2["symbol"].astype(str)
core_lfc = deseq2[deseq2["symbol"].isin(core_syms)].set_index("symbol")["log2FoldChange"]
print("\ncore genes in deseq:", len(core_lfc))

# baseMean 分箱 (log10)
deseq2["bm_bin"] = np.floor(np.log10(deseq2["baseMean"]+1)*2)/2
core_bins = deseq2[deseq2["symbol"].isin(core_syms)]["bm_bin"].value_counts()
print("core bin distribution:", dict(core_bins))

# 实际统计量
actual_frac_neg = (core_lfc < 0).mean()
actual_mean_lfc = core_lfc.mean()
# 相关系数: GM-CSF(C8) vs LN229(C9) 在 68 基因上的 rho
c8 = pd.read_csv(os.path.join(B1, "data2_C8_GMCSF_vs_ctrl.csv"))
c9 = pd.read_csv(os.path.join(B1, "data2_C9_LN229_vs_mono.csv"))
c8m = c8.drop_duplicates("symbol").set_index("symbol")["log2FoldChange"]
c9m = c9.drop_duplicates("symbol").set_index("symbol")["log2FoldChange"]
common = core_syms & set(c8m.index) & set(c9m.index)
actual_rho = stats.spearmanr(c8m[list(common)], c9m[list(common)]).statistic
print(f"actual: frac_neg={actual_frac_neg:.3f}, mean_lfc={actual_mean_lfc:.3f}, rho(C8,C9)={actual_rho:.3f}")

# 置换
rng = np.random.default_rng(42)
n_perm = 10000
perm_frac, perm_mean, perm_rho = [], [], []
# 预构建每个 bin 的候选基因
bin_pool = {b: deseq2[deseq2["bm_bin"]==b]["symbol"].values for b in core_bins.index}
for i in range(n_perm):
    chosen = []
    for b, cnt in core_bins.items():
        pool = bin_pool[b]
        chosen.extend(rng.choice(pool, size=cnt, replace=False))
    lfc = deseq2.set_index("symbol").loc[chosen, "log2FoldChange"]
    perm_frac.append((lfc < 0).mean())
    perm_mean.append(lfc.mean())
    c8s = c8m.reindex(chosen).dropna(); c9s = c9m.reindex(chosen).dropna()
    common_p = list(set(c8s.index) & set(c9s.index))
    if len(common_p) >= 5:
        perm_rho.append(stats.spearmanr(c8s[list(common_p)], c9s[list(common_p)]).statistic)
    else:
        perm_rho.append(np.nan)

perm_frac = np.array(perm_frac); perm_mean = np.array(perm_mean); perm_rho = np.array(perm_rho)
emp_frac = (perm_frac >= actual_frac_neg).mean()
emp_mean = (perm_mean <= actual_mean_lfc).mean()
emp_rho = (np.abs(perm_rho) >= abs(actual_rho)).mean()
print(f"\n=== 10,000 置换结果 ===")
print(f"方向一致率: 实际={actual_frac_neg:.3f}, 随机均值={perm_frac.mean():.3f} (95%: {np.percentile(perm_frac,2.5):.3f}-{np.percentile(perm_frac,97.5):.3f}), 经验P={emp_frac:.4f}")
print(f"平均效应: 实际={actual_mean_lfc:.3f}, 随机均值={perm_mean.mean():.3f} (95%: {np.percentile(perm_mean,2.5):.3f}-{np.percentile(perm_mean,97.5):.3f}), 经验P={emp_mean:.4f}")
print(f"rho(C8,C9): 实际={actual_rho:.3f}, 随机均值={np.nanmean(perm_rho):.3f} (95%: {np.nanpercentile(perm_rho,2.5):.3f}-{np.nanpercentile(perm_rho,97.5):.3f}), 经验P={emp_rho:.4f}")

# 保存
rob = pd.DataFrame([
    dict(test="hypergeometric_68_overlap", N=bg, IFN_ISG=len(A), CoCul_down=len(B),
         overlap=len(inter), expected=exp, enrichment=len(inter)/exp, pvalue=p_hyper),
    dict(test="permutation_frac_negative", observed=actual_frac_neg,
         perm_mean=perm_frac.mean(), perm_lo=np.percentile(perm_frac,2.5), perm_hi=np.percentile(perm_frac,97.5),
         empirical_p=emp_frac),
    dict(test="permutation_mean_log2FC", observed=actual_mean_lfc,
         perm_mean=perm_mean.mean(), perm_lo=np.percentile(perm_mean,2.5), perm_hi=np.percentile(perm_mean,97.5),
         empirical_p=emp_mean),
    dict(test="permutation_rho_GMCSF_vs_LN229", observed=actual_rho,
         perm_mean=np.nanmean(perm_rho), perm_lo=np.nanpercentile(perm_rho,2.5), perm_hi=np.nanpercentile(perm_rho,97.5),
         empirical_p=emp_rho),
])
rob.to_csv(os.path.join(OUT, "robustness_analysis.csv"), index=False)
sens.to_csv(os.path.join(OUT, "robustness_threshold_sensitivity.csv"), index=False)
pd.DataFrame(dict(perm_frac_negative=perm_frac, perm_mean_log2FC=perm_mean, perm_rho=perm_rho)).to_csv(
    os.path.join(OUT, "permutation_distributions.csv"), index=False)
print("\nDONE Part B")
