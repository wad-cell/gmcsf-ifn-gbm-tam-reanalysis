# -*- coding: utf-8 -*-
"""B5 模块三+四：核心机制链定量验证 (A-F) + 逆转分析
输出:
  - mechanism_chain_validation.csv (A-F 每步基因集层面统计)
  - rescue_analysis_full.csv
  - suppression_rescue_scatter.pdf
  - rescue_heatmap.pdf
"""
import pandas as pd, numpy as np, os
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
WS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(WS, "output", "B5")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)
rng = np.random.default_rng(42)

ev = pd.read_csv(os.path.join(OUT, "core68_evidence_matrix.csv"))
syms = ev["symbol"].tolist()

def col(cid, stat): return f"{cid}_{stat}"

def bootstrap_mean(x, n=2000, seed=42):
    r = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float); x = x[~np.isnan(x)]
    if len(x) < 3: return (np.nan, np.nan, np.nan)
    m = np.mean(x)
    bs = [np.mean(r.choice(x, len(x), replace=True)) for _ in range(n)]
    return m, np.percentile(bs, 2.5), np.percentile(bs, 97.5)

def bootstrap_spearman(x, y, n=2000, seed=42):
    r = np.random.default_rng(seed)
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y)); x, y = x[m], y[m]
    if len(x) < 5: return (np.nan, np.nan, np.nan, np.nan)
    rho, p = stats.spearmanr(x, y)
    bs = []
    for _ in range(n):
        idx = r.choice(len(x), len(x), replace=True)
        bs.append(stats.spearmanr(x[idx], y[idx])[0])
    return rho, p, np.percentile(bs, 2.5), np.percentile(bs, 97.5)

def gene_set_stats(lfc, padj, direction="down", label=""):
    """基因集层面统计: 平均效应/方向一致率/95%CI/单样本检验"""
    lfc = np.asarray(lfc, float); padj = np.asarray(padj, float)
    m = ~(np.isnan(lfc) | np.isnan(padj)); lfc, padj = lfc[m], padj[m]
    n = len(lfc)
    if n == 0: return dict(label=label, n=0)
    mean, lo, hi = bootstrap_mean(lfc)
    # 方向一致率: 与预期方向一致的基因比例 (含显著)
    if direction == "down":
        concord = np.mean(lfc < 0)
        sig = np.mean((padj < 0.05) & (lfc < 0))
    else:
        concord = np.mean(lfc > 0)
        sig = np.mean((padj < 0.05) & (lfc > 0))
    # 单样本 t 检验 (H0: mean=0)
    t, p = stats.ttest_1samp(lfc, 0)
    # 符号检验
    sign_p = stats.binomtest(int(concord * n), n, 0.5).pvalue
    return dict(label=label, n=n, mean_log2FC=mean, ci_lo=lo, ci_hi=hi,
                concordance=concord, sig_frac=sig, tstat=t, ttest_p=p, sign_p=sign_p)

rows = []

# ============ A. IFN-β 诱导 (B2_IFN_vs_MoCul) ============
rows.append(gene_set_stats(ev[col("B2_IFN_vs_MoCul","log2FC")], ev[col("B2_IFN_vs_MoCul","padj")],
                           direction="up", label="A_IFNb_induction"))

# ============ B. LN229 共培养抑制 ============
# B-direct: GSE309037 直接共培养
rows.append(gene_set_stats(ev[col("B2_CoCul_vs_MoCul","log2FC")], ev[col("B2_CoCul_vs_MoCul","padj")],
                           direction="down", label="B_LN229_direct_CoCul"))
# B-indirect: GSE309039 data2 间接共培养 (C9)
rows.append(gene_set_stats(ev[col("C9_LN229_vs_mono","log2FC")], ev[col("C9_LN229_vs_mono","padj")],
                           direction="down", label="B_LN229_indirect_CoCul"))

# ============ C. U251 阴性对照 ============
# C10 U251_vs_mono
rows.append(gene_set_stats(ev[col("C10_U251_vs_mono","log2FC")], ev[col("C10_U251_vs_mono","padj")],
                           direction="down", label="C_U251_indirect_CoCul"))
# LN229 vs U251 配对比较 (C9 vs C10, 68 基因)
l9 = ev[col("C9_LN229_vs_mono","log2FC")].values; l10 = ev[col("C10_U251_vs_mono","log2FC")].values
m = ~(np.isnan(l9) | np.isnan(l10))
w, wp = stats.wilcoxon(l9[m], l10[m])
rows.append(dict(label="C_LN229_vs_U251_wilcoxon", n=int(m.sum()),
                 mean_log2FC=np.mean(l9[m]-l10[m]), tstat=w, ttest_p=wp,
                 concordance=np.mean(l9[m] < l10[m]), sign_p=np.nan))

# ============ D. 外源 GM-CSF 复制 (C8 vs C9) ============
l8 = ev[col("C8_GMCSF_vs_ctrl","log2FC")].values; l9b = ev[col("C9_LN229_vs_mono","log2FC")].values
rho, rp, rlo, rhi = bootstrap_spearman(l8, l9b)
rows.append(dict(label="D_GMCSF_vs_LN229_spearman", n=int((~(np.isnan(l8)|np.isnan(l9b))).sum()),
                 mean_log2FC=rho, ci_lo=rlo, ci_hi=rhi, tstat=rp, ttest_p=rp,
                 concordance=np.mean(np.sign(l8)==np.sign(l9b)), sign_p=np.nan))
rows.append(gene_set_stats(l8, ev[col("C8_GMCSF_vs_ctrl","padj")], direction="down", label="D_GMCSF_effect"))

# ============ E. αGM-CSF 逆转 (C7) ============
rows.append(gene_set_stats(ev[col("C7_aGMCSF_vs_IgG_co","log2FC")], ev[col("C7_aGMCSF_vs_IgG_co","padj")],
                           direction="up", label="E_aGMCSF_rescue"))

# ============ F. 直接 vs 间接一致性 (B2_CoCul vs C9) ============
ld = ev[col("B2_CoCul_vs_MoCul","log2FC")].values; li = ev[col("C9_LN229_vs_mono","log2FC")].values
rhoF, rpF, rloF, rhiF = bootstrap_spearman(ld, li)
# 线性回归斜率
mF = ~(np.isnan(ld) | np.isnan(li))
slope, intercept, rv, pv, se = stats.linregress(ld[mF], li[mF])
# 离群基因 (残差 > 2SD)
resid = li[mF] - (slope*ld[mF] + intercept)
outliers = np.where(np.abs(resid) > 2*np.std(resid))[0]
out_syms = [syms[i] for i in np.where(mF)[0][outliers]]
rows.append(dict(label="F_direct_vs_indirect_spearman", n=int(mF.sum()),
                 mean_log2FC=rhoF, ci_lo=rloF, ci_hi=rhiF, tstat=rpF, ttest_p=rpF,
                 concordance=np.mean(np.sign(ld[mF])==np.sign(li[mF])), sign_p=np.nan,
                 slope=slope, intercept=intercept, n_outliers=len(outliers)))

mech = pd.DataFrame(rows)
mech.to_csv(os.path.join(OUT, "mechanism_chain_validation.csv"), index=False)
print("=== 机制链验证 ===")
print(mech.to_string(index=False))
print("F 离群基因:", out_syms)

# ============ 模块四: 逆转分析 ============
# suppression_effect = C4_co_vs_mono_IgG (data1, 与 rescue 同供者块)
# rescue_effect = C7_aGMCSF_vs_IgG_co (data1)
sup = ev[col("C4_co_vs_mono_IgG","log2FC")].values
res = ev[col("C7_aGMCSF_vs_IgG_co","log2FC")].values
sup_p = ev[col("C4_co_vs_mono_IgG","padj")].values
res_p = ev[col("C7_aGMCSF_vs_IgG_co","padj")].values

# 负相关检验
m = ~(np.isnan(sup) | np.isnan(res))
rho_r, p_r, rlo_r, rhi_r = bootstrap_spearman(sup, res)
print(f"\n=== 逆转分析: suppression vs rescue Spearman rho={rho_r:.3f} (95%CI {rlo_r:.3f}-{rhi_r:.3f}), p={p_r:.2e}")

# rescue index: RI = rescue_effect / (-suppression_effect), 仅对 suppression<-0.1 的基因
# 预定义阈值
SUP_TH = -0.1   # 有实质抑制
RI_FULL = 0.8   # fully rescued
RI_PART = 0.2   # partially rescued 下限
df = pd.DataFrame(dict(symbol=syms, suppression_effect=sup, rescue_effect=res,
                       suppression_padj=sup_p, rescue_padj=res_p))
df["rescue_index"] = np.where(df["suppression_effect"] < SUP_TH,
                              df["rescue_effect"] / (-df["suppression_effect"]), np.nan)
df["rescue_index_ci_lo"] = np.nan; df["rescue_index_ci_hi"] = np.nan

def classify(r):
    s, ri, rp = r["suppression_effect"], r["rescue_index"], r["rescue_padj"]
    if pd.isna(s) or pd.isna(ri):
        return "na"
    if s > 0:  # 共培养中上调 (与主方向相反)
        return "discordant"
    if ri >= RI_FULL:
        return "fully_rescued"
    if ri >= RI_PART:
        return "partially_rescued"
    if ri < RI_PART:
        return "not_rescued"
    return "na"
df["rescue_class"] = df.apply(classify, axis=1)
# 统计显著恢复 (rescue_padj<0.05 且 rescue_effect>0)
df["rescue_significant"] = (df["rescue_padj"] < 0.05) & (df["rescue_effect"] > 0)

# 敏感性分析: 不同阈值组合
sens = []
for sup_th in [-0.05, -0.1, -0.5, -1.0]:
    for ri_full in [0.5, 0.8, 1.0]:
        for ri_part in [0.1, 0.2, 0.3]:
            if ri_part >= ri_full: continue
            ri = np.where(df["suppression_effect"] < sup_th, df["rescue_effect"]/(-df["suppression_effect"]), np.nan)
            cls = []
            for s, rv in zip(df["suppression_effect"], ri):
                if pd.isna(s) or pd.isna(rv): cls.append("na")
                elif s > 0: cls.append("discordant")
                elif rv >= ri_full: cls.append("fully_rescued")
                elif rv >= ri_part: cls.append("partially_rescued")
                else: cls.append("not_rescued")
            sens.append(dict(suppression_threshold=sup_th, ri_full=ri_full, ri_part=ri_part,
                             n_analyzable=int((df["suppression_effect"]<sup_th).sum()),
                             fully=int(np.sum(np.array(cls)=="fully_rescued")),
                             partially=int(np.sum(np.array(cls)=="partially_rescued")),
                             not_rescued=int(np.sum(np.array(cls)=="not_rescued")),
                             discordant=int(np.sum(np.array(cls)=="discordant"))))
sens_df = pd.DataFrame(sens)
sens_df.to_csv(os.path.join(OUT, "rescue_sensitivity_analysis.csv"), index=False)
print("\n=== 敏感性分析 (默认阈值 sup<-0.1, full>=0.8, part>=0.2) ===")
print(df["rescue_class"].value_counts())
print("\n敏感性分析表:")
print(sens_df.to_string(index=False))

df.to_csv(os.path.join(OUT, "rescue_analysis_full.csv"), index=False)

# ============ 图 1: suppression vs rescue scatter ============
plt.figure(figsize=(7, 6))
m2 = ~(np.isnan(sup) | np.isnan(res))
colors = {"fully_rescued":"#2ca02c","partially_rescued":"#ff7f0e","not_rescued":"#d62728","discordant":"#9467bd","na":"#bbbbbb"}
for cls in colors:
    sel = (df["rescue_class"]==cls) & m2
    plt.scatter(sup[sel], res[sel], c=colors[cls], s=40, alpha=0.8, label=cls, edgecolors="none")
# 参考线: 完全恢复 (res = -sup)
xx = np.linspace(-4, 1, 100)
plt.plot(xx, -xx, "--", color="grey", lw=1, label="full rescue (res=-sup)")
plt.axhline(0, color="black", lw=0.5); plt.axvline(0, color="black", lw=0.5)
plt.xlabel("suppression_effect (log2FC, CoCul_IgG vs MoCul_IgG)")
plt.ylabel("rescue_effect (log2FC, αGM-CSF vs IgG in CoCul)")
plt.title(f"Suppression vs rescue of 68 core ISGs\nSpearman ρ={rho_r:.3f} (95%CI {rlo_r:.3f}–{rhi_r:.3f}), p={p_r:.2e}")
plt.legend(fontsize=8, loc="lower right")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "suppression_rescue_scatter.pdf")); plt.close()

# ============ 图 2: rescue heatmap ============
# 按 suppression 排序, 展示 suppression/rescue log2FC 与分类
df2 = df.dropna(subset=["suppression_effect"]).sort_values("suppression_effect")
hm = df2[["suppression_effect","rescue_effect"]].copy()
hm.columns = ["suppression","rescue"]
plt.figure(figsize=(4.5, 10))
sns.heatmap(hm.T, cmap="RdBu_r", center=0, vmin=-3, vmax=3, cbar_kws={"label":"log2FC"},
            yticklabels=True, xticklabels=df2["symbol"], linewidths=0.3)
plt.yticks(rotation=0)
plt.xticks(rotation=90, fontsize=6)
plt.title("68 core ISGs: suppression (CoCul_IgG) vs rescue (αGM-CSF)")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "rescue_heatmap.pdf")); plt.close()

print("\nDONE 模块三+四")
