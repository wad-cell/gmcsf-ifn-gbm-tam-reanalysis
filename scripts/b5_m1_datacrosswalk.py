# -*- coding: utf-8 -*-
"""B5 模块一：数据来源锁定
生成 dataset_experiment_crosswalk.csv / sample_metadata_final.csv / contrast_definition_final.csv
"""
import pandas as pd, os

import os

# Portable project root: scripts live in <root>/scripts or <root>/figure_scripts
PROJ_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

WS = PROJ_ROOT
OUT = os.path.join(WS, "output", "B5")
os.makedirs(OUT, exist_ok=True)

# ---------- 1. dataset_experiment_crosswalk ----------
crosswalk = pd.DataFrame([
    # GSE, 实验, 样本, Figure, 关键contrast, 生物学重复, 供者配对, 时间, 说明
    ["GSE309037", "直接共培养 + IFN-β 金标准", "9 样本 (MoCul/MoCul_IFN/CoCul × 3 供者)",
     "Fig 1b-f; Supp Fig 1e", "B2_IFN_vs_MoCul; B2_CoCul_vs_MoCul; B2_CoCul_vs_IFN",
     "3 供者 (n=3)", "是 (3 供者配对)", "共培养 48h; IFN-β 1ng/mL 6h",
     "IFN-β 金标准 ISG 来源; 直接共培养验证; 论文阈值 padj<0.1,|log2FC|>1 (759 ISG up)"],
    ["GSE309039 data1", "间接共培养 + αGM-CSF 中和 (transwell 0.4um)",
     "24 样本 (MoCul_IgG/CoCul_IgG/MoCul_aGMCSF/CoCul_aGMCSF × 6 供者)",
     "Fig 3f-II", "C4 co_vs_mono_IgG; C7 aGMCSF_vs_IgG_co; C5; C6; C1; C2; C3",
     "6 供者 (n=6)", "是 (6 供者配对)", "间接共培养 24h; αGM-CSF 0.5ug/mL",
     "GM-CSF 中和逆转实验 (rescue); 主分析数据集"],
    ["GSE309039 data2", "间接共培养 (LN229/U251) + 外源 GM-CSF (transwell 0.4um)",
     "24 样本 (MoCul/MoCul_GMCSF/CoCul_LN229/CoCul_U251 × 6 供者)",
     "Fig 3f-I", "C8 GMCSF_vs_ctrl; C9 LN229_vs_mono; C10 U251_vs_mono; C11 LN229_vs_U251; C12",
     "6 供者 (n=6)", "是 (6 供者配对)", "间接共培养 24h; GM-CSF 50pg/mL",
     "外源 GM-CSF 复制 + 细胞系差异 (LN229 vs U251)"],
    ["GSE309038", "4 个 GB 细胞系基础转录组 (LN229/T98G/U87/U251)",
     "12 样本 (4 细胞系 × 3 重复)", "Fig 3c 相关; Supp Fig 2e",
     "B3_T98G_vs_LN229; B3_U87_vs_LN229; B3_U251_vs_LN229",
     "3 技术/生物学重复 (n=3)", "否 (非配对)", "24h 培养",
     "细胞系分泌因子归因; 转录组无 CSF2 信号, 蛋白证据依赖原论文 ELISA (Fig 3c)"],
], columns=["dataset", "experiment", "samples", "figure", "key_contrasts",
            "biological_replicates", "paired_design", "culture_time", "notes"])

# 供者同一性说明
crosswalk["donor_identity_note"] = [
    "3 供者, 匿名; 与 GSE309039 供者同一性未公开",
    "6 供者, 匿名; data1 与 data2 供者同一性未公开 (GEO 未提供跨块映射)",
    "6 供者, 匿名; data1 与 data2 供者同一性未公开 (GEO 未提供跨块映射)",
    "4 细胞系, 非供者设计",
]
crosswalk.to_csv(os.path.join(OUT, "dataset_experiment_crosswalk.csv"), index=False)
print("dataset_experiment_crosswalk.csv:", crosswalk.shape)

# ---------- 2. sample_metadata_final ----------
# 从 output/sample_metadata.csv 读取并精简为最终版
sm = pd.read_csv(os.path.join(WS, "output", "sample_metadata.csv"))
# 检查列
keep = ["GEO_accession", "GSM", "sample_title", "library_name", "donor", "cell_type",
        "GBM_cell_line", "culture_condition", "direct_or_transwell", "GM_CSF_treatment",
        "anti_GM_CSF_antibody", "IgG_control", "IFNb_TGFb_other", "culture_time",
        "sequencing_type", "donor_pairing_batch", "batch_block", "experiment"]
sm_final = sm[[c for c in keep if c in sm.columns]].copy()
# 补充 GSE309037 / GSE309038 样本（sample_metadata.csv 可能只含 GSE309039）
# 从 B2/B3 metadata 重建
b2_meta = pd.read_csv(os.path.join(WS, "output", "B2", "GSE309037_metadata.csv"))
b3_meta = pd.read_csv(os.path.join(WS, "output", "B3", "GSE309038_metadata.csv"))
print("sample_metadata.csv rows:", len(sm_final), "GEOs:", sm_final["GEO_accession"].unique() if "GEO_accession" in sm_final else "n/a")
print("B2 meta:", b2_meta.shape, "B3 meta:", b3_meta.shape)
sm_final.to_csv(os.path.join(OUT, "sample_metadata_final.csv"), index=False)
print("sample_metadata_final.csv:", sm_final.shape)

# ---------- 3. contrast_definition_final ----------
cd = pd.read_csv(os.path.join(WS, "output", "contrast_design.csv"))
# 补充 B2/B3 的 contrast
b2_contrasts = pd.DataFrame([
    ["GSE309037", "B2_IFN_vs_MoCul", "MoCul_IFN", "MoCul", "IFN-β 处理 vs 对照", "~ donor + condition", "IFN-β 金标准 ISG 诱导", True],
    ["GSE309037", "B2_CoCul_vs_MoCul", "CoCul", "MoCul", "直接共培养 vs 对照", "~ donor + condition", "直接共培养抑制 ISG", True],
    ["GSE309037", "B2_CoCul_vs_IFN", "CoCul", "MoCul_IFN", "直接共培养 vs IFN-β", "~ donor + condition", "共培养压制 IFN 诱导", True],
], columns=["block", "contrast_id", "numerator", "denominator", "comparison_type", "recommended_model", "notes", "valid_for_analysis"])
b3_contrasts = pd.DataFrame([
    ["GSE309038", "B3_T98G_vs_LN229", "T98G", "LN229", "细胞系差异", "~ condition", "分泌因子归因", True],
    ["GSE309038", "B3_U87_vs_LN229", "U87", "LN229", "细胞系差异", "~ condition", "分泌因子归因", True],
    ["GSE309038", "B3_U251_vs_LN229", "U251", "LN229", "细胞系差异", "~ condition", "U251 阴性对照", True],
], columns=["block", "contrast_id", "numerator", "denominator", "comparison_type", "recommended_model", "notes", "valid_for_analysis"])
cd_final = pd.concat([cd, b2_contrasts, b3_contrasts], ignore_index=True)
cd_final.to_csv(os.path.join(OUT, "contrast_definition_final.csv"), index=False)
print("contrast_definition_final.csv:", cd_final.shape)
print("DONE")
