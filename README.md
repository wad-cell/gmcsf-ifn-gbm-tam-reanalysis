# GM-CSF/IFN Project - Code and Reproducibility Package v1.0.0
Archived release v1.0.0: [https://doi.org/10.5281/zenodo.22866023](https://doi.org/10.5281/zenodo.22866023)
Analysis code and reproducibility package for the GM-CSF/IFN tumour-myeloid study.
Scope is strictly limited to the five whitelisted public datasets
**GSE309037, GSE309038, GSE309039, GSE163120, GSE182109**.
All out-of-scope datasets are intentionally absent from this package.

Version history:
- v1 (C4.4E): 8-script candidate package, OUTPUT_MANIFEST verified 65/65.
- v2 (C4.4F): complete analysis-chain package, 42 Python scripts
  (37 `scripts/` + 5 `figure_scripts/`), Windows absolute paths replaced by the
  portable `PROJ_ROOT` convention, machine-name / internal workflow tokens removed.
- v3 (C4.4G, this package): closes the four remaining gaps of v2 -
  (a) the B7 barcode-aligned object now has a shipped generator script
  (`scripts/generate_gse182109_aligned_object.py`) with frozen alignment rules and
  hard failure conditions, instead of a recipe-only note;
  (b) `mygene` added and the whole third-party surface re-audited
  (`environment.lock.txt`, docs/DEPENDENCY_AUDIT_v3_FINAL.md);
  (c) the GSE309037/38/39 input chain is fully specified - per-file GEO
  supplementary names, download URLs, delimiter, sha256 and a reproducible
  fetch/verify tool (`tools/fetch_geo_supplementary.py`);
  (d) all generator/provenance footnotes have been removed from the released
  copies of the reports.
- v1.0.0 (public release): adds the MIT License, machine-readable citation and
  Zenodo metadata, and the final provenance and audit documents. No analysis
  logic or frozen scientific result was changed.

## Directory layout
```
GMCSF_IFN_CODE_RELEASE_v1.0.0/
  README.md                     this file
  LICENSE                       MIT License
  CITATION.cff                  citation metadata for GitHub and reference tools
  .zenodo.json                  Zenodo release metadata
  requirements.txt              third-party dependencies (includes mygene)
  environment.lock.txt          per-package exact / best-evidenced version pins
  RUN_ORDER.md                  actual execution order (per analysis stage)
  SOFTWARE_VERSIONS.txt         validated software environment
  INPUT_MANIFEST_v3.csv         every external / intermediate input used by the chain
  scripts/                      38 analysis scripts (B1->B7, B7 aligned object,
                                C3.1 h5ad rebuild, B4d, C2 paired tables)
  figure_scripts/               5 figure/source-data renderers (incl. C2 Figure 5)
  tools/fetch_geo_supplementary.py
                                GSE309037/38/39 raw-count fetch + sha256/parse verifier
  provenance/                   output provenance and script inclusion decisions
  docs/                         final release-readiness and reproducibility audits
```

## Path-portability convention
Every script resolves the repository root as
`PROJ_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))`
and then uses `os.path.join(PROJ_ROOT, "output", ...)` /
`os.path.join(PROJ_ROOT, "temp", ...)`. No machine-specific drive paths remain.
Expected runtime tree (created by the scripts as needed):
```
<package>/output/B1..B7   analysis result tables (the 65 frozen outputs live in
                          the lab workspace and are NOT re-distributed in this
                          zip; the same applies to C3.1 and C2 figures)
<package>/temp/           generated intermediates (see INPUT_MANIFEST_v3.csv)
```

## B7 barcode-aligned object (new in v3)
`gse182109_myeloid_expr_aligned.pkl` used to be produced by an ad-hoc workspace
script that was never shipped. It is now a first-class pipeline step:

```bash
python scripts/generate_gse182109_aligned_object.py --project-root .
python scripts/generate_gse182109_aligned_object.py --check-only     # verify only
```
Inputs: `temp/B7/gse182109_myeloid_expr.pkl` (b7_01_prep_data.py) and
`output/B6/GSE182109_myeloid_annot.csv` (b6_m8_analysis_v2.py). The frozen rules
(row order = annotation order, barcode intersection, first-come-first-served
consumption of duplicate barcodes, hard error on an unsatisfiable barcode or a
shape mismatch) are documented in the script docstring. Verified in C4.4G by a
real re-run: the rebuilt object is byte-identical to the frozen one
(sha256 e9abe6a8...577c, see INPUT_MANIFEST_v3.csv).

## GSE309037 / GSE309038 / GSE309039 input chain
Each DESeq2 stage reads one author-processed count matrix published as a GEO
series supplementary file (gzip CSV, semicolon-separated, first column = gene,
read with `pd.read_csv(..., compression="gzip", sep=";", index_col=0)`):

| accession | GEO file | read by | sha256 (prefix) |
|---|---|---|---|
| GSE309037 | GSE309037_raw_counts.csv.gz | scripts/08_b2_deseq2.py | 598c9506... |
| GSE309038 | GSE309038_raw_counts.csv.gz | scripts/12_b3_deseq2.py | b9438209... |
| GSE309039 | GSE309039_data1_raw_counts.csv.gz | scripts/02_b1_deseq2.py | b94cffe9... |
| GSE309039 | GSE309039_data2_raw_counts.csv.gz | scripts/02_b1_deseq2.py | 08aae847... |

No decompression or format-conversion step is required. Full URLs, byte sizes
and the verification tool:
```bash
python tools/fetch_geo_supplementary.py                                  # list expected local paths
python tools/fetch_geo_supplementary.py --check-project --parse          # offline sha256 + parse check
python tools/fetch_geo_supplementary.py --download /tmp/geo309           # fetch from GEO + verify
```

## Dependencies
`requirements.txt` (floors, `mygene` included) and `environment.lock.txt`
(per-package version + evidence tier; packages without original evidence are
recorded as `UNKNOWN_ORIGINAL_VERSION` and deliberately not guessed:
`mygene`, `gseapy`, `seaborn`).

## Frozen-output provenance
`provenance/COMPLETE_PROVENANCE_MATRIX_v3.csv` maps all **65 frozen
outputs** to their v3 generating script, direct inputs, producing scripts and
public upstream roots. Coverage = 65/65.

## Execution status - honest statement
- 9 scripts carry a clean-room PASS from C4.4E (b5_figures, b5_m34, b5_m6b,
  b6_permutation_audit, b7_06_rescue, b7_08_figures, make_figure5_v3_1,
  make_main_figures, make_supp_figures).
- New in C4.4G: `generate_gse182109_aligned_object.py` was really executed twice
  (rebuild + check-only) and its negative paths were exercised (missing input and
  a violated row-count precondition both abort with exit code 2).
- Every shipped Python file was syntax-checked (`py_compile`), and the import
  surface was re-audited against `requirements.txt` / `environment.lock.txt`.
- The bulk statistics stages that need the very large public downloads
  (GSE182109_RAW.tar ~2.4 GB, GSM4972211 matrices ~41 MB) were **NOT re-run from
  public data**; they remain marked `NOT RE-RUN`. Nothing in this package claims
  a full from-scratch rerun.

## Quick start
```bash
python -m venv .venv && .venv/Scripts/activate   # Windows
pip install -r requirements.txt                  # versions: environment.lock.txt
python scripts/02_b1_deseq2.py                   # then follow RUN_ORDER.md
```
No R is required; the run environment is Python 3.11.8 / PyDESeq2 0.5.4
(R-equivalent .R files exist only as paper-workflow references, not runners).

## License

Copyright (c) 2026 Jing Deng, Ying Wang, Jiateng Zeng, Lianghong Yu, and
Hongliang Ge. This software is released under the MIT License; see `LICENSE`.

## Citation and archival record

Citation metadata are provided in `CITATION.cff`. The GitHub release is intended
to be archived through Zenodo as version 1.0.0. After Zenodo publishes the
record, cite the version-specific DOI shown on the Zenodo record page.
