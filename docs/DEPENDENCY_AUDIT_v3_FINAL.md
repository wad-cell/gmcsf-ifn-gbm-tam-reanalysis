# DEPENDENCY_AUDIT_v3_FINAL.md - third-party dependency surface of the release candidate (C4.4G, sealed in C4.4H)

Scope of the audit: every `.py` file in `GMCSF_IFN_CODE_RELEASE_CANDIDATE_v3/`
(38 `scripts/`, 5 `figure_scripts/`, 1 `tools/` utility) = **44 shipped
Python files**. The single `excluded_scripts/10_b2_overlap.py` copy is
EXCLUDED_REFERENCE_ONLY, is **not** a member of the ZIP and is therefore **not
counted** as a shipped Python file. Imports were extracted with an AST scan
(`ast.Import` / `ast.ImportFrom`), so only real imports count, not mentions in
comments or docstrings.

## 1. Third-party modules actually imported

| module | shipped files importing it | consumers (examples) | in requirements.txt |
|---|---|---|---|
| pandas | 44 | all analysis and figure scripts (excluding excluded_scripts/ reference copy) | yes (`pandas>=2.1`) |
| numpy | 35 | all analysis and figure scripts | yes (`numpy>=1.26`) |
| scipy | 19 | B1-B4 statistics, B5/B6 robustness | yes (`scipy>=1.11`) |
| matplotlib | 8 | figure_scripts/*, b5_m34 | yes (`matplotlib>=3.8`) |
| pydeseq2 | 4 | 02_b1_deseq2, 08_b2_deseq2, 12_b3_deseq2, 03_b1_summary | yes (`pydeseq2==0.5.4`) |
| gseapy | 3 | b5_m5_enrichment, b7_02b_decomposition, b7_scoring | yes (`gseapy>=1.1`) |
| statsmodels | 2 | b7_04_gmcsf_ifn_corr, b7_06_rescue | yes (`statsmodels>=0.14`) |
| anndata | 1 | a3_rebuild_fulltranscriptome | yes (`anndata>=0.10`) |
| matplotlib-backed seaborn | 1 | b5_m34_mechanism_rescue | yes (`seaborn>=0.13`) |
| mygene | 1 | 04_b1_annotate | **added in C4.4G** (`mygene>=1.2.0`) |

Local package module: `b7_scoring` (shipped as `scripts/b7_scoring.py`) is
imported by `scripts/b7_02b_decomposition.py` and `scripts/b7_04_gmcsf_ifn_corr.py`
- it is not a third-party dependency.

`concurrent`, `os`, `sys`, `time`, `json`, `re`, `glob`, `subprocess`, `warnings`,
`itertools`, `collections`, `math`, `pickle` are standard-library modules and
carry no pin. `h5py` is pulled in transitively by `anndata` (h5ad backend) and is
not imported by any shipped script.

## 2. Version evidence (see `environment.lock.txt` for the machine-readable list)

| package | version | evidence tier |
|---|---|---|
| python | 3.11.8 | FROZEN_EXACT (analysis environment record) |
| pydeseq2 | 0.5.4 | FROZEN_EXACT (frozen release record) |
| matplotlib | 3.11.1 | ARTIFACT_DERIVED (metadata inside the original-run figure files) + CLEANROOM_OBSERVED |
| numpy / pandas / scipy / anndata | 2.4.6 / 2.3.3 / 1.17.1 / 0.12.19 | CLEANROOM_OBSERVED (2026-08-26 project clean-room log) |
| statsmodels | 0.14.6 (clean-room log) | CLEANROOM_OBSERVED |
| mygene | **UNKNOWN_ORIGINAL_VERSION** | no original evidence - not guessed |
| gseapy | **UNKNOWN_ORIGINAL_VERSION** | no original evidence - not guessed |
| seaborn | **UNKNOWN_ORIGINAL_VERSION** | no original evidence - not guessed (the clean-room log records seaborn as NOT_INSTALLED, so its PASS run must have used another environment) |

`UNKNOWN_ORIGINAL_VERSION` list (complete): **mygene, gseapy, seaborn**.
No version for these three is invented anywhere in this package.

## 3. Import-completeness check (executed 2026-09-10)

`importlib.util.find_spec` was run for all ten third-party modules in the
session environment: **10/10 resolvable, 0 missing**. Observed versions in that
environment (for transparency - they are NOT the original-environment versions):
numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, matplotlib 3.11.1, statsmodels 0.15.0,
gseapy 1.3.1, anndata 0.12.19, seaborn 0.13.2, mygene 3.2.2, pydeseq2 0.5.4,
h5py 3.16.0. Note that statsmodels resolves to 0.15.0 here while the project
clean-room log records 0.14.6 - further evidence that a host environment cannot
substitute for the original record.

## 4. Syntax check

`py_compile` over all 44 shipped `.py` files: **44/44 PASS, 0 FAIL**.

## 5. Install commands

```bash
pip install -r requirements.txt                              # floors, mygene included
pip install pydeseq2==0.5.4 matplotlib==3.11.1               # the two artifact-verified pins
```
For the three `UNKNOWN_ORIGINAL_VERSION` packages install any current release;
pin them yourself if you need a byte-stable environment.
