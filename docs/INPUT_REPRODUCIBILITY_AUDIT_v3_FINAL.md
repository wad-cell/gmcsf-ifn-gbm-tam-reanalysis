# INPUT_REPRODUCIBILITY_AUDIT_v3_FINAL.md - input chain of the release candidate (C4.4G, sealed in C4.4H)

This audit closes the v2 gap "series matrix supplied by authors": for every
external input it now names the exact file, the exact accession, the download
location, the parsing rules and the verified checksum. Nothing here is inferred
from a file name; all checksums were computed on the files present in the
project workspace on 2026-09-10.

## 1. GSE309037 / GSE309038 / GSE309039 - full input chain

Each DESeq2 stage reads one author-processed count matrix published as a GEO
**series supplementary file**. There is **no decompression and no format
conversion step**: the files are already plain-text CSV inside gzip.

Reader used by the scripts (`compression="gzip"`, `sep=";"`, `index_col=0`):
```python
df = pd.read_csv(<file>, compression="gzip", sep=";", index_col=0)
```

| accession | GEO supplementary file | read by | bytes | sha256 |
|---|---|---|---|---|
| GSE309037 | GSE309037_raw_counts.csv.gz | scripts/08_b2_deseq2.py | 508232 | 598c95060c49de1e0e5667cdcf2ecb9561fbc9fe08d7101bffb84595e8173bb0 |
| GSE309038 | GSE309038_raw_counts.csv.gz | scripts/12_b3_deseq2.py | 451756 | b943820976ab1954fa58235be94b7ee3505d89b9da0966320d3e1ac0cc9e4103 |
| GSE309039 | GSE309039_data1_raw_counts.csv.gz | scripts/02_b1_deseq2.py | 1005871 | b94cffe93670ff1bafd742139005763aae705f27d1c7336d3f938c92d782929e |
| GSE309039 | GSE309039_data2_raw_counts.csv.gz | scripts/02_b1_deseq2.py | 801308 | 08aae847c8e98e14fc6dacb9a0c51f40d0830c6b4580f7f2740589c60ff6c5ef |

Download URLs (GEO FTP mirror, HTTPS):
```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309037/suppl/GSE309037_raw_counts.csv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309038/suppl/GSE309038_raw_counts.csv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309039/suppl/GSE309039_data1_raw_counts.csv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309039/suppl/GSE309039_data2_raw_counts.csv.gz
```
Series pages: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309037`
(and GSE309038 / GSE309039).

Local destination used by the scripts and the verifier:
`temp/GSE309037/raw/`, `temp/GSE309038/raw/`, `temp/GSE309039/raw/`.

Verification tool added in C4.4G (`tools/fetch_geo_supplementary.py`, in the
release zip):
```bash
python tools/fetch_geo_supplementary.py                              # expected local paths
python tools/fetch_geo_supplementary.py --check-project --parse       # offline sha256 + parse
python tools/fetch_geo_supplementary.py --download <dir>              # fetch from GEO + verify
```
Executed 2026-09-10 with `--check-project --parse`: **4/4 files present,
4/4 sha256 OK, 4/4 parse OK** (first rows read successfully; observed column
counts on the first 5 rows: 9 / 12 / 24 / 24, matching the sample titles of the
original comparison designs).

## 2. GSE163120 (GSM4972211) - external scRNA

Two GEO sample supplementary files are consumed by b5_m7_scRNA_patient.py,
b6_m4/b6_m5/b6_m6 and b7_01_prep_data.py: the 10x matrix archive and its cell
annotation table (referenced in the scripts as `matrix.csv.gz` and
`annot.csv.gz`, ~41 MB combined). They are **NOT re-distributed** in this
package and were **NOT re-downloaded** in C4.4G; obtain them from the
supplementary section of
`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4972211`
(the file listing on that page is authoritative - no file name is asserted here
beyond the two names the scripts reference).

## 3. GSE182109 - bulk 10x archive

`GSE182109_RAW.tar` (~2.4 GB, 44 GSM runs), consumed by the C3.1 h5ad rebuild
and by `scripts/b6_m8_load_ndgbm.py`:
```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182109/suppl/GSE182109_RAW.tar
```
Series page: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182109`.
**NOT re-downloaded** in C4.4G (marked `NOT RE-RUN` in RUN_ORDER.md).

## 4. Generated intermediates (GENERATED_INTERMEDIATE)

`INPUT_MANIFEST_v3.csv` lists all ten generated intermediates with the producing
script, its own inputs and a sha256 of the materialised file:

| intermediate | producing script | sha256 |
|---|---|---|
| b6_bg_genes_annotated.csv | scripts/b6_prep_permutation_data.py | 12acfc9612e3c1f1c6bfe6a27eebf82089cf55b43c98965eeca403be57ae0e9e |
| gse163120_tam_expr.pkl | scripts/b7_01_prep_data.py | 512dd770dc22fbad91a3953325fe5956633e9f4cb18ab1c84052328a959bc538 |
| gse163120_tam_annot.csv | scripts/b7_01_prep_data.py | 76dfaebad953a4f18eaaf1d2edf9e29483c2c638fac996f7d6a12aaf4f58f90b |
| gse182109_myeloid_expr.pkl | scripts/b7_01_prep_data.py | dd5275df7608171ed20b615768d5280d77b386a616494701d8407b02d3a34309 |
| gse182109_myeloid_annot.csv | scripts/b7_01_prep_data.py | 3786446b229e29fa15b349782897c0d98d49d1a34d85606a9e79f57726a7a877 |
| patient_paired_GSE163120.csv | scripts/prep_patient_paired.py | a77d6e501ac9bdf8f43a548d25165077bcb20ef4c46a0e239c9ed67e32216e4b |
| patient_paired_GSE182109.csv | scripts/prep_patient_paired.py | d8d17f197fbfd4dd9eb82852307537ed2bd607765f4d5aa9042696b1541d792c |
| gse182109_myeloid_expr_aligned.pkl | **scripts/generate_gse182109_aligned_object.py (new in C4.4G)** | e9abe6a81ac7b7f0f77602de9d7c0040a1ab244d9d0a45cb662ea9b6e353577c |

## 5. The formerly unshipped alignment step (gap closed)

In v2, `temp/B7/gse182109_myeloid_expr_aligned.pkl` was described only by prose
("barcode alignment of gse182109_myeloid_expr.pkl to the b7_01 table") and its
generator was not in the package. C4.4G ships
`scripts/generate_gse182109_aligned_object.py`, which

- takes `temp/B7/gse182109_myeloid_expr.pkl` (60825 x 164, non-unique index) and
  `output/B6/GSE182109_myeloid_annot.csv` (59118 rows, 8 columns);
- fixes the rules: annotation row order is the baseline (R1), barcode
  intersection (R2), first-come-first-served consumption of duplicate barcodes
  (R3), shape and element-wise order verification (R4), hard error with exit
  code 2 on any unsatisfiable precondition (R5, no partial output);
- was **really executed** on 2026-09-10 and produced
  sha256 `e9abe6a81ac7b7f0f77602de9d7c0040a1ab244d9d0a45cb662ea9b6e353577c`,
  **byte-identical** to the frozen object in the workspace (also 59118 x 164,
  identical index/columns/values);
- was exercised on its failure paths: a violated row-count precondition and a
  missing input both aborted with exit code 2 and wrote nothing.

## 6. Honest limits

The four GSE309037/38/39 matrices and the GSM4972211 / GSE182109 downloads were
**not re-fetched** in C4.4G; their local copies were verified by checksum only.
Large public objects therefore stay `NOT RE-RUN`. All checksums above were
computed on the files actually present in the workspace, not copied from any
external description.
