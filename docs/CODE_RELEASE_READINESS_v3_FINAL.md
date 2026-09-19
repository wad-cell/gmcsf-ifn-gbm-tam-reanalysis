# CODE_RELEASE_READINESS_v3_FINAL.md - GM-CSF/IFN code release candidate (C4.4G, sealed in C4.4H)

**Verdict: GO** (for offline archiving and manuscript code-availability deposit;
license still undecided by explicit instruction, and no GitHub/Zenodo/DOI action
was taken).

## 1. What C4.4G closed (the four v2 gaps)

| # | gap in v2 | C4.4G result |
|---|---|---|
| 1 | `gse182109_myeloid_expr_aligned.pkl` had no shipped generator (recipe only) | **CLOSED** - `scripts/generate_gse182109_aligned_object.py` added (frozen rules R1-R5, hard failure, `--check-only`, `--verify-sha256`, `--expect-rows`); documented in README + RUN_ORDER step 28a. Really executed: rebuilt object **byte-identical** to the frozen one, sha256 `e9abe6a81ac7b7f0f77602de9d7c0040a1ab244d9d0a45cb662ea9b6e353577c`; negative paths exit with code 2 |
| 2 | `requirements` missing `mygene`; no full import audit | **CLOSED** - `mygene>=1.2.0` added; all 44 shipped `.py` files (38 `scripts/` + 5 `figure_scripts/` + 1 `tools/`; the `excluded_scripts/` copy is EXCLUDED_REFERENCE_ONLY and lives outside the ZIP, so it is not counted as a shipped Python file) AST-scanned: 10 third-party modules, all listed in requirements.txt; `environment.lock.txt` gives per-package evidence tiers; DEPENDENCY_AUDIT_v3.md documents the surface. 3 packages have no original evidence and are labelled **UNKNOWN_ORIGINAL_VERSION** (mygene, gseapy, seaborn) - not guessed |
| 3 | GSE309037/38/39 input chain said only "series matrix supplied by authors" | **CLOSED** - the four GEO supplementary files are named individually with accession, HTTPS URL, byte size, delimiter (`;`), reader call, gzip handling (no conversion step needed) and sha256; `tools/fetch_geo_supplementary.py` ships with download + offline verify modes; offline check on 2026-09-10 gave 4/4 sha256 OK and 4/4 parse OK |
| 4 | generator/provenance footnotes in the public reports | **CLOSED within the release copy** - a token scan (the generator/watermark token family and the AI foot-note string) over all 50 zip entries and over all sibling deliverables of this sealed set returns **0 hits**. Historical v1/v2 audit files were left untouched as instructed (two of them, the GSE309037/38/39 joint-analysis plan and the GSE309039 data-audit report, still carry such foot-notes; they are pre-existing workspace audit documents and were deliberately not edited) |

## 2. Re-executed checks (real runs, 2026-09-10)

| check | scope | result |
|---|---|---|
| `py_compile` | all 44 shipped `.py` (38 `scripts/` + 5 `figure_scripts/` + 1 `tools/`; EXCLUDED_REFERENCE_ONLY copy not counted) | **44/44 PASS, 0 FAIL** |
| import completeness | 10 third-party modules via `find_spec` | **10/10 resolvable, 0 missing** |
| aligned object rebuild | real re-run of the new generator from the small local inputs | **SUCCESS**, sha256 `e9abe6a8...577c`, byte-identical to the frozen object |
| generator failure paths | missing input / violated `--expect-rows` | **both abort with exit code 2, no partial output** |
| GSE309037/38/39 input chain | 4 files, offline | **4/4 present, 4/4 sha256 OK, 4/4 parse OK** |
| provenance coverage | 65 frozen outputs of the manuscript | **65/65 mapped**, each row carries the v3 generating script, its inputs and the execution evidence |
| portability | drive-letter / machine-path scan of the release copy | **0 hits** |

Script inventory: **51 rows** in SCRIPT_INCLUSION_EXCLUSION_v3.csv - 44 INCLUDED
(38 `scripts/` + 5 `figure_scripts/` + 1 tool), 6 EXCLUDED_SUPERSEDED, 1
EXCLUDED_REFERENCE_ONLY (R equivalents). Superseded scripts are recorded, not
deleted; the superseded `10_b2_overlap.py` reference copy stays outside the zip.

## 3. Package contents and checksums

- `GMCSF_IFN_CODE_RELEASE_CANDIDATE_v3_FINAL.zip` - **50 entries**, 123469 bytes,
  sha256 `f03590da577a389d428dfe97b51e486f4e426e822541024999fca3791cbc064c`
- the ZIP does **not** contain its own checksum file: `FINAL_CHECKSUMS_CODE_v3.sha256`
  is not a member of the archive (removed from the release copy; the archive never
  carries a hash of itself)
- sibling deliverables: COMPLETE_PROVENANCE_MATRIX_v3.csv,
  SCRIPT_INCLUSION_EXCLUSION_v3.csv, INPUT_MANIFEST_v3.csv,
  DEPENDENCY_AUDIT_v3_FINAL.md, INPUT_REPRODUCIBILITY_AUDIT_v3_FINAL.md,
  CODE_RELEASE_READINESS_v3_FINAL.md.
  See FINAL_CHECKSUMS_CODE_v3_FINAL.sha256 for their hashes.

## 4. Honest limits (unchanged from what the evidence supports)

- The very large public downloads (GSE182109_RAW.tar ~2.4 GB, GSM4972211
  matrices ~41 MB) were **not re-fetched** and the bulk statistics stages were
  **not re-run from public data**; they are marked `NOT RE-RUN` in RUN_ORDER.md
  and in the provenance matrix. No test pass is claimed for them.
- Three dependency versions remain `UNKNOWN_ORIGINAL_VERSION` (mygene, gseapy,
  seaborn) because the original environment left no evidence; the session host
  resolves them to mygene 3.2.2 / gseapy 1.3.1 / seaborn 0.13.2, which is
  recorded as host observation only, and statsmodels resolves to 0.15.0 versus
  the project clean-room log's 0.14.6.
- License: **not decided** (per instruction). LICENSE_RECOMMENDATION.md (v2)
  remains a comparison only (MIT / BSD-3-Clause / GPL-3.0; all dependencies used
  are permissive).
- No GitHub, Zenodo or journal submission action was taken.

## 5. Compliance with the write-safety boundary

Every artefact of C4.4G was created inside
`output/GMCSF_IFN_CODE_RELEASE_CANDIDATE_v3/` or in the `output/` root as a new
v3 file. **No original project file was deleted, moved, overwritten, renamed or
sent to the recycle bin**; the v1/v2 deliverables were read only; the frozen
review drafts, abstract candidates, author sign-off sheets and the isolated
superseded-output directory were never opened or touched.
