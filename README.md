# 05_replication — External replication in Greene 2025 (7,837 keloid cases)

Complete replication pipeline for the manuscript:

> **No causal association between circulating interleukin-6 and keloid risk: an adequately
> powered two-sample Mendelian randomization study with non-replication of an initial nominal
> association**
> Yong-fei Wang, Jiang-bo Kuang

## Outcome dataset

| Trait | Accession | Source | Cases / controls |
|---|---|---|---|
| Keloid (multi-ancestry meta-analysis) | GCST90652487 | Greene et al., *Nat Commun* 2025;16:7770 (PMID 40835838, PMCID PMC12368108) | 7,837 / 1,593,009 |

Not redistributed here — download the GWAS-SSF file from the GWAS Catalog by accession.
Expected file size 340 MB; **md5 `8db7428b31874b9f9a509a4e0e6ff69c`** (checked in the script).

### Contributing cohorts (from Greene 2025, Table 1)

| Cohort | EUR cases/controls | EAS cases/controls | AFR cases/controls |
|---|---|---|---|
| BioVU | 163 / 47,047 | — | 111 / 10,605 |
| eMERGE | 260 / 43,679 | — | 97 / 8,395 |
| VUMC | — | — | 122 / 356 |
| Million Veteran Program | 2,112 / 451,944 | — | 2,366 / 117,486 |
| UK Biobank | 257 / 413,923 | — | — |
| FinnGen | 1,294 / 321,903 | — | — |
| Biobank Japan | — | 1,055 / 177,671 | — |
| **Total** | **4,086 / 1,278,496** | **1,055 / 177,671** | **2,696 / 136,842** |

Sum = **7,837 cases / 1,593,009 controls**, which reconciles exactly with the paper abstract and
the GWAS Catalog summary.

> The GWAS Catalog `cohort` annotation for this study also lists `ACCC`. That label does **not**
> correspond to any cohort in the paper (0 occurrences in the full text, and the case/control
> totals are fully accounted for by the seven cohorts above), so it is not reproduced here.
> Greene et al. used **All of Us** (3,371 cases / 288,438 controls) for *independent replication*;
> that data is **not** part of GCST90652487 and is therefore not used in this analysis.

## Scripts

| File | Purpose |
|---|---|
| `IL6_Keloid_Greene2025_replication.R` | **Main script.** Five analysis sets (all 23 SNPs / IL6R only / HLA only / all minus lead / IL6R minus lead) × seven estimators; Cochran's Q + I²; MR-Egger intercept; MR-PRESSO (1,000 permutations, fixed seed); within-locus leave-one-out; Sakaue-vs-Greene side-by-side table; approximate power with a predicted-vs-observed SE self-check; **Steiger directionality test for the replication dataset**. |
| `figure_replication_forest.R` | Figure 7 (forest plot, discovery vs replication), 300 dpi TIF/PNG. |
| `figure_presso_steiger.R` | Standalone re-run of MR-PRESSO per stratum and the Steiger test (useful for auditing). |

### Three script fixes applied on 2026-09-26

1. **MR-PRESSO had never actually run.** The call passed `OUTCOME = "outcome"`. The signature of
   `MRPRESSO::mr_presso()` in MRPRESSO 1.0 is
   `mr_presso(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, OUTLIERtest, DISTORTIONtest,`
   `SignifThreshold, NbDistribution, seed)` — **there is no `OUTCOME` argument** (it exists only in
   later releases), so the call failed with "unused argument". The old
   `error = function(e) NULL` swallowed the error, so the replication MR-PRESSO step produced
   neither console output nor a results file. The argument has been removed, `seed` added, and the
   silent `NULL` replaced with `warning()`.
2. **Replication Steiger test was missing.** Added, using `get_r_from_lor()` to supply
   `r.outcome` for the binary outcome.
3. **Power table: two corrections.** (a) `n_ctrl` was 1,585,173, back-derived from the SSF `N`
   column (1,593,010), which contradicts the paper's own Table 1; it is now **1,593,009**.
   (b) Both datasets previously shared the 23-SNP cumulative R² (0.041407); the replication used
   only 22 SNPs (rs5004279 dropped as palindromic, R² = **0.039973**), and each dataset now uses
   its own R². Net effect: the replication predicted SE changes from 0.0278 to **0.0283**; all
   quoted power percentages are unchanged (discovery 48% / 99% / 52%; replication 97% / 53%).

> Note when recomputing R²: `table.harmonised_Greene.csv` has **no `samplesize.exposure` column**;
> `N_exp = 21758` must be supplied explicitly or `sum()` returns 0 and `Predicted_SE` becomes `Inf`.

## Outputs

| File | Content |
|---|---|
| `outputs/table.Greene_MR_all_methods_2026-09-24.csv` | All methods, all five strata (replication) |
| `outputs/table.Replication_Sakaue_vs_Greene_2026-09-24.csv` | Side-by-side comparison (drives Table 3 and Figure 7) |
| `outputs/table.Greene_heterogeneity_2026-09-24.csv` | Cochran's Q, df, Q-P, I² per stratum |
| `outputs/table.Greene_egger_intercept_2026-09-24.csv` | MR-Egger intercept per stratum |
| `outputs/table.Greene_MRPRESSO_v5_2026-09-26.csv` | **MR-PRESSO global test per stratum, `seed = 20260926`** — no outlier in any stratum |
| `outputs/table.Greene_steiger_v5_2026-09-26.csv` | Steiger directionality test (direction correct, P = 7.9e-171) |
| `outputs/table.Greene_LOO_IL6R_2026-09-24.csv` | Within-IL6R-cluster leave-one-out |
| `outputs/table.Greene_LOO_HLA_2026-09-24.csv` | Within-HLA-cluster leave-one-out |
| `outputs/table.Greene_coverage_audit.csv` | Instrument coverage (22/23 by rsID, 23/23 by position) |
| `outputs/table.harmonised_Greene.csv` | Harmonised SNP-level data for the replication dataset |
| `outputs/table.Sakaue_MR_all_methods_2026-09-24.csv` | **Authoritative** discovery output (supersedes `../outputs/table.MRresult.csv`) |
| `outputs/table.power_v5_2026-09-26.csv` | Approximate power, per-dataset R², corrected control count |
| `run_log_replication.txt` | Console log of the 2026-09-24 replication run (original, pre-fix) |
| `run_log_presso_steiger_v5_2026-09-26.txt` | Console log of the 2026-09-26 seeded MR-PRESSO + Steiger + power run |

### MR-PRESSO is a permutation test — the seed matters

With `NbDistribution = 1000`, the global-test P value is permutation-based and varies between
unseeded runs. After adding `seed = 20260926`:

| Analysis | RSSobs (deterministic) | unseeded | **seeded** |
|---|---|---|---|
| All 23 SNPs | 6.4090 | 1.000 | 1.000 |
| All minus rs2228145 | 6.5204 | 1.000 | 1.000 |
| HLA locus only | 0.5229 | 1.000 | 1.000 |
| IL6R locus only | 6.0611 | 0.931 | **0.934** |
| IL6R minus rs2228145 | 6.3203 | 0.864 | **0.901** |

RSSobs is unchanged and so are all conclusions: no outlier variant and no evidence of
directional pleiotropy in any stratum. The manuscript reports the seeded values.

## Environment

R 4.4.3; TwoSampleMR 0.7.9; MRPRESSO 1.0; data.table 1.18.6.1; ggplot2 4.0.3.

## Reproduction notes

- GRCh37/hg19, 1-based positions; exposure effect estimates per SD.
- `rs1129737` is absent under that rsID in the replication dataset and was matched by
  chromosome:position; alleles and effect direction were verified to be consistent.
- `rs5004279` was removed during harmonisation (palindromic with intermediate allele frequency,
  EAF 0.478).
- The bootstrap-based estimators (weighted median, weighted mode, simple mode) are **not seeded**
  in TwoSampleMR and differ slightly between runs; point estimates (β, OR) are stable. The
  archived `table.*_MR_all_methods_2026-09-24.csv` files are the ones the manuscript reports.
- The replication outcome dataset overlaps the discovery dataset in three contributing biobanks
  (UK Biobank, FinnGen, BioBank Japan). This is disclosed in the manuscript.
