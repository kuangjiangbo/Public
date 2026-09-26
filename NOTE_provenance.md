# Outputs — provenance and known inconsistencies

## 1. `table.MRresult.csv` is from an earlier run

`table.MRresult.csv` was written by an earlier execution of
`01_univariate_MR/univariate_mr_pipeline.R`. The three bootstrap-based estimators (weighted
median, weighted mode, simple mode) are stochastic — the bootstrap step is **not seeded** in
TwoSampleMR — so their SEs and P-values differ slightly between runs.

The values reported in the manuscript are those of the **later, authoritative run**, archived as
`../05_replication/outputs/table.Sakaue_MR_all_methods_2026-09-24.csv`:

| Estimator | `table.MRresult.csv` (earlier) | manuscript / authoritative file |
|---|---|---|
| Weighted median | β 0.2335, SE 0.12332, P 0.0583 | β 0.2335, SE 0.12117, P 0.0540 |
| Weighted mode | β 0.3488, SE 0.19324, P 0.0848 | β 0.3488, SE 0.18604, P 0.0742 |
| Simple mode | β 0.3058, SE 0.21799, P 0.1746 | β 0.3058, SE 0.22701, P 0.1917 |

Point estimates (β) and ORs are identical across runs; only the bootstrap-derived SEs differ.
No conclusion in the manuscript depends on the third significant figure of these estimates.

**Recommended:** state a fixed bootstrap seed and regenerate once, or keep both output sets and
retain this note. Do not delete this note while both files are present.

## 2. The power table was superseded on 2026-09-26

The earlier `table.power_2026-09-24.csv` contained two errors and has been **removed** from this
repository:

- its control count for the replication dataset was back-derived from the Greene 2025 SSF `N`
  column (1,593,010), which contradicts the paper's own Table 1 (7,837 cases / 1,593,009 controls);
- both datasets shared the 23-SNP cumulative R² (0.041407) although the replication analysis used
  only 22 SNPs (rs5004279 dropped as palindromic).

The current file is `../05_replication/outputs/table.power_v5_2026-09-26.csv`, which uses
per-dataset R² (0.041407 discovery / 0.039973 replication) and the published control count.
Net effect: the replication *predicted* SE changes from 0.0278 to 0.0283. All power percentages
quoted in the manuscript are unchanged (discovery 48% / 99% / 52%; replication 97% / 53%).

## 3. MR-PRESSO results must be read with the seed

The MR-PRESSO global test is permutation-based (`NbDistribution = 1000`). Results are only
reproducible with an explicit seed; `05_replication/IL6_Keloid_Greene2025_replication.R` uses
`seed = 20260926`. See `05_replication/README.md` for the seeded-vs-unseeded comparison.
