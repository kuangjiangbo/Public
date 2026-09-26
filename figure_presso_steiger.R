## v5 补跑（修正版）：MR-PRESSO 的参数名在 MRPRESSO 1.0 中不存在 => 去掉 OUTCOME
## 同时用 get_r_from_lor() 为二分类结局计算 r.outcome，再做 Steiger
BASE <- "F:/MR/Interleukin-6 levels and keloids"
OUT  <- file.path(BASE, "数据/Greene2025")
HARM <- file.path(OUT, "table.harmonised_Greene.csv")

suppressMessages({library(data.table); library(TwoSampleMR); library(MRPRESSO)})
cat("MRPRESSO args:", paste(names(formals(MRPRESSO::mr_presso)), collapse = ", "), "\n\n")

harm <- as.data.frame(fread(HARM))
LEAD <- "rs2228145"
sets <- list(
  "All 23 SNPs"          = harm,
  "IL6R locus only"      = harm[harm$locus == "IL6R", ],
  "HLA locus only"       = harm[harm$locus == "HLA", ],
  "All minus rs2228145"  = harm[harm$SNP != LEAD, ],
  "IL6R minus rs2228145" = harm[harm$locus == "IL6R" & harm$SNP != LEAD, ])

res <- list()
for (k in names(sets)) {
  d <- sets[[k]]
  cat("=== ", k, "  nSNP=", nrow(d), " ===\n", sep = "")
  r <- tryCatch(
    MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                        SdOutcome = "se.outcome", SdExposure = "se.exposure",
                        data = d, NbDistribution = 1000, SignifThreshold = 0.05),
    error = function(e) { cat("  !!! ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) {
    res[[k]] <- data.frame(Analysis = k, nSNP = nrow(d), RSSobs = NA_real_,
                           GlobalTest_P = NA_real_, N_outliers = NA_integer_,
                           note = "MR-PRESSO failed"); next }
  gt <- r$`MR-PRESSO results`$`Global Test`
  ot <- r$`MR-PRESSO results`$`Outlier Test`
  no <- if (is.null(ot)) NA_integer_ else sum(ot$Pvalue < 0.05, na.rm = TRUE)
  cat("  RSSobs =", gt$RSSobs, "| Global P =", gt$Pvalue, "| outliers(P<0.05) =", no, "\n")
  print(ot)
  mc <- r$`MR-PRESSO results`$`Main MR results`
  if (!is.null(mc)) { cat("  Main MR results:\n"); print(mc) }
  res[[k]] <- data.frame(Analysis = k, nSNP = nrow(d), RSSobs = gt$RSSobs,
                         GlobalTest_P = gt$Pvalue, N_outliers = no, note = NA_character_)
}
out <- rbindlist(res)
cat("\n-- MR-PRESSO 汇总（复制数据集, NbDistribution=1000）--\n"); print(out, row.names = FALSE)
fwrite(out, file.path(OUT, "table.Greene_MRPRESSO_v5_2026-09-26.csv"))
cat("  已写: table.Greene_MRPRESSO_v5_2026-09-26.csv\n\n")

## ---- Steiger：先用 get_r_from_lor 为二分类结局算 r.outcome ----
cat("== Steiger 方向性检验（Greene 2025, 全 22 SNP, 二分类结局口径）==\n")
d <- harm
ncase <- 7837; nctrl <- 1585173; prev <- ncase / (ncase + nctrl)
d$r.outcome <- get_r_from_lor(lor = d$beta.outcome, af = d$eaf.outcome,
                              ncase = ncase, ncontrol = nctrl, prevalence = prev)
d$samplesize.exposure <- 21758
d$units.exposure <- "SD"; d$units.outcome <- "log odds"
d$ncase.outcome <- ncase; d$ncontrol.outcome <- nctrl; d$prevalence.outcome <- prev
d$id.exposure <- "IL6"; d$id.outcome <- "GCST90652487"
st <- tryCatch(directionality_test(d),
  error = function(e) { cat("  !!! Steiger ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(st)) {
  print(st, row.names = FALSE)
  fwrite(as.data.table(st), file.path(OUT, "table.Greene_steiger_v5_2026-09-26.csv"))
  cat("  已写: table.Greene_steiger_v5_2026-09-26.csv\n")
}
cat("\n完成。\n")
