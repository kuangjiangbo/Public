# ================================================================================
# MR-PRESSO Analysis: IL-6 → Keloids
# ================================================================================
# MR-PRESSO detects and corrects for outlier SNPs with horizontal pleiotropy
# Recommended by reviewer as supplement/alternative to MR-Egger
# ================================================================================

cat("========================================================\n")
cat("MR-PRESSO Analysis\n")
cat("IL-6 levels → Keloids\n")
cat("========================================================\n\n")

setwd("F:/MR/Interleukin-6 levels and keloids/数据/第四次")

suppressPackageStartupMessages({
  library(MRPRESSO)
  library(dplyr)
})

# ---- Load harmonised data ----
# The harmonised data from the original analysis should contain
# beta.exposure, se.exposure, beta.outcome, se.outcome for the 23 SNPs
harm_dat <- read.csv("harmonised_data.csv", stringsAsFactors = FALSE)
cat(paste0("Loaded harmonised data: ", nrow(harm_dat), " rows\n"))

# Check which columns are available
cat("Columns:", paste(colnames(harm_dat), collapse = ", "), "\n\n")

# The original harmonised data from the 第四次 analysis
# Expected columns: SNP, beta.exposure, se.exposure, beta.outcome, se.outcome, etc.
# Let's check the structure

required_cols <- c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
missing <- required_cols[!required_cols %in% colnames(harm_dat)]
if (length(missing) > 0) {
  cat("Missing columns:", paste(missing, collapse = ", "), "\n")
  cat("Trying alternative column names...\n")
  
  # Try alternative names
  col_aliases <- list(
    beta.exposure = c("beta.exposure", "beta_exp", "beta_IL6"),
    se.exposure = c("se.exposure", "se_exp", "se_IL6"),
    beta.outcome = c("beta.outcome", "beta_out", "beta_keloid"),
    se.outcome = c("se.outcome", "se_out", "se_keloid")
  )
  
  for (req in names(col_aliases)) {
    if (!req %in% colnames(harm_dat)) {
      found <- col_aliases[[req]][col_aliases[[req]] %in% colnames(harm_dat)]
      if (length(found) > 0) {
        colnames(harm_dat)[colnames(harm_dat) == found[1]] <- req
        cat(paste0("  Remapped ", found[1], " -> ", req, "\n"))
      }
    }
  }
}

# Check again
missing <- required_cols[!required_cols %in% colnames(harm_dat)]
if (length(missing) > 0) {
  cat("\nStill missing columns:", paste(missing, collapse = ", "), "\n")
  cat("Trying to build from MVMR harmonised data...\n")
  
  # Try the MVMR harmonised data
  mvmr_harm <- read.csv(
    "F:/MR/Interleukin-6 levels and keloids/数据/Multivariable MR 数据/table.MVMR_harmonised_data.csv",
    stringsAsFactors = TRUE
  )
  
  # Extract IL-6 exposure and keloid outcome
  test_dat <- mvmr_harm[, c("SNP", "beta_IL6", "se_IL6", "beta_keloid", "se_keloid")]
  colnames(test_dat) <- c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
  harm_dat <- test_dat
  cat(paste0("Using MVMR harmonised data: ", nrow(harm_dat), " SNPs\n"))
}

# Run MR-PRESSO
cat("\n--- Running MR-PRESSO ---\n")
cat(paste0("SNPs: ", nrow(harm_dat), "\n"))
cat(paste0("Exposure: beta.exposure, SE: se.exposure\n"))
cat(paste0("Outcome: beta.outcome, SE: se.outcome\n\n"))

# MR-PRESSO with 5000 permutations
set.seed(20260618)
presso_result <- tryCatch({
  mr_presso(
    BetaOutcome = "beta.outcome",
    BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome",
    SdExposure = "se.exposure",
    data = harm_dat,
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE,
    NbDistribution = 5000
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

# ---- Print Results ----
if (!is.null(presso_result)) {
  cat("\n========================================================\n")
  cat("MR-PRESSO Results\n")
  cat("========================================================\n\n")
  
  # Global pleiotropy test
  cat("--- Global Pleiotropy Test ---\n")
  global <- presso_result$`MR-PRESSO results`$`Global Test`
  cat(paste0("  RSSobs = ", round(global$RSSobs, 3), "\n"))
  cat(paste0("  P-value = ", global$Pvalue, "\n"))
  
  if (global$Pvalue < 0.05) {
    cat("  -> Significant pleiotropy detected!\n")
  } else {
    cat("  -> No significant pleiotropy detected.\n")
  }
  
  # Outlier test
  cat("\n--- Outlier Test ---\n")
  outlier <- presso_result$`MR-PRESSO results`$`Distortion Test`
  
  if (!is.null(outlier) && length(outlier$`Outliers Indices`) > 0) {
    outliers_idx <- outlier$`Outliers Indices`
    outliers_snp <- harm_dat$SNP[outliers_idx]
    cat(paste0("  Outliers detected: ", length(outliers_idx), "\n"))
    cat(paste0("  Outlier SNPs: ", paste(outliers_snp, collapse = ", "), "\n"))
    cat(paste0("  Distortion P-value: ", outlier$`Pvalue`, "\n"))
    
    if (outlier$`Pvalue` < 0.05) {
      cat("  -> Removing outliers significantly changes the estimate!\n")
    } else {
      cat("  -> Removing outliers does NOT significantly change the estimate.\n")
    }
  } else {
    cat("  No outlier SNPs detected.\n")
  }
  
  # Causal estimates before and after
  cat("\n--- Causal Estimates ---\n")
  main <- presso_result$`MR-PRESSO results`$`Main-MR results`
  
  # Before outlier removal (Raw)
  raw_b <- main$`Raw-b`
  raw_se <- main$`Raw-se`
  raw_p <- main$`Raw-P-value`
  raw_or <- exp(raw_b)
  raw_or_l <- exp(raw_b - 1.96 * raw_se)
  raw_or_u <- exp(raw_b + 1.96 * raw_se)
  
  cat("\nBefore outlier removal (IVW):\n")
  cat(sprintf("  OR = %.4f (%.4f - %.4f)\n", raw_or, raw_or_l, raw_or_u))
  cat(sprintf("  Beta = %.4f, SE = %.4f, P = %.4f\n", raw_b, raw_se, raw_p))
  
  # After outlier removal (if applicable)
  if (!is.null(main$`b`)) {
    corr_b <- main$`b`
    corr_se <- main$`se`
    corr_p <- main$`P-value`
    corr_or <- exp(corr_b)
    corr_or_l <- exp(corr_b - 1.96 * corr_se)
    corr_or_u <- exp(corr_b + 1.96 * corr_se)
    
    cat("\nAfter outlier removal:\n")
    cat(sprintf("  OR = %.4f (%.4f - %.4f)\n", corr_or, corr_or_l, corr_or_u))
    cat(sprintf("  Beta = %.4f, SE = %.4f, P = %.4f\n", corr_b, corr_se, corr_p))
  }
  
  # Comparison with original IVW
  cat("\n\n--- Comparison with Original IVW ---\n")
  cat("Original IVW: OR = 1.214 (1.016-1.452), P = 0.033\n")
  cat(sprintf("MR-PRESSO Raw: OR = %.4f (%.4f-%.4f), P = %.4f\n", 
              raw_or, raw_or_l, raw_or_u, raw_p))
  
  cat("\n--- Interpretation for manuscript ---\n")
  if (global$Pvalue > 0.05) {
    cat("MR-PRESSO did not detect significant pleiotropy, consistent with\n")
    cat("the MR-Egger intercept test. This strengthens confidence in the IVW results.\n")
    cat("Recommendation: Report MR-PRESSO as a supplementary sensitivity analysis.\n")
  } else {
    cat("MR-PRESSO detected significant pleiotropy with outlier SNPs.\n")
    if (!is.null(outlier) && length(outlier$`Outliers Indices`) > 0) {
      cat("Outlier-corrected estimate should be compared with IVW.\n")
      cat("Discuss outliers in the context of IL6R/HLA region instruments.\n")
    }
  }
  
} else {
  cat("\nMR-PRESSO failed to run. Please check data format.\n")
}

cat(paste0("\nAnalysis completed: ", Sys.time(), "\n"))
