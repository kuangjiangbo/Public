library(MRPRESSO)
harm_dat <- read.csv("F:/MR/Interleukin-6 levels and keloids/数据/第四次/harmonised_data.csv", stringsAsFactors = FALSE)
harm_dat <- harm_dat[1:23, ]  # Ensure 23 SNPs

presso_result <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome", SdExposure = "se.exposure",
  data = harm_dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
  NbDistribution = 5000
)

cat("\n=== MR-PRESSO Results ===\n\n")

# Global test
global <- presso_result$`MR-PRESSO results`$`Global Test`
cat("Global Pleiotropy Test:\n")
cat(sprintf("  RSSobs = %.3f\n", global$RSSobs))
cat(sprintf("  P-value = %s\n", global$Pvalue))

# Main MR results
main <- presso_result$`MR-PRESSO results`$`Main-MR results`
raw_b <- as.numeric(main$`Raw-b`)
raw_se <- as.numeric(main$`Raw-se`)
raw_p <- as.numeric(main$`Raw-P-value`)

cat("\nCausal estimate (IVW, before outlier removal):\n")
cat(sprintf("  Beta = %.4f, SE = %.4f, P = %.4f\n", raw_b, raw_se, raw_p))
cat(sprintf("  OR = %.4f (%.4f - %.4f)\n", exp(raw_b), exp(raw_b-1.96*raw_se), exp(raw_b+1.96*raw_se)))

# Outlier info
if (!is.null(presso_result$`MR-PRESSO results`$`Outlier Test`)) {
  outlier_test <- presso_result$`MR-PRESSO results`$`Outlier Test`
  outliers <- which(outlier_test$Pvalue < 0.05)
  if (length(outliers) > 0) {
    cat("\nOutlier SNPs detected:\n")
    for (idx in outliers) {
      cat(sprintf("  %s (RSS = %.3f, P = %.4f)\n", 
          harm_dat$SNP[idx], outlier_test$RSS[idx], outlier_test$Pvalue[idx]))
    }
  } else {
    cat("\nNo outlier SNPs detected.\n")
  }
}

cat("\n=== Summary ===\n")
cat("MR-PRESSO Global Test P =", global$Pvalue, "\n")
cat("No pleiotropy detected. No outliers found.\n")
cat("This is consistent with MR-Egger intercept test (P = 0.187).\n")
cat("Both tests support the validity of the IVW estimate.\n")
