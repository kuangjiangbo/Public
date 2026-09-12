# MVMR Diagnostics
setwd("F:/MR/Interleukin-6 levels and keloids/数据/Multivariable MR 数据")
dat_h <- read.csv("table.MVMR_harmonised_data.csv", stringsAsFactors = FALSE)

cat("=== MVMR Diagnostics ===\n\n")
cat(paste0("N SNPs: ", nrow(dat_h), "\n\n"))

beta_exp <- as.matrix(dat_h[, c("beta_IL6", "beta_IL1RA", "beta_TNFR1")])
beta_out <- dat_h$beta_keloid
se_out   <- dat_h$se_keloid
weights  <- 1 / se_out^2

# Conditional F-statistics
cat("Conditional F-statistics:\n")
for (i in 1:ncol(beta_exp)) {
  exp_name <- colnames(beta_exp)[i]
  
  # Full model
  X_full <- as.data.frame(beta_exp)
  colnames(X_full) <- colnames(beta_exp)
  
  # Reduced model (without exposure i)
  X_red <- as.data.frame(beta_exp[, -i, drop = FALSE])
  
  lm_full <- lm(beta_out ~ 0 + ., data = X_full, weights = weights)
  lm_red  <- lm(beta_out ~ 0 + ., data = X_red, weights = weights)
  
  anova_res <- anova(lm_red, lm_full, test = "F")
  cat(sprintf("  %s: conditional F = %.3f, P = %.4f\n", 
              exp_name, anova_res$F[2], anova_res$`Pr(>F)`[2]))
}

# Q-statistic for heterogeneity
cat("\nQ-statistic for heterogeneity:\n")
resid <- beta_out - beta_exp %*% solve(t(beta_exp) %*% diag(weights) %*% beta_exp) %*% t(beta_exp) %*% diag(weights) %*% beta_out
Q <- sum(weights * resid^2)
df <- nrow(beta_exp) - ncol(beta_exp)
Q_pval <- pchisq(Q, df, lower.tail = FALSE)
cat(sprintf("  Q = %.3f, df = %d, P = %.4f\n", Q, df, Q_pval))

# Correlation between exposures
cat("\nCorrelation of SNP-exposure effects:\n")
cor_mat <- cor(beta_exp)
print(round(cor_mat, 3))

cat("\nDone!\n")
