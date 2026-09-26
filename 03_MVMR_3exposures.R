# ================================================================================
# MVMR Analysis (Final)
# Exposures: IL-6, IL-1RA, TNFR1 → Outcome: Keloids
# ================================================================================
cat("========================================================\n")
cat("MVMR: IL-6 + IL-1RA + TNFR1 → Keloids\n")
cat("========================================================\n\n")

DATA_DIR <- "F:/MR/Interleukin-6 levels and keloids/数据/Multivariable MR 数据"
setwd(DATA_DIR)

suppressPackageStartupMessages({ library(dplyr) })

# ---- Read prepared data ----
dat <- read.csv("MVMR_input_data.csv", stringsAsFactors = FALSE)
cat(paste0("Total SNPs: ", nrow(dat), "\n\n"))

# ================================================================================
# 1. Allele Harmonization
# ================================================================================
cat("--- Allele harmonization ---\n")
harmonized <- dat
harmonized$harmonized <- TRUE
harmonized$palindromic <- FALSE

for (i in 1:nrow(harmonized)) {
  ea_il6 <- toupper(harmonized$ea_IL6[i])
  oa_il6 <- toupper(harmonized$oa_IL6[i])
  
  # Palindromic check
  if (paste(sort(c(ea_il6, oa_il6)), collapse = "/") %in% c("A/T", "C/G")) {
    harmonized$palindromic[i] <- TRUE
  }
  
  # Keloid
  ea_kel <- toupper(harmonized$ea_keloid[i])
  oa_kel <- toupper(harmonized$oa_keloid[i])
  if (ea_il6 == oa_kel && oa_il6 == ea_kel) {
    harmonized$beta_keloid[i] <- -harmonized$beta_keloid[i]
  } else if (ea_il6 != ea_kel) {
    harmonized$harmonized[i] <- FALSE
  }
  
  # IL-1RA
  ea_i1 <- toupper(harmonized$ea_IL1RA[i]); oa_i1 <- toupper(harmonized$oa_IL1RA[i])
  if (ea_il6 == oa_i1 && oa_il6 == ea_i1) {
    harmonized$beta_IL1RA[i] <- -harmonized$beta_IL1RA[i]
  } else if (ea_il6 != ea_i1) {
    harmonized$harmonized[i] <- FALSE
  }
  
  # TNFR1
  ea_t1 <- toupper(harmonized$ea_TNFR1[i]); oa_t1 <- toupper(harmonized$oa_TNFR1[i])
  if (ea_il6 == oa_t1 && oa_il6 == ea_t1) {
    harmonized$beta_TNFR1[i] <- -harmonized$beta_TNFR1[i]
  } else if (ea_il6 != ea_t1) {
    harmonized$harmonized[i] <- FALSE
  }
}

dat_h <- harmonized[harmonized$harmonized == TRUE, ]
cat(paste0("Harmonized SNPs: ", nrow(dat_h), "/", nrow(dat), "\n"))
cat(paste0("Palindromic SNPs: ", sum(harmonized$palindromic), "\n\n"))

# ================================================================================
# 2. MVMR Analysis
# ================================================================================
cat("--- Running MVMR ---\n")

beta_exp <- as.matrix(dat_h[, c("beta_IL6", "beta_IL1RA", "beta_TNFR1")])
beta_out <- dat_h$beta_keloid
se_out   <- dat_h$se_keloid
weights  <- 1 / se_out^2

run_mvmr <- function(X, Y, w, label) {
  XtWX <- t(X) %*% diag(w) %*% X
  XtWY <- t(X) %*% diag(w) %*% Y
  
  if (rcond(XtWX) < 1e-14) {
    cat(paste0("  Singular matrix (rcond=", rcond(XtWX), "), trying solve(qr)...\n"))
    qr_sol <- qr.solve(XtWX, XtWY)
    beta_hat <- qr_sol
  } else {
    beta_hat <- solve(XtWX, XtWY)
  }
  
  resid <- Y - X %*% beta_hat
  sigma2 <- sum(w * resid^2) / (nrow(X) - ncol(X))
  var_beta <- sigma2 * solve(XtWX)
  se_beta <- sqrt(diag(var_beta))
  z <- beta_hat / se_beta
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  nms <- colnames(X)
  cat(paste0("\n", label, ":\n"))
  for (j in 1:length(nms)) {
    or <- exp(beta_hat[j])
    or_l <- exp(beta_hat[j] - 1.96 * se_beta[j])
    or_u <- exp(beta_hat[j] + 1.96 * se_beta[j])
    p_fmt <- ifelse(p[j] < 0.001, format(p[j], scientific=TRUE, digits=3), sprintf("%.4f", p[j]))
    cat(sprintf("  %-10s β=%.4f SE=%.4f OR=%.4f (%.4f-%.4f) Z=%.3f P=%s\n",
                nms[j], beta_hat[j], se_beta[j], or, or_l, or_u, z[j], p_fmt))
  }
  
  return(data.frame(Exposure = nms, Beta = as.numeric(beta_hat), 
                    SE = as.numeric(se_beta),
                    OR = as.numeric(exp(beta_hat)),
                    OR_lower = as.numeric(exp(beta_hat - 1.96*se_beta)),
                    OR_upper = as.numeric(exp(beta_hat + 1.96*se_beta)),
                    Z = as.numeric(z), P = as.numeric(p), 
                    stringsAsFactors = FALSE))
}

# Model 1: With intercept
X1 <- cbind(1, beta_exp)
colnames(X1) <- c("Intercept", "IL6", "IL1RA", "TNFR1")
r1 <- run_mvmr(X1, beta_out, weights, "MVMR-IVW (with intercept)")

# Model 2: Without intercept  
X2 <- beta_exp
colnames(X2) <- c("IL6", "IL1RA", "TNFR1")
r2 <- run_mvmr(X2, beta_out, weights, "MVMR-IVW (without intercept)")

# ================================================================================
# 3. Diagnostics
# ================================================================================
cat("\n\n--- Diagnostics ---\n")

# Conditional F-statistics
cat("\nConditional F-statistics:\n")
for (i in 1:ncol(beta_exp)) {
  nm <- colnames(beta_exp)[i]
  df_full <- as.data.frame(beta_exp); colnames(df_full) <- colnames(beta_exp)
  df_red  <- as.data.frame(beta_exp[, -i, drop = FALSE])
  
  lm_f <- lm(beta_out ~ 0 + ., data = df_full, weights = weights)
  lm_r <- lm(beta_out ~ 0 + ., data = df_red, weights = weights)
  a <- anova(lm_r, lm_f, test = "F")
  cat(sprintf("  %s: F=%.3f, P=%.4f\n", nm, a$F[2], a$`Pr(>F)`[2]))
}

# Q-statistic
resid_q <- beta_out - X2 %*% solve(t(X2) %*% diag(weights) %*% X2) %*% t(X2) %*% diag(weights) %*% beta_out
Q <- sum(weights * resid_q^2)
df_q <- nrow(X2) - ncol(X2)
cat(sprintf("\nCochran's Q: %.3f, df=%d, P=%.4f\n", Q, df_q, pchisq(Q, df_q, lower.tail=FALSE)))

# Exposure correlation
cat("\nCorrelation of SNP-exposure effects:\n")
print(round(cor(beta_exp), 3))

# ================================================================================
# 4. Save
# ================================================================================
cat("\n\n--- Saving results ---\n")
r1$Model <- "With intercept"
r2$Model <- "Without intercept"
all_results <- rbind(r1, r2)
write.csv(all_results, "table.MVMR_results.csv", row.names = FALSE)
write.csv(dat_h, "table.MVMR_harmonised_data.csv", row.names = FALSE)
cat("  ✓ table.MVMR_results.csv\n")
cat("  ✓ table.MVMR_harmonised_data.csv\n")

# ================================================================================
# Summary
# ================================================================================
cat("\n\n========================================================\n")
cat("Summary\n")
cat("========================================================\n")
cat(paste0("SNPs: ", nrow(dat_h), "\n"))
cat(paste0("Exposures: IL-6, IL-1RA, TNFR1\n"))
cat(paste0("Outcome: Keloids\n\n"))
cat("Key finding:\n")
il6 <- r2[r2$Exposure == "IL6", ]
cat(sprintf("IL-6 (adj. IL-1RA & TNFR1): OR=%.4f (%.4f-%.4f), P=%s\n",
            il6$OR, il6$OR_lower, il6$OR_upper,
            ifelse(il6$P<0.001,format(il6$P,scientific=TRUE,digits=3),sprintf("%.4f",il6$P))))
cat(paste0("\nFinished: ", Sys.time(), "\n"))
