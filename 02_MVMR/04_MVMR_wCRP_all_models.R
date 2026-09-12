# ================================================================================
# MVMR Analysis with CRP
# Models: (1) IL-6 + CRP → Keloid  (2) IL-6 + IL-1RA + TNFR1 + CRP → Keloid
# ================================================================================
cat("========================================================\n")
cat("MVMR Analysis with CRP\n")
cat("========================================================\n\n")

DATA_DIR <- "F:/MR/Interleukin-6 levels and keloids/数据/Multivariable MR 数据"
setwd(DATA_DIR)

run_mvmr <- function(csv_file, exp_names, model_label) {
  cat(paste0("\n\n========== ", model_label, " ==========\n"))
  
  dat <- read.csv(csv_file, stringsAsFactors = FALSE)
  cat(paste0("SNPs: ", nrow(dat), "\n"))
  
  # Harmonize alleles (align to IL-6 EA)
  for (i in 1:nrow(dat)) {
    ea_ref <- toupper(dat$ea_IL6[i])
    oa_ref <- toupper(dat$oa_IL6[i])
    
    for (exp in exp_names) {
      if (exp == "IL6") next
      ea_x <- toupper(dat[[paste0("ea_", exp)]][i])
      oa_x <- toupper(dat[[paste0("oa_", exp)]][i])
      
      if (ea_ref == oa_x && oa_ref == ea_x) {
        dat[[paste0("beta_", exp)]][i] <- -dat[[paste0("beta_", exp)]][i]
      }
    }
    
    # Outcome
    ea_k <- toupper(dat$ea_keloid[i])
    oa_k <- toupper(dat$oa_keloid[i])
    if (ea_ref == oa_k && oa_ref == ea_k) {
      dat$beta_keloid[i] <- -dat$beta_keloid[i]
    }
  }
  
  # Prepare matrices
  beta_exp <- as.matrix(dat[, paste0("beta_", exp_names), drop = FALSE])
  colnames(beta_exp) <- exp_names
  beta_out <- dat$beta_keloid
  se_out   <- dat$se_keloid
  weights  <- 1 / se_out^2
  
  # MVMR with intercept
  X <- cbind(1, beta_exp)
  colnames(X) <- c("Intercept", exp_names)
  
  XtWX <- t(X) %*% diag(weights) %*% X
  XtWY <- t(X) %*% diag(weights) %*% beta_out
  
  if (rcond(XtWX) < 1e-14) {
    beta_hat <- qr.solve(XtWX, XtWY)
  } else {
    beta_hat <- solve(XtWX, XtWY)
  }
  
  resid <- beta_out - X %*% beta_hat
  sigma2 <- sum(weights * resid^2) / (nrow(X) - ncol(X))
  se <- sqrt(diag(sigma2 * solve(XtWX)))
  z <- as.numeric(beta_hat / se)
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  # Results table
  res <- data.frame(
    Model = model_label,
    Exposure = colnames(X),
    N_SNPs = nrow(dat),
    Beta = round(as.numeric(beta_hat), 5),
    SE = round(as.numeric(se), 5),
    OR = round(exp(as.numeric(beta_hat)), 4),
    OR_lower = round(exp(as.numeric(beta_hat) - 1.96*as.numeric(se)), 4),
    OR_upper = round(exp(as.numeric(beta_hat) + 1.96*as.numeric(se)), 4),
    Z = round(as.numeric(z), 3),
    P = as.numeric(p),
    stringsAsFactors = FALSE
  )
  
  cat("\nResults (with intercept):\n")
  for (j in 1:nrow(res)) {
    p_str <- ifelse(res$P[j] < 0.001, format(res$P[j], scientific=TRUE, digits=3), sprintf("%.4f", res$P[j]))
    cat(sprintf("  %-12s OR=%.4f (%.4f-%.4f), P=%s\n",
                res$Exposure[j], res$OR[j], res$OR_lower[j], res$OR_upper[j], p_str))
  }
  
  # Also without intercept
  X2 <- beta_exp
  colnames(X2) <- exp_names
  
  XtWX2 <- t(X2) %*% diag(weights) %*% X2
  XtWY2 <- t(X2) %*% diag(weights) %*% beta_out
  
  if (rcond(XtWX2) < 1e-14) {
    beta_hat2 <- qr.solve(XtWX2, XtWY2)
  } else {
    beta_hat2 <- solve(XtWX2, XtWY2)
  }
  
  resid2 <- beta_out - X2 %*% beta_hat2
  sigma2_2 <- sum(weights * resid2^2) / (nrow(X2) - ncol(X2))
  se2 <- sqrt(diag(sigma2_2 * solve(XtWX2)))
  z2 <- as.numeric(beta_hat2 / se2)
  p2 <- 2 * pnorm(abs(z2), lower.tail = FALSE)
  
  res2 <- data.frame(
    Model = paste0(model_label, " (no intercept)"),
    Exposure = exp_names,
    N_SNPs = nrow(dat),
    Beta = round(as.numeric(beta_hat2), 5),
    SE = round(as.numeric(se2), 5),
    OR = round(exp(as.numeric(beta_hat2)), 4),
    OR_lower = round(exp(as.numeric(beta_hat2) - 1.96*as.numeric(se2)), 4),
    OR_upper = round(exp(as.numeric(beta_hat2) + 1.96*as.numeric(se2)), 4),
    Z = round(as.numeric(z2), 3),
    P = as.numeric(p2),
    stringsAsFactors = FALSE
  )
  
  cat("\nResults (no intercept):\n")
  for (j in 1:nrow(res2)) {
    p_str <- ifelse(res2$P[j] < 0.001, format(res2$P[j], scientific=TRUE, digits=3), sprintf("%.4f", res2$P[j]))
    cat(sprintf("  %-12s OR=%.4f (%.4f-%.4f), P=%s\n",
                res2$Exposure[j], res2$OR[j], res2$OR_lower[j], res2$OR_upper[j], p_str))
  }
  
  # Conditional F
  cat("\nConditional F-statistics:\n")
  for (exp in exp_names) {
    cols_all <- exp_names
    cols_red <- setdiff(exp_names, exp)
    
    df_f <- as.data.frame(beta_exp[, cols_all, drop = FALSE])
    df_r <- as.data.frame(beta_exp[, cols_red, drop = FALSE])
    colnames(df_f) <- cols_all
    colnames(df_r) <- cols_red
    
    lm_f <- lm(beta_out ~ 0 + ., data = df_f, weights = weights)
    lm_r <- lm(beta_out ~ 0 + ., data = df_r, weights = weights)
    a <- anova(lm_r, lm_f, test = "F")
    cat(sprintf("  %s: F=%.3f, P=%.4f\n", exp, a$F[2], a$`Pr(>F)`[2]))
  }
  
  # Correlation
  cat("\nExposure effect correlations:\n")
  print(round(cor(beta_exp), 3))
  
  # Q statistic
  resid_q <- beta_out - X2 %*% beta_hat2
  Q <- sum(weights * resid_q^2)
  df_q <- nrow(X2) - ncol(X2)
  cat(sprintf("\nCochran's Q: %.3f, df=%d, P=%.4f\n", Q, df_q, pchisq(Q, df_q, lower.tail=FALSE)))
  
  return(rbind(res, res2))
}

# Model 1: IL-6 + CRP
r1 <- run_mvmr("MVMR_input_IL6_CRP.csv", c("IL6", "CRP"), "IL6 + CRP")

# Model 2: All 4 exposures
r2 <- run_mvmr("MVMR_input_wCRP.csv", c("IL6", "IL1RA", "TNFR1", "CRP"), "IL6 + IL1RA + TNFR1 + CRP")

# Model 3: Compare with previous results (3 exposures without CRP)
r3 <- run_mvmr("MVMR_input_data.csv", c("IL6", "IL1RA", "TNFR1"), "IL6 + IL1RA + TNFR1")

# Save all
all_res <- rbind(r1, r2, r3)
write.csv(all_res, "table.MVMR_all_results.csv", row.names = FALSE)

cat("\n\n========================================================\n")
cat("SUMMARY: Key IL-6 Results Across Models\n")
cat("========================================================\n\n")

il6_rows <- all_res[all_res$Exposure == "IL6" & !grepl("no intercept", all_res$Model), ]
for (i in 1:nrow(il6_rows)) {
  p_str <- ifelse(il6_rows$P[i] < 0.001, format(il6_rows$P[i], scientific=TRUE, digits=3), sprintf("%.4f", il6_rows$P[i]))
  cat(sprintf("  %-45s OR=%.4f (%.4f-%.4f), P=%s\n",
              il6_rows$Model[i], il6_rows$OR[i], il6_rows$OR_lower[i], il6_rows$OR_upper[i], p_str))
}

cat(paste0("\nSaved: table.MVMR_all_results.csv\n"))
cat(paste0("Finished: ", Sys.time(), "\n"))
