# ================================================================================
# 孟德尔随机化分析脚本 (V5 - 整理版)
# 暴露：IL-6 (ebi-a-GCST90012005)
# 结局：Keloids (ebi-a-GCST90018874)
# 日期：2024-04-26
# ================================================================================

# -------------------- 0. 包加载 --------------------
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(VariantAnnotation)
  library(gwasglue)
  library(dplyr)
  library(CMplot)
  library(knitr)
  library(ggplot2)
})

cat("\n========== 包版本信息 ==========\n")
cat(paste("TwoSampleMR version:", packageVersion("TwoSampleMR"), "\n\n"))

# -------------------- 1. 工作目录与认证配置 --------------------
setwd("F:/MR/Interleukin-6 levels and keloids")
cat(paste("工作目录:", getwd(), "\n\n"))

# 配置 API 认证
api_token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJramJjbGxAMTYzLmNvbSIsImlhdCI6MTc3NjY5NjQ4OSwiZXhwIjoxNzc3OTA2MDg5fQ.bJ0zczsIGbnRcGrDE1FWTvDXsE9PZJfWJ1_RyQrIbvSToHWD_1JASJFgqcv4VEYVF5QB2C_M1hyVG6ZyBQPSw8N-f5lIKBvlSfs5nZ9ZMo3aLQ1-4Fgwuapj2V1_5FPpHknkzXbg-IeVFLxYuefz5b7CjVNF_oEP9YImKSgrTfS91fFiEiqJ4rdrE28SXqWqC33u-O2t5-pa3DlZIREujnuG9EE7qOia9iHCoSXMEQoesYyzFTNVkfJ84iQN0EpX5dIRzdCyFgS7JPg2c1wAU_dCPem_AfKPKDargnBV6GXSTKVZp5UAg2enZz6OMsQlSN5VGxQENzlEhQVyyoLONQ"
Sys.setenv(OPENGWAS_JWT = api_token)
cat("API 令牌配置完成\n\n")

# -------------------- 2. 参数设置 --------------------
clump_kb <- 10000  # Clumping 半径 (kb)
Ffilter <- 10      # F 统计量阈值

# -------------------- 3. 数据加载和处理 --------------------
cat("\033[1m========== 步骤 3: 数据加载和预处理 ==========\033[0m\n")

# 暴露数据
exposure_vcf <- "ebi-a-GCST90012005.vcf.gz"
exposure_data <- tryCatch({
  readVcf(exposure_vcf)
}, error = function(e) {
  stop(paste("无法读取暴露 VCF 文件:", exposure_vcf, "\n", e$message))
})

exposure_data <- gwasvcf_to_TwoSampleMR(vcf = exposure_data, type = "exposure")

# 结局数据
outcome_vcf <- "out-ebi-a-GCST90018874.vcf.gz"
outcome_data <- tryCatch({
  readVcf(outcome_vcf)
}, error = function(e) {
  stop(paste("无法读取结局 VCF 文件:", outcome_vcf, "\n", e$message))
})

outcome_data <- gwasvcf_to_TwoSampleMR(vcf = outcome_data, type = "outcome")

# -------------------- 4. 暴露数据 Clumping --------------------
cat("\033[1m========== 步骤 4: 暴露数据 Clumping ==========\033[0m\n")

# 提取需要的列
exp_raw <- exposure_data[, c("SNP", "chr.exposure", "pos.exposure",
                             "beta.exposure", "se.exposure",
                             "effect_allele.exposure", "other_allele.exposure",
                             "eaf.exposure", "samplesize.exposure", "pval.exposure")]
colnames(exp_raw) <- c("SNP", "CHR", "BP", "beta", "se",
                       "effect_allele", "other_allele", "eaf", "samplesize", "pval")

# 按 P 值排序
exp_raw <- exp_raw[order(exp_raw$pval), ]

# Clumping
keep_idx <- c()
last_chr <- 0
last_bp <- 0

for(i in 1:nrow(exp_raw)) {
  curr_chr <- exp_raw$CHR[i]
  curr_bp <- exp_raw$BP[i]
  
  if(curr_chr != last_chr || abs(curr_bp - last_bp) > clump_kb * 1000) {
    keep_idx <- c(keep_idx, i)
    last_chr <- curr_chr
    last_bp <- curr_bp
  }
}

exp_clumped <- exp_raw[keep_idx, ]

cat(paste("✅ Clumping 后 SNP 数量:", nrow(exp_clumped), "\n\n"))

# -------------------- 5. 数据协调 --------------------
cat("\033[1m========== 步骤 5: 数据协调 ==========\033[0m\n")

# 提取匹配的 SNP
matched_snps <- outcome_data[outcome_data$SNP %in% exp_clumped$SNP, ]

# 合并数据
merged_data <- merge(exp_clumped, matched_snps, by = "SNP", suffixes = c(".exp", ".out"))

cat(paste("✅ 合并数据后 SNP 数量:", nrow(merged_data), "\n"))

# 构建 harmonised 数据框
dat <- data.frame(
  SNP = merged_data$SNP,
  chr.exposure = merged_data$chr.exposure,
  pos.exposure = merged_data$pos.exposure,
  beta.exposure = merged_data$beta.exposure,
  se.exposure = merged_data$se.exposure,
  effect_allele.exposure = merged_data$effect_allele.exposure,
  other_allele.exposure = merged_data$other_allele.exposure,
  eaf.exposure = merged_data$eaf.exposure,
  pval.exposure = merged_data$pval.exposure,
  beta.outcome = merged_data$beta.outcome,
  se.outcome = merged_data$se.outcome,
  effect_allele.outcome = merged_data$effect_allele.outcome,
  other_allele.outcome = merged_data$other_allele.outcome,
  eaf.outcome = merged_data$eaf.outcome,
  pval.outcome = merged_data$pval.outcome,
  samplesize.exposure = merged_data$samplesize.exposure,
  samplesize.outcome = merged_data$samplesize.outcome,
  id.exposure = "IL6",
  exposure = "Interleukin-6 levels",
  id.outcome = "Keloid",
  outcome = "Keloids",
  palindromic = FALSE,
  ambiguous = FALSE,
  stringsAsFactors = FALSE
)

#  确保 mr_keep 列存在
dat$mr_keep <- TRUE
cat("✅ 确保 mr_keep 列已赋值为 TRUE\n")

# -------------------- 6. MR 分析 --------------------
cat("\033[1m========== 步骤 6: MR 分析 ==========\033[0m\n")

n_snps <- nrow(dat)
cat(paste("用于 MR 分析的 SNP 数量:", n_snps, "\n\n"))

mrResult <- tryCatch({
  mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
}, error = function(e) {
  stop(paste("MR 分析失败:", e$message))
})

# ================================================================================
# 步骤 7: 绘制图表
# ================================================================================
cat("\033[1m========== 步骤 7: 绘制图表 ==========\033[0m\n")
#绘制散点图
tryCatch({
  p <- ggplot(dat, aes(x = beta.exposure, y = beta.outcome / se.outcome)) +
    geom_point(size = 3, color = "steelblue", alpha = 0.7) +
    labs(title = "IL-6 → Keloids\nMendelian Randomization Scatter Plot") +
    theme_minimal()
  ggsave("pic.scatter_plot.png", plot = p, width = 8, height = 7, dpi = 150)
  cat("✓ 散点图已保存：pic.scatter_plot.png\n")
}, error = function(e) {
  cat("⚠️  散点图绘制失败:", e$message, "\n")
})

#森林图
tryCatch({
  res_single=mr_singlesnp(dat)
  mr_forest_plot(res_single)
  dev.off()
  cat("✓ 森林图已保存\n")
}, error = function(e) {
  cat("⚠️  森林图绘制失败:", e$message, "\n")
})

#漏斗图
tryCatch({
  mr_funnel_plot(singlesnp_results = mr_leaveoneout(dat))
  dev.off()
  cat("✓ 漏斗图已保存\n")
}, error = function(e) {
  cat("⚠️  漏斗图绘制失败:", e$message, "\n")
})

#留一法敏感性分析
tryCatch({
  mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
  dev.off()
  cat("✓ 留一法图已保存：pic.leaveoneout.png\n")
}, error = function(e) {
  cat("⚠️  留一法图绘制失败:", e$message, "\n")
})

# ================================================================================
# 步骤 8: 结果汇总
# ================================================================================
cat("\033[1m========== 步骤 8: 结果汇总 ==========\033[0m\n\n")

# 读取结果文件
mr_results <- read.csv("table.MRresult.csv")
if(file.exists("table.heterogeneity.csv")) {
  heter_results <- read.csv("table.heterogeneity.csv")
}
if(file.exists("table.pleiotropy.csv")) {
  pleio_results <- read.csv("table.pleiotropy.csv")
}

# 8.1 主要 MR 结果
cat("\033[1m========== 孟德尔随机化主要结果 ==========\033[0m\n")
ivw_result <- subset(mr_results, method == "Inverse variance weighted")

if(nrow(ivw_result) > 0) {
  cat("\nIVW 方法分析结果:\n")
  print(kable(ivw_result[, c("method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95")],
              format = "pipe"))
} else {
  cat("无 IVW 结果（SNP 数量不足）\n")
  print(kable(as.data.frame(mrResult)[, c("method", "nsnp", "b", "se", "pval")],
              format = "pipe"))
}

# 8.2 异质性检验
if(!is.null(heter_results) && nrow(heter_results) > 0) {
  cat("\n\033[1m========== 异质性检验 (Cochran's Q) ==========\033[0m\n")
  ivw_heter <- subset(heter_results, method == "Inverse variance weighted")
  if(nrow(ivw_heter) > 0) {
    print(kable(ivw_heter[, c("method", "Q", "Q_df", "Q_pval")],
                format = "pipe"))
    if(ivw_heter$Q_pval[1] < 0.05) {
      cat("\n⚠️  存在显著异质性（P < 0.05），需谨慎解释结果\n")
    } else {
      cat("\n✓ 未发现显著异质性\n")
    }
  }
}

# 8.3 多效性检验
if(!is.null(pleio_results) && nrow(pleio_results) > 0) {
  cat("\n\033[1m========== 多效性检验 (MR-Egger Intercept) ==========\033[0m\n")
  print(kable(pleio_results[, c("egger_intercept", "se", "pval")],
              format = "pipe"))
  if(pleio_results$pval[1] < 0.05) {
    cat("\n⚠️  存在显著水平多效性（P < 0.05），建议使用 MR-Egger 结果\n")
  } else {
    cat("\n✓ 未发现显著水平多效性\n")
  }
}

# 8.4 Steiger 方向性检验
if(file.exists("table.steiger.csv")) {
  cat("\n\033[1m========== Steiger 方向性检验 ==========\033[0m\n")
  steiger_results <- read.csv("table.steiger.csv")
  print(kable(steiger_results[, c("exposure", "outcome", "correct_causal_direction")],
              format = "pipe"))
}

# 8.5 最终结论
cat("\n\033[1m========== 分析结论 ==========\033[0m\n")

if(nrow(ivw_result) > 0) {
  pval <- ivw_result$pval[1]
  b <- ivw_result$b[1]
  se <- ivw_result$se[1]
  or <- ivw_result$or[1]
  ci_low <- ivw_result$or_lci95[1]
  ci_high <- ivw_result$or_uci95[1]
  
  cat(sprintf("\n主要发现 (IVW 方法):\n"))
  cat(sprintf("   β = %.4f, SE = %.4f\n", b, se))
  cat(sprintf("   OR = %.3f (95%% CI: %.3f - %.3f)\n", or, ci_low, ci_high))
  cat(sprintf("   P = %.4f\n\n", pval))
  
  if(pval < 0.05) {
    cat(sprintf("\033[32m✅ 在α=0.05 水平上，IL-6 水平对瘢痕疙瘩有显著因果效应\033[0m\n"))
    if(or > 1) {
      cat(sprintf("   解读：IL-6 水平升高增加瘢痕疙瘩风险 (OR = %.3f)\n", or))
    } else {
      cat(sprintf("   解读：IL-6 水平升高降低瘢痕疙瘩风险 (OR = %.3f)\n", or))
    }
  } else {
    cat(sprintf("\033[33m⚠️  在α=0.05 水平上，未检测到 IL-6 水平对瘢痕疙瘩的显著因果效应\033[0m\n"))
    cat(sprintf("   OR = %.3f (95%% CI: %.3f - %.3f), P = %.4f\n", or, ci_low, ci_high, pval))
    if(n_snps < 10) {
      cat("   ⚠️  注意：SNP 数量有限，统计效能可能不足\n")
    }
  }
} else {
  cat("\033[31m⚠️  无法得出结论：MR 分析结果不足\033[0m\n")
}

# 8.6 文件保存清单
cat("\n\033[1m========== 所有结果文件已保存 ==========\033[0m\n")
cat("\n表格文件:\n")
file_list <- c(
  "table.MRresult.csv",
  "table.heterogeneity.csv",
  "table.pleiotropy.csv",
  "table.steiger.csv",
  "table.leaveoneout.csv"
)
for(f in file_list) {
  if(file.exists(f)) {
    cat(sprintf("  ✓ %s\n", f))
  }
}

cat("\n图表文件:\n")
plot_list <- c(
  "pic.scatter_plot.png",
  "pic.forest.png",
  "pic.funnel_plot.png",
  "pic.leaveoneout.png"
)
for(f in plot_list) {
  if(file.exists(f)) {
    cat(sprintf("  ✓ %s\n", f))
  }
}

# ================================================================================
# 结束
# ================================================================================
cat("\n\033[1m========== 分析完成 ==========\033[0m\n")
cat(sprintf("完成时间：%s\n", Sys.time()))
cat(sprintf("研究：Interleukin-6 levels → Keloids\n"))
cat(sprintf("SNP 数量：%d\n", n_snps))
