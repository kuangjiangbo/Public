# ======================================================================
# Colocalization Analysis: IL-6 Levels ↔ Keloid Risk
# ======================================================================
# 使用 Bayesian colocalization (coloc 包) 检验 IL-6 GWAS 和 Keloid GWAS
# 在 IL6R 区域 (chr1) 和 HLA 区域 (chr6) 是否共享同一个因果变异
#
# 数据来源：
#   IL-6 GWAS:    GCST90012005 (Folkersen L et al., Nat Metab, 2020)
#                 21,758 European ancestry individuals
#                 PMID: 33067605
#   Keloid GWAS:  GCST90018874 (Sakaue S et al., Nat Genet, 2021)
#                 668 cases / 481,244 controls (European subset)
#                 PMID: 34594039
#
# 工作目录: F:/MR/Interleukin-6 levels and keloids/数据/第四次/
# 输出文件:
#   table.colocalization_summary.csv  - 共定位结果汇总表
#   pic.coloc_il6r.png               - IL6R 区域关联图
#   pic.coloc_hla.png                - HLA 区域关联图
#   colocalization_session_info.txt  - 运行日志
#
# 依赖包: coloc, VariantAnnotation, Rsamtools, GenomicRanges, 
#          dplyr, ggplot2, tidyr
# ======================================================================

# ======================================================================
# 0. 工作目录与参数设置
# ======================================================================
setwd("F:/MR/Interleukin-6 levels and keloids/数据/第四次/")
cat("=============================================\n")
cat("Colocalization Analysis: IL-6 ↔ Keloids\n")
cat("=============================================\n\n")
cat("Working directory:", getwd(), "\n")

# ---- Keloid GWAS 病例比例 (case proportion) ----
# European subset: 668 cases / 481,244 controls = 481,912 total
KELOID_CASE_PROP <- 668 / 481912  # = 0.001386
cat("Keloid case proportion (s):", 
    sprintf("%.6f (%d / %d)", KELOID_CASE_PROP, 668, 481912), "\n\n")

# ---- IL-6 GWAS 样本量 ----
IL6_N <- 21758
KELOID_N <- 481912

# ======================================================================
# 1. 安装 / 加载 R 包
# ======================================================================
required_pkgs <- c("coloc", "VariantAnnotation", "Rsamtools",
                    "GenomicRanges", "dplyr", "ggplot2", "tidyr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

cat("All packages loaded successfully.\n\n")

# ======================================================================
# 2. 定义共定位分析的两个基因组区域
# ======================================================================
# IL6R 区域: chr1, 以 IL6R 基因为中心 ±500kb
#   IL6R 基因位置: ~chr1:154,400,000-154,430,000
#   分析窗口:      chr1:154,000,000-155,000,000 (1Mb)
il6r_region <- GRanges(
  seqnames = "1",
  ranges   = IRanges(start = 154000000, end = 155000000)
)

# HLA 区域: chr6, 覆盖 MHC Class III 区域
#   分析窗口: chr6:32,000,000-33,000,000 (1Mb)
hla_region <- GRanges(
  seqnames = "6",
  ranges   = IRanges(start = 32000000, end = 33000000)
)

regions <- list(
  IL6R = il6r_region,
  HLA  = hla_region
)

cat("Analysis regions:\n")
cat("  IL6R: chr1:154,000,000-155,000,000\n")
cat("  HLA:  chr6:32,000,000-33,000,000\n\n")

# ======================================================================
# 3. 定义 VCF 文件路径
# ======================================================================
exposure_vcf <- "ebi-a-GCST90012005.vcf.gz"    # IL-6 暴露 GWAS
outcome_vcf  <- "out-ebi-a-GCST90018874.vcf.gz"  # Keloid 结局 GWAS

vcf_files <- list(
  IL6    = exposure_vcf,
  Keloid = outcome_vcf
)

cat("VCF files:\n")
cat("  IL-6 (exposure):   ", exposure_vcf, " (", 
    file.info(exposure_vcf)$size / 1e6, "MB)\n")
cat("  Keloid (outcome):  ", outcome_vcf, " (", 
    file.info(outcome_vcf)$size / 1e6, "MB)\n\n")

# ======================================================================
# 4. 创建 Tabix 索引（如果尚未创建）
# ======================================================================
for (f in vcf_files) {
  tbi_file <- paste0(f, ".tbi")
  if (!file.exists(tbi_file)) {
    cat("Creating tabix index for:", f, "...\n")
    indexTabix(f, format = "vcf")
    cat("  Done:", tbi_file, "\n")
  } else {
    cat("Tabix index exists:", tbi_file, "\n")
  }
}
cat("\n")

# ======================================================================
# 5. 定义函数: 从 VCF 中提取指定区域的 GWAS Summary 数据
# ======================================================================
# IEU OpenGWAS VCF 格式说明:
#   Genotype 字段 (不是 INFO 字段):
#     ES = effect size (beta)
#     SE = standard error
#     LP = -log10(p-value)
#     AF = effect allele frequency
#     SS = sample size (可能为空)
#     NC = number of cases (仅结局可用)

extract_gwas_region <- function(vcf_path, region_gr, trait_name = "GWAS") {
  
  # 打开 tabix 索引文件
  tab <- TabixFile(vcf_path)
  open(tab)
  
  cat("  Reading", trait_name, "from region", 
      as.character(seqnames(region_gr)), ":",
      start(region_gr), "-", end(region_gr), "...\n")
  
  vcf <- readVcf(tab, param = region_gr)
  close(tab)
  
  # 检查是否找到变异
  if (nrow(vcf) == 0) {
    cat("  WARNING: No variants found in this region!\n")
    return(NULL)
  }
  cat("  Found", nrow(vcf), "variants.\n")
  
  # 提取基本信息
  fixed <- rowRanges(vcf)
  snp_names <- names(fixed)
  
  # 染色体名称处理
  sn_rle <- seqnames(fixed)
  if (is(sn_rle, "Rle")) {
    chr <- rep(as.character(runValue(sn_rle)), times = runLength(sn_rle))
  } else {
    chr <- as.character(sn_rle)
  }
  pos <- start(fixed)
  
  # IEU VCF 的数据存储在 Genotype 字段
  geno_data <- geno(vcf)
  cat("  Genotype fields:", paste(names(geno_data), collapse = ", "), "\n")
  
  # 辅助函数: 从 matrix-of-lists 中提取数值列
  extract_col <- function(mat, col = 1) {
    sapply(seq_len(nrow(mat)), function(i) {
      val <- mat[i, col]
      if (is.null(val) || length(val) == 0) NA else as.numeric(val)
    })
  }
  
  # 提取 beta, se, p-value, allele frequency
  beta <- extract_col(geno_data$ES)
  se   <- extract_col(geno_data$SE)
  lp   <- extract_col(geno_data$LP)
  pval <- 10^(-lp)
  af   <- extract_col(geno_data$AF)
  
  cat("  Valid ES:", sum(!is.na(beta)), "/", length(beta), "\n")
  
  # 提取效应等位基因 (ALT) 和参考等位基因 (REF)
  effect_allele <- sapply(fixed$ALT, function(x) as.character(x[[1]]))
  other_allele  <- as.character(fixed$REF)
  
  # 筛选有效变异
  valid <- !is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(af)
  cat("  Valid variants (complete data):", sum(valid), "/", length(valid), "\n")
  
  if (sum(valid) < 10) {
    cat("  WARNING: Too few valid variants!\n")
    return(NULL)
  }
  
  # 构建结果数据框
  result <- data.frame(
    SNP            = snp_names[valid],
    CHR            = chr[valid],
    BP             = pos[valid],
    effect_allele  = effect_allele[valid],
    other_allele   = other_allele[valid],
    beta           = beta[valid],
    se             = se[valid],
    pval           = pval[valid],
    LP             = lp[valid],
    MAF            = pmin(af[valid], 1 - af[valid]),
    EAF            = af[valid],
    stringsAsFactors = FALSE
  )
  
  # 去重
  result <- result[!duplicated(result$SNP), ]
  
  cat("  Final variants:", nrow(result), "\n")
  return(result)
}

# ======================================================================
# 6. 提取两个 GWAS 在两个区域的数据
# ======================================================================
results_list <- list()

for (region_name in names(regions)) {
  cat("\n========= Processing region:", region_name, "=========\n\n")
  
  region_gr <- regions[[region_name]]
  
  # 提取 IL-6 暴露数据
  il6_data <- extract_gwas_region(
    vcf_path   = vcf_files$IL6,
    region_gr  = region_gr,
    trait_name = paste0("IL-6 (", region_name, ")")
  )
  
  # 提取 Keloid 结局数据
  keloid_data <- extract_gwas_region(
    vcf_path   = vcf_files$Keloid,
    region_gr  = region_gr,
    trait_name = paste0("Keloid (", region_name, ")")
  )
  
  results_list[[region_name]] <- list(
    IL6    = il6_data,
    Keloid = keloid_data
  )
}

# ======================================================================
# 7. 合并两个 GWAS 的重叠 SNP
# ======================================================================
cat("\n========= Merging overlapping SNPs =========\n\n")

for (region_name in names(results_list)) {
  cat("--- Region:", region_name, "---\n")
  
  il6_df    <- results_list[[region_name]]$IL6
  keloid_df <- results_list[[region_name]]$Keloid
  
  if (is.null(il6_df) || is.null(keloid_df)) {
    cat("  Skipping - data missing.\n\n")
    next
  }
  
  common_snps <- intersect(il6_df$SNP, keloid_df$SNP)
  cat("  IL-6 SNPs:", nrow(il6_df), "\n")
  cat("  Keloid SNPs:", nrow(keloid_df), "\n")
  cat("  Overlapping SNPs:", length(common_snps), "\n")
  
  # 合并
  merged <- il6_df %>%
    inner_join(keloid_df, by = "SNP", suffix = c(".il6", ".kel"))
  
  cat("  Merged SNPs for coloc:", nrow(merged), "\n\n")
  
  results_list[[region_name]]$merged <- merged
}

# ======================================================================
# 8. 运行 Bayesian Colocalization (coloc.abf)
# ======================================================================
cat("\n========================================\n")
cat("Running Bayesian Colocalization (coloc)\n")
cat("========================================\n\n")

coloc_results <- list()

for (region_name in names(results_list)) {
  cat("--- Region:", region_name, "---\n")
  
  merged <- results_list[[region_name]]$merged
  
  if (is.null(merged) || nrow(merged) < 50) {
    n <- ifelse(is.null(merged), 0, nrow(merged))
    cat("  Skipping - insufficient SNPs (n =", n, ", need >= 50)\n\n")
    next
  }
  
  cat("  Preparing coloc datasets...\n")
  cat("  N SNPs:", nrow(merged), "\n")
  
  # ---------------------------
  # Dataset 1: IL-6 (连续变量)
  # ---------------------------
  dataset1 <- list(
    beta     = merged$beta.il6,
    varbeta  = merged$se.il6^2,
    snp      = merged$SNP,
    position = merged$BP.il6,
    type     = "quant",
    N        = IL6_N,
    MAF      = merged$MAF.il6
  )
  
  # ---------------------------
  # Dataset 2: Keloid (二分类变量)
  #   病例比例 s = 668/481912 = 0.001386
  # ---------------------------
  dataset2 <- list(
    beta     = merged$beta.kel,
    varbeta  = merged$se.kel^2,
    snp      = merged$SNP,
    position = merged$BP.kel,
    type     = "cc",
    N        = KELOID_N,
    s        = KELOID_CASE_PROP,
    MAF      = merged$MAF.kel
  )
  
  # 删除 NA 位置
  valid_pos <- !is.na(dataset1$position) & !is.na(dataset2$position) &
               !is.na(dataset1$MAF) & !is.na(dataset2$MAF)
               
  for (nm in names(dataset1)) {
    if (length(dataset1[[nm]]) == length(valid_pos)) {
      dataset1[[nm]] <- dataset1[[nm]][valid_pos]
    }
  }
  for (nm in names(dataset2)) {
    if (length(dataset2[[nm]]) == length(valid_pos)) {
      dataset2[[nm]] <- dataset2[[nm]][valid_pos]
    }
  }
  
  cat("  Running coloc.abf with", sum(valid_pos), "SNPs...\n")
  
  # 执行贝叶斯共定位
  result <- tryCatch({
    coloc.abf(dataset1, dataset2)
  }, error = function(e) {
    cat("  ERROR in coloc.abf:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    coloc_results[[region_name]] <- result
    
    pp_h0 <- result$summary["PP.H0.abf"]
    pp_h1 <- result$summary["PP.H1.abf"]
    pp_h2 <- result$summary["PP.H2.abf"]
    pp_h3 <- result$summary["PP.H3.abf"]
    pp_h4 <- result$summary["PP.H4.abf"]
    
    cat("\n  ====== Colocalization Results ======\n")
    cat("  PP.H0 (no association):       ", sprintf("%.6f", pp_h0), "\n")
    cat("  PP.H1 (IL-6 only):            ", sprintf("%.6f", pp_h1), "\n")
    cat("  PP.H2 (Keloid only):          ", sprintf("%.6f", pp_h2), "\n")
    cat("  PP.H3 (both, distinct SNPs):  ", sprintf("%.6f", pp_h3), "\n")
    cat("  PP.H4 (shared causal variant):", sprintf("%.6f", pp_h4), "\n")
    
    if (pp_h4 > 0.8) {
      cat("\n  >>> CONCLUSION: Strong evidence for shared causal variant (PP.H4 > 0.8)\n")
    } else if (pp_h4 > 0.5) {
      cat("\n  >>> CONCLUSION: Moderate evidence for shared causal variant\n")
    } else {
      cat("\n  >>> CONCLUSION: Weak evidence (PP.H4 < 0.5).\n")
      if (pp_h3 > pp_h4) {
        cat("  >>> PP.H3 > PP.H4: suggests distinct causal variants at this locus\n")
      }
    }
    cat("\n")
  }
}

# ======================================================================
# 9. 敏感性分析: 不同先验概率组合
# ======================================================================
cat("\n========= Sensitivity Analysis =========\n\n")

if (length(coloc_results) > 0) {
  for (region_name in names(coloc_results)) {
    cat("--- Region:", region_name, "---\n")
    
    merged <- results_list[[region_name]]$merged
    if (is.null(merged) || nrow(merged) < 50) next
    
    # 5 种先验概率组合
    prior_combos <- list(
      c(1e-4, 1e-4, 1e-5),   # (1) Default coloc priors
      c(1e-4, 1e-4, 1e-4),   # (2) More permissive shared
      c(1e-4, 1e-4, 1e-6),   # (3) More stringent shared
      c(1e-3, 1e-4, 1e-5),   # (4) Higher p1 (IL-6)
      c(1e-4, 1e-3, 1e-5)    # (5) Higher p2 (Keloid)
    )
    
    prior_labels <- c(
      "Default (1e-4/1e-4/1e-5)",
      "Permissive shared",
      "Stringent shared",
      "High p1 (IL-6)",
      "High p2 (Keloid)"
    )
    
    cat(sprintf("  %-25s %s\n", "Prior setting", "PP.H4"))
    cat("  ", paste(rep("-", 40), collapse = ""), "\n")
    
    for (i in seq_along(prior_combos)) {
      pr <- prior_combos[[i]]
      
      d1 <- list(
        beta     = merged$beta.il6,
        varbeta  = merged$se.il6^2,
        snp      = merged$SNP,
        position = merged$BP.il6,
        type     = "quant",
        N        = IL6_N,
        MAF      = merged$MAF.il6
      )
      d2 <- list(
        beta     = merged$beta.kel,
        varbeta  = merged$se.kel^2,
        snp      = merged$SNP,
        position = merged$BP.kel,
        type     = "cc",
        N        = KELOID_N,
        s        = KELOID_CASE_PROP,
        MAF      = merged$MAF.kel
      )
      
      v <- !is.na(d1$position) & !is.na(d2$MAF) & !is.na(d1$MAF)
      for (nm in names(d1)) if (length(d1[[nm]]) == length(v)) d1[[nm]] <- d1[[nm]][v]
      for (nm in names(d2)) if (length(d2[[nm]]) == length(v)) d2[[nm]] <- d2[[nm]][v]
      
      rs <- coloc.abf(d1, d2, p1 = pr[1], p2 = pr[2], p12 = pr[3])
      cat(sprintf("  %-25s %.6f\n", prior_labels[i], rs$summary["PP.H4.abf"]))
    }
    cat("\n")
  }
} else {
  cat("No coloc results to analyze for sensitivity.\n\n")
}

# ======================================================================
# 10. 生成区域关联图 (Regional Association Plot)
# ======================================================================
cat("\n========= Generating Regional Plots =========\n\n")

for (region_name in names(coloc_results)) {
  merged <- results_list[[region_name]]$merged
  result <- coloc_results[[region_name]]
  
  if (is.null(merged)) next
  
  top_il6 <- merged$SNP[which.min(merged$pval.il6)]
  top_kel <- merged$SNP[which.min(merged$pval.kel)]
  
  cat("  Region:", region_name, "\n")
  cat("  Top IL-6 SNP:  ", top_il6, "(P =", format(min(merged$pval.il6), digits = 3), ")\n")
  cat("  Top Keloid SNP:", top_kel, "(P =", format(min(merged$pval.kel), digits = 3), ")\n")
  
  # 数据整理
  plot_data <- merged %>%
    dplyr::select(SNP, BP = BP.il6, 
                  LP_IL6 = LP.il6, LP_Keloid = LP.kel) %>%
    tidyr::pivot_longer(
      cols = c(LP_IL6, LP_Keloid),
      names_to = "Trait", values_to = "LP"
    ) %>%
    mutate(Trait = recode(Trait,
                          LP_IL6    = "IL-6 levels",
                          LP_Keloid = "Keloids"))
  
  pp4_val <- sprintf("%.4f", result$summary["PP.H4.abf"])
  
  p <- ggplot(plot_data, aes(x = BP / 1e6, y = LP, color = Trait)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_color_manual(
      values = c("IL-6 levels" = "#E41A1C", "Keloids" = "#377EB8")
    ) +
    labs(
      title    = paste("Regional Association Plot -", region_name, "Locus"),
      subtitle = paste0(
        "chr", gsub("chr", "", unique(merged$CHR.il6)), ":",
        round(min(merged$BP.il6) / 1e6, 3), "-",
        round(max(merged$BP.il6) / 1e6, 3), " Mb  |  PP.H4 = ", pp4_val
      ),
      x = paste0("Position on chromosome ", unique(merged$CHR.il6), " (Mb)"),
      y = expression(-log[10](italic(P)))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  plot_file <- paste0("pic.coloc_", tolower(region_name), ".png")
  ggsave(plot_file, plot = p, width = 10, height = 6, dpi = 150)
  cat("  Plot saved:", plot_file, "\n\n")
}

# ======================================================================
# 11. 结果汇总与保存
# ======================================================================
cat("\n========================================\n")
cat("Colocalization Summary Table\n")
cat("========================================\n\n")

summary_table <- data.frame(
  Region     = character(),
  N_SNPs     = integer(),
  PP.H0      = character(),
  PP.H1      = character(),
  PP.H2      = character(),
  PP.H3      = character(),
  PP.H4      = character(),
  Conclusion = character(),
  stringsAsFactors = FALSE
)

for (region_name in names(coloc_results)) {
  result <- coloc_results[[region_name]]
  merged <- results_list[[region_name]]$merged
  
  pp4 <- result$summary["PP.H4.abf"]
  pp3 <- result$summary["PP.H3.abf"]
  
  conclusion <- if (pp4 > 0.8) {
    "Shared causal variant (strong)"
  } else if (pp4 > 0.5) {
    "Moderate evidence for shared variant"
  } else if (pp3 > pp4) {
    "Likely distinct causal variants"
  } else {
    "Inconclusive"
  }
  
  summary_table <- rbind(summary_table, data.frame(
    Region     = region_name,
    N_SNPs     = nrow(merged),
    PP.H0      = sprintf("%.6f", result$summary["PP.H0.abf"]),
    PP.H1      = sprintf("%.6f", result$summary["PP.H1.abf"]),
    PP.H2      = sprintf("%.6f", result$summary["PP.H2.abf"]),
    PP.H3      = sprintf("%.6f", result$summary["PP.H3.abf"]),
    PP.H4      = sprintf("%.6f", pp4),
    Conclusion = conclusion,
    stringsAsFactors = FALSE
  ))
}

print(summary_table, row.names = FALSE)

# 保存结果表格
write.csv(summary_table, "table.colocalization_summary.csv", row.names = FALSE)
cat("\nSummary saved to: table.colocalization_summary.csv\n")

# ======================================================================
# 12. 文稿建议文本 (Interpretation Guide)
# ======================================================================
cat("\n========================================\n")
cat("Manuscript Interpretation Guide\n")
cat("========================================\n\n")

for (i in seq_len(nrow(summary_table))) {
  pp4 <- as.numeric(summary_table$PP.H4[i])
  reg <- summary_table$Region[i]
  
  cat("Region:", reg, "\n")
  cat("PP.H4 =", summary_table$PP.H4[i], "\n\n")
  
  if (pp4 > 0.8) {
    cat('  => "Bayesian colocalization at the ', reg, ' locus\n', sep = "")
    cat('      demonstrated strong evidence for a shared causal variant\n', sep = "")
    cat('      (PP.H4 = ', sprintf("%.4f", pp4), '), supporting the MR finding\n', sep = "")
    cat('      and implicating shared genetic regulation at this locus."\n\n')
  } else if (pp4 > 0.5) {
    cat('  => "Colocalization at the ', reg, ' locus provided moderate support\n', sep = "")
    cat('      for a shared causal variant (PP.H4 = ', sprintf("%.4f", pp4), ').\n', sep = "")
    cat('      Further fine-mapping studies are warranted."\n\n')
  } else {
    cat('  => "Colocalization at the ', reg, ' locus did not provide strong evidence\n', sep = "")
    cat('      for a shared causal variant (PP.H4 = ', sprintf("%.4f", pp4), ').\n', sep = "")
    cat('      The MR association at this locus may be influenced by\n', sep = "")
    cat('      distinct genetic mechanisms or residual confounding."\n\n')
  }
}

# ======================================================================
# 13. 保存运行日志
# ======================================================================
sink("colocalization_session_info.txt")
cat("Colocalization Analysis Session Info\n")
cat("====================================\n\n")
cat("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("GWAS data:\n")
cat("  IL-6 (exposure):  GCST90012005\n")
cat("    Folkersen L et al. Nat Metab. 2020. PMID: 33067605\n")
cat("    N =", IL6_N, ", European\n\n")
cat("  Keloid (outcome): GCST90018874\n")
cat("    Sakaue S et al. Nat Genet. 2021. PMID: 34594039\n")
cat("    N =", KELOID_N, ", cases = 668, controls = 481244\n")
cat("    Case proportion (s) =", KELOID_CASE_PROP, "\n\n")
cat("Analysis regions:\n")
for (nm in names(regions)) {
  cat("  ", nm, ":", as.character(seqnames(regions[[nm]])), ":",
      start(regions[[nm]]), "-", end(regions[[nm]]), "\n")
}
cat("\nR version:", R.version.string, "\n\n")
cat("R packages:\n")
for (pkg in required_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("  -", pkg, ":", as.character(packageVersion(pkg)), "\n")
  }
}
cat("\nResults:\n")
if (length(coloc_results) > 0) {
  for (nm in names(coloc_results)) {
    r <- coloc_results[[nm]]
    cat("  ", nm, ": PP.H4 =", sprintf("%.6f", r$summary["PP.H4.abf"]), "\n")
  }
} else {
  cat("  No colocalization results generated.\n")
}
sink()
cat("Session info saved to: colocalization_session_info.txt\n")

# ======================================================================
# 14. 完成
# ======================================================================
cat("\n========================================\n")
cat("Colocalization Analysis Complete!\n")
cat("========================================\n")
cat("Output files:\n")
cat("  table.colocalization_summary.csv\n")
cat("  pic.coloc_il6r.png\n")
cat("  pic.coloc_hla.png\n")
cat("  colocalization_session_info.txt\n")
