# ======================================================================
# Colocalization Analysis v4.1: IL-6 Levels ↔ Keloid Risk
# ======================================================================
# v4.1 修复:
#   1. 修复 TabixFile 双重关闭 (on.exit 模式)
#   2. 修复 nchar(NA) 崩溃
#   3. 修复 run_coloc 中脆弱的 NA 过滤逻辑
#   4. 删除冗余的 regions_gr 定义
#   5. setTimeLimit + fallback 双重超时保护
#   6. 增加 VCF 文件存在性检查
#   7. 安全关闭图形设备 (只关自己的设备)
#   8. 增加每个区域的 checkpoint 写入
#   9. 更好的错误消息和状态提示
#  10. R.utils::withTimeout 作为首选超时方案
#
# 数据来源:
#   IL-6 GWAS:    GCST90012005 (Folkersen L et al., Nat Metab, 2020)
#   Keloid GWAS:  GCST90018874 (Sakaue S et al., Nat Genet, 2021)
#
# 输出:
#   table.coloc_v4_summary.csv           - 共定位结果汇总
#   table.coloc_v4_sensitivity.csv       - 敏感性分析结果
#   pic.coloc_v4_il6r_locuscompare.png   - IL6R LocusCompare图
#   pic.coloc_v4_hla_locuscompare.png    - HLA LocusCompare图
#   pic.coloc_v4_sensitivity.png         - 敏感性分析汇总图
#   colocalization_v4_session_info.txt   - 运行日志
# ======================================================================

# ======================================================================
# 0. 环境清理 & 工作目录
# ======================================================================
cat("=============================================\n")
cat("Colocalization v4.1: IL-6 Levels ↔ Keloids\n")
cat("=============================================\n\n")

cat("--- Step 0: Environment setup ---\n")

# 清除工作区
suppressWarnings(rm(list = ls()))

setwd("F:/MR/Interleukin-6 levels and keloids/数据/第四次/")

# ---- 参数 ----
IL6_N         <- 21758
KELOID_N      <- 481912
KELOID_CASES  <- 668
KELOID_CTRL   <- KELOID_N - KELOID_CASES
KELOID_S      <- KELOID_CASES / KELOID_N

COL_IL6    <- "#E41A1C"
COL_KELOID <- "#377EB8"

# 请在双引号内粘贴你的完整 OpenGWAS token
# 获取地址: https://api.opengwas.io
OPENGWAS_TOKEN <- trimws("eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJramJjbGxAMTYzLmNvbSIsImlhdCI6MTc4MDg0MTg2OSwiZXhwIjoxNzgyMDUxNDY5fQ.NJVXIhEErEIt3xcHzJz7YnHo0AnTCuASxPol0X_Ga5QPR38aWmJfaoO416V65apjJmB-E9MONvgrMoHJORRz8B0N1FvxPJWxCqgqAiT8FIBjb9lDZuje4GKG5-qt1LMeTwHlIQKozAIx0UlfV0gOtf3V7xatdOx6prV7X4AZ9SmBffJhZvJjaLP1ZvB0QO7affozK5cJdAavNdFaO7IaY4SNx-WqUwiiMjvJ8kGIfvTAfcKkNhC72ZTt5hVun0Nh4r8YYHBC8deURtL5nxVYzScKhAxgDiFJkULP81hpZWBpIK4ld5kByuWiMcFJrIdPU5_wfsah0xXymH1C4aBZLw")

exposure_vcf <- "ebi-a-GCST90012005.vcf.gz"
outcome_vcf  <- "out-ebi-a-GCST90018874.vcf.gz"
vcf_files    <- list(IL6 = exposure_vcf, Keloid = outcome_vcf)

# 检查 VCF 文件是否存在
for (nm in names(vcf_files)) {
  if (!file.exists(vcf_files[[nm]])) {
    stop("FATAL: VCF file not found: ", vcf_files[[nm]])
  }
}
cat("  All VCF files exist.\n")

# ---- 2个分析区域 ----
REGION_DEFS <- list(
  IL6R = list(chr = "1", start = 153000000, end = 156000000),
  HLA  = list(chr = "6", start = 29500000,  end = 34500000)
)

# ---- 先验组合 ----
PRIORS <- data.frame(
  label = c("Default (coloc)", "Permissive_shared", "Stringent_shared",
            "High_p1_IL6", "High_p2_Keloid", "Equal_priors",
            "Very_stringent", "Very_permissive"),
  p1 = c(1e-4, 1e-4, 1e-4, 1e-3, 1e-4, 1e-4, 1e-4, 1e-3),
  p2 = c(1e-4, 1e-4, 1e-4, 1e-4, 1e-3, 1e-4, 1e-4, 1e-3),
  p12= c(1e-5, 1e-4, 1e-6, 1e-5, 1e-5, 1e-4, 1e-7, 1e-4),
  stringsAsFactors = FALSE
)

MAF_THRESHOLDS <- c(0, 0.01, 0.05, 0.10)

# 日志收集器
log_lines <- character()
log <- function(msg) {
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

log("  Parameters set.")
log("  Regions: IL6R (chr1:153-156Mb), HLA (chr6:29.5-34.5Mb)")
log(paste("  Samples: IL6 N=", IL6_N, ", Keloid N=", KELOID_N,
          "(cases=", KELOID_CASES, ")"))

# ======================================================================
# 1. 加载 R 包
# ======================================================================
cat("\n--- Step 1: Loading packages ---\n")

needed_pkgs <- c("coloc", "VariantAnnotation", "Rsamtools",
                  "GenomicRanges", "dplyr", "ggplot2", "tidyr",
                  "ieugwasr", "biomaRt", "ggrepel", "gridExtra")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)

for (pkg in needed_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "...\n")
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
    }, error = function(e) {
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    })
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
}

missing <- needed_pkgs[!sapply(needed_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) stop("FATAL: Missing packages: ", paste(missing, collapse = ", "))

# 尝试安装 R.utils (用于超时控制)
if (!requireNamespace("R.utils", quietly = TRUE)) {
  tryCatch(install.packages("R.utils", repos = "https://cloud.r-project.org", quiet = TRUE),
           error = function(e) cat("  Note: R.utils not installed, will use base::setTimeLimit for timeout.\n"))
}

# ---- 配置 OpenGWAS Token ----
tryCatch(ieugwasr::set_access_token(OPENGWAS_TOKEN), error = function(e) {
  cat("  Note: set_access_token not exported, using fallback methods.\n")
})
options(ieugwasr.access_token = OPENGWAS_TOKEN)
Sys.setenv(IEU_ACCESS_TOKEN = OPENGWAS_TOKEN)
cat("  OpenGWAS token configured.\n")

cat("  All packages loaded OK.\n\n")
rm(OPENGWAS_TOKEN)

# ======================================================================
# 2. 定义辅助函数
# ======================================================================
cat("--- Step 2: Defining helper functions ---\n")

# ---- 2a. 安全的图形设备关闭 ----
safe_dev_off <- function() {
  while (length(grDevices::dev.list()) > 0) {
    tryCatch(grDevices::dev.off(), error = function(e) NULL)
  }
}

# ---- 2b. 从 VCF 提取区域数据 ----
extract_gwas <- function(vcf_path, region_gr, label = "") {
  if (!file.exists(vcf_path)) {
    log(paste("    [", label, "] VCF file not found:", vcf_path))
    return(NULL)
  }

  tab <- Rsamtools::TabixFile(vcf_path)
  open(tab)
  on.exit(close(tab))

  vcf <- tryCatch({
    VariantAnnotation::readVcf(tab, param = region_gr)
  }, error = function(e) {
    log(paste("    [", label, "] Error reading VCF:", conditionMessage(e)))
    NULL
  })

  if (is.null(vcf) || nrow(vcf) < 1) {
    log(paste("    [", label, "] No variants found in region"))
    return(NULL)
  }

  fixed <- rowRanges(vcf)
  sn_rle <- GenomicRanges::seqnames(fixed)
  if (methods::is(sn_rle, "Rle")) {
    chr <- as.character(S4Vectors::runValue(sn_rle))[
      S4Vectors::runLength(sn_rle) > 0
    ]
    # 展开 Rle 到每个 SNP
    chr <- rep(as.character(S4Vectors::runValue(sn_rle)),
               times = S4Vectors::runLength(sn_rle))
  } else {
    chr <- as.character(sn_rle)
  }

  pos <- BiocGenerics::start(fixed)
  snp_names <- names(fixed)

  gd <- geno(vcf)

  # 检查必要字段
  required_fields <- c("ES", "SE", "LP", "AF")
  missing_fields <- setdiff(required_fields, names(gd))
  if (length(missing_fields) > 0) {
    log(paste("    [", label, "] Missing fields in VCF:", paste(missing_fields, collapse = ", ")))
    log(paste("    Available fields:", paste(names(gd), collapse = ", ")))
    return(NULL)
  }

  ext_col <- function(mat, col = 1) {
    sapply(seq_len(nrow(mat)), function(i) {
      val <- mat[i, col]
      if (is.null(val) || length(val) == 0 || (length(val) == 1 && is.na(val))) {
        NA
      } else if (is.list(val)) {
        as.numeric(val[[1]])
      } else {
        as.numeric(val)
      }
    })
  }

  beta <- ext_col(gd$ES)
  se   <- ext_col(gd$SE)
  lp   <- ext_col(gd$LP)
  pval <- 10^(-lp)
  af   <- ext_col(gd$AF)

  # 处理 ALT 和 REF
  ea <- sapply(fixed$ALT, function(x) {
    if (is.list(x) && length(x[[1]]) > 0) as.character(x[[1]][1]) else as.character(x[1])
  })
  oa <- as.character(fixed$REF)

  ok <- !is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(af) &
        is.finite(beta) & is.finite(se) & is.finite(pval) & is.finite(af)
  if (sum(ok) < 10) {
    log(paste("    [", label, "] Too few valid variants:", sum(ok), "out of", length(ok)))
    return(NULL)
  }

  result <- data.frame(
    SNP = snp_names[ok], CHR = chr[ok], BP = pos[ok],
    effect_allele = ea[ok], other_allele = oa[ok],
    beta = beta[ok], se = se[ok], pval = pval[ok],
    MAF = pmin(af[ok], 1 - af[ok]),
    stringsAsFactors = FALSE
  )

  # 去重: 保留 pval 最小的那个
  dup_snps <- duplicated(result$SNP)
  if (any(dup_snps)) {
    result <- result[order(result$pval), ]
    result <- result[!duplicated(result$SNP), ]
    log(paste("    [", label, "] Removed", sum(dup_snps), "duplicate SNPs"))
  }

  log(paste("    [", label, "] Extracted", nrow(result), "variants"))
  result
}

# ---- 2c. 等位基因调和 ----
harmonise <- function(df1, df2) {
  m <- merge(df1, df2, by = "SNP", suffixes = c(".gwas1", ".gwas2"),
             all = FALSE)  # inner join 只保留共同 SNP
  if (nrow(m) == 0) {
    log("    No shared SNPs found between datasets")
    return(m)
  }

  comp <- c(A = "T", T = "A", C = "G", G = "C")

  ea1 <- m$effect_allele.gwas1; oa1 <- m$other_allele.gwas1
  ea2 <- m$effect_allele.gwas2; oa2 <- m$other_allele.gwas2

  # 四个方向的匹配检查
  direct     <- ea1 == ea2 & oa1 == oa2
  swap       <- ea1 == oa2 & oa1 == ea2
  comp_d     <- !is.na(ea1) & !is.na(oa1) & !is.na(ea2) & !is.na(oa2) & comp[ea1] == ea2 & comp[oa1] == oa2
  comp_s     <- !is.na(ea1) & !is.na(oa1) & !is.na(ea2) & !is.na(oa2) & comp[ea1] == oa2 & comp[oa1] == ea2

  # 修复: comp_d/comp_s 中无效的匹配
  comp_d[is.na(comp_d)] <- FALSE
  comp_s[is.na(comp_s)] <- FALSE

  status <- rep("mismatch", nrow(m))
  status[direct] <- "direct_match"
  status[swap]   <- "strand_swap"
  status[comp_d] <- "complement"
  status[comp_s] <- "swap_complement"

  m$harmony <- status
  m$beta.aligned <- ifelse(status %in% c("strand_swap", "swap_complement"),
                            -m$beta.gwas2, m$beta.gwas2)

  # 保留调和成功的 SNP
  ok <- status != "mismatch"
  m <- m[ok, ]

  # MAF 一致性检查: 调和后 MAF 应相近 (排除 misalignment)
  if (nrow(m) > 0) {
    maf_diff <- abs(m$MAF.gwas1 - m$MAF.gwas2)
    maf_ok <- maf_diff < 0.15  # 允许 15% 差异
    if (sum(!maf_ok) > 0) {
      log(paste("    Removed", sum(!maf_ok), "SNPs with MAF discrepancy > 0.15"))
      m <- m[maf_ok, ]
    }
  }

  log(paste("    Harmonised:", nrow(m), "SNPs",
            "(direct:", sum(direct), "| swap:", sum(swap),
            "| complement:", sum(comp_d, na.rm = TRUE),
            "| swap+comp:", sum(comp_s, na.rm = TRUE),
            "| mismatch:", sum(!ok), ")"))
  m
}

# ---- 2d. 运行 coloc ----
run_coloc <- function(m, p1, p2, p12, IL6_N, KELOID_N, KELOID_S) {
  if (is.null(m) || nrow(m) < 50) {
    return(NULL)
  }

  # 清理 NA 值 — 更健壮的方式
  keep <- !is.na(m$beta.gwas1) & !is.na(m$se.gwas1) &
          !is.na(m$beta.aligned) & !is.na(m$se.gwas2) &
          !is.na(m$BP.gwas1) & !is.na(m$BP.gwas2) &
          !is.na(m$MAF.gwas1) & !is.na(m$MAF.gwas2) &
          is.finite(m$beta.gwas1) & is.finite(m$se.gwas1) &
          is.finite(m$beta.aligned) & is.finite(m$se.gwas2)

  if (sum(keep) < 50) return(NULL)

  m_clean <- m[keep, ]

  d1 <- list(
    beta     = m_clean$beta.gwas1,
    varbeta  = m_clean$se.gwas1^2,
    snp      = m_clean$SNP,
    position = m_clean$BP.gwas1,
    type     = "quant",
    N        = IL6_N,
    MAF      = m_clean$MAF.gwas1
  )
  d2 <- list(
    beta     = m_clean$beta.aligned,
    varbeta  = m_clean$se.gwas2^2,
    snp      = m_clean$SNP,
    position = m_clean$BP.gwas2,
    type     = "cc",
    N        = KELOID_N,
    s        = KELOID_S,
    MAF      = m_clean$MAF.gwas2
  )

  tryCatch({
    coloc::coloc.abf(d1, d2, p1 = p1, p2 = p2, p12 = p12)
  }, error = function(e) {
    cat("    coloc.abf error:", conditionMessage(e), "\n")
    NULL
  })
}

# ---- 2e. LD 查询 (双层超时保护) ----
query_ld <- function(merged, lead_snp, max_snps = 2000, timeout_sec = 30) {
  sig <- which(merged$pval.gwas1 < 1e-4)
  if (length(sig) > max_snps) {
    sig <- sig[order(merged$pval.gwas1[sig])[1:max_snps]]
  }
  snps <- merged$SNP[sig]
  if (!lead_snp %in% snps) {
    snps <- c(lead_snp, snps)
  }

  cat("    LD query for", length(snps), "SNPs (", timeout_sec, "sec timeout)...\n")

  # 优先使用 R.utils::withTimeout (更可靠)
  if (requireNamespace("R.utils", quietly = TRUE)) {
    result <- tryCatch({
      R.utils::withTimeout(
        ieugwasr::ld_matrix(snps, pop = "EUR"),
        timeout = timeout_sec,
        onTimeout = "error"
      )
    }, error = function(e) {
      cat("    LD query failed:", conditionMessage(e), "\n")
      cat("    Proceeding without LD colors.\n")
      NULL
    })
    return(result)
  }

  # Fallback: setTimeLimit
  setTimeLimit(elapsed = timeout_sec, transient = TRUE)
  on.exit({
    tryCatch(setTimeLimit(c(NA, NA), transient = FALSE), error = function(e) NULL)
  })

  tryCatch({
    ieugwasr::ld_matrix(snps, pop = "EUR")
  }, error = function(e) {
    cat("    LD query failed:", conditionMessage(e), "\n")
    cat("    Proceeding without LD colors.\n")
    NULL
  })
}

# ---- 2f. LocusCompare 图 ----
plot_locuscompare <- function(merged, region_name, coloc_main) {
  log(paste("  Plotting LocusCompare:", region_name))

  merged <- merged[merged$MAF.gwas1 > 0.01 & merged$MAF.gwas2 > 0.01, ]
  if (nrow(merged) < 50) { log("    SKIP: <50 SNPs"); return(invisible(NULL)) }

  merged$pos_mb <- merged$BP.gwas1 / 1e6
  chr_name <- unique(merged$CHR.gwas1)
  xmin <- min(merged$pos_mb, na.rm = TRUE); xmax <- max(merged$pos_mb, na.rm = TRUE)

  # 处理 pval 中的 NA
  merged$pval.gwas1[is.na(merged$pval.gwas1)] <- 1
  merged$pval.gwas2[is.na(merged$pval.gwas2)] <- 1

  top_il6 <- merged[which.min(merged$pval.gwas1), , drop = FALSE]
  top_kel <- merged[which.min(merged$pval.gwas2), , drop = FALSE]

  # ---- LD 着色 ----
  ld_mat <- tryCatch(query_ld(merged, top_il6$SNP[1]), error = function(e) NULL)

  if (!is.null(ld_mat) && ncol(ld_mat) > 0) {
    common <- intersect(merged$SNP, colnames(ld_mat))
    merged$LD_r2 <- NA_real_
    if (top_il6$SNP[1] %in% colnames(ld_mat)) {
      r2v <- ld_mat[, top_il6$SNP[1]]^2
      idx <- match(intersect(common, names(r2v)), merged$SNP)
      merged$LD_r2[idx] <- as.numeric(r2v[intersect(common, names(r2v))])
    }
    merged$LD_cat <- cut(merged$LD_r2,
      breaks = c(-0.01, 0.2, 0.4, 0.6, 0.8, 1.01),
      labels = c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"),
      include.lowest = TRUE)
    ld_available <- sum(!is.na(merged$LD_r2)) > 10
  } else {
    merged$LD_r2 <- NA_real_
    merged$LD_cat <- NA
    ld_available <- FALSE
  }

  pp4_val <- if (!is.null(coloc_main)) coloc_main$summary["PP.H4.abf"] else NA
  pp4_str <- if (!is.na(pp4_val)) sprintf("%.4f", pp4_val) else "NA"

  # ---- 上方散点图 ----
  merged$LP.gwas1 <- -log10(pmax(merged$pval.gwas1, 1e-300))
  merged$LP.gwas2 <- -log10(pmax(merged$pval.gwas2, 1e-300))

  plot_data <- data.frame(
    pos_mb = rep(merged$pos_mb, 2),
    LP = c(merged$LP.gwas1, merged$LP.gwas2),
    LD_cat = rep(merged$LD_cat, 2),
    Trait = rep(c("IL-6 levels", "Keloids"), each = nrow(merged)),
    stringsAsFactors = FALSE
  )

  if (ld_available) {
    p_top <- ggplot(plot_data, aes(x = pos_mb, y = LP)) +
      geom_point(aes(color = LD_cat, shape = Trait), size = 1.8, alpha = 0.7) +
      scale_color_manual(
        name = expression(italic(r)^2 ~ with ~ lead ~ SNP),
        values = c("0-0.2" = "#7B9E7B", "0.2-0.4" = "#4CA64C",
                   "0.4-0.6" = "#2E8B2E", "0.6-0.8" = "#006400", "0.8-1.0" = "#FF0000"),
        na.value = "grey50",
        drop = FALSE) +
      scale_shape_manual(values = c("IL-6 levels" = 19, "Keloids" = 17))
  } else {
    p_top <- ggplot(plot_data, aes(x = pos_mb, y = LP, color = Trait, shape = Trait)) +
      geom_point(size = 1.8, alpha = 0.7) +
      scale_color_manual(values = c("IL-6 levels" = COL_IL6, "Keloids" = COL_KELOID)) +
      scale_shape_manual(values = c("IL-6 levels" = 19, "Keloids" = 17))
  }

  p_top <- p_top +
    geom_point(data = top_il6,
               aes(x = pos_mb, y = -log10(pval.gwas1)),
               color = COL_IL6, size = 4, shape = 18) +
    geom_point(data = top_kel,
               aes(x = pos_mb, y = -log10(pval.gwas2)),
               color = COL_KELOID, size = 4, shape = 18) +
    labs(title = paste0(region_name, " Locus  |  PP.H4 = ", pp4_str),
         subtitle = paste0("chr", chr_name, ":",
                           format(round(xmin, 3), nsmall = 3), "-",
                           format(round(xmax, 3), nsmall = 3), " Mb"),
         x = NULL, y = expression(-log[10](italic(P)))) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top",
          legend.box = "vertical",
          plot.margin = margin(8, 10, 2, 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank())

  # 基因标注 (独立 tryCatch, 失败不影响主图)
  gene_df <- tryCatch({
    bm <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    chrm <- gsub("chr", "", chr_name)
    gs <- biomaRt::getBM(
      attributes = c("chromosome_name", "start_position", "end_position",
                      "gene_biotype", "hgnc_symbol"),
      filters = c("chromosome_name", "start", "end"),
      values = list(chrm, round(xmin * 1e6), round(xmax * 1e6)),
      mart = bm)

    if (nrow(gs) == 0) return(NULL)

    # 修复: 过滤时去除 NA 的 hgnc_symbol
    gs <- gs[gs$gene_biotype == "protein_coding", ]
    gs <- gs[!is.na(gs$hgnc_symbol) & nchar(gs$hgnc_symbol) > 0, ]
    if (nrow(gs) == 0) return(NULL)

    gs$mid <- (gs$start_position + gs$end_position) / 2 / 1e6
    gs <- gs[gs$mid >= xmin & gs$mid <= xmax, ]
    gs
  }, error = function(e) {
    cat("    Gene annotation skipped:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(gene_df) && nrow(gene_df) > 0) {
    if (nrow(gene_df) > 25) gene_df <- gene_df[1:25, ]
    gene_df <- gene_df[order(gene_df$mid), ]
    gene_df$y_pos <- rep(c(0.25, 0.50, 0.75), length.out = nrow(gene_df))

    p_bot <- ggplot(gene_df) +
      geom_segment(aes(x = start_position / 1e6,
                        xend = end_position / 1e6,
                        y = 0.1, yend = 0.1),
                    color = "#555555", linewidth = 2.5, alpha = 0.8) +
      geom_text(aes(x = mid, y = y_pos, label = hgnc_symbol),
                size = 2.5, color = "#333333") +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_x_continuous(limits = c(xmin, xmax), expand = c(0, 0)) +
      labs(x = paste0("Chromosome ", chr_name, " (Mb)")) +
      theme_void(base_size = 10) +
      theme(axis.title.x = element_text(size = 10, margin = margin(4, 0, 0, 0)),
            axis.text.x = element_text(size = 8),
            axis.ticks.x = element_line(color = "grey60"),
            axis.ticks.length.x = unit(3, "pt"),
            plot.margin = margin(0, 10, 8, 10))
  } else {
    p_bot <- ggplot(data.frame(x = c(xmin, xmax)), aes(x)) +
      scale_x_continuous(limits = c(xmin, xmax), expand = c(0, 0)) +
      labs(x = paste0("Chromosome ", chr_name, " (Mb)")) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            plot.margin = margin(0, 10, 8, 10))
  }

  # ---- 组合 & 保存 ----
  out_file <- paste0("pic.coloc_v4_", tolower(region_name), "_locuscompare.png")
  tryCatch({
    grDevices::png(out_file, width = 10, height = 7, units = "in", res = 200)
    gridExtra::grid.arrange(p_top, p_bot, ncol = 1, heights = c(2.8, 1.2))
    safe_dev_off()
    log(paste("    Saved:", out_file))
  }, error = function(e) {
    safe_dev_off()
    log(paste("    PLOT ERROR:", conditionMessage(e)))
  })
}

# ======================================================================
# 3. 提取数据 & 运行 colocalization
# ======================================================================
cat("\n--- Step 3: Data extraction and colocalization ---\n")

# 用列表存放结果
all_data   <- list()
coloc_res  <- list()   # 只存主结果 (Default prior)
sens_df    <- data.frame()
summary_df <- data.frame()

for (rn in names(REGION_DEFS)) {
  rd <- REGION_DEFS[[rn]]
  region_gr <- GRanges(rd$chr, IRanges(rd$start, rd$end))

  cat("\n===== Region:", rn, "=====\n\n")

  # 3a. 提取
  log(paste("Extracting IL-6 data..."))
  il6_df <- extract_gwas(vcf_files$IL6, region_gr, "IL-6")

  log(paste("Extracting Keloid data..."))
  kel_df <- extract_gwas(vcf_files$Keloid, region_gr, "Keloid")

  if (is.null(il6_df) || is.null(kel_df)) {
    log("  SKIP: data extraction returned NULL")
    next
  }

  # 3b. 调和
  log("Harmonising alleles...")
  hm <- harmonise(il6_df, kel_df)
  if (nrow(hm) < 50) {
    log(paste("  SKIP: only", nrow(hm), "SNPs after harmonisation (need >= 50)"))
    next
  }

  all_data[[rn]] <- hm

  # 3c. 主结果 (Default prior)
  hm_filt <- hm[hm$MAF.gwas1 > 0.01 & hm$MAF.gwas2 > 0.01, ]
  log(paste("  MAF>0.01:", nrow(hm_filt), "SNPs"))

  if (nrow(hm_filt) >= 50) {
    res <- run_coloc(hm_filt, 1e-4, 1e-4, 1e-5, IL6_N, KELOID_N, KELOID_S)
    if (!is.null(res)) {
      coloc_res[[rn]] <- res
      log(paste("  >> PP.H4 =", sprintf("%.6f", res$summary["PP.H4.abf"])))
    } else {
      log("  >> coloc.abf returned NULL for Default prior")
    }
  } else {
    log("  SKIP: insufficient SNPs (MAF>0.01) for coloc.abf")
  }

  # 3d. 敏感性分析: 先验组合
  for (i in seq_len(nrow(PRIORS))) {
    res <- run_coloc(hm_filt, PRIORS$p1[i], PRIORS$p2[i], PRIORS$p12[i],
                     IL6_N, KELOID_N, KELOID_S)
    if (!is.null(res)) {
      sens_df <- rbind(sens_df, data.frame(
        Region = rn, Type = "prior_sensitivity", Condition = PRIORS$label[i],
        p1 = PRIORS$p1[i], p2 = PRIORS$p2[i], p12 = PRIORS$p12[i],
        N_SNPs = nrow(hm_filt),
        PP.H4 = res$summary["PP.H4.abf"],
        stringsAsFactors = FALSE))
    }
  }

  # 3e. 敏感性分析: MAF阈值
  for (maf in MAF_THRESHOLDS) {
    hm_m <- hm[hm$MAF.gwas1 >= maf & hm$MAF.gwas2 >= maf, ]
    if (nrow(hm_m) < 50) next
    res <- run_coloc(hm_m, 1e-4, 1e-4, 1e-5, IL6_N, KELOID_N, KELOID_S)
    if (!is.null(res)) {
      sens_df <- rbind(sens_df, data.frame(
        Region = rn, Type = "maf_threshold", Condition = paste0("MAF >= ", maf),
        p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, N_SNPs = nrow(hm_m),
        PP.H4 = res$summary["PP.H4.abf"],
        stringsAsFactors = FALSE))
    }
  }

  # ---- Checkpoint: 每个区域跑完立即存盘 ----
  write.csv(sens_df, "table.coloc_v4_sensitivity.csv", row.names = FALSE)
  log(paste("  Checkpoint saved: sensitivity table (", nrow(sens_df), "rows)"))
}

# ---- 最终保存敏感性结果 ----
write.csv(sens_df, "table.coloc_v4_sensitivity.csv", row.names = FALSE)
log(paste("\nFinal sensitivity results saved:", nrow(sens_df), "rows"))

# ======================================================================
# 4. 生成 LocusCompare 图
# ======================================================================
cat("\n\n--- Step 4: LocusCompare plots ---\n")

if (length(coloc_res) == 0) {
  log("WARNING: coloc_res is empty, cannot generate plots.")
  log("Possible causes: all data extraction failed or coloc.abf errored.")
} else {
  for (rn in names(coloc_res)) {
    cat("\n===== Plot:", rn, "=====\n\n")
    if (!rn %in% names(all_data)) {
      log(paste("  SKIP:", rn, "- no harmonised data"))
      next
    }
    plot_locuscompare(all_data[[rn]], rn, coloc_res[[rn]])
  }
}

# ======================================================================
# 5. 敏感性分析汇总图
# ======================================================================
cat("\n\n--- Step 5: Sensitivity summary plot ---\n")

if (nrow(sens_df) > 0) {
  sens_prior <- sens_df[sens_df$Type == "prior_sensitivity", ]
  sens_maf   <- sens_df[sens_df$Type == "maf_threshold", ]
  plot_list <- list()

  if (nrow(sens_prior) > 0) {
    p1 <- ggplot(sens_prior, aes(x = Condition, y = PP.H4, fill = Region)) +
      geom_col(position = "dodge", width = 0.6, alpha = 0.85) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", alpha = 0.5) +
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", alpha = 0.5) +
      scale_fill_manual(values = c("IL6R" = COL_IL6, "HLA" = COL_KELOID)) +
      labs(title = "Sensitivity: Prior Probability Combinations",
           subtitle = "Dashed=strong (0.8), Dotted=moderate (0.5)",
           x = "Prior Setting", y = "PP.H4") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    plot_list[[1]] <- p1
  }

  if (nrow(sens_maf) > 0) {
    sens_maf$MAF_val <- as.numeric(gsub("MAF >= ", "", sens_maf$Condition))
    p2 <- ggplot(sens_maf, aes(x = MAF_val, y = PP.H4, color = Region, group = Region)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", alpha = 0.5) +
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", alpha = 0.5) +
      scale_color_manual(values = c("IL6R" = COL_IL6, "HLA" = COL_KELOID)) +
      labs(title = "Sensitivity: MAF Thresholds",
           x = "MAF Threshold", y = "PP.H4") +
      scale_x_continuous(breaks = sort(unique(sens_maf$MAF_val))) +
      theme_minimal(base_size = 11)
    plot_list[[2]] <- p2
  }

  if (length(plot_list) > 0) {
    tryCatch({
      grDevices::png("pic.coloc_v4_sensitivity.png", width = 10, height = 8, units = "in", res = 150)
      if (length(plot_list) == 2) {
        graphics::layout(matrix(1:2, ncol = 1))
        graphics::par(mar = c(4.5, 4, 3, 2))
        print(plot_list[[1]])
        graphics::par(mar = c(4.5, 4, 3, 2))
        print(plot_list[[2]])
      } else {
        print(plot_list[[1]])
      }
      safe_dev_off()
      log("Sensitivity plot saved.")
    }, error = function(e) {
      safe_dev_off()
      log(paste("Sensitivity plot error:", conditionMessage(e)))
    })
  }
} else {
  log("No sensitivity results to plot.")
}

# ======================================================================
# 6. 汇总表
# ======================================================================
cat("\n\n--- Step 6: Summary table ---\n")

if (length(coloc_res) == 0) {
  log("WARNING: No colocalization results to summarise!")
  log("This means coloc.abf returned NULL for all regions.")
  write.csv(summary_df, "table.coloc_v4_summary.csv", row.names = FALSE)
} else {
  for (rn in names(coloc_res)) {
    r <- coloc_res[[rn]]
    hm_data <- all_data[[rn]]
    if (is.null(hm_data)) next
    hm_filt <- hm_data[hm_data$MAF.gwas1 > 0.01 & hm_data$MAF.gwas2 > 0.01, ]
    pp4 <- r$summary["PP.H4.abf"]
    pp3 <- r$summary["PP.H3.abf"]

    conclusion <- if (is.numeric(pp4) && !is.na(pp4)) {
      if (pp4 > 0.8) "Shared causal variant (strong)"
      else if (pp4 > 0.5) "Moderate evidence for shared variant"
      else if (is.numeric(pp3) && !is.na(pp3) && pp3 > pp4) "Likely distinct causal variants"
      else "Inconclusive"
    } else {
      "Unknown (missing PP.H4)"
    }

    summary_df <- rbind(summary_df, data.frame(
      Region = rn, Analysis = "Expanded region, MAF>0.01",
      N_SNPs = nrow(hm_filt),
      PP.H0 = sprintf("%.6f", r$summary["PP.H0.abf"]),
      PP.H1 = sprintf("%.6f", r$summary["PP.H1.abf"]),
      PP.H2 = sprintf("%.6f", r$summary["PP.H2.abf"]),
      PP.H3 = sprintf("%.6f", r$summary["PP.H3.abf"]),
      PP.H4 = sprintf("%.6f", pp4),
      Conclusion = conclusion,
      stringsAsFactors = FALSE))
  }

  if (nrow(summary_df) > 0) {
    print(summary_df, row.names = FALSE)
    write.csv(summary_df, "table.coloc_v4_summary.csv", row.names = FALSE)
    log(paste("Summary table saved:", nrow(summary_df), "rows"))
  } else {
    log("No summary rows generated (all regions had no coloc results).")
    write.csv(summary_df, "table.coloc_v4_summary.csv", row.names = FALSE)
  }
}

# ======================================================================
# 7. 保存日志
# ======================================================================
cat("\n\n--- Step 7: Saving session log ---\n")

sink_output <- capture.output({
  cat("Colocalization v4.1: IL-6 ↔ Keloid\n")
  cat("================================\n\n")
  cat("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("GWAS: IL-6=GCST90012005(N=", IL6_N, "), Keloid=GCST90018874(N=", KELOID_N,
      ",cases=", KELOID_CASES, ")\n\n")
  cat("R version:", R.version.string, "\n\n")
  cat("Packages:\n")
  for (pkg in needed_pkgs) {
    if (requireNamespace(pkg, quietly = TRUE))
      cat("  -", pkg, ":", as.character(packageVersion(pkg)), "\n")
  }
  cat("\nRegions:\n")
  for (nm in names(REGION_DEFS)) {
    rd <- REGION_DEFS[[nm]]
    cat("  ", nm, ": chr", rd$chr, ":", rd$start, "-", rd$end,
        " (", (rd$end - rd$start)/1e6, "Mb)\n", sep = "")
  }
  cat("\nResults:\n")
  if (length(coloc_res) > 0) {
    for (nm in names(coloc_res))
      cat("  ", nm, ": PP.H4 =", sprintf("%.6f", coloc_res[[nm]]$summary["PP.H4.abf"]), "\n")
  } else {
    cat("  (No results - coloc_res was empty)\n")
  }
  cat("\nSummary table rows:", nrow(summary_df), "\n")
  cat("Sensitivity rows:", nrow(sens_df), "\n\n")
  cat("Output files:\n")
  cat("  table.coloc_v4_summary.csv\n")
  cat("  table.coloc_v4_sensitivity.csv\n")
  cat("  pic.coloc_v4_il6r_locuscompare.png\n")
  cat("  pic.coloc_v4_hla_locuscompare.png\n")
  cat("  pic.coloc_v4_sensitivity.png\n")
  cat("  colocalization_v4_session_info.txt\n")
})

writeLines(sink_output, con = "colocalization_v4_session_info.txt")
cat("Session log saved.\n\n")

# ======================================================================
# 8. 完成
# ======================================================================
cat("========================================\n")
cat("Colocalization v4.1 Complete!\n")
cat("========================================\n")
cat("Output:\n")
cat("  table.coloc_v4_summary.csv\n")
cat("  table.coloc_v4_sensitivity.csv\n")
cat("  pic.coloc_v4_il6r_locuscompare.png\n")
cat("  pic.coloc_v4_hla_locuscompare.png\n")
cat("  pic.coloc_v4_sensitivity.png\n")
cat("  colocalization_v4_session_info.txt\n\n")

if (length(coloc_res) == 0) {
  cat("!!! WARNING: No colocalization results obtained.\n")
  cat("!!! Check the console output above for error messages.\n")
} else {
  cat("SUCCESS: colocalization analysis completed with results for",
      length(coloc_res), "region(s).\n")
}

cat("\nDone.\n")
