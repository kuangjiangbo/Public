###############################################################################
##  IL-6 x Keloid  Mendelian randomization
##  外部重复分析 —— Greene et al. 2025 (Nat Commun 2025;16:7770, PMID 40835838)
##  结局 GWAS: GCST90652487  multi-ancestry keloid meta-analysis
##             4,086 欧洲 + 1,055 东亚 + 2,696 非裔 = 7,837 cases / 1,593,009 controls
##             GRCh37 | 1-based | GWAS-SSF v1.0 | 17,748,456 SNP | MAF >= 0.01
##
##  原稿（结局 = Sakaue GCST90018874，668 例）得到提示性结果
##  IVW OR 1.214 (P = 0.033)，且分层显示信号由 HLA 簇驱动、cis-IL6R 簇阴性。
##  本脚本用 12 倍病例数的扩展数据集重复同一套分析。
##
##  分析集合：
##        All 23 SNPs / IL6R only (cis) / HLA only
##        All minus rs2228145 / IL6R minus rs2228145
##  每个集合：IVW(固定/随机) / MR-Egger / weighted median / weighted mode /
##            simple mode + Cochran's Q + I2 + MR-Egger 截距 + MR-PRESSO
##            + 簇内 leave-one-out
##  另输出 Sakaue vs Greene 并排对照表与近似功效计算。
##
##  ── 2026-09-26 修订（v5）──
##  1) mr_presso() 调用此前误传 OUTCOME= 参数：MRPRESSO 1.0 的 mr_presso() 没有该参数
##     （那是更高版本的 API），会被当作 unused argument 直接报错，而 error=function(e) NULL
##     把错误吞掉了，导致复制期 MR-PRESSO 静默失败、既无输出也无文件。现已移除该参数、
##     加上固定 seed 与失败告警。
##  2) 新增复制期 Steiger 方向性检验（原先只做了发现期）。
##  3) 功效计算的 n_ctrl 由 1,585,173 改为 1,593,009：前者是按复制集 SSF 的 N 列（1,593,010）
##     反推的对照数，而该 N 与 Greene 2025 论文 Table 1 及 GWAS Catalog 官方汇总
##     （7,837 例 / 1,593,009 对照，合计 1,600,846）不一致；以论文口径为准。
##
##  运行： RStudio 里 Ctrl+Shift+Enter，或
##         Rscript IL6_Keloid_Greene2025_replication.R
##
##  注意：本机 TwoSampleMR 版本为 0.7.x，API 与旧版（0.5/0.6）不同：
##        - mr_ivw() 没有 model= 参数（内部即随机效应，用 sigma 做过度离散校正）
##        - mr_weighted_median() 的 bootstrap 次数由 parameters$nboot 控制，
##          没有 nbootstrap= 参数
##        本脚本已按 0.7.x 适配，并自行实现固定/随机效应 IVW 以便对照。
##
##  首次运行需分块读 340 MB 大文件（约 2-3 分钟）；结果会缓存为
##  table.Greene_extracted_records.csv，再次运行自动跳过读取。
###############################################################################

t_start <- Sys.time()

## ===========================================================================
## 0. 路径配置（目录变动只需改这一段）
## ===========================================================================
BASE_DIR    <- "F:/MR/Interleukin-6 levels and keloids"
GREENE_FILE <- file.path(BASE_DIR, "数据/Greene2025/GCST90652487.tsv.gz")
SAKAUE_FILE <- file.path(BASE_DIR, "数据/第四次/harmonised_data.csv")
OUT_DIR     <- file.path(BASE_DIR, "数据/Greene2025")
CACHE_FILE  <- file.path(OUT_DIR, "table.Greene_extracted_records.csv")

GREENE_MD5 <- "8db7428b31874b9f9a509a4e0e6ff69c"
CHUNK_ROWS <- 2e6
WINDOW_PAD <- 5000L
LEAD_SNP   <- "rs2228145"
ALPHA      <- 0.05
PRESSO_SEED   <- 20260926   # MR-PRESSO 置换检验固定种子（保证可复现）
GREENE_NCASE  <- 7837L      # Greene 2025 论文 Table 1 汇总
GREENE_NCTRL  <- 1593009L   # 注意：SSF 的 N 列（1,593,010）与论文 Table 1 不一致，以论文为准

setwd(OUT_DIR)
stopifnot(file.exists(GREENE_FILE), file.exists(SAKAUE_FILE))

## ===========================================================================
## 1. 依赖包
## ===========================================================================
need <- c("data.table", "TwoSampleMR", "MRPRESSO", "ggplot2")
for (p in need) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("缺少 R 包 '%s'，请先 install.packages('%s')", p, p))
suppressMessages({ library(data.table); library(TwoSampleMR); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a
cat("R", as.character(getRversion()),
    "| TwoSampleMR", as.character(packageVersion("TwoSampleMR")), "\n\n")

## ===========================================================================
## 2. 校验 Greene 文件
## ===========================================================================
cat("== 2. 校验文件 ==\n")
cat("  文件:", basename(GREENE_FILE), "|",
    round(file.info(GREENE_FILE)$size / 1e6, 1), "MB\n")
md5_got <- unname(tools::md5sum(GREENE_FILE))
cat("  md5 实测:", md5_got, "| 校验:",
    if (identical(tolower(md5_got), GREENE_MD5)) "通过 OK" else "不通过 FAIL", "\n\n")
if (!identical(tolower(md5_got), GREENE_MD5))
  warning("md5 与官方值不一致，文件可能不完整：", GREENE_MD5)

GREENE_COLS <- c("chromosome", "base_pair_location", "effect_allele",
                 "other_allele", "beta", "standard_error",
                 "effect_allele_frequency", "p_value", "rsid")

## ===========================================================================
## 3. 暴露侧：23 个 IL-6 工具变量
## ===========================================================================
cat("== 3. 读入暴露侧（IL-6, Folkersen GCST90012005）==\n")
sak <- as.data.frame(fread(SAKAUE_FILE, showProgress = FALSE))

exposure_dat <- data.frame(
  SNP                    = sak$SNP,
  beta.exposure          = as.numeric(sak$beta.exposure),
  se.exposure            = as.numeric(sak$se.exposure),
  effect_allele.exposure = sak$effect_allele.exposure,
  other_allele.exposure  = sak$other_allele.exposure,
  eaf.exposure           = as.numeric(sak$eaf.exposure),
  pval.exposure          = as.numeric(sak$pval.exposure),
  exposure               = "IL-6", id.exposure = "IL6",
  stringsAsFactors       = FALSE)

snp_map <- data.table(SNP = sak$SNP,
                      chr = as.integer(sak$chr.exposure),
                      pos = as.integer(sak$pos.exposure),
                      locus = ifelse(as.integer(sak$chr.exposure) == 1L, "IL6R", "HLA"))

## 累计 R2 —— 与原稿口径一致（TwoSampleMR::add_rsq 算法）：rsq = b^2 / (b^2 + N*se^2)
N_exp  <- as.numeric(sak$samplesize.exposure[1])
rsq_i  <- exposure_dat$beta.exposure^2 /
          (exposure_dat$beta.exposure^2 + N_exp * exposure_dat$se.exposure^2)
R2_all <- sum(rsq_i)
R2_unc <- sum(2 * exposure_dat$eaf.exposure * (1 - exposure_dat$eaf.exposure) *
              exposure_dat$beta.exposure^2)
cat("  工具变量:", nrow(exposure_dat), "| IL6R:", sum(snp_map$locus == "IL6R"),
    "| HLA:", sum(snp_map$locus == "HLA"), "\n")
cat("  累计 R2:", round(R2_all, 4), "(原稿口径) | 未校正口径:", round(R2_unc, 4),
    "| 最小 F:", round(min((exposure_dat$beta.exposure / exposure_dat$se.exposure)^2), 2), "\n\n")

win <- snp_map[, .(lo = min(pos) - WINDOW_PAD, hi = max(pos) + WINDOW_PAD), by = chr]
setnames(win, "chr", "chromosome")
cat("  抽取窗口:\n"); print(win)

## ===========================================================================
## 4. 从 Greene 大文件中抽取目标位点（带缓存）
##    fread 的 skip 语义：skip = N 从第 N 个数据行开始（已实测确认）
## ===========================================================================
cat("\n== 4. 抽取目标位点 ==\n")
if (file.exists(CACHE_FILE)) {
  cat("  命中缓存，跳过 340 MB 重读:", basename(CACHE_FILE), "\n")
  greene_raw <- as.data.table(fread(CACHE_FILE, showProgress = FALSE))
} else {
  hit_list <- list(); chunk_i <- 0L; n_read <- 0L; t0 <- Sys.time()
  repeat {
    chunk_i <- chunk_i + 1L
    chunk <- tryCatch(
      suppressWarnings(fread(GREENE_FILE, skip = 1L + n_read, nrows = CHUNK_ROWS,
                             header = FALSE, col.names = GREENE_COLS,
                             showProgress = FALSE, data.table = TRUE)),
      error = function(e) NULL)
    if (is.null(chunk) || nrow(chunk) == 0L) break
    keep <- chunk[chromosome %in% win$chromosome]
    if (nrow(keep) > 0L) {
      keep <- merge(keep, win, by = "chromosome")
      keep <- keep[base_pair_location >= lo & base_pair_location <= hi]
      if (nrow(keep) > 0L)
        hit_list[[length(hit_list) + 1L]] <-
          keep[, .(chromosome, base_pair_location, effect_allele, other_allele,
                   beta, standard_error, effect_allele_frequency, p_value, rsid)]
    }
    n_read <- n_read + nrow(chunk)
    cat(sprintf("   已读 %9d 行 (%.1f 分钟) 命中 %d 条\r", n_read,
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                if (length(hit_list)) nrow(rbindlist(hit_list)) else 0L))
    if (nrow(chunk) < CHUNK_ROWS) break
  }
  cat("\n  共读取", format(n_read, big.mark = ","), "行\n")
  greene_raw <- if (length(hit_list)) unique(rbindlist(hit_list)) else data.table()
  setnames(greene_raw,
           c("chromosome", "base_pair_location", "effect_allele", "other_allele",
             "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"),
           c("chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome",
             "beta.outcome", "se.outcome", "eaf.outcome", "pval.outcome", "greene_rsid"))
  fwrite(greene_raw, CACHE_FILE)
  cat("  已缓存抽取结果:", basename(CACHE_FILE), "\n")
}
cat("  窗口内记录数:", nrow(greene_raw), "\n\n")

## 覆盖率核查 + 匹配（优先 rsid，缺失则回退到 chr:pos）
gr <- copy(greene_raw)
if (!"greene_rsid" %in% names(gr)) setnames(gr, "SNP", "greene_rsid", skip_absent = TRUE)
gr[, key_pos := paste(chr.outcome, pos.outcome)]
snp_map[, key_pos := paste(chr, pos)]

m1 <- merge(snp_map, gr, by.x = "SNP", by.y = "greene_rsid")
m1[, match_by := "rsid"]
rest <- snp_map[!SNP %in% m1$SNP]
m2 <- data.table()
if (nrow(rest) > 0L) {
  m2 <- merge(rest, gr, by = "key_pos", suffixes = c("", ".y"))
  if (nrow(m2) > 0L) { m2[, SNP := SNP]; m2[, match_by := "position"] }
}
cols <- c("SNP", "locus", "beta.outcome", "se.outcome", "effect_allele.outcome",
          "other_allele.outcome", "eaf.outcome", "pval.outcome", "match_by")
fin <- rbindlist(list(m1[, ..cols], m2[, ..cols]), use.names = TRUE, fill = TRUE)

cov <- data.table(SNP = snp_map$SNP, locus = snp_map$locus,
                  found_by_rsid = snp_map$SNP %in% gr$greene_rsid,
                  found_by_position = snp_map$key_pos %in% gr$key_pos)
cat("== 覆盖率核查 ==\n")
cat("  按 rsid 命中:", sum(cov$found_by_rsid), "/", nrow(cov),
    "| 按位置命中:", sum(cov$found_by_position), "/", nrow(cov), "\n")
if (any(!cov$found_by_rsid))
  cat("  按 rsid 未命中（已尝试按位置回退）:",
      paste(cov$SNP[!cov$found_by_rsid], collapse = ", "), "\n")
print(cov)
fwrite(cov, file.path(OUT_DIR, "table.Greene_coverage_audit.csv"))

greene_out <- as.data.frame(fin[, .(
  SNP,
  beta.outcome          = as.numeric(beta.outcome),
  se.outcome            = as.numeric(se.outcome),
  effect_allele.outcome = effect_allele.outcome,
  other_allele.outcome  = other_allele.outcome,
  eaf.outcome           = as.numeric(eaf.outcome),
  pval.outcome          = as.numeric(pval.outcome),
  outcome               = "Keloid (Greene 2025)",
  id.outcome            = "GCST90652487",
  samplesize.outcome    = 1593010L,
  locus = locus, match_by = match_by)])
cat("  进入 harmonise:", nrow(greene_out), "个位点（其中按位置回退",
    sum(greene_out$match_by == "position"), "个）\n\n")

## ===========================================================================
## 5. Harmonise
## ===========================================================================
cat("== 5. Harmonise ==\n")
harm <- harmonise_data(exposure_dat, greene_out, action = 2)
harm <- harm[harm$mr_keep, ]
if (!"locus" %in% names(harm)) harm <- merge(harm, snp_map[, .(SNP, locus)], by = "SNP")
dropped <- setdiff(exposure_dat$SNP, harm$SNP)
cat("  可分析位点:", nrow(harm),
    "| 丢失:", length(dropped),
    if (length(dropped)) paste0(" (", paste(dropped, collapse = ", "), ")") else "", "\n")
cat("  不可用原因:\n")
if (length(dropped)) print(unique(greene_out[greene_out$SNP %in% dropped,
                                             c("SNP", "locus", "effect_allele.outcome",
                                               "other_allele.outcome", "eaf.outcome")]))
fwrite(as.data.table(harm), file.path(OUT_DIR, "table.harmonised_Greene.csv"))

## ===========================================================================
## 6. 分析函数
## ===========================================================================
## 6.1 显式 IVW（固定/随机），与 TwoSampleMR 内部算法等价
ivw_fun <- function(be, bo, se_e, so, model = c("random", "fixed")) {
  model <- match.arg(model)
  w   <- 1 / (so / abs(be))^2                  # = 1/se_ratio^2
  biv <- sum(w * (bo / be)) / sum(w)
  se_f <- sqrt(1 / sum(w))
  df   <- length(be) - 1
  Q    <- sum(w * (bo / be - biv)^2)
  if (model == "random") {
    scale <- if (df > 0 && Q > df) sqrt(Q / df) else 1
    se <- se_f * scale
  } else se <- se_f
  list(b = biv, se = se, pval = 2 * pnorm(abs(biv / se), lower.tail = FALSE),
       Q = Q, Q_df = df, Q_pval = pchisq(Q, df, lower.tail = FALSE), nsnp = length(be))
}

## 6.2 一套方法跑一个集合（错误不再静默丢弃，写进 error 列）
mr_panel <- function(d, label) {
  be <- d$beta.exposure; bo <- d$beta.outcome
  se_e <- d$se.exposure; so <- d$se.outcome
  p <- default_parameters(); p$nboot <- 1000        # 1000 = TwoSampleMR 默认，够用且快很多
  calls <- list(
    "IVW (random)"    = function() ivw_fun(be, bo, se_e, so, "random"),
    "IVW (fixed)"     = function() ivw_fun(be, bo, se_e, so, "fixed"),
    "IVW (TwoSampleMR)" = function() mr_ivw(be, bo, se_e, so),
    "MR Egger"        = function() mr_egger_regression(be, bo, se_e, so),
    "Weighted median" = function() mr_weighted_median(be, bo, se_e, so, parameters = p),
    "Weighted mode"   = function() mr_weighted_mode(be, bo, se_e, so, parameters = p),
    "Simple mode"     = function() mr_simple_mode(be, bo, se_e, so, parameters = p))
  out <- lapply(names(calls), function(m) {
    err <- NA_character_
    r <- tryCatch(calls[[m]](), error = function(e) { err <<- conditionMessage(e); NULL })
    if (is.null(r) || is.na(r$b %||% NA))
      return(data.frame(Analysis = label, nSNP = length(be), Method = m,
                        b = NA, se = NA, pval = NA, OR = NA, OR_lci = NA, OR_uci = NA,
                        error = err %||% "returned NA"))
    data.frame(Analysis = label, nSNP = length(be), Method = m,
               b = r$b, se = r$se, pval = r$pval, OR = exp(r$b),
               OR_lci = exp(r$b - 1.96 * r$se), OR_uci = exp(r$b + 1.96 * r$se),
               error = NA_character_)
  })
  do.call(rbind, out)
}

loo_ivw <- function(d) {
  res <- lapply(seq_len(nrow(d)), function(i) {
    dd <- d[-i, , drop = FALSE]
    if (nrow(dd) < 3) return(NULL)
    r <- ivw_fun(dd$beta.exposure, dd$beta.outcome, dd$se.exposure, dd$se.outcome, "random")
    data.frame(dropped_SNP = d$SNP[i], nSNP = nrow(dd), b = r$b, se = r$se, pval = r$pval,
               OR = exp(r$b), OR_lci = exp(r$b - 1.96 * r$se), OR_uci = exp(r$b + 1.96 * r$se))
  })
  do.call(rbind, res)
}

## 注意：MRPRESSO 1.0 的 mr_presso() **没有 OUTCOME 参数**（更高版本才有）。
## 传入 OUTCOME= 会报 unused argument；旧版本这里用 error=function(e) NULL 把错误吞掉，
## 导致复制期 MR-PRESSO 静默失败。现改为告警 + 固定 seed。
mr_presso_safe <- function(d, label, seed = PRESSO_SEED) {
  if (nrow(d) < 4) return(NULL)
  r <- tryCatch(MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        data = as.data.frame(d),
        NbDistribution = 1000, SignifThreshold = 0.05, seed = seed),
      error = function(e) {
        warning("MR-PRESSO 失败 [", label, "]: ", conditionMessage(e), call. = FALSE)
        NULL })
  if (is.null(r)) return(NULL)
  g <- r$`MR-PRESSO results`$`Global Test`
  o <- r$`MR-PRESSO results`$`Outlier Test`
  data.frame(Analysis = label, nSNP = nrow(d), RSSobs = g$RSSobs,
             GlobalTest_P = g$Pvalue,
             N_outliers = if (is.null(o)) NA_integer_ else sum(o$Pvalue < 0.05, na.rm = TRUE))
}

## 6.3 完整分层分析
analyse_outcome <- function(harm_dat, tag) {
  sets <- list(
    "All 23 SNPs"          = harm_dat,
    "IL6R locus only"      = harm_dat[harm_dat$locus == "IL6R", ],
    "HLA locus only"       = harm_dat[harm_dat$locus == "HLA", ],
    "All minus rs2228145"  = harm_dat[harm_dat$SNP != LEAD_SNP, ],
    "IL6R minus rs2228145" = harm_dat[harm_dat$locus == "IL6R" & harm_dat$SNP != LEAD_SNP, ])
  main <- do.call(rbind, lapply(names(sets), function(k) mr_panel(sets[[k]], k)))
  het  <- do.call(rbind, lapply(names(sets), function(k) {
    q <- ivw_fun(sets[[k]]$beta.exposure, sets[[k]]$beta.outcome,
                 sets[[k]]$se.exposure, sets[[k]]$se.outcome, "random")
    data.frame(Analysis = k, nSNP = nrow(sets[[k]]), Q = q$Q, Q_df = q$Q_df,
               Q_p = q$Q_pval, I2 = if (q$Q > 0) max(0, (q$Q - q$Q_df) / q$Q) * 100 else 0)
  }))
  egg <- do.call(rbind, lapply(names(sets), function(k) {
    d <- sets[[k]]
    r <- tryCatch(mr_egger_regression(d$beta.exposure, d$beta.outcome,
                                      d$se.exposure, d$se.outcome), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    data.frame(Analysis = k, nSNP = nrow(d), egger_intercept = r$b_i,
               egger_intercept_se = r$se_i, egger_intercept_p = r$pval_i)
  }))
  pres <- do.call(rbind, lapply(names(sets), function(k) mr_presso_safe(sets[[k]], k)))
  list(main = main, het = het, egger = egg, presso = pres,
       loo_il6r = loo_ivw(sets[["IL6R locus only"]]),
       loo_hla  = loo_ivw(sets[["HLA locus only"]]))
}

## ===========================================================================
## 7. Greene 2025
## ===========================================================================
cat("\n== 7. Greene 2025 结局 ===========================================\n")
res_greene <- analyse_outcome(harm, "Greene2025")
print(res_greene$main[, c("Analysis", "nSNP", "Method", "b", "se", "pval", "OR",
                          "OR_lci", "OR_uci")], row.names = FALSE)
if (any(!is.na(res_greene$main$error))) {
  cat("\n  !! 出错的方法:\n")
  print(res_greene$main[!is.na(res_greene$main$error),
        c("Analysis", "Method", "error")], row.names = FALSE)
}
cat("\n-- 异质性 --\n"); print(res_greene$het, row.names = FALSE)
cat("\n-- MR-Egger 截距 --\n"); print(res_greene$egger, row.names = FALSE)
if (!is.null(res_greene$presso)) { cat("\n-- MR-PRESSO --\n"); print(res_greene$presso, row.names = FALSE) }
cat("\n-- IL6R 簇 leave-one-out --\n"); print(res_greene$loo_il6r, row.names = FALSE)
cat("\n-- HLA 簇 leave-one-out --\n");  print(res_greene$loo_hla, row.names = FALSE)

## --- Steiger 方向性检验（复制期，2026-09-26 新增）---
## 二分类结局需先用 get_r_from_lor() 提供 r.outcome，否则 TwoSampleMR 会退回
## “所有性状均按定量处理”的近似口径。
cat("\n-- Steiger 方向性检验（复制期）--\n")
st_dat <- harm
st_dat$r.outcome <- get_r_from_lor(lor = st_dat$beta.outcome, af = st_dat$eaf.outcome,
                                   ncase = GREENE_NCASE, ncontrol = GREENE_NCTRL,
                                   prevalence = GREENE_NCASE / (GREENE_NCASE + GREENE_NCTRL))
st_dat$samplesize.exposure <- N_exp
st_dat$units.exposure <- "SD"; st_dat$units.outcome <- "log odds"
st_dat$ncase.outcome <- GREENE_NCASE; st_dat$ncontrol.outcome <- GREENE_NCTRL
st_dat$prevalence.outcome <- GREENE_NCASE / (GREENE_NCASE + GREENE_NCTRL)
res_greene$steiger <- tryCatch(
  directionality_test(st_dat),
  error = function(e) { warning("Steiger 失败: ", conditionMessage(e), call. = FALSE); NULL })
if (!is.null(res_greene$steiger)) print(res_greene$steiger, row.names = FALSE)

## ===========================================================================
## 8. Sakaue 原结局（同一套代码，保证可比）
## ===========================================================================
cat("\n== 8. Sakaue 原结局（对照）=======================================\n")
sak_dat <- as.data.frame(sak[, c("SNP", "beta.exposure", "se.exposure",
                                 "effect_allele.exposure", "other_allele.exposure",
                                 "eaf.exposure", "beta.outcome", "se.outcome",
                                 "effect_allele.outcome", "other_allele.outcome",
                                 "eaf.outcome", "pval.exposure", "pval.outcome")])
sak_dat$exposure <- "IL-6"; sak_dat$id.exposure <- "IL6"
sak_dat$outcome <- "Keloid (Sakaue 2021)"; sak_dat$id.outcome <- "GCST90018874"
sak_dat$mr_keep <- TRUE; sak_dat$palindromic <- FALSE; sak_dat$ambiguous <- FALSE
sak_dat$samplesize.outcome <- 481912L
sak_dat$locus <- snp_map$locus[match(sak_dat$SNP, snp_map$SNP)]
res_sakaue <- analyse_outcome(sak_dat, "Sakaue2021")
print(res_sakaue$main[, c("Analysis", "nSNP", "Method", "b", "se", "pval", "OR")], row.names = FALSE)

## ===========================================================================
## 9. 并排对照表（进稿件的核心表）
## ===========================================================================
cat("\n== 9. Sakaue vs Greene 并排对照 ===================================\n")
ivw_row <- function(x) {
  x <- as.data.frame(x)
  sel <- x$Method == "IVW (random)"
  if (!any(sel)) sel <- x$Method == "IVW (TwoSampleMR)"
  x[sel, c("Analysis", "nSNP", "b", "se", "pval", "OR", "OR_lci", "OR_uci"), drop = FALSE]
}
cmp <- merge(ivw_row(res_greene$main), ivw_row(res_sakaue$main), by = "Analysis",
             suffixes = c("_Greene", "_Sakaue"))
het_sel <- res_greene$het[, c("Analysis", "Q", "Q_p", "I2")]
names(het_sel) <- c("Analysis", "Q_Greene", "Q_p_Greene", "I2_Greene")
cmp <- merge(cmp, het_sel, by = "Analysis", all.x = TRUE)
cmp$Replicated <- ifelse(cmp$pval_Greene < 0.05,
                         ifelse(sign(cmp$b_Greene) == sign(cmp$b_Sakaue), "Yes (same direction)",
                                "Yes (opposite)"), "No")
print(cmp, row.names = FALSE)

## ===========================================================================
## 10. 输出
## ===========================================================================
cat("\n== 10. 写出结果 ==\n")
STAMP <- format(Sys.Date(), "%Y-%m-%d")
write_out <- function(x, nm) {
  if (is.null(x) || !nrow(x)) return(invisible(NULL))
  fp <- file.path(OUT_DIR, sprintf("%s_%s.csv", nm, STAMP))
  fwrite(as.data.table(x), fp); cat("   ->", basename(fp), "\n")
}
write_out(res_greene$main,      "table.Greene_MR_all_methods")
write_out(cmp,                  "table.Replication_Sakaue_vs_Greene")
write_out(res_greene$het,       "table.Greene_heterogeneity")
write_out(res_greene$egger,     "table.Greene_egger_intercept")
write_out(res_greene$presso,    "table.Greene_MRPRESSO")
write_out(res_greene$steiger,   "table.Greene_steiger")
write_out(res_greene$loo_il6r,  "table.Greene_LOO_IL6R")
write_out(res_greene$loo_hla,   "table.Greene_LOO_HLA")
write_out(res_sakaue$main,      "table.Sakaue_MR_all_methods")

plot_dat <- rbind(
  data.frame(Dataset = "Sakaue 2021 (668 cases)",   ivw_row(res_sakaue$main)),
  data.frame(Dataset = "Greene 2025 (7,837 cases)", ivw_row(res_greene$main)))
p <- ggplot(plot_dat, aes(x = OR, y = Analysis, colour = Dataset)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbarh(aes(xmin = OR_lci, xmax = OR_uci),
                 position = position_dodge(width = 0.6), height = 0.2) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  scale_x_log10() +
  labs(x = "OR per SD higher genetically predicted IL-6 (95% CI)", y = NULL,
       title = "IL-6 and keloid risk: Sakaue 2021 vs Greene 2025 replication") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, sprintf("figure.Replication_forest_%s.png", STAMP)),
       p, width = 9, height = 4.5, dpi = 300)
cat("   -> figure.Replication_forest_", STAMP, ".png\n", sep = "")

## ===========================================================================
## 11. 近似功效（含自检：公式预测 SE vs 实测 IVW SE）
## ===========================================================================
cat("\n== 11. 统计功效（近似）==\n")
power_mr <- function(R2, n_case, n_ctrl, OR, alpha = ALPHA) {
  N_eff <- 4 / (1 / n_case + 1 / n_ctrl)
  ncp   <- log(OR)^2 * R2 * N_eff
  pnorm(sqrt(ncp) - qnorm(1 - alpha / 2))
}
pw <- data.frame(
  Outcome     = c("Sakaue 2021", "Greene 2025"),
  n_case      = c(668, 7837),
  n_ctrl      = c(481244, GREENE_NCTRL),
  Power_OR1.2 = c(power_mr(R2_all, 668, 481244, 1.2),  power_mr(R2_all, GREENE_NCASE, GREENE_NCTRL, 1.2)),
  Power_OR1.5 = c(power_mr(R2_all, 668, 481244, 1.5),  power_mr(R2_all, GREENE_NCASE, GREENE_NCTRL, 1.5)),
  Predicted_SE = sqrt(1 / (R2_all * 4 / (1 / c(668, GREENE_NCASE) + 1 / c(481244, GREENE_NCTRL)))) )
pick_se <- function(x, set = "All 23 SNPs") {
  x <- as.data.frame(x)
  x$se[x$Method == "IVW (random)" & x$Analysis == set][1]
}
pw$Observed_SE <- c(pick_se(res_sakaue$main), pick_se(res_greene$main))
print(pw, row.names = FALSE)
cat("  注：Predicted_SE 与 Observed_SE 应接近；若相差很大，请勿直接引用 Power 数值。\n")
fwrite(pw, file.path(OUT_DIR, sprintf("table.power_%s.csv", STAMP)))

cat("\n  累计 R2 =", round(R2_all, 6),
    "| 总耗时", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "分钟\n")
cat("完成。输出目录:", OUT_DIR, "\n")

