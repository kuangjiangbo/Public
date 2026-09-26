## v5 Figure 7 重绘：逻辑分层顺序 + 300 dpi TIF
BASE <- "F:/MR/Interleukin-6 levels and keloids"
IN   <- file.path(BASE, "数据/Greene2025/table.Replication_Sakaue_vs_Greene_2026-09-24.csv")
SUB  <- file.path(BASE, "投稿期刊/新建文件夹")
FIGD <- file.path(BASE, "图表/Figures")

suppressMessages({library(data.table); library(ggplot2)})
cmp <- as.data.frame(fread(IN))
cat("读入:", nrow(cmp), "行\n"); print(cmp[, c("Analysis", "OR_Sakaue", "OR_Greene")])

lv <- c("HLA locus only", "IL6R minus rs2228145", "IL6R locus only",
        "All minus rs2228145", "All 23 SNPs")   # ggplot: 第一个在最下 -> 这样顶端是 All 23 SNPs
plot_dat <- rbind(
  data.frame(Dataset = "Discovery: Sakaue 2021 (668 cases)",    Analysis = cmp$Analysis,
             OR = cmp$OR_Sakaue, lci = cmp$OR_lci_Sakaue, uci = cmp$OR_uci_Sakaue),
  data.frame(Dataset = "Replication: Greene 2025 (7,837 cases)", Analysis = cmp$Analysis,
             OR = cmp$OR_Greene, lci = cmp$OR_lci_Greene, uci = cmp$OR_uci_Greene))
plot_dat$Analysis <- factor(plot_dat$Analysis, levels = lv)

p <- ggplot(plot_dat, aes(y = Analysis, x = OR, colour = Dataset)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbar(aes(xmin = lci, xmax = uci), orientation = "y",
                width = 0.18, linewidth = 0.5,
                position = position_dodge(width = 0.55)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.55)) +
  scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6)) +
  scale_colour_manual(values = c("Discovery: Sakaue 2021 (668 cases)" = "#1B7F79",
                                 "Replication: Greene 2025 (7,837 cases)" = "#C0392B")) +
  labs(x = "OR per SD higher genetically predicted IL-6 (95% CI)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        axis.text.y = element_text(size = 10))
ggsave(file.path(SUB,  "figure.Replication_forest_v5_2026-09-26.tif"),
       p, width = 7.2, height = 3.6, dpi = 300, compression = "lzw")
ggsave(file.path(SUB,  "figure.Replication_forest_v5_2026-09-26.png"),
       p, width = 7.2, height = 3.6, dpi = 300)
ggsave(file.path(FIGD, "Figure_7_Replication_Forest.tif"),
       p, width = 7.2, height = 3.6, dpi = 300, compression = "lzw")
cat("已输出 TIF/PNG\n")
