# MVMR Forest Plot
library(ggplot2)

# Data for forest plot
dat <- data.frame(
  Model = c(
    "Univariate IVW (ref)",
    "MVMR: IL-6 + CRP",
    "MVMR: IL-6 + CRP",
    "MVMR: IL-6 + IL-1RA + TNFR1",
    "MVMR: IL-6 + IL-1RA + TNFR1 + CRP",
    "MVMR: IL-6 + IL-1RA + TNFR1 + CRP",
    "MVMR: IL-6 + IL-1RA + TNFR1 + CRP",
    "MVMR: IL-6 + IL-1RA + TNFR1 + CRP"
  ),
  Exposure = c(
    "IL-6",
    "IL-6", "CRP",
    "IL-6",
    "IL-6", "IL-1RA", "TNFR1", "CRP"
  ),
  OR = c(1.214, 1.150, 0.879, 1.057, 1.053, 0.979, 1.093, 0.874),
  Lower = c(1.016, 1.029, 0.819, 0.914, 0.939, 0.928, 0.873, 0.810),
  Upper = c(1.452, 1.285, 0.943, 1.223, 1.181, 1.033, 1.368, 0.943),
  Pval = c(0.033, 0.014, 0.0003, 0.457, 0.375, 0.436, 0.439, 0.0006)
)

# Create label column
dat$Label <- paste0(dat$Model, " | ", dat$Exposure)
dat$Label <- factor(dat$Label, levels = rev(dat$Label))

# Create proper model grouping labels
dat$Group <- ""
for (i in 1:nrow(dat)) {
  if (dat$Model[i] == "Univariate IVW (ref)") {
    dat$Group[i] <- "Univariate MR"
  } else if (dat$Model[i] == "MVMR: IL-6 + CRP") {
    dat$Group[i] <- "MVMR: IL-6 + CRP"
  } else if (dat$Model[i] == "MVMR: IL-6 + IL-1RA + TNFR1") {
    dat$Group[i] <- "MVMR: IL-6 + IL-1RA + TNFR1"
  } else {
    dat$Group[i] <- "MVMR: All four exposures"
  }
}

# Manual y positions for better grouping
dat$y <- nrow(dat):1

# Reference line for OR=1
ref_line <- 1

p <- ggplot(dat, aes(x = OR, y = y)) +
  geom_vline(xintercept = ref_line, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, color = Group), height = 0.15, linewidth = 0.8) +
  geom_point(aes(color = Group, fill = Group), size = 3, shape = 21, stroke = 0.8) +
  scale_y_continuous(
    breaks = dat$y,
    labels = paste0(dat$Exposure, "  (", dat$Pval, ")"),
    name = ""
  ) +
  scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.7),
                name = "Odds Ratio (95% CI)") +
  scale_color_manual(values = c(
    "Univariate MR" = "#E41A1C",
    "MVMR: IL-6 + CRP" = "#377EB8",
    "MVMR: IL-6 + IL-1RA + TNFR1" = "#4DAF4A",
    "MVMR: All four exposures" = "#984EA3"
  )) +
  scale_fill_manual(values = c(
    "Univariate MR" = "#E41A1C",
    "MVMR: IL-6 + CRP" = "#377EB8",
    "MVMR: IL-6 + IL-1RA + TNFR1" = "#4DAF4A",
    "MVMR: All four exposures" = "#984EA3"
  )) +
  labs(
    title = "Multivariable Mendelian Randomization Results",
    subtitle = "Causal effects of inflammatory biomarkers on keloid risk",
    caption = "P values shown in parentheses. Univariate IVW result shown as reference."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.y = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  coord_cartesian(xlim = c(0.7, 1.7))

# Save
setwd("F:/MR/Interleukin-6 levels and keloids/数据/Multivariable MR 数据")
ggsave("pic.MVMR_forest_plot.png", plot = p, width = 9, height = 5, dpi = 300)
cat("✓ Saved: pic.MVMR_forest_plot.png\n")
