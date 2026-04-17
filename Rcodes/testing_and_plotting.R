library(dplyr)
library(ggplot2)
results <- list()
process_cycle <- function(prefix) {
  # LOAD
  data_clean <- readRDS(paste0(prefix, "data_clean.rds"))
  ksd_exact  <- readRDS(paste0(prefix, "ksd_exact.rds"))
  out_md     <- readRDS(paste0(prefix, "out_md.rds"))
  # KSD
  threshold_ksd <- quantile(ksd_exact, 0.05)
  out_ksd <- ksd_exact < threshold_ksd
  # INDICES
  idx_md  <- which(out_md)
  idx_ksd <- which(out_ksd)
  ksd_only <- setdiff(idx_ksd, idx_md)
  md_only  <- setdiff(idx_md, idx_ksd)
  common   <- intersect(idx_md, idx_ksd)
  ksd_values <- data_clean$LBXGH[ksd_only]
  md_values  <- data_clean$LBXGH[md_only]
  wilcox_test <- wilcox.test(
    ksd_values,
    md_values,
    alternative = "greater",
    paired = FALSE
  )
  summarize_group <- function(idx) {
    data_clean[idx, ] %>%
      summarise(
        glucose = mean(LBXGLU),
        hba1c   = mean(LBXGH),
        trig    = mean(LBXTLG),
        hdl     = mean(LBDHDD)
      )
  }
  list(
    n = nrow(data_clean),
    ksd_only_n = length(ksd_only),
    md_only_n  = length(md_only),
    common_n   = length(common),
    ksd_stats = summarize_group(ksd_only),
    md_stats  = summarize_group(md_only),
    wilcox_pvalue = wilcox_test$p.value,
    
    ksd_only_idx = ksd_only,
    md_only_idx  = md_only
  )
}
results[["15_16"]] <- process_cycle("15_16_")
results[["17_18"]] <- process_cycle("17_18_")
results[["21_23"]] <- process_cycle("21_23_")

data.frame(
  Cycle = c("15-16","17-18","21-23"),
  KSD_HbA1c = c(
    results[["15_16"]]$ksd_stats$hba1c,
    results[["17_18"]]$ksd_stats$hba1c,
    results[["21_23"]]$ksd_stats$hba1c
  ),
  MD_HbA1c = c(
    results[["15_16"]]$md_stats$hba1c,
    results[["17_18"]]$md_stats$hba1c,
    results[["21_23"]]$md_stats$hba1c
  ),
  Wilcox_pvalue = c(
    results[["15_16"]]$wilcox_pvalue,
    results[["17_18"]]$wilcox_pvalue,
    results[["21_23"]]$wilcox_pvalue
  )
)

results[["15_16"]]$wilcox_pvalue
results[["17_18"]]$wilcox_pvalue
results[["21_23"]]$wilcox_pvalue

prefix <- "21_23_"   # change to other prefixes as well
data_clean <- readRDS(paste0(prefix, "data_clean.rds"))
ksd_exact  <- readRDS(paste0(prefix, "ksd_exact.rds"))
out_md     <- readRDS(paste0(prefix, "out_md.rds"))
threshold_ksd <- quantile(ksd_exact, 0.05)
out_ksd <- ksd_exact < threshold_ksd
idx_md  <- which(out_md)
idx_ksd <- which(out_ksd)
ksd_only <- setdiff(idx_ksd, idx_md)
md_only  <- setdiff(idx_md, idx_ksd)
common   <- intersect(idx_md, idx_ksd)
X <- data_clean %>%
  select(LBXGLU, LBXGH, LBXTLG, LBDHDD, RIDAGEYR) %>%
  scale()
pca <- prcomp(X)
df <- as.data.frame(pca$x[,1:2])
colnames(df) <- c("PC1", "PC2")
df$group <- "Normal"
df$group[ksd_only] <- "KSD_only"
df$group[md_only]  <- "MD_only"
df$group[common]   <- "Both"
ggplot(df, aes(PC1, PC2, color = group)) +
  geom_point(
    data = subset(df, group == "Normal"),
    alpha = 0.4, size = 1.2
  ) +
  geom_point(
    data = subset(df, group != "Normal"),
    size = 2.5
  ) +
  scale_color_manual(values = c(
    "Normal" = "grey70",
    "KSD_only" = "red",
    "MD_only" = "blue",
    "Both" = "black"
  )) +
  labs(
    title = paste("Outlier Detection Year: ", sep = ""),
    subtitle = "Comparison across NHANES cycles (2015–16, 2017–18, 2021–23)",
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Group"
  ) +
  theme_minimal(base_size = 14)

