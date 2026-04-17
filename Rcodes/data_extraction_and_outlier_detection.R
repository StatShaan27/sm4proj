library(haven)
library(dplyr)
library(tidyr)
library(dbscan)
library(beepr)

prefix <- "17_18_" # change to other prefixes as well
demo <- read_xpt(paste0(prefix, "DEMO_L.xpt"))
bmx  <- read_xpt(paste0(prefix, "BMX_L.xpt"))
glu  <- read_xpt(paste0(prefix, "GLU_L.xpt"))
ghb  <- read_xpt(paste0(prefix, "GHB_L.xpt"))
trig <- read_xpt(paste0(prefix, "TRIGLY_L.xpt"))
hdl  <- read_xpt(paste0(prefix, "HDL_L.xpt"))
dim(demo)
dim(bmx)
dim(glu)
dim(ghb)
dim(trig)
dim(hdl)
data <- demo %>%
  select(SEQN, RIDAGEYR) %>%
  left_join(bmx %>% select(SEQN, BMXBMI), by = "SEQN") %>%
  left_join(glu %>% select(SEQN, LBXGLU), by = "SEQN") %>%
  left_join(ghb %>% select(SEQN, LBXGH), by = "SEQN") %>%
  left_join(
    trig %>% 
      select(SEQN, LBXTLG) %>% 
      rename(LBXTLG = LBXTLG),
    by = "SEQN"
  ) %>%
  left_join(hdl %>% select(SEQN, LBDHDD), by = "SEQN")
data_clean <- data %>%
  filter(RIDAGEYR >= 18,
         BMXBMI >= 18.5, BMXBMI < 25) %>%
  drop_na()
dim(data_clean)
X <- data_clean %>%
  select(LBXGLU, LBXGH, LBXTLG, LBDHDD, RIDAGEYR)
X <- scale(X)
X_mat <- as.matrix(X)
dim(X_mat)


#md
md <- mahalanobis(X_mat, colMeans(X_mat), cov(X_mat))
threshold_md <- quantile(md, 0.95)
out_md <- md > threshold_md
sum(out_md)

#ksd
dists <- as.matrix(dist(X_mat))
sigma <- median(dists)
K <- exp(-(dists^2) / (2 * sigma^2))

n <- nrow(X_mat)
ksd_exact <- numeric(n)
cat("\nStarting EXACT KSD...\n")
start_time <- Sys.time()
for (i in 1:n) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("i=%d/%d (%.1f%%) | %.2f min\r", i, n, 100*i/n, elapsed))
  s_val <- 0
  for (j in 1:n) {
    for (k in 1:n) {
      num <- K[i,i] + K[j,k] - K[i,j] - K[i,k]
      denom1 <- sqrt(K[i,i] + K[j,j] - 2*K[i,j])
      denom2 <- sqrt(K[i,i] + K[k,k] - 2*K[i,k])
      if (denom1 > 0 && denom2 > 0) {
        s_val <- s_val + num / (denom1 * denom2)
      }
    }
  }
  
  ksd_exact[i] <- 1 - (s_val / (n^2))
}
cat("\nKSD completed in:",
    difftime(Sys.time(), start_time, units = "mins"), "minutes\n")
beepr::beep()

# KSD OUTLIERS
threshold_ksd <- quantile(ksd_exact, 0.05)
out_ksd <- ksd_exact < threshold_ksd

# INDEX SETS
idx_md  <- which(out_md)
idx_ksd <- which(out_ksd)
ksd_only <- setdiff(idx_ksd, idx_md)
md_only  <- setdiff(idx_md, idx_ksd)
common   <- intersect(idx_md, idx_ksd)


#verification:
cat("MD:", length(idx_md), "\n")
cat("KSD:", length(idx_ksd), "\n")
cat("Overlap:", length(common), "\n")
cat("KSD only:", length(ksd_only), "\n")
cat("MD only:", length(md_only), "\n")


#analysis
analyze <- function(idx, name, df) {
  cat(paste("\n---", name, "---\n"))
  print(
    df[idx, ] %>%
      summarise(
        glucose = mean(LBXGLU),
        hba1c   = mean(LBXGH),
        trig    = mean(LBXTLG),
        hdl     = mean(LBDHDD)
      )
  )
}
analyze(ksd_only, "KSD only", data_clean)
analyze(md_only, "MD only", data_clean)

