## 87_merge_ppc.R : merges the per-combination predictive checks of 53_ext_checks.R
## into out/ppc_ext_attempts.csv, out/ppc_ext_jumps.csv and out/ppc_ext_openers.csv
suppressPackageStartupMessages({library(dplyr); library(readr)})
for (k in c("attempts", "jumps", "openers")) {
  fs <- Sys.glob(sprintf("out/ppc/*_%s.csv", k))
  if (length(fs)) write_csv(bind_rows(lapply(fs, read_csv, show_col_types = FALSE)), sprintf("out/ppc_ext_%s.csv", k))
}
cat("merged\n")
