## 00_timing.R : time of one evaluation of the log posterior of the final model (men's snatch)
suppressMessages(source("R/50_ext_lib.R"))
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE, col_types = cols(category = col_character(), .default = col_guess()))
d <- build_data_ext(build_data(att, "Men", "snatch", G = 24), jac = 0)
p <- fill_par2(c(readRDS("R/attempt_starts.rds")$men_snatch, log_tau = log(0.02),
  b0 = 160, b1 = 20, b2 = -3, b3 = 0, g0 = 2.7, g1 = 0.07, g2 = -0.01, xi = -0.3, mx_m2 = -2, mx_m3 = -1.4, mx_ls2 = 1, mx_ls3 = 0.4,
  mx_a2 = -3, mx_a3 = -0.6, mx_u = 0.7, mx_lsig = -0.3))
v <- log_post_ext(p, d); t0 <- Sys.time(); for (i in 1:5) v <- log_post_ext(p, d)
writeLines(as.character(as.numeric(Sys.time() - t0, units = "secs") / 5), "out/timing_eval.txt")
