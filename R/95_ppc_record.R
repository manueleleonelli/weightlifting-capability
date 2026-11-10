## 95_ppc_record.R
## -----------------------------------------------------------------------
## Posterior predictive check of the record book. For WL_NDRAW draws of the
## final model, full attempt sequences are simulated for the athlete-day-lifts
## of the analysis set (93_sim_lib.R), and the best made lifts of the
## simulated data are compared with the observed ones:
##   - the local record since 2019 within five kilograms of each bodyweight
##     of Table 6 of the paper;
##   - the 0.99 and 0.999 quantiles of the best made lift on the standardized
##     scale (best - mu) / omega of each draw;
##   - the share of the top one per cent of best lifts (on that scale) made on
##     a day that ended with a successful attempt.
## Outputs: out/ppc_record_bw.csv (records by bodyweight) and
##          out/ppc_record_summary.csv (the other statistics)
source("R/50_ext_lib.R"); source("R/93_sim_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
SUF <- Sys.getenv("WL_TAG_SUFFIX", "_rtau_sm3_mix_fix_log_tau_-3.91202")
N_DRAW <- as.integer(Sys.getenv("WL_NDRAW", "200"))
rec_rows <- list(); sum_rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  f <- sprintf("out/mcmcx_e2_%s_%s%s.rds", tolower(sex), lift, SUF)
  if (!file.exists(f)) { cat("no chain for", sex, lift, "\n"); next }
  D <- readRDS(f)$draws
  d <- build_data_ext(build_data(att, sex, lift, G = 24), jac = 0)
  p_rep <- repeat_rate(d)
  bws <- if (sex == "Men") c(61, 67, 73, 81, 89, 96, 102, 110, 130, 150) else c(49, 55, 59, 64, 71, 76, 81, 86, 100, 120)
  recent <- d$u$year >= 2019
  win <- sapply(bws, function(b) abs(d$u$bodyweight - b) <= 5 & recent)
  records <- function(best) apply(win, 2, function(h) max(best[h], na.rm = TRUE))
  stats <- function(best, ended, mu, om) {
    z <- (best - mu) / om; ok <- !is.na(z)
    top <- ok & z >= quantile(z[ok], 0.99)
    c(q99 = unname(quantile(z[ok], 0.99)), q999 = unname(quantile(z[ok], 0.999)), top_success = mean(ended[top]))
  }
  best_obs <- best_made(d$W, d$S); end_obs <- ended_success(d$W, d$S)
  rec_obs <- records(best_obs)
  idx <- round(seq(1, nrow(D), length.out = min(N_DRAW, nrow(D))))
  R <- matrix(NA_real_, length(idx), length(bws)); ST <- matrix(NA_real_, length(idx), 6)
  for (k in seq_along(idx)) {
    par <- fill_par2(D[idx[k], ]); p <- unpack2(par)
    mu <- as.numeric(d$X %*% p$beta); om <- omega_of(p, d$X)
    sim <- sim_attempts(d, par, p_rep, seed = 5000 + k)
    best <- best_made(sim$W, sim$S)
    R[k, ] <- records(best)
    ST[k, ] <- c(stats(best, ended_success(sim$W, sim$S), mu, om), stats(best_obs, end_obs, mu, om))
    if (k %% 20 == 0) cat(sprintf("%s %s: %d of %d draws\n", sex, lift, k, length(idx)))
  }
  rec_rows[[length(rec_rows) + 1]] <- tibble(sex = sex, lift = lift, bodyweight = bws, record = rec_obs,
    sim_median = apply(R, 2, median), sim_lo = apply(R, 2, quantile, 0.05), sim_hi = apply(R, 2, quantile, 0.95),
    p_sim_ge_obs = colMeans(R >= matrix(rec_obs, nrow(R), length(bws), byrow = TRUE)))
  nm <- c("q99", "q999", "top_success")
  for (s in seq_along(nm)) sum_rows[[length(sum_rows) + 1]] <- tibble(sex = sex, lift = lift, statistic = nm[s],
    observed = median(ST[, 3 + s]), sim_median = median(ST[, s]), sim_lo = quantile(ST[, s], 0.05), sim_hi = quantile(ST[, s], 0.95),
    p_sim_ge_obs = mean(ST[, s] >= ST[, 3 + s]))
}
rb <- bind_rows(rec_rows); sb <- bind_rows(sum_rows)
write_csv(rb, "out/ppc_record_bw.csv"); write_csv(sb, "out/ppc_record_summary.csv")
print(as.data.frame(rb |> mutate(across(where(is.numeric), ~ round(.x, 3)))))
print(as.data.frame(sb |> mutate(across(where(is.numeric), ~ round(.x, 3)))))
