## 58_ext_summaries.R
## -----------------------------------------------------------------------
## Posterior summaries of the final model (extended kernel: log margin,
## situation-specific lower scales, relative attempt scale fixed, athlete
## factor) from the chains out/mcmcx_<tag>.rds, one per sex and lift:
##   out/postx_params.csv     parameter medians, intervals, ESS
##   out/postx_quantiles.csv  capability quantiles by bodyweight, 2026
##   out/postx_ceiling.csv    endpoint against the local record of the current era (median over all
##                            draws, infinite endpoints included, as in 68_; interval given finite)
##   out/postx_latent.csv     capability against the best made lift, by type of day
## WL_TAG_SUFFIX gives the common suffix of the chain files (default
## "_rtau_sm3_fix_log_tau_-3.91202"); WL_NDRAW draws for the derived quantities.
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
N_DRAW <- as.integer(Sys.getenv("WL_NDRAW", "300"))
SUF <- Sys.getenv("WL_TAG_SUFFIX", "_rtau_sm3_mix_fix_log_tau_-3.91202")
local_record <- function(u, best, bw, half = 5) {
  h <- abs(u$bodyweight - bw) <= half & !is.na(best)
  tibble(record = max(best[h & u$year >= 2019]), record_all = max(best[h]), n = sum(h & u$year >= 2019), n_all = sum(h))
}
ess <- function(v) { a <- acf(v, lag.max = 300, plot = FALSE)$acf[-1]; n <- length(v); n / (1 + 2 * sum(a[cumsum(a < 0) == 0])) }
par_rows <- list(); q_rows <- list(); c_rows <- list(); l_rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  f <- sprintf("out/mcmcx_e2_%s_%s%s.rds", tolower(sex), lift, SUF)
  if (!file.exists(f)) { cat("no chain for", sex, lift, "\n"); next }
  m <- readRDS(f); D <- m$draws; ref <- m$ref
  cat(sprintf("== %s %s: %d draws, acceptance %.3f\n", sex, lift, nrow(D), m$acc))
  d <- build_data_ext(build_data(att, sex, lift, G = 24), jac = m$jac)
  best_all <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  nat <- cbind(D, tau = exp(D[, "log_tau"]), sL = exp(D[, "lam"]), nu = exp(D[, "lnu"]),
               sp1 = exp(D[, "lsp1"]), sp2 = exp(D[, "lsp2"]), sp3 = exp(D[, "lsp3"]),
               sm1 = exp(D[, "lsm"]), sm2 = exp(D[, "lsm2"]), sm3 = exp(D[, "lsm3"]))
  if ("mx_u" %in% colnames(D) && all(is.finite(D[, "mx_u"]))) {   # normal-mixture body: threshold probability and tail scale
    Hu <- apply(D, 1, function(x) mix_Hu(unpack2(fill_par2(x))))
    nat <- cbind(nat, Hu = Hu, sigma_tail = exp(D[, "mx_lsig"]))
  }
  for (nm in colnames(nat)) {
    v <- nat[, nm]; if (all(is.na(v))) next
    dgi <- if (!is.null(m$diag)) match(nm, m$diag$parameter) else NA   # across-chain diagnostics of 86_combine_chains.R
    par_rows[[length(par_rows) + 1]] <- tibble(sex = sex, lift = lift, parameter = nm, median = median(v),
      lo = quantile(v, .025), hi = quantile(v, .975), sd = sd(v),
      ess = if (!is.na(dgi)) m$diag$ess_bulk[dgi] else if (nm %in% m$free) ess(v) else NA,
      rhat = if (!is.na(dgi)) m$diag$rhat[dgi] else NA)
  }
  par_rows[[length(par_rows) + 1]] <- tibble(sex = sex, lift = lift, parameter = "P(xi<0)", median = mean(D[, "xi"] < 0), lo = NA, hi = NA, sd = NA, ess = NA, rhat = NA)
  bws <- if (sex == "Men") c(61, 67, 73, 81, 89, 96, 102, 110, 130, 150) else c(49, 55, 59, 64, 71, 76, 81, 86, 100, 120)
  idx <- round(seq(1, nrow(D), length.out = min(N_DRAW, nrow(D))))
  for (bw in bws) {
    xd <- design_at(ref, bw, 2026)
    Qm <- sapply(idx, function(i) {
      p <- unpack2(D[i, ]); mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
      c(mu + om * std_quantile(p, c(0.9, 0.99, 0.999, 0.9999)), mu + om * std_endpoint(p), mu, om)
    })
    for (k in 1:4) q_rows[[length(q_rows) + 1]] <- tibble(sex = sex, lift = lift, bodyweight = bw,
      level = c(0.9, 0.99, 0.999, 0.9999)[k], median = median(Qm[k, ]), lo = quantile(Qm[k, ], .025), hi = quantile(Qm[k, ], .975))
    q_rows[[length(q_rows) + 1]] <- tibble(sex = sex, lift = lift, bodyweight = bw, level = NA,
      median = median(Qm[6, ]), lo = quantile(Qm[6, ], .025), hi = quantile(Qm[6, ], .975))
    rec <- local_record(d$u, best_all, bw); ep <- Qm[5, ]
    c_rows[[length(c_rows) + 1]] <- tibble(sex = sex, lift = lift, bodyweight = bw,
      record = rec$record, record_all = rec$record_all, n_lifts = rec$n, n_lifts_all = rec$n_all,
      p_finite = mean(is.finite(ep)), p_within_5 = mean(ep <= 1.05 * rec$record), p_within_10 = mean(ep <= 1.10 * rec$record),
      p_within_25 = mean(ep <= 1.25 * rec$record), p_within_50 = mean(ep <= 1.50 * rec$record),
      ceiling_median = median(ep), ceiling_lo = quantile(ep[is.finite(ep)], .025), ceiling_hi = quantile(ep[is.finite(ep)], .975),
      q999_median = median(Qm[3, ]), p_record_above_q999 = mean(rec$record > Qm[3, ]))
  }
  ## capability against the best made lift, at twenty draws, athlete factor off
  ended_success <- !is.na(best_all) & apply(d$Sint, 1, function(s) { s <- s[s >= 0]; s[length(s)] == 1 })
  bracketed <- !is.na(best_all) & d$has_miss & d$has_make
  type <- ifelse(is.na(best_all), "no made lift", ifelse(ended_success, "ended on a success", "bracketed"))
  gaps <- list()
  for (i in round(seq(1, nrow(D), length.out = 20))) {
    L <- latent_summary_ext(D[i, ], d)
    gaps[[length(gaps) + 1]] <- tibble(type = type, gap = L[, 1] - best_all, rel = (L[, 1] - best_all) / best_all, zsd = L[, 2])
  }
  gl <- bind_rows(gaps, .id = "draw")
  for (tp in c("bracketed", "ended on a success", "all")) {
    sub <- if (tp == "all") gl |> filter(type != "no made lift") else gl |> filter(type == tp)
    byd <- sub |> group_by(draw) |> summarise(gap = mean(gap), rel = mean(rel), zsd = mean(zsd), .groups = "drop")
    n_tp <- if (tp == "all") sum(type != "no made lift") else sum(type == tp)
    l_rows[[length(l_rows) + 1]] <- tibble(sex = sex, lift = lift, type = tp, n = n_tp,
      gap_mean = median(byd$gap), gap_lo = quantile(byd$gap, .025), gap_hi = quantile(byd$gap, .975),
      rel_mean = median(byd$rel), latent_sd = median(byd$zsd))
  }
}
write_csv(bind_rows(par_rows), "out/postx_params.csv"); write_csv(bind_rows(q_rows), "out/postx_quantiles.csv")
write_csv(bind_rows(c_rows), "out/postx_ceiling.csv"); write_csv(bind_rows(l_rows), "out/postx_latent.csv")
print(as.data.frame(bind_rows(par_rows) |> filter(parameter %in% c("xi", "kappa", "rhoa", "delta2", "delta3", "sL", "nu", "m1", "m2", "m3", "sp1", "sp2", "sp3", "sm1", "sm2", "sm3", "g0", "b3", "P(xi<0)")) |> mutate(across(where(is.numeric), ~round(.x, 3)))))
print(as.data.frame(bind_rows(c_rows) |> filter(bodyweight %in% c(73, 89, 110, 150, 59, 71, 86, 120)) |> mutate(across(where(is.numeric), ~round(.x, 2)))))
print(as.data.frame(bind_rows(l_rows) |> mutate(across(where(is.numeric), ~round(.x, 3)))))
