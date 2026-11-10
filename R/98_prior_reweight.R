## 98_prior_reweight.R
## -----------------------------------------------------------------------
## Prior sensitivity of the tail by importance reweighting of the posterior
## draws of the final model: each draw gets the weight
##   p_new(theta) / p_reference(theta),
## which changes only the prior of xi or of the threshold u. Variants:
##   xi_wide   xi ~ N(0, 0.5^2) instead of N(0, 0.3^2)
##   u_wide    u ~ N(1.5, 2^2) instead of N(1.5, 1), same support (0.2, 4)
##   u_max3    the threshold bounded below 3 standardized units
## Reported for each sex and lift at the reference bodyweight in 2026: the
## posterior probability of a finite endpoint, the posterior median and 90 per
## cent interval of the endpoint, the median of the 0.999 quantile and of u,
## and the effective sample size of the weights.
## Settings: WL_NDRAW (draws used, default all).
## Output: out/prior_reweight.csv
source("R/50_ext_lib.R"); source("R/93_sim_lib.R")
SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
N_DRAW <- as.integer(Sys.getenv("WL_NDRAW", "1000000"))
wq <- function(x, w, p) { o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w); x[which(cw >= p)[1]] }
rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  f <- sprintf("out/mcmcx_e2_%s_%s%s.rds", tolower(sex), lift, SUF)
  if (!file.exists(f)) next
  x <- readRDS(f); D <- x$draws
  idx <- round(seq(1, nrow(D), length.out = min(N_DRAW, nrow(D))))
  Q <- t(sapply(idx, function(i) ref_quantities(fill_par2(D[i, ]), x$ref, sex)))
  xi <- Q[, "xi"]; u <- Q[, "u"]
  W <- list(reference = rep(1, length(idx)),
            xi_wide = dnorm(xi, 0, 0.5) / dnorm(xi, 0, 0.3),
            u_wide = dnorm(u, 1.5, 2) / dnorm(u, 1.5, 1),
            u_max3 = as.numeric(u < 3))
  for (v in names(W)) { w <- W[[v]]
    rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, variant = v,
      p_finite = sum(w * (xi < 0)) / sum(w), xi_median = wq(xi, w, 0.5), u_median = wq(u, w, 0.5),
      q999_median = wq(Q[, "q999"], w, 0.5), endpoint_median = wq(Q[, "endpoint"], w, 0.5),
      endpoint_lo = wq(Q[, "endpoint"], w, 0.05), endpoint_hi = wq(Q[, "endpoint"], w, 0.95),
      ess = sum(w)^2 / sum(w^2), n_draws = length(idx)) }
  cat(sex, lift, "done\n")
}
out <- bind_rows(rows); write_csv(out, "out/prior_reweight.csv")
print(as.data.frame(out |> mutate(across(where(is.numeric), ~ round(.x, 3)))))
