## 63_bodyfree_pot.R
## -----------------------------------------------------------------------
## Body-free check of the tail: no model for the bulk at all. The threshold
## u(x) is the q-quantile of the best made lift, estimated by linear quantile
## regression on the design of the paper (intercept, standardized log
## bodyweight, its orthogonal square, centred year). The exceedances
## y - u(x) > 0 follow a generalized Pareto distribution with log scale
## linear in log bodyweight and a constant shape. The endpoint at covariate
## x is u(x) + sigma(x) / (-xi). Point estimates use all athlete-days
## (independence likelihood); intervals come from a bootstrap that resamples
## athletes, which carries the dependence between the days of one athlete.
## Output: out/bodyfree_pot.csv
source("R/30_marginal_lib.R")
suppressPackageStartupMessages(library(quantreg))
set.seed(12)
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
B <- as.integer(Sys.getenv("WL_BOOT", "1000"))
fit_pot <- function(y, X, q, start = NULL) {
  qr <- rq.fit(X, y, tau = q, method = "br")
  u <- as.numeric(X %*% qr$coefficients)
  h <- y > u; r <- y[h] - u[h]; Z <- X[h, 1:2, drop = FALSE]
  nll <- function(th) {
    s <- exp(as.numeric(Z %*% th[1:2])); k <- th[3]
    a <- 1 + k * r / s
    if (any(a <= 0)) return(1e10)
    sum(log(s)) + (1 / k + 1) * sum(log(a))
  }
  st <- if (is.null(start)) c(log(mean(r)), 0, -0.1) else start
  o <- optim(st, nll, control = list(maxit = 3000, reltol = 1e-12))
  o <- optim(o$par, nll, method = "BFGS", control = list(maxit = 1000))
  list(beta_u = qr$coefficients, th = o$par, n_exc = sum(h))
}
endpoint_at <- function(f, xd) {
  u <- sum(f$beta_u * xd); s <- exp(sum(f$th[1:2] * xd[1:2])); k <- f$th[3]
  if (k >= 0) Inf else u + s / (-k)
}
final <- read_csv("out/postx_ceiling.csv", show_col_types = FALSE)
rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  d <- build_data(att, sex, lift, G = 8)
  best <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  ok <- !is.na(best); y <- best[ok]; X <- d$X[ok, ]; ath <- d$u$athlete[ok]
  idx_by_ath <- split(seq_along(y), ath); A <- length(idx_by_ath)
  bws <- if (sex == "Men") c(73, 89, 110, 150) else c(59, 71, 86, 120)
  for (q in c(0.90, 0.95)) {
    f <- fit_pot(y, X, q)
    bs <- replicate(B, {
      ii <- unlist(idx_by_ath[sample.int(A, A, replace = TRUE)], use.names = FALSE)
      fb <- tryCatch(fit_pot(y[ii], X[ii, ], q, start = f$th), error = function(e) NULL)
      if (is.null(fb)) rep(NA, 1 + length(bws)) else c(fb$th[3], sapply(bws, function(bw) endpoint_at(fb, design_at(d$ref, bw, 2026))))
    })
    cat(sprintf("%s %s q = %.2f: %d exceedances, xi %.3f (%.3f, %.3f)\n", sex, lift, q, f$n_exc, f$th[3],
                quantile(bs[1, ], .025, na.rm = TRUE), quantile(bs[1, ], .975, na.rm = TRUE)))
    for (k in seq_along(bws)) {
      ep <- bs[1 + k, ]; fr <- final[final$sex == sex & final$lift == lift & final$bodyweight == bws[k], ]
      rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, q = q, bodyweight = bws[k], record = fr$record,
        n_exc = f$n_exc, xi = f$th[3], xi_lo = quantile(bs[1, ], .025, na.rm = TRUE), xi_hi = quantile(bs[1, ], .975, na.rm = TRUE),
        endpoint = endpoint_at(f, design_at(d$ref, bws[k], 2026)),
        ep_lo = quantile(ep, .025, na.rm = TRUE), ep_hi = quantile(ep, .975, na.rm = TRUE), p_finite = mean(is.finite(ep), na.rm = TRUE),
        latent = fr$ceiling_median, latent_lo = fr$ceiling_lo, latent_hi = fr$ceiling_hi)
    }
  }
}
out <- bind_rows(rows); write_csv(out, "out/bodyfree_pot.csv")
print(as.data.frame(out |> mutate(across(where(is.numeric), ~round(.x, 3)))), width = 220)
