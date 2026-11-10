## 62_bestlift_evt.R
## -----------------------------------------------------------------------
## Benchmark without the attempt and declaration layers: the capability
## distribution of the paper (location and log scale in bodyweight and year,
## half-t / half-normal body, generalized Pareto tail above the fixed splice,
## jump parameter kappa) fitted directly to the best made lift of each
## athlete-day, treating the best lift as the capability.
##   bl_exch   athlete-days exchangeable
##   bl_ath    days of one athlete sharing a Gaussian factor on the copula
##             scale (equicorrelation rho_a), as in the final model
##   bl_max    one value per athlete: the largest standardized residual of the
##             athlete over his or her days (residuals from bl_ath), generalized
##             Pareto fitted above its 0.9 quantile; the endpoint of the
##             athlete maxima is the endpoint of the day distribution
## Same priors as the paper for the shared parameters; posterior modes and
## Laplace intervals (2,000 draws on the transformed scale).
## Output: out/bestlift_evt.csv (shape and endpoint by sex, lift, bodyweight,
## with the final latent model and the local record for comparison).
source("R/50_ext_lib.R")   # also sources 30_marginal_lib.R
set.seed(11)
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
NB <- c("b0", "b1", "b2", "b3", "g0", "g1", "g2", "xi", "kappa", "lam", "lnu", "rhoa")
BND <- list(xi = c(XI_LO, XI_HI), kappa = c(-1.5, 1.5), lam = c(-2, 2), lnu = c(log(2), log(200)), rhoa = c(0, 0.95))
to_r <- function(p) { for (nm in intersect(names(BND), names(p))) { b <- BND[[nm]]; p[nm] <- qlogis((p[nm] - b[1]) / (b[2] - b[1])) }; p }
from_r <- function(t) { for (nm in intersect(names(BND), names(t))) { b <- BND[[nm]]; t[nm] <- b[1] + (b[2] - b[1]) * plogis(t[nm]) }; t }
full21 <- function(p) {
  v <- rep(NA_real_, 21); names(v) <- PAR_NAMES
  v[c("b0", "b1", "b2", "b3", "g0", "g1", "g2", "xi", "kappa", "lam", "lnu")] <- p[c("b0", "b1", "b2", "b3", "g0", "g1", "g2", "xi", "kappa", "lam", "lnu")]
  v
}
logpost_bl <- function(p, y, X, aid, athlete = TRUE) {
  beta <- p[1:4]; gam <- p[5:7]; xi <- p[["xi"]]; ka <- p[["kappa"]]; sL <- exp(p[["lam"]]); nu <- exp(p[["lnu"]])
  mu <- as.numeric(X %*% beta); om <- as.numeric(exp(X[, 1:3] %*% gam))
  e <- (y - mu) / om
  lf <- log_f_std(e, xi, ka, sL, nu) - log(om)
  if (any(!is.finite(lf))) return(-Inf)
  ll <- sum(lf)
  if (athlete) {
    rho <- p[["rhoa"]]
    cz <- qnorm(pmin(pmax(1 - surv_std(e, xi, ka, sL, nu), 1e-12), 1 - 1e-12))
    n <- as.numeric(rowsum(rep(1, length(cz)), aid)); s1 <- as.numeric(rowsum(cz, aid)); s2 <- as.numeric(rowsum(cz^2, aid))
    logdet <- (n - 1) * log(1 - rho) + log(1 + (n - 1) * rho)
    Q <- (s2 - rho / (1 + (n - 1) * rho) * s1^2) / (1 - rho)
    ll <- ll + sum(-0.5 * logdet - 0.5 * Q + 0.5 * s2)
  }
  lp <- sum(dnorm(beta, 0, 300, log = TRUE)) + dnorm(gam[1], log(15), 0.7, log = TRUE) + sum(dnorm(gam[2:3], 0, 1, log = TRUE)) +
    dnorm(xi, 0, 0.3, log = TRUE) + dnorm(ka, 0, 0.35, log = TRUE) + dnorm(p[["lam"]], 0, 1, log = TRUE) + dnorm(p[["lnu"]], log(10), 1.5, log = TRUE)
  ll + lp
}
fit_bl <- function(start, y, X, aid, athlete) {
  nms <- if (athlete) NB else NB[-12]
  t0 <- to_r(start[nms])
  obj <- function(t) { names(t) <- nms; v <- logpost_bl(from_r(t), y, X, aid, athlete); if (is.finite(v)) -v else 1e12 }
  o <- optim(t0, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12, parscale = c(rep(1, 4), rep(0.1, 3), rep(0.3, length(nms) - 7))))
  o <- optim(o$par, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12, parscale = c(rep(1, 4), rep(0.1, 3), rep(0.3, length(nms) - 7))))
  H <- optimHess(o$par, obj)
  names(o$par) <- nms
  list(t = o$par, par = from_r(o$par), value = -o$value, H = H, nms = nms, conv = o$convergence)
}
draws_bl <- function(f, n = 2000) {
  V <- tryCatch(solve(f$H), error = function(e) MASS::ginv(f$H)); V <- (V + t(V)) / 2
  ev <- eigen(V, symmetric = TRUE); ev$values <- pmax(ev$values, 1e-10)
  Z <- matrix(rnorm(n * length(f$t)), n) %*% t(ev$vectors %*% diag(sqrt(ev$values)))
  T <- sweep(Z, 2, f$t, "+"); colnames(T) <- f$nms
  t(apply(T, 1, function(t) from_r(t)))
}
gpd_fit <- function(x) {   # MLE of the generalized Pareto for exceedances x > 0
  nll <- function(th) { s <- exp(th[1]); k <- th[2]; a <- 1 + k * x / s; if (any(a <= 0)) return(1e10); length(x) * log(s) + (1 / k + 1) * sum(log(a)) }
  o <- optim(c(log(mean(x)), -0.1), nll, hessian = TRUE)
  list(sigma = exp(o$par[1]), xi = o$par[2], V = solve(o$hessian), th = o$par)
}

rows <- list()
final <- read_csv("out/postx_ceiling.csv", show_col_types = FALSE)
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  d <- build_data(att, sex, lift, G = 8)
  best <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  ok <- !is.na(best); y <- best[ok]; X <- d$X[ok, ]; aid <- as.integer(factor(d$u$athlete[ok]))
  fmf <- sprintf("out/extmap_e2_%s_%s_rtau_sm3_mix_fix_log_tau_-3.91202.rds", tolower(sex), lift)
  fm <- fill_par2(readRDS(fmf)$par)
  ## Starting values. The location and scale of the mixture fit refer to its reference component, which
  ## for men is the lower one, so they are a poor start for a unimodal body and can lead the optimizer to
  ## a local optimum with the lower scale at its bound; a least-squares fit of the best lifts on the design
  ## is the second source of starts. Each source is tried with a non-negative and two negative shapes (with
  ## a negative shape some best lifts can lie beyond the endpoint at the start, in which case that start is
  ## dropped), and the best of the fits is kept.
  ls0 <- lm(y ~ X - 1); r0 <- residuals(ls0)
  loc_mix <- fm[c("b0", "b1", "b2", "b3", "g0", "g1", "g2")]
  loc_ls <- c(coef(ls0), log(sd(r0)), 0, 0)
  fit_multi <- function(starts, athlete) {
    fits <- lapply(starts, function(s) tryCatch(fit_bl(s, y, X, aid, athlete), error = function(e) NULL))
    fits <- Filter(function(f) !is.null(f) && is.finite(f$value) && f$value > -1e11, fits)
    if (!length(fits)) stop("no feasible start for ", sex, " ", lift)
    cat(sprintf("   %s fits: log posterior %s\n", if (athlete) "athlete" else "exchangeable",
                paste(sprintf("%.1f", sort(sapply(fits, `[[`, "value"), decreasing = TRUE)), collapse = ", ")))
    fits[[which.max(sapply(fits, `[[`, "value"))]]
  }
  starts <- list()
  for (loc in list(loc_mix, loc_ls)) for (x0 in c(0.05, -0.05, -0.15)) {
    s <- c(loc, xi = x0, kappa = 0, lam = log(1.6), lnu = log(8), rhoa = 0.8); names(s) <- NB
    starts[[length(starts) + 1]] <- s
  }
  fe <- fit_multi(starts, FALSE)
  fa <- fit_multi(c(list(c(fe$par[NB[-12]], rhoa = 0.8)), starts[c(1, 4)]), TRUE)
  cat(sprintf("== %s %s: n = %d days, %d athletes | exch: xi %.3f kappa %.3f lp %.1f | athlete: xi %.3f kappa %.3f rho_a %.3f lp %.1f\n",
              sex, lift, length(y), max(aid), fe$par["xi"], fe$par["kappa"], fe$value, fa$par["xi"], fa$par["kappa"], fa$par["rhoa"], fa$value))
  De <- draws_bl(fe); Da <- draws_bl(fa)
  ## athlete maxima of standardized residuals under the athlete fit
  pa <- fa$par; mu <- as.numeric(X %*% pa[1:4]); om <- as.numeric(exp(X[, 1:3] %*% pa[5:7])); e <- (y - mu) / om
  emax <- as.numeric(tapply(e, aid, max)); u0 <- quantile(emax, 0.9); g <- gpd_fit(emax[emax > u0] - u0)
  gd <- MASS::mvrnorm(2000, g$th, g$V)
  ep_std <- ifelse(gd[, 2] < 0, u0 + exp(gd[, 1]) / (-gd[, 2]), Inf)
  cat(sprintf("   athlete maxima: %d athletes, %d exceedances, xi %.3f (se %.3f)\n", length(emax), sum(emax > u0), g$xi, sqrt(g$V[2, 2])))
  bws <- if (sex == "Men") c(73, 89, 110, 150) else c(59, 71, 86, 120)
  for (bw in bws) {
    xd <- design_at(d$ref, bw, 2026)
    epf <- function(D) apply(D, 1, function(p) { if (p[["xi"]] >= 0) return(Inf); cap_endpoint(full21(p), d$ref, bw, 2026) })
    q9f <- function(D) apply(D, 1, function(p) { pp <- unpack(full21(p)); sum(pp$beta * xd) + omega_of(pp, matrix(xd, 1)) * q_std(0.999, pp$xi, pp$kappa, pp$sL, pp$nuL) })
    ee <- epf(De); ea <- epf(Da)
    mu_b <- sum(pa[1:4] * xd); om_b <- as.numeric(exp(xd[1:3] %*% pa[5:7]))
    em <- mu_b + om_b * ep_std
    fr <- final[final$sex == sex & final$lift == lift & final$bodyweight == bw, ]
    s <- function(v) c(median(v[is.finite(v)]), quantile(v[is.finite(v)], c(.025, .975)), mean(is.finite(v)))
    rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, bodyweight = bw, record = fr$record,
      xi_exch = fe$par[["xi"]], xi_ath = fa$par[["xi"]], xi_max = g$xi, rho_a = fa$par[["rhoa"]],
      q999_exch = median(q9f(De)), q999_ath = median(q9f(Da)), q999_latent = fr$q999_median,
      ceil_exch = s(ee)[1], ceil_exch_lo = s(ee)[2], ceil_exch_hi = s(ee)[3],
      ceil_ath = s(ea)[1], ceil_ath_lo = s(ea)[2], ceil_ath_hi = s(ea)[3],
      ceil_max = s(em)[1], ceil_max_lo = s(em)[2], ceil_max_hi = s(em)[3], p_finite_max = s(em)[4],
      ceil_latent = fr$ceiling_median, ceil_latent_lo = fr$ceiling_lo, ceil_latent_hi = fr$ceiling_hi)
  }
}
SUFQ <- if (abs(Q_SPLICE - 0.9) > 1e-9) sprintf("_q%02d", round(100 * Q_SPLICE)) else ""
out <- bind_rows(rows); write_csv(out, sprintf("out/bestlift_evt%s.csv", SUFQ))
print(as.data.frame(out |> mutate(across(where(is.numeric), ~round(.x, 2)))), width = 250)
