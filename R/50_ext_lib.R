## 50_ext_lib.R
## -----------------------------------------------------------------------
## Two extensions of the marginal-likelihood model (see marginal_kernel.cpp):
##   (1) a shared self-assessment error per athlete-day-lift in the
##       declaration model, parameter lse = log s_eta;
##   (2) an athlete level: the days of one athlete share a Gaussian factor on
##       the copula scale, parameter rhoa in [0, 1).
## Either is switched off by setting its parameter to NA. The parameter
## vector is the base vector plus (lse, rhoa). Units are sorted by athlete
## and the kernel returns one log likelihood per athlete.
source("R/30_marginal_lib.R")

PAR_NAMES2 <- c(PAR_NAMES, "lse", "rhoa", "lsm2", "lsm3", "lpi",
                "mx_m2", "mx_m3", "mx_ls2", "mx_ls3", "mx_a2", "mx_a3", "mx_u", "mx_lsig")
## mx_*: normal-mixture body with a generalized Pareto tail above the estimated
## threshold mx_u (standardized units); first component N(0, 1); active when
## mx_u is finite, in which case kappa, lam and lnu are not used
MIX_NAMES <- c("mx_m2", "mx_m3", "mx_ls2", "mx_ls3", "mx_a2", "mx_a3", "mx_u", "mx_lsig")
## full extended vector from a named vector (missing names are NA)
fill_par2 <- function(x) { v <- rep(NA_real_, length(PAR_NAMES2)); names(v) <- PAR_NAMES2; nm <- intersect(names(x), PAR_NAMES2); v[nm] <- x[nm]; v }
## lpi: logit of the technical-miss probability (NA: no technical misses)
## lsm2, lsm3: situation-specific lower scales of the margin density (NA: equal to lsm, the opener's)
N_THREADS <- as.integer(Sys.getenv("WL_THREADS", "1"))
JAC_DEFAULT <- as.integer(Sys.getenv("WL_JAC", "0"))   # margin parametrisation, see marginal_kernel.cpp
RTAU <- as.integer(Sys.getenv("WL_RTAU", "1"))         # 1 (default): relative attempt scale, P(S) = F(log(Z/w)/tau + delta)
unpack2 <- function(par) {
  p <- unpack(par[PAR_NAMES])
  p$s_eta <- if (is.finite(par["lse"])) exp(par[["lse"]]) else -1
  p$rho_a <- if (is.finite(par["rhoa"])) par[["rhoa"]] else -1
  sm1 <- p$dec[7]
  p$dec <- c(p$dec[1:7], if (is.finite(par["lsm2"])) exp(par[["lsm2"]]) else sm1, if (is.finite(par["lsm3"])) exp(par[["lsm3"]]) else sm1)
  p$pi_miss <- if (is.finite(par["lpi"])) plogis(par[["lpi"]]) else -1
  p$mix <- isTRUE(is.finite(par["mx_u"]))
  if (p$mix) {
    a <- c(0, par[["mx_a2"]], par[["mx_a3"]]); w <- exp(a - max(a)); w <- w / sum(w)
    p$mx <- list(w = w, m = c(0, par[["mx_m2"]], par[["mx_m3"]]), s = c(1, exp(par[["mx_ls2"]]), exp(par[["mx_ls3"]])),
                 u = par[["mx_u"]], sig = exp(par[["mx_lsig"]]))
    p$bmix <- c(3, p$mx$w, p$mx$m, p$mx$s, p$mx$u, p$mx$sig, p$xi)
    p$kappa <- 0; p$sL <- 1; p$nuL <- 10
  } else p$bmix <- numeric(0)
  p
}

## Gauss-Hermite nodes and weights for the weight function exp(-x^2)
gauss_hermite <- function(K) {
  k <- seq_len(K - 1)
  J <- matrix(0, K, K)
  J[cbind(k, k + 1)] <- sqrt(k / 2); J[cbind(k + 1, k)] <- sqrt(k / 2)
  e <- eigen(J, symmetric = TRUE)
  x <- e$values; w <- sqrt(pi) * e$vectors[1, ]^2
  o <- order(x); list(x = x[o], w = w[o])
}

## sort the units of a data object by athlete and add the athlete index
build_data_ext <- function(d, K_gh = 12, K_eta = 8, halfwidth = 6, trunc = TRUE, jac = JAC_DEFAULT) {
  o <- order(d$u$athlete)
  d$u <- d$u[o, ]; d$X <- d$X[o, , drop = FALSE]; d$xstd <- d$xstd[o]
  d$W <- d$W[o, , drop = FALSE]; d$S <- d$S[o, , drop = FALSE]; d$Sint <- d$Sint[o, , drop = FALSE]
  d$D <- d$D[o, , drop = FALSE]; d$has_miss <- d$has_miss[o]; d$has_make <- d$has_make[o]
  d$wmade_min <- d$wmade_min[o]; d$wmiss_max <- d$wmiss_max[o]
  r <- rle(d$u$athlete)
  d$alen <- r$lengths
  d$astart <- c(0L, cumsum(r$lengths)[-length(r$lengths)])
  d$athletes <- r$values
  d$gh <- gauss_hermite(K_gh)
  d$ghe <- gauss_legendre(K_eta); d$inflate <- halfwidth   # eta nodes (Gauss-Legendre per segment) and window half-width in sd
  d$trunc <- as.integer(trunc)   # rule truncation w >= w_prev + 1 of the declaration densities
  d$jac <- as.integer(jac)       # 1: relative margin Z/w - 1 (density Z/w^2 f); 0: log margin log(Z/w) (density f/w)
  d
}

## per-athlete log likelihoods
unit_loglik_ext <- function(par, d) {
  p <- unpack2(par)
  mu <- as.numeric(d$X %*% p$beta); om <- omega_of(p, d$X)
  sc <- splice_consts(p$sL, p$nuL)
  unit_ll_ext_cpp(d$W, d$Sint, d$D, mu, om, p$xi, p$kappa, p$tau, p$delta,
                  sc$zq, sc$ptail, sc$sig0, p$sL, p$nuL, d$gl$x, d$gl$w, REACH_MULT,
                  p$dec, p$s_eta, p$rho_a, as.integer(d$astart), as.integer(d$alen),
                  d$gh$x, d$gh$w, d$ghe$x, d$ghe$w, d$inflate, d$trunc, N_THREADS, d$jac, RTAU, p$pi_miss, p$bmix)
}

log_prior_ext <- function(par, xi_sd = 0.3) {
  pb <- par[PAR_NAMES]
  if (isTRUE(is.finite(par["mx_u"]))) pb[c("kappa", "lam", "lnu")] <- c(0, 0, log(10))
  lp <- log_prior(pb, xi_sd)
  if (!is.finite(lp)) return(-Inf)
  if (RTAU == 1) lp <- lp - dnorm(par[["log_tau"]], log(2), 0.7, log = TRUE) + dnorm(par[["log_tau"]], log(0.025), 0.7, log = TRUE)
  if (is.finite(par["lse"])) {
    if (par[["lse"]] < log(0.002) || par[["lse"]] > log(0.3)) return(-Inf)
    lp <- lp + dnorm(par[["lse"]], log(0.02), 1, log = TRUE)
  }
  if (is.finite(par["rhoa"])) {
    if (par[["rhoa"]] <= 0 || par[["rhoa"]] >= 0.95) return(-Inf)
  }
  for (nm in c("lsm2", "lsm3")) if (is.finite(par[nm])) {
    if (par[[nm]] <= log(0.003) || par[[nm]] >= 0) return(-Inf)
    lp <- lp + dnorm(par[[nm]], log(0.03), 1, log = TRUE)
  }
  if (is.finite(par["lpi"])) {
    if (par[["lpi"]] <= qlogis(0.001) || par[["lpi"]] >= qlogis(0.6)) return(-Inf)
    lp <- lp + dnorm(par[["lpi"]], qlogis(0.1), 1.5, log = TRUE)
  }
  if (isTRUE(is.finite(par["mx_u"]))) {
    m2 <- par[["mx_m2"]]; m3 <- par[["mx_m3"]]; ls <- c(par[["mx_ls2"]], par[["mx_ls3"]])
    ## components 2 and 3 are exchangeable (no ordering): every reported quantity is label-invariant
    if (any(ls <= log(0.05)) || any(ls >= log(20)) || par[["mx_u"]] <= 0.2 || par[["mx_u"]] >= 4 ||
        par[["mx_lsig"]] <= log(0.01) || par[["mx_lsig"]] >= log(10)) return(-Inf)
    lp <- lp + sum(dnorm(c(m2, m3), 0, 3, log = TRUE)) + sum(dnorm(ls, 0, 1, log = TRUE)) +
      sum(dnorm(c(par[["mx_a2"]], par[["mx_a3"]]), 0, 1.5, log = TRUE)) +
      dnorm(par[["mx_u"]], 1.5, 1, log = TRUE) + dnorm(par[["mx_lsig"]], log(0.5), 1, log = TRUE)
  }
  lp
}
log_post_ext <- function(par, d, xi_sd = 0.3) {
  lp <- log_prior_ext(par, xi_sd)
  if (!is.finite(lp)) return(-Inf)
  lp + sum(unit_loglik_ext(par, d))
}

## optimiser on the transformed scale, as fit_map() but for the extended vector
BOUNDS2 <- c(BOUNDS, list(lse = c(log(0.002), log(0.3)), rhoa = c(0, 0.95), lsm2 = c(log(0.003), 0), lsm3 = c(log(0.003), 0), lpi = c(qlogis(0.001), qlogis(0.6)),
                          mx_ls2 = c(log(0.05), log(20)), mx_ls3 = c(log(0.05), log(20)), mx_u = c(0.2, 4), mx_lsig = c(log(0.01), log(10))))
PSCALE2 <- c(PSCALE, lse = 0.3, rhoa = 0.1, lsm2 = PSCALE[["lsm"]], lsm3 = PSCALE[["lsm"]], lpi = 0.3,
             mx_m2 = 0.2, mx_m3 = 0.2, mx_ls2 = 0.1, mx_ls3 = 0.1, mx_a2 = 0.2, mx_a3 = 0.2, mx_u = 0.05, mx_lsig = 0.1)
to_real2 <- function(par) {
  t <- par
  for (nm in intersect(names(BOUNDS2), names(par)[is.finite(par)])) {
    b <- BOUNDS2[[nm]]
    t[nm] <- qlogis(pmin(pmax((par[nm] - b[1]) / (b[2] - b[1]), 1e-9), 1 - 1e-9))
  }
  t
}
from_real2 <- function(t) {
  par <- t
  for (nm in intersect(names(BOUNDS2), names(t)[is.finite(t)])) {
    b <- BOUNDS2[[nm]]
    par[nm] <- b[1] + (b[2] - b[1]) * plogis(t[nm])
  }
  par
}
fit_map_ext <- function(d, start, fixed = NULL, maxit = 400) {
  free <- setdiff(PAR_NAMES2, names(fixed))
  full <- start[PAR_NAMES2]; full[names(fixed)] <- fixed
  free <- free[is.finite(full[free])]
  tfull <- to_real2(full)
  obj <- function(th) {
    t <- tfull; t[free] <- th
    v <- log_post_ext(from_real2(t), d)
    if (!is.finite(v)) 1e12 else -v
  }
  ps <- PSCALE2[free]; ps[names(ps) %in% names(BOUNDS2)] <- 1
  o <- optim(tfull[free], obj, method = "BFGS",
             control = list(parscale = ps, maxit = maxit, reltol = 1e-10))
  t <- tfull; t[free] <- o$par
  par <- from_real2(t)
  list(par = par, value = -o$value, conv = o$convergence,
       loglik = sum(unit_loglik_ext(par, d)), evals = o$counts)
}

## Laplace covariance with fixed steps, extended vector
laplace_cov_ext <- function(par, d, free, step_div = 50) {
  f <- function(p) { v <- log_post_ext(p, d); if (is.finite(v)) -v else NA }
  k <- length(free); h <- PSCALE2[free] / step_div
  H <- matrix(NA_real_, k, k); f0 <- f(par); fp <- fm <- numeric(k)
  for (i in 1:k) {
    pp <- par; pp[free[i]] <- pp[free[i]] + h[i]; fp[i] <- f(pp)
    pm <- par; pm[free[i]] <- pm[free[i]] - h[i]; fm[i] <- f(pm)
    H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h[i]^2
  }
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    p1 <- par; p1[free[i]] <- p1[free[i]] + h[i]; p1[free[j]] <- p1[free[j]] + h[j]
    p2 <- par; p2[free[i]] <- p2[free[i]] - h[i]; p2[free[j]] <- p2[free[j]] - h[j]
    H[i, j] <- H[j, i] <- (f(p1) - fp[i] - fp[j] + 2 * f0 - fm[i] - fm[j] + f(p2)) / (2 * h[i] * h[j])
  }
  H[is.na(H)] <- 0
  Hs <- diag(PSCALE2[free]) %*% ((H + t(H)) / 2) %*% diag(PSCALE2[free])
  e <- eigen(Hs, symmetric = TRUE); lam <- e$values; fl <- 1e-4 * max(lam)
  nfl <- sum(lam < fl); lam[lam < fl] <- fl
  V <- diag(PSCALE2[free]) %*% (e$vectors %*% diag(1 / lam, k) %*% t(e$vectors)) %*% diag(PSCALE2[free])
  dimnames(V) <- list(free, free); attr(V, "n_floored") <- nfl; V
}

## ---- R reference implementations for checking the kernel ---------------
## shared-error declaration factor by brute force on an eta grid, for one unit
ref_eta_integral <- function(z, w, sit, dec, s_eta, wprev, trunc = TRUE, jac = 1, n_grid = 40001, half = 8) {
  eta <- seq(-half * s_eta, half * s_eta, length.out = n_grid); d_eta <- eta[2] - eta[1]
  lg <- rep(0, n_grid)
  for (j in seq_along(w)) {
    m <- dec[sit[j]]; sp <- dec[3 + sit[j]]; sm <- if (length(dec) >= 9) dec[6 + sit[j]] else dec[7]
    tp <- function(r, m, sp) { cc <- 2 / (sm + sp); ifelse(r < m, cc * sm * pnorm((r - m) / sm), cc * sm / 2 + cc * sp * (pnorm((r - m) / sp) - 0.5)) }
    if (jac == 1) { zest <- z * (1 + eta); r <- zest / w[j] - 1; lj <- log(zest) - 2 * log(w[j]); rmax <- zest / (wprev[j] + 1) - 1 }
    else { r <- log(z / w[j]) + eta; lj <- -log(w[j]); rmax <- log(z / (wprev[j] + 1)) + eta }
    u <- ifelse(r < m, (r - m) / sm, (r - m) / sp)
    lg <- lg + log(2 / (sm + sp)) - 0.9189385332046727 - 0.5 * u^2 + lj
    if (trunc && sit[j] >= 2) { lg <- lg - log(pmax(tp(rmax, m, sp), 1e-300)); lg[r > rmax + 1e-9] <- -Inf }
  }
  lg <- lg + dnorm(eta, 0, s_eta, log = TRUE)
  mx <- max(lg); mx + log(sum(exp(lg - mx)) * d_eta)
}


## ---- two-lift copula stage through the extended marginal model ----------
## par_s, par_c: extended parameter vectors of the two lifts (PAR_NAMES2);
## ds, dc: data objects with the same units in the same order (pair_data);
## returns unit joint log likelihoods under the copula (rho, psi) and,
## optionally, one draw of the pair of capabilities per unit.
lift_pack_ext <- function(par) {
  p <- unpack2(par); sc <- splice_consts(p$sL, p$nuL)
  list(vec = c(p$xi, p$kappa, p$tau, p$delta[2], p$delta[3], sc$zq, sc$sig0, p$sL, sc$ptail, p$nuL, p$pi_miss), dec = p$dec, p = p, bmix = p$bmix)
}
joint_loglik_ext <- function(par_s, par_c, ds, dc, rho, psi, draw = FALSE, G = 24, jac = JAC_DEFAULT) {
  ls <- lift_pack_ext(par_s); lc <- lift_pack_ext(par_c)
  mus <- as.numeric(ds$X %*% ls$p$beta); oms <- omega_of(ls$p, ds$X)
  muc <- as.numeric(dc$X %*% lc$p$beta); omc <- omega_of(lc$p, dc$X)
  joint_ll_ext_cpp(ds$W, ds$Sint, ds$D, dc$W, dc$Sint, dc$D, mus, oms, muc, omc, ls$vec, lc$vec, ls$dec, lc$dec,
                   as.integer(G), REACH_MULT, as.integer(jac), RTAU, 1L, rho, psi, draw, N_THREADS, ls$bmix, lc$bmix)
}

## posterior mean, sd and exceedance probability of the capability of each unit (athlete factor off)
latent_summary_ext <- function(par, d, thr = rep(-Inf, d$n), G = d$G, jac = d$jac) {
  lp <- lift_pack_ext(par)
  mu <- as.numeric(d$X %*% lp$p$beta); om <- omega_of(lp$p, d$X)
  unit_latent_ext_cpp(d$W, d$Sint, d$D, mu, om, lp$vec, lp$dec, as.integer(G), REACH_MULT, as.integer(jac), RTAU, 1L, as.numeric(thr), N_THREADS, lp$bmix)
}


## ---- quantiles and endpoint of the standardized error, either body ------
mix_cdf <- function(e, mx, xi) {
  Su <- sum(mx$w * pnorm(mx$u, mx$m, mx$s, lower.tail = FALSE))
  sapply(e, function(x) if (x <= mx$u) 1 - sum(mx$w * pnorm(x, mx$m, mx$s, lower.tail = FALSE))
         else 1 - Su * max(1 + xi * (x - mx$u) / mx$sig, 0)^(-1 / xi))
}
std_quantile <- function(p, u) {
  if (!isTRUE(p$mix)) return(q_std(u, p$xi, p$kappa, p$sL, p$nuL))
  mx <- p$mx; Su <- sum(mx$w * pnorm(mx$u, mx$m, mx$s, lower.tail = FALSE))
  out <- numeric(length(u)); hi <- u > 1 - Su
  if (any(hi)) out[hi] <- mx$u + mx$sig / p$xi * (((1 - u[hi]) / Su)^(-p$xi) - 1)
  if (any(!hi)) {   # body: invert the mixture cdf on a fine grid, then refine by one Newton step
    eg <- seq(min(mx$m - 10 * mx$s), mx$u, length.out = 40001)
    Fg <- 1 - colSums(mx$w * sapply(eg, function(x) pnorm(x, mx$m, mx$s, lower.tail = FALSE)))
    x <- approx(Fg, eg, xout = u[!hi], ties = "ordered", rule = 2)$y
    fx <- colSums(mx$w * sapply(x, function(z) dnorm(z, mx$m, mx$s)))
    Fx <- 1 - colSums(mx$w * sapply(x, function(z) pnorm(z, mx$m, mx$s, lower.tail = FALSE)))
    out[!hi] <- x - (Fx - u[!hi]) / pmax(fx, 1e-300)
  }
  out
}
std_endpoint <- function(p) {
  if (p$xi >= 0) return(Inf)
  if (isTRUE(p$mix)) return(p$mx$u + p$mx$sig / (-p$xi))
  sc <- splice_consts(p$sL, p$nuL); sc$zq + sc$sig0 * exp(p$kappa) / (-p$xi)
}
## threshold probability H(u) of the mixture body (the share of athlete-days below the tail)
mix_Hu <- function(p) if (isTRUE(p$mix)) 1 - sum(p$mx$w * pnorm(p$mx$u, p$mx$m, p$mx$s, lower.tail = FALSE)) else NA
