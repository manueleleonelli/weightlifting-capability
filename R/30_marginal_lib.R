## 30_marginal_lib.R
## -----------------------------------------------------------------------
## Marginal likelihood of the latent capability model, with the capabilities
## integrated out by Gauss-Legendre quadrature, unit by unit. No data
## augmentation, so the global parameters can be sampled or optimised
## directly. One lift at a time (the copula is handled in a second stage).
##
## Model (per lift l, unit i = athlete-day):
##   Z_i = mu(x_i) + omega(x_i) e_i,  e_i ~ spliced N(0,1) / GPD(xi, sigma0 e^kappa)
##   P(S_ij = 1 | Z_i) = plogis((Z_i - w_ij)/tau + delta_j),  delta_1 = 0
##
## Functions
##   prep_units(att, sex, lift)     unit table with weights, outcomes, design
##   gauss_legendre(G)              nodes and weights on [-1, 1]
##   unit_loglik_marg(par, d)       vector of unit log marginal likelihoods
##   log_post(par, d)               log posterior (sum + prior)
##   unpack(par)                    named list from the parameter vector

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})

Q_SPLICE <- as.numeric(Sys.getenv("WL_Q", "0.90"))
ZQ <- qnorm(Q_SPLICE)
SIG0 <- (1 - Q_SPLICE) / dnorm(ZQ)

## ---- quadrature ------------------------------------------------------
gauss_legendre <- function(G) {
  k <- seq_len(G - 1)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, G, G)
  J[cbind(k, k + 1)] <- b; J[cbind(k + 1, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  x <- e$values; w <- 2 * e$vectors[1, ]^2
  o <- order(x)
  list(x = x[o], w = w[o])
}

## ---- data ------------------------------------------------------------
## Cleaning: units whose attempt weights are not non-decreasing, whose weight
## does not increase after a successful attempt (which the rules require), or
## whose largest declared weight exceeds 1.5 times the smallest, are
## transcription errors in the source sheets and are dropped (a few dozen in
## 78,412).
prep_units <- function(att, sex, lift_name, design = NULL) {
  base <- att |>
    filter(!is_youth, outcome %in% c("good", "fail"), weight > 0,
           !is.na(bodyweight), bodyweight > 30, !is.na(year),
           gender == sex, lift == lift_name) |>
    mutate(success = as.integer(outcome == "good"),
           unit = paste(event_id, athlete, category, sep = "|")) |>
    arrange(unit, attempt_no)
  u <- base |>
    group_by(unit) |>
    summarise(athlete = first(athlete), event_id = first(event_id),
              bodyweight = first(bodyweight), year = first(year),
              m = n(),
              w1 = weight[1], w2 = weight[2], w3 = weight[3],
              s1 = success[1], s2 = success[2], s3 = success[3],
              wmin = min(weight), wmax = max(weight),
              mono = all(diff(weight) >= 0),
              rule = all(diff(weight) > 0 | success[-n()] == 0),
              .groups = "drop")
  ## rule: after a successful attempt the next declared weight must increase
  bad <- !u$mono | !u$rule | u$wmax / u$wmin > 1.5
  u <- u[!bad, ]
  attr(u, "n_dropped") <- sum(bad)
  u
}

## Design matrix as in the paper: standardised log bodyweight, its
## orthogonalised square, and centred year. The centring constants are
## stored so that predictions at new covariate values use the same mapping.
make_design <- function(u, ref = NULL) {
  lb <- log(u$bodyweight)
  if (is.null(ref)) {
    centre <- mean(lb); scale <- sd(lb)
    x <- (lb - centre) / scale
    a2 <- mean(x^2); b2 <- sum((x^2 - a2) * x) / sum(x * x)
    raw2 <- x^2 - a2 - b2 * x; c2 <- sd(raw2)
    tbar <- mean(u$year - 2012)
    ref <- list(centre = centre, scale = scale, a2 = a2, b2 = b2, c2 = c2,
                tbar = tbar)
  }
  x <- (lb - ref$centre) / ref$scale
  X <- cbind(1, x, (x^2 - ref$a2 - ref$b2 * x) / ref$c2,
             (u$year - 2012) - ref$tbar)
  list(X = X, xstd = x, ref = ref)
}

design_at <- function(ref, bw, year) {
  x <- (log(bw) - ref$centre) / ref$scale
  c(1, x, (x^2 - ref$a2 - ref$b2 * x) / ref$c2, (year - 2012) - ref$tbar)
}

## Assemble everything the likelihood needs for one sex and lift.
build_data <- function(att, sex, lift_name, G = 24) {
  u <- prep_units(att, sex, lift_name)
  des <- make_design(u)
  gl <- gauss_legendre(G)
  W <- as.matrix(u[, c("w1", "w2", "w3")])
  S <- as.matrix(u[, c("s1", "s2", "s3")])
  has_miss <- rowSums(S == 0, na.rm = TRUE) > 0
  has_make <- rowSums(S == 1, na.rm = TRUE) > 0
  wmade_min <- suppressWarnings(apply(ifelse(S == 1, W, NA), 1, min, na.rm = TRUE))
  wmiss_max <- suppressWarnings(apply(ifelse(S == 0, W, NA), 1, max, na.rm = TRUE))
  Sint <- S; storage.mode(Sint) <- "integer"; Sint[is.na(Sint)] <- -1L
  ## declaration situations: 1 opener; 2 raise after a make; 3 raise after a
  ## miss; 0 repeat of a missed weight (no new declaration)
  Dm <- matrix(0L, nrow(W), 3)
  Dm[, 1] <- 1L
  for (j in 2:3) {
    h <- !is.na(W[, j])
    raise <- h & W[, j] > W[, j - 1]
    Dm[h & raise & S[, j - 1] == 1, j] <- 2L
    Dm[h & raise & S[, j - 1] == 0, j] <- 3L
  }
  list(n = nrow(u), u = u, X = des$X, xstd = des$xstd, ref = des$ref,
       W = W, S = S, Sint = Sint, D = Dm, has_miss = has_miss, has_make = has_make,
       wmade_min = wmade_min, wmiss_max = wmiss_max,
       gl = gl, G = G, sex = sex, lift = lift_name,
       n_dropped = attr(u, "n_dropped"))
}

## ---- parameter vector ------------------------------------------------
## par = (beta[1:4], gamma[1:2], xi, kappa, delta2, delta3, log_tau, lam,
##        m1, m2, m3, lsp1, lsp2, lsp3, lsm)
## lam = log of the left-hand scale of the two-piece normal bulk (right-hand
## scale 1); lam = 0 is the symmetric normal bulk of the original model.
## m_s, lsp_s: mode and log upper scale of the two-piece normal distribution
## of the relative margin (Z - w)/w of a declared weight in situation s
## (1 opener, 2 raise after a make, 3 raise after a miss); lsm: log of the
## shared lower scale. The declaration model is off when lsm is NA.
## lnu = log of the degrees of freedom of the half-t left side of the bulk.
PAR_NAMES <- c("b0", "b1", "b2", "b3", "g0", "g1", "g2", "xi", "kappa",
               "delta2", "delta3", "log_tau", "lam",
               "m1", "m2", "m3", "lsp1", "lsp2", "lsp3", "lsm", "lnu")
DEC_NAMES <- c("m1", "m2", "m3", "lsp1", "lsp2", "lsp3", "lsm")
unpack <- function(par) {
  if (length(par) < 13) par <- c(par, 0)
  op <- length(par) >= 20 && is.finite(par[20])
  dec <- if (op) c(par[14:16], exp(par[17:19]), exp(par[20])) else c(rep(0, 6), -1)
  nuL <- if (length(par) >= 21 && is.finite(par[21])) exp(par[21]) else 1e6
  list(beta = par[1:4], gamma = par[5:7], xi = par[8], kappa = par[9],
       delta = c(0, par[10], par[11]), tau = exp(par[12]), sL = exp(par[13]),
       dec = dec, nuL = nuL)
}
## scale of the upper half of the body: log-quadratic in bodyweight, using
## the same centred and orthogonalised terms as the location
omega_of <- function(p, X) as.numeric(exp(X[, 1:3, drop = FALSE] %*% p$gamma))

## Splice constants. The generalized Pareto starts at the Q_SPLICE quantile of
## the right half-normal, a fixed number of right-hand scales above the mode,
## so that the threshold does not move with the left side of the bulk. Its
## mass is c (1 - Q_SPLICE) / 2 with c the density level at the mode, and the
## continuity scale follows.
bulk_const <- function(sL, nuL) {
  t0 <- exp(lgamma((nuL + 1) / 2) - lgamma(nuL / 2)) / sqrt(nuL * pi)
  2 / (1 + sL * dnorm(0) / t0)
}
splice_consts <- function(sL, nuL = 1e6, q = Q_SPLICE) {
  zq <- qnorm((1 + q) / 2)
  cc <- bulk_const(sL, nuL)
  ptail <- cc * (1 - q) / 2
  f_zq <- cc * dnorm(zq)
  list(zq = zq, sig0 = ptail / f_zq, ptail = ptail, cc = cc)
}

## ---- marginal density pieces ----------------------------------------
## log density of the standardised error, vectorised over a matrix e
log_f_std <- function(e, xi, kappa, sL = 1, nuL = 1e6) {
  sc <- splice_consts(sL, nuL); ZQ <- sc$zq; sig <- sc$sig0 * exp(kappa)
  lc <- log(sc$cc)
  out <- ifelse(e < 0, lc + dnorm(0, log = TRUE) - ((nuL + 1) / 2) * log1p((e / sL)^2 / nuL),
                lc + dnorm(e, log = TRUE))
  hi <- e > ZQ
  if (any(hi)) {
    y <- (e[hi] - ZQ) / sig
    if (abs(xi) < 1e-6) {
      v <- log(sc$ptail) - log(sig) - y
    } else {
      a <- 1 + xi * y
      v <- log(sc$ptail) - log(sig) - (1 / xi + 1) * log(pmax(a, 1e-300))
      v[a <= 0] <- -Inf
    }
    out[hi] <- v
  }
  out
}

## survival function of the standardised error
surv_std <- function(e, xi, kappa, sL = 1, nuL = 1e6) {
  sc <- splice_consts(sL, nuL); ZQ <- sc$zq; sig <- sc$sig0 * exp(kappa)
  t0 <- exp(lgamma((nuL + 1) / 2) - lgamma(nuL / 2)) / sqrt(nuL * pi)
  left_mass <- sc$cc * dnorm(0) * sL / t0
  out <- ifelse(e < 0, 1 - left_mass * pt(e / sL, nuL), sc$cc * pnorm(e, lower.tail = FALSE))
  hi <- e > ZQ
  if (any(hi)) {
    y <- (e[hi] - ZQ) / sig
    if (abs(xi) < 1e-6) {
      v <- sc$ptail * exp(-y)
    } else {
      a <- 1 + xi * y
      v <- sc$ptail * pmax(a, 0)^(-1 / xi)
    }
    out[hi] <- v
  }
  out
}

## quantile of the standardised error for u > Q_SPLICE (and below, normal)
q_std <- function(u, xi, kappa, sL = 1, nuL = 1e6) {
  sc <- splice_consts(sL, nuL); ZQ <- sc$zq; sig <- sc$sig0 * exp(kappa)
  t0 <- exp(lgamma((nuL + 1) / 2) - lgamma(nuL / 2)) / sqrt(nuL * pi)
  left_mass <- sc$cc * dnorm(0) * sL / t0
  u_tail <- 1 - sc$ptail
  uu <- pmin(u, u_tail)
  out <- ifelse(uu < left_mass, sL * qt(uu / left_mass, nuL),
                qnorm(1 - (1 - uu) / sc$cc))
  hi <- u > u_tail
  if (any(hi)) {
    p <- (1 - u[hi]) / sc$ptail
    if (abs(xi) < 1e-6) out[hi] <- ZQ - sig * log(p)
    else out[hi] <- ZQ + sig * (p^(-xi) - 1) / xi
  }
  out
}

## ---- attempt likelihood on a matrix of capabilities ------------------
## Zm: n x K matrix of candidate capabilities; returns n x K log A
log_attempt <- function(Zm, W, S, tau, delta) {
  out <- matrix(0, nrow(Zm), ncol(Zm))
  for (j in 1:3) {
    h <- !is.na(W[, j])
    if (!any(h)) next
    eta <- (Zm[h, , drop = FALSE] - W[h, j]) / tau + delta[j]
    sgn <- ifelse(S[h, j] == 1, 1, -1)
    out[h, ] <- out[h, ] + plogis(sgn * eta, log.p = TRUE)
  }
  out
}

## ---- unit log marginal likelihoods ----------------------------------
## For each unit the integration range is [a, b]:
##   a = lightest made weight - 12 tau   (or mu - 7 omega if nothing was made)
##   b = heaviest missed weight + 12 tau (or heaviest made weight + 12 tau,
##       with the analytic survival beyond b added, if nothing was missed)
## and the range is split at the splice point whenever it falls inside,
## with G Gauss-Legendre nodes on each piece.
unit_loglik_marg <- function(par, d, return_parts = FALSE) {
  p <- unpack(par)
  mu <- as.numeric(d$X %*% p$beta)
  om <- omega_of(p, d$X)
  tau <- p$tau
  sc <- splice_consts(p$sL, p$nuL); ZQ <- sc$zq
  reach <- 12 * tau
  a <- ifelse(d$has_make, d$wmade_min - reach, mu - 7 * om * p$sL)
  b <- ifelse(d$has_miss, d$wmiss_max + reach, d$W[cbind(seq_len(d$n), d$u$m)] + reach)
  if (p$dec[7] > 0) {
    any_dec <- rep(FALSE, d$n)
    for (j in 1:3) {
      h <- !is.na(d$W[, j]) & d$D[, j] > 0
      sit <- d$D[, j]; m <- p$dec[pmax(sit, 1)]; sp <- p$dec[3 + pmax(sit, 1)]; sm <- p$dec[7]
      lo_h <- d$W[, j] * (1 + m - 8 * sm); hi_h <- d$W[, j] * (1 + m + 8 * sp)
      a <- ifelse(h, pmax(a, lo_h), a)
      b <- ifelse(h & !d$has_miss & !any_dec, hi_h, ifelse(h, pmin(b, hi_h), b))
      any_dec <- any_dec | h
    }
  }
  b <- pmax(b, a + 1e-3)
  zs <- mu + om * ZQ
  ## split at the splice point when it lies inside the range, otherwise at
  ## the midpoint so that both sets of nodes are used
  m <- ifelse(zs > a & zs < b, zs, (a + b) / 2)
  nodes <- d$gl$x; wts <- d$gl$w
  ## piece 1: [a, m]; piece 2: [m, b]
  h1 <- (m - a) / 2; c1 <- (m + a) / 2
  h2 <- (b - m) / 2; c2 <- (b + m) / 2
  Z1 <- outer(c1, rep(1, d$G)) + outer(h1, nodes)
  Z2 <- outer(c2, rep(1, d$G)) + outer(h2, nodes)
  Zm <- cbind(Z1, Z2)
  E <- (Zm - mu) / om
  lf <- log_f_std(E, p$xi, p$kappa, p$sL, p$nuL) - log(om)
  la <- log_attempt(Zm, d$W, d$S, tau, p$delta)
  if (p$dec[7] > 0) {
    sm <- p$dec[7]
    tpn_cdf <- function(r, m, sp) {
      cc <- 2 / (sm + sp)
      ifelse(r < m, cc * sm * pnorm((r - m) / sm), cc * sm / 2 + cc * sp * (pnorm((r - m) / sp) - 0.5))
    }
    for (j in 1:3) {
      h <- which(!is.na(d$W[, j]) & d$D[, j] > 0)
      if (!length(h)) next
      sit <- d$D[h, j]; m <- p$dec[sit]; sp <- p$dec[3 + sit]; w <- d$W[h, j]
      r <- (Zm[h, , drop = FALSE] - w) / w
      u <- ifelse(r < m, (r - m) / sm, (r - m) / sp)
      v <- log(2 / (sm + sp)) - 0.9189385332046727 - log(w) - 0.5 * u^2
      if (j > 1) {
        wprev <- d$W[h, j - 1]
        rmax <- (Zm[h, , drop = FALSE] - wprev - 1) / (wprev + 1)
        Gm <- tpn_cdf(rmax, m, sp)
        v <- v - log(pmax(Gm, 1e-300))
        v[r > rmax] <- -Inf
      }
      la[h, ] <- la[h, ] + v
    }
  }
  lw <- cbind(outer(h1, wts), outer(h2, wts))
  lw <- log(pmax(lw, 1e-300))
  M <- lf + la + lw
  mx <- apply(M, 1, max)
  mx[!is.finite(mx)] <- 0
  integ <- exp(mx) * rowSums(exp(M - mx))
  ## tail mass beyond b for units with no miss (attempt likelihood ~ 1 there)
  tail_mass <- ifelse(d$has_miss | p$dec[7] > 0, 0, surv_std((b - mu) / om, p$xi, p$kappa, p$sL, p$nuL))
  L <- integ + tail_mass
  ll <- log(pmax(L, 1e-300))
  if (return_parts) return(list(ll = ll, integ = integ, tail = tail_mass, a = a, b = b))
  ll
}

## ---- prior and posterior --------------------------------------------
XI_LO <- -1.5; XI_HI <- 0.4
log_prior <- function(par, xi_sd = 0.3, tau_fixed = NULL) {
  p <- unpack(par)
  if (p$xi <= XI_LO || p$xi >= XI_HI || abs(p$kappa) > 1.5 || abs(log(p$sL)) > 2) return(-Inf)
  lp <- sum(dnorm(p$beta, 0, 300, log = TRUE)) +
    dnorm(p$gamma[1], log(15), 0.7, log = TRUE) +
    dnorm(p$gamma[2], 0, 1, log = TRUE) + dnorm(p$gamma[3], 0, 1, log = TRUE) +
    dnorm(p$xi, 0, xi_sd, log = TRUE) +
    dnorm(p$kappa, 0, 0.35, log = TRUE) +
    sum(dnorm(p$delta[2:3], 0, 1.5, log = TRUE))
  if (is.null(tau_fixed)) lp <- lp + dnorm(par[12], log(2), 0.7, log = TRUE)
  if (length(par) >= 13) lp <- lp + dnorm(par[13], 0, 1, log = TRUE)
  if (length(par) >= 21 && is.finite(par[21])) {
    if (par[21] < log(2) || par[21] > log(200)) return(-Inf)
    lp <- lp + dnorm(par[21], log(10), 1.5, log = TRUE)
  }
  if (length(par) >= 20 && is.finite(par[20])) {
    if (any(par[14:16] <= -0.2) || any(par[14:16] >= 0.3) ||
        any(par[17:20] <= log(0.003)) || any(par[17:20] >= 0)) return(-Inf)
    lp <- lp + sum(dnorm(par[14:16], 0.03, 0.1, log = TRUE)) +
      sum(dnorm(par[17:20], log(0.03), 1, log = TRUE))
  }
  lp
}

log_post <- function(par, d, xi_sd = 0.3) {
  lp <- log_prior(par, xi_sd)
  if (!is.finite(lp)) return(-Inf)
  lp + sum(unit_loglik_marg(par, d))
}

## ---- derived quantities ---------------------------------------------
## capability quantile at covariate value xd for probability u
cap_quantile <- function(par, ref, bw, year, u) {
  p <- unpack(par)
  xd <- design_at(ref, bw, year)
  mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
  mu + om * q_std(u, p$xi, p$kappa, p$sL, p$nuL)
}
cap_endpoint <- function(par, ref, bw, year) {
  p <- unpack(par)
  if (p$xi >= 0) return(Inf)
  xd <- design_at(ref, bw, year)
  mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
  sc <- splice_consts(p$sL, p$nuL)
  mu + om * (sc$zq + sc$sig0 * exp(p$kappa) / (-p$xi))
}

## ---- compiled kernel --------------------------------------------------
## The C++ kernel in marginal_kernel.cpp implements the same quadrature as
## unit_loglik_marg() about twenty times faster. Build it once per session.
if (requireNamespace("Rcpp", quietly = TRUE)) {
  ## compiled once, reused by later sessions. Every call of sourceCpp rewrites the index of the
  ## cache (rcpp_cache/.../cache.rds), so two scripts started at the same moment can find it
  ## empty or half written; the call is then repeated after a short pause.
  for (k in 1:6) {
    err <- tryCatch({ Rcpp::sourceCpp("R/marginal_kernel.cpp", cacheDir = Sys.getenv("WL_RCPP_CACHE", "rcpp_cache")); NULL },
                    error = function(e) e)
    if (is.null(err)) break
    if (k == 6 || !grepl("zero-byte|error reading|cache|restore file|load", conditionMessage(err))) stop(err)
    Sys.sleep(runif(1, 1, 3))
  }
  REACH_MULT <- 12
  unit_loglik_fast <- function(par, d) {
    p <- unpack(par)
    mu <- as.numeric(d$X %*% p$beta)
    om <- omega_of(p, d$X)
    sc <- splice_consts(p$sL, p$nuL)
    unit_ll_cpp(d$W, d$Sint, mu, om, p$xi, p$kappa, p$tau, p$delta,
                sc$zq, sc$ptail, sc$sig0, p$sL, p$nuL, d$gl$x, d$gl$w, REACH_MULT,
                d$D, p$dec)
  }
  log_post_fast <- function(par, d, xi_sd = 0.3) {
    lp <- log_prior(par, xi_sd)
    if (!is.finite(lp)) return(-Inf)
    lp + sum(unit_loglik_fast(par, d))
  }
  latent_summary <- function(par, d, probs = c(0.025, 0.5, 0.975)) {
    p <- unpack(par)
    mu <- as.numeric(d$X %*% p$beta)
    om <- omega_of(p, d$X)
    sc <- splice_consts(p$sL, p$nuL)
    unit_latent_cpp(d$W, d$Sint, mu, om, p$xi, p$kappa, p$tau, p$delta,
                    sc$zq, sc$ptail, sc$sig0, p$sL, p$nuL, d$gl$x, d$gl$w, REACH_MULT, probs,
                    d$D, p$dec)
  }
}

## ---- bounded reparametrisation for optimisation ----------------------
## xi, kappa and lam are mapped to the real line through a logistic map on
## their prior supports, so that an optimiser never steps outside them. The
## MAP is unchanged (the map is monotone); the posterior density is evaluated
## in the original parametrisation, without a Jacobian.
BOUNDS <- list(xi = c(XI_LO, XI_HI), kappa = c(-1.5, 1.5), lam = c(-2, 2),
               m1 = c(-0.2, 0.3), m2 = c(-0.2, 0.3), m3 = c(-0.2, 0.3),
               lsp1 = c(log(0.003), 0), lsp2 = c(log(0.003), 0), lsp3 = c(log(0.003), 0),
               lsm = c(log(0.003), 0), lnu = c(log(2), log(200)))
to_real <- function(par) {
  t <- par
  for (nm in intersect(names(BOUNDS), names(par)[is.finite(par)])) {
    b <- BOUNDS[[nm]]
    t[nm] <- qlogis(pmin(pmax((par[nm] - b[1]) / (b[2] - b[1]), 1e-9), 1 - 1e-9))
  }
  t
}
from_real <- function(t) {
  par <- t
  for (nm in intersect(names(BOUNDS), names(t)[is.finite(t)])) {
    b <- BOUNDS[[nm]]
    par[nm] <- b[1] + (b[2] - b[1]) * plogis(t[nm])
  }
  par
}

## ---- bivariate stage wrapper -----------------------------------------
## par_s, par_c: full parameter vectors of the two lifts; ds, dc: their data
## objects (same units in the same order). Returns the vector of unit joint
## log likelihoods under the t copula (rho, psi), and optionally one
## posterior draw of the two capabilities per unit.
lift_pack <- function(par) {
  p <- unpack(par); sc <- splice_consts(p$sL, p$nuL)
  c(p$xi, p$kappa, p$tau, p$delta[2], p$delta[3], sc$zq, sc$sig0, p$sL, sc$ptail, p$nuL)
}
joint_loglik <- function(par_s, par_c, ds, dc, rho, psi, draw = FALSE) {
  p_s <- unpack(par_s); p_c <- unpack(par_c)
  mus <- as.numeric(ds$X %*% p_s$beta); oms <- omega_of(p_s, ds$X)
  muc <- as.numeric(dc$X %*% p_c$beta); omc <- omega_of(p_c, dc$X)
  joint_ll_cpp(ds$W, ds$Sint, dc$W, dc$Sint, ds$D, dc$D, mus, oms, muc, omc,
               lift_pack(par_s), lift_pack(par_c), p_s$dec, p_c$dec,
               0, ds$gl$x, ds$gl$w, REACH_MULT, rho, psi, draw)
}
## pair the two lifts of each athlete-day
pair_data <- function(ds, dc) {
  key <- intersect(ds$u$unit, dc$u$unit)
  sub <- function(d, idx) {
    d$u <- d$u[idx, ]; d$X <- d$X[idx, , drop = FALSE]; d$xstd <- d$xstd[idx]
    d$W <- d$W[idx, , drop = FALSE]; d$S <- d$S[idx, , drop = FALSE]; d$Sint <- d$Sint[idx, , drop = FALSE]
    d$D <- d$D[idx, , drop = FALSE]
    d$has_miss <- d$has_miss[idx]; d$has_make <- d$has_make[idx]
    d$wmade_min <- d$wmade_min[idx]; d$wmiss_max <- d$wmiss_max[idx]; d$n <- length(idx); d
  }
  is <- match(key, ds$u$unit); ic <- match(key, dc$u$unit)
  list(s = sub(ds, is), c = sub(dc, ic))
}

## ---- MAP optimiser ---------------------------------------------------
## Box-constrained quasi-Newton (PORT) on the natural parameters; bounds
## slightly inside the prior supports. 'fixed' is a named vector of
## parameters held at the given values (NA entries switch a component off).
PSCALE <- c(100, 10, 5, 1, 1, 0.3, 0.3, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, rep(0.02, 3), rep(0.3, 4), 0.5)
names(PSCALE) <- PAR_NAMES
LOWER <- c(rep(-Inf, 7), XI_LO + 1e-3, -1.5 + 1e-3, -6, -6, log(0.3), -2 + 1e-3,
           rep(-0.2 + 1e-3, 3), rep(log(0.003) + 1e-3, 4), log(2) + 1e-3)
UPPER <- c(rep(Inf, 7), XI_HI - 1e-3, 1.5 - 1e-3, 6, 6, log(8), 2 - 1e-3,
           rep(0.3 - 1e-3, 3), rep(0 - 1e-3, 4), log(200) - 1e-3)
names(LOWER) <- names(UPPER) <- PAR_NAMES
fit_map <- function(d, start, fixed = NULL, xi_sd = 0.3, maxit = 500, method = "BFGS") {
  free <- setdiff(PAR_NAMES, names(fixed))
  full <- start[PAR_NAMES]; full[names(fixed)] <- fixed
  if (method == "BFGS") {
    ## quasi-Newton on the parameters mapped to the real line; the finite
    ## difference step of optim (one thousandth of the scale) is large enough
    ## to average over the small discontinuities of the quadrature rule
    tfull <- to_real(full)
    obj <- function(th) {
      t <- tfull; t[free] <- th
      v <- log_post_fast(from_real(t), d, xi_sd)
      if (!is.finite(v)) 1e12 else -v
    }
    ps <- PSCALE[free]; ps[names(ps) %in% names(BOUNDS)] <- 1
    o <- optim(tfull[free], obj, method = "BFGS",
               control = list(parscale = ps, maxit = maxit, reltol = 1e-10))
    t <- tfull; t[free] <- o$par
    par <- from_real(t)
    return(list(par = par, value = -o$value, conv = o$convergence, message = "",
                loglik = sum(unit_loglik_fast(par, d)), evals = o$counts))
  }
  obj <- function(th) {
    par <- full; par[free] <- th
    v <- log_post_fast(par, d, xi_sd)
    if (!is.finite(v)) 1e12 else -v
  }
  st <- pmin(pmax(full[free], LOWER[free] + 1e-6), UPPER[free] - 1e-6)
  o <- nlminb(st, obj, lower = LOWER[free], upper = UPPER[free],
              scale = 1 / PSCALE[free],
              control = list(iter.max = maxit, eval.max = 20 * maxit, rel.tol = 1e-10))
  par <- full; par[free] <- o$par
  list(par = par, value = -o$objective, conv = o$convergence, message = o$message,
       loglik = sum(unit_loglik_fast(par, d)), evals = o$evaluations)
}

## ---- Laplace covariance with fixed finite-difference steps ------------
## Central second differences with absolute steps of one fiftieth of the
## parameter scale, which are large enough to average over the small
## discontinuities of the quadrature rule. Eigenvalues of the curvature that
## are negative or below one thousandth of the largest are floored, so that
## the result is a usable proposal covariance even when the mode is on a
## flat direction.
laplace_cov <- function(par, d, free, step_div = 50) {
  f <- function(p) { v <- log_post_fast(p, d); if (is.finite(v)) -v else NA }
  k <- length(free); h <- PSCALE[free] / step_div
  H <- matrix(NA_real_, k, k)
  f0 <- f(par)
  fp <- fm <- numeric(k)
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
  ## work in scaled coordinates so that the eigenvalue floor is meaningful
  Hs <- diag(PSCALE[free]) %*% ((H + t(H)) / 2) %*% diag(PSCALE[free])
  e <- eigen(Hs, symmetric = TRUE)
  lam <- e$values; floor_val <- 1e-4 * max(lam)
  n_floor <- sum(lam < floor_val)
  lam[lam < floor_val] <- floor_val
  Vs <- e$vectors %*% diag(1 / lam, k) %*% t(e$vectors)
  V <- diag(PSCALE[free]) %*% Vs %*% diag(PSCALE[free])
  dimnames(V) <- list(free, free)
  attr(V, "n_floored") <- n_floor
  V
}
