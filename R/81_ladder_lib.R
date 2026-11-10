## 81_ladder_lib.R
## -----------------------------------------------------------------------
## Log posterior of the first two rungs of the ladder (Table 5 of the paper):
## the best made lift y of each athlete-day taken as the capability, with the
## capability distribution of the final model (location and log scale in
## bodyweight and year, three-component normal-mixture body, generalized
## Pareto tail above the threshold u) and the same priors.
##   L1  days exchangeable
##   L2  the days of one athlete share a Gaussian factor on the copula scale
##       with equicorrelation rho_a (closed form for each athlete)
## The mixture components 2 and 3 are exchangeable: the density does not
## depend on their labels, and the prior is symmetric in them.
LADDER_NAMES <- c("b0", "b1", "b2", "b3", "g0", "g1", "g2", "xi", MIX_NAMES)
LADDER_NEUTRAL <- c(delta2 = 0, delta3 = 0, log_tau = log(0.02), m1 = 0.03, m2 = 0.03, m3 = 0.03,
                    lsp1 = log(0.03), lsp2 = log(0.03), lsp3 = log(0.03), lsm = log(0.03))
ladder_data <- function(att, sex, lift) {
  d <- build_data(att, sex, lift, G = 8)
  best <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  ok <- !is.na(best)
  list(y = best[ok], X = d$X[ok, ], aid = as.integer(factor(d$u$athlete[ok])), ref = d$ref, n = sum(ok))
}
.lmix <- function(e, w, m, s) {
  l1 <- log(w[1]) + dnorm(e, m[1], s[1], log = TRUE); l2 <- log(w[2]) + dnorm(e, m[2], s[2], log = TRUE)
  l3 <- log(w[3]) + dnorm(e, m[3], s[3], log = TRUE); mx <- pmax(l1, l2, l3)
  mx + log(exp(l1 - mx) + exp(l2 - mx) + exp(l3 - mx))
}
ladder_logpost <- function(p, L, athlete) {
  pp <- fill_par2(c(p[intersect(names(p), c(LADDER_NAMES, "rhoa"))], LADDER_NEUTRAL))
  pp[c("kappa", "lam", "lnu")] <- c(0, 0, log(10))
  lp <- log_prior_ext(pp[setdiff(names(pp), "rhoa")], 0.3); if (!is.finite(lp)) return(-Inf)
  q <- unpack2(pp); mx <- q$mx; xi <- q$xi
  mu <- as.numeric(L$X %*% q$beta); om <- as.numeric(exp(L$X[, 1:3] %*% q$gamma)); e <- (L$y - mu) / om
  Su <- sum(mx$w * pnorm(mx$u, mx$m, mx$s, lower.tail = FALSE)); lo <- e <= mx$u
  a <- 1 + xi * (e[!lo] - mx$u) / mx$sig; if (any(a <= 0)) return(-Inf)
  lf <- numeric(length(e)); lf[lo] <- .lmix(e[lo], mx$w, mx$m, mx$s); lf[!lo] <- log(Su) - log(mx$sig) - (1 / xi + 1) * log(a)
  ll <- sum(lf - log(om))
  if (athlete) {
    rho <- p[["rhoa"]]; if (!(rho > 0 && rho < 0.95)) return(-Inf)
    S <- numeric(length(e)); el <- e[lo]
    S[lo] <- mx$w[1] * pnorm(el, mx$m[1], mx$s[1], lower.tail = FALSE) + mx$w[2] * pnorm(el, mx$m[2], mx$s[2], lower.tail = FALSE) +
      mx$w[3] * pnorm(el, mx$m[3], mx$s[3], lower.tail = FALSE)
    S[!lo] <- Su * a^(-1 / xi)
    cz <- qnorm(pmin(pmax(1 - S, 1e-12), 1 - 1e-12))
    n <- as.numeric(rowsum(rep(1, length(cz)), L$aid)); s1 <- as.numeric(rowsum(cz, L$aid)); s2 <- as.numeric(rowsum(cz^2, L$aid))
    ll <- ll + sum(-0.5 * ((n - 1) * log(1 - rho) + log(1 + (n - 1) * rho)) - 0.5 * (s2 - rho / (1 + (n - 1) * rho) * s1^2) / (1 - rho) + 0.5 * s2)
  }
  ll + lp
}

## a start with an empty or degenerate component (weight below e^-4, or centred
## above the threshold, where the body is not used) gets that component moved to
## the left tail with a small weight, as in 52_ext_map.R
sanitize_mix <- function(par) {
  for (k in c("2", "3")) {
    a <- paste0("mx_a", k); m <- paste0("mx_m", k); s <- paste0("mx_ls", k)
    if (par[[a]] < -4 || par[[m]] > par[["mx_u"]] || par[[s]] < log(0.1)) { par[[a]] <- -3; par[[m]] <- -4; par[[s]] <- log(2.5) }
  }
  par
}

## The best made lift is a whole number of kilograms, so its likelihood is the
## probability of the interval [y - 0.5, y + 0.5) rather than a density. This
## is the correct likelihood for the recorded data, and unlike the density it
## is a smooth function of the threshold u: with a density, the log likelihood
## jumps every time u crosses an observation, and a random-walk sampler cannot
## move the threshold across thousands of such steps.
LADDER_IC <- Sys.getenv("WL_LADDER_IC", "1") == "1"
.surv_mix <- function(e, mx, xi) {   # survival function of the standardized error
  Su <- sum(mx$w * pnorm(mx$u, mx$m, mx$s, lower.tail = FALSE)); out <- numeric(length(e)); lo <- e <= mx$u
  if (any(lo)) out[lo] <- mx$w[1] * pnorm(e[lo], mx$m[1], mx$s[1], lower.tail = FALSE) + mx$w[2] * pnorm(e[lo], mx$m[2], mx$s[2], lower.tail = FALSE) +
    mx$w[3] * pnorm(e[lo], mx$m[3], mx$s[3], lower.tail = FALSE)
  if (any(!lo)) { a <- pmax(1 + xi * (e[!lo] - mx$u) / mx$sig, 0); out[!lo] <- Su * a^(-1 / xi) }
  out
}
.cdf_mix_body <- function(e, mx) mx$w[1] * pnorm(e, mx$m[1], mx$s[1]) + mx$w[2] * pnorm(e, mx$m[2], mx$s[2]) + mx$w[3] * pnorm(e, mx$m[3], mx$s[3])
ladder_logpost_ic <- function(p, L, athlete) {
  pp <- fill_par2(c(p[intersect(names(p), c(LADDER_NAMES, "rhoa"))], LADDER_NEUTRAL))
  pp[c("kappa", "lam", "lnu")] <- c(0, 0, log(10))
  lp <- log_prior_ext(pp[setdiff(names(pp), "rhoa")], 0.3); if (!is.finite(lp)) return(-Inf)
  q <- unpack2(pp); mx <- q$mx; xi <- q$xi
  mu <- as.numeric(L$X %*% q$beta); om <- as.numeric(exp(L$X[, 1:3] %*% q$gamma))
  elo <- (L$y - 0.5 - mu) / om; ehi <- (L$y + 0.5 - mu) / om
  if (xi < 0 && any(elo >= mx$u + mx$sig / (-xi))) return(-Inf)       # an observation beyond the endpoint
  pr <- numeric(length(elo)); low <- ehi <= 0                          # lower body: differences of the cdf; elsewhere of the survival
  if (any(low)) pr[low] <- .cdf_mix_body(ehi[low], mx) - .cdf_mix_body(elo[low], mx)
  if (any(!low)) pr[!low] <- .surv_mix(elo[!low], mx, xi) - .surv_mix(ehi[!low], mx, xi)
  if (any(!(pr > 0))) return(-Inf)
  ll <- sum(log(pr))
  if (athlete) {
    rho <- p[["rhoa"]]; if (!(rho > 0 && rho < 0.95)) return(-Inf)
    e <- (L$y - mu) / om; S <- .surv_mix(e, mx, xi)
    cz <- qnorm(pmin(pmax(1 - S, 1e-12), 1 - 1e-12))
    n <- as.numeric(rowsum(rep(1, length(cz)), L$aid)); s1 <- as.numeric(rowsum(cz, L$aid)); s2 <- as.numeric(rowsum(cz^2, L$aid))
    ll <- ll + sum(-0.5 * ((n - 1) * log(1 - rho) + log(1 + (n - 1) * rho)) - 0.5 * (s2 - rho / (1 + (n - 1) * rho) * s1^2) / (1 - rho) + 0.5 * s2)
  }
  ll + lp
}
if (LADDER_IC) ladder_logpost <- ladder_logpost_ic
