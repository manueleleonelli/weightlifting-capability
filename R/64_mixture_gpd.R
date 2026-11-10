## 64_mixture_gpd.R
## -----------------------------------------------------------------------
## Mixture bodies with a generalized Pareto tail above a threshold that is
## profiled over, fitted to the best made lift of each athlete-day (days
## treated as exchangeable, as in the body-free check of 63_):
##   MNGPD  y = mu(x) + omega(x) e; e has a three-component normal mixture
##          body (first component N(0, 1), which fixes location and scale)
##          and, above u, density (1 - H(u)) g(e - u; sigma, xi);
##          mu linear and log omega quadratic in the design of the paper
##   MGGPD  y = mu(x) V with log mu linear in the design; V > 0 has a
##          three-component gamma mixture body (first component mean 1) and
##          the generalized Pareto tail above u (Nascimento et al. 2012 with a
##          multiplicative regression)
## For each model the threshold is held at a grid of values (empirical
## quantiles of the residual or ratio from a fit without a tail) and every
## other parameter is maximised; the profile log likelihood over the grid is
## comparable across thresholds because the full-range density is fitted.
## Output: out/mixture_gpd_<sex>.csv
source("R/30_marginal_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
sex <- Sys.getenv("WL_SEX", "Men")
PGRID <- c(0.80, 0.85, 0.90, 0.93, 0.95, 0.97)
softw <- function(a) { e <- exp(c(0, a)); e / sum(e) }
lgpd <- function(r, lsig, xi) { s <- exp(lsig); a <- 1 + xi * r / s; ifelse(a > 0, -lsig - (1 / xi + 1) * log(pmax(a, 1e-300)), -Inf) }

## ---- MNGPD on the residual ------------------------------------------------
nm_A <- c("b0", "b1", "b2", "b3", "g0", "g1", "g2", "m2", "m3", "ls2", "ls3", "a2", "a3", "lsig", "xi")
ll_A <- function(p, y, X, u, tail = TRUE) {
  mu <- as.numeric(X %*% p[1:4]); om <- as.numeric(exp(X[, 1:3] %*% p[5:7])); e <- (y - mu) / om
  m <- c(0, p[["m2"]], p[["m3"]]); s <- c(1, exp(p[["ls2"]]), exp(p[["ls3"]])); w <- softw(c(p[["a2"]], p[["a3"]]))
  if (any(s < 0.05) || any(s > 20) || abs(p[["xi"]]) > 1.5) return(-Inf)
  hb <- w[1] * dnorm(e, m[1], s[1]) + w[2] * dnorm(e, m[2], s[2]) + w[3] * dnorm(e, m[3], s[3])
  lf <- log(pmax(hb, 1e-300))
  if (tail) {
    Hu <- sum(w * pnorm(u, m, s)); hi <- e > u
    lf[hi] <- log(1 - Hu) + lgpd(e[hi] - u, p[["lsig"]], p[["xi"]])
  }
  v <- sum(lf - log(om)); if (is.finite(v)) v else -Inf
}
ep_A <- function(p, xd, u) { if (p[["xi"]] >= 0) return(Inf); sum(p[1:4] * xd) + exp(sum(p[5:7] * xd[1:3])) * (u + exp(p[["lsig"]]) / (-p[["xi"]])) }

## ---- MGGPD on the ratio -----------------------------------------------------
nm_B <- c("b0", "b1", "b2", "b3", "lnu2", "lnu3", "lal1", "lal2", "lal3", "a2", "a3", "lsig", "xi")
ll_B <- function(p, y, X, u, tail = TRUE) {
  mu <- as.numeric(exp(X %*% p[1:4])); v <- y / mu
  nu <- c(1, exp(p[["lnu2"]]), exp(p[["lnu3"]])); al <- exp(c(p[["lal1"]], p[["lal2"]], p[["lal3"]])); w <- softw(c(p[["a2"]], p[["a3"]]))
  if (any(al > 1e5) || any(al < 0.5) || abs(p[["xi"]]) > 1.5) return(-Inf)
  hb <- w[1] * dgamma(v, al[1], al[1] / nu[1]) + w[2] * dgamma(v, al[2], al[2] / nu[2]) + w[3] * dgamma(v, al[3], al[3] / nu[3])
  lf <- log(pmax(hb, 1e-300))
  if (tail) {
    Hu <- sum(w * pgamma(u, al, al / nu)); hi <- v > u
    lf[hi] <- log(1 - Hu) + lgpd(v[hi] - u, p[["lsig"]], p[["xi"]])
  }
  r <- sum(lf - log(mu)); if (is.finite(r)) r else -Inf
}
ep_B <- function(p, xd, u) { if (p[["xi"]] >= 0) return(Inf); exp(sum(p[1:4] * xd)) * (u + exp(p[["lsig"]]) / (-p[["xi"]])) }

fit <- function(ll, start, y, X, u, tail, ps) {
  nms <- names(start)
  obj <- function(t) { names(t) <- nms; v <- ll(t, y, X, u, tail); if (is.finite(v)) -v else 1e12 }
  o <- optim(start, obj, method = "BFGS", control = list(maxit = 3000, reltol = 1e-12, parscale = ps))
  o <- optim(o$par, obj, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12, parscale = ps))
  o <- optim(o$par, obj, method = "BFGS", control = list(maxit = 3000, reltol = 1e-12, parscale = ps))
  names(o$par) <- nms; list(par = o$par, ll = -o$value)
}

rows <- list()
for (lift in c("snatch", "cleanjerk")) {
  d <- build_data(att, sex, lift, G = 8)
  best <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  ok <- !is.na(best); y <- best[ok]; X <- d$X[ok, ]
  bw_ref <- if (sex == "Men") 110 else 86; xd <- design_at(d$ref, bw_ref, 2026)
  ## MNGPD: body-only fit for starting values and the threshold grid
  ls0 <- lm(y ~ X - 1); r0 <- residuals(ls0)
  sA <- c(coef(ls0), log(sd(r0)) , 0, 0, -1, 1, 0.5, 0.3, -2, -1, 0, 0); names(sA) <- nm_A
  sA[5] <- log(sd(r0) / 2); sA["m2"] <- -1.5; sA["m3"] <- 0.5
  psA <- c(1, 1, 1, 0.1, 0.05, 0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 0.2, 0.2, 0.1, 0.05)
  b0 <- fit(ll_A, sA, y, X, Inf, FALSE, psA)
  e0 <- (y - as.numeric(X %*% b0$par[1:4])) / as.numeric(exp(X[, 1:3] %*% b0$par[5:7]))
  cur <- b0$par
  for (pg in PGRID) {
    u <- as.numeric(quantile(e0, pg)); ex <- e0[e0 > u] - u
    cur[["lsig"]] <- log(mean(ex)); cur[["xi"]] <- -0.2
    f <- fit(ll_A, cur, y, X, u, TRUE, psA); cur <- f$par
    rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, model = "MNGPD", p_u = pg, u = u, n_exc = sum(e0 > u),
      loglik = f$ll, xi = f$par[["xi"]], endpoint = ep_A(f$par, xd, u), bw_ref = bw_ref)
    cat(sprintf("%s %s MNGPD p_u %.2f: loglik %.1f xi %.3f endpoint %.1f\n", sex, lift, pg, f$ll, f$par[["xi"]], ep_A(f$par, xd, u)))
  }
  ## MGGPD: body-only fit for starting values and the threshold grid
  lr <- lm(log(y) ~ X - 1)
  sB <- c(coef(lr), -0.3, 0.1, log(200), log(30), log(100), 0, 0, -3, -0.2); names(sB) <- nm_B
  psB <- c(0.01, 0.01, 0.01, 0.001, 0.1, 0.1, 0.3, 0.3, 0.3, 0.2, 0.2, 0.1, 0.05)
  g0 <- fit(ll_B, sB, y, X, Inf, FALSE, psB)
  v0 <- y / as.numeric(exp(X %*% g0$par[1:4]))
  cur <- g0$par
  for (pg in PGRID) {
    u <- as.numeric(quantile(v0, pg)); ex <- v0[v0 > u] - u
    cur[["lsig"]] <- log(mean(ex)); cur[["xi"]] <- -0.2
    f <- fit(ll_B, cur, y, X, u, TRUE, psB); cur <- f$par
    rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, model = "MGGPD", p_u = pg, u = u, n_exc = sum(v0 > u),
      loglik = f$ll, xi = f$par[["xi"]], endpoint = ep_B(f$par, xd, u), bw_ref = bw_ref)
    cat(sprintf("%s %s MGGPD p_u %.2f: loglik %.1f xi %.3f endpoint %.1f\n", sex, lift, pg, f$ll, f$par[["xi"]], ep_B(f$par, xd, u)))
  }
  cat(sprintf("%s %s body only: MN loglik %.1f, MG loglik %.1f\n", sex, lift, b0$ll, g0$ll))
}
write_csv(bind_rows(rows), sprintf("out/mixture_gpd_%s.csv", tolower(sex)))
