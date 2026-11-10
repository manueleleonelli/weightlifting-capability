## 66_mix_start.R
## -----------------------------------------------------------------------
## Starting values for the latent model with the normal-mixture body: the
## MNGPD of 64_mixture_gpd.R fitted to the best made lifts at the threshold
## quantile with the highest profile likelihood (WL_PU, default from
## out/mixture_gpd_<sex>.csv), saved as out/mixstart_<sex>_<lift>.rds with the
## names of the extended parameter vector.
source("R/30_marginal_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
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


prof <- bind_rows(read_csv("out/mixture_gpd_men.csv", show_col_types = FALSE), read_csv("out/mixture_gpd_women.csv", show_col_types = FALSE))
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  pr <- prof |> filter(sex == !!sex, lift == !!lift, model == "MNGPD")
  pu <- pr$p_u[which.max(pr$loglik)]
  d <- build_data(att, sex, lift, G = 8)
  best <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  ok <- !is.na(best); y <- best[ok]; X <- d$X[ok, ]
  ls0 <- lm(y ~ X - 1); r0 <- residuals(ls0)
  sA <- c(coef(ls0), log(sd(r0)), 0, 0, -1, 1, 0.5, 0.3, -2, -1, 0, 0); names(sA) <- nm_A
  sA[5] <- log(sd(r0) / 2); sA["m2"] <- -1.5; sA["m3"] <- 0.5
  psA <- c(1, 1, 1, 0.1, 0.05, 0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 0.2, 0.2, 0.1, 0.05)
  b0 <- fit(ll_A, sA, y, X, Inf, FALSE, psA)
  e0 <- (y - as.numeric(X %*% b0$par[1:4])) / as.numeric(exp(X[, 1:3] %*% b0$par[5:7]))
  u <- as.numeric(quantile(e0, pu)); cur <- b0$par; cur[["lsig"]] <- log(mean(e0[e0 > u] - u)); cur[["xi"]] <- -0.2
  f <- fit(ll_A, cur, y, X, u, TRUE, psA)
  pp <- f$par
  out <- c(b0 = pp[["b0"]], b1 = pp[["b1"]], b2 = pp[["b2"]], b3 = pp[["b3"]], g0 = pp[["g0"]], g1 = pp[["g1"]], g2 = pp[["g2"]],
           xi = pp[["xi"]], mx_m2 = pp[["m2"]], mx_m3 = pp[["m3"]], mx_ls2 = pp[["ls2"]], mx_ls3 = pp[["ls3"]],
           mx_a2 = pp[["a2"]], mx_a3 = pp[["a3"]], mx_u = u, mx_lsig = pp[["lsig"]])
  if (out[["mx_m2"]] > out[["mx_m3"]]) out[c("mx_m2", "mx_m3", "mx_ls2", "mx_ls3", "mx_a2", "mx_a3")] <- out[c("mx_m3", "mx_m2", "mx_ls3", "mx_ls2", "mx_a3", "mx_a2")]
  saveRDS(list(par = out, loglik = f$ll, p_u = pu), sprintf("out/mixstart_%s_%s.rds", tolower(sex), lift))
  cat(sprintf("%s %s: p_u %.2f u %.3f loglik %.1f xi %.3f\n", sex, lift, pu, u, f$ll, pp[["xi"]])); print(round(out, 4))
}
