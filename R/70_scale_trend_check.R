## 70_scale_trend_check.R
## -----------------------------------------------------------------------
## Does the spread of capability change with time? First rung of the ladder
## (best made lift of each athlete-day, days exchangeable, normal-mixture body,
## GPD tail above an estimated threshold), posterior mode with log omega
## quadratic in bodyweight as in the paper, and with a linear term in the
## centred year added to log omega. Reports the log likelihood, the year
## coefficient of log omega, the tail shape and the endpoint and
## one-in-a-thousand level at the reference bodyweight in 2026 and in 2008,
## against the local records of the current era and of all years.
## Output: out/scale_trend_check.csv
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
rec <- read_csv("out/postx_ceiling.csv", show_col_types = FALSE)
softw <- function(a) { e <- exp(c(0, a)); e / sum(e) }
nll <- function(th, y, X, trend) {
  b <- th[1:4]; g <- if (trend) th[5:8] else c(th[5:7], 0); k <- if (trend) 8 else 7
  m <- c(0, th[k + 1], th[k + 2]); s <- c(1, exp(th[k + 3]), exp(th[k + 4])); w <- softw(th[k + 5:6])
  u <- th[k + 7]; sig <- exp(th[k + 8]); xi <- th[k + 9]
  if (m[2] > m[3] || any(s < 0.05 | s > 20) || u < 0.2 || u > 4 || xi < -1.5 || xi > 0.4) return(1e12)
  mu <- as.numeric(X %*% b); om <- as.numeric(exp(X %*% g)); e <- (y - mu) / om
  lo <- e <= u
  hb <- w[1] * dnorm(e[lo], m[1], s[1]) + w[2] * dnorm(e[lo], m[2], s[2]) + w[3] * dnorm(e[lo], m[3], s[3])
  Su <- sum(w * pnorm(u, m, s, lower.tail = FALSE)); a <- 1 + xi * (e[!lo] - u) / sig
  if (any(a <= 0)) return(1e12)
  ll <- sum(log(pmax(hb, 1e-300))) + sum(log(Su) - log(sig) - (1 / xi + 1) * log(a)) - sum(log(om))
  if (!is.finite(ll)) 1e12 else -ll
}
endpoint <- function(th, xd, trend) {
  k <- if (trend) 8 else 7; g <- if (trend) th[5:8] else c(th[5:7], 0)
  xi <- th[k + 9]; if (xi >= 0) return(Inf)
  sum(th[1:4] * xd) + exp(sum(g * xd)) * (th[k + 7] + exp(th[k + 8]) / (-xi))
}
q999 <- function(th, xd, trend) {
  k <- if (trend) 8 else 7; g <- if (trend) th[5:8] else c(th[5:7], 0)
  m <- c(0, th[k + 1], th[k + 2]); s <- c(1, exp(th[k + 3]), exp(th[k + 4])); w <- softw(th[k + 5:6])
  u <- th[k + 7]; Su <- sum(w * pnorm(u, m, s, lower.tail = FALSE)); xi <- th[k + 9]; sig <- exp(th[k + 8])
  e <- u + sig / xi * (((0.001) / Su)^(-xi) - 1)
  sum(th[1:4] * xd) + exp(sum(g * xd)) * e
}
rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  L1 <- readRDS(sprintf("out/ladder_L1_%s_%s.rds", tolower(sex), lift))
  st <- apply(L1$draws, 2, median)
  d <- build_data(att, sex, lift, G = 8)
  best <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  ok <- !is.na(best); y <- best[ok]; X <- d$X[ok, ]
  th0 <- unname(c(st[c("b0", "b1", "b2", "b3", "g0", "g1", "g2")], st[c("mx_m2", "mx_m3", "mx_ls2", "mx_ls3", "mx_a2", "mx_a3", "mx_u", "mx_lsig", "xi")]))
  ps0 <- c(1, 1, 0.5, 0.05, 0.02, 0.02, 0.02, rep(0.1, 6), 0.05, 0.05, 0.02)
  f0 <- optim(th0, nll, y = y, X = X, trend = FALSE, method = "BFGS", control = list(maxit = 2000, parscale = ps0, reltol = 1e-12))
  f0 <- optim(f0$par, nll, y = y, X = X, trend = FALSE, method = "Nelder-Mead", control = list(maxit = 6000, parscale = ps0, reltol = 1e-12))
  f0 <- optim(f0$par, nll, y = y, X = X, trend = FALSE, method = "BFGS", control = list(maxit = 2000, parscale = ps0, reltol = 1e-12))
  th1 <- c(f0$par[1:7], 0, f0$par[8:16]); ps1 <- c(ps0[1:7], 0.002, ps0[8:16])
  f1 <- optim(th1, nll, y = y, X = X, trend = TRUE, method = "BFGS", control = list(maxit = 2000, parscale = ps1, reltol = 1e-12))
  f1 <- optim(f1$par, nll, y = y, X = X, trend = TRUE, method = "Nelder-Mead", control = list(maxit = 8000, parscale = ps1, reltol = 1e-12))
  f1 <- optim(f1$par, nll, y = y, X = X, trend = TRUE, method = "BFGS", control = list(maxit = 2000, parscale = ps1, reltol = 1e-12))
  bw <- if (sex == "Men") 110 else 86; r <- rec |> filter(sex == !!sex, lift == !!lift, bodyweight == bw)
  for (yr in c(2008, 2026)) {
    xd <- design_at(d$ref, bw, yr)
    rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, year = yr, bodyweight = bw, record = r$record, record_all = r$record_all,
      loglik0 = -f0$value, loglik1 = -f1$value, g3 = f1$par[8], xi0 = f0$par[16], xi1 = f1$par[17],
      endpoint0 = endpoint(f0$par, xd, FALSE), endpoint1 = endpoint(f1$par, xd, TRUE), q999_0 = q999(f0$par, xd, FALSE), q999_1 = q999(f1$par, xd, TRUE))
  }
  print(as.data.frame(bind_rows(rows) |> filter(sex == !!sex, lift == !!lift) |> mutate(across(where(is.numeric), ~round(.x, 3)))))
}
write_csv(bind_rows(rows), "out/scale_trend_check.csv")
