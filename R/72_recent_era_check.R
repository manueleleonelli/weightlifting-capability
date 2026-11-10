## 72_recent_era_check.R
## -----------------------------------------------------------------------
## Does the tail come from the early years? The record-book fit of the first
## rung of the ladder (best made lift, days exchangeable, normal-mixture body,
## GPD tail above an estimated threshold), posterior mode, on all athlete-days
## and on those from 2013 and from 2017 onwards, with the endpoint and the
## one-in-a-thousand level at the reference bodyweight in 2026.
## Output: out/recent_era_check.csv
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
softw <- function(a) { e <- exp(c(0, a)); e / sum(e) }
nll <- function(th, y, X) {
  b <- th[1:4]; g <- th[5:7]
  m <- c(0, th[8], th[9]); s <- c(1, exp(th[10]), exp(th[11])); w <- softw(th[12:13])
  u <- th[14]; sig <- exp(th[15]); xi <- th[16]
  if (m[2] > m[3] || any(s < 0.05 | s > 20) || u < 0.2 || u > 4 || xi < -1.5 || xi > 0.4) return(1e12)
  mu <- as.numeric(X %*% b); om <- as.numeric(exp(X[, 1:3] %*% g)); e <- (y - mu) / om
  lo <- e <= u
  hb <- w[1] * dnorm(e[lo], m[1], s[1]) + w[2] * dnorm(e[lo], m[2], s[2]) + w[3] * dnorm(e[lo], m[3], s[3])
  Su <- sum(w * pnorm(u, m, s, lower.tail = FALSE)); a <- 1 + xi * (e[!lo] - u) / sig
  if (any(a <= 0)) return(1e12)
  ll <- sum(log(pmax(hb, 1e-300))) + sum(log(Su) - log(sig) - (1 / xi + 1) * log(a)) - sum(log(om))
  if (!is.finite(ll)) 1e12 else -ll
}
summ <- function(th, xd) {
  m <- c(0, th[8], th[9]); s <- c(1, exp(th[10]), exp(th[11])); w <- softw(th[12:13])
  u <- th[14]; sig <- exp(th[15]); xi <- th[16]; Su <- sum(w * pnorm(u, m, s, lower.tail = FALSE))
  mu <- sum(th[1:4] * xd); om <- exp(sum(th[5:7] * xd[1:3]))
  c(xi = xi, q999 = mu + om * (u + sig / xi * ((0.001 / Su)^(-xi) - 1)), endpoint = if (xi < 0) mu + om * (u + sig / (-xi)) else Inf)
}
rec <- read_csv("out/postx_ceiling.csv", show_col_types = FALSE)
rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  L1 <- readRDS(sprintf("out/ladder_L1_%s_%s.rds", tolower(sex), lift)); st <- apply(L1$draws, 2, median)
  th0 <- unname(c(st[c("b0", "b1", "b2", "b3", "g0", "g1", "g2")], st[c("mx_m2", "mx_m3", "mx_ls2", "mx_ls3", "mx_a2", "mx_a3", "mx_u", "mx_lsig", "xi")]))
  ps <- c(1, 1, 0.5, 0.05, 0.02, 0.02, 0.02, rep(0.1, 6), 0.05, 0.05, 0.02)
  d <- build_data(att, sex, lift, G = 8)
  best <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  bw <- if (sex == "Men") 110 else 86; xd <- design_at(d$ref, bw, 2026)
  r <- rec |> filter(sex == !!sex, lift == !!lift, bodyweight == bw)
  for (from in c(1998, 2013, 2017)) {
    ok <- !is.na(best) & d$u$year >= from; y <- best[ok]; X <- d$X[ok, ]
    f <- optim(th0, nll, y = y, X = X, method = "BFGS", control = list(maxit = 2000, parscale = ps, reltol = 1e-12))
    f <- optim(f$par, nll, y = y, X = X, method = "Nelder-Mead", control = list(maxit = 6000, parscale = ps, reltol = 1e-12))
    f <- optim(f$par, nll, y = y, X = X, method = "BFGS", control = list(maxit = 2000, parscale = ps, reltol = 1e-12))
    sm <- summ(f$par, xd)
    rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, from = from, n = length(y), bodyweight = bw, record = r$record,
                                       xi = sm[["xi"]], q999 = sm[["q999"]], endpoint = sm[["endpoint"]], trend = f$par[4])
    print(as.data.frame(tail(rows, 1)[[1]] |> mutate(across(where(is.numeric), ~round(.x, 3)))))
  }
}
write_csv(bind_rows(rows), "out/recent_era_check.csv")
