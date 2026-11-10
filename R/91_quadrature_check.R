## 91_quadrature_check.R
## -----------------------------------------------------------------------
## Accuracy of the quadrature rules (Table S1 of the supplement): each rule
## against a more accurate reference, at the posterior modes of the final
## model at tau = 0.02 (and of the shared-error variant for the last two).
##   capability      the rule of Section S1.1 against a grid of 400,001
##                   points, for WL_NUNIT athletes with a single competition
##                   per sex and lift, athlete factor off
##   athlete factor  12 Gauss-Hermite nodes against 20 nodes and against a
##                   grid of 4,001 points, every athlete
##   shared error    K_eta = 4 nodes per segment against a grid of 200,001
##                   points, for WL_NUNIT athlete-day-lifts with at least two
##                   declared weights per sex and lift, at three capabilities
##                   around the heaviest declared weight; and K_eta = 4
##                   against 64 in the total log likelihood (from
##                   out/shared_error.csv of 90_shared_error.R)
## Output: out/quadrature_check.csv (largest absolute difference in log
## likelihood for each check)
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
N_UNIT <- as.integer(Sys.getenv("WL_NUNIT", "25"))
res <- list(cap = 0, gh20 = 0, grid = 0, eta = 0)
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  cat("==", sex, lift, "\n")
  par <- fill_par2(readRDS(sprintf("out/extmap_e2_%s_%s%s.rds", tolower(sex), lift, SUF))$par)
  p <- unpack2(par)
  d <- build_data_ext(build_data(att, sex, lift, G = 24), jac = 0)
  mu <- as.numeric(d$X %*% p$beta); om <- omega_of(p, d$X)

  ## 1. capability, athlete factor off
  ll <- unit_loglik_ext(fill_par2(c(par[setdiff(PAR_NAMES2, "rhoa")], rhoa = NA)), d)
  one <- which(d$alen == 1)
  for (ia in one[unique(round(seq(1, length(one), length.out = N_UNIT)))]) {
    i <- d$astart[ia] + 1; W <- d$W[i, ]; S <- d$Sint[i, ]; D <- d$D[i, ]
    z <- seq(mu[i] + om[i] * min(p$mx$m - 9 * p$mx$s), max(W, na.rm = TRUE) * 1.3, length.out = 400001); dz <- z[2] - z[1]
    z <- z[z > 0]
    lv <- log_f_mix_cpp((z - mu[i]) / om[i], p$bmix) - log(om[i])
    for (j in 1:3) if (!is.na(W[j])) {
      eta <- log(z / W[j]) / p$tau + p$delta[j]
      lv <- lv + if (S[j] == 1) plogis(eta, log.p = TRUE) else plogis(-eta, log.p = TRUE)
    }
    okd <- !is.na(W) & D > 0
    if (any(okd)) lv <- lv + sapply(z, function(zz) log_declaration_cpp(zz, W[okd], c(0, W[1], W[2])[okd], D[okd], p$dec, 1L, 0L))
    ok <- is.finite(lv); mx <- max(lv[ok])
    res$cap <- max(res$cap, abs(ll[ia] - (mx + log(sum(exp(lv[ok] - mx)) * dz))))
  }
  cat(sprintf("  capability: %.2e\n", res$cap))

  ## 2. athlete factor
  a12 <- unit_loglik_ext(par, d)
  d20 <- d; d20$gh <- gauss_hermite(20); a20 <- unit_loglik_ext(par, d20)
  dg <- d; x <- seq(-8, 8, length.out = 4001); dg$gh <- list(x = x, w = (x[2] - x[1]) * exp(-x^2)); ag <- unit_loglik_ext(par, dg)
  res$gh20 <- max(res$gh20, abs(a12 - a20)); res$grid <- max(res$grid, abs(a12 - ag))
  cat(sprintf("  athlete factor: 20 nodes %.2e, grid %.2e\n", max(abs(a12 - a20)), max(abs(a12 - ag))))

  ## 3. shared error
  f12 <- sprintf("out/extmap_e12_%s_%s%s.rds", tolower(sex), lift, SUF)
  if (file.exists(f12)) {
    p12 <- unpack2(fill_par2(readRDS(f12)$par)); gl <- gauss_legendre(4)
    two <- which(rowSums(!is.na(d$W) & d$D > 0) >= 2)
    for (i in two[unique(round(seq(1, length(two), length.out = N_UNIT)))]) {
      W <- d$W[i, ]; D <- d$D[i, ]; okd <- !is.na(W) & D > 0
      wd <- W[okd]; wp <- c(0, W[1], W[2])[okd]; sd_ <- D[okd]
      for (z in max(wd) * c(1.00, 1.03, 1.06)) {
        k <- log_eta_integral_pt_cpp(z, wd, wp, sd_, p12$dec, p12$s_eta, gl$x, gl$w, 0L)
        r <- ref_eta_integral(z, wd, sd_, p12$dec, p12$s_eta, wp, trunc = TRUE, jac = 0, n_grid = 200001)
        if (is.finite(k) && is.finite(r)) res$eta <- max(res$eta, abs(k - r))
      }
    }
    cat(sprintf("  shared error: %.2e\n", res$eta))
  }
}
se <- if (file.exists("out/shared_error.csv")) max(read_csv("out/shared_error.csv", show_col_types = FALSE)$quad_error) else NA
out <- tibble(integral = c("capability", "athlete factor", "athlete factor", "shared error", "shared error"),
              rule = c("rule of S1.1", "12 Gauss-Hermite nodes", "12 Gauss-Hermite nodes", "K_eta = 4", "K_eta = 4"),
              reference = c("grid of 400,001 points", "20 nodes", "grid of 4,001 points", "grid of 200,001 points", "K_eta = 64"),
              largest_difference = c(res$cap, res$gh20, res$grid, res$eta, se),
              per = c("athlete-day-lift", "athlete", "athlete", "athlete-day-lift", "total"))
write_csv(out, "out/quadrature_check.csv"); print(as.data.frame(out))
