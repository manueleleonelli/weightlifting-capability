## 65_mix_check.R : checks of the normal-mixture body in the extended kernel
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
f <- readRDS("out/extmap_e2_men_snatch_rtau_sm3_mix_fix_log_tau_-3.91202.rds")
par <- fill_par2(f$par)
par[c("kappa", "lam", "lnu")] <- NA
p <- unpack2(par)

## 1. density integrates to one; survival matches the integral of the density
ig <- integrate(function(e) exp(log_f_mix_cpp(e, p$bmix)), -40, std_endpoint(p), subdivisions = 2000)$value
cat(sprintf("integral of the density: %.8f\n", ig))
for (e0 in c(-3, 0, 1.4, 1.6, 2.5)) {
  num <- integrate(function(e) exp(log_f_mix_cpp(e, p$bmix)), e0, std_endpoint(p))$value
  cat(sprintf("  survival at %.1f: kernel %.8f  integral %.8f\n", e0, surv_mix_cpp(e0, p$bmix), num))
}
q <- std_quantile(p, c(0.5, 0.9, 0.99, 0.999))
cat("quantiles 0.5 0.9 0.99 0.999:", round(q, 4), " cdf back:", round(1 - surv_mix_cpp(q, p$bmix), 6), "\n")

## 2. unit likelihoods against brute-force integration in R (single-day athletes)
d <- build_data_ext(build_data(att, "Men", "snatch", G = 24), jac = 0)
ll <- unit_loglik_ext(fill_par2(c(par[setdiff(PAR_NAMES2, "rhoa")], rhoa = NA)), d)
one <- which(d$alen == 1)
mu <- as.numeric(d$X %*% p$beta); om <- omega_of(p, d$X)
for (k in pmin(c(3, 300, 1500, 3000), length(one))) {
  ia <- one[k]; i <- d$astart[ia] + 1
  W <- d$W[i, ]; S <- d$Sint[i, ]; D <- d$D[i, ]
  z <- seq(mu[i] + om[i] * min(p$mx$m - 9 * p$mx$s), max(W, na.rm = TRUE) * 1.3, length.out = 400001); dz <- z[2] - z[1]
  e <- (z - mu[i]) / om[i]
  lv <- log_f_mix_cpp(e, p$bmix) - log(om[i])
  for (j in 1:3) if (!is.na(W[j])) {
    eta <- log(z / W[j]) / p$tau + p$delta[j]
    lv <- lv + if (S[j] == 1) plogis(eta, log.p = TRUE) else plogis(-eta, log.p = TRUE)
  }
  okd <- !is.na(W) & D > 0; wd <- W[okd]; sd_ <- D[okd]; wp <- c(0, W[1], W[2])[okd]
  if (length(wd)) lv <- lv + sapply(z, function(zz) log_declaration_cpp(zz, wd, wp, sd_, p$dec, 1L, 0L))
  mx <- max(lv[is.finite(lv)]); ref <- mx + log(sum(exp(lv[is.finite(lv)] - mx)) * dz)
  cat(sprintf("unit %5d: kernel %.5f  brute force %.5f\n", i, ll[ia], ref))
}

## 3. quadrature convergence at G = 24 against G = 96, athlete factor on
par_c <- par; par_c["rhoa"] <- 0.9
d96 <- build_data_ext(build_data(att, "Men", "snatch", G = 96), jac = 0)
t0 <- Sys.time(); a <- unit_loglik_ext(par_c, d); t1 <- Sys.time(); b <- unit_loglik_ext(par_c, d96)
cat(sprintf("athlete factor on: G=24 total %.3f (%.2f s), G=96 total %.3f, max abs diff per athlete %.5f\n",
            sum(a), as.numeric(t1 - t0, units = "secs"), sum(b), max(abs(a - b))))
## 4. threshold derivative is smooth: log likelihood along mx_u
for (uu in c(1.40, 1.45, 1.50, 1.55, 1.60)) { pu <- par_c; pu["mx_u"] <- uu; cat(sprintf("  u = %.2f: %.3f\n", uu, sum(unit_loglik_ext(pu, d)))) }
