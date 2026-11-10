## 57_ext_copula.R
## -----------------------------------------------------------------------
## Second stage through the extended marginal model: dependence between the
## two lifts. For each sex, the joint marginal likelihood of the paired
## athlete-days under a t copula (rho, psi) is evaluated on a grid with the
## marginal parameters at posterior draws (or at the mode) of the first
## stage, the posterior over the grid averaged over those draws. Athlete-days
## are treated as exchangeable in this stage (the athlete factor of the
## marginal model is switched off). Also writes latent pair draws, chi(u) on
## the latent capabilities and the return level of the total.
## Inputs: WL_FS_<SEX>, WL_FC_<SEX> (paths of the snatch and clean-and-jerk
## first-stage files, mcmcx_*.rds with $draws or extmap_*.rds with $par);
## WL_NMARG marginal draws (default 3), WL_NLAT latent draws (default 40).
## Outputs: out/copulax_<sex>.rds, out/copulax_summary.csv,
##          out/latentx_draws_<sex>.rds, out/chix_latent.csv, out/totalx_returnlevel.csv
source("R/50_ext_lib.R")
set.seed(7)
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
N_MARG <- as.integer(Sys.getenv("WL_NMARG", "3")); N_LAT <- as.integer(Sys.getenv("WL_NLAT", "40"))
rho_grid <- c(0.88, 0.90, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98)
psi_grid <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
get_draws <- function(f) {
  x <- readRDS(f)
  D <- if (!is.null(x$draws)) x$draws else matrix(x$par, 1, dimnames = list(NULL, names(x$par)))
  M <- matrix(NA_real_, nrow(D), length(PAR_NAMES2), dimnames = list(NULL, PAR_NAMES2))
  for (nm in intersect(colnames(D), PAR_NAMES2)) M[, nm] <- D[, nm]
  D <- M; D[, "rhoa"] <- NA   # exchangeable days in this stage
  list(D = D, ref = x$ref)
}
summ_rows <- list(); chi_rows <- list(); ret_rows <- list()
for (sex in c("Men", "Women")) {
  fs <- Sys.getenv(paste0("WL_FS_", toupper(sex)), ""); fc <- Sys.getenv(paste0("WL_FC_", toupper(sex)), "")
  if (fs == "" || !file.exists(fs) || !file.exists(fc)) { cat("no first-stage files for", sex, "\n"); next }
  ms <- get_draws(fs); mc <- get_draws(fc)
  ds <- build_data(att, sex, "snatch", G = 24); dc <- build_data(att, sex, "cleanjerk", G = 24)
  pr <- pair_data(ds, dc)
  cat(sprintf("== %s: %d paired athlete-days\n", sex, pr$s$n))
  nm <- min(N_MARG, nrow(ms$D)); is <- round(seq(1, nrow(ms$D), length.out = nm)); ic <- round(seq(1, nrow(mc$D), length.out = nm))
  post <- matrix(0, length(rho_grid), length(psi_grid)); t0 <- Sys.time()
  for (k in seq_len(nm)) {
    par_s <- ms$D[is[k], ]; par_c <- mc$D[ic[k], ]
    ll <- matrix(NA_real_, length(rho_grid), length(psi_grid))
    for (i in seq_along(rho_grid)) for (j in seq_along(psi_grid))
      ll[i, j] <- sum(joint_loglik_ext(par_s, par_c, pr$s, pr$c, rho_grid[i], psi_grid[j])$ll)
    w <- exp(ll - max(ll)); post <- post + w / sum(w) / nm
    cat(sprintf("  marginal draw %d: best rho %.3f psi %.3f (%s)\n", k, rho_grid[which.max(apply(ll, 1, max))],
                psi_grid[which.max(apply(ll, 2, max))], format(Sys.time() - t0)))
  }
  p_rho <- rowSums(post); p_psi <- colSums(post)
  cdf <- function(p, g) function(q) g[which(cumsum(p) >= q)[1]]
  summ_rows[[sex]] <- tibble(
    sex = sex, n_pairs = pr$s$n,
    rho_median = cdf(p_rho, rho_grid)(0.5), rho_lo = cdf(p_rho, rho_grid)(0.025), rho_hi = cdf(p_rho, rho_grid)(0.975),
    psi_median = cdf(p_psi, psi_grid)(0.5), psi_lo = cdf(p_psi, psi_grid)(0.025), psi_hi = cdf(p_psi, psi_grid)(0.975),
    p_psi_zero = p_psi[1], p_psi_le_002 = sum(p_psi[psi_grid <= 0.02]),
    eta_gauss_median = (1 + cdf(p_rho, rho_grid)(0.5)) / 2)
  print(as.data.frame(summ_rows[[sex]]))
  rho_hat <- summ_rows[[sex]]$rho_median; psi_hat <- summ_rows[[sex]]$psi_median
  saveRDS(list(sex = sex, rho_grid = rho_grid, psi_grid = psi_grid, post = post, rho = rho_hat, psi = psi_hat),
          sprintf("out/copulax_%s.rds", tolower(sex)))

  Zarr <- array(NA_real_, c(pr$s$n, 2, N_LAT)); ks <- sample(nrow(ms$D), N_LAT, replace = TRUE); kc <- sample(nrow(mc$D), N_LAT, replace = TRUE)
  for (k in seq_len(N_LAT)) {
    jd <- joint_loglik_ext(ms$D[ks[k], ], mc$D[kc[k], ], pr$s, pr$c, rho_hat, psi_hat, draw = TRUE)
    Zarr[, 1, k] <- jd$zs; Zarr[, 2, k] <- jd$zc
  }
  saveRDS(list(sex = sex, Z = Zarr, units = pr$s$u$unit, athlete = pr$s$u$athlete,
               bodyweight = pr$s$u$bodyweight, year = pr$s$u$year, rho = rho_hat, psi = psi_hat),
          sprintf("out/latentx_draws_%s.rds", tolower(sex)))
  for (k in seq_len(N_LAT)) {
    p_s <- unpack2(ms$D[ks[k], ]); p_c <- unpack2(mc$D[kc[k], ])
    es <- (Zarr[, 1, k] - as.numeric(pr$s$X %*% p_s$beta)) / omega_of(p_s, pr$s$X)
    ec <- (Zarr[, 2, k] - as.numeric(pr$c$X %*% p_c$beta)) / omega_of(p_c, pr$c$X)
    us <- rank(es) / (length(es) + 1); uc <- rank(ec) / (length(ec) + 1)
    for (u in c(0.5, 0.75, 0.9, 0.95, 0.975, 0.99))
      chi_rows[[length(chi_rows) + 1]] <- tibble(sex = sex, draw = k, u = u, chi = mean(uc[us > u] > u))
  }
  nu <- if (psi_hat > 0) 1 / psi_hat else Inf
  M <- 400000
  z1 <- rnorm(M); z2 <- rho_hat * z1 + sqrt(1 - rho_hat^2) * rnorm(M)
  if (is.finite(nu)) { wsc <- sqrt(nu / rchisq(M, nu)); u1 <- pt(z1 * wsc, nu); u2 <- pt(z2 * wsc, nu) } else { u1 <- pnorm(z1); u2 <- pnorm(z2) }
  for (u in c(0.5, 0.75, 0.9, 0.95, 0.975, 0.99))
    chi_rows[[length(chi_rows) + 1]] <- tibble(sex = sex, draw = 0, u = u, chi = mean(u2[u1 > u] > u))

  bw <- if (sex == "Men") 110 else 86
  for (N_RET in c(1000, 10000)) {
    joint <- numeric(40); marg <- numeric(40)
    for (k in 1:40) {
      p_s <- unpack2(ms$D[sample(nrow(ms$D), 1), ]); p_c <- unpack2(mc$D[sample(nrow(mc$D), 1), ])
      xd <- design_at(ms$ref, bw, 2026); xdc <- design_at(mc$ref, bw, 2026)
      mus <- sum(p_s$beta * xd); oms <- omega_of(p_s, matrix(xd, 1)); muc <- sum(p_c$beta * xdc); omc <- omega_of(p_c, matrix(xdc, 1))
      Mm <- 200000
      z1 <- rnorm(Mm); z2 <- rho_hat * z1 + sqrt(1 - rho_hat^2) * rnorm(Mm)
      if (is.finite(nu)) { wsc <- sqrt(nu / rchisq(Mm, nu)); u1 <- pt(z1 * wsc, nu); u2 <- pt(z2 * wsc, nu) } else { u1 <- pnorm(z1); u2 <- pnorm(z2) }
      zs <- mus + oms * std_quantile(p_s, u1); zc <- muc + omc * std_quantile(p_c, u2)
      pr_ret <- 1 - 1 / N_RET
      joint[k] <- quantile(zs + zc, pr_ret)
      marg[k] <- mus + oms * std_quantile(p_s, pr_ret) + muc + omc * std_quantile(p_c, pr_ret)
    }
    ret_rows[[length(ret_rows) + 1]] <- tibble(sex = sex, bodyweight = bw, horizon = N_RET,
      total_joint = median(joint), joint_lo = quantile(joint, .025), joint_hi = quantile(joint, .975),
      sum_marginal = median(marg), deficit = median(marg - joint),
      deficit_lo = quantile(marg - joint, .025), deficit_hi = quantile(marg - joint, .975))
  }
}
write_csv(bind_rows(summ_rows), "out/copulax_summary.csv")
chi <- bind_rows(chi_rows) |> group_by(sex, u, model = draw == 0) |>
  summarise(lo = quantile(chi, .025), hi = quantile(chi, .975), chi = median(chi), .groups = "drop")
write_csv(chi, "out/chix_latent.csv"); print(as.data.frame(chi))
write_csv(bind_rows(ret_rows), "out/totalx_returnlevel.csv"); print(as.data.frame(bind_rows(ret_rows)))
