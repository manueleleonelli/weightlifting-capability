## 52_ext_map.R
## -----------------------------------------------------------------------
## Posterior modes of the base model and of its two extensions, for one sex
## and lift, through the extended kernel (exact declaration Jacobian, units
## grouped by athlete):
##   base  independent declaration errors, exchangeable athlete-days
##   e1    shared self-assessment error per athlete-day (parameter lse)
##   e2    athlete-level Gaussian factor on the copula scale (parameter rhoa)
##   e12   both
## Settings through environment variables:
##   WL_SEX=Men|Women  WL_LIFT=snatch|cleanjerk  WL_MODEL=base|e1|e2|e12
##   WL_G quadrature nodes per piece (default 40 for base and e2, 24 for e1
##        and e12, whose integrand is smoother and whose kernel is slower)
##   WL_KETA nodes per segment of the eta integral (default 8)
##   WL_MAXIT BFGS iterations (default 300); WL_START an .rds fit to start from
##   WL_THREADS OpenMP threads for the kernel (default 2)
##   WL_LAPLACE=1 also computes the Laplace covariance at the mode
## Output: out/extmap_<model>_<sex>_<lift>.rds and a row appended to
## out/ext_map_summary.csv (log likelihood at the fitting grid and at G = 64,
## principal parameters, capability quantiles at the reference bodyweight).
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
sex <- Sys.getenv("WL_SEX", "Men"); lift <- Sys.getenv("WL_LIFT", "snatch")
model <- Sys.getenv("WL_MODEL", "e1")
G <- as.integer(Sys.getenv("WL_G", "24"))
K_ETA <- as.integer(Sys.getenv("WL_KETA", "8"))
MAXIT <- as.integer(Sys.getenv("WL_MAXIT", "300"))
JAC <- as.integer(Sys.getenv("WL_JAC", "0"))
tag <- paste0(model, "_", tolower(sex), "_", lift, if (JAC == 1) "_rel" else "", if (RTAU == 1) "_rtau" else "")
d <- build_data_ext(build_data(att, sex, lift, G = G), K_eta = K_ETA, jac = JAC)
cat(sprintf("== %s | n = %d units, %d athletes | G = %d, K_eta = %d, threads = %d, margin %s\n",
            tag, d$n, length(d$athletes), G, K_ETA, N_THREADS, if (JAC == 1) "relative" else "log"))

## ---- start ----------------------------------------------------------
start_file <- Sys.getenv("WL_START", "")
if (start_file == "") {
  sfx <- paste0(if (JAC == 1) "_rel" else "", if (RTAU == 1) "_rtau" else "")
  cands <- c(if (model == "e12") sprintf("out/extmap_e1_%s_%s%s.rds", tolower(sex), lift, sfx),
             sprintf("out/extmap_base_%s_%s%s.rds", tolower(sex), lift, sfx))
  start_file <- c(cands[file.exists(cands)], "")[1]
}
if (start_file != "" && file.exists(start_file)) {
  prev <- readRDS(start_file); start <- fill_par2(prev$par); cat("starting from", start_file, "\n")
  ## optional: the margin block from another fit (WL_MARGINS_FROM), used for the anchors of a posterior
  ## with two modes, whose margin parameters lie on a flat ridge and can stop at different points
  mf <- Sys.getenv("WL_MARGINS_FROM", "")
  if (mf != "") { mb <- c("m1", "m2", "m3", "lsp1", "lsp2", "lsp3", "lsm", "lsm2", "lsm3"); start[mb] <- fill_par2(readRDS(mf)$par)[mb]; cat("margins from", mf, "\n") }
} else {
  m <- readRDS(sprintf("out/mcmc_%s_%s.rds", tolower(sex), lift))
  start <- fill_par2(apply(m$draws, 2, median))
  cat("starting from the posterior medians of the base chain\n")
}
if (model %in% c("e1", "e12")) { if (!is.finite(start["lse"])) start["lse"] <- log(0.01) } else start["lse"] <- NA
## situation-specific lower scales (WL_SM3=1): start at the shared value
if (Sys.getenv("WL_SM3", "0") == "1") { for (nm in c("lsm2", "lsm3")) if (!is.finite(start[nm])) start[nm] <- start["lsm"]; tag <- paste0(tag, "_sm3") } else start[c("lsm2", "lsm3")] <- NA
## technical-miss floor (WL_FLOOR=1)
if (Sys.getenv("WL_FLOOR", "0") == "1") { if (!is.finite(start["lpi"])) start["lpi"] <- qlogis(0.1); tag <- paste0(tag, "_floor") } else start["lpi"] <- NA
## normal-mixture body with estimated threshold (WL_MIX=1): location, scale, shape and
## mixture parameters from out/mixstart_<sex>_<lift>.rds unless the start file has them
if (Sys.getenv("WL_MIX", "0") == "1") {
  if (!is.finite(start["mx_u"])) {
    ms <- readRDS(sprintf("out/mixstart_%s_%s.rds", tolower(sex), lift))$par
    start[names(ms)] <- ms
    for (k in c("2", "3")) if (start[[paste0("mx_a", k)]] < -6) {   # empty component: move it to the left tail with a small weight
      start[paste0("mx_a", k)] <- -3; start[paste0("mx_m", k)] <- -4; start[paste0("mx_ls", k)] <- log(2.5) }
    if (start[["mx_m2"]] > start[["mx_m3"]]) start[c("mx_m2", "mx_m3", "mx_ls2", "mx_ls3", "mx_a2", "mx_a3")] <- start[c("mx_m3", "mx_m2", "mx_ls3", "mx_ls2", "mx_a3", "mx_a2")]
  }
  start[c("kappa", "lam", "lnu")] <- NA; tag <- paste0(tag, "_mix")
} else start[MIX_NAMES] <- NA
if (model %in% c("e2", "e12")) {
  if (!is.finite(start["rhoa"])) {
    e2f <- sprintf("out/extmap_e2_%s_%s%s.rds", tolower(sex), lift, paste0(if (JAC == 1) "_rel" else "", if (RTAU == 1) "_rtau" else ""))
    start["rhoa"] <- if (file.exists(e2f)) readRDS(e2f)$par["rhoa"] else 0.3
  }
} else start["rhoa"] <- NA
start["rhoa"] <- if (is.finite(start["rhoa"])) min(max(start["rhoa"], 0.02), 0.93) else NA
## relative attempt scale: convert an absolute start value (kg) to a proportion of the typical weight
if (RTAU == 1 && exp(start["log_tau"]) > 0.5) start["log_tau"] <- log(exp(start["log_tau"]) / median(d$W[, 1], na.rm = TRUE))
if (RTAU == 0 && exp(start["log_tau"]) < 0.5) start["log_tau"] <- log(exp(start["log_tau"]) * median(d$W[, 1], na.rm = TRUE))

## fixed parameters, e.g. WL_FIX="xi=-0.3,kappa=0"
FIX <- Sys.getenv("WL_FIX", "")
fixed <- c(); if (FIX != "") { kv <- strsplit(strsplit(FIX, ",")[[1]], "="); fixed <- sapply(kv, function(x) as.numeric(x[2])); names(fixed) <- sapply(kv, function(x) x[1]) }
if (length(fixed)) { start[names(fixed)] <- fixed; tag <- paste0(tag, "_fix_", gsub("[=,]", "_", FIX)) }
tag <- paste0(tag, Sys.getenv("WL_TAGX", ""))   # optional suffix, e.g. for a fit from another start
## an interrupted fit resumes from its checkpoint (the best point it had visited)
ckf <- sprintf("out/extmap_%s_ckpt.rds", tag)
if (file.exists(ckf)) {
  ck <- tryCatch(readRDS(ckf), error = function(e) NULL)
  if (!is.null(ck) && is.finite(ck$value)) {
    start[names(ck$par)] <- ck$par
    cat(sprintf("resuming from the checkpoint: logpost %.2f after %d evaluations\n", ck$value, ck$n_eval))
  }
}

## ---- fit --------------------------------------------------------------
n_eval <- 0; t_start <- Sys.time(); best_v <- -Inf; best_par <- NULL
obj_wrap <- function(par, d) {
  v <- log_post_ext(par, d)
  n_eval <<- n_eval + 1
  if (is.finite(v) && v > best_v) { best_v <<- v; best_par <<- par }
  if (n_eval %% 200 == 0 && !is.null(best_par)) { tmp <- paste0(ckf, ".tmp"); saveRDS(list(par = best_par, value = best_v, n_eval = n_eval), tmp); file.rename(tmp, ckf) }
  if (n_eval %% 50 == 0) cat(sprintf("  eval %5d  logpost %.2f  elapsed %.1f min\n", n_eval, v,
                                     as.numeric(Sys.time() - t_start, units = "mins")))
  v
}
fit_map_ext2 <- function(d, start, maxit) {
  free <- setdiff(PAR_NAMES2[is.finite(start[PAR_NAMES2])], names(fixed))
  full <- start[PAR_NAMES2]
  tfull <- to_real2(full)
  obj <- function(th) { t <- tfull; t[free] <- th; v <- obj_wrap(from_real2(t), d); if (!is.finite(v)) 1e12 else -v }
  ps <- PSCALE2[free]; ps[names(ps) %in% names(BOUNDS2)] <- 1
  o <- optim(tfull[free], obj, method = "BFGS",
             control = list(parscale = ps, maxit = maxit, reltol = as.numeric(Sys.getenv("WL_RELTOL", "1e-10")), ndeps = rep(1e-3, length(free))))
  t <- tfull; t[free] <- o$par
  par <- from_real2(t)
  list(par = par, value = -o$value, conv = o$convergence, free = free, evals = o$counts)
}
cat("start logpost:", log_post_ext(start, d), "\n")
fit <- fit_map_ext2(d, start, MAXIT)
fit$loglik <- sum(unit_loglik_ext(fit$par, d))
fit$time_min <- as.numeric(Sys.time() - t_start, units = "mins")
fit$G <- G; fit$K_eta <- K_ETA; fit$model <- model; fit$sex <- sex; fit$lift <- lift; fit$ref <- d$ref
cat(sprintf("done: logpost %.2f  loglik %.2f  conv %d  evals %s  %.1f min\n", fit$value, fit$loglik, fit$conv,
            paste(fit$evals, collapse = "/"), fit$time_min))
print(round(fit$par, 4))

## log likelihood on a fine grid, comparable across models
d64 <- build_data_ext(build_data(att, sex, lift, G = 64), K_eta = K_ETA, jac = JAC)
fit$loglik64 <- sum(unit_loglik_ext(fit$par, d64))
cat(sprintf("loglik at G = 64: %.2f\n", fit$loglik64))

## Laplace covariance (optional)
if (Sys.getenv("WL_LAPLACE", "0") == "1") {
  cat("Laplace covariance...\n")
  fit$V <- laplace_cov_ext(fit$par, d, fit$free)
  fit$sd <- sqrt(diag(fit$V))
  cat("approximate sd:\n"); print(round(fit$sd, 4))
}
saveRDS(fit, sprintf("out/extmap_%s.rds", tag)); unlink(ckf)

## ---- summary row ----------------------------------------------------
p <- unpack2(fit$par)
bw_ref <- if (sex == "Men") 110 else 86
xd <- design_at(d$ref, bw_ref, 2026)
mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
qs <- mu + om * std_quantile(p, c(0.99, 0.999))
ep <- mu + om * std_endpoint(p)
rec <- if (file.exists("out/post_ceiling.csv")) {
  pc <- read_csv("out/post_ceiling.csv", show_col_types = FALSE)
  pc$record[pc$sex == sex & pc$lift == lift & pc$bodyweight == bw_ref][1]
} else NA
row <- tibble(model = model, margin = if (JAC == 1) "relative" else "log", tau_scale = if (RTAU == 1) "relative" else "absolute", sex = sex, lift = lift, G = G, K_eta = K_ETA, n_par = length(fit$free),
              loglik = fit$loglik, loglik64 = fit$loglik64, logpost = fit$value, conv = fit$conv,
              xi = p$xi, kappa = p$kappa, tau = p$tau, delta2 = p$delta[2], delta3 = p$delta[3],
              s_eta = if (p$s_eta > 0) p$s_eta else NA, rho_a = if (p$rho_a > 0) p$rho_a else NA, pi_miss = if (p$pi_miss > 0) p$pi_miss else NA,
              sL = p$sL, nuL = p$nuL,
              m1 = p$dec[1], m2 = p$dec[2], m3 = p$dec[3], sp1 = p$dec[4], sp2 = p$dec[5], sp3 = p$dec[6], sm = p$dec[7], sm2 = p$dec[8], sm3 = p$dec[9],
              bw_ref = bw_ref, mu_ref = mu, omega_ref = om, q99 = qs[1], q999 = qs[2], endpoint = ep, record = rec,
              xi_sd = if (!is.null(fit$sd)) fit$sd["xi"] else NA,
              lse_sd = if (!is.null(fit$sd) && "lse" %in% names(fit$sd)) fit$sd["lse"] else NA,
              rhoa_sd = if (!is.null(fit$sd) && "rhoa" %in% names(fit$sd)) fit$sd["rhoa"] else NA,
              minutes = fit$time_min, finished = format(Sys.time(), "%Y-%m-%d %H:%M"))
row$fixed <- FIX
row <- row |> mutate(across(everything(), as.character)); row$file <- sprintf("out/extmap_%s.rds", tag)
## one summary file per fit (parallel runs never write to the same file)
dir.create("out/mapsummary", showWarnings = FALSE); write_csv(row, sprintf("out/mapsummary/%s.csv", tag))
f <- "out/ext_map_summary.csv"
if (Sys.getenv("WL_SHARED_SUMMARY", "0") == "1") if (file.exists(f)) { old <- read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (!"tau_scale" %in% names(old)) old$tau_scale <- "absolute"
  if (!"fixed" %in% names(old)) old$fixed <- NA
  old <- old[!(old$model == model & old$margin == row$margin & old$tau_scale == row$tau_scale & old$sex == sex & old$lift == lift & (old$fixed == FIX | (is.na(old$fixed) & FIX == ""))), ]
  write_csv(bind_rows(old, row), f) } else write_csv(row, f)
print(as.data.frame(row[, c("model", "sex", "lift", "loglik", "loglik64", "xi", "tau", "s_eta", "rho_a", "q99", "q999", "endpoint", "record")]))
