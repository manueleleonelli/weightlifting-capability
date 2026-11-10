## 90_shared_error.R
## -----------------------------------------------------------------------
## The final model with a shared self-assessment error added (model e12 of
## 52_ext_map.R: athlete factor and shared error together) against the final
## model (e2), posterior modes at the reference anchor tau = 0.02, for the four
## combinations of sex and lift (Supplement S7).
## The fits compute the shared-error integral with K_eta = 4 nodes per
## segment. As a check, the log likelihood at the e12 mode is recomputed with
## 16 and 64 nodes; a difference of more than a few units flags the fit.
## Where the posterior has a second mode (84_final_select.R) and run_all.R
## fitted the shared error from it as well, a second row compares the two
## modes with and without the shared error (column mode = "second").
## Inputs: out/extmap_e2_<sex>_<lift><SUF>.rds, out/extmap_e12_<sex>_<lift><SUF>.rds
##         (and the same with _second)
## Output: out/shared_error.csv, out/shared_error_rows.tex
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
TEST <- Sys.getenv("RUN_ALL_TEST", "0") == "1"
K_CHECK <- if (TEST) c(4, 8) else c(4, 16, 64)
## summary of a mode: tail shape, order shifts, opener margin, athlete correlation, and the
## one-in-a-thousand level and the endpoint at the reference bodyweight in 2026
summ <- function(m, ref, sex) {
  p <- unpack2(fill_par2(m$par)); bw_ref <- if (sex == "Men") 110 else 86; xd <- design_at(ref, bw_ref, 2026)
  mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
  list(xi = p$xi, delta2 = p$delta[2], delta3 = p$delta[3], m1 = p$dec[1], rho_a = if (p$rho_a > 0) p$rho_a else NA, bw_ref = bw_ref,
       q999 = unname(mu + om * std_quantile(p, 0.999)), endpoint = unname(mu + om * std_endpoint(p)))
}
rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) for (mode in c("kept", "second")) {
  sfx <- if (mode == "second") "_second" else ""
  f2 <- sprintf("out/extmap_e2_%s_%s%s%s.rds", tolower(sex), lift, SUF, sfx); f12 <- sprintf("out/extmap_e12_%s_%s%s%s.rds", tolower(sex), lift, SUF, sfx)
  if (!file.exists(f2) || !file.exists(f12)) { if (mode == "kept") cat("missing fit for", sex, lift, "\n"); next }
  m2 <- readRDS(f2); m12 <- readRDS(f12)
  d4 <- build_data_ext(build_data(att, sex, lift, G = 24), K_eta = K_CHECK[1], jac = 0)
  s2 <- summ(m2, d4$ref, sex); s12 <- summ(m12, d4$ref, sex)
  ## accuracy of the shared-error integral at the e12 mode
  llk <- sapply(K_CHECK, function(K) {
    d <- if (K == K_CHECK[1]) d4 else build_data_ext(build_data(att, sex, lift, G = 24), K_eta = K, jac = 0)
    t0 <- Sys.time(); v <- sum(unit_loglik_ext(fill_par2(m12$par), d))
    cat(sprintf("%s %s (%s mode): e12 log likelihood with K_eta = %d: %.2f (%.0f s)\n", sex, lift, mode, K, v, as.numeric(Sys.time() - t0, units = "secs")))
    v })
  rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, mode = mode,
    loglik_e2 = m2$loglik, loglik_e12 = m12$loglik, gain = m12$loglik - m2$loglik,
    loglik_e12_K16 = if (length(llk) >= 2) llk[2] else NA, loglik_e12_K64 = if (length(llk) >= 3) llk[3] else NA,
    quad_error = max(abs(llk - llk[1])), conv_e12 = m12$conv, minutes_e12 = m12$time_min,
    s_eta = exp(m12$par[["lse"]]), rho_a_e2 = s2$rho_a, rho_a_e12 = s12$rho_a,
    xi_e2 = s2$xi, xi_e12 = s12$xi, delta2_e2 = s2$delta2, delta2_e12 = s12$delta2, delta3_e2 = s2$delta3, delta3_e12 = s12$delta3,
    m1_e2 = s2$m1, m1_e12 = s12$m1, bw_ref = s2$bw_ref,
    q999_e2 = s2$q999, q999_e12 = s12$q999, endpoint_e2 = s2$endpoint, endpoint_e12 = s12$endpoint)
}
out <- bind_rows(rows)
write_csv(out, "out/shared_error.csv")
print(as.data.frame(out |> mutate(across(where(is.numeric), ~ round(.x, 3)))), width = 250)
if (any(out$quad_error > 5, na.rm = TRUE)) cat("WARNING: the shared-error integral is inaccurate at the e12 mode for",
  paste(out$sex[out$quad_error > 5], out$lift[out$quad_error > 5]), "\n")
## LaTeX rows: sex, lift, gain in log likelihood, s_eta, xi, q999 and endpoint without and with the shared error
neg <- function(s) gsub("-", "$-$", s)
LIFT <- c(snatch = "Snatch", cleanjerk = "Clean and jerk")
tex <- sprintf("%s & %s & %s & %.3f & %s & %s & %.1f & %.1f & %.1f & %.1f \\\\", out$sex, paste0(LIFT[out$lift], ifelse(out$mode == "second", ", second mode", "")),
               neg(formatC(round(out$gain), format = "d", big.mark = "{,}")), out$s_eta, neg(sprintf("%.2f", out$xi_e2)), neg(sprintf("%.2f", out$xi_e12)),
               out$q999_e2, out$q999_e12, out$endpoint_e2, out$endpoint_e12)
writeLines(tex, "out/shared_error_rows.tex")
