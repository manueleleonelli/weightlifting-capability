## 69_anchor_mix.R
## -----------------------------------------------------------------------
## The final model (normal-mixture body, GPD tail above an estimated
## threshold, athlete factor) at the three anchors tau = 0.01, 0.02, 0.03:
## posterior modes from 52_ext_map.R (out/extmap_e2_<sex>_<lift>_rtau_sm3_mix_
## fix_log_tau_<value>.rds), and the fit with tau free when it exists. For each fit: log likelihood and log posterior,
## attempt-order shifts, opener and raise modes, location, tail shape,
## threshold, endpoint and one-in-a-thousand level at the reference
## bodyweight in 2026, and the mean of capability minus the best made lift
## (all athlete-days with a made lift, bracketed days, days that ended on a
## success), the capabilities evaluated at their conditional means given the
## day with the athlete factor off.
## Where the posterior has a second mode and run_all.R fitted the anchors from
## it as well (files with the suffix _second), the same rows are added with
## mode = "second". Each of those anchors is also fitted from the margins of the
## other mode (suffix _swap), and the fit with the higher log posterior is used.
## Output: out/anchor_mix.csv (out/anchor_mix_free.csv for the fit with tau free)
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
TAUS <- c("0.01" = "-4.60517", "0.02" = "-3.91202", "0.03" = "-3.50656")
rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  dd <- NULL
  for (mode in c("kept", "second")) for (k in c(names(TAUS), "free")) {
    if (mode == "second" && k == "free") next
    sfx <- if (mode == "second") "_second" else ""
    f <- if (k == "free") sprintf("out/extmap_e2_%s_%s_rtau_sm3_mix.rds", tolower(sex), lift) else
      sprintf("out/extmap_e2_%s_%s_rtau_sm3_mix_fix_log_tau_%s%s.rds", tolower(sex), lift, TAUS[[k]], sfx)
    if (!file.exists(f)) { if (mode == "kept") cat("missing", f, "\n"); next }
    fit <- readRDS(f)
    ## the same anchor started from the margins of the other mode (run_all.R, stage 4k): the higher of the two is kept
    fs <- sub("\\.rds$", "_swap.rds", f)
    if (file.exists(fs)) { alt <- readRDS(fs); if (alt$value > fit$value) { cat(sprintf("%s %s %s tau %s: the start from the other mode's margins is higher (%.2f against %.2f)\n", sex, lift, mode, k, alt$value, fit$value)); fit <- alt } }
    if (is.null(dd)) {
      dd <- build_data_ext(build_data(att, sex, lift, G = 24), jac = 0)
      best <- apply(ifelse(dd$S == 1, dd$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
      ended <- !is.na(best) & apply(dd$Sint, 1, function(s) { s <- s[s >= 0]; s[length(s)] == 1 })
      brk <- !is.na(best) & dd$has_miss & dd$has_make
    }
    par <- fill_par2(fit$par); p <- unpack2(par)
    if (is.null(fit$loglik)) fit$loglik <- sum(unit_loglik_ext(par, dd))   # a fit stopped before the end (tau free)
    bw <- if (sex == "Men") 110 else 86; xd <- design_at(dd$ref, bw, 2026)
    mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
    L <- latent_summary_ext(par, dd); g <- L[, 1] - best
    rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, mode = mode, free = k == "free", tau = exp(par[["log_tau"]]), loglik = fit$loglik, logpost = fit$value,
      delta2 = p$delta[2], delta3 = p$delta[3], opener_mode = par[["m1"]], raise_mode = par[["m2"]],
      xi = p$xi, u = p$mx$u, Hu = mix_Hu(p), location = mu, scale = om,
      q999 = mu + om * std_quantile(p, 0.999), endpoint = mu + om * std_endpoint(p),
      gap_all = mean(g, na.rm = TRUE), gap_bracketed = mean(g[brk]), gap_ended_success = mean(g[ended]))
    print(as.data.frame(tail(rows, 1)[[1]] |> mutate(across(where(is.numeric), ~round(.x, 3)))))
  }
}
out <- bind_rows(rows)
write_csv(out |> filter(!free), "out/anchor_mix.csv")
if (any(out$free)) write_csv(out |> filter(free), "out/anchor_mix_free.csv")
