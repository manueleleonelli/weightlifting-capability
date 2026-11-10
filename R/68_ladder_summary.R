## 68_ladder_summary.R
## -----------------------------------------------------------------------
## Summary of the ladder (Table 5 of the paper, Table S15 of the supplement): tail shape and endpoint at
## the reference bodyweight in 2026 under
##   L1  best lifts, days exchangeable (67_ladder.R)
##   L2  best lifts, athlete factor (67_ladder.R)
##   L3  the full latent model (54_ext_mcmc.R with WL_MIX=1)
##   POT the body-free check (63_bodyfree_pot.R, threshold quantile 0.9)
## Output: out/ladder_summary.csv and LaTeX rows in out/ladder_rows.tex
source("R/50_ext_lib.R")
SUF3 <- Sys.getenv("WL_TAG_SUFFIX", "_rtau_sm3_mix_fix_log_tau_-3.91202")
NDRAW <- as.integer(Sys.getenv("WL_NDRAW", "2000"))
rec_tab <- read_csv("out/postx_ceiling.csv", show_col_types = FALSE)
pot <- read_csv("out/bodyfree_pot.csv", show_col_types = FALSE)
summ_draws <- function(D, ref, bw, record) {
  idx <- round(seq(1, nrow(D), length.out = min(NDRAW, nrow(D))))
  xd <- design_at(ref, bw, 2026)
  vals <- t(sapply(idx, function(i) {
    p <- unpack2(fill_par2(D[i, ]))
    mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
    c(xi = unname(p$xi), ep = unname(mu + om * std_endpoint(p)), q999 = unname(mu + om * std_quantile(p, 0.999)))
  }))
  ep <- vals[, "ep"]
  tibble(xi = median(vals[, "xi"]), xi_lo = quantile(vals[, "xi"], .025), xi_hi = quantile(vals[, "xi"], .975),
         p_finite = mean(is.finite(ep)), ceiling = median(ep), q999 = median(vals[, "q999"]),
         ceiling_lo = quantile(ep[is.finite(ep)], .025), ceiling_hi = quantile(ep[is.finite(ep)], .975),
         p10 = mean(ep <= 1.10 * record), p25 = mean(ep <= 1.25 * record), p50 = mean(ep <= 1.50 * record))
}
rows <- list()
ROWS_ONLY <- Sys.getenv("WL_ROWS_ONLY", "0") == "1"   # rewrite the LaTeX rows from out/ladder_summary.csv
if (!ROWS_ONLY) for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  bw <- if (sex == "Men") 110 else 86
  record <- rec_tab$record[rec_tab$sex == sex & rec_tab$lift == lift & rec_tab$bodyweight == bw][1]
  for (rung in c("L1", "L2", "L3")) {
    f <- if (rung == "L3") sprintf("out/mcmcx_e2_%s_%s%s.rds", tolower(sex), lift, SUF3) else sprintf("out/ladder_%s_%s_%s.rds", rung, tolower(sex), lift)
    if (!file.exists(f)) { cat("missing", f, "\n"); next }
    x <- readRDS(f)
    rows[[length(rows) + 1]] <- bind_cols(tibble(sex = sex, lift = lift, rung = rung, bodyweight = bw, record = record),
                                          summ_draws(x$draws, x$ref, bw, record))
  }
  pr <- pot |> filter(sex == !!sex, lift == !!lift, q == 0.9, bodyweight == bw)
  rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, rung = "POT", bodyweight = bw, record = record,
    xi = pr$xi, xi_lo = pr$xi_lo, xi_hi = pr$xi_hi, p_finite = pr$p_finite, ceiling = pr$endpoint, q999 = NA,
    ceiling_lo = pr$ep_lo, ceiling_hi = pr$ep_hi,
    p10 = NA, p25 = NA, p50 = NA)
}
if (ROWS_ONLY) out <- read_csv("out/ladder_summary.csv", show_col_types = FALSE) else { out <- bind_rows(rows); write_csv(out, "out/ladder_summary.csv") }
print(as.data.frame(out |> mutate(across(where(is.numeric), ~round(.x, 3)))))
## LaTeX rows
lab <- c(L1 = "Record book", L2 = "Record book, athlete factor", L3 = "Full model", POT = "No body (check)")
fm <- function(x) gsub("-", "$-$", formatC(x, format = "f", digits = 2))
cmb <- list(c("Men", "snatch"), c("Men", "cleanjerk"), c("Women", "snatch"), c("Women", "cleanjerk"))
lines <- c()
for (rg in c("L1", "L2", "L3", "POT")) {
  cx <- sapply(cmb, function(cb) { r <- out |> filter(sex == cb[1], lift == cb[2], rung == rg); if (!nrow(r)) "--" else sprintf("%s (%s, %s)", fm(r$xi), fm(r$xi_lo), fm(r$xi_hi)) })
  ce <- sapply(cmb, function(cb) { r <- out |> filter(sex == cb[1], lift == cb[2], rung == rg); if (!nrow(r)) "--" else if (is.na(r$p10)) sprintf("%.0f", r$ceiling) else if (!is.finite(r$ceiling)) sprintf("$\\infty$ (%.2f)", r$p10) else sprintf("%.0f (%.2f)", r$ceiling, r$p10) })
  pf <- sapply(cmb, function(cb) { r <- out |> filter(sex == cb[1], lift == cb[2], rung == rg); if (!nrow(r) || is.na(r$p_finite)) "--" else sprintf("%.2f", r$p_finite) })
  lines <- c(lines, sprintf("%s & $\\xi_l$ & %s \\\\", lab[[rg]], paste(cx, collapse = " & ")),
             sprintf(" & $P(\\xi_l < 0)$ & %s \\\\", paste(pf, collapse = " & ")),
             sprintf(" & ceiling & %s \\\\", paste(ce, collapse = " & ")))
}
writeLines(lines, "out/ladder_rows.tex"); cat(lines, sep = "\n")
