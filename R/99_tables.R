## 99_tables.R
## -----------------------------------------------------------------------
## Writes the body of every table of the paper and of the supplement, as the
## LaTeX rows between \midrule and \bottomrule, to one file per table in the
## folder WL_TABLES (default tables/). Reads the outputs of the main data set
## (WL_RUN_MAIN, default runs/since2017) and of the whole archive
## (WL_RUN_FULL, default runs/full). Run from the project folder after
## run_all.R; it is also the last step of run_all.R. A table whose inputs are
## missing is skipped with a message.
suppressPackageStartupMessages({ library(dplyr); library(readr) })
MAIN <- Sys.getenv("WL_RUN_MAIN", "runs/since2017"); FULL <- Sys.getenv("WL_RUN_FULL", "runs/full")
OUT <- Sys.getenv("WL_TABLES", "tables"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- helpers ---------------------------------------------------------------
rd <- function(run, file) { p <- file.path(run, "out", file); if (!file.exists(p)) stop("missing ", p); read_csv(p, show_col_types = FALSE) }
minus <- function(s) gsub("(?<![\\w$])-(?=[0-9])", "$-$", s, perl = TRUE)
fx <- function(x, d) {                      # fixed decimals, LaTeX minus, infinity
  out <- character(length(x)); r <- round(x, d); r[is.finite(r) & r == 0] <- 0
  out[is.finite(x)] <- minus(formatC(r[is.finite(x)], format = "f", digits = d))
  out[is.infinite(x)] <- "$\\infty$"; out[is.na(x)] <- ""; out
}
big <- function(x, d = 0) minus(formatC(round(x, d), format = "f", digits = d, big.mark = "{,}"))
ci <- function(m, lo, hi, d) paste0(fx(m, d), " (", fx(lo, d), ", ", fx(hi, d), ")")
row <- function(...) paste0(paste(..., sep = " & "), " \\\\")
## the checks at the mode: where the posterior has two regions, those at the mode of the region with the larger mass
dominant <- function(s) { b <- rd(MAIN, "basin_weights.csv"); b <- b[b$sex == s[1] & b$lift == s[2], ]
  if (nrow(b) == 2 && isTRUE(all(b$two_regions)) && all(is.finite(b$weight))) b$basin[which.max(b$weight)] else "kept" }
CB <- list(c("Men", "snatch"), c("Men", "cleanjerk"), c("Women", "snatch"), c("Women", "cleanjerk"))
LIFT <- c(snatch = "Snatch", cleanjerk = "Clean and jerk"); LIFT_SHORT <- c(snatch = "Snatch", cleanjerk = "C\\&J")
REF <- c(Men = 110, Women = 86)
sci1 <- function(x) {                       # one significant figure, as a LaTeX power of ten
  if (!is.finite(x)) return("")
  e <- floor(log10(x)); m <- round(x / 10^e); if (m == 10) { m <- 1; e <- e + 1 }
  if (e >= 0) return(as.character(round(x))); if (m == 1) sprintf("$10^{%d}$", e) else sprintf("$%d \\times 10^{%d}$", m, e)
}
write_tab <- function(name, f) {
  r <- tryCatch(f(), error = function(e) { message("skipped ", name, ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) { writeLines(r, file.path(OUT, paste0(name, ".tex"))); cat("written", name, "\n") }
}
par_get <- function(pp, s, nm) { r <- pp[pp$sex == s[1] & pp$lift == s[2] & pp$parameter == nm, ]; if (!nrow(r)) stop("no ", nm, " for ", paste(s, collapse = " ")); r }
par_row <- function(pp, label, nm, d, scale = 1, interval = TRUE)
  row(label, paste(sapply(CB, function(s) { r <- par_get(pp, s, nm)
    if (interval) ci(scale * r$median, scale * r$lo, scale * r$hi, d) else fx(scale * r$median, d) }), collapse = " & "))

## ---- main text -------------------------------------------------------------
write_tab("table1_coverage", function() {
  cv <- rd(FULL, "desc_coverage.csv")
  row(sub("-", "--", cv$period), cv$competitions, big(cv$athlete_day_lifts), big(cv$athletes))
})

write_tab("table2_openers", function() {
  op <- rd(MAIN, "desc_opener_prev.csv"); bands <- unique(op$band)
  lab <- function(b) if (startsWith(b, "<")) sprintf("below $%s$", sub("<", "", b)) else if (startsWith(b, ">")) sprintf("above $%s$", sub(">", "", b)) else
    do.call(sprintf, c(list("$%s$ to $%s$"), as.list(strsplit(b, "-")[[1]])))
  sapply(bands, function(b) row(lab(b), paste(sapply(CB, function(s) { r <- op[op$sex == s[1] & op$lift == s[2] & op$band == b, ]
    paste(big(r$n), fx(r$opener_success, 2), sep = " & ") }), collapse = " & ")))
})

write_tab("table3_attempts", function() {
  pp <- rd(MAIN, "postx_params.csv")
  c(par_row(pp, "$\\delta_2$ & ", "delta2", 2), par_row(pp, "$\\delta_3$ & ", "delta3", 2),
    par_row(pp, "Opener & mode", "m1", 1, 100), par_row(pp, " & upper scale", "sp1", 1, 100), par_row(pp, " & lower scale", "sm1", 1, 100),
    par_row(pp, "Raise after a success & mode", "m2", 1, 100), par_row(pp, " & upper scale", "sp2", 1, 100), par_row(pp, " & lower scale", "sm2", 1, 100))
})

write_tab("table4_tail", function() {
  pp <- rd(MAIN, "postx_params.csv")
  c(par_row(pp, "$u_l$", "mx_u", 2), par_row(pp, "$H_l(u_l)$", "Hu", 3), par_row(pp, "$\\sigma_l$", "sigma_tail", 2),
    par_row(pp, "$\\xi_l$", "xi", 2), par_row(pp, "$\\rho_{a,l}$", "rhoa", 3), par_row(pp, "$P(\\xi_l < 0)$", "P(xi<0)", 2, interval = FALSE),
    par_row(pp, "Location trend (kg per year)", "b3", 2))
})

ladder_rows <- function(run, with_p) {
  ls <- rd(run, "ladder_summary.csv"); pp <- rd(run, "postx_params.csv")
  at <- function(s, rung) { r <- ls[ls$sex == s[1] & ls$lift == s[2] & ls$rung == rung & ls$bodyweight == REF[[s[1]]], ]; if (nrow(r) != 1) stop("ladder row ", rung); r }
  cell <- function(rung, what) paste(sapply(CB, function(s) {
    r <- at(s, rung)
    if (rung == "L3" && what == "xi") { x <- par_get(pp, s, "xi"); return(ci(x$median, x$lo, x$hi, 2)) }
    if (rung == "L3" && what == "p") return(fx(par_get(pp, s, "P(xi<0)")$median, 2))
    switch(what, xi = ci(r$xi, r$xi_lo, r$xi_hi, 2), p = fx(r$p_finite, 2),
           ceiling = if (rung == "POT") fx(r$ceiling, 0) else sprintf("%s (%s)", fx(r$ceiling, 0), fx(r$p10, 2))) }), collapse = " & ")
  out <- c()
  for (k in list(c("L1", "Record book"), c("L2", "Record book, athlete factor"), c("L3", "Full model"), c("POT", "No body (check)"))) {
    out <- c(out, row(k[2], "$\\xi_l$", cell(k[1], "xi")))
    if (with_p) out <- c(out, row("", "$P(\\xi_l < 0)$", cell(k[1], "p")))
    out <- c(out, row("", "ceiling", cell(k[1], "ceiling")))
  }
  out
}
write_tab("table5_ladder", function() ladder_rows(MAIN, TRUE))

write_tab("table6_quantiles", function() {
  qq <- rd(MAIN, "postx_quantiles.csv"); cc <- rd(MAIN, "postx_ceiling.csv"); out <- c()
  for (sx in c("Men", "Women")) {
    if (sx == "Women") out <- c(out, "\\addlinespace")
    out <- c(out, sprintf("\\multicolumn{9}{l}{\\textit{%s}} \\\\", sx))
    for (bw in sort(unique(cc$bodyweight[cc$sex == sx]))) {
      cells <- c(bw)
      for (lf in c("snatch", "cleanjerk")) {
        q <- function(l) { r <- qq[qq$sex == sx & qq$lift == lf & qq$bodyweight == bw & !is.na(qq$level) & abs(qq$level - l) < 1e-9, ]; ci(r$median, r$lo, r$hi, 1) }
        c0 <- cc[cc$sex == sx & cc$lift == lf & cc$bodyweight == bw, ]
        cells <- c(cells, q(0.99), q(0.999), fx(c0$ceiling_median, 0), c0$record)
      }
      out <- c(out, do.call(row, as.list(cells)))
    }
  }
  out
})

write_tab("table7_ceiling", function() {
  cc <- rd(MAIN, "postx_ceiling.csv"); pp <- rd(MAIN, "postx_params.csv"); am <- rd(MAIN, "anchor_mix.csv")
  ## the modes at the three anchors: those of the region with the larger posterior mass, where there are two
  if (!"mode" %in% names(am)) am$mode <- "kept"; dom <- dominant
  BW <- list(Men = c(73, 89, 110, 150), Women = c(59, 71, 86, 120)); out <- c()
  for (s in CB) for (bw in BW[[s[1]]]) {
    r <- cc[cc$sex == s[1] & cc$lift == s[2] & cc$bodyweight == bw, ]
    modes <- if (bw == REF[[s[1]]]) sapply(c(0.01, 0.02, 0.03), function(t) fx(am$endpoint[am$sex == s[1] & am$lift == s[2] & am$mode == dom(s) & abs(am$tau - t) < 1e-4], 0)) else c("", "", "")
    out <- c(out, row(s[1], LIFT_SHORT[[s[2]]], bw, r$record, fx(par_get(pp, s, "P(xi<0)")$median, 2), fx(r$p_within_10, 2), fx(r$p_within_25, 2),
                      fx(r$p_within_50, 2), fx(r$ceiling_median, 0), modes[1], modes[2], modes[3]))
  }
  out
})

## ---- supplement -------------------------------------------------------------
write_tab("tableS1_quadrature", function() {
  qc <- rd(MAIN, "quadrature_check.csv")
  rule <- c("Rule of S1.1", "$12$ Gauss--Hermite nodes", "$12$ Gauss--Hermite nodes", "$K_\\eta = 4$", "$K_\\eta = 4$")
  ref <- c("Grid of $400{,}001$ points", "$20$ nodes", "Grid of $4{,}001$ points", "Grid of $200{,}001$ points", "$K_\\eta = 64$")
  int <- c("Capability", "Athlete factor", "Athlete factor", "Shared error", "Shared error")
  per <- ifelse(qc$per == "total", "in total", paste("per", qc$per))
  row(int, rule, ref, paste(sapply(qc$largest_difference, sci1), per))
})

write_tab("tableS2_ridge", function() {
  rg <- rd(MAIN, "ridge_men_snatch.csv") |> arrange(free, tau)
  l0 <- rg$loglik[!rg$free & abs(rg$tau - 0.02) < 1e-4]
  lab <- ifelse(rg$free, sprintf("free (%.3f)", rg$tau), sprintf("%.3f", rg$tau))
  row(lab, fx(rg$loglik - l0, 0), fx(100 * rg$m1, 1), fx(100 * rg$m2, 1), fx(rg$delta2, 2), fx(rg$delta3, 2), fx(rg$xi, 2), fx(rg$u, 2), fx(rg$q999, 1), fx(rg$endpoint, 1))
})

write_tab("tableS3_simulation", function() {
  sf <- rd(MAIN, "sim_fits.csv")
  q <- list(c("xi", "Tail shape $\\xi$", 3), c("u", "Threshold $u$", 2), c("rhoa", "Athlete correlation $\\rho_a$", 3),
            c("delta2", "Shift $\\delta_2$", 2), c("delta3", "Shift $\\delta_3$", 2), c("opener_mode", "Opener mode $m_1$ (\\%)", 2),
            c("raise_mode", "Raise mode $m_2$ (\\%)", 2), c("q999", "$q_{0.999}$ (kg)", 1), c("endpoint", "Endpoint (kg)", 1))
  out <- c()
  for (s in list(c("Men", "snatch", "Men's snatch"), c("Women", "cleanjerk", "Women's clean and jerk"))) {
    g <- sf[sf$sex == s[1] & sf$lift == s[2], ]; if (!nrow(g)) next
    tr <- g[g$what == "truth", ][1, ]; fi <- g[g$what == "fit", ]; a <- fi[abs(fi$fit_tau - 0.02) < 1e-6, ]
    out <- c(out, sprintf("\\multicolumn{5}{l}{\\emph{%s}} \\\\", s[3]))
    for (k in q) {
      nm <- k[1]; d <- as.integer(k[3]); x <- a[[nm]]
      fit <- if (any(!is.finite(x))) paste0(fx(median(x), d), "$^{*}$") else sprintf("%s (%s)", fx(mean(x), d), fx(sd(x), d))
      ch <- sapply(c(0.01, 0.03), function(t) {
        b <- fi[abs(fi$fit_tau - t) < 1e-6, ]; base <- a[[nm]][match(b$rep, a$rep)]; dd <- b[[nm]] - base; ok <- is.finite(dd)
        if (!length(dd)) return("") else paste0(fx(mean(dd[ok]), d), if (sum(ok) < length(dd)) "$^{\\dagger}$" else "") })
      out <- c(out, row(k[2], fx(tr[[nm]], d), fit, ch[1], ch[2]))
    }
  }
  out
})

write_tab("tableS4_bodyfree", function() {
  bp <- rd(MAIN, "bodyfree_pot.csv"); out <- c()
  for (s in CB) for (q in c(0.90, 0.95)) {
    r <- bp[bp$sex == s[1] & bp$lift == s[2] & abs(bp$q - q) < 1e-9 & bp$bodyweight == REF[[s[1]]], ]
    out <- c(out, row(if (q == 0.90) s[1] else "", if (q == 0.90) LIFT[[s[2]]] else "", sprintf("%.2f", q), big(r$n_exc),
                      ci(r$xi, r$xi_lo, r$xi_hi, 2), fx(r$endpoint, 0), if (q == 0.90) r$record else ""))
  }
  out
})

write_tab("tableS5_fixed_threshold", function() {
  mg <- bind_rows(rd(MAIN, "mixture_gpd_men.csv"), rd(MAIN, "mixture_gpd_women.csv")) |> filter(model == "MNGPD")
  out <- c()
  for (q in c(0.85, 0.90, 0.95)) {
    v <- lapply(CB, function(s) mg[mg$sex == s[1] & mg$lift == s[2] & abs(mg$p_u - q) < 1e-9, ])
    out <- c(out, row(sprintf("Threshold at %.2f", q), "$\\hat\\xi$", paste(sapply(v, function(r) fx(r$xi, 2)), collapse = " & ")),
             row("", "ceiling", paste(sapply(v, function(r) fx(r$endpoint, 0)), collapse = " & ")))
  }
  out
})

write_tab("tableS6_ppc_attempts", function() {
  pa <- rd(MAIN, "ppc_ext_attempts.csv"); out <- c()
  for (s in CB) for (j in 1:3) { lab <- if (dominant(s) == "second" && any(pa$model == "final_second" & pa$sex == s[1] & pa$lift == s[2])) "final_second" else "final"
    r <- pa[pa$sex == s[1] & pa$lift == s[2] & pa$attempt == j & pa$model == lab, ]
    out <- c(out, row(s[1], LIFT[[s[2]]], j, sprintf("%.3f", r$observed), sprintf("%.3f", r$predicted))) }
  out
})

write_tab("tableS7_ppc_records", function() {
  b <- rd(MAIN, "ppc_record_bw.csv"); out <- c()
  for (sx in c("Men", "Women")) {
    out <- c(out, sprintf("\\multicolumn{5}{l}{\\emph{%s}} \\\\", sx))
    for (bw in sort(unique(b$bodyweight[b$sex == sx]))) {
      cells <- c(bw)
      for (lf in c("snatch", "cleanjerk")) { r <- b[b$sex == sx & b$lift == lf & b$bodyweight == bw, ]
        cells <- c(cells, r$record, sprintf("%s (%s, %s)", fx(r$sim_median, 0), fx(r$sim_lo, 0), fx(r$sim_hi, 0))) }
      out <- c(out, do.call(row, as.list(cells)))
    }
  }
  out
})

write_tab("tableS8_raise_after_miss", function() {
  pp <- rd(MAIN, "postx_params.csv")
  c(par_row(pp, "Raise after a miss & mode", "m3", 1, 100), par_row(pp, " & upper scale", "sp3", 1, 100), par_row(pp, " & lower scale", "sm3", 1, 100))
})

write_tab("tableS9_anchors", function() {
  am <- rd(MAIN, "anchor_mix.csv"); fr <- tryCatch(rd(MAIN, "anchor_mix_free.csv"), error = function(e) NULL)
  if (!"mode" %in% names(am)) am$mode <- "kept"
  line <- function(r, lab) row(r$sex, paste0(LIFT[[r$lift]], if (!is.null(r$mode) && r$mode == "second") ", second mode" else ""), lab,
                               big(r$loglik, 1), fx(100 * r$opener_mode, 1), fx(r$delta2, 2), fx(r$delta3, 2),
                               fx(r$xi, 2), fx(r$q999, 1), fx(r$endpoint, 1), fx(r$gap_all, 1))
  out <- c()
  for (s in CB) for (md in c("kept", "second")) { a <- am[am$sex == s[1] & am$lift == s[2] & am$mode == md, ] |> arrange(tau)
    for (i in seq_len(nrow(a))) out <- c(out, line(a[i, ], sprintf("$\\tau = %.2f$", a$tau[i]))) }
  if (!is.null(fr)) for (i in seq_len(nrow(fr))) out <- c(out, line(fr[i, ], sprintf("$\\tau$ free ($\\hat\\tau = %.3f$)", fr$tau[i])))
  out
})

write_tab("tableS10_priors", function() {
  pr <- rd(MAIN, "prior_reweight.csv")
  lab <- c(reference = "Reference", xi_wide = "$\\xi_l \\sim N(0, 0.5^2)$", u_wide = "$u_l \\sim N(1.5, 2^2)$", u_max3 = "$u_l < 3$")
  out <- c()
  for (s in CB) { r <- pr[pr$sex == s[1] & pr$lift == s[2], ]; r <- r[match(names(lab), r$variant), ]
    for (i in seq_len(nrow(r))) out <- c(out, row(if (i == 1) s[1] else "", if (i == 1) LIFT[[s[2]]] else "", lab[[r$variant[i]]], big(r$ess[i]),
      fx(r$xi_median[i], 2), fx(r$u_median[i], 2), fx(r$q999_median[i], 1), sprintf("%.3f", r$p_finite[i]),
      paste0(fx(r$endpoint_median[i], 1), " (", fx(r$endpoint_lo[i], 1), ", ", fx(r$endpoint_hi[i], 1), ")"))) }
  out
})

write_tab("tableS11_scale_trend", function() {
  st <- rd(MAIN, "scale_trend_check.csv") |> filter(year == 2026)
  sapply(CB, function(s) { r <- st[st$sex == s[1] & st$lift == s[2], ]
    row(s[1], LIFT[[s[2]]], fx(r$loglik1 - r$loglik0, 1), fx(r$g3, 4), fx(r$xi0, 2), fx(r$xi1, 2), fx(r$q999_0, 1), fx(r$q999_1, 1), fx(r$endpoint0, 1), fx(r$endpoint1, 1)) })
})

write_tab("tableS12_shared_error", function() {
  se <- rd(MAIN, "shared_error.csv"); if (!"mode" %in% names(se)) se$mode <- "kept"
  unlist(lapply(CB, function(s) { rr <- se[se$sex == s[1] & se$lift == s[2], ]; rr <- rr[order(rr$mode != "kept"), ]
    sapply(seq_len(nrow(rr)), function(j) { r <- rr[j, ]
      row(s[1], paste0(LIFT[[s[2]]], if (r$mode == "second") ", second mode" else ""), fx(r$gain, 0), fx(r$s_eta, 3), fx(r$xi_e2, 2), fx(r$xi_e12, 2),
          fx(r$q999_e2, 1), fx(r$q999_e12, 1), fx(r$endpoint_e2, 1), fx(r$endpoint_e12, 1)) }) }))
})

write_tab("tableS13_convergence", function() {
  dg <- rd(MAIN, "diag_summary.csv"); out <- c()
  for (s in CB) for (k in list(c("final", "Final model"), c("L1", "Record book"), c("L2", "Record book, athlete factor"))) {
    nm <- paste(k[1], s[1], s[2])
    rr <- dg[dg$fit == nm | startsWith(dg$fit, paste0(nm, ", ")), ]
    for (j in seq_len(nrow(rr))) { r <- rr[j, ]
      lab <- if (r$fit == nm) k[2] else paste0(k[2], ", ", sub(paste0(nm, ", "), "", r$fit, fixed = TRUE), if (is.finite(r$weight)) paste0(" (", fx(r$weight, 2), ")") else "")
      acc <- mean(as.numeric(strsplit(r$acceptance, " ")[[1]]))
      out <- c(out, row(s[1], LIFT[[s[2]]], lab, fx(acc, 2), fx(r$max_rhat, 2), fx(r$min_ess_bulk, 0), fx(r$median_ess_bulk, 0),
                        fx(r$ess_xi, 0), fx(r$ess_u, 0), fx(r$ess_endpoint, 0), fx(r$rhat_endpoint, 2))) }
  }
  out
})

write_tab("tableS14_archive", function() {
  out <- c()
  for (s in CB) for (k in list(c(MAIN, "since 2017"), c(FULL, "since 1998"))) {
    pp <- rd(k[1], "postx_params.csv"); ls <- rd(k[1], "ladder_summary.csv"); ut <- rd(k[1], "desc_unit_types.csv")
    cr <- ls[ls$sex == s[1] & ls$lift == s[2] & ls$rung == "L3" & ls$bodyweight == REF[[s[1]]], ]; xi <- par_get(pp, s, "xi")
    out <- c(out, row(if (k[2] == "since 2017") s[1] else "", if (k[2] == "since 2017") LIFT[[s[2]]] else "", k[2],
                      big(sum(ut$n[ut$sex == s[1] & ut$lift == s[2]])), ci(xi$median, xi$lo, xi$hi, 2), fx(par_get(pp, s, "P(xi<0)")$median, 2),
                      fx(cr$q999, 1), fx(cr$ceiling, 0), fx(par_get(pp, s, "rhoa")$median, 3), fx(par_get(pp, s, "b3")$median, 2),
                      fx(par_get(pp, s, "delta2")$median, 2), fx(par_get(pp, s, "delta3")$median, 2)))
  }
  out
})

write_tab("tableS15_archive_ladder", function() ladder_rows(FULL, FALSE))
