## 59_ext_figures.R
## -----------------------------------------------------------------------
## Figures of the paper (1 to 3) and of the supplement (S1), from the outputs of the final model.
##   fig-jump.pdf      Figure 1: observed and predicted success of the second
##                     attempt after a successful opener, by size of the raise
##   fig-ridge.pdf     Figure 2: the simplified model along the anchor
##   fig-chimodel.pdf  Figure 3: chi(u) on the latent capabilities against the copula
##   fig-trace.pdf     Figure S1: trace plots of the one-in-a-thousand level and the endpoint
suppressPackageStartupMessages({library(dplyr); library(readr)})
source("R/50_ext_lib.R")
dir.create("figs", showWarnings = FALSE, recursive = TRUE)
FIRE <- "firebrick"; GREY <- "grey35"; BLUE <- "steelblue4"
have <- function(f) file.exists(f)
SUF <- Sys.getenv("WL_TAG_SUFFIX", "_rtau_sm3_mix_fix_log_tau_-3.91202")

## ---- success by raise -------------------------------------------------
if (have("out/desc_jump_rates.csv") && have("out/ppc_ext_jumps.csv")) {
  jr <- read_csv("out/desc_jump_rates.csv", show_col_types = FALSE) |>
    filter(situation == "after a make", attempt == 2, jump_band != "0")
  pp <- read_csv("out/ppc_ext_jumps.csv", show_col_types = FALSE) |> filter(situation == "second after a make")
  ## the checks at the mode: where the posterior has two regions, those at the mode of the region with the larger mass
  dom_label <- function(sx, lf) { if (!have("out/basin_weights.csv")) return("final")
    b <- read_csv("out/basin_weights.csv", show_col_types = FALSE); b <- b[b$sex == sx & b$lift == lf, ]
    if (nrow(b) == 2 && isTRUE(all(b$two_regions)) && all(is.finite(b$weight)) && b$basin[which.max(b$weight)] == "second" &&
        any(pp$model == "final_second" & pp$sex == sx & pp$lift == lf)) "final_second" else "final" }
  bands <- c("1", "2", "3", "4", "5", "6-7", "8-10", ">10")
  pdf("figs/fig-jump.pdf", width = 8, height = 3.4)
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.3, 0.7, 0))
  for (sx in c("Men", "Women")) {
    plot(NA, xlim = c(1, length(bands)), ylim = c(0.15, 0.8), xaxt = "n", xlab = "Raise after a successful opener (kg)",
         ylab = "Proportion of second attempts made", main = sx)
    axis(1, at = seq_along(bands), labels = bands, cex.axis = 0.9, gap.axis = 0)   # keep every band label
    for (lf in c("snatch", "cleanjerk")) {
      s <- jr |> filter(sex == sx, lift == lf) |> mutate(x = match(jump_band, bands)) |> arrange(x)
      col <- if (lf == "snatch") FIRE else GREY
      lines(s$x, s$rate, col = col, lwd = 2); points(s$x, s$rate, pch = 19, col = col, cex = 0.5 + 1.2 * sqrt(s$n / max(jr$n)))
      q <- pp |> filter(sex == sx, lift == lf, model == dom_label(sx, lf)) |> mutate(x = match(jump, bands)) |> arrange(x)
      if (nrow(q)) lines(q$x, q$predicted, col = col, lwd = 1.5, lty = 2)
    }
    if (sx == "Men") legend("bottomleft", c("snatch, observed", "clean and jerk, observed", "model, from the declarations"),
                            col = c(FIRE, GREY, GREY), lwd = c(2, 2, 1.5), lty = c(1, 1, 2), pch = c(19, 19, NA), bty = "n", cex = 0.85)
  }
  dev.off()
}

## ---- the ridge: exchangeable model along the anchor -----------------
## the exchangeable model (no athlete factor, one lower margin scale) with the
## normal-mixture body, men's snatch, at fixed tau and with tau free
fs <- Sys.glob("out/extmap_base_men_snatch_rtau_mix_fix_log_tau_*.rds"); fs <- fs[!grepl("ckpt|_[mc]\\.rds$", fs)]   # not the two starts of 89_ridge.R
free_f <- "out/extmap_base_men_snatch_rtau_mix.rds"
if (length(fs) >= 3) {
  rows <- lapply(c(fs, if (have(free_f)) free_f), function(f) {
    x <- readRDS(f); par <- fill_par2(x$par); p <- unpack2(par)
    xd <- design_at(x$ref, 110, 2026); mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
    tibble(tau = p$tau, loglik = x$loglik, m1 = p$dec[1], m2 = p$dec[2], delta2 = p$delta[2], delta3 = p$delta[3], xi = p$xi, u = par[["mx_u"]],
           mu110 = mu, q999 = mu + om * std_quantile(p, 0.999), endpoint = mu + om * std_endpoint(p), free = f == free_f,
           stopped = isTRUE(x$stopped))
  }) |> bind_rows() |> arrange(tau)
  write_csv(rows, "out/ridge_men_snatch.csv")
  rec <- 200
  pdf("figs/fig-ridge.pdf", width = 9, height = 3.2)
  par(mfrow = c(1, 3), mar = c(4, 4.2, 1.5, 1), mgp = c(2.4, 0.7, 0))
  plot(rows$tau, rows$loglik - max(rows$loglik), log = "x", type = "b", pch = ifelse(rows$free, 17, 19), col = FIRE, lwd = 2,
       xlab = expression(paste("Relative attempt scale ", tau)), ylab = "Log likelihood, relative to maximum")
  abline(v = 0.02, col = "grey70", lty = 3)
  plot(rows$tau, 100 * rows$m1, log = "x", type = "b", pch = ifelse(rows$free, 17, 19), col = FIRE, lwd = 2, ylim = c(0, 100 * max(rows$m1) * 1.05),
       xlab = expression(paste("Relative attempt scale ", tau)), ylab = "Margin mode (% of the bar)")
  lines(rows$tau, 100 * rows$m2, type = "b", pch = ifelse(rows$free, 17, 19), col = GREY, lwd = 2)
  abline(v = 0.02, col = "grey70", lty = 3)
  legend("topleft", c("opener", "raise after a success"), col = c(FIRE, GREY), lwd = 2, bty = "n", cex = 0.9)
  yl <- range(c(rows$q999, rows$endpoint[is.finite(rows$endpoint)], rec)) + c(-3, 3)
  plot(rows$tau, rows$endpoint, log = "x", type = "b", pch = ifelse(rows$free, 17, 19), col = FIRE, lwd = 2, ylim = yl,
       xlab = expression(paste("Relative attempt scale ", tau)), ylab = "Capability at 110 kg, 2026 (kg)")
  lines(rows$tau, rows$q999, type = "b", pch = ifelse(rows$free, 17, 19), col = GREY, lwd = 2)
  abline(h = rec, col = "grey50", lty = 2); abline(v = 0.02, col = "grey70", lty = 3)
  legend("topleft", c("endpoint", "one in a thousand", "record of the current era"), col = c(FIRE, GREY, "grey50"), lwd = c(2, 2, 1), lty = c(1, 1, 2), bty = "n", cex = 0.9)
  dev.off()
}


## ---- chi(u) -----------------------------------------------------------
if (have("out/chix_latent.csv")) {
  ch <- read_csv("out/chix_latent.csv", show_col_types = FALSE)
  pdf("figs/fig-chimodel.pdf", width = 8, height = 3.4)
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.3, 0.7, 0))
  for (sx in c("Men", "Women")) {
    a <- ch |> filter(sex == sx, !model) |> arrange(u); b <- ch |> filter(sex == sx, model) |> arrange(u)
    if (!nrow(a)) next
    plot(a$u, a$chi, type = "b", pch = 19, col = FIRE, lwd = 2, ylim = c(0.3, 1), xlab = "quantile level u", ylab = expression(chi(u)), main = sx)
    arrows(a$u, a$lo, a$u, a$hi, angle = 90, code = 3, length = 0.03, col = FIRE)
    lines(b$u, b$chi, type = "b", pch = 17, col = GREY, lwd = 1.5, lty = 2)
    if (sx == "Men") legend("bottomleft", c("latent capabilities", "fitted copula"), col = c(FIRE, GREY), pch = c(19, 17), lty = c(1, 2), bty = "n", cex = 0.85)
  }
  dev.off()
}

## ---- traces ------------------------------------------------------------
## Figure S1: the one-in-a-thousand level and the endpoint of the capability
## distribution at the reference bodyweight in 2026, draw by draw, for every
## chain of the final model (the joint set alone where one was run, since it
## is the set that moves between the regions of the posterior). The endpoint
## panel is cut at 1.5 times the local record of the current era; draws above
## the cut or with an unbounded tail are marked at the top of the panel.
fs <- Sys.glob(sprintf("out/mcmcx_e2_*%s.rds", SUF))
if (length(fs)) {
  ord <- c("men_snatch", "men_cleanjerk", "women_snatch", "women_cleanjerk")
  key <- sub(sprintf("^out/mcmcx_e2_(.*)%s\\.rds$", SUF), "\\1", fs); fs <- fs[order(match(key, ord))]
  rec_tab <- if (have("out/postx_ceiling.csv")) read_csv("out/postx_ceiling.csv", show_col_types = FALSE) else NULL
  MAXD <- 4500   # draws shown per chain (every draw for the men's chains and the joint set)
  cols <- c(FIRE, BLUE, GREY, "darkgoldenrod3", "darkorchid4", "darkgreen", "deeppink3", "black")
  pdf("figs/fig-trace.pdf", width = 8, height = 2.3 * length(fs))
  par(mfrow = c(length(fs), 2), mar = c(3, 3.8, 2, 1), mgp = c(2.2, 0.6, 0))
  for (f in fs) {
    m <- readRDS(f); chs <- if (!is.null(m$chains)) m$chains else list(m$draws)
    if (!is.null(m$chain_from) && any(m$chain_from == "joint")) chs <- chs[m$chain_from == "joint"]
    bw <- if (m$sex == "Men") 110 else 86; xd <- design_at(m$ref, bw, 2026)
    rec <- if (!is.null(rec_tab)) rec_tab$record[rec_tab$sex == m$sex & rec_tab$lift == m$lift & rec_tab$bodyweight == bw][1] else NA
    tr <- lapply(chs, function(D) {
      idx <- unique(round(seq(1, nrow(D), length.out = min(MAXD, nrow(D)))))
      v <- t(sapply(idx, function(i) { p <- unpack2(fill_par2(D[i, ])); mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
        c(q999 = unname(mu + om * std_quantile(p, 0.999)), endpoint = unname(mu + om * std_endpoint(p))) }))
      list(idx = idx, q999 = v[, "q999"], endpoint = v[, "endpoint"]) })
    xl <- c(1, max(sapply(chs, nrow))); lab <- paste(m$sex, if (m$lift == "snatch") "snatch" else "clean and jerk")
    ## the one-in-a-thousand level
    yl <- range(unlist(lapply(tr, `[[`, "q999")))
    plot(NA, xlim = xl, ylim = yl, xlab = "retained draw", ylab = "one-in-a-thousand level (kg)", main = lab, cex.main = 0.9)
    for (k in seq_along(tr)) lines(tr[[k]]$idx, tr[[k]]$q999, col = adjustcolor(cols[(k - 1) %% length(cols) + 1], 0.7), lwd = 0.6)
    ## the endpoint, cut at 1.5 times the local record
    ep <- unlist(lapply(tr, `[[`, "endpoint")); fin <- ep[is.finite(ep)]
    cap <- if (is.finite(rec)) 1.5 * rec else quantile(fin, 0.99)
    top <- min(cap, max(fin)); lo <- min(fin, rec, na.rm = TRUE); yl <- c(lo, top + 0.04 * (top - lo))   # the record is always in the panel
    plot(NA, xlim = xl, ylim = yl, xlab = "retained draw", ylab = "endpoint (kg)", main = lab, cex.main = 0.9)
    for (k in seq_along(tr)) {
      e <- tr[[k]]$endpoint; cc <- adjustcolor(cols[(k - 1) %% length(cols) + 1], 0.7)
      over <- !is.finite(e) | e > top; e[over] <- NA
      lines(tr[[k]]$idx, e, col = cc, lwd = 0.6)
      if (any(over)) points(tr[[k]]$idx[over], rep(yl[2], sum(over)), pch = "|", col = cc, cex = 0.5)
    }
    if (is.finite(rec)) abline(h = rec, col = "grey50", lty = 2)
  }
  dev.off()
}
cat("figures written\n")
