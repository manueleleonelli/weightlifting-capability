## 84_final_select.R
## -----------------------------------------------------------------------
## For each sex and lift, compares the posterior modes of the final model
## found by 52_ext_map.R from the starts of 82_final_starts.R and keeps the
## better of main and alt under the name the other scripts read, with the
## Laplace covariance used to start the chains. This kept mode starts the
## chains and is the mode of the checks at the mode (anchors, shared error,
## predictive checks).
## The other modes (scl, x, and whichever of main and alt was not kept) are a
## search for further basins of the posterior. The second mode is the highest
## of them whose threshold differs from that of the kept mode by more than
## 0.4 standardized units, if any. It is written as ..._second.rds and is used
## by 96_basin_weights.R to tell the basins apart; for the combinations in
## WL_SECOND (e.g. "Women snatch") it also gets its Laplace covariance,
## because a second set of chains starts from it (88_final_demcz.R with
## WL_FROM=second). A mode above the kept one is reported in the column
## higher_elsewhere of out/mode_compare.csv. run_all.R runs this script twice,
## before the chains and again before the weights of the modes, once the
## searches from the third start (which run beside the chains) are done.
## Output: out/mode_compare.csv; out/extmap_e2_<sex>_<lift>_rtau_sm3_mix_fix_log_tau_-3.91202.rds
## (and _main.rds, the main start's mode, when the alt start gives the higher one)
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
only <- Sys.getenv("WL_SEXLIFT", "")   # e.g. "Men snatch" to do one combination
second_for <- trimws(strsplit(Sys.getenv("WL_SECOND", ""), ",")[[1]])
SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
rows <- list()
## Laplace covariance, computed once and stored in the mode's own file as well
with_cov <- function(m, sex, lift, file) {
  if (is.null(m$V)) {
    cat(sex, lift, ": Laplace covariance\n")
    d <- build_data_ext(build_data(att, sex, lift, G = 24), jac = 0)
    m$V <- laplace_cov_ext(m$par, d, m$free); m$sd <- sqrt(diag(m$V))
    saveRDS(m, file)
  }
  m
}
val_of <- function(m) unname(m$value)
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  if (only != "" && only != paste(sex, lift)) next
  base <- sprintf("out/extmap_e2_%s_%s%s", tolower(sex), lift, SUF); fm <- paste0(base, ".rds")
  cur <- readRDS(fm); if (is.null(cur$kind)) cur$kind <- "main"
  ## the candidates; the main start's mode is in _main.rds once the alt one was kept
  cand <- list(main = if (cur$kind == "main") cur else readRDS(paste0(base, "_main.rds")))
  for (k in c("alt", "scl", "x")) { f <- sprintf("%s_%s.rds", base, k); if (file.exists(f)) { m <- readRDS(f); m$kind <- k; cand[[k]] <- m } }
  cand$main$kind <- "main"
  std <- intersect(c("main", "alt"), names(cand))
  best <- std[which.max(sapply(cand[std], val_of))]
  if (best != cur$kind) {
    if (cur$kind == "main") file.copy(fm, paste0(base, "_main.rds"), overwrite = TRUE)
    cur <- cand[[best]]; cat(sex, lift, ": the", best, "start gives the higher mode\n")
  }
  cur <- with_cov(cur, sex, lift, fm); saveRDS(cur, fm)
  ## second mode
  oth <- cand[setdiff(names(cand), best)]
  oth <- oth[vapply(oth, function(m) abs(m$par[["mx_u"]] - cur$par[["mx_u"]]) > 0.4, logical(1))]
  sec <- if (length(oth)) oth[[which.max(sapply(oth, val_of))]] else NULL
  fs <- paste0(base, "_second.rds")
  if (!is.null(sec)) {
    if (paste(sex, lift) %in% second_for) sec <- with_cov(sec, sex, lift, sprintf("%s_%s.rds", base, sec$kind))
    saveRDS(sec, fs)
  } else if (file.exists(fs)) file.remove(fs)
  top <- max(sapply(cand, val_of))
  if (top > val_of(cur) + 0.5) cat(sex, lift, ": a mode", round(top - val_of(cur), 2), "units of log posterior above the kept one, from the",
                                   names(cand)[which.max(sapply(cand, val_of))], "start\n")
  g <- function(k, nm) if (!is.null(cand[[k]])) if (nm == "value") val_of(cand[[k]]) else cand[[k]]$par[[nm]] else NA
  rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, kept = best, second = if (is.null(sec)) NA else sec$kind,
    higher_elsewhere = round(top - val_of(cur), 2),
    logpost_main = g("main", "value"), logpost_alt = g("alt", "value"), logpost_scl = g("scl", "value"), logpost_x = g("x", "value"),
    u_main = g("main", "mx_u"), u_alt = g("alt", "mx_u"), u_scl = g("scl", "mx_u"), u_x = g("x", "mx_u"),
    xi_main = g("main", "xi"), xi_alt = g("alt", "xi"), xi_scl = g("scl", "xi"), xi_x = g("x", "xi"),
    g0_main = g("main", "g0"), g0_alt = g("alt", "g0"), g0_scl = g("scl", "g0"), g0_x = g("x", "g0"))
}
out <- bind_rows(rows)
f <- "out/mode_compare.csv"
if (file.exists(f) && only != "") { old <- read_csv(f, show_col_types = FALSE); old <- old[!(paste(old$sex, old$lift) %in% paste(out$sex, out$lift)), ]; out <- bind_rows(old, out) }
write_csv(out, f); print(as.data.frame(out))
