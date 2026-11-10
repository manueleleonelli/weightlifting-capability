## 82_final_starts.R
## -----------------------------------------------------------------------
## Starting points for the posterior mode of the final model (athlete factor,
## situation-specific lower margin scales, normal-mixture body with GPD tail
## above an estimated threshold, tau = 0.02), one per sex and lift and up to
## four kinds:
##   main  location, scale, mixture and tail from the best-lift fit of
##         66_mix_start.R on this data set; attempt and declaration
##         parameters from R/attempt_starts.rds (the modes of the full
##         archive; they are proportions of the bar and do not depend on the
##         design); rho_a from the same file
##   alt   as main with the threshold at 2.2 and a positive tail shape, to
##         check whether the posterior has a second mode at a high threshold
##   scl   as main moved along the ridge between the scale of capability and
##         the free mixture components: the log scale g0 is lowered by 0.4 and
##         the free components, the threshold and the tail scale are widened
##         by the factor exp(0.4) in standardized units, so that they describe
##         about the same capabilities in kilograms. It searches for modes with
##         a narrower scale of capability that the other two starts miss.
##   x     a mode found in our runs that the optimization from scl reaches in
##         some runs and not in others (the path of the optimizer depends on
##         rounding in the likelihood), stored in R/mode_starts.rds so that
##         every run finds it: the women's snatch since 2017, where the
##         chains around it give the second set of 88_final_demcz.R.
##         Written only for the data set named in WL_DATASET (run_all.R sets it).
## The mode used for the chains and for the checks at the mode is the better
## of main and alt; scl and x only look for further modes (84_final_select.R).
## Output: out/start_final_<sex>_<lift>_<kind>.rds (list with $par)
source("R/50_ext_lib.R"); source("R/81_ladder_lib.R")
ast <- readRDS("R/attempt_starts.rds")
stored <- if (file.exists("R/mode_starts.rds")) readRDS("R/mode_starts.rds")[[Sys.getenv("WL_DATASET", "")]] else NULL
along_scale <- function(p, f = -0.4) {
  p[["g0"]] <- p[["g0"]] + f; s <- exp(-f)
  for (k in 2:3) { p[[paste0("mx_m", k)]] <- p[[paste0("mx_m", k)]] * s; p[[paste0("mx_ls", k)]] <- p[[paste0("mx_ls", k)]] - f }
  p[["mx_u"]] <- p[["mx_u"]] * s; p[["mx_lsig"]] <- p[["mx_lsig"]] - f
  p
}
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  ms <- sanitize_mix(readRDS(sprintf("out/mixstart_%s_%s.rds", tolower(sex), lift))$par)
  a <- ast[[paste(tolower(sex), lift, sep = "_")]]
  par <- fill_par2(c(ms, a, log_tau = log(0.02)))
  par[c("kappa", "lam", "lnu", "lse", "lpi")] <- NA
  saveRDS(list(par = par), sprintf("out/start_final_%s_%s_main.rds", tolower(sex), lift))
  alt <- par; alt[["mx_u"]] <- 2.2; alt[["xi"]] <- 0.05; alt[["mx_lsig"]] <- log(0.4)
  saveRDS(list(par = alt), sprintf("out/start_final_%s_%s_alt.rds", tolower(sex), lift))
  saveRDS(list(par = along_scale(par)), sprintf("out/start_final_%s_%s_scl.rds", tolower(sex), lift))
  x <- stored[[paste(sex, lift)]]
  if (!is.null(x)) { x <- fill_par2(x); x[c("kappa", "lam", "lnu", "lse", "lpi")] <- NA
    saveRDS(list(par = x), sprintf("out/start_final_%s_%s_x.rds", tolower(sex), lift)); cat(sex, lift, ": stored start written\n") }
}
