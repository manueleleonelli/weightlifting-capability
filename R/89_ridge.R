## 89_ridge.R
## -----------------------------------------------------------------------
## The ridge of Figure 2: the exchangeable model (no athlete factor, one lower
## margin scale) with the normal-mixture body, men's snatch, at fixed values
## of the relative attempt scale tau and with tau free.
## The fits are run by 52_ext_map.R. Each fixed-tau fit is started twice, from
## the posterior mode of the final model and from the fit at the neighbouring
## tau (continuation along the ridge, outwards from tau = 0.02), and the
## higher of the two is kept: along a ridge a single start can stop at a
## local mode (on the full archive the fit at tau = 0.05 started from the
## final model ended at xi = -1.02). The fit with tau free starts from the
## kept fit at the largest tau.
## Output, in the names 59_ext_figures.R reads:
##   out/extmap_base_men_snatch_rtau_mix_fix_log_tau_<log tau>.rds, and
##   out/extmap_base_men_snatch_rtau_mix.rds (tau free); the summary rows in
##   out/mapsummary/; the two starts of each fit in out/ridge_starts.csv
RSCRIPT <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
SUF2 <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
mode_file <- sprintf("out/extmap_e2_men_snatch%s.rds", SUF2)
stopifnot(file.exists(mode_file))
base_env <- c(WL_MODEL = "base", WL_SM3 = "0", WL_MIX = "1", WL_G = "24", WL_SEX = "Men", WL_LIFT = "snatch",
              WL_RELTOL = "1e-9", WL_LAPLACE = "0")
dir.create("logs", showWarnings = FALSE)
fit52 <- function(lt, start, tagx, maxit = 300) {
  if (Sys.getenv("RUN_ALL_TEST", "0") == "1") maxit <- 3   # quick test of the pipeline
  tag <- paste0("base_men_snatch_rtau_mix", if (is.na(lt)) "" else paste0("_fix_log_tau_", lt), tagx)
  out <- sprintf("out/extmap_%s.rds", tag)
  if (!file.exists(out)) {
    env <- c(base_env, WL_MAXIT = as.character(maxit), WL_START = start, WL_TAGX = tagx,
             WL_FIX = if (is.na(lt)) "" else paste0("log_tau=", lt))
    do.call(Sys.setenv, as.list(env))
    if (is.na(lt)) Sys.unsetenv("WL_FIX")
    st <- system2(RSCRIPT, "R/52_ext_map.R", stdout = sprintf("logs/ridge_%s.log", tag), stderr = sprintf("logs/ridge_%s.log", tag))
    Sys.unsetenv(names(env))
    if (st != 0 || !file.exists(out)) stop("the fit ", tag, " failed, see logs/ridge_", tag, ".log")
  }
  out
}
keep_best <- function(lt, files) {
  v <- sapply(files, function(f) readRDS(f)$value)
  best <- files[which.max(v)]
  canon <- sprintf("out/extmap_base_men_snatch_rtau_mix%s.rds", if (is.na(lt)) "" else paste0("_fix_log_tau_", lt))
  file.copy(best, canon, overwrite = TRUE)
  ## summary row under the canonical name
  sb <- sub("^out/extmap_", "out/mapsummary/", sub("\\.rds$", ".csv", best)); sc <- sub("^out/extmap_", "out/mapsummary/", sub("\\.rds$", ".csv", canon))
  if (file.exists(sb)) file.copy(sb, sc, overwrite = TRUE)
  rows <- data.frame(tau = if (is.na(lt)) NA else exp(as.numeric(lt)), start = basename(files), logpost = v,
                     xi = sapply(files, function(f) readRDS(f)$par[["xi"]]), kept = files == best, row.names = NULL)
  print(rows)
  rows
}
grid_up <- c("-3.91202", "-3.50656", "-2.99573", "-2.52573")   # 0.02, 0.03, 0.05, 0.08
grid_down <- c("-4.19971", "-4.60517")                          # 0.015, 0.01
log_rows <- list(); prev <- mode_file
for (lt in grid_up) {
  fs <- unique(c(fit52(lt, mode_file, "_m"), if (prev != mode_file) fit52(lt, prev, "_c")))
  log_rows[[length(log_rows) + 1]] <- keep_best(lt, fs)
  prev <- sprintf("out/extmap_base_men_snatch_rtau_mix_fix_log_tau_%s.rds", lt)
}
top <- prev
prev <- "out/extmap_base_men_snatch_rtau_mix_fix_log_tau_-3.91202.rds"
for (lt in grid_down) {
  fs <- c(fit52(lt, mode_file, "_m"), fit52(lt, prev, "_c"))
  log_rows[[length(log_rows) + 1]] <- keep_best(lt, fs)
  prev <- sprintf("out/extmap_base_men_snatch_rtau_mix_fix_log_tau_%s.rds", lt)
}
## tau free, from the kept fit at the largest tau
log_rows[[length(log_rows) + 1]] <- keep_best(NA, fit52(NA, top, "_c", maxit = 600))
write.csv(do.call(rbind, log_rows), "out/ridge_starts.csv", row.names = FALSE)
## the two-start files are kept under their own names; remove their summary rows so that the tables see one row per fit
unlink(Sys.glob("out/mapsummary/base_men_snatch_rtau_mix*_[mc].csv"))
cat("ridge done\n")
