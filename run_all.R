## =======================================================================
## run_all.R
## Every result of "How much can a human being lift? A Bayesian latent
## capability model for Olympic weightlifting" and of its supplement, for
## the athlete-days from 2017 onwards (the analysis of the paper) and for
## the whole archive from 1998 (Section S7 of the supplement).
##
## Run from this folder, in a terminal:   Rscript run_all.R
## or in RStudio with this folder as the working directory: source("run_all.R")
##
## The run is long (see README.md and the time estimate printed at the
## start). It can be stopped at any moment and restarted with the same
## command: finished steps are skipped and chains resume from their last
## checkpoint. Each data set gets its own folder, runs/since2017 and
## runs/full, with its out/ (results, chains), logs/ and figs/. The tables
## of the paper are written to tables/ (see README.md).
## =======================================================================

## ---------------- settings ------------------------------------------------
N_CORES       <- max(1, parallel::detectCores(logical = FALSE))   # parallel jobs: one per physical core
DATASETS      <- c(full = 0, since2017 = 2017)   # name = first year included (0: all)
MAIN          <- "since2017"  # the data set of the paper: the ridge (Figure 2) and the attempt-model variants are fitted on it

## final model: 4 chains per sex and lift, run together by differential-evolution
## Markov chains (DE-MCzs, R/80_sampler.R); one generation updates every
## chain once, so it costs CHAINS likelihood evaluations.
## ITER and BURN are given per data set. They can be changed between restarts:
## a resumed sampler continues to the new ITER (or stops at once if it is
## already past it) and the burn-in is applied at the end, so nothing is lost.
SAMPLER       <- "demcz"    # "demcz", or "amwg" (adaptive Metropolis, one job per chain)
COORD         <- "sstar"    # sampler coordinates (see R/80_sampler.R)
CHAINS        <- 4
ITER          <- c(full = 120000, since2017 = 120000)   # generations (iterations per chain)
## the women's lifts since 2017 put part of their mass on shapes at or above zero and the
## chains move slowly between finite and infinite endpoints (Section S6): they run longer
LONG          <- list(full = character(0), since2017 = c("Women snatch", "Women cleanjerk"))
ITER_LONG     <- c(full = 120000, since2017 = 360000)
## the women's snatch since 2017 has a second mode (Section S6): a second set of CHAINS chains
## starts from it
TWO_MODES     <- list(full = character(0), since2017 = "Women snatch")
ITER_SECOND   <- c(full = 120000, since2017 = 120000)
## where the chains of a combination sit in two regions of the posterior and move between them only
## slowly (R/96_basin_weights.R finds this from the chains above), a joint set of chains starts from
## their last states with an archive of draws from both regions and frequent jumps of full length,
## so that it moves between the regions every few hundred generations; the share of its draws in
## each region is the weight of the region in the pooled posterior (Section S6)
ITER_JOINT    <- c(full = 40000, since2017 = 40000)
BURN          <- c(full = 30000, since2017 = 30000)     # discarded
THIN          <- 20         # every 20th retained state of each chain is stored

## record-book rungs of the ladder (cheap): 4 chains per rung, sex and lift
LADDER_CHAINS <- 4
LADDER_ITER   <- 200000
LADDER_BURN   <- 50000
LADDER_THIN   <- 50

BOOT          <- 1000       # athlete bootstrap of the body-free check
COPULA_DRAWS  <- 5          # marginal posterior draws in the copula stage
NREP_SIM      <- 10         # simulation study (Section S2.2): data sets per sex and lift, fitted at tau = 0.02,
NREP_ANCHOR   <- 5          # of which the first NREP_ANCHOR are also fitted at tau = 0.01 and 0.03
NDRAW_PPC     <- 200        # posterior draws in the predictive check of the record book (Section S5)
RUN_E1        <- FALSE      # a shared self-assessment error in the simplified model, men's snatch; not
                            # reported in the paper (about 10 s per likelihood evaluation)
RUN_E12       <- TRUE       # the final model with the shared self-assessment error added, posterior
                            # modes for the four combinations (Section S5, Table S12); about 5 s per
                            # likelihood evaluation, 6 to 12 hours per job, the four jobs run in parallel
## ---------------------------------------------------------------------------

## a quick end-to-end test of the whole pipeline on a small part of the data:
##   RUN_ALL_TEST=1 Rscript run_all.R      (results meaningless; checks that every step runs)
TEST <- Sys.getenv("RUN_ALL_TEST", "0") == "1"
if (TEST) {
  DATASETS <- c(full = 2022, since2017 = 2024); CHAINS <- 3; ITER <- c(full = 200, since2017 = 300); BURN <- c(full = 100, since2017 = 100); THIN <- 10
  ITER_LONG <- c(full = 200, since2017 = 400); ITER_SECOND <- c(full = 200, since2017 = 300); ITER_JOINT <- c(full = 200, since2017 = 200); NREP_SIM <- 1; NREP_ANCHOR <- 1; NDRAW_PPC <- 3
  LADDER_CHAINS <- 3; LADDER_ITER <- 1500; LADDER_BURN <- 500; LADDER_THIN <- 10; BOOT <- 5; COPULA_DRAWS <- 1
}

PROJECT <- normalizePath(".")
stopifnot(file.exists(file.path(PROJECT, "R", "50_ext_lib.R")))
RSCRIPT <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
need <- c("Rcpp", "dplyr", "tidyr", "readr", "quantreg", "numDeriv", "MASS", "posterior")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss)) { cat("installing:", miss, "\n"); install.packages(miss, repos = "https://cloud.r-project.org") }

log_msg <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%a %H:%M"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = file.path(PROJECT, "runs", "run_all.log"), append = TRUE) }
dir.create(file.path(PROJECT, "runs"), showWarnings = FALSE)

## ---------------- data sets ---------------------------------------------
setup_run <- function(name, from_year) {
  rd <- file.path(PROJECT, "runs", name)
  for (sub in c("data", "out", "logs", "figs")) dir.create(file.path(rd, sub), recursive = TRUE, showWarnings = FALSE)
  ## code: a copy, so that each run is self-contained
  dir.create(file.path(rd, "R"), showWarnings = FALSE)
  file.copy(list.files(file.path(PROJECT, "R"), full.names = TRUE), file.path(rd, "R"), overwrite = TRUE)
  df <- file.path(rd, "data", "attempts_clean.csv")
  if (!file.exists(df)) {
    src <- c(file.path(PROJECT, "data", "attempts_clean.csv"), file.path(PROJECT, "data", "attempts_clean.csv.gz"))
    src <- src[file.exists(src)][1]; if (is.na(src)) stop("data/attempts_clean.csv(.gz) not found")
    att <- readr::read_csv(src, show_col_types = FALSE, col_types = readr::cols(category = readr::col_character(), .default = readr::col_guess()))
    att <- att[att$year >= from_year, ]
    readr::write_csv(att, df)
  }
  rd
}
RUNS <- sapply(names(DATASETS), function(n) setup_run(if (TEST) paste0("test_", n) else n, DATASETS[[n]]))
names(RUNS) <- names(DATASETS)

## ---------------- jobs ----------------------------------------------------
## a job runs one script with Rscript in a run folder, with settings passed as
## environment variables; 'done' is the file whose existence means it is finished
job <- function(run, script, env = list(), done = NULL, name = script, weight = 1)
  list(run = run, script = script, env = env, done = done, name = name, weight = weight)
COMMON_ENV <- list(WL_RTAU = "1", WL_JAC = "0", WL_THREADS = "1", WL_RCPP_CACHE = "rcpp_cache")
run_one <- function(j) {
  owd <- setwd(j$run); on.exit(setwd(owd))
  if (!is.null(j$done) && all(file.exists(j$done))) return(list(name = j$name, status = "skipped"))
  env <- c(COMMON_ENV, j$env); env <- env[!duplicated(names(env), fromLast = TRUE)]   # a later setting overrides an earlier one
  if (TEST && !is.null(env$WL_MAXIT)) env$WL_MAXIT <- "3"
  env <- lapply(env, function(v) if (is.numeric(v)) format(v, scientific = FALSE, trim = TRUE) else as.character(v))
  do.call(Sys.setenv, env)
  lf <- file.path("logs", paste0(gsub("[^A-Za-z0-9_.-]", "_", j$name), ".log"))
  t0 <- Sys.time()
  st <- system2(RSCRIPT, j$script, stdout = lf, stderr = lf)
  ok <- st == 0 && (is.null(j$done) || all(file.exists(j$done)))
  if (!ok && any(grepl("sourceCppDynlib|zero-byte", readLines(lf, warn = FALSE)))) {   # two scripts read the index of the kernel cache while a third rewrote it
    Sys.sleep(5); st <- system2(RSCRIPT, j$script, stdout = lf, stderr = lf)
    ok <- st == 0 && (is.null(j$done) || all(file.exists(j$done)))
  }
  Sys.unsetenv(names(env))
  list(name = j$name, run = basename(j$run), status = if (ok) "ok" else paste("FAILED, see", file.path(j$run, lf)),
       minutes = round(as.numeric(Sys.time() - t0, units = "mins"), 1))
}
run_jobs <- function(jobs, label) {
  jobs <- Filter(function(j) is.null(j$done) || !all(file.exists(file.path(j$run, j$done))), jobs)
  if (!length(jobs)) { log_msg(label, ": nothing to do"); return(invisible()) }
  jobs <- jobs[order(-sapply(jobs, `[[`, "weight"))]           # longest first
  log_msg(label, ": ", length(jobs), " jobs on ", min(N_CORES, length(jobs)), " cores")
  if (N_CORES > 1 && length(jobs) > 1) {
    cl <- parallel::makeCluster(min(N_CORES, length(jobs))); on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("RSCRIPT", "COMMON_ENV", "TEST"), envir = globalenv())
    res <- parallel::clusterApplyLB(cl, jobs, run_one)
  } else res <- lapply(jobs, run_one)
  for (r in res) log_msg("  ", r$status, "  ", r$run, " ", r$name, if (!is.null(r$minutes)) paste0(" (", r$minutes, " min)") else "")
  bad <- Filter(function(r) grepl("FAILED", r$status), res)
  if (length(bad)) stop(length(bad), " job(s) failed in '", label, "'; fix and rerun run_all.R (finished work is kept)")
}

## ---------------- compile the kernel once per run ----------------------
for (rd in RUNS) {
  r <- run_one(job(rd, "R/00_compile.R", name = "compile"))
  if (grepl("FAILED", r$status)) {
    lg <- readLines(file.path(rd, "logs", "compile.log"))
    cat("\n---- end of", file.path(rd, "logs", "compile.log"), "----\n", tail(lg, 25), sep = "\n")
    if (any(grepl("xcrun|developer path|command line tools", lg, ignore.case = TRUE)))
      cat("\nThe Apple command line tools are missing: run   xcode-select --install   in Terminal, then rerun.\n")
    stop("the C++ kernel did not compile (log above)")
  }
}

SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
SL <- list(c("Men", "snatch"), c("Men", "cleanjerk"), c("Women", "snatch"), c("Women", "cleanjerk"))
iter_of <- function(n, s, second = FALSE, joint = FALSE) if (joint) ITER_JOINT[[n]] else if (second) ITER_SECOND[[n]] else if (paste(s[1], s[2]) %in% LONG[[n]]) ITER_LONG[[n]] else ITER[[n]]

## ---------------- time estimate ----------------------------------------
r <- run_one(job(RUNS[[1]], "R/00_timing.R", done = "out/timing_eval.txt", name = "timing"))
tf <- file.path(RUNS[[1]], "out", "timing_eval.txt")
t_eval <- if (file.exists(tf)) suppressWarnings(as.numeric(readLines(tf))[1]) else NA
if (length(t_eval) && is.finite(t_eval)) {
  ## the women's and the 2017 data are cheaper than the men's snatch on all years
  ## with every core busy an evaluation takes about 1.8 times as long as alone (measured on an 8-core laptop)
  rel <- c(full = 1, since2017 = 0.45)[names(DATASETS)]
  done_gen <- unlist(lapply(names(DATASETS), function(n) sapply(c("men_snatch", "men_cleanjerk", "women_snatch", "women_cleanjerk"), function(sl) {
    if (file.exists(file.path(RUNS[[n]], "out", "chains", paste0("final_", sl, "_c1.rds")))) return(Inf)
    ck <- file.path(RUNS[[n]], "out", "chains", paste0("final_", sl, "_demcz_ckpt.rds"))
    if (file.exists(ck)) tryCatch(readRDS(ck)$g, error = function(e) 0) else 0 })))
  left <- pmax(as.numeric(sapply(names(DATASETS), function(n) sapply(SL, function(s) iter_of(n, s)))) - done_gen, 0)
  w <- rep(rel, each = 4) * rep(c(1, 1, 0.7, 0.7), length(rel)) * 1.8
  for (n in names(DATASETS)) for (sl in TWO_MODES[[n]]) {      # second sets of chains
    f <- tolower(sub(" ", "_", sl))
    if (file.exists(file.path(RUNS[[n]], "out", "chains", paste0("final_", f, "_second_c1.rds")))) next
    ck <- file.path(RUNS[[n]], "out", "chains", paste0("final_", f, "_second_demcz_ckpt.rds"))
    left <- c(left, max(ITER_SECOND[[n]] - (if (file.exists(ck)) tryCatch(readRDS(ck)$g, error = function(e) 0) else 0), 0))
    w <- c(w, rel[[n]] * (if (grepl("Women", sl)) 0.7 else 1) * 1.8)
  }
  job_h <- if (SAMPLER == "demcz") CHAINS * left * t_eval / 3600 * w else rep(2 * left * t_eval / 3600 * w, each = CHAINS)
  log_msg(sprintf("one likelihood evaluation (men's snatch, all years): %.2f s; the %d final-model jobs need %.0f to %.0f more hours each",
                  t_eval, length(job_h), min(job_h), max(job_h)))
  ## the ladder likelihood costs about a twentieth of the final model's; about 30 more hours for the other jobs of the stage
  other_h <- sum(8 * LADDER_CHAINS * LADDER_ITER * t_eval / 20 / 3600 * rel * 1.8) + if (TEST) 0 else 30
  log_msg(sprintf("with %d cores the chain stage takes about %.0f hours; the other stages add about half a day",
                  N_CORES, max(max(job_h), (sum(job_h) + other_h) / N_CORES)))
}

low <- function(s) tolower(s)
J <- function(f) lapply(names(RUNS), f)   # build jobs for every data set
flat <- function(x) unlist(x, recursive = FALSE)

## ---------------- stage 1: descriptives and best-lift fits ---------------
run_jobs(flat(J(function(n) { rd <- RUNS[[n]]; list(
  job(rd, "R/31_descriptives.R", done = "out/desc_jump_rates.csv", name = "descriptives"),
  job(rd, "R/56_ext_descriptives.R", done = "out/desc_opener_prev.csv", name = "openers"),
  job(rd, "R/64_mixture_gpd.R", list(WL_SEX = "Men"), done = "out/mixture_gpd_men.csv", name = "mixture_gpd_men", weight = 3),
  job(rd, "R/64_mixture_gpd.R", list(WL_SEX = "Women"), done = "out/mixture_gpd_women.csv", name = "mixture_gpd_women", weight = 3)) })),
  "stage 1: descriptives and best-lift mixture fits")
run_jobs(flat(J(function(n) { rd <- RUNS[[n]]; list(
  job(rd, "R/66_mix_start.R", done = "out/mixstart_women_cleanjerk.rds", name = "mixture_starts")) })), "stage 1b: starting values")
## modes found in our runs that the optimizer reaches only in some runs (R/mode_starts.rds, see R/82_final_starts.R)
stored_x <- function(n) sprintf("out/start_final_%s_x.rds", sub(" ", "_", low(names(readRDS("R/mode_starts.rds")[[n]]))))
run_jobs(flat(J(function(n) { rd <- RUNS[[n]]; list(
  job(rd, "R/82_final_starts.R", list(WL_DATASET = n), done = c("out/start_final_women_cleanjerk_scl.rds", stored_x(n)), name = "final_starts")) })),
  "stage 1c: starts for the final model")

## ---------------- stage 2: modes of the final model, ladder chains ------
MODE_ENV <- list(WL_SM3 = "1", WL_MIX = "1", WL_MODEL = "e2", WL_G = "24", WL_FIX = "log_tau=-3.91202", WL_MAXIT = "300", WL_RELTOL = "1e-9")
run_jobs(flat(J(function(n) { rd <- RUNS[[n]]; c(
  lapply(SL, function(s) job(rd, "R/52_ext_map.R", c(MODE_ENV, WL_SEX = s[1], WL_LIFT = s[2], WL_LAPLACE = "1",
    WL_START = sprintf("out/start_final_%s_%s_main.rds", low(s[1]), s[2])), done = sprintf("out/extmap_e2_%s_%s%s.rds", low(s[1]), s[2], SUF),
    name = paste("mode", s[1], s[2]), weight = 20)),
  lapply(SL, function(s) job(rd, "R/52_ext_map.R", c(MODE_ENV, WL_SEX = s[1], WL_LIFT = s[2], WL_TAGX = "_alt",
    WL_START = sprintf("out/start_final_%s_%s_alt.rds", low(s[1]), s[2])), done = sprintf("out/extmap_e2_%s_%s%s_alt.rds", low(s[1]), s[2], SUF),
    name = paste("mode_alt", s[1], s[2]), weight = 15)),
  lapply(Filter(function(s) file.exists(file.path(rd, sprintf("out/start_final_%s_%s_x.rds", low(s[1]), s[2]))), SL), function(s)
    job(rd, "R/52_ext_map.R", c(MODE_ENV, WL_SEX = s[1], WL_LIFT = s[2], WL_TAGX = "_x",
    WL_START = sprintf("out/start_final_%s_%s_x.rds", low(s[1]), s[2])), done = sprintf("out/extmap_e2_%s_%s%s_x.rds", low(s[1]), s[2], SUF),
    name = paste("mode_x", s[1], s[2]), weight = 5))) })),
  "stage 2: posterior modes of the final model")
select_modes <- function() flat(J(function(n) list(job(RUNS[[n]], "R/84_final_select.R", list(WL_SECOND = paste(TWO_MODES[[n]], collapse = ",")),
  name = "select modes", weight = 2))))
run_jobs(select_modes(), "stage 2b: choice of the mode and Laplace covariance")

## ---------------- stage 3: chains, anchors, ridge, variants, checks -----
COST <- c(full = 1, since2017 = 0.45)   # relative cost of a likelihood evaluation, for the order of the jobs
run_jobs(flat(J(function(n) { rd <- RUNS[[n]]; main_only <- n == MAIN; c(
  ## final-model chains, longest first
  if (SAMPLER == "demcz") lapply(SL, function(s)       # one job per sex and lift, its chains run together
    job(rd, "R/88_final_demcz.R", list(WL_SEX = s[1], WL_LIFT = s[2], WL_CHAINS = CHAINS, WL_ITER = iter_of(n, s), WL_BURN = BURN[[n]], WL_THIN = THIN,
        WL_COORD = COORD), done = sprintf("out/chains/final_%s_%s_c%d.rds", low(s[1]), s[2], seq_len(CHAINS)),
        name = sprintf("chains %s %s", s[1], s[2]), weight = 1000 * COST[[n]] * (if (s[1] == "Men") 1 else 0.7) * iter_of(n, s) / ITER[[n]])),
  ## second set of chains from the second mode (TWO_MODES)
  if (SAMPLER == "demcz") lapply(Filter(function(s) paste(s[1], s[2]) %in% TWO_MODES[[n]] &&
      file.exists(file.path(rd, "out", sprintf("extmap_e2_%s_%s%s_second.rds", low(s[1]), s[2], SUF))), SL), function(s)
    job(rd, "R/88_final_demcz.R", list(WL_SEX = s[1], WL_LIFT = s[2], WL_CHAINS = CHAINS, WL_ITER = iter_of(n, s, TRUE), WL_BURN = BURN[[n]], WL_THIN = THIN,
        WL_COORD = COORD, WL_FROM = "second"), done = sprintf("out/chains/final_%s_%s_second_c%d.rds", low(s[1]), s[2], seq_len(CHAINS)),
        name = sprintf("chains %s %s second mode", s[1], s[2]), weight = 1000 * COST[[n]] * (if (s[1] == "Men") 1 else 0.7) * iter_of(n, s, TRUE) / ITER[[n]])),
  if (SAMPLER != "demcz") flat(lapply(SL, function(s) lapply(seq_len(CHAINS), function(k)   # one job per chain
    job(rd, "R/83_final_chain.R", list(WL_SEX = s[1], WL_LIFT = s[2], WL_CHAIN = k, WL_ITER = iter_of(n, s), WL_BURN = BURN[[n]], WL_THIN = THIN,
        WL_COORD = COORD), done = sprintf("out/chains/final_%s_%s_c%d.rds", low(s[1]), s[2], k),
        name = sprintf("chain %s %s c%d", s[1], s[2], k), weight = 1000 * COST[[n]] * (if (s[1] == "Men") 1 else 0.7) * iter_of(n, s) / ITER[[n]])))),
  ## record-book rungs of the ladder (4 chains each, run together)
  flat(lapply(SL, function(s) lapply(c("L1", "L2"), function(r)
    job(rd, "R/85_ladder_chain.R", list(WL_SEX = s[1], WL_LIFT = s[2], WL_RUNG = r, WL_CHAINS = LADDER_CHAINS, WL_ITER = LADDER_ITER,
      WL_BURN = LADDER_BURN, WL_THIN = LADDER_THIN, WL_COORD = "sstar"),
      done = sprintf("out/chains/ladder_%s_%s_%s_c%d.rds", r, low(s[1]), s[2], seq_len(LADDER_CHAINS)),
      name = sprintf("ladder %s %s %s", r, s[1], s[2]), weight = 10)))),
  ## the search for further modes from the third start (82_final_starts.R); 84_final_select.R
  ## takes them into account when it runs again before the weights of the modes
  lapply(SL, function(s) job(rd, "R/52_ext_map.R", c(MODE_ENV, WL_SEX = s[1], WL_LIFT = s[2], WL_TAGX = "_scl",
    WL_START = sprintf("out/start_final_%s_%s_scl.rds", low(s[1]), s[2])), done = sprintf("out/extmap_e2_%s_%s%s_scl.rds", low(s[1]), s[2], SUF),
    name = paste("mode_scl", s[1], s[2]), weight = 15)),
  ## the two other anchors, from the mode at tau = 0.02
  flat(lapply(SL, function(s) lapply(c("-4.60517", "-3.50656"), function(lt)
    job(rd, "R/52_ext_map.R", c(MODE_ENV[names(MODE_ENV) != "WL_FIX"], WL_SEX = s[1], WL_LIFT = s[2], WL_FIX = paste0("log_tau=", lt), WL_MAXIT = "150",
      WL_START = sprintf("out/extmap_e2_%s_%s%s.rds", low(s[1]), s[2], SUF)), done = sprintf("out/extmap_e2_%s_%s_rtau_sm3_mix_fix_log_tau_%s.rds", low(s[1]), s[2], lt),
      name = sprintf("anchor %s %s %s", s[1], s[2], lt), weight = 15)))),
  ## the final model with tau free (men's snatch): it drifts along the direction of Proposition 1
  list(job(rd, "R/52_ext_map.R", c(MODE_ENV[names(MODE_ENV) != "WL_FIX"], WL_SEX = "Men", WL_LIFT = "snatch", WL_MAXIT = "600",
      WL_START = sprintf("out/extmap_e2_men_snatch%s.rds", SUF)), done = "out/extmap_e2_men_snatch_rtau_sm3_mix.rds", name = "tau free", weight = 60)),
  ## predictive checks at the mode
  lapply(SL, function(s) job(rd, "R/53_ext_checks.R", list(WL_SEX = s[1], WL_LIFT = s[2], WL_G = "24",
      WL_FITS = sprintf("final=out/extmap_e2_%s_%s%s.rds", low(s[1]), s[2], SUF)), done = sprintf("out/ppc_%s_%s.done", low(s[1]), s[2]),
      name = paste("ppc", s[1], s[2]), weight = 5)),
  ## the ridge (Figure 2) and the attempt-model variants: the main data set only
  if (main_only) c(
    list(job(rd, "R/89_ridge.R", done = c(sprintf("out/extmap_base_men_snatch_rtau_mix_fix_log_tau_%s.rds",
        c("-4.60517", "-4.19971", "-3.91202", "-3.50656", "-2.99573", "-2.52573")), "out/extmap_base_men_snatch_rtau_mix.rds"),
        name = "ridge", weight = 150),
      job(rd, "R/52_ext_map.R", list(WL_MODEL = "base", WL_SM3 = "1", WL_MIX = "1", WL_G = "24", WL_SEX = "Men", WL_LIFT = "snatch",
        WL_FIX = "log_tau=-3.91202", WL_MAXIT = "300", WL_RELTOL = "1e-9", WL_START = sprintf("out/extmap_e2_men_snatch%s.rds", SUF)),
        done = "out/extmap_base_men_snatch_rtau_sm3_mix_fix_log_tau_-3.91202.rds", name = "variant reference", weight = 12),
      job(rd, "R/52_ext_map.R", list(WL_MODEL = "base", WL_SM3 = "1", WL_MIX = "1", WL_FLOOR = "1", WL_G = "24", WL_SEX = "Men", WL_LIFT = "snatch",
        WL_FIX = "log_tau=-3.91202", WL_MAXIT = "300", WL_RELTOL = "1e-9", WL_START = sprintf("out/extmap_e2_men_snatch%s.rds", SUF)),
        done = "out/extmap_base_men_snatch_rtau_sm3_floor_mix_fix_log_tau_-3.91202.rds", name = "variant floor", weight = 15)),
    if (RUN_E1) list(job(rd, "R/52_ext_map.R", list(WL_MODEL = "e1", WL_SM3 = "1", WL_MIX = "1", WL_G = "24", WL_SEX = "Men", WL_LIFT = "snatch",
        WL_FIX = "log_tau=-3.91202", WL_MAXIT = "150", WL_RELTOL = "1e-8", WL_START = sprintf("out/extmap_e2_men_snatch%s.rds", SUF)),
        done = "out/extmap_e1_men_snatch_rtau_sm3_mix_fix_log_tau_-3.91202.rds", name = "variant shared error", weight = 200)),
    if (RUN_E12) lapply(SL, function(s) job(rd, "R/52_ext_map.R", list(WL_MODEL = "e12", WL_SM3 = "1", WL_MIX = "1", WL_G = "24", WL_KETA = "4",
        WL_SEX = s[1], WL_LIFT = s[2], WL_FIX = "log_tau=-3.91202", WL_MAXIT = "200", WL_RELTOL = "1e-8",
        WL_START = sprintf("out/extmap_e2_%s_%s%s.rds", low(s[1]), s[2], SUF)),
        done = sprintf("out/extmap_e12_%s_%s%s.rds", low(s[1]), s[2], SUF), name = paste("shared error", s[1], s[2]), weight = 300)))) })),
  "stage 3: chains of the final model and of the ladder, anchors, ridge, variants, predictive checks")

## ---------------- stage 4: pooling, summaries, copula, checks ---------------
run_jobs(select_modes(), "stage 4: the modes again, with those of the third start")
regions <- function(done = "out/basin_weights.csv") flat(J(function(n) list(job(RUNS[[n]], "R/96_basin_weights.R", done = done, name = "regions of the posterior"))))
run_jobs(regions(), "stage 4: regions of the posterior")
## the joint sets, for the combinations whose chains sit in two regions
joint_jobs <- flat(J(function(n) { rd <- RUNS[[n]]
  bw <- read.csv(file.path(rd, "out/basin_weights.csv")); need <- unique(bw[bw$two_regions, c("sex", "lift")])
  lapply(seq_len(nrow(need)), function(i) job(rd, "R/88_final_demcz.R", list(WL_SEX = need$sex[i], WL_LIFT = need$lift[i], WL_FROM = "joint",
    WL_ITER = iter_of(n, NULL, joint = TRUE), WL_BURN = round(iter_of(n, NULL, joint = TRUE) / 10), WL_THIN = THIN, WL_PJUMP = "0.5", WL_COORD = COORD),
    done = sprintf("out/chains/final_%s_%s_joint_c1.rds", low(need$sex[i]), need$lift[i]), name = paste("joint chains", need$sex[i], need$lift[i]),
    weight = 1000 * COST[[n]])) }))
joint_jobs <- lapply(joint_jobs, function(j) { j$env$WL_PAR <- max(2, N_CORES %/% max(1, length(joint_jobs))); j })
## the shared error fitted from the second mode as well, where there is one (Section S5), so that the
## check does not depend on which mode the fit started from; it is the longest job of this stage
e12_second <- if (RUN_E12) lapply(Filter(function(s) file.exists(file.path(RUNS[[MAIN]], sprintf("out/extmap_e2_%s_%s%s_second.rds", low(s[1]), s[2], SUF))),
  lapply(TWO_MODES[[MAIN]], function(x) strsplit(x, " ")[[1]])), function(s)
  job(RUNS[[MAIN]], "R/52_ext_map.R", list(WL_MODEL = "e12", WL_SM3 = "1", WL_MIX = "1", WL_G = "24", WL_KETA = "4", WL_SEX = s[1], WL_LIFT = s[2],
    WL_FIX = "log_tau=-3.91202", WL_MAXIT = "200", WL_RELTOL = "1e-8", WL_TAGX = "_second",
    WL_START = sprintf("out/extmap_e2_%s_%s%s_second.rds", low(s[1]), s[2], SUF)),
    done = sprintf("out/extmap_e12_%s_%s%s_second.rds", low(s[1]), s[2], SUF), name = paste("shared error", s[1], s[2], "second mode"), weight = 300)) else list()
## the two other anchors from the second mode as well (Table S9)
anchors_second <- lapply(Filter(function(s) file.exists(file.path(RUNS[[MAIN]], sprintf("out/extmap_e2_%s_%s%s_second.rds", low(s[1]), s[2], SUF))),
  lapply(TWO_MODES[[MAIN]], function(x) strsplit(x, " ")[[1]])), function(s) lapply(c("-4.60517", "-3.50656"), function(lt)
  job(RUNS[[MAIN]], "R/52_ext_map.R", c(MODE_ENV[names(MODE_ENV) != "WL_FIX"], WL_SEX = s[1], WL_LIFT = s[2], WL_FIX = paste0("log_tau=", lt), WL_MAXIT = "150",
    WL_TAGX = "_second", WL_START = sprintf("out/extmap_e2_%s_%s%s_second.rds", low(s[1]), s[2], SUF)),
    done = sprintf("out/extmap_e2_%s_%s_rtau_sm3_mix_fix_log_tau_%s_second.rds", low(s[1]), s[2], lt),
    name = sprintf("anchor %s %s %s second mode", s[1], s[2], lt), weight = 15)))
## the predictive checks at the second mode (model label final_second in out/ppc_ext_*.csv)
ppc_second <- lapply(Filter(function(s) file.exists(file.path(RUNS[[MAIN]], sprintf("out/extmap_e2_%s_%s%s_second.rds", low(s[1]), s[2], SUF))),
  lapply(TWO_MODES[[MAIN]], function(x) strsplit(x, " ")[[1]])), function(s)
  job(RUNS[[MAIN]], "R/53_ext_checks.R", list(WL_SEX = s[1], WL_LIFT = s[2], WL_G = "24", WL_PPC_TAG = "_second",
    WL_FITS = sprintf("final_second=out/extmap_e2_%s_%s%s_second.rds", low(s[1]), s[2], SUF)), done = sprintf("out/ppc_%s_%s_second.done", low(s[1]), s[2]),
    name = paste("ppc", s[1], s[2], "second mode"), weight = 5))
run_jobs(c(joint_jobs, e12_second, flat(anchors_second), ppc_second), "stage 4j: joint chains between the two regions; anchors, shared error and predictive checks at the second mode")
## each anchor of a posterior with two modes refitted from the margins of the other mode: the margin
## parameters lie on a flat ridge at tau = 0.03, and 69_anchor_mix.R keeps the higher of the two fits
anchors_swap <- lapply(Filter(function(s) file.exists(file.path(RUNS[[MAIN]], sprintf("out/extmap_e2_%s_%s%s_second.rds", low(s[1]), s[2], SUF))),
  lapply(TWO_MODES[[MAIN]], function(x) strsplit(x, " ")[[1]])), function(s) flat(lapply(c("-4.60517", "-3.50656"), function(lt) {
  fk <- sprintf("out/extmap_e2_%s_%s_rtau_sm3_mix_fix_log_tau_%s", low(s[1]), s[2], lt)
  lapply(list(c("", "_second"), c("_second", "")), function(m)
    job(RUNS[[MAIN]], "R/52_ext_map.R", c(MODE_ENV[names(MODE_ENV) != "WL_FIX"], WL_SEX = s[1], WL_LIFT = s[2], WL_FIX = paste0("log_tau=", lt), WL_MAXIT = "300", WL_RELTOL = "1e-10",
      WL_TAGX = paste0(m[1], "_swap"), WL_START = paste0(fk, m[1], ".rds"), WL_MARGINS_FROM = paste0(fk, m[2], ".rds")),
      done = paste0(fk, m[1], "_swap.rds"), name = sprintf("anchor %s %s %s%s from the other mode's margins", s[1], s[2], lt, sub("_", ", ", m[1])), weight = 15)) })))
run_jobs(flat(anchors_swap), "stage 4k: anchors of the two modes refitted from each other's margins")
run_jobs(regions(done = NULL), "stage 4: weights of the regions")
run_jobs(flat(J(function(n) list(job(RUNS[[n]], "R/86_combine_chains.R", done = "out/diag_summary.csv", name = "combine chains"),
                                  job(RUNS[[n]], "R/87_merge_ppc.R", done = "out/ppc_ext_jumps.csv", name = "merge predictive checks")))),
  "stage 4a: pooling the chains")
run_jobs(flat(J(function(n) list(job(RUNS[[n]], "R/58_ext_summaries.R", list(WL_NDRAW = "4000"), done = "out/postx_ceiling.csv", name = "posterior summaries"),
                                  job(RUNS[[n]], "R/92_record_check.R", list(WL_NDRAW = "4000"), done = "out/record_check.csv", name = "records against the ceiling")))),
  "stage 4b: posterior summaries")
run_jobs(flat(J(function(n) { rd <- RUNS[[n]]; c(list(
  job(rd, "R/57_ext_copula.R", list(WL_NMARG = COPULA_DRAWS, WL_NLAT = "40",
      WL_FS_MEN = sprintf("out/mcmcx_e2_men_snatch%s.rds", SUF), WL_FC_MEN = sprintf("out/mcmcx_e2_men_cleanjerk%s.rds", SUF),
      WL_FS_WOMEN = sprintf("out/mcmcx_e2_women_snatch%s.rds", SUF), WL_FC_WOMEN = sprintf("out/mcmcx_e2_women_cleanjerk%s.rds", SUF)),
      done = "out/totalx_returnlevel.csv", name = "copula", weight = 30),
  job(rd, "R/63_bodyfree_pot.R", list(WL_BOOT = BOOT), done = "out/bodyfree_pot.csv", name = "body-free check", weight = 10),
  job(rd, "R/62_bestlift_evt.R", list(WL_Q = "0.80"), done = "out/bestlift_evt_q80.csv", name = "fixed splice 0.80", weight = 3),
  job(rd, "R/62_bestlift_evt.R", list(WL_Q = "0.90"), done = "out/bestlift_evt.csv", name = "fixed splice 0.90", weight = 3),
  job(rd, "R/62_bestlift_evt.R", list(WL_Q = "0.95"), done = "out/bestlift_evt_q95.csv", name = "fixed splice 0.95", weight = 3),
  job(rd, "R/70_scale_trend_check.R", done = "out/scale_trend_check.csv", name = "time trend in the scale", weight = 2),
  job(rd, "R/65_mix_check.R", done = NULL, name = "kernel checks", weight = 2)),
  if (n == names(RUNS)[1]) list(job(rd, "R/72_recent_era_check.R", done = "out/recent_era_check.csv", name = "recent era check", weight = 2))) })),
  "stage 4c: copula, checks and comparisons")
run_jobs(flat(J(function(n) list(job(RUNS[[n]], "R/68_ladder_summary.R", list(WL_NDRAW = "4000"), done = "out/ladder_summary.csv", name = "ladder summary")))),
  "stage 4d: ladder table")
run_jobs(flat(J(function(n) list(job(RUNS[[n]], "R/59_ext_figures.R", name = "figures")))), "stage 5: figures")

## ---------------- stage 6: simulation study, record book, priors, quadrature (main data set) ----
rd <- RUNS[[MAIN]]
sim_jobs <- flat(lapply(list(c("Men", "snatch"), c("Women", "cleanjerk")), function(s) c(
  lapply(seq_len(NREP_SIM), function(r) job(rd, "R/94_simulation.R", list(WL_SEX = s[1], WL_LIFT = s[2], WL_REP = r, WL_FITTAU = "0.02", WL_MAXIT = "300"),
    done = sprintf("out/sim/fit_%s_%s_r%02d_tau0.020.rds", low(s[1]), s[2], r), name = sprintf("simulation %s %s r%02d tau 0.02", s[1], s[2], r), weight = 20)),
  flat(lapply(c("0.01", "0.03"), function(tt) lapply(seq_len(NREP_ANCHOR), function(r)
    job(rd, "R/94_simulation.R", list(WL_SEX = s[1], WL_LIFT = s[2], WL_REP = r, WL_FITTAU = tt, WL_MAXIT = "300"),
      done = sprintf("out/sim/fit_%s_%s_r%02d_tau%.3f.rds", low(s[1]), s[2], r, as.numeric(tt)),
      name = sprintf("simulation %s %s r%02d tau %s", s[1], s[2], r, tt), weight = 20)))))))
run_jobs(c(sim_jobs, list(
  job(rd, "R/95_ppc_record.R", list(WL_NDRAW = NDRAW_PPC), done = "out/ppc_record_summary.csv", name = "predictive check of the record book", weight = 10),
  job(rd, "R/98_prior_reweight.R", done = "out/prior_reweight.csv", name = "prior reweighting", weight = 5),
  job(rd, "R/91_quadrature_check.R", list(WL_NUNIT = if (TEST) "2" else "25"), done = "out/quadrature_check.csv", name = "quadrature check", weight = 10))),
  "stage 6: simulation study, record-book check, prior sensitivity, quadrature")
run_jobs(c(list(job(rd, "R/97_review_summary.R", done = "out/sim_summary.csv", name = "simulation summary"),
  job(rd, "R/87_merge_ppc.R", name = "merge predictive checks")),
  flat(J(function(n) list(job(RUNS[[n]], "R/69_anchor_mix.R", done = "out/anchor_mix.csv", name = "anchors summary", weight = 3)))),
  if (RUN_E12) list(job(rd, "R/90_shared_error.R", done = "out/shared_error.csv", name = "shared error summary", weight = 5))),
  "stage 6b: simulation summary, predictive checks merged, anchors summary, shared error summary")

## ---------------- stage 7: the tables of the paper ------------------------
Sys.setenv(WL_RUN_MAIN = RUNS[[MAIN]], WL_RUN_FULL = RUNS[["full"]], WL_TABLES = file.path(PROJECT, if (TEST) "tables_test" else "tables"))
st <- system2(RSCRIPT, "R/99_tables.R", stdout = file.path(PROJECT, "runs", "tables.log"), stderr = file.path(PROJECT, "runs", "tables.log"))
if (st != 0) stop("the tables failed, see runs/tables.log")
## ---------------- one archive of the results ----------------------------------
owd <- setwd(PROJECT)
tb <- if (TEST) "tables_test" else "tables"
res <- c(list.files(tb, full.names = TRUE), file.path("runs", c("run_all.log", "tables.log")),
         unlist(lapply(RUNS, function(rd) { r <- file.path("runs", basename(rd))
           c(list.files(file.path(r, "out"), pattern = "\\.(csv|txt|tex)$", full.names = TRUE),
             list.files(file.path(r, "out"), pattern = "^(mcmcx|extmap)_e2_.*\\.rds$", full.names = TRUE),
             list.files(file.path(r, "out", "sim"), pattern = "\\.rds$", full.names = TRUE),
             list.files(file.path(r, "figs"), full.names = TRUE), list.files(file.path(r, "logs"), full.names = TRUE)) })))
res <- res[file.exists(res)]
zf <- if (TEST) "results_test.zip" else "results.zip"
if (file.exists(zf)) file.remove(zf)
utils::zip(zf, res, flags = "-q")
setwd(owd)
log_msg("all done. Tables in ", tb, "/, figures in ", file.path("runs", basename(RUNS[[MAIN]]), "figs"),
        ", everything else in runs/*/out (see README.md). The results are collected in ", zf)
