## 94_simulation.R
## -----------------------------------------------------------------------
## Simulation study: one replicate and one fit per job. Data are simulated from
## the posterior mode of the final model at tau = 0.02 (93_sim_lib.R), on the
## design of the analysis set (bodyweights, years, athletes, number of
## attempts), and the final model is refitted by its posterior mode with the
## attempt scale fixed at WL_FITTAU, started with the tail shape at zero and the
## threshold at 1.5. A fit at 0.02 checks recovery; fits at
## 0.01 and 0.03 show what a wrong anchor does. The three fits of a replicate
## see the same simulated data.
## Settings: WL_SEX, WL_LIFT, WL_REP (replicate number), WL_FITTAU, WL_MAXIT.
## Output: out/sim/fit_<sex>_<lift>_r<rep>_tau<fittau>.rds
source("R/50_ext_lib.R"); source("R/93_sim_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
sex <- Sys.getenv("WL_SEX", "Men"); lift <- Sys.getenv("WL_LIFT", "snatch")
REP <- as.integer(Sys.getenv("WL_REP", "1")); FIT_TAU <- as.numeric(Sys.getenv("WL_FITTAU", "0.02"))
MAXIT <- as.integer(Sys.getenv("WL_MAXIT", "300"))
SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
dir.create("out/sim", showWarnings = FALSE, recursive = TRUE)
tag <- sprintf("fit_%s_%s_r%02d_tau%.3f", tolower(sex), lift, REP, FIT_TAU)
outf <- file.path("out/sim", paste0(tag, ".rds"))
if (file.exists(outf)) { cat("done already\n"); quit(save = "no") }

d <- build_data_ext(build_data(att, sex, lift, G = 24), jac = 0)
truth <- fill_par2(readRDS(sprintf("out/extmap_e2_%s_%s%s.rds", tolower(sex), lift, SUF))$par)
truth[["log_tau"]] <- log(0.02)
seed <- 1000 * match(paste(sex, lift), c("Men snatch", "Men cleanjerk", "Women snatch", "Women cleanjerk")) + REP
sim <- sim_attempts(d, truth, repeat_rate(d), seed)
ds <- with_attempts(d, sim$W, sim$S)
cat(sprintf("== %s %s, replicate %d (seed %d), fit at tau = %.3f: %d athlete-day-lifts, success rate %.3f\n",
            sex, lift, REP, seed, FIT_TAU, ds$n, mean(ds$S == 1, na.rm = TRUE)))

## the fit starts from the true attempt, declaration, body and regression parameters, with the
## tail neutral: shape zero and threshold at the centre of its prior
start <- truth; start[["log_tau"]] <- log(FIT_TAU); start[["xi"]] <- 0; start[["mx_u"]] <- 1.5
fit <- fit_mode(ds, start, fixed = c(log_tau = log(FIT_TAU)), maxit = MAXIT,
                ckpt = file.path("out/sim", paste0(tag, "_ckpt.rds")))
cat(sprintf("fit: logpost %.2f, convergence %d, %d evaluations, %.1f min\n", fit$value, fit$conv, fit$evals, fit$minutes))
q_true <- ref_quantities(truth, d$ref, sex); q_fit <- ref_quantities(fit$par, d$ref, sex)
best <- best_made(ds$W, ds$S)
saveRDS(list(sex = sex, lift = lift, rep = REP, seed = seed, fit_tau = FIT_TAU, truth = truth, par = fit$par,
             value = fit$value, loglik = fit$loglik, conv = fit$conv, evals = fit$evals, minutes = fit$minutes,
             q_true = q_true, q_fit = q_fit, success_rate = mean(ds$S == 1, na.rm = TRUE),
             gap_true = mean(sim$Z - best, na.rm = TRUE)), outf)
unlink(file.path("out/sim", paste0(tag, "_ckpt.rds")))
print(rbind(truth = q_true, fit = q_fit))
