## 83_final_chain.R
## -----------------------------------------------------------------------
## One chain of the final model (tau = 0.02) for one sex and lift, with the
## sampler of 80_sampler.R. Chain 1 starts at the posterior mode; the others
## start from overdispersed draws around it (twice the Laplace standard
## deviations, and the threshold moved by up to half a scale unit), so that
## the convergence diagnostics compare chains that began apart.
## Settings: WL_SEX, WL_LIFT, WL_CHAIN, WL_ITER, WL_BURN, WL_THIN, WL_THREADS.
## Output: out/chains/final_<sex>_<lift>_c<chain>.rds (checkpoint *_ckpt.rds
## while running; rerunning the script resumes from it)
source("R/50_ext_lib.R"); source("R/80_sampler.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
sex <- Sys.getenv("WL_SEX", "Men"); lift <- Sys.getenv("WL_LIFT", "snatch"); chain <- as.integer(Sys.getenv("WL_CHAIN", "1"))
N_ITER <- as.integer(Sys.getenv("WL_ITER", "60000")); N_BURN <- as.integer(Sys.getenv("WL_BURN", "20000")); THIN <- as.integer(Sys.getenv("WL_THIN", "10"))
SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
dir.create("out/chains", showWarnings = FALSE)
tag <- sprintf("final_%s_%s_c%d", tolower(sex), lift, chain); outf <- sprintf("out/chains/%s.rds", tag)
if (file.exists(outf)) { cat("done already:", outf, "\n"); quit(save = "no") }
d <- build_data_ext(build_data(att, sex, lift, G = 24), jac = 0)
mode <- readRDS(sprintf("out/extmap_e2_%s_%s%s.rds", tolower(sex), lift, SUF))
full <- fill_par2(mode$par); full[["log_tau"]] <- log(0.02)
free <- setdiff(names(full)[is.finite(full)], "log_tau")
logpost <- function(p) { x <- full; x[names(p)] <- p; log_post_ext(x, d) }

## proposal covariance on the sampler coordinates from the Laplace covariance
V <- mode$V[free, free]
nt <- names(par_to_theta(full[free]))
J <- numDeriv::jacobian(function(v) { names(v) <- free; as.numeric(par_to_theta(v)) }, full[free])
Vt <- J %*% V %*% t(J); dimnames(Vt) <- list(nt, nt); Vt <- (Vt + t(Vt)) / 2
## the curvature at the mode understates the spread of the threshold: widen it
tu <- if ("kg_u" %in% nt) "kg_u" else "mx_u"; su <- if (tu == "kg_u") exp(full[["g0"]]) else 1
Vt[tu, tu] <- max(Vt[tu, tu], (0.05 * su)^2)

## start
set.seed(1000 + chain)
start <- full[free]
if (chain > 1) {
  L <- chol_safe(Vt)
  for (tries in 1:200) {
    th <- par_to_theta(full[free]); th <- th + 2 * as.numeric(rnorm(length(th)) %*% L)
    th[[tu]] <- th[[tu]] + su * runif(1, -0.5, 0.5)
    p <- theta_to_par(th); if (is.null(p)) next
    if (is.finite(logpost(p))) { start <- p; break }
  }
}
cat(sprintf("%s: start log posterior %.1f (mode %.1f)\n", tag, logpost(start), mode$value))
blocks <- list(all = nt, tail = tail_block(nt))
res <- run_amwg(logpost, start, Vt, blocks, N_ITER, N_BURN, THIN, ckpt = sprintf("out/chains/%s_ckpt.rds", tag),
                seed = chain, ckpt_every = 500, verbose_every = 1000, label = tag)
## full parameter vectors (fixed values filled in) for the summary scripts
D <- matrix(full, nrow(res$draws), length(full), byrow = TRUE, dimnames = list(NULL, names(full)))
D[, colnames(res$draws)] <- res$draws
res$draws <- D; res$ref <- d$ref; res$sex <- sex; res$lift <- lift; res$chain <- chain; res$start <- start
saveRDS(res, outf); unlink(sprintf("out/chains/%s_ckpt.rds", tag))
cat("written", outf, "\n")
