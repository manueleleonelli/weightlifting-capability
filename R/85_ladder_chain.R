## 85_ladder_chain.R
## -----------------------------------------------------------------------
## A record-book rung of the ladder (L1 or L2) for one sex and lift, with the
## population sampler run_demcz of 80_sampler.R: WL_CHAINS chains run
## together and share an archive of past states. The initial archive is
## drawn around the posterior mode with four times the Laplace variance and
## with the threshold spread uniformly between 0.6 and 2.6 scale units, so
## that a posterior with more than one mode in the threshold is explored
## from the start and shows up in R-hat if the chains do not agree.
## Settings: WL_SEX, WL_LIFT, WL_RUNG (L1|L2), WL_CHAINS, WL_ITER (generations),
##   WL_BURN, WL_THIN.
## Output: out/chains/ladder_<rung>_<sex>_<lift>_c<k>.rds, one per chain
source("R/50_ext_lib.R"); source("R/80_sampler.R"); source("R/81_ladder_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
sex <- Sys.getenv("WL_SEX", "Men"); lift <- Sys.getenv("WL_LIFT", "snatch"); rung <- Sys.getenv("WL_RUNG", "L1")
NCH <- as.integer(Sys.getenv("WL_CHAINS", "4"))
N_GEN <- as.integer(Sys.getenv("WL_ITER", "200000")); N_BURN <- as.integer(Sys.getenv("WL_BURN", "50000")); THIN <- as.integer(Sys.getenv("WL_THIN", "20"))
dir.create("out/chains", showWarnings = FALSE)
tag <- sprintf("ladder_%s_%s_%s", rung, tolower(sex), lift)
outf <- sprintf("out/chains/%s_c%d.rds", tag, seq_len(NCH))
if (all(file.exists(outf))) { cat("done already\n"); quit(save = "no") }
L <- ladder_data(att, sex, lift); athlete <- rung == "L2"
logpost <- function(p) ladder_logpost(p, L, athlete)

## posterior mode from the best-lift mixture fit, on the sampler coordinates
ms <- readRDS(sprintf("out/mixstart_%s_%s.rds", tolower(sex), lift))$par
st <- sanitize_mix(ms[LADDER_NAMES]); if (athlete) st[["rhoa"]] <- 0.85
th <- par_to_theta(st)
f <- function(t) { names(t) <- names(th); p <- theta_to_par(t); if (is.null(p)) return(1e300); v <- logpost(p); if (is.finite(v)) -(v + log_jacobian(p)) else 1e300 }
s0 <- exp(st[["g0"]])
ps <- c(b0 = 1, b1 = 1, b2 = 0.5, b3 = 0.05, g0 = 0.02, g1 = 0.02, g2 = 0.02, xi = 0.02, mx_m2 = 0.1, mx_m3 = 0.1, mx_ls2 = 0.05, mx_ls3 = 0.05,
        mx_a2 = 0.1, mx_a3 = 0.1, mx_u = 0.05, mx_sstar = 0.02, rhoa = 0.01,
        kg_m2 = 0.1 * s0, kg_m3 = 0.1 * s0, kg_ls2 = 0.05, kg_ls3 = 0.05, kg_u = 0.05 * s0, kg_ss = 0.02 * s0)[names(th)]
o <- optim(th, f, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12, parscale = ps))
o <- optim(o$par, f, method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-12, parscale = ps))
o <- optim(o$par, f, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12, parscale = ps))
th0 <- o$par; names(th0) <- names(th); k <- length(th0)
H <- tryCatch(optimHess(th0, f), error = function(e) NULL)
V0 <- if (!is.null(H) && all(is.finite(H))) tryCatch(solve(H), error = function(e) NULL) else NULL
if (is.null(V0) || any(!is.finite(V0)) || any(diag(V0) <= 0)) V0 <- diag(ps^2)
V0 <- (V0 + t(V0)) / 2; dimnames(V0) <- list(names(th0), names(th0))
cat(sprintf("%s: mode u %.3f xi %.3f logpost %.2f\n", tag, theta_to_par(th0)[["mx_u"]], th0[["xi"]], -o$value))

## initial archive: around the mode, with the threshold spread over (0.6, 2.6)
set.seed(31 + match(rung, c("L1", "L2"))); Lc <- chol_safe(4 * V0)
Z0 <- matrix(NA_real_, 0, k, dimnames = list(NULL, names(th0))); tries <- 0
while (nrow(Z0) < 10 * k && tries < 200 * k) {
  tries <- tries + 1
  p <- theta_to_par(th0 + as.numeric(rnorm(k) %*% Lc)); if (is.null(p)) next
  p[["mx_u"]] <- runif(1, 0.6, 2.6); t1 <- par_to_theta(p)
  if (f(t1) < 1e299) Z0 <- rbind(Z0, t1)
}
if (nrow(Z0) < max(NCH, 3)) stop("could not build the initial archive")
X0 <- Z0[sample.int(nrow(Z0), NCH), , drop = FALSE]
res <- run_demcz(logpost, Z0, X0, N_GEN, N_BURN, THIN, ckpt = sprintf("out/chains/%s_ckpt.rds", tag), seed = 13,
                 ckpt_every = 2000, verbose_every = 20000, label = tag)
for (i in seq_len(NCH))
  saveRDS(list(draws = res$chains[[i]], lp = res$lp[[i]], acc = res$acc[i], n_iter = res$n_gen, n_burn = res$n_burn, thin = res$thin,
               hours = res$hours, ref = L$ref, sex = sex, lift = lift, rung = rung, chain = i, n = L$n,
               mode = theta_to_par(th0), mode_logpost = -o$value, sampler = "DE-MCzs"), outf[i])
cat("written", outf, "\n")   # the checkpoint is kept: deleting the outputs and rerunning with more generations extends the run
