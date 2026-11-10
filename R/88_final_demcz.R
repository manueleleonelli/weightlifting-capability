## 88_final_demcz.R
## -----------------------------------------------------------------------
## The final model (tau = 0.02) for one sex and lift with the population
## sampler run_demcz of 80_sampler.R: WL_CHAINS chains run together in one
## process and share an archive of past states. The initial archive is drawn
## around the posterior mode from the Laplace approximation with four times
## its variance and the threshold spread over half a scale unit either side;
## the chains start from members of it.
## Settings: WL_SEX, WL_LIFT, WL_CHAINS, WL_ITER (generations), WL_BURN, WL_THIN;
## WL_FROM=second starts a second set of chains around the second mode of
## 84_final_select.R (..._second.rds) instead of the kept mode.
## WL_FROM=joint starts one chain from the last state of every chain of the
## first two sets, with an archive made of draws of both sets, so that the
## proposals carry the differences between the two regions of the posterior
## and the chains move between them (Section S6). WL_PAR=n evaluates the
## chains' proposals on n worker processes (parallel::makeCluster).
## Output: out/chains/final_<sex>_<lift>_c<k>.rds (final_<sex>_<lift>_second_c<k>.rds
## for the second set, final_<sex>_<lift>_joint_c<k>.rds for the joint set), one
## file per chain, in the format of 83_final_chain.R
## (checkpoint out/chains/final_<sex>_<lift>[_second|_joint]_demcz_ckpt.rds)
source("R/50_ext_lib.R"); source("R/80_sampler.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
sex <- Sys.getenv("WL_SEX", "Men"); lift <- Sys.getenv("WL_LIFT", "snatch")
NCH <- as.integer(Sys.getenv("WL_CHAINS", "4"))
N_GEN <- as.integer(Sys.getenv("WL_ITER", "60000")); N_BURN <- as.integer(Sys.getenv("WL_BURN", "20000")); THIN <- as.integer(Sys.getenv("WL_THIN", "10"))
SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
dir.create("out/chains", showWarnings = FALSE)
FROM <- Sys.getenv("WL_FROM", "")   # "", "second" or "joint"
PAR <- as.integer(Sys.getenv("WL_PAR", "0"))
sfx <- switch(FROM, second = "_second", joint = "_joint", "")
tag <- sprintf("final_%s_%s%s", tolower(sex), lift, sfx)
prev <- if (FROM == "joint") c(Sys.glob(sprintf("out/chains/final_%s_%s_c[0-9]*.rds", tolower(sex), lift)),
                               Sys.glob(sprintf("out/chains/final_%s_%s_second_c[0-9]*.rds", tolower(sex), lift))) else character(0)
if (FROM == "joint") NCH <- length(prev)
outf <- sprintf("out/chains/%s_c%d.rds", tag, seq_len(NCH))
if (all(file.exists(outf))) { cat("done already\n"); quit(save = "no") }
d <- build_data_ext(build_data(att, sex, lift, G = 24), jac = 0)
mode <- readRDS(sprintf("out/extmap_e2_%s_%s%s%s.rds", tolower(sex), lift, SUF, if (FROM == "second") "_second" else ""))
full <- fill_par2(mode$par); full[["log_tau"]] <- log(0.02)
free <- setdiff(names(full)[is.finite(full)], "log_tau")
logpost <- function(p) { x <- full; x[names(p)] <- p; log_post_ext(x, d) }
## worker processes for the proposals of a generation
cl <- NULL
if (PAR > 1) {
  cl <- parallel::makeCluster(PAR)
  parallel::clusterExport(cl, c("sex", "lift", "full"))
  for (w in seq_along(cl)) invisible(parallel::clusterEvalQ(cl[w], {   # one worker at a time: they share the compiled kernel's cache
    source("R/50_ext_lib.R"); source("R/80_sampler.R")
    att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE, col_types = cols(category = col_character(), .default = col_guess()))
    d <- build_data_ext(build_data(att, sex, lift, G = 24), jac = 0)
    logpost <- function(p) { x <- full; x[names(p)] <- p; log_post_ext(x, d) }
    target <- function(th) { p <- theta_to_par(th); if (is.null(p)) return(-Inf); v <- logpost(p); if (!is.finite(v)) return(-Inf); v + log_jacobian(p) }
    NULL }))
}

th_mode <- par_to_theta(full[free]); nt <- names(th_mode); k <- length(nt)
if (FROM == "joint") {
  ## archive: 400 draws of each earlier chain, both regions represented; starts: their last states
  ch <- lapply(prev, readRDS)
  Z0 <- do.call(rbind, lapply(ch, function(x) { i <- unique(round(seq(1, nrow(x$draws), length.out = min(400, nrow(x$draws)))))
    t(apply(x$draws[i, free, drop = FALSE], 1, par_to_theta)) }))
  colnames(Z0) <- nt
  X0 <- t(sapply(ch, function(x) par_to_theta(x$draws[nrow(x$draws), free]))); colnames(X0) <- nt
  cat(sprintf("joint set: %d chains from %s, archive of %d states\n", NCH, paste(basename(prev), collapse = " "), nrow(Z0)))
} else {
  ## initial archive on the sampler coordinates, from the Laplace approximation at the mode
  J <- numDeriv::jacobian(function(v) { names(v) <- free; as.numeric(par_to_theta(v)) }, full[free])
  Vt <- J %*% mode$V[free, free] %*% t(J); dimnames(Vt) <- list(nt, nt); Vt <- (Vt + t(Vt)) / 2
  tu <- if ("kg_u" %in% nt) "kg_u" else "mx_u"; su <- if (tu == "kg_u") exp(full[["g0"]]) else 1
  Vt[tu, tu] <- max(Vt[tu, tu], (0.05 * su)^2)
  set.seed(77); L <- chol_safe(4 * Vt)
  Z0 <- t(replicate(10 * k, { th <- th_mode + as.numeric(rnorm(k) %*% L); th[[tu]] <- th[[tu]] + su * runif(1, -0.5, 0.5); th }))
  colnames(Z0) <- nt
  ## starting states: archive members with finite posterior density
  X0 <- matrix(NA_real_, NCH, k, dimnames = list(NULL, nt)); j <- 0
  for (r in sample.int(nrow(Z0))) {
    p <- theta_to_par(Z0[r, ]); if (is.null(p) || !is.finite(logpost(p))) next
    j <- j + 1; X0[j, ] <- Z0[r, ]; if (j == NCH) break
  }
  if (j < NCH) stop("could not find ", NCH, " starting states")
}
P_JUMP <- if (FROM == "joint") as.numeric(Sys.getenv("WL_PJUMP", "0.4")) else NULL   # the joint set jumps between the regions often
res <- run_demcz(logpost, Z0, X0, N_GEN, N_BURN, THIN, ckpt = sprintf("out/chains/%s_demcz_ckpt.rds", tag), seed = switch(FROM, second = 12, joint = 13, 11),
                 ckpt_every = 100, verbose_every = 500, label = tag, cl = cl, p_jump = P_JUMP)
for (i in seq_len(NCH)) {
  P <- res$chains[[i]]
  D <- matrix(full, nrow(P), length(full), byrow = TRUE, dimnames = list(NULL, names(full))); D[, colnames(P)] <- P
  saveRDS(list(draws = D, lp = res$lp[[i]], acc = res$acc[i], n_iter = res$n_gen, n_burn = res$n_burn, thin = res$thin,
               hours = res$hours, ref = d$ref, sex = sex, lift = lift, chain = i, sampler = "DE-MCzs", from = if (FROM %in% c("second", "joint")) FROM else "kept"), outf[i])
}
if (!is.null(cl)) parallel::stopCluster(cl)
cat("written", outf, "\n")   # the checkpoint is kept: deleting the outputs and rerunning with more generations extends the run
