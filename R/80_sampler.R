## 80_sampler.R
## -----------------------------------------------------------------------
## Samplers for the posteriors of the paper, with checkpoints so that a long
## run can be interrupted and resumed.
##
## run_demcz   differential-evolution Markov chains with sampling from the
##             past and snooker updates (ter Braak and Vrugt 2008), several
##             chains run together. Used by run_all.R for the final model and
##             for the record-book rungs of the ladder.
## run_amwg    adaptive random-walk Metropolis within Gibbs (Haario et al.
##             2001; Roberts and Rosenthal 2009): each iteration updates all
##             coordinates jointly and then the body and tail on their own,
##             with Gaussian proposals whose covariance is the empirical
##             covariance of the chain during burn-in, scaled to an acceptance
##             rate of 0.234; adaptation stops at the end of burn-in.
##
## Why not a plain random walk: the threshold, the generalized Pareto scale,
## the mixture components (one of them with a weight of a few per cent) and
## the location and scale coefficients move together along curved ridges of
## the posterior, and a Gaussian random walk with one fixed proposal crawls
## along them (effective sample sizes of about twenty in 17,000 draws in the
## third revision). The population sampler takes the shape of the posterior
## from the archive of past states of all chains.
##
## Both work on coordinates theta in which the generalized Pareto scale sigma
## is replaced by the threshold-stable scale sigma* = sigma - xi u (Coles 2001,
## Section 4.3.4): if the tail above u is GPD(sigma, xi), the tail above any
## higher u' is GPD(sigma + xi (u' - u), xi), so (sigma*, xi) does not change
## when the threshold moves, and the endpoint on the standardized scale is
## u + sigma / (-xi) = sigma* / (-xi). The change of coordinates is exact: the
## target density carries the Jacobian, so the posterior is unchanged.

## ---- coordinates -------------------------------------------------------
## par: named parameter vector with mx_lsig = log sigma. theta: the sampler's
## coordinates. Two changes of coordinates, both exact (the target carries the
## log Jacobian, so the posterior is unchanged):
##   sigma* = sigma - xi u (threshold-stable GPD scale, stored as mx_sstar)
##   with WL_COORD = "kg", the body and tail are also expressed
##   in kilograms at the centre of the design (b = 0, t = 0), where the
##   location is b0 and the scale exp(g0): component centres kg_m = b0 +
##   exp(g0) m, log component scales kg_ls = g0 + ls, threshold kg_u = b0 +
##   exp(g0) u and GPD scale kg_ss = exp(g0) sigma*. The data fix the
##   components and the threshold in kilograms, so in these coordinates b0 and
##   g0 move the reference component alone instead of dragging the rest of
##   the distribution with them.
COORD <- Sys.getenv("WL_COORD", "sstar")
par_to_theta <- function(par) {
  th <- par
  if (isTRUE(is.finite(par["mx_u"]))) {
    sstar <- exp(par[["mx_lsig"]]) - par[["xi"]] * par[["mx_u"]]
    th[["mx_lsig"]] <- sstar; names(th)[names(th) == "mx_lsig"] <- "mx_sstar"
    if (COORD == "kg") {
      c0 <- par[["b0"]]; s0 <- exp(par[["g0"]])
      th[["mx_m2"]] <- c0 + s0 * par[["mx_m2"]]; th[["mx_m3"]] <- c0 + s0 * par[["mx_m3"]]
      th[["mx_ls2"]] <- par[["g0"]] + par[["mx_ls2"]]; th[["mx_ls3"]] <- par[["g0"]] + par[["mx_ls3"]]
      th[["mx_u"]] <- c0 + s0 * par[["mx_u"]]; th[["mx_sstar"]] <- s0 * sstar
      nm <- names(th); nm[nm %in% c("mx_m2", "mx_m3", "mx_ls2", "mx_ls3", "mx_u", "mx_sstar")] <-
        c(mx_m2 = "kg_m2", mx_m3 = "kg_m3", mx_ls2 = "kg_ls2", mx_ls3 = "kg_ls3", mx_u = "kg_u", mx_sstar = "kg_ss")[nm[nm %in% c("mx_m2", "mx_m3", "mx_ls2", "mx_ls3", "mx_u", "mx_sstar")]]
      names(th) <- nm
    }
  }
  th
}
theta_to_par <- function(th) {
  par <- th; nm <- names(th)
  if ("kg_u" %in% nm) {
    c0 <- th[["b0"]]; s0 <- exp(th[["g0"]])
    par[["kg_m2"]] <- (th[["kg_m2"]] - c0) / s0; par[["kg_m3"]] <- (th[["kg_m3"]] - c0) / s0
    par[["kg_ls2"]] <- th[["kg_ls2"]] - th[["g0"]]; par[["kg_ls3"]] <- th[["kg_ls3"]] - th[["g0"]]
    par[["kg_u"]] <- (th[["kg_u"]] - c0) / s0; par[["kg_ss"]] <- th[["kg_ss"]] / s0
    nm[match(c("kg_m2", "kg_m3", "kg_ls2", "kg_ls3", "kg_u", "kg_ss"), nm)] <- c("mx_m2", "mx_m3", "mx_ls2", "mx_ls3", "mx_u", "mx_sstar")
    names(par) <- nm
  }
  if ("mx_sstar" %in% names(par)) {
    sig <- par[["mx_sstar"]] + par[["xi"]] * par[["mx_u"]]
    if (!is.finite(sig) || sig <= 0) return(NULL)
    par[["mx_sstar"]] <- log(sig); names(par)[names(par) == "mx_sstar"] <- "mx_lsig"
  }
  par
}
## log |d par / d theta|
log_jacobian <- function(par) {
  if (!isTRUE(is.finite(par["mx_u"]))) return(0)
  -par[["mx_lsig"]] - if (COORD == "kg") 4 * par[["g0"]] else 0
}
## names of the body-and-tail block in the sampler coordinates
tail_block <- function(nms) nms[nms == "xi" | grepl("^(mx_|kg_)", nms)]

## ---- the sampler -------------------------------------------------------
## logpost(par) returns the log posterior of a named parameter vector;
## start: named par vector (free coordinates only, plus any fixed ones in
## 'fixed'); V0: starting proposal covariance on the theta scale for the
## free coordinates; blocks: list of name vectors (theta names);
## n_iter, n_burn, thin; ckpt: checkpoint file; seed.
run_amwg <- function(logpost, start, V0, blocks, n_iter, n_burn, thin = 10, ckpt = NULL, seed = 1,
                     ckpt_every = 1000, adapt_every = 200, verbose_every = 1000, label = "") {
  target <- function(th) {
    p <- theta_to_par(th); if (is.null(p)) return(-Inf)
    v <- logpost(p); if (!is.finite(v)) return(-Inf)
    v + log_jacobian(p)
  }
  th0 <- par_to_theta(start); free <- names(th0)
  if (!is.null(ckpt) && file.exists(ckpt)) {
    st <- readRDS(ckpt); cat(sprintf("%s resuming at iteration %d\n", label, st$it))
    assign(".Random.seed", st$seed_state, envir = .GlobalEnv)
  } else {
    set.seed(seed)
    k <- length(free)
    st <- list(it = 0, cur = th0, cur_lp = target(th0), lsc = sapply(blocks, function(b) log(2.38^2 / length(b))),
               Sig = lapply(blocks, function(b) V0[b, b, drop = FALSE]), acc = setNames(numeric(length(blocks)), names(blocks)),
               hist = matrix(NA_real_, 0, k, dimnames = list(NULL, free)), draws = matrix(NA_real_, 0, k, dimnames = list(NULL, free)),
               lp = numeric(0), n_iter = n_iter, n_burn = n_burn, thin = thin, blocks = blocks, t_used = 0)
    if (!is.finite(st$cur_lp)) stop("the start has zero posterior density")
  }
  ## adaptation history: every iteration during burn-in, kept as a growing buffer
  hist_buf <- matrix(NA_real_, st$n_burn, length(free), dimnames = list(NULL, free))
  if (nrow(st$hist)) hist_buf[seq_len(nrow(st$hist)), ] <- st$hist
  draws_buf <- matrix(NA_real_, ceiling(max(st$n_iter - st$n_burn, 0) / st$thin), length(free), dimnames = list(NULL, free))
  lp_buf <- numeric(nrow(draws_buf)); nd <- nrow(st$draws)
  if (nd) { draws_buf[seq_len(nd), ] <- st$draws; lp_buf[seq_len(nd)] <- st$lp }
  chols <- lapply(seq_along(blocks), function(bi) chol_safe(exp(st$lsc[bi]) * st$Sig[[bi]]))
  cur <- st$cur; cur_lp <- st$cur_lp; t0 <- Sys.time()
  save_state <- function(it) {
    st$t_used <<- st$t_used + as.numeric(Sys.time() - t0, units = "secs"); t0 <<- Sys.time()
    if (is.null(ckpt)) return(invisible())
    s2 <- st; s2$it <- it; s2$cur <- cur; s2$cur_lp <- cur_lp
    s2$hist <- hist_buf[seq_len(min(it, s2$n_burn)), , drop = FALSE]
    s2$draws <- draws_buf[seq_len(nd), , drop = FALSE]; s2$lp <- lp_buf[seq_len(nd)]
    s2$seed_state <- get(".Random.seed", envir = .GlobalEnv)
    tmp <- paste0(ckpt, ".tmp"); saveRDS(s2, tmp); file.rename(tmp, ckpt)
  }
  if (st$it < st$n_iter) for (it in (st$it + 1):st$n_iter) {
    for (bi in seq_along(blocks)) {
      b <- blocks[[bi]]; prop <- cur
      prop[b] <- cur[b] + as.numeric(rnorm(length(b)) %*% chols[[bi]])
      lp <- target(prop)
      a <- if (is.finite(lp)) min(1, exp(lp - cur_lp)) else 0
      if (runif(1) < a) { cur <- prop; cur_lp <- lp; st$acc[bi] <- st$acc[bi] + 1 }
      if (it <= st$n_burn) st$lsc[bi] <- st$lsc[bi] + (a - 0.234) / max(100, it)^0.6
    }
    if (it <= st$n_burn) {
      hist_buf[it, ] <- cur
      if (it %% adapt_every == 0 || it == st$n_burn) {
        for (bi in seq_along(blocks)) {
          b <- blocks[[bi]]
          if (it >= 1000) {   # empirical covariance from the second half of the burn-in so far
            H <- hist_buf[max(1, floor(it / 2)):it, b, drop = FALSE]
            S <- cov(H); if (all(is.finite(S))) st$Sig[[bi]] <- S + diag(1e-10 + 1e-6 * diag(S), length(b))
          }
          chols[[bi]] <- chol_safe(exp(st$lsc[bi]) * st$Sig[[bi]])
        }
      }
    } else if ((it - st$n_burn) %% st$thin == 0) {
      nd <- nd + 1; draws_buf[nd, ] <- cur; lp_buf[nd] <- cur_lp
    }
    if (it %% ckpt_every == 0) save_state(it)
    if (it %% verbose_every == 0) { pp <- theta_to_par(cur)
      cat(sprintf("%s it %7d  acc %s  lp %.1f  xi %.3f u %.3f  (%.1f h)\n", label, it,
                  paste(sprintf("%.2f", st$acc / it), collapse = "/"), cur_lp,
                  if ("xi" %in% names(pp)) pp[["xi"]] else NA, if ("mx_u" %in% names(pp)) pp[["mx_u"]] else NA,
                  (st$t_used + as.numeric(Sys.time() - t0, units = "secs")) / 3600)) }
  }
  save_state(st$n_iter)
  ## draws back on the parameter scale
  D <- draws_buf[seq_len(nd), , drop = FALSE]
  P <- t(apply(D, 1, function(th) theta_to_par(th)))
  if (nd == 1) P <- matrix(P, 1, dimnames = list(NULL, names(theta_to_par(D[1, ]))))
  list(draws = P, lp = lp_buf[seq_len(nd)], acc = st$acc / st$n_iter, Sig = st$Sig, lsc = st$lsc,
       n_iter = st$n_iter, n_burn = st$n_burn, thin = st$thin, hours = st$t_used / 3600)
}
chol_safe <- function(S) {
  S <- (S + t(S)) / 2
  out <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(out)) { ev <- eigen(S, symmetric = TRUE); ev$values <- pmax(ev$values, 1e-10 * max(ev$values, 1e-12))
    out <- chol(ev$vectors %*% diag(ev$values, length(ev$values)) %*% t(ev$vectors)) }
  out
}

## ---- convergence diagnostics -------------------------------------------
## rank-normalized split R-hat and bulk effective sample size (Vehtari et al.
## 2021), from the posterior package when it is installed, otherwise the
## classical split R-hat and an autocorrelation-based ESS summed over chains
diag_chains <- function(chains) {   # chains: list of draws matrices with the same columns
  nms <- colnames(chains[[1]]); n <- min(sapply(chains, nrow))
  arr <- array(NA_real_, c(n, length(chains), length(nms)), dimnames = list(NULL, NULL, nms))
  for (c in seq_along(chains)) arr[, c, ] <- chains[[c]][seq_len(n), nms]
  if (requireNamespace("posterior", quietly = TRUE)) {
    out <- t(sapply(nms, function(v) { x <- arr[, , v]; if (sd(x) == 0) return(c(rhat = NA, ess_bulk = NA, ess_tail = NA))
      c(rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x), ess_tail = posterior::ess_tail(x)) }))
  } else {
    ess1 <- function(v) { a <- acf(v, lag.max = min(length(v) - 1, 2000), plot = FALSE)$acf[-1]; length(v) / (1 + 2 * sum(a[cumsum(a < 0) == 0])) }
    out <- t(sapply(nms, function(v) { x <- arr[, , v, drop = FALSE][, , 1]; if (is.null(dim(x))) x <- matrix(x)
      if (sd(x) == 0) return(c(rhat = NA, ess_bulk = NA, ess_tail = NA))
      h <- floor(nrow(x) / 2); xs <- cbind(x[1:h, , drop = FALSE], x[(h + 1):(2 * h), , drop = FALSE])
      W <- mean(apply(xs, 2, var)); B <- h * var(colMeans(xs)); rh <- sqrt(((h - 1) / h * W + B / h) / W)
      c(rhat = rh, ess_bulk = sum(apply(x, 2, ess1)), ess_tail = NA) }))
  }
  data.frame(parameter = nms, out, row.names = NULL)
}

## ---- population sampler: DE-MCzs -------------------------------------------
## Differential-evolution Markov chains with sampling from the past (ter Braak
## and Vrugt 2008): n_chains chains run together; each proposes a move along
## the difference of two past states drawn from a common archive, on a random
## subset of the coordinates, with an occasional jump of full length (gamma =
## 1) between modes and a snooker update one time in ten. The proposals take
## the shape of the posterior from the archive, so strong and curved
## dependence between parameters costs much less than with a Gaussian
## random walk. Each chain is a valid Markov chain; R-hat compares them.
## Z0: initial archive (theta scale, at least 10 d rows); X0: starting states
## (n_chains x d); n_gen generations, the first n_burn discarded, every
## thin-th kept (thin is rounded up to a multiple of K).
## Every K generations the current state of each chain is appended to the
## archive, which is saved in the checkpoint; the returned draws are read off
## the archive, so a run resumed with a different n_gen or n_burn loses
## nothing. If the checkpoint is already past n_gen, the run stops there and
## returns everything sampled.
## cl: an optional cluster of parallel::makeCluster whose workers hold a function
## target(th) equal to the one below (log posterior with the Jacobian, -Inf where
## the state is invalid). The proposals of a generation depend only on the archive
## and on the chain's own state, so they are evaluated together on the workers;
## the random numbers are drawn in the main process and the run is reproducible.
## p_jump: probability of a jump of full length along the difference of two archive
## states on all coordinates (gamma = 1), the move that carries a chain from one
## region of the posterior to another; 0.03 in the default scheme (one proposal in
## ten has gamma = 1, one in three of those uses all coordinates), raised for the
## joint set of 88_final_demcz.R.
run_demcz <- function(logpost, Z0, X0, n_gen, n_burn, thin = 10, ckpt = NULL, seed = 1, K = 10,
                      ckpt_every = 200, verbose_every = 500, label = "", cl = NULL, p_jump = NULL) {
  target <- function(th) {
    p <- theta_to_par(th); if (is.null(p)) return(-Inf)
    v <- logpost(p); if (!is.finite(v)) return(-Inf)
    v + log_jacobian(p)
  }
  evaluate <- function(xps) if (is.null(cl)) vapply(xps, target, numeric(1)) else unlist(parallel::parLapply(cl, xps, "target"))
  nms <- colnames(Z0); d <- length(nms); NCH <- nrow(X0)
  if (!is.null(ckpt) && file.exists(ckpt)) {
    st <- readRDS(ckpt); cat(sprintf("%s resuming at generation %d\n", label, st$g))
    assign(".Random.seed", st$seed_state, envir = .GlobalEnv)
    if (is.null(st$nz0)) st$nz0 <- nrow(st$Z) - NCH * (st$g %/% K)          # checkpoints of the first version
    if (is.null(st$Zlp)) st$Zlp <- rep(NA_real_, nrow(st$Z))
  } else {
    set.seed(seed)
    lp0 <- apply(X0, 1, target); if (any(!is.finite(lp0))) stop("a start has zero posterior density")
    st <- list(g = 0, Z = Z0, Zlp = rep(NA_real_, nrow(Z0)), nz0 = nrow(Z0), X = X0, lpX = lp0,
               acc = numeric(NCH), nev = numeric(NCH), t_used = 0)
  }
  n_more <- max(n_gen - st$g, 0)
  Zb <- matrix(NA_real_, nrow(st$Z) + ceiling(n_more / K) * NCH + NCH, d, dimnames = list(NULL, nms)); nz <- nrow(st$Z); Zb[seq_len(nz), ] <- st$Z
  Zlpb <- c(st$Zlp, rep(NA_real_, nrow(Zb) - nz))
  X <- st$X; lpX <- st$lpX; t0 <- Sys.time(); CR <- c(1 / 3, 2 / 3, 1); g_done <- st$g
  scale_eps <- apply(st$Z, 2, sd) * 1e-4
  save_state <- function(g) {
    st$t_used <<- st$t_used + as.numeric(Sys.time() - t0, units = "secs"); t0 <<- Sys.time()
    if (is.null(ckpt)) return(invisible())
    s2 <- st; s2$g <- g; s2$X <- X; s2$lpX <- lpX; s2$Z <- Zb[seq_len(nz), , drop = FALSE]; s2$Zlp <- Zlpb[seq_len(nz)]
    s2$draws <- NULL; s2$lp <- NULL
    s2$seed_state <- get(".Random.seed", envir = .GlobalEnv)
    tmp <- paste0(ckpt, ".tmp"); saveRDS(s2, tmp); file.rename(tmp, ckpt)
  }
  if (n_more > 0) for (g in (st$g + 1):n_gen) {
    ## proposals of every chain, from the archive and the chain's own state
    xps <- vector("list", NCH); extra <- numeric(NCH); ok <- logical(NCH)
    for (i in seq_len(NCH)) {
      x <- X[i, ]; r <- sample.int(nz, 3)
      if (runif(1) < 0.1) {               # snooker update
        z <- Zb[r[1], ]; dirv <- x - z; nn <- sum(dirv^2)
        if (nn > 0) {
          pr <- function(v) sum(v * dirv) / nn * dirv
          xp <- x + runif(1, 1.2, 2.2) * (pr(Zb[r[2], ]) - pr(Zb[r[3], ]))
          extra[i] <- (d - 1) * (0.5 * log(sum((xp - z)^2)) - 0.5 * log(nn)); ok[i] <- TRUE
        } else xp <- x
      } else if (!is.null(p_jump) && runif(1) < p_jump) {   # jump of full length on all coordinates
        xp <- x + (Zb[r[1], ] - Zb[r[2], ]) + rnorm(d, 0, scale_eps); ok[i] <- TRUE
      } else {                              # parallel direction on a random subspace
        cr <- CR[sample.int(3, 1)]; sub <- which(runif(d) < cr); if (!length(sub)) sub <- sample.int(d, 1)
        gam <- if (runif(1) < 0.1) 1 else 2.38 / sqrt(2 * length(sub))
        xp <- x; xp[sub] <- x[sub] + gam * (Zb[r[1], sub] - Zb[r[2], sub]) + rnorm(length(sub), 0, scale_eps[sub]); ok[i] <- TRUE
      }
      xps[[i]] <- xp
    }
    lps <- rep(-Inf, NCH); if (any(ok)) lps[ok] <- evaluate(xps[ok])
    for (i in seq_len(NCH)) {
      a <- if (ok[i]) lps[i] - lpX[i] + extra[i] else -Inf
      st$nev[i] <- st$nev[i] + 1
      if (is.finite(a) && log(runif(1)) < a) { X[i, ] <- xps[[i]]; lpX[i] <- lps[i]; st$acc[i] <- st$acc[i] + 1 }
    }
    if (g %% K == 0) { Zb[nz + seq_len(NCH), ] <- X; Zlpb[nz + seq_len(NCH)] <- lpX; nz <- nz + NCH }
    g_done <- g
    if (g %% ckpt_every == 0) save_state(g)
    if (g %% verbose_every == 0) { pp <- theta_to_par(X[1, ])
      cat(sprintf("%s gen %7d  acc %s  lp %s  xi %.3f u %.3f  (%.1f h)\n", label, g, paste(sprintf("%.2f", st$acc / st$nev), collapse = "/"),
                  paste(sprintf("%.1f", lpX), collapse = "/"), if ("xi" %in% names(pp)) pp[["xi"]] else NA, if ("mx_u" %in% names(pp)) pp[["mx_u"]] else NA,
                  (st$t_used + as.numeric(Sys.time() - t0, units = "secs")) / 3600)) }
  }
  save_state(g_done)
  ## draws: archive rows of each chain after burn-in, every thin-th generation
  thin_k <- K * max(1, ceiling(thin / K))
  gens <- seq(K, K * ((nz - st$nz0) %/% NCH), by = K)
  keep <- gens[gens > n_burn & gens %% thin_k == 0]
  if (!length(keep)) stop("no generations after burn-in: increase the number of generations")
  chains <- lapply(seq_len(NCH), function(i) {
    rows <- st$nz0 + (keep / K - 1) * NCH + i
    P <- t(apply(Zb[rows, , drop = FALSE], 1, theta_to_par))
    if (length(rows) == 1) P <- matrix(P, 1, dimnames = list(NULL, names(theta_to_par(Zb[rows, ]))))
    P })
  lps <- lapply(seq_len(NCH), function(i) Zlpb[st$nz0 + (keep / K - 1) * NCH + i])
  list(chains = chains, lp = lps, acc = st$acc / pmax(st$nev, 1), n_gen = g_done, n_burn = n_burn, thin = thin_k,
       hours = st$t_used / 3600, archive_size = nz)
}
