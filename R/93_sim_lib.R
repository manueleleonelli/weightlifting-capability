## 93_sim_lib.R
## -----------------------------------------------------------------------
## Simulation of attempt sequences from the final model of Section 4, used by
## 94_simulation.R (simulation study) and 95_ppc_record.R (predictive check of
## the record book). For each athlete-day-lift of a data object d (from
## build_data_ext), with parameter vector par:
##   capability    C = sqrt(rho_a) V_a + sqrt(1 - rho_a) eps,  Z = mu + omega G^{-1}(Phi(C))
##   opener        R_1 ~ two-piece normal(m_1, s_1-, s_1+),    w_1 = round(Z exp(-R_1))
##   attempt j     after a success: a raise, margin from situation 2 truncated at
##                 log{Z / (w_{j-1} + 1)}; after a miss: the same weight with the
##                 observed repeat rate, otherwise a raise from situation 3
##   outcome       S_j ~ Bernoulli(F(log(Z / w_j) / tau + delta_j))
## Declared weights are whole kilograms, as in the data. Each athlete-day-lift
## keeps its observed number of attempts.

## two-piece normal: distribution function and quantile function (log margin)
tpn_cdf <- function(r, m, sm, sp) {
  ifelse(r < m, 2 * sm / (sm + sp) * pnorm((r - m) / sm),
         sm / (sm + sp) + 2 * sp / (sm + sp) * (pnorm((r - m) / sp) - 0.5))
}
tpn_q <- function(u, m, sm, sp) {
  a <- sm / (sm + sp)
  lo <- pmin(pmax(u * (sm + sp) / (2 * sm), 1e-300), 0.5)                  # lower half, u < a
  hi <- pmin(pmax(0.5 + (u - a) * (sm + sp) / (2 * sp), 0.5), 1 - 1e-16)   # upper half, u >= a
  ifelse(u < a, m + sm * qnorm(lo), m + sp * qnorm(hi))
}

## observed share of repeated weights after a miss, for attempts 2 and 3
repeat_rate <- function(d) {
  sapply(2:3, function(j) {
    h <- !is.na(d$W[, j]) & !is.na(d$S[, j - 1]) & d$S[, j - 1] == 0
    if (!any(h)) return(0.75)
    mean(d$W[h, j] == d$W[h, j - 1])
  })
}

## simulate weights and outcomes for every athlete-day-lift of d
sim_attempts <- function(d, par, p_rep, seed) {
  set.seed(seed)
  p <- unpack2(fill_par2(par)); n <- d$n
  mu <- as.numeric(d$X %*% p$beta); om <- omega_of(p, d$X)
  ath <- match(d$u$athlete, unique(d$u$athlete))
  rho <- if (p$rho_a > 0) p$rho_a else 0
  V <- rnorm(max(ath))[ath]
  C <- sqrt(rho) * V + sqrt(1 - rho) * rnorm(n)
  Z <- mu + om * std_quantile(p, pnorm(C))
  ## the wide lower components of the body put a little mass (about 0.2 per cent
  ## of athlete-days for men) below any weight ever lifted, even below zero:
  ## such capabilities are redrawn, i.e. capability is conditioned to exceed half
  ## the lightest best lift of the data
  zmin <- 0.5 * min(d$W[d$S == 1], na.rm = TRUE)
  for (it in 1:50) {
    low <- which(Z < zmin); if (!length(low)) break
    C[low] <- sqrt(rho) * V[low] + sqrt(1 - rho) * rnorm(length(low))
    Z[low] <- mu[low] + om[low] * std_quantile(p, pnorm(C[low]))
  }
  Z <- pmax(Z, zmin)
  m <- p$dec[1:3]; sp <- p$dec[4:6]; sm <- p$dec[7:9]
  natt <- rowSums(!is.na(d$W))
  W <- matrix(NA_real_, n, 3); S <- matrix(NA_integer_, n, 3)
  r1 <- tpn_q(runif(n), m[1], sm[1], sp[1])
  W[, 1] <- pmax(round(Z * exp(-r1)), 1)
  S[, 1] <- rbinom(n, 1, plogis(log(Z / W[, 1]) / p$tau + p$delta[1]))
  for (j in 2:3) {
    h <- natt >= j
    prev <- W[, j - 1]; made <- S[, j - 1] == 1
    rep_j <- h & !made & runif(n) < p_rep[j - 1]
    raise <- h & !rep_j
    k <- ifelse(made, 2, 3)
    t <- log(Z / (prev + 1))
    Ft <- tpn_cdf(t, m[k], sm[k], sp[k])
    r <- tpn_q(runif(n) * Ft, m[k], sm[k], sp[k])
    r <- ifelse(is.finite(r), pmin(r, t), t)
    w <- pmax(round(Z * exp(-r)), prev + 1)
    W[raise, j] <- w[raise]; W[rep_j, j] <- prev[rep_j]
    S[h, j] <- rbinom(sum(h), 1, plogis(log(Z[h] / W[h, j]) / p$tau + p$delta[j]))
  }
  list(W = W, S = S, Z = Z)
}

## a copy of d with new weights and outcomes and the derived fields of build_data
with_attempts <- function(d, W, S) {
  d$W <- W; d$S <- S
  d$has_miss <- rowSums(S == 0, na.rm = TRUE) > 0
  d$has_make <- rowSums(S == 1, na.rm = TRUE) > 0
  d$wmade_min <- suppressWarnings(apply(ifelse(S == 1, W, NA), 1, min, na.rm = TRUE))
  d$wmiss_max <- suppressWarnings(apply(ifelse(S == 0, W, NA), 1, max, na.rm = TRUE))
  Sint <- S; storage.mode(Sint) <- "integer"; Sint[is.na(Sint)] <- -1L; d$Sint <- Sint
  Dm <- matrix(0L, nrow(W), 3); Dm[, 1] <- 1L
  for (j in 2:3) {
    h <- !is.na(W[, j]); raise <- h & W[, j] > W[, j - 1]
    Dm[h & raise & S[, j - 1] == 1, j] <- 2L
    Dm[h & raise & S[, j - 1] == 0, j] <- 3L
  }
  d$D <- Dm
  d
}

## best made lift and whether the day ended on a success
best_made <- function(W, S) {
  M <- ifelse(!is.na(S) & S == 1, W, -Inf)
  b <- pmax(M[, 1], M[, 2], M[, 3]); b[!is.finite(b)] <- NA; b
}
ended_success <- function(W, S) {
  last <- ifelse(!is.na(S[, 3]), S[, 3], ifelse(!is.na(S[, 2]), S[, 2], S[, 1]))
  !is.na(last) & last == 1
}

## quantities at the reference bodyweight in 2026 for a parameter vector
ref_quantities <- function(par, ref, sex) {
  p <- unpack2(fill_par2(par)); bw <- if (sex == "Men") 110 else 86
  xd <- design_at(ref, bw, 2026); mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
  v <- c(xi = p$xi, u = if (isTRUE(p$mix)) p$mx$u else NA, Hu = mix_Hu(p), rhoa = p$rho_a,
         delta2 = p$delta[2], delta3 = p$delta[3], opener_mode = 100 * p$dec[1], raise_mode = 100 * p$dec[2],
         q999 = mu + om * std_quantile(p, 0.999), endpoint = mu + om * std_endpoint(p))
  setNames(as.numeric(v), c("xi", "u", "Hu", "rhoa", "delta2", "delta3", "opener_mode", "raise_mode", "q999", "endpoint"))
}

## posterior mode with fixed parameters, as fit_map_ext2 of 52_ext_map.R, with an
## optional log posterior and a checkpoint so that an interrupted fit resumes
fit_mode <- function(d, start, fixed, maxit = 300, reltol = 1e-9, logpost = log_post_ext, ckpt = NULL) {
  start <- fill_par2(start); start[names(fixed)] <- fixed
  if (!is.null(ckpt) && file.exists(ckpt)) {
    ck <- tryCatch(readRDS(ckpt), error = function(e) NULL)
    if (!is.null(ck) && is.finite(ck$value)) start[names(ck$par)] <- ck$par
  }
  free <- setdiff(PAR_NAMES2[is.finite(start[PAR_NAMES2])], names(fixed))
  tfull <- to_real2(start[PAR_NAMES2])
  n_eval <- 0; best_v <- -Inf; best_par <- NULL; t0 <- Sys.time()
  obj <- function(th) {
    t <- tfull; t[free] <- th; par <- from_real2(t)
    v <- logpost(par, d); n_eval <<- n_eval + 1
    if (is.finite(v) && v > best_v) { best_v <<- v; best_par <<- par }
    if (!is.null(ckpt) && n_eval %% 200 == 0 && !is.null(best_par)) { tmp <- paste0(ckpt, ".tmp"); saveRDS(list(par = best_par, value = best_v), tmp); file.rename(tmp, ckpt) }
    if (!is.finite(v)) 1e12 else -v
  }
  ps <- PSCALE2[free]; ps[names(ps) %in% names(BOUNDS2)] <- 1
  o <- optim(tfull[free], obj, method = "BFGS", control = list(parscale = ps, maxit = maxit, reltol = reltol, ndeps = rep(1e-3, length(free))))
  t <- tfull; t[free] <- o$par; par <- from_real2(t)
  if (is.finite(best_v) && best_v > -o$value) par <- best_par
  list(par = par, value = logpost(par, d), conv = o$convergence, evals = n_eval,
       loglik = sum(unit_loglik_ext(par, d)), minutes = as.numeric(Sys.time() - t0, units = "mins"))
}
