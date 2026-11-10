## 86_combine_chains.R
## -----------------------------------------------------------------------
## Pools the chains written by 83_final_chain.R and 85_ladder_chain.R, with
## convergence diagnostics (rank-normalized split R-hat, bulk and tail
## effective sample sizes) for every parameter and for the endpoint and the
## one-in-a-thousand level at the reference bodyweight in 2026.
## When the chains of the final model move only slowly between two regions
## of the posterior (out/basin_weights.csv and out/basin_draws_<sex>_<lift>.rds
## of 96_basin_weights.R), the draws of each region are pooled over all the
## chains of the three sets (kept, second and joint), the diagnostics are
## computed within each region (R-hat over the chains with at least 1,000
## draws in the region, cut to a common length, since a shorter visit holds
## too few effective draws for a comparison between chains; effective sample
## sizes summed over those chains) and the pooled draws give each region the weight estimated
## from the joint set: all the draws of the region that is shortest relative
## to its weight are kept, and an evenly spaced subset of the other.
## Outputs, in the names the summary scripts read:
##   out/mcmcx_e2_<sex>_<lift>_rtau_sm3_mix_fix_log_tau_-3.91202.rds  (final model)
##   out/ladder_<rung>_<sex>_<lift>.rds                              (ladder)
##   out/diag_chains.csv   one row per fit and parameter
##   out/diag_summary.csv  one row per fit: chains, draws, max R-hat, min ESS
source("R/50_ext_lib.R"); source("R/80_sampler.R")
SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
diag_rows <- list(); summ_rows <- list()
derived <- function(D, ref, sex) {
  bw <- if (sex == "Men") 110 else 86; xd <- design_at(ref, bw, 2026)
  t(apply(D, 1, function(x) { p <- unpack2(fill_par2(x)); mu <- sum(p$beta * xd); om <- omega_of(p, matrix(xd, 1))
    ep <- unname(mu + om * std_endpoint(p)); c(endpoint = if (is.finite(ep)) ep else 1e6, q999 = unname(mu + om * std_quantile(p, 0.999))) }))
}
diagnose <- function(ch, label, extra_cols = list()) {
  n <- min(sapply(ch, function(x) nrow(x$draws)))
  draws <- lapply(ch, function(x) x$draws[seq_len(n), , drop = FALSE])
  free <- colnames(draws[[1]])[apply(do.call(rbind, draws), 2, function(v) all(is.finite(v)) && sd(v) > 0)]
  dg <- diag_chains(lapply(draws, function(D) D[, free, drop = FALSE]))
  ## derived quantities, on a thinned subset of each chain
  idx <- unique(round(seq(1, n, length.out = min(n, 1000))))
  dv <- lapply(draws, function(D) derived(D[idx, , drop = FALSE], ch[[1]]$ref, ch[[1]]$sex))
  dg2 <- diag_chains(dv)
  dg <- rbind(dg, dg2)
  diag_rows[[length(diag_rows) + 1]] <<- data.frame(fit = label, dg)
  summ_rows[[length(summ_rows) + 1]] <<- data.frame(fit = label, chains = length(ch), draws_per_chain = n,
    iterations = ch[[1]]$n_iter, burn = ch[[1]]$n_burn, thin = ch[[1]]$thin,
    acceptance = paste(sapply(ch, function(x) paste(sprintf("%.2f", x$acc), collapse = "/")), collapse = " "),
    max_rhat = max(dg$rhat, na.rm = TRUE), min_ess_bulk = min(dg$ess_bulk, na.rm = TRUE),
    median_ess_bulk = median(dg$ess_bulk, na.rm = TRUE), ess_xi = dg$ess_bulk[dg$parameter == "xi"],
    ess_u = if ("mx_u" %in% dg$parameter) dg$ess_bulk[dg$parameter == "mx_u"] else NA,
    ess_endpoint = dg$ess_bulk[dg$parameter == "endpoint"], rhat_endpoint = dg$rhat[dg$parameter == "endpoint"],
    hours_per_chain = mean(sapply(ch, function(x) x$hours)))
  for (nm in names(extra_cols)) summ_rows[[length(summ_rows)]][[nm]] <<- extra_cols[[nm]]
  list(draws = draws, lp = lapply(ch, function(x) x$lp[seq_len(n)]), dg = dg, free = free, n = n)
}
## diagnostics within one region: pieces are the draws of each chain in the region
diagnose_pieces <- function(pieces, lps, label, ch, extra_cols = list(), n_min = 1000) {
  ref <- ch[[1]]$ref; sex <- ch[[1]]$sex
  all_draws <- do.call(rbind, pieces); all_lp <- unlist(lps)      # every draw of the region is pooled
  sizes <- sapply(pieces, nrow)
  if (sum(sizes >= n_min) < 2) n_min <- max(1, sort(sizes, decreasing = TRUE)[min(2, length(sizes))])   # short test runs: the two longest visits
  keep <- sizes >= n_min; pieces <- pieces[keep]                   # the diagnostics use the chains with enough draws in the region
  free <- colnames(pieces[[1]])[apply(do.call(rbind, pieces), 2, function(v) all(is.finite(v)) && sd(v) > 0)]
  one <- function(P) {   # R-hat over the pieces cut to a common length, ESS summed over the whole pieces
    rh <- diag_chains(P)$rhat
    es <- Reduce(`+`, lapply(P, function(D) { g <- diag_chains(list(D)); cbind(ess_bulk = g$ess_bulk, ess_tail = g$ess_tail) }))
    data.frame(parameter = colnames(P[[1]]), rhat = rh, es, row.names = NULL) }
  dg <- one(lapply(pieces, function(D) D[, free, drop = FALSE]))
  dv <- lapply(pieces, function(D) { idx <- unique(round(seq(1, nrow(D), length.out = min(nrow(D), 1000)))); derived(D[idx, , drop = FALSE], ref, sex) })
  dg <- rbind(dg, one(dv))
  diag_rows[[length(diag_rows) + 1]] <<- data.frame(fit = label, dg)
  summ_rows[[length(summ_rows) + 1]] <<- data.frame(fit = label, chains = length(pieces), draws_per_chain = round(mean(sapply(pieces, nrow))),
    iterations = ch[[1]]$n_iter, burn = ch[[1]]$n_burn, thin = ch[[1]]$thin,
    acceptance = paste(sapply(ch[keep], function(x) paste(sprintf("%.2f", x$acc), collapse = "/")), collapse = " "),
    max_rhat = max(dg$rhat, na.rm = TRUE), min_ess_bulk = min(dg$ess_bulk, na.rm = TRUE),
    median_ess_bulk = median(dg$ess_bulk, na.rm = TRUE), ess_xi = dg$ess_bulk[dg$parameter == "xi"],
    ess_u = if ("mx_u" %in% dg$parameter) dg$ess_bulk[dg$parameter == "mx_u"] else NA,
    ess_endpoint = dg$ess_bulk[dg$parameter == "endpoint"], rhat_endpoint = dg$rhat[dg$parameter == "endpoint"],
    hours_per_chain = mean(sapply(ch[keep], function(x) x$hours)))
  for (nm in names(extra_cols)) summ_rows[[length(summ_rows)]][[nm]] <<- extra_cols[[nm]]
  list(draws = all_draws, lp = all_lp, dg = dg, free = free)
}
combine <- function(files, out_file, label, extra = list(), basins = NULL, assignment = NULL) {
  ch <- lapply(files, readRDS)
  from <- sapply(ch, function(x) if (is.null(x$from)) "kept" else x$from)
  two <- !is.null(basins) && nrow(basins) == 2 && isTRUE(all(basins$two_regions)) && !is.null(assignment)
  if (two && any(is.na(basins$weight))) { cat(label, ": two regions but no joint set yet; the draws are pooled as they are\n"); two <- FALSE }
  if (!two) {
    g <- diagnose(ch, label); draws <- g$draws
    out <- c(list(draws = do.call(rbind, draws), chains = draws, chain_from = from, lp = unlist(g$lp),
                  acc = mean(sapply(ch, function(x) mean(x$acc))), free = g$free, ref = ch[[1]]$ref, sex = ch[[1]]$sex,
                  lift = ch[[1]]$lift, diag = g$dg, n_chains = length(ch)), extra)
  } else {
    ## two regions: every draw assigned to the nearer mode (96_basin_weights.R), diagnostics within
    ## each region, pooled draws weighted by the masses of the regions
    ord <- match(basename(files), assignment$files); reg <- assignment$region[ord]
    lab <- ifelse(basins$mode_u == min(basins$mode_u), "lower threshold", "higher threshold")
    ## the joint set on its own, as a check of the mixing between the regions
    if (any(from == "joint")) diagnose(ch[from == "joint"], paste0(label, ", joint set"))
    gs <- lapply(seq_len(nrow(basins)), function(b) {
      pieces <- lapply(seq_along(ch), function(i) ch[[i]]$draws[reg[[i]] == basins$basin[b], , drop = FALSE])
      lps <- lapply(seq_along(ch), function(i) ch[[i]]$lp[reg[[i]] == basins$basin[b]])
      diagnose_pieces(pieces, lps, paste0(label, ", ", lab[b]), ch, list(weight = basins$weight[b])) })
    nb <- sapply(gs, function(g) nrow(g$draws)); w <- basins$weight
    N <- min(sum(nb), floor(min(nb / pmax(w, 1e-12))))
    pick <- lapply(seq_along(gs), function(b) { D <- gs[[b]]$draws; L <- gs[[b]]$lp
      m <- min(nrow(D), round(w[b] * N)); i <- if (m > 0) unique(round(seq(1, nrow(D), length.out = m))) else integer(0)
      list(D = D[i, , drop = FALSE], L = L[i]) })
    top <- which.max(w)
    out <- c(list(draws = do.call(rbind, lapply(pick, `[[`, "D")), chains = lapply(ch, `[[`, "draws"), chain_from = from,
                  lp = unlist(lapply(pick, `[[`, "L")), acc = mean(sapply(ch, function(x) mean(x$acc))), free = gs[[top]]$free,
                  ref = ch[[1]]$ref, sex = ch[[1]]$sex, lift = ch[[1]]$lift, diag = gs[[top]]$dg, n_chains = length(ch),
                  regions = reg, basins = data.frame(basins, label = lab, draws_kept = sapply(pick, function(p) nrow(p$D)))), extra)
    cat(label, ": regions", paste(sprintf("%s %.3f", lab, w), collapse = ", "), "\n")
  }
  saveRDS(out, out_file); cat("written", out_file, "\n")
}
basin_tab <- if (file.exists("out/basin_weights.csv")) read_csv("out/basin_weights.csv", show_col_types = FALSE) else NULL
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  pre <- sprintf("out/chains/final_%s_%s", tolower(sex), lift)
  fs <- c(Sys.glob(paste0(pre, "_c[0-9]*.rds")), Sys.glob(paste0(pre, "_second_c[0-9]*.rds")), Sys.glob(paste0(pre, "_joint_c[0-9]*.rds")))
  if (length(fs)) {
    mode <- readRDS(sprintf("out/extmap_e2_%s_%s%s.rds", tolower(sex), lift, SUF))
    b <- if (!is.null(basin_tab)) as.data.frame(basin_tab[basin_tab$sex == sex & basin_tab$lift == lift, ]) else NULL
    fa <- sprintf("out/basin_draws_%s_%s.rds", tolower(sex), lift)
    combine(fs, sprintf("out/mcmcx_e2_%s_%s%s.rds", tolower(sex), lift, SUF), sprintf("final %s %s", sex, lift),
            list(jac = 0, map = mode, model = "e2", rtau = 1, G = 24), basins = b, assignment = if (file.exists(fa)) readRDS(fa) else NULL)
  }
  for (rung in c("L1", "L2")) {
    fs <- Sys.glob(sprintf("out/chains/ladder_%s_%s_%s_c*.rds", rung, tolower(sex), lift)); fs <- fs[!grepl("ckpt", fs)]
    if (length(fs)) combine(fs, sprintf("out/ladder_%s_%s_%s.rds", rung, tolower(sex), lift), sprintf("%s %s %s", rung, sex, lift), list(rung = rung))
  }
}
if (length(diag_rows)) {
  write_csv(bind_rows(diag_rows), "out/diag_chains.csv"); write_csv(bind_rows(summ_rows), "out/diag_summary.csv")
  print(as.data.frame(bind_rows(summ_rows) |> mutate(across(where(is.numeric), ~round(.x, 3)))))
}
