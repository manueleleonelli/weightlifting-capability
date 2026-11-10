## 96_basin_weights.R
## -----------------------------------------------------------------------
## The two regions of the posterior of the final model, for the combinations
## of sex and lift whose chains move only slowly between them, and the weight
## of each region in the pooled posterior.
## The chains of the final model are those started around the kept mode
## (final_<sex>_<lift>_c<k>), a second set started around the second mode of
## 84_final_select.R where run_all.R asks for one (final_<sex>_<lift>_second_c<k>),
## and, where the two regions are found, a joint set started from the last
## states of those chains with an archive holding draws of both regions
## (final_<sex>_<lift>_joint_c<k>, 88_final_demcz.R with WL_FROM=joint).
## Every retained draw is assigned to the nearer of the two modes in the
## threshold, the log scale g0 and the shape, each divided by its typical
## standard deviation within a chain (the median over the chains of the first
## two sets). When the chains of the first two sets agree on the share of
## their draws near each mode (the largest and the smallest share differ by at
## most 0.5), they move freely between the two regions and the draws are
## pooled as they are. When they disagree by more, the chains move between
## the regions too slowly for the pooled draws to carry the right share of
## each, and run_all.R runs the joint set, whose proposals follow the
## differences between the regions and which moves between them every few
## hundred generations. The share of the joint set's draws in each region is
## the weight of the region in the pooled posterior of 86_combine_chains.R,
## and its standard error comes from the shares of batches of the joint
## chains. Until the joint set exists the weights are missing and the file
## only records which combinations need it.
## Output: out/basin_weights.csv, one row per sex, lift and mode
##         out/basin_draws_<sex>_<lift>.rds, the assignment of every draw of
##         every chain to the two regions, where the regions are found
source("R/50_ext_lib.R"); source("R/80_sampler.R")
SUF <- "_rtau_sm3_mix_fix_log_tau_-3.91202"
V <- c("mx_u", "g0", "xi")
SPREAD <- 0.5          # chains whose shares near the second mode differ by more are treated as two regions
N_MIN <- 500           # a region needs at least this many draws (or a fifth of all draws) to count
N_BATCH <- 10          # batches per joint chain for the standard error of the weights
rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  pre <- sprintf("out/chains/final_%s_%s", tolower(sex), lift)
  fs <- list(kept = Sys.glob(paste0(pre, "_c[0-9]*.rds")), second = Sys.glob(paste0(pre, "_second_c[0-9]*.rds")),
             joint = Sys.glob(paste0(pre, "_joint_c[0-9]*.rds")))
  if (!length(fs$kept)) next
  base <- sprintf("out/extmap_e2_%s_%s%s", tolower(sex), lift, SUF)
  modes <- list(kept = readRDS(paste0(base, ".rds")))
  if (file.exists(paste0(base, "_second.rds"))) modes$second <- readRDS(paste0(base, "_second.rds"))
  ch <- lapply(unlist(fs), readRDS); set_of <- rep(names(fs), sapply(fs, length)); names(ch) <- basename(unlist(fs))
  n_all <- sum(sapply(ch, function(x) nrow(x$draws)))
  fa <- file.path("out", sprintf("basin_draws_%s_%s.rds", tolower(sex), lift))
  row_of <- function(b, n, share, shares = NA, spread = NA, two = FALSE, w = if (b == "kept") 1 else 0, se = NA, joint = NA, trans = NA)
    tibble(sex = sex, lift = lift, basin = b, mode_u = modes[[b]]$par[["mx_u"]], mode_xi = modes[[b]]$par[["xi"]],
           mode_logpost = unname(modes[[b]]$value), n_draws = n, share_draws = share, chain_shares = shares, spread = spread,
           two_regions = two, weight = w, weight_se = se, joint_chain_shares = joint, joint_transitions = trans)
  if (length(modes) < 2) { rows[[length(rows) + 1]] <- row_of("kept", n_all, 1); if (file.exists(fa)) file.remove(fa); next }
  ## assignment of every draw to the nearer mode, on the scale of the first two sets
  M <- t(sapply(modes, function(m) m$par[V]))
  S <- apply(t(sapply(ch[set_of != "joint"], function(x) apply(x$draws[, V, drop = FALSE], 2, sd))), 2, median)
  dist <- function(Y) matrix(sapply(rownames(M), function(m) sqrt(rowSums(sweep(sweep(Y, 2, M[m, ]), 2, S, "/")^2))), ncol = nrow(M))
  nearest <- function(Y) rownames(M)[apply(dist(Y), 1, which.min)]
  reg <- lapply(ch, function(x) nearest(x$draws[, V, drop = FALSE]))
  first <- set_of != "joint"
  shares <- sapply(reg[first], function(r) mean(r == "second")); spread <- max(shares) - min(shares)
  n_reg <- sapply(rownames(M), function(b) sum(unlist(reg) == b))
  sh_txt <- function(b, which = first) paste(sprintf("%.2f", sapply(reg[which], function(r) mean(r == b))), collapse = " ")
  cat(sprintf("%s %s: shares of the chains near the second mode %s (spread %.2f)\n", sex, lift, sh_txt("second"), spread))
  two <- spread > SPREAD && all(n_reg >= min(N_MIN, floor(0.2 * n_all)))
  if (!two) {
    for (b in rownames(M)) rows[[length(rows) + 1]] <- row_of(b, n_reg[[b]], n_reg[[b]] / n_all, sh_txt(b), spread, FALSE)
    if (file.exists(fa)) file.remove(fa); next
  }
  ## two regions: the weights are the shares of the joint set's draws in each
  w <- setNames(rep(NA_real_, nrow(M)), rownames(M)); se <- w; jsh <- NA; trans <- NA
  if (any(set_of == "joint")) {
    rj <- reg[set_of == "joint"]
    for (b in rownames(M)) {
      w[[b]] <- mean(unlist(rj) == b)
      bm <- unlist(lapply(rj, function(r) { q <- cut(seq_along(r), N_BATCH, labels = FALSE); tapply(r == b, q, mean) }))
      se[[b]] <- sd(bm) / sqrt(length(bm))
    }
    jsh <- sh_txt("second", set_of == "joint")
    trans <- sum(sapply(rj, function(r) sum(r[-1] != r[-length(r)])))
    cat(sprintf("%s %s: the joint set of %d chains moved between the regions %d times and spent %.3f (se %.3f) of its draws near the second mode\n",
                sex, lift, length(rj), trans, w[["second"]], se[["second"]]))
  } else cat(sex, lift, ": two regions; the joint set is needed for the weights\n")
  for (b in rownames(M)) rows[[length(rows) + 1]] <- row_of(b, n_reg[[b]], n_reg[[b]] / n_all, sh_txt(b), spread, TRUE, w[[b]], se[[b]], jsh, trans)
  saveRDS(list(files = names(ch), set = set_of, region = reg, modes = M, scale = S, coords = V, weight = w), fa)
}
out <- bind_rows(rows); write_csv(out, "out/basin_weights.csv"); print(as.data.frame(out))
