## 92_record_check.R
## -----------------------------------------------------------------------
## The local records of Tables 6 and 7 against the ceiling at the lifter's
## own bodyweight. For every bodyweight of Table 6, the heaviest best lift
## since 2019 within five kilograms (the local record of the paper) and the
## heaviest since the start of the analysis set, with the bodyweight and year
## of the lifter, and the posterior median of the ceiling in 2026 at the
## nominal bodyweight and at the lifter's bodyweight, and for the heaviest
## lift since the start also at the lifter's bodyweight in the year it was made.
## Inputs: out/mcmcx_e2_<sex>_<lift><SUF>.rds, data/attempts_clean.csv
## Output: out/record_check.csv
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
SUF <- Sys.getenv("WL_TAG_SUFFIX", "_rtau_sm3_mix_fix_log_tau_-3.91202")
N_DRAW <- as.integer(Sys.getenv("WL_NDRAW", "4000"))
rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  f <- sprintf("out/mcmcx_e2_%s_%s%s.rds", tolower(sex), lift, SUF)
  if (!file.exists(f)) { cat("no chain for", sex, lift, "\n"); next }
  m <- readRDS(f); D <- m$draws; ref <- m$ref
  d <- build_data(att, sex, lift, G = 24)
  best <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  idx <- round(seq(1, nrow(D), length.out = min(N_DRAW, nrow(D))))
  P <- lapply(idx, function(i) unpack2(D[i, ]))
  ceil_at <- function(bw, year = 2026) median(sapply(P, function(p) { xd <- design_at(ref, bw, year)
    sum(p$beta * xd) + omega_of(p, matrix(xd, 1)) * std_endpoint(p) }))
  bws <- if (sex == "Men") c(61, 67, 73, 81, 89, 96, 102, 110, 130, 150) else c(49, 55, 59, 64, 71, 76, 81, 86, 100, 120)
  for (bw in bws) {
    h <- which(abs(d$u$bodyweight - bw) <= 5 & !is.na(best))
    h19 <- h[d$u$year[h] >= 2019]
    r19 <- h19[which.max(best[h19])]; rall <- h[which.max(best[h])]
    rows[[length(rows) + 1]] <- tibble(sex = sex, lift = lift, bodyweight = bw,
      record = best[r19], record_bw = d$u$bodyweight[r19], record_year = d$u$year[r19],
      record_all = best[rall], record_all_bw = d$u$bodyweight[rall], record_all_year = d$u$year[rall],
      ceiling_nominal = ceil_at(bw), ceiling_at_record_bw = ceil_at(d$u$bodyweight[r19]),
      ceiling_at_record_all_bw = ceil_at(d$u$bodyweight[rall]),
      ceiling_at_record_all_bw_year = ceil_at(d$u$bodyweight[rall], d$u$year[rall]))
    cat(sprintf("%s %s %d: record %g at %.2f (%d), ceiling %.1f nominal, %.1f at the lifter's bodyweight\n", sex, lift, bw,
                best[r19], d$u$bodyweight[r19], d$u$year[r19], rows[[length(rows)]]$ceiling_nominal, rows[[length(rows)]]$ceiling_at_record_bw))
  }
}
out <- bind_rows(rows)
write_csv(out, "out/record_check.csv")
