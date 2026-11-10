## 53_ext_checks.R
## -----------------------------------------------------------------------
## Predictive checks for the base model and the shared-error extension
## (e1), at their posterior modes from 52_ext_map.R: the probability that
## attempt j succeeds given everything observed before its outcome (the
## declared weights up to and including w_j, the earlier outcomes), computed
## as the ratio of marginal likelihoods with S_j set to 1 and 0 and the later
## attempts dropped. The athlete factor is switched off so that the
## prediction uses the day alone. Output:
##   out/ppc_ext_attempts.csv   success by attempt number, observed and predicted
##   out/ppc_ext_jumps.csv      success of the second and third attempts after a
##                              success, by size of the raise
##   out/ppc_ext_openers.csv    opener success by band of the opener relative to
##                              the athlete's best lift of the day
## Settings: WL_SEX, WL_LIFT, WL_G (default 32), WL_FITS (label=file,...), WL_PPC_TAG (a
## suffix for the output files, e.g. "_second" for the checks at a second mode).
source("R/50_ext_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
sex <- Sys.getenv("WL_SEX", "Men"); lift <- Sys.getenv("WL_LIFT", "snatch")
G <- as.integer(Sys.getenv("WL_G", "32"))
d <- build_data_ext(build_data(att, sex, lift, G = G), jac = as.integer(Sys.getenv("WL_JAC", "0")))
## one "athlete" per unit so that the kernel returns one value per unit
d$alen <- rep(1L, d$n); d$astart <- 0:(d$n - 1); d$athletes <- seq_len(d$n)

## fits to check: WL_FITS="label=path,label=path" (default: the base and e1 fits of this sex and lift)
fits <- list()
spec <- Sys.getenv("WL_FITS", "")
if (spec != "") {
  for (kv in strsplit(strsplit(spec, ",")[[1]], "=")) if (file.exists(kv[2])) fits[[kv[1]]] <- readRDS(kv[2])$par
} else for (mname in c("base", "e1")) {
  f <- sprintf("out/extmap_%s_%s_%s.rds", mname, tolower(sex), lift)
  if (file.exists(f)) fits[[mname]] <- readRDS(f)$par
}
stopifnot(length(fits) > 0)
d$jac <- as.integer(Sys.getenv("WL_JAC", "0"))

refresh <- function(x) {
  x$has_miss <- rowSums(x$Sint == 0) > 0; x$has_make <- rowSums(x$Sint == 1) > 0
  x$wmade_min <- suppressWarnings(apply(ifelse(x$Sint == 1, x$W, NA), 1, min, na.rm = TRUE))
  x$wmiss_max <- suppressWarnings(apply(ifelse(x$Sint == 0, x$W, NA), 1, max, na.rm = TRUE))
  x
}
pred_attempt <- function(par, j) {
  par["rhoa"] <- NA
  dA <- d; dB <- d
  keep <- !is.na(d$W[, j])
  if (j < 3) for (jj in (j + 1):3) { dA$W[, jj] <- NA; dB$W[, jj] <- NA; dA$Sint[, jj] <- -1L; dB$Sint[, jj] <- -1L; dA$D[, jj] <- 0L; dB$D[, jj] <- 0L }
  dA$Sint[, j] <- 1L; dB$Sint[, j] <- 0L
  dA <- refresh(dA); dB <- refresh(dB)
  lA <- unit_loglik_ext(par, dA); lB <- unit_loglik_ext(par, dB)
  p <- 1 / (1 + exp(lB - lA)); p[!keep] <- NA; p
}

ppa <- list(); ppj <- list(); ppo <- list()
best <- suppressWarnings(apply(ifelse(d$Sint == 1, d$W, NA), 1, max, na.rm = TRUE)); best[!is.finite(best)] <- NA
for (mname in names(fits)) {
  cat("model", mname, "\n")
  P <- sapply(1:3, function(j) pred_attempt(fits[[mname]], j))
  for (j in 1:3) ppa[[length(ppa) + 1]] <- tibble(model = mname, sex = sex, lift = lift, attempt = j,
    n = sum(!is.na(d$S[, j])), observed = mean(d$S[, j], na.rm = TRUE), predicted = mean(P[, j], na.rm = TRUE))
  for (j in 2:3) {
    jump <- d$W[, j] - d$W[, j - 1]; ok <- !is.na(d$S[, j]) & d$S[, j - 1] == 1
    jb <- cut(jump, c(-1, 0, 1, 2, 3, 4, 5, 7, 10, 1000), labels = c("0", "1", "2", "3", "4", "5", "6-7", "8-10", ">10"))
    for (b in levels(jb)) { h <- ok & jb == b & !is.na(jb)
      if (sum(h) > 0) ppj[[length(ppj) + 1]] <- tibble(model = mname, sex = sex, lift = lift,
        situation = if (j == 2) "second after a make" else "third after a make", jump = b,
        n = sum(h), observed = mean(d$S[h, j]), predicted = mean(P[h, j])) }
  }
  ## opener success by the opener's position relative to the day's best made lift (as a check of the margin model)
  rel <- (best - d$W[, 1]) / d$W[, 1]
  rb <- cut(rel, c(-Inf, 0, 0.02, 0.04, 0.06, 0.08, 0.12, Inf), labels = c("0 (opener is the best)", "0-2%", "2-4%", "4-6%", "6-8%", "8-12%", ">12%"))
  for (b in levels(rb)) { h <- !is.na(rb) & rb == b
    if (sum(h) > 0) ppo[[length(ppo) + 1]] <- tibble(model = mname, sex = sex, lift = lift, band = b, n = sum(h),
      observed = mean(d$S[h, 1]), predicted = mean(P[h, 1])) }
  ## days without a made lift
  h <- is.na(best)
  ppo[[length(ppo) + 1]] <- tibble(model = mname, sex = sex, lift = lift, band = "no made lift", n = sum(h),
    observed = mean(d$S[h, 1]), predicted = mean(P[h, 1]))
}
upd <- function(f, new) {
  if (file.exists(f)) { old <- read_csv(f, show_col_types = FALSE); old <- old[!(old$sex == sex & old$lift == lift), ]; new <- bind_rows(old, new) }
  write_csv(new, f)
}
## one set of files per sex and lift (parallel runs never write to the same file); 87_merge_ppc.R merges them
dir.create("out/ppc", showWarnings = FALSE)
PTAG <- Sys.getenv("WL_PPC_TAG", "")
for (k in c("attempts", "jumps", "openers")) write_csv(bind_rows(list(attempts = ppa, jumps = ppj, openers = ppo)[[k]]),
                                                      sprintf("out/ppc/%s_%s%s_%s.csv", tolower(sex), lift, PTAG, k))
writeLines("done", sprintf("out/ppc_%s_%s%s.done", tolower(sex), lift, PTAG))
print(as.data.frame(bind_rows(ppa) |> mutate(across(where(is.numeric), ~round(.x, 3)))))
print(as.data.frame(bind_rows(ppj) |> filter(situation == "second after a make") |> mutate(across(where(is.numeric), ~round(.x, 3)))))
print(as.data.frame(bind_rows(ppo) |> mutate(across(where(is.numeric), ~round(.x, 3)))))
