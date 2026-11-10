## 97_review_summary.R
## -----------------------------------------------------------------------
## Collects the simulation study of Section S2.2 (out/sim/*.rds) into tables.
## Outputs: out/sim_fits.csv        one row per replicate and fit
##          out/sim_summary.csv     by sex, lift and fitted anchor: truth, mean
##                                  and sd of the estimates, bias and RMSE
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })
sf <- Sys.glob("out/sim/fit_*.rds"); sf <- sf[!grepl("_ckpt", sf)]
if (length(sf)) {
  fits <- bind_rows(lapply(sf, function(f) { x <- readRDS(f)
    bind_rows(tibble(sex = x$sex, lift = x$lift, rep = x$rep, fit_tau = x$fit_tau, what = "truth", !!!as.list(x$q_true)),
              tibble(sex = x$sex, lift = x$lift, rep = x$rep, fit_tau = x$fit_tau, what = "fit", !!!as.list(x$q_fit),
                     conv = x$conv, evals = x$evals, minutes = x$minutes, loglik = x$loglik)) }))
  write_csv(fits, "out/sim_fits.csv")
  qn <- c("xi", "u", "Hu", "rhoa", "delta2", "delta3", "opener_mode", "raise_mode", "q999", "endpoint")
  long <- fits |> select(sex, lift, rep, fit_tau, what, all_of(qn)) |>
    pivot_longer(all_of(qn), names_to = "quantity") |> pivot_wider(names_from = what, values_from = value)
  summ <- long |> group_by(sex, lift, fit_tau, quantity) |>
    summarise(n_rep = n(), truth = first(truth), mean_fit = mean(fit), sd_fit = sd(fit),
              bias = mean(fit - truth), rmse = sqrt(mean((fit - truth)^2)), .groups = "drop")
  write_csv(summ, "out/sim_summary.csv")
  print(as.data.frame(summ |> mutate(across(where(is.numeric), ~ round(.x, 3)))))
}
