## 56_ext_descriptives.R
## -----------------------------------------------------------------------
## Descriptives used by the identification argument: the opener relative to
## the athlete's best made lift of the previous 400 days, and the success of
## the opener by band of that ratio; by sex and lift. Also the between-day
## standard deviation of the log best lift of the same athlete.
## Output: out/desc_opener_prev.csv, out/desc_between_day.csv
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))
rows <- list(); bd <- list()
for (sx in c("Men", "Women")) for (lf in c("snatch", "cleanjerk")) {
  u <- att |> filter(gender == sx, lift == lf, !is_youth) |>
    mutate(success = as.integer(outcome == "good")) |>
    group_by(athlete, event_id, date) |>
    summarise(w1 = weight[attempt_no == 1][1], s1 = success[attempt_no == 1][1],
              best = suppressWarnings(max(weight[success == 1], na.rm = TRUE)), .groups = "drop") |>
    mutate(best = ifelse(is.finite(best), best, NA)) |> arrange(athlete, date)
  u <- u |> group_by(athlete) |>
    mutate(prev_best = sapply(seq_along(date), function(i) {
      h <- date < date[i] & date >= date[i] - 400
      if (any(h) && any(!is.na(best[h]))) max(best[h], na.rm = TRUE) else NA })) |> ungroup()
  r <- u |> filter(!is.na(prev_best), !is.na(w1)) |> mutate(ratio = w1 / prev_best)
  r <- r |> mutate(band = cut(ratio, c(0, 0.88, 0.91, 0.94, 0.97, 1.0, 1.03, 10),
                              labels = c("<0.88", "0.88-0.91", "0.91-0.94", "0.94-0.97", "0.97-1.00", "1.00-1.03", ">1.03")))
  rows[[length(rows) + 1]] <- r |> group_by(band) |>
    summarise(n = n(), opener_success = mean(s1), .groups = "drop") |>
    mutate(sex = sx, lift = lf, median_ratio_all = median(r$ratio), q10 = quantile(r$ratio, .1), q90 = quantile(r$ratio, .9), n_days = nrow(r))
  v <- u |> filter(!is.na(best)) |> group_by(athlete) |>
    filter(n() >= 4, as.numeric(max(date) - min(date)) <= 730) |>
    mutate(lb = log(best) - mean(log(best))) |> ungroup()
  bd[[length(bd) + 1]] <- tibble(sex = sx, lift = lf, n_athletes = n_distinct(v$athlete), n_days = nrow(v), sd_log_best = sd(v$lb))
}
write_csv(bind_rows(rows), "out/desc_opener_prev.csv"); write_csv(bind_rows(bd), "out/desc_between_day.csv")
print(as.data.frame(bind_rows(rows) |> select(sex, lift, band, n, opener_success) |> mutate(opener_success = round(opener_success, 3))))
print(as.data.frame(bind_rows(bd)))
