## 31_descriptives.R
## -----------------------------------------------------------------------
## Descriptive tables for the paper, computed on the analysis sample (senior
## competition, recorded bodyweight, cleaned units). Outputs in out/desc_*.csv:
##   coverage        competitions, athlete-day-lifts, athletes by period
##   unit_types      bracketed / clean sweep / bomb-out shares by sex and lift
##   attempt_rates   success rate by attempt number
##   jump_rates      success rate of the next attempt by situation and jump
##   climb           distribution of best made minus opener
##   exponents       allometric exponents of best made lift on bodyweight
##   concentration   distinct athletes among the top athlete-days
source("R/30_marginal_lib.R")
att <- read_csv("data/attempts_clean.csv", show_col_types = FALSE,
                col_types = cols(category = col_character(), .default = col_guess()))

units <- list(); rows <- list()
for (sex in c("Men", "Women")) for (lift in c("snatch", "cleanjerk")) {
  d <- build_data(att, sex, lift, G = 8)
  u <- d$u |> mutate(sex = sex, lift = lift)
  best <- apply(ifelse(d$S == 1, d$W, NA), 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  u$best <- best
  u$type <- case_when(is.na(best) ~ "all missed",
                      d$has_miss & d$wmiss_max > best ~ "bracketed",
                      TRUE ~ "ended on a success")
  units[[paste(sex, lift)]] <- u
}
U <- bind_rows(units)

## coverage (athlete-day-lifts; competitions and athletes within period)
cov <- U |> mutate(period = cut(year, c(1997, 2004, 2011, 2016, 2018, 2026),
                                labels = c("1998-2004", "2005-2011", "2012-2016", "2017-2018", "2019-2026"))) |>
  group_by(period) |>
  summarise(competitions = n_distinct(event_id), athlete_day_lifts = n(),
            athletes = n_distinct(paste(sex, athlete)), .groups = "drop")
write_csv(cov, "out/desc_coverage.csv"); print(as.data.frame(cov))
cat(sprintf("total: %d athlete-day-lifts, %d attempts, %d athletes, %d competitions\n",
            nrow(U), sum(U$m), n_distinct(U$athlete), n_distinct(U$event_id)))

## unit types
ut <- U |> count(sex, lift, type) |> group_by(sex, lift) |> mutate(share = n / sum(n)) |> ungroup()
write_csv(ut, "out/desc_unit_types.csv"); print(as.data.frame(ut))

## success by attempt number
ar <- U |> group_by(sex, lift) |>
  summarise(a1 = mean(s1), a2 = mean(s2, na.rm = TRUE), a3 = mean(s3, na.rm = TRUE), .groups = "drop")
write_csv(ar, "out/desc_attempt_rates.csv"); print(as.data.frame(ar))

## success of the next attempt by situation and jump size
jr <- bind_rows(
  U |> filter(!is.na(s2)) |> transmute(sex, lift, attempt = 2, prev = s1, jump = w2 - w1, success = s2),
  U |> filter(!is.na(s3)) |> transmute(sex, lift, attempt = 3, prev = s2, jump = w3 - w2, success = s3)) |>
  mutate(situation = ifelse(prev == 1, "after a make", ifelse(jump == 0, "repeat after a miss", "raise after a miss")),
         jump_band = cut(jump, c(-1, 0, 1, 2, 3, 4, 5, 7, 10, 1000),
                         labels = c("0", "1", "2", "3", "4", "5", "6-7", "8-10", ">10"))) |>
  group_by(sex, lift, attempt, situation, jump_band) |>
  summarise(n = n(), rate = mean(success), .groups = "drop")
write_csv(jr, "out/desc_jump_rates.csv")
print(as.data.frame(jr |> filter(situation == "after a make", attempt == 2) |> dplyr::select(-attempt)))

## climb from the opener to the best made lift, given a successful opener
cl <- U |> filter(s1 == 1) |> mutate(climb = best - w1) |>
  group_by(sex, lift) |>
  summarise(n = n(), p0 = mean(climb == 0), p1_3 = mean(climb >= 1 & climb <= 3),
            p4_6 = mean(climb >= 4 & climb <= 6), p7_10 = mean(climb >= 7 & climb <= 10),
            p11plus = mean(climb > 10), median = median(climb), rel_median = median(climb / w1), .groups = "drop")
write_csv(cl, "out/desc_climb.csv"); print(as.data.frame(cl))

## allometric exponents (best made lift on bodyweight, with a linear year trend)
ex <- U |> filter(!is.na(best)) |> group_by(sex, lift) |>
  group_modify(~{ f <- lm(log(best) ~ log(bodyweight) + year, data = .x)
    tibble(n = nrow(.x), exponent = coef(f)[2], se = sqrt(vcov(f)[2, 2])) }) |> ungroup()
write_csv(ex, "out/desc_exponents.csv"); print(as.data.frame(ex))

## concentration of the top athlete-days among athletes, by residual rank
conc <- U |> filter(!is.na(best)) |> group_by(sex, lift) |>
  group_modify(~{ f <- lm(best ~ poly(log(bodyweight), 2) + year, data = .x)
    r <- residuals(f); o <- order(-r)
    tibble(top20 = n_distinct(.x$athlete[o[1:20]]), top50 = n_distinct(.x$athlete[o[1:50]]),
           top100 = n_distinct(.x$athlete[o[1:100]]), top1pct = n_distinct(.x$athlete[o[1:round(0.01 * nrow(.x))]]),
           n_top1pct = round(0.01 * nrow(.x)),
           max_share_top100 = max(table(.x$athlete[o[1:100]])) / 100) }) |> ungroup()
write_csv(conc, "out/desc_concentration.csv"); print(as.data.frame(conc))
