# How much can a human being lift?

Code and data for

> Leonelli, M. (2026). How much can a human being lift? A Bayesian latent capability model for Olympic weightlifting. Submitted.

One command reproduces every table and figure of the paper and of its supplement.

## Contents

```
run_all.R                 the driver: runs every step in order and writes the tables
R/                        the scripts called by run_all.R and the C++ likelihood kernel
data/attempts_clean.csv.gz  the archive of attempts
```

The data are the attempts of the international competitions sanctioned by the International Weightlifting Federation from 1998 to 2026, one row per attempt, processed from the result sheets compiled by the OpenWeightlifting project (https://github.com/euanwm/OpenWeightlifting). The columns give the competition (`event_id`, `event_name`, `date`, `year`, `is_youth`), the athlete (`athlete`, `nation`, `gender`, `bodyweight`, `category`), the lift (`lift`, `attempt_no`), the declared weight (`weight`), the outcome (`outcome`, good or fail) and the best lift reported in the result sheet (`reported_best`). The analysis keeps the senior competitions (`is_youth` false) and drops athlete-day-lifts whose weights break the rules of the sport (`prep_units` in `R/30_marginal_lib.R`); Section 2 of the paper describes the data.

Two features of the source are kept as they are. Athletes are identified by sex and name, and one women's athlete-day at the 2024 European Championships (64 kg) carries the name of a male lifter, so it counts as a separate athlete. At the 2023 Pan-American Games 86 per cent of the recorded bodyweights equal the category limit, which suggests that the result sheets give the limit in place of the weigh-in.

## Requirements

- R 4.1 or later. The packages Rcpp, dplyr, tidyr, readr, quantreg, numDeriv, MASS and posterior are installed by `run_all.R` if missing.
- A C++ compiler for Rcpp. On macOS run `xcode-select --install` once, and on Windows install Rtools.

## Running

From this folder, in a terminal:

```
Rscript run_all.R
```

On a Mac, `caffeinate -i Rscript run_all.R` keeps the computer awake while the lid is open. The run uses one process per physical core and takes about ten days on an eight-core laptop. Most of that time goes on the Markov chains of the final model. The longest are those of the women's lifts since 2017, which run three times as long as the others. The women's snatch gets a second set of chains from a second mode of its posterior and, because the chains move between the two regions of the posterior only slowly, a joint set of chains with an archive of draws from both regions, whose share of draws in each region weights the pooled posterior (Section S6). The fits with the shared self-assessment error (Section S5) take about ten hours each. The time estimate printed at the start is based on one likelihood evaluation on your machine. The progress is written to `runs/run_all.log`.

The run can be stopped at any time and restarted with the same command. Finished steps are skipped and the chains resume from their last checkpoint. A run that is interrupted and resumed follows a different path from an uninterrupted one, so the posterior summaries agree with the paper up to Monte Carlo error.

To check that every step runs, on a small part of the data and in about an hour:

```
RUN_ALL_TEST=1 Rscript run_all.R
```

The results of the test are meaningless and go to `runs/test_*` and `tables_test/`.

## Output

| Where | What |
|---|---|
| `tables/` | the body of each table, as LaTeX rows |
| `runs/since2017/figs/` | Figures 1 to 3 and the trace plots |
| `runs/since2017/out/` | every result for the athlete-days since 2017, the data of the paper |
| `runs/full/out/` | the same for the whole archive since 1998 (Section S7) |
| `results.zip` | all of the above, with the logs |

### Tables

| Table | File in `tables/` | Computed by |
|---|---|---|
| 1 Coverage of the archive | `table1_coverage.tex` | `31_descriptives.R` |
| 2 Openers against the previous best | `table2_openers.tex` | `56_ext_descriptives.R` |
| 3 Attempt and declaration parameters | `table3_attempts.tex` | `88_final_demcz.R`, `86_combine_chains.R`, `58_ext_summaries.R` |
| 4 Capability distribution | `table4_tail.tex` | as Table 3 |
| 5 The tail under the layers of the model | `table5_ladder.tex` | `85_ladder_chain.R`, `63_bodyfree_pot.R`, `68_ladder_summary.R` |
| 6 Capability quantiles | `table6_quantiles.tex` | `58_ext_summaries.R` |
| 7 The ceiling | `table7_ceiling.tex` | `58_ext_summaries.R`, `69_anchor_mix.R` |
| S1 Accuracy of the quadrature | `tableS1_quadrature.tex` | `91_quadrature_check.R`, `90_shared_error.R` |
| S2 The simplified model along the anchor | `tableS2_ridge.tex` | `89_ridge.R`, `59_ext_figures.R` |
| S3 Simulation study | `tableS3_simulation.tex` | `94_simulation.R`, `97_review_summary.R` |
| S4 The tail without a body | `tableS4_bodyfree.tex` | `63_bodyfree_pot.R` |
| S5 Fixed and estimated threshold | `tableS5_fixed_threshold.tex` | `64_mixture_gpd.R` |
| S6 Predictive check of the attempts | `tableS6_ppc_attempts.tex` | `53_ext_checks.R`, `87_merge_ppc.R` |
| S7 Predictive check of the record book | `tableS7_ppc_records.tex` | `95_ppc_record.R` |
| S8 Raise after a miss | `tableS8_raise_after_miss.tex` | as Table 3 |
| S9 The final model at the three anchors | `tableS9_anchors.tex` | `52_ext_map.R`, `69_anchor_mix.R` |
| S10 Prior sensitivity | `tableS10_priors.tex` | `98_prior_reweight.R` |
| S11 Time trend in the scale | `tableS11_scale_trend.tex` | `70_scale_trend_check.R` |
| S12 Shared self-assessment error | `tableS12_shared_error.tex` | `52_ext_map.R`, `90_shared_error.R` |
| S13 Convergence of the chains | `tableS13_convergence.tex` | `86_combine_chains.R`, `96_basin_weights.R` |
| S14 The whole archive | `tableS14_archive.tex` | as Table 3, on both data sets |
| S15 The whole archive, layers of the model | `tableS15_archive_ladder.tex` | as Table 5, on the whole archive |

All tables are written by `R/99_tables.R`, which can be run again on its own from this folder after `run_all.R` has finished.

### Figures

| Figure | File in `runs/since2017/figs/` |
|---|---|
| 1 Success of the second attempt by size of the raise | `fig-jump.pdf` |
| 2 The simplified model along the anchor | `fig-ridge.pdf` |
| 3 Extremal dependence between the lifts | `fig-chimodel.pdf` |
| (not in the paper) Trace plots of the one-in-a-thousand level and the endpoint | `fig-trace.pdf` |

### Numbers quoted in the text

Numbers that appear in the text but in no table are in `runs/since2017/out/`:

| File | What |
|---|---|
| `postx_latent.csv` | capability against the best made lift, by how the day ended (Section 5.5) |
| `desc_unit_types.csv` | share of days ending on a made attempt (Section 5.5) |
| `record_check.csv` | the records against the ceiling at the lifter's bodyweight and year (Section 5.7) |
| `ppc_record_summary.csv` | the top of the best lifts in simulated competitions (Section 5.7, Discussion) |
| `copulax_summary.csv`, `chix_latent.csv`, `totalx_returnlevel.csv` | the copula, extremal dependence and the competition total (Section 5.6, Section S8) |
| `anchor_mix.csv` | difference between capability and the best made lift at the three anchors (Section 5.7) |
| `mapsummary/` | the posterior modes of every fit, including the floor on the miss probability (Section S5) |
| `mode_compare.csv` | the posterior modes of the final model from each start (Section S6) |
| `basin_weights.csv` | the two modes of the women's snatch, the shares of the chains near each, and the weights of the two regions from the joint set of chains (Section S6) |
| `desc_*.csv` | the descriptive statistics of Section 5.1 |

## The scripts

The numbers in the file names give the order in which the steps were written, not the order in which they run; `run_all.R` runs them in the right order.

| Script | Purpose |
|---|---|
| `30_marginal_lib.R`, `50_ext_lib.R`, `marginal_kernel.cpp` | the model: data preparation, the likelihood with capability integrated out by quadrature (Section S1), priors and helpers |
| `80_sampler.R` | the differential-evolution sampler with checkpoints |
| `81_ladder_lib.R` | the record-book fits of Table 5 |
| `93_sim_lib.R` | simulation of complete competitions from the model |
| `00_compile.R`, `00_timing.R` | compile the kernel; time one likelihood evaluation |
| `31_*`, `56_*` | descriptive statistics |
| `64_*`, `66_*`, `82_*`, `52_*`, `84_*` | starting values and posterior modes of the final model (`R/attempt_starts.rds` holds the starting values of the attempt parameters, `R/mode_starts.rds` a mode of the women's snatch since 2017 that the optimizer reaches only in some runs) |
| `88_*` (`83_*` as an alternative sampler), `85_*`, `96_*`, `86_*` | the chains of the final model and of the record-book fits, the two regions of the women's snatch posterior and their weights, the pooling of the chains |
| `58_*`, `68_*`, `69_*`, `92_*` | posterior summaries |
| `57_*` | the copula stage |
| `53_*`, `87_*`, `95_*` | predictive checks |
| `62_*`, `63_*`, `70_*`, `72_*` | comparisons and checks on the best lifts |
| `89_*`, `90_*`, `91_*`, `94_*`, `97_*`, `98_*`, `65_*` | the anchor, the shared error, the quadrature, the simulation study, the priors, kernel checks |
| `59_*`, `99_*` | figures and tables |
