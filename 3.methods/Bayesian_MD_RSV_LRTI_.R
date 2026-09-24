# ============================================================
# BAYESIAN WANING EFFICACY ANALYSIS
# Endpoint: Nirsevimab - medically attended RSV-LRTI
# Data: MD_RSV_LRTI.rds
#
# FINAL thesis:
#   - follow-up: 0-150 days post-dose
#   - reconstructed IPD from published Kaplan-Meier curves
#   - published risk-table intervals: 0, 30, 60, 90, 120, 150 days
#   - risk-scale VE
#   - interval-level BINOMIAL likelihood
#   - pi_v,k = pi_p,k * (1 - VE_r(t_k))
#   - candidate waning forms:
#       1) Exponential
#       2) Erlang-3
#       3) Power-law
#   - model choice by LOOIC
#   - HMC/NUTS in Stan via cmdstanr
#   - 4 independent chains
#   - 20,000 iterations per chain, 1,000 warm-up 
#   - convergence target: R-hat < 1.01 and bulk ESS > 400
#   - posterior median + equal-tailed 95% CrI
#   - VE output every 15 days
#
#
# ============================================================


# ============================================================
# 0. USER SETTINGS
# ============================================================

DATA_DIR <- "D:/Desktop/Tiwonge/Tiwonge"
DATA_FILE <- file.path(DATA_DIR, "MD_RSV_LRTI.rds")
OUT_DIR <- file.path(
  DATA_DIR,
  "analysis_outputs_MD_RSV_LRTI_BAYESIAN"
)

# specific endpoint settings
ANALYSIS_T_MAX <- 150
RISK_TABLE_BREAKS <- c(0, 30, 60, 90, 120, 150)
OUTPUT_GRID <- seq(0, 150, by = 15)

# four chains, 20,000 iterations each, 1,000 warm-up.
CHAINS <- 4
ITER_TOTAL <- 20000
ITER_WARMUP <- 1000
ITER_SAMPLING <- ITER_TOTAL - ITER_WARMUP

SEED <- 2024
ADAPT_DELTA <- 0.95
MAX_TREEDEPTH <- 12

# ------------------------------------------------------------
# PRIOR SETTINGS
# ------------------------------------------------------------
#
#   VE0 ~ Beta(a0,b0)
#   kappa ~ Half-Normal(0, sigma_kappa^2)
#   beta_Gamma ~ Gamma(2,1)
#
# The prior mean for VE0 was aligned with the earliest
# published trial efficacy estimate. For this endpoint the published cumulative
# VE was 79.5%, so the prior mean is centred at 0.795.
#
# 
VE0_PRIOR_MEAN <- 0.795
VE0_PRIOR_CONCENTRATION <- 4
A0 <- VE0_PRIOR_MEAN * VE0_PRIOR_CONCENTRATION
B0 <- (1 - VE0_PRIOR_MEAN) * VE0_PRIOR_CONCENTRATION
#

SIGMA_KAPPA <- 0.02

# specified Erlang/Gamma prior
BETA_GAMMA_SHAPE <- 2
BETA_GAMMA_RATE <- 1

RHO_SHAPE <- 2
RHO_RATE <- 1

# Beta(1,1) is used transparently.
BASELINE_RISK_ALPHA <- 1
BASELINE_RISK_BETA <- 1


# ============================================================
# 1. PACKAGE SETUP
# ============================================================

required_packages <- c(
  "cmdstanr",
  "posterior",
  "loo",
  "ggplot2",
  "dplyr",
  "tidyr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required package '", pkg, "' is not installed.\n",
      "Install it before running this script."
    )
  }
}

library(cmdstanr)
library(posterior)
library(loo)
library(ggplot2)
library(dplyr)
library(tidyr)

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

cat("Data file:", DATA_FILE, "\n")
cat("Output directory:", OUT_DIR, "\n")

if (!file.exists(DATA_FILE)) {
  stop("Cannot find: ", DATA_FILE)
}


# ============================================================
# 2. LOAD DATA
# ============================================================

data_list <- readRDS(DATA_FILE)

cat("\nTop-level names in RDS object:\n")
print(names(data_list))

if (is.null(data_list$ipd)) {
  stop(
    "No $ipd element was found. The Bayesian analysis in the report uses ",
    "reconstructed individual event/censoring times."
  )
}

ipd_raw <- as.data.frame(data_list$ipd)

cat("\nIPD columns:\n")
print(names(ipd_raw))


# ============================================================
# 3. ROBUST VARIABLE IDENTIFICATION
# ============================================================

find_col <- function(df, candidates, label) {
  nms <- names(df)
  low <- tolower(nms)
  candidates_low <- tolower(candidates)

  exact <- match(candidates_low, low, nomatch = 0)
  exact <- exact[exact > 0]

  if (length(exact) > 0) {
    return(nms[exact[1]])
  }

  for (candidate in candidates_low) {
    hit <- grep(candidate, low, fixed = TRUE)
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }

  stop(
    "Could not identify ", label, ".\n",
    "Tried: ", paste(candidates, collapse = ", "), "\n",
    "Available columns: ", paste(nms, collapse = ", ")
  )
}

time_col <- find_col(
  ipd_raw,
  c("time", "t", "day", "days", "followup", "follow_up", "fu", "ftime"),
  "follow-up time"
)

event_col <- find_col(
  ipd_raw,
  c("event", "status", "delta", "d", "case", "infection", "infected", "outcome"),
  "event indicator"
)

arm_col <- find_col(
  ipd_raw,
  c("arm", "group", "treatment", "treat", "trt", "vacc", "vaccine"),
  "treatment arm"
)


# ============================================================
# 4. CLEAN IPD WITHOUT ALTERING THE SCIENTIFIC DATA
# ============================================================

standardize_arm <- function(x) {
  raw <- trimws(as.character(x))
  z <- tolower(raw)
  z_compact <- gsub("[_ -]+", "", z)

  out <- rep(NA_integer_, length(z))

  # Text labels
  out[grepl("plac|control|standard", z_compact)] <- 0L
  out[grepl("nirse|nir|active|treat|mab", z_compact)] <- 1L

  
  num <- suppressWarnings(as.numeric(raw))
  out[is.na(out) & !is.na(num) & num == 0] <- 0L
  out[is.na(out) & !is.na(num) & num == 1] <- 1L

  out
}

ipd <- data.frame(
  time = suppressWarnings(as.numeric(ipd_raw[[time_col]])),
  status = suppressWarnings(as.numeric(ipd_raw[[event_col]])),
  arm = standardize_arm(ipd_raw[[arm_col]])
)

bad_missing <- !is.finite(ipd$time) | is.na(ipd$status) | is.na(ipd$arm)

if (any(bad_missing)) {
  warning(
    sum(bad_missing),
    " IPD rows have missing/non-numeric time, status, or arm and will be removed."
  )
  ipd <- ipd[!bad_missing, , drop = FALSE]
}


#if (any(ipd$time < 0)) {
  #print(head(ipd[ipd$time < 0, , drop = FALSE], 20))
 # stop(
   # "Negative follow-up times were found. The old code converted them using ",
   # "abs(time), which changes the data. Correct the reconstruction/source data ",
    #"instead of reflecting negative times into positive follow-up."
 # )
#}


ipd$status <- ifelse(ipd$status > 0, 1L, 0L)

if (!all(c(0L, 1L) %in% unique(ipd$arm))) {
  stop(
    "Could not identify both placebo and nirsevimab arms after standardisation."
  )
}

# This endpoint has administrative follow-up to day 150.
# Records beyond 150 are administratively censored at 150.
beyond_150 <- ipd$time > ANALYSIS_T_MAX
if (any(beyond_150)) {
  ipd$time[beyond_150] <- ANALYSIS_T_MAX
  ipd$status[beyond_150] <- 0L
}

cat("\nClean IPD summary:\n")
print(table(Arm = ipd$arm, Event = ipd$status))
cat("N total:", nrow(ipd), "\n")
cat("Maximum follow-up:", max(ipd$time), "\n")
cat("Total reconstructed events:", sum(ipd$status), "\n")

# validation for this endpoint:
# N = 2350, total events = 70.
if (nrow(ipd) != 2350) {
  warning(
    "You reported 2,350 infants for this endpoint, but this RDS ",
    "contains ", nrow(ipd), " usable records."
  )
}

if (sum(ipd$status) != 70) {
  warning(
    "You reported 70 total events for this endpoint, but this RDS ",
    "contains ", sum(ipd$status), " events."
  )
}


# ============================================================
# 5. VERIFY PUBLISHED/RECONSTRUCTED NUMBERS AT RISK
# Day:         0   30   60   90  120  150
# Placebo:   786  772  756  737  729  724
# Nirsevimab 1564 1553 1546 1538 1527 1519
# ============================================================

reported_risk <- data.frame(
  day = c(0, 30, 60, 90, 120, 150),
  placebo_reported = c(786, 772, 756, 737, 729, 724),
  nirsevimab_reported = c(1564, 1553, 1546, 1538, 1527, 1519)
)

# "At risk" immediately before time t: observed time >= t.
reported_risk$placebo_ipd <- sapply(
  reported_risk$day,
  function(tt) sum(ipd$arm == 0L & ipd$time >= tt)
)

reported_risk$nirsevimab_ipd <- sapply(
  reported_risk$day,
  function(tt) sum(ipd$arm == 1L & ipd$time >= tt)
)

reported_risk$placebo_match <- (
  reported_risk$placebo_reported == reported_risk$placebo_ipd
)

reported_risk$nirsevimab_match <- (
  reported_risk$nirsevimab_reported == reported_risk$nirsevimab_ipd
)

cat("\nRisk-table validation:\n")
print(reported_risk)

write.csv(
  reported_risk,
  file.path(OUT_DIR, "risk_table_validation.csv"),
  row.names = FALSE
)

if (!all(reported_risk$placebo_match) ||
    !all(reported_risk$nirsevimab_match)) {
  warning(
    "The IPD risk table does not exactly reproduce the reported risk table. ",
    "Review the RDS before interpreting any Bayesian results as a replication."
  )
}


# ============================================================
# 6. BUILD THE FIVE SPECIFIC LIKELIHOOD INTERVALS
#
# Published risk table defines:
#   (0,30], (30,60], (60,90], (90,120], (120,150]
#
# For interval k:
#   n_p,k = placebo number at risk at interval start
#   n_v,k = nirsevimab number at risk at interval start
#   c_p,k = placebo events in interval
#   c_v,k = nirsevimab events in interval
#
# The Bayesian likelihood is:
#   c_p,k ~ Binomial(n_p,k, pi_p,k)
#   c_v,k ~ Binomial(n_v,k, pi_v,k)
#   pi_v,k = pi_p,k * [1 - VE_r(t_k)]
#
# t_k is represented by the midpoint of the interval for the parametric VE curve.
# ============================================================

make_interval_data <- function(dat, breaks) {

  stopifnot(length(breaks) >= 2)

  K <- length(breaks) - 1L

  out <- lapply(seq_len(K), function(k) {

    lo <- breaks[k]
    hi <- breaks[k + 1]

    # At risk immediately before the start of the interval.
    n_p <- sum(dat$arm == 0L & dat$time >= lo)
    n_v <- sum(dat$arm == 1L & dat$time >= lo)

    # Use (lo, hi] consistently. Exact day-zero events are not expected.
    c_p <- sum(
      dat$arm == 0L &
      dat$status == 1L &
      dat$time > lo &
      dat$time <= hi
    )

    c_v <- sum(
      dat$arm == 1L &
      dat$status == 1L &
      dat$time > lo &
      dat$time <= hi
    )

    data.frame(
      k = k,
      lo = lo,
      hi = hi,
      t_mid = (lo + hi) / 2,
      n_p = n_p,
      c_p = c_p,
      n_v = n_v,
      c_v = c_v
    )
  })

  out <- bind_rows(out)

  if (any(out$c_p > out$n_p) || any(out$c_v > out$n_v)) {
    stop("An interval has more events than participants at risk.")
  }

  out
}

interval_dat <- make_interval_data(ipd, RISK_TABLE_BREAKS)

cat("\nInterval-level Bayesian data:\n")
print(interval_dat)

cat("\nEvents represented in intervals:", sum(interval_dat$c_p + interval_dat$c_v), "\n")

if (sum(interval_dat$c_p + interval_dat$c_v) != sum(ipd$status)) {
  warning(
    "The interval event total differs from the IPD event total. ",
    "Check for events exactly at time zero or boundary-time coding."
  )
}

write.csv(
  interval_dat,
  file.path(OUT_DIR, "bayesian_interval_data.csv"),
  row.names = FALSE
)


# ============================================================
# 7. STAN MODEL CODE
# ============================================================

# One placebo baseline attack risk is estimated for each interval.
# Conditional on the parametric VE trajectory:
#   pi_v[k] = pi_p[k] * (1 - VE_k)
#
# Separate event likelihood contributions are retained in generated quantities
# for LOOIC calculation.

stan_common_data <- '
data {
  int<lower=1> K;
  array[K] int<lower=0> n_p;
  array[K] int<lower=0> c_p;
  array[K] int<lower=0> n_v;
  array[K] int<lower=0> c_v;
  vector<lower=0>[K] t_mid;

  real<lower=0> a0;
  real<lower=0> b0;

  real<lower=0> baseline_alpha;
  real<lower=0> baseline_beta;
}
'

stan_exponential <- paste0(
stan_common_data,
'
data {
  real<lower=0> sigma_kappa;
}
parameters {
  vector<lower=0, upper=1>[K] pi_p;
  real<lower=0, upper=1> VE0;
  real<lower=0> kappa;
}
model {
  VE0 ~ beta(a0, b0);
  kappa ~ normal(0, sigma_kappa);
  pi_p ~ beta(baseline_alpha, baseline_beta);

  for (k in 1:K) {
    real VE_k = VE0 * exp(-kappa * t_mid[k]);
    real pi_v = pi_p[k] * (1 - VE_k);

    c_p[k] ~ binomial(n_p[k], pi_p[k]);
    c_v[k] ~ binomial(n_v[k], pi_v);
  }
}
generated quantities {
  vector[2 * K] log_lik;

  for (k in 1:K) {
    real VE_k = VE0 * exp(-kappa * t_mid[k]);
    real pi_v = pi_p[k] * (1 - VE_k);

    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k], pi_p[k]);
    log_lik[K + k] = binomial_lpmf(c_v[k] | n_v[k], pi_v);
  }
}
'
)

stan_erlang3 <- paste0(
stan_common_data,
'
parameters {
  vector<lower=0, upper=1>[K] pi_p;
  real<lower=0, upper=1> VE0;
  real<lower=0> beta_Gamma;
}
model {
  VE0 ~ beta(a0, b0);
  beta_Gamma ~ gamma(', BETA_GAMMA_SHAPE, ', ', BETA_GAMMA_RATE, ');
  pi_p ~ beta(baseline_alpha, baseline_beta);

  for (k in 1:K) {
    real VE_k =
      VE0 * (1 - gamma_cdf(t_mid[k] | 3, beta_Gamma));

    real pi_v = pi_p[k] * (1 - VE_k);

    c_p[k] ~ binomial(n_p[k], pi_p[k]);
    c_v[k] ~ binomial(n_v[k], pi_v);
  }
}
generated quantities {
  vector[2 * K] log_lik;

  for (k in 1:K) {
    real VE_k =
      VE0 * (1 - gamma_cdf(t_mid[k] | 3, beta_Gamma));

    real pi_v = pi_p[k] * (1 - VE_k);

    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k], pi_p[k]);
    log_lik[K + k] = binomial_lpmf(c_v[k] | n_v[k], pi_v);
  }
}
'
)

stan_powerlaw <- paste0(
stan_common_data,
'
data {
  real<lower=0> sigma_kappa;
  real<lower=0> rho_shape;
  real<lower=0> rho_rate;
}
parameters {
  vector<lower=0, upper=1>[K] pi_p;
  real<lower=0, upper=1> VE0;
  real<lower=0> kappa;
  real<lower=0> rho;
}
model {
  VE0 ~ beta(a0, b0);
  kappa ~ normal(0, sigma_kappa);
  rho ~ gamma(rho_shape, rho_rate);
  pi_p ~ beta(baseline_alpha, baseline_beta);

  for (k in 1:K) {
    real VE_k =
      VE0 / (1 + pow(kappa * t_mid[k], rho));

    real pi_v = pi_p[k] * (1 - VE_k);

    c_p[k] ~ binomial(n_p[k], pi_p[k]);
    c_v[k] ~ binomial(n_v[k], pi_v);
  }
}
generated quantities {
  vector[2 * K] log_lik;

  for (k in 1:K) {
    real VE_k =
      VE0 / (1 + pow(kappa * t_mid[k], rho));

    real pi_v = pi_p[k] * (1 - VE_k);

    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k], pi_p[k]);
    log_lik[K + k] = binomial_lpmf(c_v[k] | n_v[k], pi_v);
  }
}
'
)

stan_files <- c(
  exponential = file.path(OUT_DIR, "nirsevimab_exponential.stan"),
  erlang3 = file.path(OUT_DIR, "nirsevimab_erlang3.stan"),
  powerlaw = file.path(OUT_DIR, "nirsevimab_powerlaw.stan")
)

writeLines(stan_exponential, stan_files["exponential"])
writeLines(stan_erlang3, stan_files["erlang3"])
writeLines(stan_powerlaw, stan_files["powerlaw"])


# ============================================================
# 8. STAN DATA LISTS
# ============================================================

base_stan_data <- list(
  K = nrow(interval_dat),
  n_p = as.integer(interval_dat$n_p),
  c_p = as.integer(interval_dat$c_p),
  n_v = as.integer(interval_dat$n_v),
  c_v = as.integer(interval_dat$c_v),
  t_mid = as.numeric(interval_dat$t_mid),

  a0 = A0,
  b0 = B0,

  baseline_alpha = BASELINE_RISK_ALPHA,
  baseline_beta = BASELINE_RISK_BETA
)

stan_data_exp <- c(
  base_stan_data,
  list(sigma_kappa = SIGMA_KAPPA)
)

stan_data_erl <- base_stan_data

stan_data_pow <- c(
  base_stan_data,
  list(
    sigma_kappa = SIGMA_KAPPA,
    rho_shape = RHO_SHAPE,
    rho_rate = RHO_RATE
  )
)


# ============================================================
# 9. COMPILE
# ============================================================

cat("\nCompiling Stan models...\n")

mod_exp <- cmdstan_model(stan_files["exponential"])
mod_erl <- cmdstan_model(stan_files["erlang3"])
mod_pow <- cmdstan_model(stan_files["powerlaw"])


# ============================================================
# 10. FIT WITH HMC/NUTS
# ============================================================

sample_model <- function(model, data, seed, name) {

  cat("\n============================================================\n")
  cat("Sampling:", name, "\n")
  cat("============================================================\n")

  model$sample(
    data = data,
    seed = seed,
    chains = CHAINS,
    parallel_chains = CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    adapt_delta = ADAPT_DELTA,
    max_treedepth = MAX_TREEDEPTH,
    refresh = 1000
  )
}

fit_exp <- sample_model(
  mod_exp,
  stan_data_exp,
  SEED + 1,
  "Exponential"
)

fit_erl <- sample_model(
  mod_erl,
  stan_data_erl,
  SEED + 2,
  "Erlang-3"
)

fit_pow <- sample_model(
  mod_pow,
  stan_data_pow,
  SEED + 3,
  "Power-law"
)

fits <- list(
  exponential = fit_exp,
  erlang3 = fit_erl,
  powerlaw = fit_pow
)


# ============================================================
# 11. CONVERGENCE DIAGNOSTICS
# ============================================================

diagnose_model <- function(fit, model_name, parameter_names) {

  s <- fit$summary(variables = parameter_names)

  diagnostic <- fit$diagnostic_summary()

  data.frame(
    model = model_name,
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_bulk_ess = min(s$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(s$ess_tail, na.rm = TRUE),
    divergences = sum(diagnostic$num_divergent),
    max_treedepth_hits = sum(diagnostic$num_max_treedepth),
    passes_report_rule =
      max(s$rhat, na.rm = TRUE) < 1.01 &&
      min(s$ess_bulk, na.rm = TRUE) > 400
  )
}

conv <- bind_rows(
  diagnose_model(
    fit_exp,
    "Exponential",
    c("VE0", "kappa")
  ),
  diagnose_model(
    fit_erl,
    "Erlang-3",
    c("VE0", "beta_Gamma")
  ),
  diagnose_model(
    fit_pow,
    "Power-law",
    c("VE0", "kappa", "rho")
  )
)

cat("\nConvergence summary:\n")
print(conv)

write.csv(
  conv,
  file.path(OUT_DIR, "convergence_summary.csv"),
  row.names = FALSE
)


# ============================================================
# 12. TRACE PLOTS - FOUR CHAINS
# ============================================================

save_trace_plot <- function(fit, parameters, title, filename) {

  draws <- fit$draws(
    variables = parameters,
    format = "df"
  )

  trace_data <- as.data.frame(draws) |>
    select(.chain, .iteration, all_of(parameters)) |>
    pivot_longer(
      cols = all_of(parameters),
      names_to = "parameter",
      values_to = "value"
    )

  p <- ggplot(
    trace_data,
    aes(
      x = .iteration,
      y = value,
      colour = factor(.chain),
      group = .chain
    )
  ) +
    geom_line(alpha = 0.55, linewidth = 0.22) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
    labs(
      title = title,
      x = "Post-warm-up iteration",
      y = "Posterior draw",
      colour = "Chain"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(
    file.path(OUT_DIR, filename),
    p,
    width = 9,
    height = max(5, length(parameters) * 2.4),
    dpi = 300
  )
}

save_trace_plot(
  fit_exp,
  c("VE0", "kappa"),
  "Nirsevimab medically attended RSV-LRTI - Exponential trace plots",
  "trace_exponential.png"
)

save_trace_plot(
  fit_erl,
  c("VE0", "beta_Gamma"),
  "Nirsevimab medically attended RSV-LRTI - Erlang-3 trace plots",
  "trace_erlang3.png"
)

save_trace_plot(
  fit_pow,
  c("VE0", "kappa", "rho"),
  "Nirsevimab medically attended RSV-LRTI - Power-law trace plots",
  "trace_powerlaw.png"
)


# ============================================================
# 13. LOOIC MODEL SELECTION
# ============================================================

calculate_loo <- function(fit) {

  log_lik_matrix <- fit$draws(
    variables = "log_lik",
    format = "matrix"
  )

  loo(log_lik_matrix)
}

loo_exp <- calculate_loo(fit_exp)
loo_erl <- calculate_loo(fit_erl)
loo_pow <- calculate_loo(fit_pow)

loo_objects <- list(
  exponential = loo_exp,
  erlang3 = loo_erl,
  powerlaw = loo_pow
)

loo_table <- bind_rows(
  lapply(names(loo_objects), function(name) {

    x <- loo_objects[[name]]

    data.frame(
      model = name,
      elpd_loo = x$estimates["elpd_loo", "Estimate"],
      se_elpd_loo = x$estimates["elpd_loo", "SE"],
      p_loo = x$estimates["p_loo", "Estimate"],
      looic = -2 * x$estimates["elpd_loo", "Estimate"],
      se_looic = 2 * x$estimates["elpd_loo", "SE"]
    )
  })
) |>
  arrange(looic)

cat("\nLOOIC model selection (lower is better):\n")
print(loo_table)

write.csv(
  loo_table,
  file.path(OUT_DIR, "LOOIC_model_selection.csv"),
  row.names = FALSE
)

selected_model_name <- loo_table$model[1]
selected_fit <- fits[[selected_model_name]]

cat("\nSelected model:", selected_model_name, "\n")

# Final report says nirsevimab selected a single-stage exponential decay.
if (selected_model_name != "exponential") {
  warning(
    "The final report states that exponential/single-stage decay had the ",
    "lowest LOOIC for nirsevimab. This run selected ", selected_model_name, ".\n",
    "Do NOT force the answer to exponential. Instead check whether the RDS, ",
    "interval construction, or numerical prior hyperparameters differ from the ",
    "original analysis."
  )
}


# ============================================================
# 14. PARETO-k LOO DIAGNOSTICS
# ============================================================

pareto_summary <- bind_rows(
  lapply(names(loo_objects), function(name) {

    k <- pareto_k_values(loo_objects[[name]])

    data.frame(
      model = name,
      max_pareto_k = max(k, na.rm = TRUE),
      n_k_gt_0_7 = sum(k > 0.7, na.rm = TRUE),
      n_k_gt_1 = sum(k > 1, na.rm = TRUE)
    )
  })
)

cat("\nPareto-k diagnostics:\n")
print(pareto_summary)

write.csv(
  pareto_summary,
  file.path(OUT_DIR, "LOO_pareto_k_diagnostics.csv"),
  row.names = FALSE
)


# ============================================================
# 15. COMPUTE POSTERIOR VE(t)
# ============================================================

posterior_ve <- function(fit, form, time_grid) {

  draws <- as_draws_df(fit$draws())

  if (form == "exponential") {

    ve_draws <- sapply(
      time_grid,
      function(t) draws$VE0 * exp(-draws$kappa * t)
    )

  } else if (form == "erlang3") {

    ve_draws <- sapply(
      time_grid,
      function(t) {
        draws$VE0 *
          (1 - pgamma(t, shape = 3, rate = draws$beta_Gamma))
      }
    )

  } else if (form == "powerlaw") {

    ve_draws <- sapply(
      time_grid,
      function(t) {
        draws$VE0 /
          (1 + (draws$kappa * t)^draws$rho)
      }
    )

  } else {
    stop("Unknown model form: ", form)
  }

  if (is.null(dim(ve_draws))) {
    ve_draws <- matrix(ve_draws, ncol = length(time_grid))
  }

  data.frame(
    day = time_grid,
    posterior_median = apply(ve_draws, 2, median),
    lower_95_CrI = apply(
      ve_draws,
      2,
      quantile,
      probs = 0.025,
      names = FALSE
    ),
    upper_95_CrI = apply(
      ve_draws,
      2,
      quantile,
      probs = 0.975,
      names = FALSE
    ),
    posterior_mean = colMeans(ve_draws)
  )
}

ve_selected <- posterior_ve(
  selected_fit,
  selected_model_name,
  OUTPUT_GRID
)

cat("\nPosterior VE at 15-day intervals:\n")
print(ve_selected)

write.csv(
  ve_selected,
  file.path(OUT_DIR, "Bayesian_VE_15_day_grid.csv"),
  row.names = FALSE
)


# ============================================================
# 16. COMPARE OUTPUT WITH THE FINAL-REPORT BAYESIAN TABLE
#
# Final report Table 5.2, M4 Bayesian:
# Day 0   0.85 (0.81, 0.93)
# ...
# Day 150 0.65 (0.53, 0.77)
#
# This is a validation target, NOT hard-coded into estimation.
# ============================================================

report_bayesian <- data.frame(
  day = seq(0, 150, by = 15),
  report_median = c(
    0.85, 0.83, 0.81, 0.79, 0.77, 0.75,
    0.73, 0.71, 0.69, 0.67, 0.65
  ),
  report_lower = c(
    0.81, 0.76, 0.75, 0.72, 0.70, 0.67,
    0.65, 0.63, 0.60, 0.56, 0.53
  ),
  report_upper = c(
    0.93, 0.90, 0.88, 0.86, 0.84, 0.82,
    0.81, 0.80, 0.79, 0.78, 0.77
  )
)

validation <- left_join(
  ve_selected,
  report_bayesian,
  by = "day"
) |>
  mutate(
    median_difference =
      posterior_median - report_median,
    lower_difference =
      lower_95_CrI - report_lower,
    upper_difference =
      upper_95_CrI - report_upper
  )

cat("\nComparison with Bayesian values:\n")
print(validation)

write.csv(
  validation,
  file.path(OUT_DIR, "comparison_with_Table_5_2.csv"),
  row.names = FALSE
)


# ============================================================
# 17. PLOT SELECTED VE CURVE
# ============================================================

p_ve <- ggplot(
  ve_selected,
  aes(x = day)
) +
  geom_ribbon(
    aes(
      ymin = lower_95_CrI,
      ymax = upper_95_CrI
    ),
    alpha = 0.20
  ) +
  geom_line(
    aes(y = posterior_median),
    linewidth = 0.9
  ) +
  geom_point(
    aes(y = posterior_median),
    size = 1.8
  ) +
  coord_cartesian(
    xlim = c(0, 150),
    ylim = c(0, 1)
  ) +
  scale_x_continuous(
    breaks = seq(0, 150, by = 15)
  ) +
  labs(
    title = paste0(
      "Bayesian VE(t): Nirsevimab against medically attended RSV-LRTI\n",
      "LOOIC-selected form: ", selected_model_name
    ),
    x = "Days post-dose",
    y = "Risk-scale vaccine efficacy",
    caption = "Line: posterior median; band: equal-tailed 95% credible interval"
  ) +
  theme_bw()

ggsave(
  file.path(OUT_DIR, "Bayesian_waning_curve_selected_model.png"),
  p_ve,
  width = 9,
  height = 6,
  dpi = 300
)


# ============================================================
# 18. PARAMETER SUMMARIES
# ============================================================

selected_parameters <- switch(
  selected_model_name,
  exponential = c("VE0", "kappa"),
  erlang3 = c("VE0", "beta_Gamma"),
  powerlaw = c("VE0", "kappa", "rho")
)

parameter_summary <- selected_fit$summary(
  variables = selected_parameters,
  probs = c(0.025, 0.5, 0.975)
)

cat("\nSelected-model parameter summary:\n")
print(parameter_summary)

write.csv(
  parameter_summary,
  file.path(OUT_DIR, "selected_model_parameter_summary.csv"),
  row.names = FALSE
)


# ============================================================
# 19. PRIOR SENSITIVITY
#
# ============================================================

prior_sensitivity_plan <- data.frame(
  scenario = c(
    "base",
    "VE0_more_diffuse",
    "VE0_more_concentrated",
    "kappa_tighter",
    "kappa_wider"
  ),
  VE0_mean = VE0_PRIOR_MEAN,
  VE0_concentration = c(
    VE0_PRIOR_CONCENTRATION,
    2,
    8,
    VE0_PRIOR_CONCENTRATION,
    VE0_PRIOR_CONCENTRATION
  ),
  sigma_kappa = c(
    SIGMA_KAPPA,
    SIGMA_KAPPA,
    SIGMA_KAPPA,
    SIGMA_KAPPA / 2,
    SIGMA_KAPPA * 2
  )
)

write.csv(
  prior_sensitivity_plan,
  file.path(OUT_DIR, "prior_sensitivity_plan.csv"),
  row.names = FALSE
)


# ============================================================
# 20. SAVE ANALYSIS METADATA
# ============================================================

analysis_summary <- list(
  endpoint = "Nirsevimab - medically attended RSV-LRTI",
  horizon_days = ANALYSIS_T_MAX,
  risk_table_breaks = RISK_TABLE_BREAKS,
  output_grid = OUTPUT_GRID,

  report_alignment = list(
    risk_scale = TRUE,
    binomial_likelihood = TRUE,
    candidate_forms = c(
      "exponential",
      "erlang3",
      "powerlaw"
    ),
    model_selection = "LOOIC",
    sampler = "HMC/NUTS via cmdstanr",
    chains = CHAINS,
    iterations_per_chain = ITER_TOTAL,
    warmup_per_chain = ITER_WARMUP,
    rhat_target = 1.01,
    bulk_ess_target = 400
  ),

  prior_settings = list(
    VE0_prior_mean = VE0_PRIOR_MEAN,
    VE0_prior_concentration = VE0_PRIOR_CONCENTRATION,
    a0 = A0,
    b0 = B0,
    sigma_kappa = SIGMA_KAPPA,
    beta_Gamma_shape = BETA_GAMMA_SHAPE,
    beta_Gamma_rate = BETA_GAMMA_RATE,
    rho_shape = RHO_SHAPE,
    rho_rate = RHO_RATE,
    baseline_risk_alpha = BASELINE_RISK_ALPHA,
    baseline_risk_beta = BASELINE_RISK_BETA
  ),

  interval_data = interval_dat,
  risk_table_validation = reported_risk,
  convergence = conv,
  looic = loo_table,
  pareto_k = pareto_summary,
  selected_model = selected_model_name,
  VE = ve_selected,
  final_report_validation = validation
)

saveRDS(
  analysis_summary,
  file.path(OUT_DIR, "Bayesian_MD_RSV_LRTI_analysis_summary.rds")
)

# CmdStan's posterior draws are already stored in CSV output files.
# Record their paths for reproducibility.
cmdstan_files <- c(
  paste(
    "Exponential:",
    paste(fit_exp$output_files(), collapse = "; ")
  ),
  paste(
    "Erlang-3:",
    paste(fit_erl$output_files(), collapse = "; ")
  ),
  paste(
    "Power-law:",
    paste(fit_pow$output_files(), collapse = "; ")
  )
)

writeLines(
  cmdstan_files,
  file.path(OUT_DIR, "cmdstan_output_file_locations.txt")
)


# ============================================================
# 21. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("Endpoint: nirsevimab - medically attended RSV-LRTI\n")
cat("Follow-up: 0-150 days\n")
cat("Events in IPD:", sum(ipd$status), "\n")
cat("Selected model by LOOIC:", selected_model_name, "\n")
cat("Expected : exponential\n")
cat("Output directory:", OUT_DIR, "\n")
cat("============================================================\n")

# End of script
