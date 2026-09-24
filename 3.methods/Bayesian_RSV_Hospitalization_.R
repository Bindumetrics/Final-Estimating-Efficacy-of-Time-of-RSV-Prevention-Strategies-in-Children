# ============================================================
# BAYESIAN WANING EFFICACY ANALYSIS
# Endpoint: Nirsevimab - RSV LRTI Hospital Admission
# Data: RSV_Hospitalization.rds
#
#   - follow-up: 0-150 days post-dose
#   - reconstructed IPD from published Kaplan-Meier curves
#   - published risk-table boundaries: 0, 30, 60, 90, 120, 150 days
#   - risk-scale VE
#   - interval-level BINOMIAL likelihood
#   - pi_v,k = pi_p,k * (1 - VE_r(t_k))
#   - candidate waning forms:
#       1) Exponential
#       2) Erlang-3
#       3) Power-law
#   - model choice by LOOIC
#   - HMC/NUTS in Stan via cmdstanr
#   - four independent chains
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
DATA_FILE <- file.path(DATA_DIR, "RSV_Hospitalization.rds")
OUT_DIR <- file.path(
  DATA_DIR,
  "analysis_outputs_RSV_Hospitalization_BAYES"
)

ANALYSIS_T_MAX <- 150
RISK_TABLE_BREAKS <- c(0, 30, 60, 90, 120, 150)
OUTPUT_GRID <- seq(0, 150, by = 15)

CHAINS <- 4
ITER_TOTAL <- 20000
ITER_WARMUP <- 1000
ITER_SAMPLING <- ITER_TOTAL - ITER_WARMUP

SEED <- 2024
ADAPT_DELTA <- 0.95
MAX_TREEDEPTH <- 12


# ============================================================
# 0A. PRIOR SETTINGS
# ============================================================

#   VE0 ~ Beta(a0,b0)
#   kappa ~ Half-Normal(0, sigma_kappa^2)
#   beta_Gamma ~ Gamma(2,1)
#
# The prior mean for VE0 was aligned with the earliest
# published trial estimate. For this endpoint, the published cumulative VE over
# 150 days was 77.3%, so the prior is centred at 0.773.
#
VE0_PRIOR_MEAN <- 0.773
VE0_PRIOR_CONCENTRATION <- 4
A0 <- VE0_PRIOR_MEAN * VE0_PRIOR_CONCENTRATION
B0 <- (1 - VE0_PRIOR_MEAN) * VE0_PRIOR_CONCENTRATION

SIGMA_KAPPA <- 0.02

# specified Erlang/Gamma prior
BETA_GAMMA_SHAPE <- 2
BETA_GAMMA_RATE <- 1


RHO_SHAPE <- 2
RHO_RATE <- 1

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

if (!file.exists(DATA_FILE)) {
  stop("Cannot find data file: ", DATA_FILE)
}

cat("Data file:", DATA_FILE, "\n")
cat("Output directory:", OUT_DIR, "\n")


# ============================================================
# 2. LOAD RDS
# ============================================================

data_list <- readRDS(DATA_FILE)

cat("\nTop-level names in RDS object:\n")
print(names(data_list))

if (is.null(data_list$ipd)) {
  stop(
    "No $ipd element was found. The Bayesian analysis is based ",
    "on reconstructed individual event/censoring times."
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
  cand_low <- tolower(candidates)

  exact <- match(cand_low, low, nomatch = 0)
  exact <- exact[exact > 0]

  if (length(exact) > 0) {
    return(nms[exact[1]])
  }

  for (candidate in cand_low) {
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
# 4. CLEAN IPD WITHOUT ALTERING ITS SCIENTIFIC CONTENT
# ============================================================

standardize_arm <- function(x) {
  raw <- trimws(as.character(x))
  z <- tolower(raw)
  z_compact <- gsub("[_ -]+", "", z)

  out <- rep(NA_integer_, length(z))

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

bad <- !is.finite(ipd$time) | is.na(ipd$status) | is.na(ipd$arm)

if (any(bad)) {
  warning(
    sum(bad),
    " rows have missing/non-numeric time, status, or arm and will be removed."
  )
  ipd <- ipd[!bad, , drop = FALSE]
}

# DO NOT transform negative times with abs().
if (any(ipd$time < 0)) {
  print(head(ipd[ipd$time < 0, , drop = FALSE], 20))
  stop(
    "Negative follow-up times were found. The earlier code used abs(time), ",
    "which changes the data. Correct the source/reconstruction instead."
  )
}

# Preserve reconstructed timing precision: do not round times.
ipd$status <- ifelse(ipd$status > 0, 1L, 0L)

if (!all(c(0L, 1L) %in% unique(ipd$arm))) {
  stop("Both placebo and nirsevimab arms could not be identified.")
}

# Administrative censoring at the endpoint-specific maximum follow-up.
beyond_horizon <- ipd$time > ANALYSIS_T_MAX
if (any(beyond_horizon)) {
  ipd$time[beyond_horizon] <- ANALYSIS_T_MAX
  ipd$status[beyond_horizon] <- 0L
}

cat("\nClean IPD summary:\n")
print(table(Arm = ipd$arm, Event = ipd$status))
cat("N total:", nrow(ipd), "\n")
cat("Maximum follow-up:", max(ipd$time), "\n")
cat("Total events:", sum(ipd$status), "\n")


# ============================================================
# 5. VALIDATE 
# ============================================================

EXPECTED_N_TOTAL <- 2350
EXPECTED_N_PLACEBO <- 786
EXPECTED_N_NIRSEVIMAB <- 1564
EXPECTED_EVENTS <- 30

if (nrow(ipd) != EXPECTED_N_TOTAL) {
  warning(
    " N = 2350. This RDS has ", nrow(ipd), " usable records."
  )
}

if (sum(ipd$arm == 0L) != EXPECTED_N_PLACEBO) {
  warning(
    "placebo N = 786. This RDS has ",
    sum(ipd$arm == 0L), "."
  )
}

if (sum(ipd$arm == 1L) != EXPECTED_N_NIRSEVIMAB) {
  warning(
    "nirsevimab N = 1564. This RDS has ",
    sum(ipd$arm == 1L), "."
  )
}

if (sum(ipd$status) != EXPECTED_EVENTS) {
  warning(
    " total events = 30. This RDS has ",
    sum(ipd$status), "."
  )
}


# ============================================================
# 6. RISK-TABLE VALIDATION
# ============================================================

reported_risk <- data.frame(
  day = c(0, 30, 60, 90, 120, 150),
  placebo_reported = c(786, 778, 769, 761, 757, 753),
  nirsevimab_reported = c(1564, 1554, 1547, 1540, 1535, 1529)
)

reported_risk$placebo_ipd <- sapply(
  reported_risk$day,
  function(tt) sum(ipd$arm == 0L & ipd$time >= tt)
)

reported_risk$nirsevimab_ipd <- sapply(
  reported_risk$day,
  function(tt) sum(ipd$arm == 1L & ipd$time >= tt)
)

reported_risk$placebo_match <-
  reported_risk$placebo_reported == reported_risk$placebo_ipd

reported_risk$nirsevimab_match <-
  reported_risk$nirsevimab_reported == reported_risk$nirsevimab_ipd

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
    "The RDS does not exactly reproduce the reported hospitalisation ",
    "risk table. Review the RDS before calling this a replication."
  )
}


# ============================================================
# 7. BUILD REPORT-SPECIFIC INTERVAL BINOMIAL DATA
#
# Intervals:
#   (0,30], (30,60], (60,90], (90,120], (120,150]
#
# Likelihood:
#   c_p,k ~ Binomial(n_p,k, pi_p,k)
#   c_v,k ~ Binomial(n_v,k, pi_v,k)
#   pi_v,k = pi_p,k * [1 - VE_r(t_k)]
# ============================================================

make_interval_data <- function(dat, breaks) {

  K <- length(breaks) - 1L

  out <- lapply(seq_len(K), function(k) {

    lo <- breaks[k]
    hi <- breaks[k + 1]

    n_p <- sum(dat$arm == 0L & dat$time >= lo)
    n_v <- sum(dat$arm == 1L & dat$time >= lo)

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
    stop("At least one interval has more events than participants at risk.")
  }

  out
}

interval_dat <- make_interval_data(ipd, RISK_TABLE_BREAKS)

cat("\nInterval data used in Bayesian likelihood:\n")
print(interval_dat)

cat(
  "Total interval events:",
  sum(interval_dat$c_p + interval_dat$c_v),
  "\n"
)

if (sum(interval_dat$c_p + interval_dat$c_v) != sum(ipd$status)) {
  warning(
    "Interval-event total differs from the IPD event total. ",
    "Check events exactly at interval boundaries or day 0."
  )
}

write.csv(
  interval_dat,
  file.path(OUT_DIR, "bayesian_interval_data.csv"),
  row.names = FALSE
)


# ============================================================
# 8. STAN MODEL DEFINITIONS
# ============================================================

stan_common <- '
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
stan_common,
'
data {
  real<lower=0> sigma_kappa;
}
parameters {
  vector<lower=0,upper=1>[K] pi_p;
  real<lower=0,upper=1> VE0;
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
  vector[2*K] log_lik;

  for (k in 1:K) {
    real VE_k = VE0 * exp(-kappa * t_mid[k]);
    real pi_v = pi_p[k] * (1 - VE_k);

    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k], pi_p[k]);
    log_lik[K+k] = binomial_lpmf(c_v[k] | n_v[k], pi_v);
  }
}
'
)

stan_erlang3 <- paste0(
stan_common,
'
parameters {
  vector<lower=0,upper=1>[K] pi_p;
  real<lower=0,upper=1> VE0;
  real<lower=0> beta_Gamma;
}
model {
  VE0 ~ beta(a0, b0);
  beta_Gamma ~ gamma(', BETA_GAMMA_SHAPE, ', ', BETA_GAMMA_RATE, ');
  pi_p ~ beta(baseline_alpha, baseline_beta);

  for (k in 1:K) {
    real VE_k = VE0 * (1 - gamma_cdf(t_mid[k] | 3, beta_Gamma));
    real pi_v = pi_p[k] * (1 - VE_k);

    c_p[k] ~ binomial(n_p[k], pi_p[k]);
    c_v[k] ~ binomial(n_v[k], pi_v);
  }
}
generated quantities {
  vector[2*K] log_lik;

  for (k in 1:K) {
    real VE_k = VE0 * (1 - gamma_cdf(t_mid[k] | 3, beta_Gamma));
    real pi_v = pi_p[k] * (1 - VE_k);

    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k], pi_p[k]);
    log_lik[K+k] = binomial_lpmf(c_v[k] | n_v[k], pi_v);
  }
}
'
)

stan_powerlaw <- paste0(
stan_common,
'
data {
  real<lower=0> sigma_kappa;
  real<lower=0> rho_shape;
  real<lower=0> rho_rate;
}
parameters {
  vector<lower=0,upper=1>[K] pi_p;
  real<lower=0,upper=1> VE0;
  real<lower=0> kappa;
  real<lower=0> rho;
}
model {
  VE0 ~ beta(a0, b0);
  kappa ~ normal(0, sigma_kappa);
  rho ~ gamma(rho_shape, rho_rate);
  pi_p ~ beta(baseline_alpha, baseline_beta);

  for (k in 1:K) {
    real VE_k = VE0 / (1 + pow(kappa * t_mid[k], rho));
    real pi_v = pi_p[k] * (1 - VE_k);

    c_p[k] ~ binomial(n_p[k], pi_p[k]);
    c_v[k] ~ binomial(n_v[k], pi_v);
  }
}
generated quantities {
  vector[2*K] log_lik;

  for (k in 1:K) {
    real VE_k = VE0 / (1 + pow(kappa * t_mid[k], rho));
    real pi_v = pi_p[k] * (1 - VE_k);

    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k], pi_p[k]);
    log_lik[K+k] = binomial_lpmf(c_v[k] | n_v[k], pi_v);
  }
}
'
)

stan_files <- c(
  exponential = file.path(OUT_DIR, "hospitalisation_exponential.stan"),
  erlang3 = file.path(OUT_DIR, "hospitalisation_erlang3.stan"),
  powerlaw = file.path(OUT_DIR, "hospitalisation_powerlaw.stan")
)

writeLines(stan_exponential, stan_files["exponential"])
writeLines(stan_erlang3, stan_files["erlang3"])
writeLines(stan_powerlaw, stan_files["powerlaw"])


# ============================================================
# 9. STAN DATA
# ============================================================

stan_data_base <- list(
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
  stan_data_base,
  list(sigma_kappa = SIGMA_KAPPA)
)

stan_data_erl <- stan_data_base

stan_data_pow <- c(
  stan_data_base,
  list(
    sigma_kappa = SIGMA_KAPPA,
    rho_shape = RHO_SHAPE,
    rho_rate = RHO_RATE
  )
)


# ============================================================
# 10. COMPILE AND FIT WITH HMC/NUTS
# ============================================================

cat("\nCompiling Stan models...\n")

mod_exp <- cmdstan_model(stan_files["exponential"])
mod_erl <- cmdstan_model(stan_files["erlang3"])
mod_pow <- cmdstan_model(stan_files["powerlaw"])

sample_model <- function(model, data, seed, label) {

  cat("\n============================================================\n")
  cat("Sampling:", label, "\n")
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
  mod_exp, stan_data_exp, SEED + 1, "Exponential"
)

fit_erl <- sample_model(
  mod_erl, stan_data_erl, SEED + 2, "Erlang-3"
)

fit_pow <- sample_model(
  mod_pow, stan_data_pow, SEED + 3, "Power-law"
)

fits <- list(
  exponential = fit_exp,
  erlang3 = fit_erl,
  powerlaw = fit_pow
)


# ============================================================
# 11. CONVERGENCE DIAGNOSTICS
# ============================================================

diagnose_fit <- function(fit, label, parameters) {

  s <- fit$summary(variables = parameters)
  d <- fit$diagnostic_summary()

  data.frame(
    model = label,
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_bulk_ess = min(s$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(s$ess_tail, na.rm = TRUE),
    divergences = sum(d$num_divergent),
    max_treedepth_hits = sum(d$num_max_treedepth),
    passes_report_rule =
      max(s$rhat, na.rm = TRUE) < 1.01 &&
      min(s$ess_bulk, na.rm = TRUE) > 400
  )
}

convergence <- bind_rows(
  diagnose_fit(fit_exp, "Exponential", c("VE0", "kappa")),
  diagnose_fit(fit_erl, "Erlang-3", c("VE0", "beta_Gamma")),
  diagnose_fit(fit_pow, "Power-law", c("VE0", "kappa", "rho"))
)

cat("\nConvergence diagnostics:\n")
print(convergence)

write.csv(
  convergence,
  file.path(OUT_DIR, "convergence_diagnostics.csv"),
  row.names = FALSE
)


# ============================================================
# 12. FOUR-CHAIN TRACE PLOTS
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
  "Nirsevimab RSV LRTI hospital admission - Exponential",
  "trace_exponential.png"
)

save_trace_plot(
  fit_erl,
  c("VE0", "beta_Gamma"),
  "Nirsevimab RSV LRTI hospital admission - Erlang-3",
  "trace_erlang3.png"
)

save_trace_plot(
  fit_pow,
  c("VE0", "kappa", "rho"),
  "Nirsevimab RSV LRTI hospital admission - Power-law",
  "trace_powerlaw.png"
)


# ============================================================
# 13. LOOIC MODEL SELECTION
# ============================================================

get_loo <- function(fit) {

  log_lik <- fit$draws(
    variables = "log_lik",
    format = "matrix"
  )

  loo(log_lik)
}

loo_exp <- get_loo(fit_exp)
loo_erl <- get_loo(fit_erl)
loo_pow <- get_loo(fit_pow)

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

selected_model <- loo_table$model[1]
selected_fit <- fits[[selected_model]]

cat("\nSelected waning form:", selected_model, "\n")

if (selected_model != "exponential") {
  warning(
    "The exponential/single-stage decay had the ",
    "lowest LOOIC for this endpoint. This run selected ", selected_model, ".\n"
  )
}


# ============================================================
# 14. PARETO-k DIAGNOSTICS
# ============================================================

pareto_table <- bind_rows(
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
print(pareto_table)

write.csv(
  pareto_table,
  file.path(OUT_DIR, "LOO_pareto_k_diagnostics.csv"),
  row.names = FALSE
)


# ============================================================
# 15. POSTERIOR VE(t)
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
    stop("Unknown waning form: ", form)
  }

  if (is.null(dim(ve_draws))) {
    ve_draws <- matrix(ve_draws, ncol = length(time_grid))
  }

  data.frame(
    day = time_grid,
    posterior_median = apply(ve_draws, 2, median),
    lower_95_CrI = apply(
      ve_draws, 2, quantile, probs = 0.025, names = FALSE
    ),
    upper_95_CrI = apply(
      ve_draws, 2, quantile, probs = 0.975, names = FALSE
    ),
    posterior_mean = colMeans(ve_draws)
  )
}

ve_selected <- posterior_ve(
  selected_fit,
  selected_model,
  OUTPUT_GRID
)

cat("\nSelected-model VE estimates:\n")
print(ve_selected)

write.csv(
  ve_selected,
  file.path(OUT_DIR, "Bayesian_VE_15_day_grid.csv"),
  row.names = FALSE
)


# ============================================================
# 16. VALIDATE AGAINST TABLE 5.1
# ============================================================

report_table_5_1 <- data.frame(
  day = seq(0, 150, by = 15),
  report_median = c(
    0.78, 0.76, 0.74, 0.72, 0.70, 0.68,
    0.66, 0.64, 0.62, 0.61, 0.60
  ),
  report_lower = c(
    0.64, 0.63, 0.61, 0.59, 0.57, 0.54,
    0.51, 0.48, 0.44, 0.41, 0.39
  ),
  report_upper = c(
    0.86, 0.84, 0.82, 0.80, 0.78, 0.77,
    0.76, 0.75, 0.74, 0.73, 0.73
  )
)

validation <- left_join(
  ve_selected,
  report_table_5_1,
  by = "day"
) |>
  mutate(
    median_difference = posterior_median - report_median,
    lower_difference = lower_95_CrI - report_lower,
    upper_difference = upper_95_CrI - report_upper
  )

cat("\nComparison with Table 5.1:\n")
print(validation)

write.csv(
  validation,
  file.path(OUT_DIR, "comparison_with_Table_5_1.csv"),
  row.names = FALSE
)


# ============================================================
# 17. VE CURVE PLOT
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
      "Bayesian VE(t): Nirsevimab against RSV LRTI hospital admission\n",
      "LOOIC-selected form: ", selected_model
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
# 18. SELECTED-MODEL PARAMETER SUMMARY
# ============================================================

selected_parameters <- switch(
  selected_model,
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
# 19. PRIOR SENSITIVITY PLAN
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
# 20. SAVE ANALYSIS SUMMARY
# ============================================================

analysis_summary <- list(
  endpoint = "Nirsevimab - RSV LRTI hospital admission",
  data_file = DATA_FILE,
  horizon_days = ANALYSIS_T_MAX,
  risk_table_breaks = RISK_TABLE_BREAKS,
  output_grid = OUTPUT_GRID,

  final_report_targets = list(
    total_N = EXPECTED_N_TOTAL,
    placebo_N = EXPECTED_N_PLACEBO,
    nirsevimab_N = EXPECTED_N_NIRSEVIMAB,
    total_events = EXPECTED_EVENTS,
    selected_waning_form = "exponential",
    Bayesian_day0 = c(median = 0.78, lower = 0.64, upper = 0.86),
    Bayesian_day150 = c(median = 0.60, lower = 0.39, upper = 0.73)
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

  risk_table_validation = reported_risk,
  interval_data = interval_dat,
  convergence = convergence,
  looic = loo_table,
  pareto_k = pareto_table,
  selected_model = selected_model,
  VE = ve_selected,
  report_validation = validation
)

saveRDS(
  analysis_summary,
  file.path(
    OUT_DIR,
    "Bayesian_RSV_Hospitalization_analysis_summary.rds"
  )
)

writeLines(
  c(
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
  ),
  file.path(OUT_DIR, "cmdstan_output_file_locations.txt")
)


# ============================================================
# 21. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("Endpoint: Nirsevimab - RSV LRTI hospital admission\n")
cat("Follow-up: 0-150 days\n")
cat("N in cleaned IPD:", nrow(ipd), "\n")
cat("Events in cleaned IPD:", sum(ipd$status), "\n")
cat("LOOIC-selected model:", selected_model, "\n")
cat("Final expected model: exponential\n")
cat("Output directory:", OUT_DIR, "\n")
cat("============================================================\n")
