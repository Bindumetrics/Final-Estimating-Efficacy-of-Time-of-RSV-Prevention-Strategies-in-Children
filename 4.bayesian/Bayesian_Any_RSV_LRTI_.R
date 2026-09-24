# ============================================================
# BAYESIAN WANING-EFFICACY ANALYSIS 
# Endpoint: Maternal RSVpreF vaccine - medically attended RSV-LRTI (all cases)
# Data file: Any_RSV_LRTI.rds
#
# This script is reworked to match Chapter 4.6 of the final report:
#   * risk-scale VE
#   * interval-level BINOMIAL likelihood
#   * pi_vk = pi_pk * (1 - VE_r,k)
#   * candidate waning functions:
#       1) exponential
#       2) Erlang-3
#       3) power-law
#   * model selection by LOOIC
#   * HMC/NUTS through cmdstanr
#   * four chains
#   * 20,000 iterations per chain, 1,000 warm-up
#   * convergence checks: R-hat < 1.01 and bulk ESS > 400
#   * pointwise posterior median and equal-tailed 95% CrI
#   * VE evaluated every 15 days through day 360
#
# IMPORTANT IMPLEMENTATION NOTES
# -------------------------------------
#
#   VE0 ~ Beta(a0,b0)
#   kappa ~ Half-Normal(0,sigma_kappa^2)
#   beta_Gamma ~ Gamma(2,1)
#
# We do not explicitly state a prior for the interval-specific baseline
# risks pi_pk or for rho in the power-law model. Therefore:
#   * pi_pk is given Beta(1,1), i.e. the explicit uniform prior implied by an
#     otherwise unconstrained probability parameter.
#   * rho is given Gamma(2,1) as a weak positive prior.
#   * VE0 prior is centred at the earliest published maternal all-cases VE
#     reported in the thesis (57.1% at 90 days), with a weak concentration of 4.
# These choices are clearly isolated below so they can be replaced if the
# original exact hyperparameter values are available.
#
# The maternal-vaccine endpoints selected Erlang-3 by
# LOOIC. This script nevertheless fits all three forms and selects by LOOIC,
# matching the stated model-selection procedure.
#
# 
# ============================================================


# ============================================================
# 0. USER SETTINGS
# ============================================================

DATA_DIR  <- "D:/Desktop/Tiwonge/Tiwonge"
DATA_FILE <- file.path(DATA_DIR, "Any_RSV_LRTI.rds")
OUT_DIR   <- file.path(DATA_DIR, "analysis_outputs_Any_RSV_LRTI_360_DAYS_BAYES_REPORT_ALIGNED")

ANALYSIS_T_MAX <- 360

# Interval width for the likelihood.
# The thesis reconstruction uses 30-day risk-table intervals; VE itself is
# evaluated every 15 days.
LIKELIHOOD_INTERVAL_DAYS <- 30
OUTPUT_GRID_DAYS <- 15

# Final-report sampler settings
CHAINS <- 4
ITER_TOTAL <- 20000
ITER_WARMUP <- 1000
ITER_SAMPLING <- ITER_TOTAL - ITER_WARMUP
ADAPT_DELTA <- 0.95
MAX_TREEDEPTH <- 12
SEED <- 2024

# ---------------------------
# Prior hyperparameters
# ---------------------------

# VE0 ~ Beta(a0,b0), prior mean aligned with earliest published estimate.
VE0_PRIOR_MEAN <- 0.571
VE0_PRIOR_CONCENTRATION <- 4
A0 <- VE0_PRIOR_MEAN * VE0_PRIOR_CONCENTRATION
B0 <- (1 - VE0_PRIOR_MEAN) * VE0_PRIOR_CONCENTRATION

# kappa ~ Half-Normal(0, sigma_kappa^2), chosen so the implied
# half-life is not implausibly longer than one RSV season.
# 
SIGMA_KAPPA <- 0.03

#  beta_Gamma ~ Gamma(2,1)
BETA_GAMMA_SHAPE <- 2
BETA_GAMMA_RATE  <- 1

# Power-law rho prior.
RHO_SHAPE <- 2
RHO_RATE  <- 1


# ============================================================
# 1. PACKAGES AND OUTPUT DIRECTORY
# ============================================================

required <- c("cmdstanr", "posterior", "loo", "ggplot2", "dplyr", "tidyr")

for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(
      "Package '", p, "' is required but is not installed.\n",
      "Install it before rerunning this script."
    )
  }
}

library(cmdstanr)
library(posterior)
library(loo)
library(ggplot2)
library(dplyr)
library(tidyr)

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

cat("Data file:", DATA_FILE, "\n")
cat("Output directory:", OUT_DIR, "\n")


# ============================================================
# 2. LOAD THE RDS AND EXTRACT THE RECONSTRUCTED / EXTENDED IPD
# ============================================================

if (!file.exists(DATA_FILE)) stop("Data file not found: ", DATA_FILE)

obj <- readRDS(DATA_FILE)

cat("\nTop-level RDS elements:\n")
print(names(obj))

if (is.null(obj$ipd)) {
  stop(
    "The final-report Bayesian model should be built from reconstructed IPD.\n",
    "No obj$ipd element was found in Any_RSV_LRTI.rds."
  )
}

ipd_raw <- as.data.frame(obj$ipd)


# ============================================================
# 3. ROBUST VARIABLE DETECTION
# ============================================================

find_col <- function(df, candidates, label) {
  nm <- names(df)
  low <- tolower(nm)

  # exact
  hit <- match(tolower(candidates), low, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit) > 0) return(nm[hit[1]])

  # partial
  for (cc in tolower(candidates)) {
    h <- grep(cc, low, fixed = TRUE)
    if (length(h) > 0) return(nm[h[1]])
  }

  stop(
    "Could not identify ", label, ".\n",
    "Tried: ", paste(candidates, collapse = ", "), "\n",
    "Available columns: ", paste(nm, collapse = ", ")
  )
}

time_col <- find_col(
  ipd_raw,
  c("time", "t", "day", "days", "followup", "follow_up", "ftime"),
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


standardize_arm <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z2 <- gsub("[_ -]+", "", z)

  # Handle common text labels
  out <- ifelse(
    grepl("plac|control|standard", z2),
    0L,
    ifelse(grepl("maternal|vacc|active|treat|rsvpref", z2), 1L, NA_integer_)
  )

  # Handle 0/1 numeric arm coding if present
  num <- suppressWarnings(as.numeric(z))
  out[is.na(out) & num == 0] <- 0L
  out[is.na(out) & num == 1] <- 1L

  out
}

ipd <- data.frame(
  time = suppressWarnings(as.numeric(ipd_raw[[time_col]])),
  status = suppressWarnings(as.numeric(ipd_raw[[event_col]])),
  arm = standardize_arm(ipd_raw[[arm_col]])
)

# Keep valid records only
ipd <- ipd[
  is.finite(ipd$time) &
  !is.na(ipd$status) &
  !is.na(ipd$arm),
]

# Do NOT use abs(time). Negative follow-up times indicate a data problem.
if (any(ipd$time < 0)) {
  bad <- head(ipd[ipd$time < 0, ], 10)
  print(bad)
  stop(
    "Negative follow-up times were found. The earlier script used abs(time), ",
    "which changes the data and is not appropriate. Correct the source data first."
  )
}

ipd$status <- ifelse(ipd$status > 0, 1L, 0L)

if (!all(c(0L, 1L) %in% unique(ipd$arm))) {
  stop("Both placebo/control (0) and maternal-vaccine (1) arms were not identified.")
}

cat("\nIPD summary:\n")
print(table(ipd$arm, ipd$status))
cat("Observed maximum follow-up:", max(ipd$time), "days\n")

# The final report's maternal all-cases analysis is 360 days.
#
if (max(ipd$time) < ANALYSIS_T_MAX - 1e-8) {
  stop(
    "This RDS only reaches ", max(ipd$time), " days, but the report's all-cases ",
    "analysis uses a 360-day dataset with a Weibull-extended 180-360 day segment.\n",
    "Use the already-extended 360-day RDS or perform the Chapter 3 Weibull extension first."
  )
}

# Administrative restriction to the report horizon
# Anyone whose original reconstructed/extended time is beyond day 360 is
# censored administratively at day 360.
beyond_horizon <- ipd$time > ANALYSIS_T_MAX
ipd$time <- pmin(ipd$time, ANALYSIS_T_MAX)
ipd$status[beyond_horizon] <- 0L


# ============================================================
# 4. CONSTRUCT INTERVAL-LEVEL BINOMIAL DATA
#
# Report likelihood:
#   c_vk ~ Binomial(n_vk, pi_vk)
#   c_pk ~ Binomial(n_pk, pi_pk)
#   pi_vk = pi_pk * (1 - VE_r,k)
#
# n_*k is the number at risk at the beginning of interval k.
# c_*k is the number of events during interval k.
# ============================================================

make_interval_data <- function(dat, horizon = 360, width = 30) {

  breaks <- seq(0, horizon, by = width)
  if (tail(breaks, 1) < horizon) breaks <- c(breaks, horizon)

  K <- length(breaks) - 1L

  rows <- lapply(seq_len(K), function(k) {

    lo <- breaks[k]
    hi <- breaks[k + 1]

    # At risk immediately after lo.
    # Events/censorings exactly at lo belong to the previous closed-right interval,
    # except at baseline.
    at_risk <- if (lo == 0) {
      dat$time >= lo
    } else {
      dat$time > lo
    }

    # Interval convention: (lo, hi], with day 0 included in first interval if needed.
    in_event <- dat$status == 1L &
      dat$time > lo &
      dat$time <= hi

    data.frame(
      k = k,
      lo = lo,
      hi = hi,
      t_mid = (lo + hi) / 2,
      n_p = sum(at_risk & dat$arm == 0L),
      c_p = sum(in_event & dat$arm == 0L),
      n_v = sum(at_risk & dat$arm == 1L),
      c_v = sum(in_event & dat$arm == 1L)
    )
  })

  out <- bind_rows(rows)

  if (any(out$c_p > out$n_p) || any(out$c_v > out$n_v)) {
    stop("At least one interval has event count greater than number at risk.")
  }

  out
}

interval_dat <- make_interval_data(
  ipd,
  horizon = ANALYSIS_T_MAX,
  width = LIKELIHOOD_INTERVAL_DAYS
)

cat("\nInterval-level data used by the Bayesian likelihood:\n")
print(interval_dat)

write.csv(
  interval_dat,
  file.path(OUT_DIR, "interval_binomial_data.csv"),
  row.names = FALSE
)


# ============================================================
# 5. WRITE THREE STAN MODELS
# ============================================================

# All three models estimate interval-specific placebo baseline risk pi_p[k].
# This makes the Bayesian model operate on the risk scale exactly as described
# in Eq. 4.20 of the report.

common_data_block <- '
data {
  int<lower=1> K;
  array[K] int<lower=0> n_p;
  array[K] int<lower=0> c_p;
  array[K] int<lower=0> n_v;
  array[K] int<lower=0> c_v;
  vector<lower=0>[K] t_mid;

  real<lower=0> a0;
  real<lower=0> b0;
}
'

stan_exp <- paste0(
common_data_block,
'
parameters {
  vector<lower=0,upper=1>[K] pi_p;
  real<lower=0,upper=1> VE0;
  real<lower=0> kappa;
}
model {
  VE0 ~ beta(a0, b0);
  kappa ~ normal(0, ', SIGMA_KAPPA, ');  // half-normal due to lower bound
  pi_p ~ beta(1, 1);

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

# In Stan gamma_cdf(x | alpha, beta), beta is a RATE.
# This follows the report notation literally: 1 - F_Gamma(t; 3, beta_Gamma).
stan_erlang <- paste0(
common_data_block,
'
parameters {
  vector<lower=0,upper=1>[K] pi_p;
  real<lower=0,upper=1> VE0;
  real<lower=0> beta_Gamma;
}
model {
  VE0 ~ beta(a0, b0);
  beta_Gamma ~ gamma(', BETA_GAMMA_SHAPE, ', ', BETA_GAMMA_RATE, ');
  pi_p ~ beta(1, 1);

  for (k in 1:K) {
    real surv3 = 1 - gamma_cdf(t_mid[k] | 3, beta_Gamma);
    real VE_k = VE0 * surv3;
    real pi_v = pi_p[k] * (1 - VE_k);

    c_p[k] ~ binomial(n_p[k], pi_p[k]);
    c_v[k] ~ binomial(n_v[k], pi_v);
  }
}
generated quantities {
  vector[2*K] log_lik;
  for (k in 1:K) {
    real surv3 = 1 - gamma_cdf(t_mid[k] | 3, beta_Gamma);
    real VE_k = VE0 * surv3;
    real pi_v = pi_p[k] * (1 - VE_k);

    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k], pi_p[k]);
    log_lik[K+k] = binomial_lpmf(c_v[k] | n_v[k], pi_v);
  }
}
'
)

stan_power <- paste0(
common_data_block,
'
parameters {
  vector<lower=0,upper=1>[K] pi_p;
  real<lower=0,upper=1> VE0;
  real<lower=0> kappa;
  real<lower=0> rho;
}
model {
  VE0 ~ beta(a0, b0);
  kappa ~ normal(0, ', SIGMA_KAPPA, ');  // half-normal
  rho ~ gamma(', RHO_SHAPE, ', ', RHO_RATE, ');
  pi_p ~ beta(1, 1);

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

stan_paths <- c(
  exponential = file.path(OUT_DIR, "waning_exponential.stan"),
  erlang3     = file.path(OUT_DIR, "waning_erlang3.stan"),
  powerlaw    = file.path(OUT_DIR, "waning_powerlaw.stan")
)

writeLines(stan_exp, stan_paths["exponential"])
writeLines(stan_erlang, stan_paths["erlang3"])
writeLines(stan_power, stan_paths["powerlaw"])


# ============================================================
# 6. STAN DATA
# ============================================================

stan_data <- list(
  K = nrow(interval_dat),
  n_p = as.integer(interval_dat$n_p),
  c_p = as.integer(interval_dat$c_p),
  n_v = as.integer(interval_dat$n_v),
  c_v = as.integer(interval_dat$c_v),
  t_mid = as.numeric(interval_dat$t_mid),
  a0 = A0,
  b0 = B0
)


# ============================================================
# 7. COMPILE AND SAMPLE USING HMC/NUTS
# ============================================================

compile_and_fit <- function(path, model_name, seed_offset = 0) {

  cat("\n============================================================\n")
  cat("Fitting:", model_name, "\n")
  cat("============================================================\n")

  mod <- cmdstan_model(path)

  fit <- mod$sample(
    data = stan_data,
    seed = SEED + seed_offset,
    chains = CHAINS,
    parallel_chains = CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    adapt_delta = ADAPT_DELTA,
    max_treedepth = MAX_TREEDEPTH,
    refresh = 500
  )

  fit
}

fit_exp <- compile_and_fit(stan_paths["exponential"], "Exponential", 1)
fit_erl <- compile_and_fit(stan_paths["erlang3"], "Erlang-3", 2)
fit_pow <- compile_and_fit(stan_paths["powerlaw"], "Power-law", 3)

fits <- list(
  exponential = fit_exp,
  erlang3 = fit_erl,
  powerlaw = fit_pow
)


# ============================================================
# 8. CONVERGENCE DIAGNOSTICS
# ============================================================

diagnose_fit <- function(fit, name) {

  sm <- fit$summary()

  # Exclude generated log_lik from the headline parameter diagnostic
  pars <- sm[!grepl("^log_lik", sm$variable), ]

  rhat_max <- max(pars$rhat, na.rm = TRUE)
  ess_bulk_min <- min(pars$ess_bulk, na.rm = TRUE)
  ess_tail_min <- min(pars$ess_tail, na.rm = TRUE)

  diag <- fit$diagnostic_summary()
  divergences <- sum(diag$num_divergent)
  max_td <- sum(diag$num_max_treedepth)

  data.frame(
    model = name,
    rhat_max = rhat_max,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    divergences = divergences,
    max_treedepth_hits = max_td,
    report_convergence_rule =
      rhat_max < 1.01 && ess_bulk_min > 400
  )
}

conv <- bind_rows(
  diagnose_fit(fit_exp, "Exponential"),
  diagnose_fit(fit_erl, "Erlang-3"),
  diagnose_fit(fit_pow, "Power-law")
)

cat("\nConvergence diagnostics:\n")
print(conv)

write.csv(
  conv,
  file.path(OUT_DIR, "convergence_diagnostics.csv"),
  row.names = FALSE
)


# ============================================================
# 9. TRACE PLOTS
# ============================================================

save_trace <- function(fit, pars, filename, title) {

  draws <- fit$draws(variables = pars, format = "df")

  long <- draws |>
    as.data.frame() |>
    select(.chain, .iteration, all_of(pars)) |>
    pivot_longer(
      cols = all_of(pars),
      names_to = "parameter",
      values_to = "value"
    )

  p <- ggplot(
    long,
    aes(x = .iteration, y = value, group = .chain, colour = factor(.chain))
  ) +
    geom_line(alpha = 0.65, linewidth = 0.25) +
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
    height = max(5, 2.3 * length(pars)),
    dpi = 300
  )
}

save_trace(
  fit_exp,
  c("VE0", "kappa"),
  "trace_exponential.png",
  "Bayesian trace plots - Exponential"
)

save_trace(
  fit_erl,
  c("VE0", "beta_Gamma"),
  "trace_erlang3.png",
  "Bayesian trace plots - Erlang-3"
)

save_trace(
  fit_pow,
  c("VE0", "kappa", "rho"),
  "trace_powerlaw.png",
  "Bayesian trace plots - Power-law"
)


# ============================================================
# 10. LOOIC MODEL SELECTION
# ============================================================

get_loo <- function(fit) {
  ll <- fit$draws("log_lik", format = "matrix")
  loo(ll)
}

loo_exp <- get_loo(fit_exp)
loo_erl <- get_loo(fit_erl)
loo_pow <- get_loo(fit_pow)

loos <- list(
  exponential = loo_exp,
  erlang3 = loo_erl,
  powerlaw = loo_pow
)

loo_compare_table <- loo_compare(loos)

cat("\nLOO comparison (higher elpd_loo is better):\n")
print(loo_compare_table)

# Convert to a simple table including LOOIC = -2 * elpd_loo
loo_summary <- bind_rows(lapply(names(loos), function(nm) {
  x <- loos[[nm]]
  data.frame(
    model = nm,
    elpd_loo = x$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = x$estimates["elpd_loo", "SE"],
    p_loo = x$estimates["p_loo", "Estimate"],
    looic = -2 * x$estimates["elpd_loo", "Estimate"]
  )
})) |>
  arrange(looic)

cat("\nLOOIC summary (lower is better):\n")
print(loo_summary)

write.csv(
  loo_summary,
  file.path(OUT_DIR, "LOOIC_model_selection.csv"),
  row.names = FALSE
)

selected_name <- loo_summary$model[1]
selected_fit <- fits[[selected_name]]

cat("\nSelected waning form by LOOIC:", selected_name, "\n")

# The maternal-vaccine endpoints selected Erlang-3.
if (selected_name != "erlang3") {
  warning(
    "This run did not select Erlang-3"
  )
}


# ============================================================
# 11. POSTERIOR VE(t) ON THE REPORT'S 15-DAY GRID
# ============================================================

grid <- seq(0, ANALYSIS_T_MAX, by = OUTPUT_GRID_DAYS)

posterior_ve <- function(fit, form, grid) {

  dr <- as_draws_df(fit$draws())

  if (form == "exponential") {

    ve_mat <- sapply(grid, function(t) {
      dr$VE0 * exp(-dr$kappa * t)
    })

  } else if (form == "erlang3") {

    ve_mat <- sapply(grid, function(t) {
      # R pgamma uses rate= to match Stan gamma_cdf(..., beta_Gamma)
      dr$VE0 * (1 - pgamma(t, shape = 3, rate = dr$beta_Gamma))
    })

  } else if (form == "powerlaw") {

    ve_mat <- sapply(grid, function(t) {
      dr$VE0 / (1 + (dr$kappa * t)^dr$rho)
    })

  } else {
    stop("Unknown form")
  }

  # draws x grid
  if (is.null(dim(ve_mat))) ve_mat <- matrix(ve_mat, ncol = length(grid))

  out <- data.frame(
    day = grid,
    median = apply(ve_mat, 2, median),
    lower_95 = apply(ve_mat, 2, quantile, probs = 0.025),
    upper_95 = apply(ve_mat, 2, quantile, probs = 0.975),
    mean = colMeans(ve_mat)
  )

  out
}

ve_selected <- posterior_ve(selected_fit, selected_name, grid)

write.csv(
  ve_selected,
  file.path(OUT_DIR, "Bayesian_VE_selected_model_15day_grid.csv"),
  row.names = FALSE
)

cat("\nSelected-model VE summary:\n")
print(ve_selected)


# ============================================================
# 12. PLOT SELECTED BAYESIAN WANING CURVE
# ============================================================

p_ve <- ggplot(ve_selected, aes(x = day)) +
  geom_ribbon(
    aes(ymin = lower_95, ymax = upper_95),
    alpha = 0.20
  ) +
  geom_line(
    aes(y = median),
    linewidth = 0.9
  ) +
  geom_vline(
    xintercept = 180,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  coord_cartesian(
    xlim = c(0, ANALYSIS_T_MAX),
    ylim = c(0, 1)
  ) +
  labs(
    title = paste0(
      "Bayesian time-varying efficacy - maternal RSVpreF (",
      selected_name, ")"
    ),
    subtitle = "Dashed line: boundary between reconstructed 0-180 d and Weibull-extended 180-360 d segment",
    x = "Days after birth",
    y = "Risk-scale vaccine efficacy",
    caption = "Solid line: posterior median; shaded band: equal-tailed 95% credible interval"
  ) +
  theme_bw()

ggsave(
  file.path(OUT_DIR, "Bayesian_waning_selected_model.png"),
  p_ve,
  width = 9,
  height = 6,
  dpi = 300
)


# ============================================================
# 13. REPORT-RELEVANT SUMMARY DAYS
# ============================================================

days_to_report <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135,
                    150, 165, 180, 195, 210, 225, 240, 255, 270,
                    285, 300, 315, 330, 345, 360)

report_table <- ve_selected |>
  filter(day %in% days_to_report) |>
  mutate(
    VE_CrI = sprintf(
      "%.2f (%.2f, %.2f)",
      median, lower_95, upper_95
    )
  )

write.csv(
  report_table,
  file.path(OUT_DIR, "Bayesian_report_table_15day.csv"),
  row.names = FALSE
)

cat("\n15-day report table:\n")
print(report_table)


# ============================================================
# 14. PARAMETER SUMMARIES
# ============================================================

selected_parameters <- switch(
  selected_name,
  exponential = c("VE0", "kappa"),
  erlang3 = c("VE0", "beta_Gamma"),
  powerlaw = c("VE0", "kappa", "rho")
)

param_summary <- selected_fit$summary(selected_parameters) |>
  select(variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

write.csv(
  param_summary,
  file.path(OUT_DIR, "selected_model_parameter_summary.csv"),
  row.names = FALSE
)

cat("\nSelected-model parameter summary:\n")
print(param_summary)


# ============================================================
# 15. SAVE ALL FIT OBJECTS AND METADATA
# ============================================================

saveRDS(
  list(
    interval_data = interval_dat,
    prior_settings = list(
      a0 = A0,
      b0 = B0,
      VE0_prior_mean = VE0_PRIOR_MEAN,
      VE0_prior_concentration = VE0_PRIOR_CONCENTRATION,
      sigma_kappa = SIGMA_KAPPA,
      beta_Gamma_shape = BETA_GAMMA_SHAPE,
      beta_Gamma_rate = BETA_GAMMA_RATE,
      rho_shape = RHO_SHAPE,
      rho_rate = RHO_RATE
    ),
    loo_summary = loo_summary,
    convergence = conv,
    selected_model = selected_name,
    VE = ve_selected,
    report_table = report_table
  ),
  file.path(OUT_DIR, "Bayesian_analysis_summary.rds")
)

# Save CmdStan CSV file locations rather than attempting to embed the external
# CmdStan processes in a conventional RData object.
writeLines(
  c(
    paste("Exponential CSV:", paste(fit_exp$output_files(), collapse = "; ")),
    paste("Erlang-3 CSV:", paste(fit_erl$output_files(), collapse = "; ")),
    paste("Power-law CSV:", paste(fit_pow$output_files(), collapse = "; "))
  ),
  file.path(OUT_DIR, "cmdstan_output_files.txt")
)

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("Selected model:", selected_name, "\n")
cat("Outputs saved to:", OUT_DIR, "\n")
cat("============================================================\n")
