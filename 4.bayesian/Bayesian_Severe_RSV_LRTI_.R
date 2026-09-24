# ============================================================
# Bayesian Waning Efficacy Analysis
# Maternal RSVpreF: Medically Attended Severe RSV-LRTI
#
# Endpoint:
#   Maternal vaccine (RSVpreF), medically attended severe RSV-LRTI
# Time scale:
#   Days since birth
# Follow-up:
#   0-180 days
#
# The methodology implemented:
#   * reconstructed IPD is the analysis source
#   * interval-binomial likelihood on the risk scale
#   * VE_r(t) = 1 - p_v(t)/p_p(t)
#   * candidate waning forms: exponential, Erlang-3, power law
#   * model selection by LOOIC
#   * HMC/NUTS via cmdstanr
#   * 4 chains, 20,000 iterations/chain, 1,000 warm-up
#   * posterior median + equal-tailed 95% CrI every 15 days
#
# Endpoint validation targets :
#   N placebo = 3495; N vaccine = 3480
#   total reconstructed events at day 180 = 81
#   cumulative reconstructed events:
#       day       0  30  60  90 120 150 180
#       vaccine   0   1   4   6  11  16  19
#       placebo   0   9  28  33  47  55  62
#   Bayesian selected Erlang-3 by LOOIC.
#
# ============================================================

# ---------------- USER SETTINGS ----------------
DATA_DIR  <- "D:/Desktop/Tiwonge/Tiwonge"
DATA_FILE <- file.path(DATA_DIR, "Severe_RSV_LRTI.rds")
OUT_DIR   <- file.path(DATA_DIR, "analysis_outputs_Severe_RSV_LRTI_BAYES")

FOLLOWUP_DAYS <- 180L
INTERVAL_DAYS <- 30L
OUTPUT_STEP   <- 15L

CHAINS        <- 4L
ITER_TOTAL    <- 20000L
ITER_WARMUP   <- 1000L
ITER_SAMPLING <- ITER_TOTAL - ITER_WARMUP
SEED          <- 2024L
ADAPT_DELTA   <- 0.95
MAX_TREEDEPTH <- 12L

# ---------------- PRIORS ----------------
# VE0 ~ Beta(a0,b0), with mean aligned to earliest published estimate.
# Earliest published severe endpoint estimate reported in the thesis = 0.818 at 90 d.
# Concentration is NOT reported, so 4 is an explicit weakly-informative completion.
VE0_PRIOR_MEAN <- 0.818
VE0_PRIOR_CONCENTRATION <- 4
A0 <- VE0_PRIOR_MEAN * VE0_PRIOR_CONCENTRATION
B0 <- (1 - VE0_PRIOR_MEAN) * VE0_PRIOR_CONCENTRATION

# kappa ~ Half-Normal(0, sigma_kappa^2)
SIGMA_KAPPA <- 0.03

# beta_Gamma ~ Gamma(2,1).
BETA_GAMMA_SHAPE <- 2
BETA_GAMMA_RATE  <- 1

RHO_SHAPE <- 2
RHO_RATE  <- 1

# Explicit weak prior:
BASELINE_RISK_A <- 1
BASELINE_RISK_B <- 1

# ---------------- SETUP ----------------
needed <- c("cmdstanr", "posterior", "loo", "ggplot2", "dplyr", "tidyr")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required packages first: ", paste(missing, collapse = ", "))
}
if (!cmdstanr::cmdstan_version(error_on_NA = FALSE) |> length()) {
  stop("CmdStan is not installed. Run cmdstanr::install_cmdstan() once, then rerun.")
}

library(cmdstanr)
library(posterior)
library(loo)
library(ggplot2)
library(dplyr)
library(tidyr)

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
save_plot <- function(p, name, width=9, height=6) {
  print(p)
  ggplot2::ggsave(file.path(OUT_DIR, name), p, width=width, height=height, dpi=300)
}

# ---------------- LOAD RECONSTRUCTED IPD ----------------
obj <- readRDS(DATA_FILE)
#if (is.null(obj$ipd)) {
  #stop("The analysis uses reconstructed individual event/censoring times. ",
   #    "This RDS has no $ipd component; do not silently re-create IPD from KM drops.")
#}
ipd0 <- as.data.frame(obj$ipd)

find_col <- function(df, candidates, label) {
  nm <- names(df); lo <- tolower(nm)
  exact <- match(tolower(candidates), lo, nomatch=0L)
  exact <- exact[exact > 0L]
  if (length(exact)) return(nm[exact[1]])
  for (cc in tolower(candidates)) {
    h <- grep(cc, lo, fixed=TRUE)
    if (length(h)) return(nm[h[1]])
  }
  stop("Cannot identify ", label, ". Available columns: ", paste(nm, collapse=", "))
}

time_col  <- find_col(ipd0, c("time","t","day","days","followup","follow_up","fu","ftime"), "time")
event_col <- find_col(ipd0, c("event","status","d","delta","case","outcome"), "event")
arm_col   <- find_col(ipd0, c("arm","group","treatment","treat","trt","vacc","vaccine"), "arm")

standardize_arm <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    ux <- sort(unique(x[!is.na(x)]))
    if (all(ux %in% c(0,1))) return(ifelse(x == 1, 1L, 0L))
  }
  z <- tolower(trimws(as.character(x)))
  z2 <- gsub("[_ -]+", "", z)
  out <- rep(NA_integer_, length(z2))
  out[grepl("plac|control|standard", z2)] <- 0L
  out[grepl("vacc|vaccine|maternal|rsvpref|abrysvo|active|treat", z2)] <- 1L
  if (anyNA(out)) {
    stop("Could not standardize every treatment arm. Values: ",
         paste(unique(as.character(x)), collapse=", "))
  }
  out
}

ipd <- data.frame(
  time   = suppressWarnings(as.numeric(ipd0[[time_col]])),
  status = suppressWarnings(as.numeric(ipd0[[event_col]])),
  trt    = standardize_arm(ipd0[[arm_col]])
)

if (any(!is.finite(ipd$time))) stop("Non-finite follow-up times detected.")
if (any(ipd$time < 0)) stop("Negative follow-up times detected. Do NOT repair them with abs(); fix upstream reconstruction.")
ipd$status <- ifelse(ipd$status > 0, 1L, 0L)

# Administrative censoring at the endpoint's published 180-day horizon.
after_end <- ipd$time > FOLLOWUP_DAYS
ipd$time[after_end] <- FOLLOWUP_DAYS
ipd$status[after_end] <- 0L

# ---------------- ENDPOINT VALIDATION ----------------
expected_n <- c(placebo=3495L, vaccine=3480L)
observed_n <- c(placebo=sum(ipd$trt==0), vaccine=sum(ipd$trt==1))
cat("\nSample sizes from RDS:\n"); print(observed_n)
if (!all(observed_n == expected_n)) {
  warning("Sample sizes differ from thesis validation table. Expected P/T = 3495/3480.")
}

expected_days <- c(0,30,60,90,120,150,180)
expected_cum_v <- c(0,1,4,6,11,16,19)
expected_cum_p <- c(0,9,28,33,47,55,62)

cum_events <- function(trt_value, day) {
  sum(ipd$trt == trt_value & ipd$status == 1 & ipd$time <= day)
}
validation <- data.frame(
  day = expected_days,
  vaccine_observed = sapply(expected_days, \(d) cum_events(1,d)),
  vaccine_thesis = expected_cum_v,
  placebo_observed = sapply(expected_days, \(d) cum_events(0,d)),
  placebo_thesis = expected_cum_p
)
validation$vaccine_match <- validation$vaccine_observed == validation$vaccine_thesis
validation$placebo_match <- validation$placebo_observed == validation$placebo_thesis
print(validation)
write.csv(validation, file.path(OUT_DIR, "00_endpoint_validation.csv"), row.names=FALSE)

total_events <- sum(ipd$status)
cat("\nTotal reconstructed events through day 180:", total_events, "\n")
if (total_events != 81L) warning(" reports 81 total events; RDS currently gives ", total_events, ".")

# ---------------- INTERVAL BINOMIAL DATA ----------------
# Risk table boundaries in the thesis are every 30 days.
breaks <- seq(0, FOLLOWUP_DAYS, by=INTERVAL_DAYS)
K <- length(breaks)-1L

at_risk_start <- function(trt_value, t0) {
  # At risk immediately after baseline at t=0; for later intervals, observed time >= t0.
  sum(ipd$trt == trt_value & ipd$time >= t0)
}
events_interval <- function(trt_value, lo, hi, first=FALSE) {
  if (first) {
    sum(ipd$trt == trt_value & ipd$status==1 & ipd$time >= lo & ipd$time <= hi)
  } else {
    sum(ipd$trt == trt_value & ipd$status==1 & ipd$time > lo & ipd$time <= hi)
  }
}

interval_data <- do.call(rbind, lapply(seq_len(K), function(k) {
  lo <- breaks[k]; hi <- breaks[k+1]
  data.frame(
    k=k, lo=lo, hi=hi, t_mid=(lo+hi)/2,
    n_p=at_risk_start(0,lo),
    c_p=events_interval(0,lo,hi,k==1),
    n_v=at_risk_start(1,lo),
    c_v=events_interval(1,lo,hi,k==1)
  )
}))
print(interval_data)
write.csv(interval_data, file.path(OUT_DIR, "01_interval_binomial_data.csv"), row.names=FALSE)

if (sum(interval_data$c_p) != 62L || sum(interval_data$c_v) != 19L) {
  warning("Interval counts do not reproduce totals P=62 and vaccine=19. Check IPD boundary/event times.")
}

# ---------------- STAN MODEL FACTORY ----------------
stan_common_data <- list(
  K=K,
  t_mid=interval_data$t_mid,
  n_p=interval_data$n_p,
  c_p=interval_data$c_p,
  n_v=interval_data$n_v,
  c_v=interval_data$c_v,
  a0=A0, b0=B0,
  sigma_kappa=SIGMA_KAPPA,
  beta_gamma_shape=BETA_GAMMA_SHAPE,
  beta_gamma_rate=BETA_GAMMA_RATE,
  rho_shape=RHO_SHAPE,
  rho_rate=RHO_RATE,
  base_a=BASELINE_RISK_A,
  base_b=BASELINE_RISK_B
)

stan_header <- '
data {
  int<lower=1> K;
  vector<lower=0>[K] t_mid;
  array[K] int<lower=0> n_p;
  array[K] int<lower=0> c_p;
  array[K] int<lower=0> n_v;
  array[K] int<lower=0> c_v;
  real<lower=0> a0;
  real<lower=0> b0;
  real<lower=0> sigma_kappa;
  real<lower=0> beta_gamma_shape;
  real<lower=0> beta_gamma_rate;
  real<lower=0> rho_shape;
  real<lower=0> rho_rate;
  real<lower=0> base_a;
  real<lower=0> base_b;
}
'

stan_exp <- paste0(stan_header, '
parameters {
  real<lower=0,upper=1> VE0;
  real<lower=0> kappa;
  vector<lower=0,upper=1>[K] pi_p;
}
transformed parameters {
  vector[K] VE;
  vector[K] pi_v;
  for (k in 1:K) {
    VE[k] = VE0 * exp(-kappa * t_mid[k]);
    pi_v[k] = pi_p[k] * (1 - VE[k]);
  }
}
model {
  VE0 ~ beta(a0,b0);
  kappa ~ normal(0,sigma_kappa);
  pi_p ~ beta(base_a,base_b);
  c_p ~ binomial(n_p,pi_p);
  c_v ~ binomial(n_v,pi_v);
}
generated quantities {
  vector[2*K] log_lik;
  for (k in 1:K) {
    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k],pi_p[k]);
    log_lik[K+k] = binomial_lpmf(c_v[k] | n_v[k],pi_v[k]);
  }
}
')

stan_erlang <- paste0(stan_header, '
parameters {
  real<lower=0,upper=1> VE0;
  real<lower=0> beta_Gamma;
  vector<lower=0,upper=1>[K] pi_p;
}
transformed parameters {
  vector[K] VE;
  vector[K] pi_v;
  for (k in 1:K) {
    VE[k] = VE0 * (1 - gamma_cdf(t_mid[k] | 3, beta_Gamma));
    pi_v[k] = pi_p[k] * (1 - VE[k]);
  }
}
model {
  VE0 ~ beta(a0,b0);
  beta_Gamma ~ gamma(beta_gamma_shape,beta_gamma_rate);
  pi_p ~ beta(base_a,base_b);
  c_p ~ binomial(n_p,pi_p);
  c_v ~ binomial(n_v,pi_v);
}
generated quantities {
  vector[2*K] log_lik;
  for (k in 1:K) {
    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k],pi_p[k]);
    log_lik[K+k] = binomial_lpmf(c_v[k] | n_v[k],pi_v[k]);
  }
}
')

stan_power <- paste0(stan_header, '
parameters {
  real<lower=0,upper=1> VE0;
  real<lower=0> kappa;
  real<lower=0> rho;
  vector<lower=0,upper=1>[K] pi_p;
}
transformed parameters {
  vector[K] VE;
  vector[K] pi_v;
  for (k in 1:K) {
    VE[k] = VE0 / (1 + pow(kappa*t_mid[k],rho));
    pi_v[k] = pi_p[k] * (1 - VE[k]);
  }
}
model {
  VE0 ~ beta(a0,b0);
  kappa ~ normal(0,sigma_kappa);
  rho ~ gamma(rho_shape,rho_rate);
  pi_p ~ beta(base_a,base_b);
  c_p ~ binomial(n_p,pi_p);
  c_v ~ binomial(n_v,pi_v);
}
generated quantities {
  vector[2*K] log_lik;
  for (k in 1:K) {
    log_lik[k] = binomial_lpmf(c_p[k] | n_p[k],pi_p[k]);
    log_lik[K+k] = binomial_lpmf(c_v[k] | n_v[k],pi_v[k]);
  }
}
')

stan_files <- c(
  exponential=file.path(OUT_DIR,"model_exponential.stan"),
  erlang3=file.path(OUT_DIR,"model_erlang3.stan"),
  powerlaw=file.path(OUT_DIR,"model_powerlaw.stan")
)
writeLines(stan_exp, stan_files["exponential"])
writeLines(stan_erlang, stan_files["erlang3"])
writeLines(stan_power, stan_files["powerlaw"])

# ---------------- FIT ALL CANDIDATE FORMS ----------------
fits <- list()
for (nm in names(stan_files)) {
  cat("\nCompiling/fitting:", nm, "\n")
  mod <- cmdstanr::cmdstan_model(stan_files[[nm]])
  fits[[nm]] <- mod$sample(
    data=stan_common_data,
    seed=SEED,
    chains=CHAINS,
    parallel_chains=CHAINS,
    iter_warmup=ITER_WARMUP,
    iter_sampling=ITER_SAMPLING,
    adapt_delta=ADAPT_DELTA,
    max_treedepth=MAX_TREEDEPTH,
    refresh=500
  )
}

# ---------------- CONVERGENCE DIAGNOSTICS ----------------
diag_rows <- list()
for (nm in names(fits)) {
  s <- fits[[nm]]$summary()
  monitored <- s[!grepl("^log_lik",s$variable),]
  diag_rows[[nm]] <- data.frame(
    model=nm,
    max_rhat=max(monitored$rhat,na.rm=TRUE),
    min_ess_bulk=min(monitored$ess_bulk,na.rm=TRUE),
    min_ess_tail=min(monitored$ess_tail,na.rm=TRUE)
  )
}
diagnostics <- bind_rows(diag_rows)
print(diagnostics)
write.csv(diagnostics,file.path(OUT_DIR,"02_convergence_diagnostics.csv"),row.names=FALSE)

if (any(diagnostics$max_rhat >= 1.01,na.rm=TRUE))
  warning("At least one candidate has R-hat >= 1.01.")
if (any(diagnostics$min_ess_bulk <= 400,na.rm=TRUE))
  warning("At least one candidate has bulk ESS <= 400.")

# ---------------- LOOIC MODEL SELECTION ----------------
loos <- lapply(fits, function(fit) {
  ll <- fit$draws("log_lik", format="matrix")
  loo::loo(ll)
})
loo_cmp <- loo::loo_compare(loos)

loo_table <- data.frame(
  model=rownames(loo_cmp),
  elpd_diff=loo_cmp[,"elpd_diff"],
  se_diff=loo_cmp[,"se_diff"],
  row.names=NULL
)
# LOOIC = -2 * elpd_loo; report absolute LOOIC separately.
loo_abs <- data.frame(
  model=names(loos),
  LOOIC=sapply(loos, function(x) -2*x$estimates["elpd_loo","Estimate"]),
  SE=sapply(loos, function(x) 2*x$estimates["elpd_loo","SE"])
)
loo_abs <- loo_abs[order(loo_abs$LOOIC),]
selected <- loo_abs$model[1]
print(loo_abs)
write.csv(loo_abs,file.path(OUT_DIR,"03_LOOIC_model_comparison.csv"),row.names=FALSE)

cat("\nLOOIC-selected model:", selected, "\n")
if (selected != "erlang3") {
  warning("The thesis reports Erlang-3 as selected for maternal severe RSV-LRTI, ",
          "but this run selected ", selected,
          ". Check exact prior hyperparameters, interval construction, and Erlang parameterisation.")
}

# ---------------- POSTERIOR VE CURVES ----------------
grid <- seq(0,FOLLOWUP_DAYS,by=OUTPUT_STEP)

ve_draw_matrix <- function(fit, model_name, grid) {
  dr <- posterior::as_draws_df(fit$draws())
  if (model_name=="exponential") {
    out <- sapply(grid, function(t) dr$VE0*exp(-dr$kappa*t))
  } else if (model_name=="erlang3") {
    out <- sapply(grid, function(t) dr$VE0*(1-pgamma(t,shape=3,rate=dr$beta_Gamma)))
  } else {
    out <- sapply(grid, function(t) dr$VE0/(1+(dr$kappa*t)^dr$rho))
  }
  as.matrix(out)
}

curve_list <- list()
for (nm in names(fits)) {
  m <- ve_draw_matrix(fits[[nm]],nm,grid)
  curve_list[[nm]] <- data.frame(
    model=nm,
    day=grid,
    median=apply(m,2,median),
    lower=apply(m,2,quantile,0.025),
    upper=apply(m,2,quantile,0.975)
  )
}
curves <- bind_rows(curve_list)
write.csv(curves,file.path(OUT_DIR,"04_all_candidate_VE_curves.csv"),row.names=FALSE)

selected_curve <- curves[curves$model==selected,]
write.csv(selected_curve,file.path(OUT_DIR,"05_selected_model_VE_every_15_days.csv"),row.names=FALSE)

# Table 5.3 M4 values: validation targets only, never fitted as data.
thesis_m4 <- data.frame(
  day=c(0,15,30,45,60,75,90,105,120,135,150,165,180),
  thesis_median=c(.82,.77,.72,.66,.61,.56,.52,.48,.45,.41,.39,.37,.35),
  thesis_lower=c(.68,.65,.60,.55,.50,.44,.38,.33,.28,.24,.20,.18,.17),
  thesis_upper=c(.89,.84,.81,.76,.74,.70,.66,.62,.65,.62,.60,.59,.58)
)
comparison <- merge(selected_curve,thesis_m4,by="day",all.x=TRUE)
comparison$median_difference <- comparison$median-comparison$thesis_median
write.csv(comparison,file.path(OUT_DIR,"06_selected_vs_thesis_Table5_3.csv"),row.names=FALSE)

p <- ggplot(selected_curve,aes(day,median)) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.2) +
  geom_line(linewidth=.9) +
  geom_point(data=thesis_m4,aes(day,thesis_median),inherit.aes=FALSE,size=1.7) +
  coord_cartesian(xlim=c(0,180),ylim=c(0,1)) +
  labs(
    title=paste0("Maternal RSVpreF: severe RSV-LRTI — Bayesian ",selected),
    subtitle="Line/ribbon: current posterior; points: thesis Table 5.3 posterior medians",
    x="Days since birth", y="Risk-scale vaccine efficacy"
  ) +
  theme_bw()
save_plot(p,"07_selected_Bayesian_waning_curve.png")

# Trace plots for scientifically central parameters only.
for (nm in names(fits)) {
  pars <- if (nm=="exponential") c("VE0","kappa") else
          if (nm=="erlang3") c("VE0","beta_Gamma") else c("VE0","kappa","rho")
  d <- posterior::as_draws_df(fits[[nm]]$draws(pars))
  d$.draw_id <- seq_len(nrow(d))
  dl <- tidyr::pivot_longer(as.data.frame(d)[,c(".draw_id",pars)], - .draw_id,
                            names_to="parameter",values_to="value")
  pt <- ggplot(dl,aes(.draw_id,value)) +
    geom_line(alpha=.45) + facet_wrap(~parameter,scales="free_y") +
    labs(title=paste("Trace-like combined draws:",nm),x="Stored draw",y=NULL) +
    theme_bw()
  save_plot(pt,paste0("trace_",nm,".png"),10,6)
}

# Save summaries and CmdStan file references.
param_summary <- bind_rows(lapply(names(fits), function(nm) {
  x <- fits[[nm]]$summary()
  x$model <- nm
  x
}))
write.csv(param_summary,file.path(OUT_DIR,"08_parameter_summary.csv"),row.names=FALSE)

saveRDS(
  list(
    endpoint="Maternal RSVpreF - medically attended severe RSV-LRTI",
    followup_days=FOLLOWUP_DAYS,
    interval_data=interval_data,
    endpoint_validation=validation,
    looic=loo_abs,
    selected_model=selected,
    selected_curve=selected_curve,
    thesis_validation_targets=thesis_m4,
    assumptions=list(
      VE0_prior_mean=VE0_PRIOR_MEAN,
      VE0_prior_concentration=VE0_PRIOR_CONCENTRATION,
      sigma_kappa=SIGMA_KAPPA,
      baseline_risk_prior=c(BASELINE_RISK_A,BASELINE_RISK_B),
      rho_prior=c(RHO_SHAPE,RHO_RATE),
      erlang_beta_interpretation="Stan gamma rate"
    ),
    cmdstan_output_files=lapply(fits, \(x) x$output_files())
  ),
  file.path(OUT_DIR,"Severe_RSV_LRTI_Bayesian_summary.rds")
)

cat("\n============================================================\n")
cat("Analysis complete.\n")
cat("Selected model:",selected,"\n")
cat("The expected selected model: erlang3\n")
cat("Outputs:",OUT_DIR,"\n")
cat("IMPORTANT: exact numerical reproduction still depends on the thesis-unreported\n")
cat("prior hyperparameters and the original Erlang beta/time parameterisation.\n")
cat("============================================================\n")
