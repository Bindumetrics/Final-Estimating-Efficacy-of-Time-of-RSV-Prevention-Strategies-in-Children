################################################################################
# TDC COX MODEL WITH LOG(TIME) INTERACTION
#
#
# specification:
#   lambda(t | Z_i) = lambda_0(t) exp{ beta Z_i + gamma Z_i log(t) }
#
#   eta(t)   = beta + gamma log(t)
#   HR(t)    = exp{eta(t)}
#   VE_h(t)  = 1 - exp{eta(t)}
#
#   beta and gamma are estimated jointly by Cox partial likelihood.
#   Pointwise 95% CIs are obtained from the variance-covariance matrix
#   of (beta, gamma) using:
#
#     Var[eta(t)] =
#       Var(beta) + log(t)^2 Var(gamma)
#       + 2 log(t) Cov(beta,gamma)
#
#   CI_VE(t) =
#     [1-exp{eta(t)+1.96 SE_eta(t)},
#      1-exp{eta(t)-1.96 SE_eta(t)}]
#
#
################################################################################

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

#===============================================================================
# 1) Prepare data
#===============================================================================
tdc_prep_data <- function(data,
                          time_col   = "time",
                          status_col = "status",
                          treat_col  = "treat",
                          ref_level  = "placebo",
                          vax_level  = "nirsevimab") {

  stopifnot(is.data.frame(data))

  if (!all(c(time_col, status_col, treat_col) %in% names(data))) {
    stop(
      "Data must contain columns: ",
      paste(c(time_col, status_col, treat_col), collapse = ", ")
    )
  }

  dat <- data %>%
    transmute(
      time   = as.numeric(.data[[time_col]]),
      status = as.integer(.data[[status_col]]),
      treat  = factor(
        as.character(.data[[treat_col]]),
        levels = c(ref_level, vax_level)
      )
    )

  if (anyNA(dat$time) || anyNA(dat$status) || anyNA(dat$treat)) {
    stop("Missing or invalid values found in time, status or treatment.")
  }

  if (!all(dat$status %in% c(0L, 1L))) {
    stop("status must be coded 0/1.")
  }

  if (any(dat$time <= 0)) {
    stop(
      "The report specifies log(t), so all Cox analysis times must be > 0. ",
      "If time=0 observations exist, check the reconstruction/time origin."
    )
  }

  dat %>%
    mutate(vax = as.integer(treat == vax_level))
}


#===============================================================================
# 2) Fit TDC model: beta*Z + gamma*Z*log(t)
#===============================================================================
tdc_fit_logt <- function(data,
                         time_col   = "time",
                         status_col = "status",
                         treat_col  = "treat",
                         ref_level  = "placebo",
                         vax_level  = "nirsevimab",
                         ties       = "efron") {

  dat <- tdc_prep_data(
    data = data,
    time_col = time_col,
    status_col = status_col,
    treat_col = treat_col,
    ref_level = ref_level,
    vax_level = vax_level
  )

  fit <- coxph(
    Surv(time, status) ~ vax + tt(vax),
    data = dat,
    ties = ties,
    x = TRUE,
    model = TRUE,
    tt = function(x, t, ...) x * log(t)
  )

  cf <- coef(fit)

  beta_name  <- "vax"
  gamma_name <- grep("^tt\\(vax\\)", names(cf), value = TRUE)

  if (length(gamma_name) != 1L) {
    stop("Could not uniquely identify the time-interaction coefficient.")
  }

  beta_hat  <- unname(cf[beta_name])
  gamma_hat <- unname(cf[gamma_name])

  V <- vcov(fit)[c(beta_name, gamma_name),
                 c(beta_name, gamma_name),
                 drop = FALSE]

  list(
    data = dat,
    fit = fit,
    beta_hat = beta_hat,
    gamma_hat = gamma_hat,
    vcov = V,
    ref_level = ref_level,
    vax_level = vax_level,
    ties = ties
  )
}


#===============================================================================
# 3) Pointwise VE(t) and delta-method 95% CIs
#===============================================================================
tdc_predict_logt <- function(fit_obj,
                             eval_times,
                             level = 0.95) {

  if (any(eval_times <= 0)) {
    stop(
      "The fitted report model uses log(t), which is undefined at t=0. ",
      "Use strictly positive evaluation times."
    )
  }

  beta_hat  <- fit_obj$beta_hat
  gamma_hat <- fit_obj$gamma_hat
  V         <- fit_obj$vcov

  lt <- log(eval_times)

  # eta(t) = beta + gamma log(t)
  eta <- beta_hat + gamma_hat * lt

  # Delta-method variance of eta(t):
  # [1, log(t)] V [1, log(t)]'
  var_eta <- (
    V[1, 1] +
    (lt^2) * V[2, 2] +
    2 * lt * V[1, 2]
  )

  # Guard against tiny negative numerical values
  var_eta <- pmax(var_eta, 0)
  se_eta  <- sqrt(var_eta)

  zcrit <- qnorm(1 - (1 - level) / 2)

  eta_lower <- eta - zcrit * se_eta
  eta_upper <- eta + zcrit * se_eta

  HR_t <- exp(eta)
  VE_t <- 1 - HR_t

  # Monotone reversal on VE scale:
  # larger eta -> larger HR -> smaller VE
  VE_lower <- 1 - exp(eta_upper)
  VE_upper <- 1 - exp(eta_lower)

  tibble(
    time = eval_times,
    eta_t = eta,
    SE_eta_t = se_eta,
    HR_t = HR_t,
    VE_t = VE_t,
    VE_lower = VE_lower,
    VE_upper = VE_upper
  )
}


#===============================================================================
# 4) AIC comparison of candidate time functions
#    The log(t) was selected because
#    it gave the lowest AIC among candidate specifications.
#
#    Candidate models here:
#      - constant Cox effect
#      - linear time interaction
#      - log(time) interaction
#      - square-root time interaction
#
#    All are fitted to the same data and use the same Cox partial-likelihood
#    framework, so their AIC values are directly comparable.
#===============================================================================
tdc_compare_aic <- function(data,
                            time_col   = "time",
                            status_col = "status",
                            treat_col  = "treat",
                            ref_level  = "placebo",
                            vax_level  = "nirsevimab",
                            ties       = "efron") {

  dat <- tdc_prep_data(
    data = data,
    time_col = time_col,
    status_col = status_col,
    treat_col = treat_col,
    ref_level = ref_level,
    vax_level = vax_level
  )

  fit_constant <- coxph(
    Surv(time, status) ~ vax,
    data = dat,
    ties = ties
  )

  fit_linear <- coxph(
    Surv(time, status) ~ vax + tt(vax),
    data = dat,
    ties = ties,
    tt = function(x, t, ...) x * t
  )

  fit_log <- coxph(
    Surv(time, status) ~ vax + tt(vax),
    data = dat,
    ties = ties,
    tt = function(x, t, ...) x * log(t)
  )

  fit_sqrt <- coxph(
    Surv(time, status) ~ vax + tt(vax),
    data = dat,
    ties = ties,
    tt = function(x, t, ...) x * sqrt(t)
  )

  out <- tibble(
    specification = c(
      "Constant treatment effect",
      "Linear time: Z*t",
      "Log time: Z*log(t)",
      "Square-root time: Z*sqrt(t)"
    ),
    AIC = c(
      AIC(fit_constant),
      AIC(fit_linear),
      AIC(fit_log),
      AIC(fit_sqrt)
    )
  ) %>%
    arrange(AIC) %>%
    mutate(
      delta_AIC = AIC - min(AIC),
      selected = row_number() == 1L
    )

  attr(out, "fits") <- list(
    constant = fit_constant,
    linear = fit_linear,
    log = fit_log,
    sqrt = fit_sqrt
  )

  out
}


#===============================================================================
# 5) One-shot analysis evaluated every 15 days
#===============================================================================
tdc_logt_analysis <- function(data,
                              time_col   = "time",
                              status_col = "status",
                              treat_col  = "treat",
                              ref_level  = "placebo",
                              vax_level  = "nirsevimab",
                              eval_every = 15,
                              max_followup = NULL,
                              first_eval = 1,
                              level = 0.95,
                              ties = "efron") {

  fit_obj <- tdc_fit_logt(
    data = data,
    time_col = time_col,
    status_col = status_col,
    treat_col = treat_col,
    ref_level = ref_level,
    vax_level = vax_level,
    ties = ties
  )

  if (is.null(max_followup)) {
    max_followup <- max(fit_obj$data$time, na.rm = TRUE)
  }

  if (first_eval <= 0) {
    stop("first_eval must be > 0 because the report model uses log(t).")
  }

  # 15-day reporting grid, while respecting log(t)>-Inf.
  #
  # Example with first_eval=15:
  #   15, 30, 45, ..., 150
  #
  # Example with first_eval=1:
  #   1, 16, 31, ...
  #
  # But we have 15,30,45,... so we use first_eval=15.
  eval_times <- seq(first_eval, max_followup, by = eval_every)

  # Ensure end-of-follow-up can be included if it is not exactly on grid
  if (tail(eval_times, 1) < max_followup) {
    eval_times <- c(eval_times, max_followup)
  }

  curve <- tdc_predict_logt(
    fit_obj = fit_obj,
    eval_times = eval_times,
    level = level
  )

  list(
    fit_obj = fit_obj,
    curve = curve,
    settings = list(
      model = "Cox TDC with Z*log(t)",
      eval_every = eval_every,
      first_eval = first_eval,
      max_followup = max_followup,
      confidence = level,
      CI_method = "delta method from vcov(beta, gamma)",
      ties = ties
    )
  )
}


#===============================================================================
# 6) Plot helper
#===============================================================================
plot_tdc_logt <- function(res,
                          ylim = c(0, 1),
                          percent = FALSE) {

  b0 <- res$fit_obj$beta_hat
  b1 <- res$fit_obj$gamma_hat

  pdat <- res$curve

  if (percent) {
    pdat <- pdat %>%
      mutate(
        VE_t = 100 * VE_t,
        VE_lower = 100 * VE_lower,
        VE_upper = 100 * VE_upper
      )

    ylim <- 100 * ylim
    ylab <- "Vaccine efficacy (%)"
  } else {
    ylab <- "Vaccine efficacy"
  }

  ggplot(pdat, aes(x = time, y = VE_t)) +
    geom_ribbon(
      aes(ymin = VE_lower, ymax = VE_upper),
      alpha = 0.20
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    coord_cartesian(ylim = ylim) +
    labs(
      title = "TDC Cox model: time-varying vaccine efficacy",
      subtitle = paste0(
        "log HR(t) = beta + gamma log(t); beta = ",
        round(b0, 3),
        ", gamma = ",
        round(b1, 3),
        "; 95% pointwise delta-method CIs"
      ),
      x = "Days since dose / birth",
      y = ylab
    ) +
    theme_bw()
}


#===============================================================================
# 7) Example usage
#===============================================================================

# A) Confirm that log(t) is the lowest-AIC candidate:
#
# aic_table <- tdc_compare_aic(
#   data = mydata,
#   time_col = "time",
#   status_col = "status",
#   treat_col = "treat",
#   ref_level = "placebo",
#   vax_level = "nirsevimab"
# )
# print(aic_table)


# B) Final log(t) TDC model:
#
#
# res <- tdc_logt_analysis(
#   data = mydata,
#   time_col = "time",
#   status_col = "status",
#   treat_col = "treat",
#   ref_level = "placebo",
#   vax_level = "nirsevimab",
#   eval_every = 15,
#   first_eval = 15,
#   max_followup = 150
# )
#
# res$curve
# plot_tdc_logt(res, percent = TRUE)


