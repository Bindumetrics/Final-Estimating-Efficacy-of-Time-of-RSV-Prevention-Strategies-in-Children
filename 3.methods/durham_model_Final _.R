################################################################################
# DURHAM'S SMOOTHED SCHOENFELD RESIDUAL METHOD
# 
#
# Method specification:
#   1. Fit a standard Cox PH model to the full dataset.
#   2. Obtain scaled Schoenfeld residuals at event times.
#   3. Use E[r_i^* | T_i=t] ≈ beta(t) - beta_hat.
#   4. Fit a cubic smoothing spline to r_i^* versus event time.
#   5. Select the smoothing parameter by GCV.
#   6. beta_hat(t) = beta_hat + s_hat(t).
#   7. VE_h(t) = 1 - exp{beta_hat(t)}.
#   8. Evaluate VE at 15-day intervals.
#   9. Report pointwise 95% CIs
################################################################################

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tibble)
  library(mgcv)
  library(ggplot2)
})

#-----------------------------#
# 1) Input validation / prep  #
#-----------------------------#
durham_prep_data <- function(data,
                             time_col   = "time",
                             status_col = "status",
                             treat_col  = "treat",
                             ref_level  = "placebo",
                             vax_level  = "nirsevimab") {

  stopifnot(is.data.frame(data))

  if (!all(c(time_col, status_col, treat_col) %in% names(data))) {
    stop("Data must contain columns: ",
         paste(c(time_col, status_col, treat_col), collapse = ", "))
  }

  dat <- data %>%
    transmute(
      time   = as.numeric(.data[[time_col]]),
      status = as.integer(.data[[status_col]]),
      treat  = factor(as.character(.data[[treat_col]]),
                      levels = c(ref_level, vax_level))
    )

  if (anyNA(dat$time) || anyNA(dat$status) || anyNA(dat$treat)) {
    stop("Missing/invalid values found in time, status or treatment.")
  }

  if (!all(dat$status %in% c(0L, 1L))) {
    stop("status must be coded 0/1.")
  }

  if (any(dat$time <= 0)) {
    stop("All survival times must be strictly positive.")
  }

  if (sum(dat$status) < 4) {
    warning("Very few events: the Durham curve and pointwise SEs may be unstable.")
  }

  dat
}

#------------------------------------#
# 2) Standard Cox PH model           #
#------------------------------------#
durham_fit_cox <- function(dat, ties = "efron") {

  fit <- coxph(
    Surv(time, status) ~ treat,
    data = dat,
    ties = ties,
    x = TRUE,
    model = TRUE
  )

  beta_hat <- unname(coef(fit)[1])
  se_beta  <- sqrt(vcov(fit)[1, 1])

  list(
    fit      = fit,
    beta_hat = beta_hat,
    se_beta  = se_beta,
    HR       = exp(beta_hat),
    VE       = 1 - exp(beta_hat),
    n        = nrow(dat),
    events   = sum(dat$status)
  )
}

#--------------------------------------------------#
# 3) Grambsch–Therneau test for time-varying beta #
#--------------------------------------------------#
durham_ph_test <- function(cox_fit, transform = "identity") {

  zph <- cox.zph(cox_fit, transform = transform)

  list(
    zph      = zph,
    p_treat  = unname(zph$table[1, "p"]),
    p_global = unname(zph$table[nrow(zph$table), "p"])
  )
}

#------------------------------------------------------------------#
# 4) Durham time-varying coefficient using cubic smoothing spline  #
#------------------------------------------------------------------#
durham_curve <- function(dat,
                         eval_times = NULL,
                         eval_every = 15,
                         max_followup = NULL,
                         zph_transform = "identity",
                         ties = "efron",
                         spline_k = NULL) {

  # Standard Cox model
  cox_obj <- durham_fit_cox(dat, ties = ties)
  fit <- cox_obj$fit
  beta_hat <- cox_obj$beta_hat

  # Scaled Schoenfeld residual process
  # With a single treatment coefficient:
  # E[r_i^* | T_i=t] ≈ beta(t) - beta_hat
  zph <- cox.zph(fit, transform = zph_transform)

  residual_process <- tibble(
    time = as.numeric(zph$x),
    r_scaled = as.numeric(zph$y[, 1])
  ) %>%
    filter(is.finite(time), is.finite(r_scaled)) %>%
    arrange(time)

  if (nrow(residual_process) < 4) {
    stop("Too few usable event-time residuals to fit the Durham smoother.")
  }

  # Evaluation grid: report estimates every 15 days
  if (is.null(max_followup)) {
    max_followup <- max(dat$time, na.rm = TRUE)
  }

  if (is.null(eval_times)) {
    eval_times <- seq(0, max_followup, by = eval_every)
  }

  # A cubic penalised spline is fitted to the scaled Schoenfeld residuals.
  # method = "GCV.Cp" selects the smoothing parameter by GCV.
  #
  # bs = "cr" gives a cubic regression spline with the usual curvature
  # penalty on the second derivative.
  if (is.null(spline_k)) {
    # Keep the basis moderate relative to the number of events.
    spline_k <- min(20L, max(5L, floor(nrow(residual_process) / 3)))
  }

  spline_k <- min(spline_k, nrow(residual_process) - 1L)
  spline_k <- max(spline_k, 4L)

  spline_fit <- gam(
    r_scaled ~ s(time, bs = "cr", k = spline_k),
    data = residual_process,
    method = "GCV.Cp"
  )

  # Predict the smoothed residual component and its pointwise SE
  pred <- predict(
    spline_fit,
    newdata = data.frame(time = eval_times),
    type = "response",
    se.fit = TRUE
  )

  s_hat <- as.numeric(pred$fit)
  se_s  <- as.numeric(pred$se.fit)

  # Report specification:
  # beta_hat(t) = beta_hat + s_hat(t)
  beta_t <- beta_hat + s_hat

  # Pointwise SE for beta(t).
  # describes uncertainty as arising from the fitted
  # time-varying coefficient/smoother. We combine the global Cox
  # coefficient variance with the smoother's pointwise variance.
  #
  #
  # This matches the pointwise-SE procedure.
  se_beta_t <- sqrt(cox_obj$se_beta^2 + se_s^2)

  zcrit <- qnorm(0.975)

  beta_lo <- beta_t - zcrit * se_beta_t
  beta_hi <- beta_t + zcrit * se_beta_t

  # VE_h(t) = 1 - exp{beta(t)}
  VE_t <- 1 - exp(beta_t)

  # Because VE = 1-exp(beta) is monotonically decreasing in beta,
  # beta upper -> VE lower; beta lower -> VE upper.
  VE_lower <- 1 - exp(beta_hi)
  VE_upper <- 1 - exp(beta_lo)

  curve <- tibble(
    time       = eval_times,
    beta_t     = beta_t,
    SE_beta_t  = se_beta_t,
    HR_t       = exp(beta_t),
    VE_t       = VE_t,
    VE_lower   = VE_lower,
    VE_upper   = VE_upper
  )

  list(
    cox              = cox_obj,
    zph              = zph,
    residual_process = residual_process,
    spline_fit       = spline_fit,
    curve            = curve,
    settings         = list(
      eval_every     = eval_every,
      max_followup   = max_followup,
      zph_transform  = zph_transform,
      ties           = ties,
      spline_basis   = "cubic regression spline",
      spline_k       = spline_k,
      smoothing_rule = "GCV"
    )
  )
}

#---------------------------------------------#
# 5) One-shot analysis                        #
#---------------------------------------------#
durham_waning_analysis <- function(data,
                                   time_col = "time",
                                   status_col = "status",
                                   treat_col = "treat",
                                   ref_level = "placebo",
                                   vax_level = "nirsevimab",
                                   eval_times = NULL,
                                   eval_every = 15,
                                   max_followup = NULL,
                                   zph_transform = "identity",
                                   ties = "efron",
                                   spline_k = NULL) {

  dat <- durham_prep_data(
    data = data,
    time_col = time_col,
    status_col = status_col,
    treat_col = treat_col,
    ref_level = ref_level,
    vax_level = vax_level
  )

  main <- durham_curve(
    dat = dat,
    eval_times = eval_times,
    eval_every = eval_every,
    max_followup = max_followup,
    zph_transform = zph_transform,
    ties = ties,
    spline_k = spline_k
  )

  ph <- durham_ph_test(
    main$cox$fit,
    transform = zph_transform
  )

  list(
    data = dat,
    cox = main$cox,
    ph = ph,
    residual_process = main$residual_process,
    spline_fit = main$spline_fit,
    curve = main$curve,
    settings = main$settings
  )
}

#---------------------------------------------#
# 6) Plot helper                              #
#---------------------------------------------#
plot_durham_ve <- function(res,
                           ylim = c(0, 1),
                           percent = FALSE) {

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
      title = "Durham smoothed Schoenfeld-residual VE(t)",
      subtitle = paste0(
        "Cubic smoothing spline; GCV-selected smoothing; ",
        "95% pointwise confidence intervals"
      ),
      x = "Days since dose/birth",
      y = ylab
    ) +
    theme_bw()
}

#---------------------------------------------#
# 7) Example usage                            #
#---------------------------------------------#
# res <- durham_waning_analysis(
#   data = mydata,
#   time_col = "time",
#   status_col = "status",
#   treat_col = "treat",
#   ref_level = "placebo",
#   vax_level = "nirsevimab",
#   eval_every = 15,
#   max_followup = 150,       # use 180 or 360 for maternal-vaccine analyses
#   zph_transform = "identity"
# )
#
# res$ph$p_treat
# res$curve
# plot_durham_ve(res, ylim = c(0, 1))
#
# For maternal vaccine:
#   vax_level = "maternal_vaccine"
#   max_followup = 180 or 360
################################################################################
