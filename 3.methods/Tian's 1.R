################################################################################
# TIAN'S KERNEL-WEIGHTED PARTIAL LIKELIHOOD METHOD
#
#   1. Estimate beta(t0) locally at each evaluation time t0.
#   2. Use kernel-weighted Cox partial likelihood.
#   3. Epanechnikov kernel:
#          K(u) = 3/4 * (1-u^2) * I(|u| <= 1)
#          K_b(u) = K(u/b) / b
#   4. Select bandwidth b by leave-one-event-out cross-validation.
#   5. VE_h(t0) = 1 - exp{beta_hat(t0)}.
#   6. Evaluate VE at 15-day intervals.
#   7. Construct a simultaneous confidence band across the evaluation grid.
#
################################################################################

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

#===============================================================================
# 1. DATA PREPARATION
#===============================================================================

tian_prep_data <- function(data,
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
  
  if (anyNA(dat$time) ||
      anyNA(dat$status) ||
      anyNA(dat$treat)) {
    stop("Missing/invalid values found in time, status or treatment.")
  }
  
  if (!all(dat$status %in% c(0L, 1L))) {
    stop("status must be coded 0/1.")
  }
  
  if (any(dat$time <= 0)) {
    stop("All survival times must be strictly positive.")
  }
  
  dat$Z <- as.integer(dat$treat == vax_level)
  
  if (sum(dat$status) < 4) {
    warning(
      "Very few events: Tian's local estimates may be unstable."
    )
  }
  
  dat
}


#===============================================================================
# 2. EPANECHNIKOV KERNEL
#===============================================================================

epanechnikov_kernel <- function(u) {
  
  out <- numeric(length(u))
  
  inside <- abs(u) <= 1
  
  out[inside] <- 0.75 * (1 - u[inside]^2)
  
  out
}


kernel_weight <- function(event_time, t0, bandwidth) {
  
  if (bandwidth <= 0) {
    stop("Bandwidth must be positive.")
  }
  
  u <- (event_time - t0) / bandwidth
  
  epanechnikov_kernel(u) / bandwidth
}


#===============================================================================
# 3. LOCAL KERNEL-WEIGHTED COX PARTIAL LOG-LIKELIHOOD
#===============================================================================

tian_local_loglik <- function(beta,
                              dat,
                              t0,
                              bandwidth,
                              leave_out_event = NULL) {
  
  event_rows <- which(dat$status == 1L)
  
  if (!is.null(leave_out_event)) {
    event_rows <- setdiff(event_rows, leave_out_event)
  }
  
  loglik <- 0
  
  for (i in event_rows) {
    
    ti <- dat$time[i]
    
    w_i <- kernel_weight(
      event_time = ti,
      t0         = t0,
      bandwidth  = bandwidth
    )
    
    # Event contributes nothing outside the kernel window
    if (w_i <= 0) next
    
    # Cox risk set at event time ti
    risk <- which(dat$time >= ti)
    
    eta <- beta * dat$Z[risk]
    
    # Stable log(sum(exp(eta)))
    m <- max(eta)
    
    log_denom <- m + log(sum(exp(eta - m)))
    
    contribution <-
      beta * dat$Z[i] -
      log_denom
    
    loglik <- loglik + w_i * contribution
  }
  
  loglik
}


#===============================================================================
# 4. FIT beta(t0) AT ONE EVALUATION TIME
#===============================================================================

tian_fit_local <- function(dat,
                           t0,
                           bandwidth,
                           beta_interval = c(-8, 4)) {
  
  event_times <- dat$time[dat$status == 1L]
  
  local_events <- event_times[
    abs(event_times - t0) <= bandwidth
  ]
  
  if (length(local_events) < 2) {
    
    return(
      list(
        beta = NA_real_,
        se   = NA_real_,
        loglik = NA_real_,
        n_local_events = length(local_events)
      )
    )
  }
  
  objective <- function(beta) {
    
    -tian_local_loglik(
      beta      = beta,
      dat       = dat,
      t0        = t0,
      bandwidth = bandwidth
    )
  }
  
  fit <- tryCatch(
    optimize(
      f        = objective,
      interval = beta_interval
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    
    return(
      list(
        beta = NA_real_,
        se   = NA_real_,
        loglik = NA_real_,
        n_local_events = length(local_events)
      )
    )
  }
  
  beta_hat <- fit$minimum
  
  # Numerical observed information
  h <- 1e-4
  
  l0 <- tian_local_loglik(
    beta_hat,
    dat,
    t0,
    bandwidth
  )
  
  lp <- tian_local_loglik(
    beta_hat + h,
    dat,
    t0,
    bandwidth
  )
  
  lm <- tian_local_loglik(
    beta_hat - h,
    dat,
    t0,
    bandwidth
  )
  
  second_derivative <-
    (lp - 2 * l0 + lm) / h^2
  
  information <- -second_derivative
  
  se_beta <-
    if (is.finite(information) && information > 0) {
      sqrt(1 / information)
    } else {
      NA_real_
    }
  
  list(
    beta = beta_hat,
    se   = se_beta,
    loglik = l0,
    n_local_events = length(local_events)
  )
}


#===============================================================================
# 5. LEAVE-ONE-EVENT-OUT CROSS-VALIDATION
#===============================================================================

tian_cv_bandwidth <- function(dat,
                              bandwidth_grid = NULL,
                              beta_interval = c(-8, 4),
                              verbose = TRUE) {
  
  event_rows <- which(dat$status == 1L)
  event_times <- dat$time[event_rows]
  
  if (length(event_rows) < 4) {
    stop("Too few events for leave-one-event-out bandwidth selection.")
  }
  
  followup <- max(dat$time)
  
  if (is.null(bandwidth_grid)) {
    
    bandwidth_grid <- unique(
      round(
        seq(
          max(15, 0.05 * followup),
          max(30, 0.40 * followup),
          length.out = 15
        )
      )
    )
  }
  
  cv_results <- tibble(
    bandwidth = bandwidth_grid,
    cv_score  = NA_real_
  )
  
  for (b_index in seq_along(bandwidth_grid)) {
    
    b <- bandwidth_grid[b_index]
    
    if (verbose) {
      message(
        "Evaluating bandwidth ",
        b,
        " (",
        b_index,
        "/",
        length(bandwidth_grid),
        ")"
      )
    }
    
    score <- 0
    usable <- 0
    
    for (k in seq_along(event_rows)) {
      
      omitted_row <- event_rows[k]
      t0 <- event_times[k]
      
      # Fit local beta excluding this event contribution
      objective <- function(beta) {
        
        -tian_local_loglik(
          beta            = beta,
          dat             = dat,
          t0              = t0,
          bandwidth       = b,
          leave_out_event = omitted_row
        )
      }
      
      fit <- tryCatch(
        optimize(
          objective,
          interval = beta_interval
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit)) next
      
      beta_minus_i <- fit$minimum
      
      # Predictive Cox contribution of omitted event
      ti <- dat$time[omitted_row]
      
      risk <- which(dat$time >= ti)
      
      eta <- beta_minus_i * dat$Z[risk]
      
      m <- max(eta)
      
      log_denom <-
        m + log(sum(exp(eta - m)))
      
      predictive_loglik <-
        beta_minus_i * dat$Z[omitted_row] -
        log_denom
      
      score <- score + predictive_loglik
      usable <- usable + 1
    }
    
    cv_results$cv_score[b_index] <-
      if (usable > 0) score else NA_real_
  }
  
  if (all(is.na(cv_results$cv_score))) {
    stop("Bandwidth cross-validation failed.")
  }
  
  best_index <-
    which.max(cv_results$cv_score)
  
  best_bandwidth <-
    cv_results$bandwidth[best_index]
  
  list(
    bandwidth = best_bandwidth,
    table     = cv_results
  )
}


#===============================================================================
# 6. ESTIMATE THE COMPLETE TIAN VE CURVE
#===============================================================================

tian_curve <- function(dat,
                       eval_times = NULL,
                       eval_every = 15,
                       max_followup = NULL,
                       bandwidth = NULL,
                       bandwidth_grid = NULL,
                       beta_interval = c(-8, 4),
                       verbose = TRUE) {
  
  if (is.null(max_followup)) {
    max_followup <- max(dat$time)
  }
  
  if (is.null(eval_times)) {
    
    eval_times <- seq(
      0,
      max_followup,
      by = eval_every
    )
  }
  
  #---------------------------------------------------------------
  # Bandwidth selection
  #---------------------------------------------------------------
  
  if (is.null(bandwidth)) {
    
    cv <- tian_cv_bandwidth(
      dat            = dat,
      bandwidth_grid = bandwidth_grid,
      beta_interval  = beta_interval,
      verbose        = verbose
    )
    
    bandwidth <- cv$bandwidth
    
  } else {
    
    cv <- NULL
  }
  
  if (verbose) {
    message(
      "Selected Tian bandwidth = ",
      bandwidth,
      " days"
    )
  }
  
  #---------------------------------------------------------------
  # Local estimates
  #---------------------------------------------------------------
  
  fits <- lapply(
    eval_times,
    function(t0) {
      
      tian_fit_local(
        dat           = dat,
        t0            = t0,
        bandwidth     = bandwidth,
        beta_interval = beta_interval
      )
    }
  )
  
  beta_hat <-
    vapply(
      fits,
      function(x) x$beta,
      numeric(1)
    )
  
  se_beta <-
    vapply(
      fits,
      function(x) x$se,
      numeric(1)
    )
  
  n_local_events <-
    vapply(
      fits,
      function(x) x$n_local_events,
      numeric(1)
    )
  
  #---------------------------------------------------------------
  # Convert log-HR to vaccine efficacy
  #---------------------------------------------------------------
  
  VE <- 1 - exp(beta_hat)
  
  # Pointwise intervals retained internally
  beta_lo_point <- beta_hat - 1.96 * se_beta
  beta_hi_point <- beta_hat + 1.96 * se_beta
  
  # Note reversal when converting beta CI to VE
  VE_lo_point <- 1 - exp(beta_hi_point)
  VE_hi_point <- 1 - exp(beta_lo_point)
  
  curve <- tibble(
    time            = eval_times,
    beta            = beta_hat,
    se_beta         = se_beta,
    VE              = VE,
    VE_lo_pointwise = VE_lo_point,
    VE_hi_pointwise = VE_hi_point,
    n_local_events  = n_local_events
  )
  
  list(
    curve      = curve,
    bandwidth  = bandwidth,
    cv         = cv,
    eval_times = eval_times
  )
}


#===============================================================================
# 7. SIMULTANEOUS CONFIDENCE BAND
#===============================================================================

tian_simultaneous_band <- function(dat,
                                   fitted_object,
                                   B = 500,
                                   alpha = 0.05,
                                   seed = 12345,
                                   beta_interval = c(-8, 4),
                                   verbose = TRUE) {
  
  set.seed(seed)
  
  original <- fitted_object$curve
  eval_times <- fitted_object$eval_times
  bandwidth <- fitted_object$bandwidth
  
  beta_original <- original$beta
  se_original <- original$se_beta
  
  valid_original <-
    is.finite(beta_original) &
    is.finite(se_original) &
    se_original > 0
  
  max_statistics <- rep(NA_real_, B)
  
  n <- nrow(dat)
  
  for (b in seq_len(B)) {
    
    if (verbose && (b %% 25 == 0 || b == 1)) {
      message("Simultaneous-band bootstrap: ", b, "/", B)
    }
    
    boot_index <- sample(
      seq_len(n),
      size = n,
      replace = TRUE
    )
    
    boot_dat <- dat[boot_index, , drop = FALSE]
    
    boot_fits <- lapply(
      eval_times,
      function(t0) {
        
        tian_fit_local(
          dat           = boot_dat,
          t0            = t0,
          bandwidth     = bandwidth,
          beta_interval = beta_interval
        )
      }
    )
    
    beta_boot <-
      vapply(
        boot_fits,
        function(x) x$beta,
        numeric(1)
      )
    
    valid <-
      valid_original &
      is.finite(beta_boot)
    
    if (sum(valid) < 2) next
    
    standardized_difference <-
      abs(
        (beta_boot[valid] -
           beta_original[valid]) /
          se_original[valid]
      )
    
    max_statistics[b] <-
      max(
        standardized_difference,
        na.rm = TRUE
      )
  }
  
  max_statistics <-
    max_statistics[is.finite(max_statistics)]
  
  if (length(max_statistics) < 20) {
    stop(
      "Too few successful bootstrap replicates for simultaneous band."
    )
  }
  
  critical_value <-
    unname(
      quantile(
        max_statistics,
        probs = 1 - alpha,
        na.rm = TRUE
      )
    )
  
  beta_lower <-
    beta_original -
    critical_value * se_original
  
  beta_upper <-
    beta_original +
    critical_value * se_original
  
  # Reverse endpoints under VE = 1-exp(beta)
  VE_lower <-
    1 - exp(beta_upper)
  
  VE_upper <-
    1 - exp(beta_lower)
  
  fitted_object$curve <-
    fitted_object$curve %>%
    mutate(
      beta_lo_simultaneous = beta_lower,
      beta_hi_simultaneous = beta_upper,
      VE_lo_simultaneous   = VE_lower,
      VE_hi_simultaneous   = VE_upper
    )
  
  fitted_object$simultaneous_band <-
    list(
      alpha          = alpha,
      critical_value = critical_value,
      successful_B   = length(max_statistics),
      requested_B    = B
    )
  
  fitted_object
}


#===============================================================================
# 8. COMPLETE WRAPPER
#===============================================================================

fit_tian <- function(data,
                     time_col   = "time",
                     status_col = "status",
                     treat_col  = "treat",
                     ref_level  = "placebo",
                     vax_level  = "nirsevimab",
                     eval_every = 15,
                     max_followup = NULL,
                     bandwidth_grid = NULL,
                     B = 500,
                     alpha = 0.05,
                     seed = 12345,
                     verbose = TRUE) {
  
  dat <- tian_prep_data(
    data        = data,
    time_col    = time_col,
    status_col  = status_col,
    treat_col   = treat_col,
    ref_level   = ref_level,
    vax_level   = vax_level
  )
  
  fit <- tian_curve(
    dat            = dat,
    eval_every     = eval_every,
    max_followup   = max_followup,
    bandwidth_grid = bandwidth_grid,
    verbose        = verbose
  )
  
  fit <- tian_simultaneous_band(
    dat           = dat,
    fitted_object = fit,
    B             = B,
    alpha         = alpha,
    seed          = seed,
    verbose       = verbose
  )
  
  fit$data <- dat
  
  fit
}


#===============================================================================
# 9. PLOT
#===============================================================================

plot_tian <- function(fit) {
  
  pdat <- fit$curve
  
  ggplot(
    pdat,
    aes(x = time, y = VE)
  ) +
    geom_ribbon(
      aes(
        ymin = VE_lo_simultaneous,
        ymax = VE_hi_simultaneous
      ),
      alpha = 0.20
    ) +
    geom_line(
      linewidth = 0.9
    ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    labs(
      x = "Time (days)",
      y = "Vaccine efficacy",
      title = "Tian kernel-weighted partial likelihood",
      subtitle = paste0(
        "Epanechnikov kernel; LOEO-CV bandwidth = ",
        fit$bandwidth,
        " days"
      )
    ) +
    theme_minimal()
}


################################################################################
# 
################################################################################
#
#   time    = observed follow-up time
#   status  = 1 event, 0 censored
#   treat   = "placebo" or "nirsevimab"
#
# fit_tian_nirs <- fit_tian(
#     data          = ipd_nirsevimab,
#     time_col      = "time",
#     status_col    = "status",
#     treat_col     = "treat",
#     ref_level     = "placebo",
#     vax_level     = "nirsevimab",
#     eval_every    = 15,
#     max_followup  = 150,
#     bandwidth_grid = seq(20, 90, by = 5),
#     B             = 500,
#     seed          = 12345
# )
#
# fit_tian_nirs$bandwidth
# fit_tian_nirs$cv$table
# fit_tian_nirs$curve
#
# plot_tian(fit_tian_nirs)
################################################################################