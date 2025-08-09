# ===============================================================================
# MaxLik_fit_functions_reduced_model.R
# 
# Clean maximum likelihood functions for epidemic modeling analysis
# 
# This file contains all the essential maximum likelihood estimation functions
# for fitting reduced SEIR models with heterogeneity parameter v and 
# non-pharmaceutical interventions to epidemic data.
#
# Key functions:
# - Single epidemic fitting: fit4_reducedm_loglik.NPI()
# - Two epidemics fitting: fit4_2epic_reduced.loglik.NPI()
# - Profile likelihood analysis for both scenarios
# - Supporting likelihood and objective functions
# 
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# Date: August 2025
# ===============================================================================

# Load required libraries
library(deSolve)     # For ODE solving
library(tidyverse)   # For data manipulation

# ===============================================================================
# SECTION 1: CORE LIKELIHOOD FUNCTIONS
# ===============================================================================

#' Poisson log-likelihood function for reduced SEIR model with NPI
#' 
#' This function calculates the log-likelihood of observed epidemic data 
#' given model parameters. It integrates the reduced SEIR model and compares
#' the predicted daily incidence with observed cases using Poisson likelihood.
#' 
#' Note: This function expects the global variable 'times' to be set
#' 
#' @param params Named vector of model parameters including R0, v, intervention parameters
#' @param sim.data Data frame with 'time' and 'reports' columns
#' @param initial_state Initial conditions [S, E, I, R, C]
#' @param times Time vector (optional, uses global 'times' if NULL)
#' @return Log-likelihood value
poisson.loglik.NPI.reduced <- function(params, sim.data, initial_state, times = NULL) {
  
  # Use global times variable if not provided (matches original code pattern)
  if (is.null(times)) {
    if (!exists("times", envir = .GlobalEnv)) {
      stop("The 'times' variable must be defined in the global environment or passed as argument")
    }
    times <- get("times", envir = .GlobalEnv)
  }
  
  # Solve the ODE system using our utility function
  out <- as.data.frame(ode(
    y = initial_state, 
    times = times, 
    func = Reduced.m_intervene, 
    parms = params
  ))
  
  # Calculate daily incidence from cumulative cases
  # Note: We exclude the first time point (t=0) to match data structure
  Daily_incidence <- c(0, diff(out[, "C"]))
  df <- out %>% mutate(Inc = Daily_incidence)
  
  # Ensure lambda values are positive for Poisson likelihood
  # Small positive value prevents log(0) errors
  lambda_ <- pmax(df[,"Inc"], 0.0001)
  
  # Validate input data structure
  if (!is.data.frame(sim.data)) {
    sim.data <- as.data.frame(sim.data)
  }
  
  if (!"reports" %in% names(sim.data)) {
    stop("The 'reports' column is missing in sim.data")
  }
  
  if (nrow(sim.data) != nrow(df)) {
    stop("The lengths of sim.data and model output do not match")
  }
  
  # Calculate Poisson log-likelihood
  # sum(log(P(observed | expected))) where P is Poisson pmf
  loglik <- sum(dpois(
    x = sim.data[,"reports"], 
    lambda = lambda_, 
    log = TRUE
  ))
  
  return(loglik)
}

# ===============================================================================
# SECTION 2: SINGLE EPIDEMIC FUNCTIONS
# ===============================================================================

#' Objective function for single epidemic optimization
#' 
#' This function transforms parameters from optimization space to natural space
#' and returns the negative log-likelihood for minimization.
#' 
#' Parameter transformations:
#' - R0: exp(par[1]) to ensure positivity
#' - v: exp(par[2]) to ensure positivity  
#' - t0: exp(par[3]) to ensure positivity
#' - c_value2: expit(par[4]) to ensure 0 < c_value2 < 1
#' 
#' @param par Vector of transformed parameters [log(R0), log(v), log(t0), logit(c_value2)]
#' @param sim.data Observed epidemic data
#' @return Negative log-likelihood for minimization
f4_optim_reducedm.NPI <- function(par, sim.data) {
  
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),           # Basic reproduction number (positive)
    v = exp(par[2]),            # Heterogeneity parameter (positive)
    t0 = exp(par[3]),           # Intervention start time (positive)
    t1 = t1_spec,               # Fixed intervention timing
    t2 = t2_spec,               # Fixed intervention timing  
    t3 = t3_spec,               # Fixed intervention timing
    c_value1 = c_value1_spec,   # Fixed baseline transmission
    c_value2 = expit(par[4]),   # Intervention strength (0-1)
    c_value3 = c_value3_spec,   # Fixed final transmission
    rho = rho_spec,             # Fixed relative infectiousness
    delta = delta_spec,         # Fixed incubation rate
    gamma = gamma_spec,         # Fixed recovery rate
    N = N,                      # Fixed population size
    tfinal = tfinal_spec        # Fixed simulation time
  )
  
  # Calculate log-likelihood and return negative for minimization
  loglik <- poisson.loglik.NPI.reduced(
    params = params, 
    sim.data = sim.data, 
    initial_state = initial_state
  )
  
  return(-loglik)
}

#' Fit reduced SEIR model to single epidemic data
#' 
#' This is the main function for fitting the reduced SEIR model with NPI
#' to a single epidemic dataset. It estimates R0, v, t0, and c_value2.
#' 
#' @param dat Data frame with 'time' and 'reports' columns
#' @return List containing:
#'   - parms: Parameter estimates on natural scale
#'   - trans_parms: Parameter estimates on transformed scale  
#'   - trans_hessian: Hessian matrix for confidence intervals
fit4_reducedm_loglik.NPI <- function(dat) {
  
  # Starting values for optimization (on transformed scale)
  # These are reasonable defaults: R0~2, v~2, t0~12 days, c_value2~0.4
  start_par <- c(log(2), log(2), log(12), logit(0.4))
  
  # Run optimization using Nelder-Mead algorithm
  fit1 <- optim(
    par = start_par,
    fn = f4_optim_reducedm.NPI,
    sim.data = dat,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 1000),
    hessian = TRUE
  )
  
  # Transform fitted parameters back to natural scale
  fittedparams <- c(
    R0 = exp(fit1$par[1]),
    v = exp(fit1$par[2]),
    t0 = exp(fit1$par[3]),
    c_value2 = expit(fit1$par[4]),
    AIC = 2 * length(fit1$par) + 2 * fit1$value,  # Akaike Information Criterion
    negloglik = fit1$value,                        # Minimized negative log-likelihood
    loglik = -fit1$value,                         # Maximized log-likelihood
    convergence = fit1$convergence                 # Convergence status (0 = success)
  )
  
  return(list(
    parms = fittedparams,        # Natural scale parameters
    trans_parms = fit1$par,      # Transformed scale parameters
    trans_hessian = fit1$hessian # Hessian matrix
  ))
}

# ===============================================================================
# SECTION 3: TWO EPIDEMICS FUNCTIONS  
# ===============================================================================

#' Objective function for two epidemics optimization
#' 
#' This function fits the same model parameters to two concurrent epidemics
#' with different initial conditions. The total likelihood is the sum of
#' individual epidemic likelihoods (assuming independence).
#' 
#' @param par Vector of transformed parameters
#' @param sim.data_1 First epidemic dataset
#' @param sim.data_2 Second epidemic dataset  
#' @return Combined negative log-likelihood
f4_2epi_optim_reduced.NPI <- function(par, sim.data_1, sim.data_2) {
  
  # Transform parameters to natural scale (same as single epidemic)
  params <- c(
    R0 = exp(par[1]),
    v = exp(par[2]),
    t0 = exp(par[3]),
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[4]),
    c_value3 = c_value3_spec,
    rho = rho_spec,
    delta = delta_spec,
    gamma = gamma_spec,
    N = N,
    tfinal = tfinal_spec
  )
  
  # Calculate log-likelihood for first epidemic
  loglik1 <- poisson.loglik.NPI.reduced(
    params = params, 
    sim.data = sim.data_1, 
    initial_state = initial_state_1
  )
  
  # Calculate log-likelihood for second epidemic  
  loglik2 <- poisson.loglik.NPI.reduced(
    params = params, 
    sim.data = sim.data_2, 
    initial_state = initial_state_2
  )
  
  # Sum log-likelihoods (valid for independent epidemics)
  combined_loglik = loglik1 + loglik2
  
  # Return negative combined log-likelihood for minimization
  return(-combined_loglik)
}

#' Fit reduced SEIR model to two concurrent epidemics
#' 
#' This function simultaneously fits the reduced SEIR model to two epidemic
#' datasets, improving parameter identifiability by leveraging information
#' from multiple epidemics with different initial conditions.
#' 
#' @param dat1 First epidemic dataset
#' @param dat2 Second epidemic dataset
#' @return List containing parameter estimates and diagnostics
fit4_2epic_reduced.loglik.NPI <- function(dat1, dat2) {
  
  # Starting values (same as single epidemic)
  start_par <- c(log(2), log(2), log(12), logit(0.4))
  
  # Run optimization
  fit1 <- optim(
    par = start_par,
    fn = f4_2epi_optim_reduced.NPI,
    sim.data_1 = dat1,
    sim.data_2 = dat2,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 1500),  # More iterations for complexity
    hessian = TRUE
  )
  
  # Transform fitted parameters to natural scale
  fittedparams <- c(
    R0 = exp(fit1$par[1]),
    v = exp(fit1$par[2]),
    t0 = exp(fit1$par[3]),
    c_value2 = expit(fit1$par[4]),
    AIC = 2 * length(fit1$par) + 2 * fit1$value,
    negloglik = fit1$value,
    loglik = -fit1$value,
    convergence = fit1$convergence
  )
  
  return(list(
    parms = fittedparams,
    trans_parms = fit1$par,
    trans_hessian = fit1$hessian
  ))
}


# ===============================================================================
# OPTIMIZATION FUNCTIONS FOR HOMOGENEOUS MODEL
# ===============================================================================
poisson.loglik.withNPI.reduced <- function(params, sim.data, initial_state) {
  # Get the time vector from sim.data
  times_data <- sim.data$time
  
  # Check if we need to add time 0 for integration
  if (min(times_data) > 0) {
    # Data starts at time 1, so add time 0 for ODE integration
    times_integration <- c(0, times_data)
    needs_trimming <- TRUE
  } else {
    # Data already includes time 0
    times_integration <- times_data
    needs_trimming <- FALSE
  }
  
  # Integrate the model equations
  out <- as.data.frame(ode(
    y = initial_state,
    times = times_integration,
    func = Reduced.m_intervene,
    parms = params
  ))
  
  # If we added time 0, remove it to match sim.data
  if (needs_trimming) {
    out <- out[-1, ]  # Remove first row (time 0)
  }
  
  # Calculate daily incidence from cumulative cases
  # IMPORTANT: Don't prepend 0 - just use diff directly
  Daily_incidence <- diff(c(0, out[, "C"]))  # This ensures first value is out$C[1] - 0
  
  # Add incidence to dataframe
  df <- out %>% mutate(Inc = Daily_incidence)
  
  # Ensure values are valid for Poisson likelihood
  lambda_ <- ifelse(df[,"Inc"] <= 0, 0.0001, df[,"Inc"])
  
  # Validate sim.data
  if (!is.data.frame(sim.data)) {
    sim.data <- as.data.frame(sim.data)
  }
  
  # Check if the "reports" column exists
  if (!"reports" %in% names(sim.data)) {
    stop("The 'reports' column is missing in sim.data")
  }
  
  # Check if the lengths match
  if (nrow(sim.data) != length(lambda_)) {
    cat("Debug: sim.data rows =", nrow(sim.data), ", lambda length =", length(lambda_), "\n")
    cat("Debug: df rows =", nrow(df), "\n")
    cat("Debug: sim.data time range =", range(sim.data$time), "\n")
    cat("Debug: df time range =", range(df$time), "\n")
    stop("The lengths of sim.data and df do not match")
  }
  
  # Compute log likelihood
  loglik <- sum(dpois(
    x = sim.data[,"reports"],
    lambda = lambda_,
    log = TRUE
  ))
  
  return(loglik)
}

# ===============================================================================
# Homogeneous model functions (single epidemic)
# ===============================================================================

#' Objective function for homogeneous model with NPI
#'
#' @param par Vector of transformed parameters to be estimated
#' @param sim.data Observed epidemic data
#' @return Negative log-likelihood value
f_optim_reducedm.poisloglikwithNPI <- function(par, sim.data) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),
    v = v_spec,  # Fixed at 0 for homogeneous model
    t0 = exp(par[2]), 
    t1 = t1_spec, 
    t2 = t2_spec, 
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[3]),
    c_value3 = c_value3_spec,
    rho = rho_spec,
    delta = delta_spec,
    gamma = gamma_spec,
    N = N, 
    tfinal = tfinal_spec
  )
  
  # Calculate the negative log-likelihood
  loglik <- -poisson.loglik.withNPI.reduced(
    params, 
    sim.data = sim.data, 
    initial_state = initial_state
  )
  
  return(loglik)
}

#' Function to fit homogeneous model to a single epidemic with NPI
#'
#' @param dat Data frame with time and reports columns
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit3_hom_1epic_loglikwithNPI <- function(dat) {
  fit <- optim(
    par = c(log(2), log(12), logit(0.3)), 
    fn = f_optim_reducedm.poisloglikwithNPI,
    sim.data = dat,
    method = "Nelder-Mead", 
    control = list(trace = 0, maxit = 1500),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit$par[1]),
    t0 = exp(fit$par[2]),
    c_value2 = expit(fit$par[3]),
    AIC = 2 * length(fit$par) - 2 * (-fit$value),
    value = fit$value,
    convergence = fit$convergence
  )
  
  return(list(
    parms = fittedparams,
    trans_parms = fit$par,
    trans_hessian = fit$hessian
  ))
}

# ===============================================================================
# Homogeneous model functions (dual epidemic)
# ===============================================================================

#' Objective function for homogeneous model with NPI for two epidemics
#'
#' @param par Vector of transformed parameters to be estimated
#' @param sim.data_1 First epidemic dataset
#' @param sim.data_2 Second epidemic dataset
#' @return Combined negative log-likelihood
f3_optim_reducedm.poisloglikwithNPI <- function(par, sim.data_1, sim.data_2) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),
    v = v_spec,  # Fixed at 0 for homogeneous model
    t0 = exp(par[2]), 
    t1 = t1_spec, 
    t2 = t2_spec, 
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[3]),
    c_value3 = c_value3_spec,
    rho = rho_spec,
    delta = delta_spec,
    gamma = gamma_spec,
    N = N, 
    tfinal = tfinal_spec
  )
  
  # Calculate log-likelihood for first epidemic
  loglik1 <- -poisson.loglik.withNPI.reduced(
    params, 
    sim.data = sim.data_1,
    initial_state = initial_state_1
  )
  
  # Calculate log-likelihood for second epidemic
  loglik2 <- -poisson.loglik.withNPI.reduced(
    params, 
    sim.data = sim.data_2,
    initial_state = initial_state_2
  )
  
  # Sum negative log-likelihoods
  loglik = loglik1 + loglik2
  
  return(loglik)
}

#' Function to fit homogeneous model to two epidemics with NPI
#'
#' @param dat1 First epidemic dataset
#' @param dat2 Second epidemic dataset
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit3_hom_2epic_loglikwithNPI <- function(dat1, dat2) {
  fit <- optim(
    par = c(log(2), log(12), logit(0.2)), 
    fn = f3_optim_reducedm.poisloglikwithNPI,
    sim.data_1 = dat1,
    sim.data_2 = dat2,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 1600),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit$par[1]),
    t0 = exp(fit$par[2]),
    c_value2 = expit(fit$par[3]),
    AIC = 2 * length(fit$par) - 2 * (-fit$value),
    value = fit$value,
    convergence = fit$convergence
  )
  
  return(list(
    parms = fittedparams,
    trans_parms = fit$par,
    trans_hessian = fit$hessian
  ))
}

# ===============================================================================
# Additional utility functions for likelihood calculations
# ===============================================================================

#' Log-likelihood function for homogeneous model with NPI
#'
#' @param params Model parameters
#' @param sim.data Observed epidemic data
#' @param initial_state Initial state vector for the SEIR model
#' @return Log-likelihood value




# ===============================================================================
# Additional MaxLik Functions for No-NPI Models
# 
# These functions handle fitting models without non-pharmaceutical interventions
# Note: We reuse existing simulation functions by setting c_value2_spec = 1
# ===============================================================================

# ===============================================================================
# Homogeneous model without NPIs (single epidemic)
# ===============================================================================

#' Objective function for homogeneous model without NPIs
#'
#' @param par Vector of transformed parameters to be estimated (only R0)
#' @param sim.data Observed epidemic data
#' @return Negative log-likelihood value
f_optim_hom_noNPI <- function(par, sim.data) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),           # R0 (positive)
    v = v_spec,                 # Fixed coefficient of variation (0 for homogeneous)
    rho = rho_spec,             # Relative infectiousness in E compartment (fixed)
    delta = delta_spec,         # Rate of transition from E to I (fixed)
    gamma = gamma_spec,         # Recovery rate (fixed)
    N = N                       # Population size (fixed)
  )
  
  # Calculate the log-likelihood using existing function but with c_value2 = 1 (no NPI)
  loglik <- poisson.loglik.withNPI.reduced(
    params = c(params, 
               t0 = t0_spec, t1 = tfinal_spec, t2 = tfinal_spec, t3 = tfinal_spec,
               c_value1 = 1, c_value2 = 1, c_value3 = 1, tfinal = tfinal_spec), 
    sim.data = sim.data, 
    initial_state = initial_state
  )
  
  # Return negative log-likelihood for minimization
  return(-loglik)
}

#' Function to fit homogeneous model to single epidemic without NPIs
#'
#' @param dat Data frame with time and reports columns
#' @return List containing parameter estimates, transformed parameters, and Hessian

# ===============================================================================
# Additional MaxLik Functions for No-NPI Models
# 
# These functions should be added to MaxLik_fit_functions.R
# They handle fitting models without non-pharmaceutical interventions
# Note: We reuse existing simulation functions by setting c_value2_spec = 1
# ===============================================================================

# ===============================================================================
# Homogeneous model without NPIs (single epidemic)
# ===============================================================================

#' Objective function for homogeneous model without NPIs
#'
#' @param par Vector of transformed parameters to be estimated (only R0)
#' @param sim.data Observed epidemic data
#' @return Negative log-likelihood value
f_optim_hom_noNPI <- function(par, sim.data) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),           # R0 (positive)
    v = v_spec,                 # Fixed coefficient of variation (0 for homogeneous)
    rho = rho_spec,             # Relative infectiousness in E compartment (fixed)
    delta = delta_spec,         # Rate of transition from E to I (fixed)
    gamma = gamma_spec,         # Recovery rate (fixed)
    N = N                       # Population size (fixed)
  )
  
  # Calculate the log-likelihood using existing function but with c_value2 = 1 (no NPI)
  loglik <- poisson.loglik.withNPI.reduced(
    params = c(params, 
               t0 = t0_spec, t1 = tfinal_spec, t2 = tfinal_spec, t3 = tfinal_spec,
               c_value1 = 1, c_value2 = 1, c_value3 = 1, tfinal = tfinal_spec), 
    sim.data = sim.data, 
    initial_state = initial_state
  )
  
  # Return negative log-likelihood for minimization
  return(-loglik)
}

#' Function to fit homogeneous model to single epidemic without NPIs
#'
#' @param dat Data frame with time and reports columns
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit1_hom_1epic_loglik_noNPI <- function(dat) {
  # Starting values for optimization (on transformed scale)
  start_par <- c(log(2))  # Only R0
  
  # Run optimization to find maximum likelihood estimates
  fit1 <- optim(
    par = start_par,
    fn = f_optim_hom_noNPI,
    sim.data = dat,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 1500),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit1$par[1]),
    # Calculate AIC: 2*k + 2*minimized_negative_loglik
    AIC = 2 * length(fit1$par) + 2 * fit1$value,
    # Store the minimized negative log-likelihood
    negloglik = fit1$value,
    # Store the maximized log-likelihood
    loglik = -fit1$value,
    # Convergence code
    convergence = fit1$convergence
  )
  
  # Return results
  return(list(
    parms = fittedparams,        # Parameter estimates on natural scale
    trans_parms = fit1$par,      # Parameter estimates on transformed scale
    trans_hessian = fit1$hessian # Hessian matrix (for confidence intervals)
  ))
}

# ===============================================================================
# Heterogeneous model without NPIs (single epidemic) - Using Reduced Model
# ===============================================================================

#' Objective function for heterogeneous reduced model without NPIs
#'
#' @param par Vector of transformed parameters (R0, v)
#' @param sim.data Observed epidemic data
#' @return Negative log-likelihood value
f2_optim_het_reduced_noNPI <- function(par, sim.data) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),          # Basic reproduction number
    v = exp(par[2]),           # Coefficient of variation
    rho = rho_spec,            # Relative infectiousness in E compartment (fixed)
    delta = delta_spec,        # Rate of transition from E to I (fixed)
    gamma = gamma_spec,        # Recovery rate (fixed)
    N = N,                     # Population size (fixed)
    t0 = t0_spec,                    # No intervention timing needed
    t1 = tfinal_spec,          # Set intervention times beyond simulation period
    t2 = tfinal_spec,
    t3 = tfinal_spec,
    c_value1 = 1,              # No intervention effect
    c_value2 = 1,              # No intervention effect (key: = 1 means no NPI)
    c_value3 = 1,              # No intervention effect
    tfinal = tfinal_spec
  )
  
  # Calculate the log-likelihood using existing reduced model function
  loglik <- poisson.loglik.NPI.reduced(
    params, 
    sim.data = sim.data, 
    initial_state = initial_state,
    times = times
  )
  
  # Return negative log-likelihood for minimization
  return(-loglik)
}

#' Function to fit heterogeneous reduced model to single epidemic without NPIs
#'
#' @param dat Data frame with time and reports columns
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit2_het_1epic_loglik_noNPI <- function(dat) {
  # Starting values for optimization (on transformed scale)
  start_par <- c(log(2), log(2))  # R0 and v
  
  # Run optimization to find maximum likelihood estimates
  fit1 <- optim(
    par = start_par,
    fn = f2_optim_het_reduced_noNPI,
    sim.data = dat,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 1800),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit1$par[1]),
    v = exp(fit1$par[2]),
    # Calculate AIC: 2*k + 2*minimized_negative_loglik
    AIC = 2 * length(fit1$par) + 2 * fit1$value,
    # Store the minimized negative log-likelihood
    negloglik = fit1$value,
    # Store the maximized log-likelihood
    loglik = -fit1$value,
    # Convergence code
    convergence = fit1$convergence
  )
  
  # Return results
  return(list(
    parms = fittedparams,        # Parameter estimates on natural scale
    trans_parms = fit1$par,      # Parameter estimates on transformed scale
    trans_hessian = fit1$hessian # Hessian matrix (for confidence intervals)
  ))
}



# ===============================================================================
# SECTION 4: PROFILE LIKELIHOOD FUNCTIONS
# ===============================================================================

#' Profile likelihood analysis for single epidemic
#' 
#' This function performs profile likelihood analysis for a specified parameter,
#' providing confidence intervals and assessing parameter identifiability.
#' 
#' @param sim_data Epidemic dataset
#' @param param_to_profile Parameter name: "R0", "v", "t0", or "c_value2"
#' @param value_range Range of parameter values to profile (auto-computed if NULL)
#' @param n_points Number of points for profile (default: 15)
#' @param plot Whether to create profile likelihood plot (default: TRUE)
#' @return List containing profile results, plot, MLE, and confidence intervals
profile_likelihood_reducedm <- function(sim_data, param_to_profile, value_range = NULL, 
                                        n_points = 15, plot = TRUE, custom_title = NULL) {
  
  # First fit the model to get MLE
  mle_fit <- fit4_reducedm_loglik.NPI(dat = sim_data)
  mle_params <- mle_fit$parms
  mle_trans_params <- mle_fit$trans_parms
  
  # Process Hessian matrix to set appropriate parameter range
  if (!is.null(mle_fit$trans_hessian)) {
    tryCatch({
      z_variance <- solve(mle_fit$trans_hessian)
      z_se <- sqrt(diag(z_variance))
      
      # Identify parameter index
      if(param_to_profile == "R0") param_index <- 1
      else if(param_to_profile == "v") param_index <- 2
      else if(param_to_profile == "t0") param_index <- 3
      else if(param_to_profile == "c_value2") param_index <- 4
      
      # Set range based on Hessian if not provided
      if(is.null(value_range)) {
        se <- z_se[param_index]
        trans_lower <- mle_trans_params[param_index] - 2.5 * se
        trans_upper <- mle_trans_params[param_index] + 2.5 * se
        
        if(param_to_profile %in% c("R0", "v", "t0")) {
          value_range <- seq(exp(trans_lower), exp(trans_upper), length.out = n_points)
        } else if(param_to_profile == "c_value2") {
          value_range <- seq(expit(trans_lower), expit(trans_upper), length.out = n_points)
        }
      }
    }, error = function(e) {
      # Fallback to default range if Hessian is not invertible
      if(is.null(value_range)) {
        if(param_to_profile == "R0") {
          value_range <- seq(max(0.5, mle_params["R0"] * 0.7), mle_params["R0"] * 1.3, length.out = n_points)
        } else if(param_to_profile == "v") {
          value_range <- seq(max(0.1, mle_params["v"] * 0.7), mle_params["v"] * 1.3, length.out = n_points)
        } else if(param_to_profile == "t0") {
          value_range <- seq(max(1, mle_params["t0"] * 0.8), mle_params["t0"] * 1.2, length.out = n_points)
        } else if(param_to_profile == "c_value2") {
          value_range <- seq(max(0.1, mle_params["c_value2"] * 0.7), min(0.9, mle_params["c_value2"] * 1.3), length.out = n_points)
        }
      }
    })
  }
  
  # Initialize results dataframe
  profile_results <- data.frame(
    param_value = value_range,
    neg_loglik = NA,
    R0 = NA,
    v = NA,
    t0 = NA,
    c_value2 = NA
  )
  
  # Get the true values from global environment if available
  true_values <- list(
    R0 = if(exists("R0_spec")) R0_spec else mle_params["R0"],
    v = if(exists("CV_true")) CV_true else mle_params["v"],
    t0 = if(exists("t0_spec")) t0_spec else mle_params["t0"],
    c_value2 = if(exists("c_value2_spec")) c_value2_spec else mle_params["c_value2"]
  )
  
  # Set initial state if not already set globally
  if(!exists("initial_state", envir = .GlobalEnv)) {
    assign("initial_state", c(S = 99920, E = 60, I = 24, R = 0, C = 0), envir = .GlobalEnv)
  }
  
  # Profile over parameter values
  for(i in 1:length(value_range)) {
    profile_value <- value_range[i]
    
    # Transform to the appropriate scale
    if(param_to_profile == "R0") {
      profile_trans_value <- log(profile_value)
      param_index <- 1
    } else if(param_to_profile == "v") {
      profile_trans_value <- log(profile_value)
      param_index <- 2
    } else if(param_to_profile == "t0") {
      profile_trans_value <- log(profile_value)
      param_index <- 3
    } else if(param_to_profile == "c_value2") {
      profile_trans_value <- logit(profile_value)
      param_index <- 4
    }
    
    # Get other parameter indices
    other_indices <- setdiff(1:4, param_index)
    start_par <- mle_trans_params[other_indices]
    
    # Optimize the other parameters while fixing the profiled one
    profile_obj_fn <- function(par, fixed_param_index, fixed_value) {
      full_par <- numeric(4)
      full_par[fixed_param_index] <- fixed_value
      full_par[other_indices] <- par
      return(f4_optim_reducedm.NPI(full_par, sim_data))
    }
    
    opt_result <- tryCatch({
      optim(
        par = start_par,
        fn = profile_obj_fn,
        fixed_param_index = param_index,
        fixed_value = profile_trans_value,
        method = "Nelder-Mead",
        control = list(maxit = 1000)
      )
    }, error = function(e) {
      list(value = NA, par = rep(NA, length(start_par)))
    })
    
    # Store results
    profile_results$neg_loglik[i] <- opt_result$value
    
    # Store optimized parameters
    opt_full_par <- numeric(4)
    opt_full_par[param_index] <- profile_trans_value
    opt_full_par[other_indices] <- opt_result$par
    
    # Transform back to original scale
    profile_results$R0[i] <- exp(opt_full_par[1])
    profile_results$v[i] <- exp(opt_full_par[2])
    profile_results$t0[i] <- exp(opt_full_par[3])
    profile_results$c_value2[i] <- expit(opt_full_par[4])
  }
  
  # Remove NA values
  profile_results <- profile_results[!is.na(profile_results$neg_loglik), ]
  
  # Calculate confidence intervals
  min_neg_loglik <- min(profile_results$neg_loglik)
  profile_results$LR_stat <- 2 * (profile_results$neg_loglik - min_neg_loglik)
  conf_threshold <- qchisq(0.95, df = 1)
  
  ci_data <- profile_results[profile_results$LR_stat <= conf_threshold, ]
  if(nrow(ci_data) > 0) {
    ci_lower <- min(ci_data$param_value)
    ci_upper <- max(ci_data$param_value)
  } else {
    ci_lower <- min(profile_results$param_value)
    ci_upper <- max(profile_results$param_value)
  }
  
  # Extract values for plotting
  mle_x <- mle_params[param_to_profile]
  true_x <- true_values[[param_to_profile]]
  
  # Create plots if requested
  if(plot) {
    # Set up parameter names and titles
   # param_labels <- list(
      #R0 = expression(R[0]),
     # v = "CV",
    #  t0 = expression("t0"),
   #   c_value2 = "c"
   # )
    param_labels <- list(
    R0 = expression(R[0]),           # R with subscript 0
    v = expression(nu),              # Greek letter nu
    t0 = expression(t[0]),           # t with subscript 0
    c_value2 = expression(c[1])      # c with subscript 1
    )
    
    # Create customizable title
    if(!is.null(custom_title)) {
      profile_title <- custom_title
    } else {
      epidemic_type <- "Single Epidemic"
      npi_text <- if(exists("c_value2_spec")) paste0("NPI = ", c_value2_spec) else ""
      profile_title <- paste0("Profile for ", param_labels[[param_to_profile]], 
                              " (", epidemic_type, 
                              if(npi_text != "") paste0(", ", npi_text), ")")
    }
    
    # Create the profile likelihood plot with enhanced visibility
    p <- ggplot(profile_results, aes(x = param_value, y = -neg_loglik)) +
      geom_line(color = "steelblue", size = 1.2) +
      #geom_point(size = 2.5, color = "darkblue", alpha = 0.8) +
      geom_point(color = "darkblue", size = 2) +
      geom_hline(yintercept = -min_neg_loglik - conf_threshold/2, 
                 linetype = "dashed", color = "red", size = 1.2, alpha = 0.7) +
      geom_vline(xintercept = mle_x, color = "purple", linetype = "solid", size = 1.1) +
      geom_vline(xintercept = true_x, color = "green4", linetype = "dotted", size = 1.5) +
      geom_vline(xintercept = ci_lower, color = "blue", linetype = "dotted", size = 1.2) +
      geom_vline(xintercept = ci_upper, color = "blue", linetype = "dotted", size = 1.2) +
      annotate("text", x = mle_x, y = max(-profile_results$neg_loglik), 
               label = "MLE", color = "purple", hjust = -0.3, size = 4, fontface = "bold") +
      #annotate("text", x = true_x, y = max(-profile_results$neg_loglik) - 0.8, 
      # label = "True", color = "green4", hjust = -0.3, size = 4, fontface = "bold") +
      labs(#title = profile_title,
        #subtitle = paste0("95% CI: [", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
        x = param_labels[[param_to_profile]],
        y = "Log Likelihood") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
      )
    
    
    
    
    
    print(p)
  }
  
  # Return results
  return(list(
    profile_results = profile_results,
    profile_plot = if(plot) p else NULL,
    mle = mle_params,
    ci = c(ci_lower, ci_upper),
    ci_width = ci_upper - ci_lower
  ))
}




# Custom two-epidemic objective function wrapper for profile likelihood
two_epi_profile_obj_fn <- function(par, fixed_param_index, fixed_value, sim_data_1, sim_data_2) {
  # Insert the fixed parameter value
  full_par <- numeric(4)
  full_par[fixed_param_index] <- fixed_value
  full_par[setdiff(1:4, fixed_param_index)] <- par
  
  # Call the original objective function
  return(f4_2epi_optim_reduced.NPI(full_par, sim_data_1, sim_data_2))
}


#' Profile likelihood analysis for two epidemics
#' 
#' Similar to single epidemic profile likelihood but using two concurrent
#' epidemics for improved parameter identifiability.
#' 
#' @param sim_data_large First (larger) epidemic dataset
#' @param sim_data_small Second (smaller) epidemic dataset  
#' @param param_to_profile Parameter to profile
#' @param value_range Parameter range (auto-computed if NULL)
#' @param n_points Number of profile points
#' @param plot Whether to create plot
#' @return Profile likelihood results
profile_likelihood_two_epidemics <- function(sim_data_large, sim_data_small, param_to_profile, 
                                              value_range = NULL, n_points = 15, plot = F, 
                                              custom_title = NULL) {
  
  # First fit the model to get MLE
  mle_fit <- fit4_2epic_reduced.loglik.NPI(dat1 = sim_data_large, dat2 = sim_data_small)
  mle_params <- mle_fit$parms
  mle_trans_params <- mle_fit$trans_parms
  
  # Process Hessian matrix to set appropriate parameter range
  if (!is.null(mle_fit$trans_hessian)) {
    tryCatch({
      z_variance <- solve(mle_fit$trans_hessian)
      z_se <- sqrt(diag(z_variance))
      
      # Identify parameter index
      if(param_to_profile == "R0") param_index <- 1
      else if(param_to_profile == "v") param_index <- 2
      else if(param_to_profile == "t0") param_index <- 3
      else if(param_to_profile == "c_value2") param_index <- 4
      
      # Set range based on Hessian if not provided
      if(is.null(value_range)) {
        se <- z_se[param_index]
        trans_lower <- mle_trans_params[param_index] - 2.5 * se
        trans_upper <- mle_trans_params[param_index] + 2.5 * se
        
        if(param_to_profile %in% c("R0", "v", "t0")) {
          value_range <- seq(exp(trans_lower), exp(trans_upper), length.out = n_points)
        } else if(param_to_profile == "c_value2") {
          value_range <- seq(expit(trans_lower), expit(trans_upper), length.out = n_points)
        }
      }
    }, error = function(e) {
      # Fallback to default range if Hessian is not invertible
      if(is.null(value_range)) {
        if(param_to_profile == "R0") {
          value_range <- seq(max(0.5, mle_params["R0"] * 0.7), mle_params["R0"] * 1.3, length.out = n_points)
        } else if(param_to_profile == "v") {
          value_range <- seq(max(0.1, mle_params["v"] * 0.7), mle_params["v"] * 1.3, length.out = n_points)
        } else if(param_to_profile == "t0") {
          value_range <- seq(max(1, mle_params["t0"] * 0.8), mle_params["t0"] * 1.2, length.out = n_points)
        } else if(param_to_profile == "c_value2") {
          value_range <- seq(max(0.1, mle_params["c_value2"] * 0.7), min(0.9, mle_params["c_value2"] * 1.3), length.out = n_points)
        }
      }
    })
  }
  
  # Initialize results dataframe
  profile_results <- data.frame(
    param_value = value_range,
    neg_loglik = NA,
    R0 = NA,
    v = NA,
    t0 = NA,
    c_value2 = NA
  )
  
  # Get the true values from global environment if available
  true_values <- list(
    R0 = if(exists("R0_spec")) R0_spec else mle_params["R0"],
    v = if(exists("CV_true")) CV_true else mle_params["v"],
    t0 = if(exists("t0_spec")) t0_spec else mle_params["t0"],
    c_value2 = if(exists("c_value2_spec")) c_value2_spec else mle_params["c_value2"]
  )
  
  
  
  # Profile over parameter values
  for(i in 1:length(value_range)) {
    profile_value <- value_range[i]
    
    # Transform to the appropriate scale
    if(param_to_profile == "R0") {
      profile_trans_value <- log(profile_value)
      param_index <- 1
    } else if(param_to_profile == "v") {
      profile_trans_value <- log(profile_value)
      param_index <- 2
    } else if(param_to_profile == "t0") {
      profile_trans_value <- log(profile_value)
      param_index <- 3
    } else if(param_to_profile == "c_value2") {
      profile_trans_value <- logit(profile_value)
      param_index <- 4
    }
    
    # Get other parameter indices
    other_indices <- setdiff(1:4, param_index)
    start_par <- mle_trans_params[other_indices]
    
    # Optimize the other parameters while fixing the profiled one
    profile_obj_fn <- function(par, fixed_param_index, fixed_value) {
      full_par <- numeric(4)
      full_par[fixed_param_index] <- fixed_value
      full_par[other_indices] <- par
      return(f4_2epi_optim_reduced.NPI(full_par, sim_data_large, sim_data_small))
    }
    
    opt_result <- tryCatch({
      optim(
        par = start_par,
        fn = profile_obj_fn,
        fixed_param_index = param_index,
        fixed_value = profile_trans_value,
        method = "Nelder-Mead",
        control = list(maxit = 1000)
      )
    }, error = function(e) {
      list(value = NA, par = rep(NA, length(start_par)))
    })
    
    # Store results
    profile_results$neg_loglik[i] <- opt_result$value
    
    # Store optimized parameters
    opt_full_par <- numeric(4)
    opt_full_par[param_index] <- profile_trans_value
    opt_full_par[other_indices] <- opt_result$par
    
    # Transform back to original scale
    profile_results$R0[i] <- exp(opt_full_par[1])
    profile_results$v[i] <- exp(opt_full_par[2])
    profile_results$t0[i] <- exp(opt_full_par[3])
    profile_results$c_value2[i] <- expit(opt_full_par[4])
  }
  
  # Remove NA values
  profile_results <- profile_results[!is.na(profile_results$neg_loglik), ]
  
  # Calculate confidence intervals
  min_neg_loglik <- min(profile_results$neg_loglik)
  profile_results$LR_stat <- 2 * (profile_results$neg_loglik - min_neg_loglik)
  conf_threshold <- qchisq(0.95, df = 1)
  
  ci_data <- profile_results[profile_results$LR_stat <= conf_threshold, ]
  if(nrow(ci_data) > 0) {
    ci_lower <- min(ci_data$param_value)
    ci_upper <- max(ci_data$param_value)
  } else {
    ci_lower <- min(profile_results$param_value)
    ci_upper <- max(profile_results$param_value)
  }
  
  # Extract values for plotting
  mle_x <- mle_params[param_to_profile]
  true_x <- true_values[[param_to_profile]]
  
  # Create plots if requested
  if(plot) {
    
    
    
    param_labels <- list(
      R0 = expression(R[0]),           # R with subscript 0
      v = expression(nu),              # Greek letter nu
      t0 = expression(t[0]),           # t with subscript 0
      c_value2 = expression(c[1])      # c with subscript 1
    )
    
    # Create customizable title
    if(!is.null(custom_title)) {
      profile_title <- custom_title
    } else {
      epidemic_type <- "Two-Epidemic"
      npi_text <- if(exists("c_value2_spec")) paste0("NPI = ", c_value2_spec) else ""
      profile_title <- paste("Profile for", param_labels[[param_to_profile]], 
                             paste0("(", epidemic_type, if(npi_text != "") paste0(", ", npi_text), ")"))
    }
    
    # Create the profile likelihood plot with enhanced visibility
    p <- ggplot(profile_results, aes(x = param_value, y = -neg_loglik)) +
      geom_line(color = "steelblue", size = 1.2) +
      #geom_point(size = 2.5, color = "darkblue", alpha = 0.8) +
      geom_point(color = "darkblue", size = 2) +
      geom_hline(yintercept = -min_neg_loglik - conf_threshold/2, 
                 linetype = "dashed", color = "red", size = 1.2, alpha = 0.7) +
      geom_vline(xintercept = mle_x, color = "purple", linetype = "solid", size = 1.1) +
      geom_vline(xintercept = true_x, color = "green4", linetype = "dotted", size = 1.5) +
      geom_vline(xintercept = ci_lower, color = "blue", linetype = "dotted", size = 1.2) +
      geom_vline(xintercept = ci_upper, color = "blue", linetype = "dotted", size = 1.2) +
      annotate("text", x = mle_x, y = max(-profile_results$neg_loglik), 
               label = "MLE", color = "purple", hjust = -0.3, size = 4, fontface = "bold") +
      #annotate("text", x = true_x, y = max(-profile_results$neg_loglik) - 0.8, 
      # label = "True", color = "green4", hjust = -0.3, size = 4, fontface = "bold") +
      labs(#title = profile_title,
        #subtitle = paste0("95% CI: [", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
        x = param_labels[[param_to_profile]],
        y = "Log Likelihood") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
      )
    
    
    print(p)
  }
  
  # Return results
  return(list(
    profile_results = profile_results,
    profile_plot = if(plot) p else NULL,
    mle = mle_params,
    ci = c(ci_lower, ci_upper),
    ci_width = ci_upper - ci_lower
  ))
}























# ===============================================================================
# END OF MAXIMUM LIKELIHOOD FUNCTIONS
# ===============================================================================
#
# Summary of key functions provided:
#
# CORE FUNCTIONS:
# 1. poisson.loglik.NPI.reduced() - Calculate Poisson log-likelihood
# 2. f4_optim_reducedm.NPI() - Objective function for single epidemic
# 3. fit4_reducedm_loglik.NPI() - Fit model to single epidemic 
# 4. f4_2epi_optim_reduced.NPI() - Objective function for two epidemics
# 5. fit4_2epic_reduced.loglik.NPI() - Fit model to two epidemics
# 6. profile_likelihood_reducedm() - Profile likelihood for single epidemic
# 7. profile_likelihood_two_epidemics() - Profile likelihood for two epidemics
#
# GLOBAL VARIABLES EXPECTED:
# - times: Time vector for integration
# - initial_state: Initial conditions for single epidemic
# - initial_state_1, initial_state_2: Initial conditions for two epidemics
# - Various _spec parameters (t1_spec, t2_spec, rho_spec, etc.)
# - N: Population size
# 
# These functions support the reduced SEIR model with heterogeneity parameter v
# and provide comprehensive maximum likelihood estimation capabilities.
# ===============================================================================
