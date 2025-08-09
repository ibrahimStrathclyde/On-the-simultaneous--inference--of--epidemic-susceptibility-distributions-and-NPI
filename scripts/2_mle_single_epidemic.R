# ===============================================================================
# 2_mle_single_epidemic.R
# 
# Maximum Likelihood Estimation and Profile Likelihood Analysis for Single Epidemic
# 
# This script demonstrates parameter estimation for the reduced SEIR model with
# heterogeneity parameter v using data from a single epidemic. It performs
# comprehensive profile likelihood analysis to assess parameter identifiability
# and characterize parameter correlations.
#
# Key Analysis:
# - Maximum likelihood estimation of R_0, \nu, t_0, and c_1
# - Profile likelihood analysis for all parameters
# - Confidence interval calculation and coverage assessment
# - Parameter bias and precision analysis
# - Foundation for comparison with two-epidemic approach
# 
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# Date: August 2025
# ===============================================================================

# Load required libraries
library(tidyverse)
library(deSolve)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

# Source required functions
source("utility_functions.R")
source("MaxLik_fit_functions_reduced_model.R")

# ===============================================================================
# SECTION 1: SETUP AND PARAMETERS
# ===============================================================================

cat("===============================================================================\n")
cat("SINGLE EPIDEMIC MLE PROFILING ANALYSIS\n")
cat("===============================================================================\n")
cat("This script analyzes parameter identifiability using data from a single\n")
cat("epidemic to establish baseline parameter estimation performance.\n\n")

# Create output directories
if (!dir.exists("results")) {
  dir.create("results")
}
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Model parameters
N <- 100000                     # Population size
alpha_gamma_shape <- 0.5        # Shape parameter for gamma distribution
CV_true <- 1/sqrt(alpha_gamma_shape)  # Coefficient of variation (approximately 1.414)

# SEIR model parameters
R0_spec <- 3.0                  # Basic reproduction number
delta_spec <- 1/5.5             # Incubation rate (1/incubation period)
rho_spec <- 0.5                 # Relative infectiousness in E compartment
gamma_spec <- 1/4               # Recovery rate (1/infectious period)

# Intervention parameters
t0_spec <- 15                   # Time when behavioral changes begin
t1_spec <- 20                   # Time when full intervention starts
t2_spec <- 99                   # Time when intervention ends
t3_spec <- t2_spec + 1          # Time when transmission returns to baseline
tfinal_spec <- t3_spec          # Final simulation time
c_value1_spec <- 1              # Baseline transmission factor
c_value2_spec <- 0.3            # Intervention strength (70% reduction)
c_value3_spec <- 1              # Post-intervention transmission factor

# Single epidemic initial conditions
I0_fixed <- 40                  # Initial infected individuals
E0_fixed <- I0_fixed * 2.5      # Initial exposed individuals

cat("Model Parameters:\n")
cat("- Population size (N):", N, "\n")
cat("- True R0:", R0_spec, "\n")
cat("- True CV (v):", round(CV_true, 3), "\n")
cat("- Intervention strength:", c_value2_spec, "(", (1-c_value2_spec)*100, "% reduction)\n")
cat("- Initial infected (I0):", I0_fixed, "\n")
cat("- Initial exposed (E0):", E0_fixed, "\n\n")

# Set global variables needed by MaxLik functions
assign("N", N, envir = .GlobalEnv)
assign("rho_spec", rho_spec, envir = .GlobalEnv)
assign("delta_spec", delta_spec, envir = .GlobalEnv)
assign("gamma_spec", gamma_spec, envir = .GlobalEnv)
assign("t1_spec", t1_spec, envir = .GlobalEnv)
assign("t2_spec", t2_spec, envir = .GlobalEnv)
assign("t3_spec", t3_spec, envir = .GlobalEnv)
assign("c_value1_spec", c_value1_spec, envir = .GlobalEnv)
assign("c_value3_spec", c_value3_spec, envir = .GlobalEnv)
assign("tfinal_spec", tfinal_spec, envir = .GlobalEnv)

# Create initial state for single epidemic
initial_state <- c(S = N - E0_fixed - I0_fixed, E = E0_fixed, I = I0_fixed, R = 0, C = 0)

# Set global initial state
assign("initial_state", initial_state, envir = .GlobalEnv)

# Analysis settings
n_datasets <- 10                # Number of datasets to analyze for robustness
n_profile_points <- 15          # Number of points for profile likelihood

# Storage for results across datasets
all_ci_tables <- list()
all_profiles <- list()
successful_datasets <- 0

cat("Analysis Settings:\n")
cat("- Number of datasets:", n_datasets, "\n")
cat("- Profile points per parameter:", n_profile_points, "\n\n")

# ===============================================================================
# SECTION 2: MAIN ANALYSIS LOOP
# ===============================================================================

cat("Starting analysis of", n_datasets, "simulated datasets...\n\n")

for (dataset_id in 1:n_datasets) {
  cat("Analyzing dataset", dataset_id, "of", n_datasets, "\n")
  
  # Generate epidemic dataset with different random seed
  set.seed(1000 + dataset_id)
  
  # Simulate single epidemic
  sim_result <- simulate_cases_reduced_model(
    R0 = R0_spec,
    delta = delta_spec,
    rho = rho_spec,
    gamma = gamma_spec,
    v = CV_true,
    N = N,
    E0 = E0_fixed,
    I0 = I0_fixed,
    t0 = t0_spec,
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = c_value2_spec,
    c_value3 = c_value3_spec,
    tfinal = tfinal_spec
  )
  
  sim_data <- sim_result$sim_data
  
  cat("  Generated epidemic data (", sum(sim_data$reports), "total cases)\n")
  
  # Set global times variable (required by MaxLik functions)
  times <- sim_data$time
  assign("times", times, envir = .GlobalEnv)
  
  # ===============================================================================
  # SUBSECTION 2.1: Profile Likelihood Analysis
  # ===============================================================================
  
  cat("  Running profile likelihood analysis...\n")
  
  # Profile likelihood for coefficient of variation (v)
  cat("    - Profiling coefficient of variation (v)...\n")
  cv_profile <- tryCatch({
    profile_likelihood_reducedm(
      sim_data = sim_data,
      param_to_profile = "v",
      n_points = n_profile_points,
      plot = T  #  plot individual datasets
    )
  }, error = function(e) {
    cat("    Error in CV profile:", e$message, "\n")
    return(NULL)
  })
  
  # Profile likelihood for intervention strength (c_value2)
  cat("    - Profiling intervention strength (c_value2)...\n")
  c_profile <- tryCatch({
    profile_likelihood_reducedm(
      sim_data = sim_data,
      param_to_profile = "c_value2",
      n_points = n_profile_points,
      plot = T
    )
  }, error = function(e) {
    cat("    Error in c_value2 profile:", e$message, "\n")
    return(NULL)
  })
  
  # Profile likelihood for basic reproduction number (R0)
  cat("    - Profiling basic reproduction number (R0)...\n")
  r0_profile <- tryCatch({
    profile_likelihood_reducedm(
      sim_data = sim_data,
      param_to_profile = "R0",
      n_points = n_profile_points,
      plot = T
    )
  }, error = function(e) {
    cat("    Error in R0 profile:", e$message, "\n")
    return(NULL)
  })
  
  # Profile likelihood for intervention timing (t0)
  cat("    - Profiling intervention timing (t0)...\n")
  t0_profile <- tryCatch({
    profile_likelihood_reducedm(
      sim_data = sim_data,
      param_to_profile = "t0",
      n_points = n_profile_points,
      plot = T
    )
  }, error = function(e) {
    cat("    Error in t0 profile:", e$message, "\n")
    return(NULL)
  })
  
  # ===============================================================================
  # SUBSECTION 2.2: Store Results
  # ===============================================================================
  
  # Check if all profiles were successful
  profiles_successful <- !is.null(cv_profile) && !is.null(c_profile) && 
    !is.null(r0_profile) && !is.null(t0_profile)
  
  if (profiles_successful) {
    successful_datasets <- successful_datasets + 1
    
    # Create confidence interval table
    ci_table <- data.frame(
      Dataset = dataset_id,
      Parameter = c("R0", "V", "t0", "c"),
      True_Value = c(R0_spec, CV_true, t0_spec, c_value2_spec),
      MLE = c(
        r0_profile$mle,
        cv_profile$mle,
        t0_profile$mle,
        c_profile$mle
      ),
      Lower_CI = c(
        r0_profile$ci[1],
        cv_profile$ci[1],
        t0_profile$ci[1],
        c_profile$ci[1]
      ),
      Upper_CI = c(
        r0_profile$ci[2],
        cv_profile$ci[2],
        t0_profile$ci[2],
        c_profile$ci[2]
      )
    )
    
    # Calculate CI metrics
    ci_table$CI_Width <- ci_table$Upper_CI - ci_table$Lower_CI
    ci_table$Relative_Width <- (ci_table$CI_Width / ci_table$True_Value) * 100
    ci_table$Contains_True <- (ci_table$True_Value >= ci_table$Lower_CI) & 
      (ci_table$True_Value <= ci_table$Upper_CI)
    
    # Store results
    all_ci_tables[[dataset_id]] <- ci_table
    all_profiles[[dataset_id]] <- list(
      cv = cv_profile$profile_results,
      c = c_profile$profile_results,
      r0 = r0_profile$profile_results,
      t0 = t0_profile$profile_results
    )
    
    cat("  -> Dataset", dataset_id, "completed successfully\n")
  } else {
    cat("  -> Dataset", dataset_id, "failed - skipping\n")
  }
}

# ===============================================================================
# SECTION 3: SUMMARY ANALYSIS AND RESULTS
# ===============================================================================

cat("\n===============================================================================\n")
cat("ANALYSIS SUMMARY\n")
cat("===============================================================================\n")
cat("Successfully analyzed", successful_datasets, "out of", n_datasets, "datasets\n\n")

if (successful_datasets > 0) {
  # Combine all CI tables
  all_ci_df <- do.call(rbind, all_ci_tables)
  
  # Calculate summary statistics
  ci_summary <- all_ci_df %>%
    group_by(Parameter) %>%
    summarise(
      Mean_MLE = mean(MLE, na.rm = TRUE),
      SD_MLE = sd(MLE, na.rm = TRUE),
      Bias = mean(MLE - True_Value, na.rm = TRUE),
      RMSE = sqrt(mean((MLE - True_Value)^2, na.rm = TRUE)),
      Mean_CI_Width = mean(CI_Width, na.rm = TRUE),
      Coverage = mean(Contains_True, na.rm = TRUE) * 100,
      .groups = 'drop'
    )
  
  cat("Parameter Estimation Summary:\n")
  print(ci_summary)
  cat("\n")
  
  # ===============================================================================
  # SECTION 4: SAVE RESULTS
  # ===============================================================================
  
  cat("Saving results...\n")
  
  # Save detailed CI results
  write.csv(all_ci_df, "results/single_epidemic_profile_ci_results.csv", row.names = FALSE)
  
  # Save summary statistics
  write.csv(ci_summary, "results/single_epidemic_summary_statistics.csv", row.names = FALSE)
  
  cat("Results saved to:\n")
  cat("- results/single_epidemic_profile_ci_results.csv\n")
  cat("- results/single_epidemic_summary_statistics.csv\n\n")
  
  # ===============================================================================
  # SECTION 5: VISUALIZATION
  # ===============================================================================
  
  cat("Creating visualizations...\n")
  
  # 5.1: Parameter Estimation Accuracy Plot
  #estimation_plot <- ggplot(all_ci_df, aes(x = Parameter, y = MLE - True_Value)) +
   # geom_boxplot(aes(fill = Parameter), alpha = 0.7) +
   # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    #labs(
    #  title = "Parameter Estimation Bias (Single Epidemic)",
    #  subtitle = paste("Based on", successful_datasets, "simulated datasets"),
      #x = "Parameter",
     # y = "Bias (MLE - True Value)",
   #  # fill = "Parameter"
    #) +
    #theme_minimal() +
   # theme(
    #  plot.title = element_text(size = 14, face = "bold"),
   #   legend.position = "none"
 #   )
  
  #print(estimation_plot)
 # ggsave("figures/single_epidemic_parameter_bias.png", estimation_plot, 
     #    width = 10, height = 6, dpi = 300)
  
  # 5.2: Confidence Interval Coverage Plot
  #coverage_plot <- ggplot(ci_summary, aes(x = Parameter, y = Coverage)) +
  #  geom_col(aes(fill = Parameter), alpha = 0.8) +
   # geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    #labs(
    #  title = "95% Confidence Interval Coverage (Single Epidemic)",
    #  x = "Parameter",
     # y = "Coverage (%)",
    #  fill = "Parameter"
    #) +
    #ylim(0, 100) +
    #theme_minimal() +
   # theme(
    #  plot.title = element_text(size = 14, face = "bold"),
    #  legend.position = "none"
   # )
  #
 # print(coverage_plot)
  #ggsave("figures/single_epidemic_ci_coverage.png", coverage_plot, 
     #    width = 10, height = 6, dpi = 300)
  
  # 5.3: CI Width Analysis Plot
  #width_plot <- ggplot(all_ci_df, aes(x = Parameter, y = Relative_Width)) +
   # geom_boxplot(aes(fill = Parameter), alpha = 0.7) +
   # labs(
    #  title = "Confidence Interval Width (Single Epidemic)",
     # subtitle = "Relative to true parameter value",
    #  x = "Parameter",
      #y = "CI Width (% of true value)",
     # fill = "Parameter"
   # ) +
   # theme_minimal() +
  #  theme(
   #   plot.title = element_text(size = 14, face = "bold"),
     # legend.position = "none"
   # )
  
  #print(width_plot)
 # ggsave("figures/single_epidemic_ci_width.png", width_plot, 
   #      width = 10, height = 6, dpi = 300)
  
  # 5.4: MLE Distribution Plot
 # mle_distribution_plot <- ggplot(all_ci_df, aes(x = MLE, fill = Parameter)) +
 #   geom_histogram(alpha = 0.7, bins = 8) +
  #  geom_vline(aes(xintercept = True_Value), color = "red", linetype = "dashed", size = 1) +
  #  facet_wrap(~ Parameter, scales = "free") +
   # labs(
    #  title = "Distribution of MLE Estimates (Single Epidemic)",
    #  subtitle = "Red dashed line shows true value",
     # x = "MLE Estimate",
     # y = "Frequency",
    #  fill = "Parameter"
    #) +
   # theme_minimal() +
    #theme(
     # plot.title = element_text(size = 14, face = "bold"),
    #  legend.position = "none"
   # )
  
  #print(mle_distribution_plot)
 # ggsave("figures/single_epidemic_mle_distribution.png", mle_distribution_plot, 
       #  width = 12, height = 8, dpi = 300)
  
 # cat("Plots saved to figures/ directory\n\n")
  
  
} else {
  cat("No successful datasets to analyze!\n")
}



cat("===============================================================================\n")
cat("SINGLE EPIDEMIC ANALYSIS COMPLETED\n")
cat("===============================================================================\n")
cat("Key findings:\n")
if (successful_datasets > 0) {
  cat("- Successfully analyzed", successful_datasets, "datasets\n")
  cat("- Mean CI coverage across parameters:", round(mean(ci_summary$Coverage), 1), "%\n")
  cat("- Parameter identifiability varies across parameters\n")
  cat("- Results provide baseline for comparison with two-epidemic approach\n")
} else {
  cat("- Analysis failed - check parameter settings and functions\n")
}
cat("- See results/ and figures/ directories for detailed outputs\n")
cat("- Run 3_mle_two_epidemics.R next to compare approaches\n")

# ===============================================================================
# END OF SCRIPT
# ===============================================================================
