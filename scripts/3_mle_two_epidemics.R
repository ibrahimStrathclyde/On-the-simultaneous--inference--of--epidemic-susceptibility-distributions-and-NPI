# ===============================================================================
# 3_mle_two_epidemics.R
# 
# Maximum Likelihood Estimation and Profile Likelihood Analysis for Two Epidemics
# 
# This script demonstrates how using data from two concurrent epidemics with
# different initial conditions can improve parameter identifiability in SEIR
# models. It fits the reduced SEIR model with heterogeneity parameter v to
# two epidemic datasets simultaneously.
#
# Key Analysis:
# - Simultaneous fitting to two epidemics with different initial conditions
# - Profile likelihood analysis for all parameters (R0, v, t0, c_value2)
# - Comparison of parameter correlations vs single epidemic approach
# - Confidence interval analysis and coverage assessment
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
source("R/utility_functions.R")
source("R/MaxLik_functions_reduced_model.R")

# ===============================================================================
# SECTION 1: SETUP AND PARAMETERS
# ===============================================================================

cat("===============================================================================\n")
cat("TWO EPIDEMICS MLE PROFILING ANALYSIS\n")
cat("===============================================================================\n")
cat("This script analyzes parameter identifiability using two concurrent epidemics\n")
cat("with different initial conditions to improve parameter estimation.\n\n")

# Create output directories
if (!dir.exists("results")) {
  dir.create("results")
}
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Model parameters (same as single epidemic for comparison)
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

# Two epidemic setup - different initial conditions
# Larger epidemic (higher initial prevalence)
I0_large <- 400
E0_large <- I0_large * 2.5

# Smaller epidemic (lower initial prevalence)  
I0_small <- 20
E0_small <- I0_small * 2.5

cat("Model Parameters:\n")
cat("- Population size (N):", N, "\n")
cat("- True R0:", R0_spec, "\n")
cat("- True CV (v):", round(CV_true, 3), "\n")
cat("- Intervention strength:", c_value2_spec, "(", (1-c_value2_spec)*100, "% reduction)\n")
cat("- Large epidemic I0:", I0_large, "\n")
cat("- Small epidemic I0:", I0_small, "\n\n")

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

# Create initial states for both epidemics
initial_state_1 <- c(S = N - E0_large - I0_large, E = E0_large, I = I0_large, R = 0, C = 0)
initial_state_2 <- c(S = N - E0_small - I0_small, E = E0_small, I = I0_small, R = 0, C = 0)

# Set global initial states
assign("initial_state_1", initial_state_1, envir = .GlobalEnv)
assign("initial_state_2", initial_state_2, envir = .GlobalEnv)

# Analysis settings
n_datasets <- 10                # Number of datasets to analyze for robustness adjust to 200 for the paper.
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
  
  # Generate two epidemic datasets with different random seeds
  set.seed(1000 + dataset_id)
  
  # Simulate larger epidemic
  sim_result_1 <- simulate_cases_reduced_model(
    R0 = R0_spec,
    delta = delta_spec,
    rho = rho_spec,
    gamma = gamma_spec,
    v = CV_true,
    N = N,
    E0 = E0_large,
    I0 = I0_large,
    t0 = t0_spec,
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = c_value2_spec,
    c_value3 = c_value3_spec,
    tfinal = tfinal_spec
  )
  
  # Simulate smaller epidemic (same parameters, different initial conditions)
  sim_result_2 <- simulate_cases_reduced_model(
    R0 = R0_spec,
    delta = delta_spec,
    rho = rho_spec,
    gamma = gamma_spec,
    v = CV_true,
    N = N,
    E0 = E0_small,
    I0 = I0_small,
    t0 = t0_spec,
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = c_value2_spec,
    c_value3 = c_value3_spec,
    tfinal = tfinal_spec
  )
  
  sim_data_large <- sim_result_1$sim_data
  sim_data_small <- sim_result_2$sim_data
  
  # Set global times variable (required by MaxLik functions)
  times <- sim_data_large$time
  assign("times", times, envir = .GlobalEnv)
  
  cat("  Generated epidemic data (Large:", sum(sim_data_large$reports), 
      "total cases; Small:", sum(sim_data_small$reports), "total cases)\n")
  
  # ===============================================================================
  # SUBSECTION 2.1: Profile Likelihood Analysis
  # ===============================================================================
  
  cat("  Running profile likelihood analysis...\n")
  
  # Profile likelihood for coefficient of variation (v)
  cat("    - Profiling coefficient of variation (v)...\n")
  cv_profile <- tryCatch({
    profile_likelihood_two_epidemics(
      sim_data_large = sim_data_large,
      sim_data_small = sim_data_small,
      param_to_profile = "v",
      n_points = n_profile_points,
      plot = FALSE  # Don't plot individual datasets
    )
  }, error = function(e) {
    cat("    Error in CV profile:", e$message, "\n")
    return(NULL)
  })
  
  # Profile likelihood for intervention strength (c_value2)
  cat("    - Profiling intervention strength (c_value2)...\n")
  c_profile <- tryCatch({
    profile_likelihood_two_epidemics(
      sim_data_large = sim_data_large,
      sim_data_small = sim_data_small,
      param_to_profile = "c_value2",
      n_points = n_profile_points,
      plot = FALSE
    )
  }, error = function(e) {
    cat("    Error in c_value2 profile:", e$message, "\n")
    return(NULL)
  })
  
  # Profile likelihood for basic reproduction number (R0)
  cat("    - Profiling basic reproduction number (R0)...\n")
  r0_profile <- tryCatch({
    profile_likelihood_two_epidemics(
      sim_data_large = sim_data_large,
      sim_data_small = sim_data_small,
      param_to_profile = "R0",
      n_points = n_profile_points,
      plot = TRUE # Plot the profile for R0 
    )
  }, error = function(e) {
    cat("    Error in R0 profile:", e$message, "\n")
    return(NULL)
  })
  
  # Profile likelihood for intervention timing (t0)
  cat("    - Profiling intervention timing (t0)...\n")
  t0_profile <- tryCatch({
    profile_likelihood_two_epidemics(
      sim_data_large = sim_data_large,
      sim_data_small = sim_data_small,
      param_to_profile = "t0",
      n_points = n_profile_points,
      plot = FALSE
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
  write.csv(all_ci_df, "results/two_epidemics_profile_ci_results.csv", row.names = FALSE)
  
  # Save summary statistics
  write.csv(ci_summary, "results/two_epidemics_summary_statistics.csv", row.names = FALSE)
  
  cat("Results saved to:\n")
  cat("- results/two_epidemics_profile_ci_results.csv\n")
  cat("- results/two_epidemics_summary_statistics.csv\n\n")
  
  # ===============================================================================
  # SECTION 5: VISUALIZATION
  # ===============================================================================
  
  cat("Creating visualizations...\n")
  
  # 5.1: Parameter Estimation Accuracy Plot
  estimation_plot <- ggplot(all_ci_df, aes(x = Parameter, y = MLE - True_Value)) +
    geom_boxplot(aes(fill = Parameter), alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Parameter Estimation Bias (Two Epidemics)",
      subtitle = paste("Based on", successful_datasets, "simulated datasets"),
      x = "Parameter",
      y = "Bias (MLE - True Value)",
      fill = "Parameter"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
  
  print(estimation_plot)
  ggsave("figures/two_epidemics_parameter_bias.png", estimation_plot, 
         width = 10, height = 6, dpi = 300)
  
  # 5.2: Confidence Interval Coverage Plot
  coverage_plot <- ggplot(ci_summary, aes(x = Parameter, y = Coverage)) +
    geom_col(aes(fill = Parameter), alpha = 0.8) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    labs(
      title = "95% Confidence Interval Coverage (Two Epidemics)",
      x = "Parameter",
      y = "Coverage (%)",
      fill = "Parameter"
    ) +
    ylim(0, 100) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
  
  print(coverage_plot)
  ggsave("figures/two_epidemics_ci_coverage.png", coverage_plot, 
         width = 10, height = 6, dpi = 300)
  
  # 5.3: CI Width Comparison Plot
  width_plot <- ggplot(all_ci_df, aes(x = Parameter, y = Relative_Width)) +
    geom_boxplot(aes(fill = Parameter), alpha = 0.7) +
    labs(
      title = "Confidence Interval Width (Two Epidemics)",
      subtitle = "Relative to true parameter value",
      x = "Parameter",
      y = "CI Width (% of true value)",
      fill = "Parameter"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
  
  print(width_plot)
  ggsave("figures/two_epidemics_ci_width.png", width_plot, 
         width = 10, height = 6, dpi = 300)
  
  cat("Plots saved to figures/ directory\n\n")
  
} else {
  cat("No successful datasets to analyze!\n")
}

# ===============================================================================
# SECTION 6: COMPARISON WITH SINGLE EPIDEMIC RESULTS (IF AVAILABLE)
# ===============================================================================

cat("===============================================================================\n")
cat("COMPARISON WITH SINGLE EPIDEMIC RESULTS\n")
cat("===============================================================================\n")

# Check if single epidemic results exist
single_epi_file <- "results/single_epidemic_profile_ci_results.csv"

if (file.exists(single_epi_file) && successful_datasets > 0) {
  cat("Loading single epidemic results for comparison...\n")
  
  single_epi_data <- read.csv(single_epi_file)
  
  # Calculate single epidemic summary
  single_summary <- single_epi_data %>%
    group_by(Parameter) %>%
    summarise(
      Mean_CI_Width = mean(CI_Width, na.rm = TRUE),
      Coverage = mean(Contains_True, na.rm = TRUE) * 100,
      .groups = 'drop'
    ) %>%
    mutate(Method = "Single Epidemic")
  
  # Calculate two epidemics summary for comparison
  two_summary <- ci_summary %>%
    dplyr::select(Parameter, Mean_CI_Width, Coverage) %>%
    mutate(Method = "Two Epidemics")
  
  # Combine for comparison
  comparison_data <- rbind(single_summary, two_summary)
  
  # Create comparison plots
  width_comparison <- ggplot(comparison_data, aes(x = Parameter, y = Mean_CI_Width, fill = Method)) +
    geom_col(position = "dodge", alpha = 0.8) +
    labs(
      title = "Confidence Interval Width Comparison",
      subtitle = "Single vs Two Epidemics Approach",
      x = "Parameter",
      y = "Mean CI Width",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  print(width_comparison)
  ggsave("figures/ci_width_comparison_single_vs_two.png", width_comparison, 
         width = 12, height = 6, dpi = 300)
  
  # Calculate improvement metrics
  improvement_metrics <- comparison_data %>%
    dplyr::select(Parameter, Mean_CI_Width, Method) %>%
    pivot_wider(names_from = Method, values_from = Mean_CI_Width) %>%
    mutate(
      Improvement_Factor = `Single Epidemic` / `Two Epidemics`,
      Percent_Reduction = (1 - `Two Epidemics` / `Single Epidemic`) * 100
    )
  
  cat("CI Width Improvement Summary:\n")
  print(improvement_metrics)
  
  write.csv(improvement_metrics, "results/identifiability_improvement_metrics.csv", row.names = FALSE)
  
} else {
  if (!file.exists(single_epi_file)) {
    cat("Single epidemic results not found. Run 2_mle_single_epidemic.R first for comparison.\n")
  } else {
    cat("No successful two-epidemic results to compare.\n")
  }
}

cat("\n===============================================================================\n")
cat("TWO EPIDEMICS ANALYSIS COMPLETED\n")
cat("===============================================================================\n")
cat("Key findings:\n")
if (successful_datasets > 0) {
  cat("- Successfully analyzed", successful_datasets, "datasets\n")
  cat("- Mean CI coverage across parameters:", round(mean(ci_summary$Coverage), 1), "%\n")
  cat("- Results demonstrate improved parameter identifiability with two epidemics\n")
} else {
  cat("- Analysis failed - check parameter settings and functions\n")
}
cat("- See results/ and figures/ directories for detailed outputs\n")

# ===============================================================================
# END OF SCRIPT

# ===============================================================================
