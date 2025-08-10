# Simple Single Epidemic Correlation Investigation (cleaned up)


library(tidyverse)
library(GGally)
library(deSolve)
library(gridExtra)
library(grid)
library(knitr)

# Source required functions
source("R/MaxLik_functions_reduced_model.R")
source("R/utility_functions.R")

# =============================================================================
# PARAMETERS AND SETTINGS
# =============================================================================

# Analysis settings
n_replicates <- 200           # Number of replicates per initial condition
K <- 20                      # Number of groups for gamma discretization
N <- 100000                  # Population size
alpha_gamma_shape <- 0.5     # Shape and rate of gamma distribution

# True parameter values
R0_spec <- 3
delta_spec <- 1/5.5
rho_spec <- 0.5
gamma_spec <- 1/4
CV_true <- sqrt(1/alpha_gamma_shape)  # ≈ 1.414

# Intervention parameters
t0_spec <- 15        # When people start distancing
t1_spec <- 20        # When lockdown begins  
t2_spec <- 99        # When lockdown ends
t3_spec <- t2_spec + 1
tfinal_spec <- t3_spec
c_value1_spec <- 1
c_value2_spec <- 0.3 # Main NPI strength (change this to test different strengths)
c_value3_spec <- 1

# Initial condition values to test
I0_values <- c(20,40, 80, 160, 320, 400)

# Create output folder with descriptive name
model_name <- "single_epidemic"
npi_label <- paste0("npi", c_value2_spec)
output_folder <- paste0(model_name, "_", npi_label)
dir.create(output_folder, showWarnings = FALSE)

cat("=============================================================================\n")
cat("SINGLE EPIDEMIC CORRELATION ANALYSIS\n")
cat("=============================================================================\n")
cat("Model:", model_name, "\n")
cat("NPI strength:", c_value2_spec, "\n")
cat("Initial conditions:", paste(I0_values, collapse = ", "), "\n")
cat("Replicates per condition:", n_replicates, "\n")
cat("Output folder:", output_folder, "\n")
cat("=============================================================================\n")

# =============================================================================
# HELPER FUNCTIONS (simple versions)
# =============================================================================

# Safe correlation function (handles NaN/empty data)
safe_cor <- function(x, y) {
  if(length(x) < 3 || length(y) < 3 || all(is.na(x)) || all(is.na(y))) {
    return(NA)
  }
  tryCatch(cor(x, y, use = "pairwise.complete.obs"), error = function(e) NA)
}

# Generate consistent filenames
make_filename <- function(plot_type, extension = "png") {
  paste0(model_name, "_", npi_label, "_", plot_type, ".", extension)
}

# =============================================================================
# MAIN ANALYSIS LOOP
# =============================================================================

# Initialize empty results dataframe
results <- data.frame()

# Loop through different initial conditions
for (i in 1:length(I0_values)) {
  
  I0_1 <- I0_values[i]
  E0_1 <- I0_1 * 2.5
  
  # Create initial state for the current initial conditions
  initial_state <- c(S = N - E0_1 - I0_1, E = E0_1, I = I0_1, R = 0, C = 0)
  assign("initial_state", initial_state, envir = .GlobalEnv)
  
  cat(sprintf("\nSimulating with initial conditions: I0 = %d, E0 = %d\n", I0_1, E0_1))
  
  # Set seed for reproducibility
  set.seed(12375 + i)
  
  # Simulate data for multiple datasets with current initial conditions
  sim_data <- NULL
  
  for (j in 1:n_replicates) {
    if (j %% 50 == 0) {
      cat(sprintf("Simulating dataset %d of %d for I0 = %d\n", j, n_replicates, I0_1))
    }
    
    # Simulate cases using heterogeneous model
    dsimdiscre <- simulate_cases_reduced_model(
      R0 = R0_spec,
      delta = delta_spec,
      rho = rho_spec, 
      gamma = gamma_spec,
      v = CV_true,
      N = N,
      E0 = E0_1,
      I0 = I0_1,
      t0 = t0_spec,
      t1 = t1_spec,
      t2 = t2_spec,
      t3 = t3_spec,
      c_value1 = c_value1_spec,
      c_value2 = c_value2_spec,
      c_value3 = c_value3_spec,
      tfinal = tfinal_spec
    ) 
    
    sim.data_1 <- dsimdiscre$sim_data
    sim.data_1$data_set <- j
    sim.data_1$I0 <- I0_1
    sim.data_1$E0 <- E0_1
    
    sim_data <- if (exists("sim_data") && !is.null(sim_data)) 
      bind_rows(sim_data, sim.data_1) else sim.data_1
  }
  
  # Set up global times variable
  times <- sim.data_1$time
  assign("times", times, envir = .GlobalEnv)
  
  # Fit the reduced model to each dataset for current initial conditions
  for (j in 1:n_replicates) { 
    if (j %% 50 == 0) {
      cat(sprintf("Fitting dataset %d for I0 = %d\n", j, I0_1))
    }
    
    # Extract current dataset
    sim.data_1 <- sim_data %>% 
      filter(data_set == j, I0 == I0_1) %>% 
      dplyr::select(-data_set, -I0, -E0)
    
    # Fit the reduced model with error handling
    z_mle <- tryCatch({
      fit4_reducedm_loglik.NPI(dat = sim.data_1)
    }, error = function(e) {
      cat("Error in fitting dataset", j, ":", e$message, "\n")
      return(NULL)
    })
    
    # Skip this dataset if fitting failed
    if (is.null(z_mle)) {
      cat("Skipping dataset", j, "due to fitting error\n")
      next
    }
    
    # Set up results dataframe for this iteration
    i_results <- as.data.frame(matrix(z_mle$parms, nrow = 1))
    colnames(i_results) <- names(z_mle$parms)
    
    # Add I0 and E0 columns
    i_results$I0 <- I0_1
    i_results$E0 <- E0_1
    
    # Initialize values for Hessian results
    z_se <- numeric(length(z_mle$trans_parms))
    z_cor <- c(0, 0, 0)
    z_hess <- 0
    z_pd <- 0
    z_ratio <- 0
    
    # Process Hessian matrix if available
    if (!is.null(z_mle$trans_hessian)) {
      tryCatch({
        z_hess <- 1
        z_eigen <- eigen(z_mle$trans_hessian)
        z_ratios <- z_eigen$values[1] / z_eigen$values
        z_ratio <- z_ratios[length(z_ratios)]
        
        if (all(z_eigen$values > 0)) {
          z_pd <- 1
          z_variance <- solve(z_mle$trans_hessian)
          z_d <- diag(1 / sqrt(diag(z_variance)), nrow = nrow(z_variance))
          z_correlation <- z_d %*% (z_variance %*% z_d)
          z_se <- sqrt(diag(z_variance))
          
          # Extract key correlations (R0-v, R0-t0, v-c_value2)
          z_cor <- c(
            z_correlation[2, 1],  # R0_v_trans_cor
            z_correlation[3, 1],  # R0_t0_trans_cor
            z_correlation[4, 2]   # v_c_value2_trans_cor
          )
          
          # Calculate confidence intervals
          par_ucl <- z_mle$trans_parms + 1.96 * sqrt(diag(z_variance))
          par_lcl <- z_mle$trans_parms - 1.96 * sqrt(diag(z_variance))
          
          C_intervals <- as.data.frame(matrix(c(
            exp(par_lcl[1]), exp(par_ucl[1]),          # R0
            exp(par_lcl[2]), exp(par_ucl[2]),          # v
            exp(par_lcl[3]), exp(par_ucl[3]),          # t0
            expit(par_lcl[4]), expit(par_ucl[4])       # c_value2
          ), nrow = 1, byrow = TRUE))
          
          colnames(C_intervals) <- c(
            "R0_lcl", "R0_ucl", 
            "v_lcl", "v_ucl", 
            "t0_lcl", "t0_ucl", 
            "c_value2_lcl", "c_value2_ucl"
          )
        }
      }, error = function(e) {
        # Keep default values if Hessian processing fails
      })
    }
    
    # Create dataframes for transformed parameters
    z1 <- as.data.frame(matrix(z_mle$trans_parms, nrow = 1))
    colnames(z1) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans", sep = "_")
    
    z_se <- as.data.frame(matrix(z_se, nrow = 1))
    colnames(z_se) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans_se", sep = "_")
    
    z_cor <- as.data.frame(matrix(z_cor, nrow = 1))
    colnames(z_cor) <- c("R0_v_trans_cor", "R0_t0_trans_cor", "v_c_value2_trans_cor")
    
    # Combine all results for this iteration
    i_results <- i_results %>% bind_cols(z1, z_se, z_cor)
    i_results <- if (exists("C_intervals")) bind_cols(i_results, C_intervals) else i_results
    
    # Add diagnostic values and dataset ID
    i_results$hess_exists <- z_hess
    i_results$hess_pd <- z_pd
    i_results$ratio_max_min_evalue <- z_ratio
    i_results$dataset_id <- j
    
    # Append to results
    results <- if (exists("results") && nrow(results) > 0) bind_rows(results, i_results) else i_results
  }
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save results with descriptive filename
results_file <- make_filename("results", "csv")
write.csv(results, file.path(output_folder, results_file), row.names = FALSE)
cat("\nResults saved to:", file.path(output_folder, results_file), "\n")

# Calculate and save summary statistics
if (nrow(results) > 0) {
  result_summary <- results %>% 
    group_by(I0) %>%
    summarize(
      n_total = n(),
      n_converged = sum(convergence == 0, na.rm = TRUE),
      n_valid_hessian = sum(convergence == 0 & hess_pd == 1, na.rm = TRUE),
      success_rate = n_converged / n_total,
      
      # Parameter estimates (only converged results)
      R0_mean = mean(R0[convergence == 0], na.rm = TRUE),
      R0_sd = sd(R0[convergence == 0], na.rm = TRUE),
      R0_bias = mean(R0[convergence == 0] - R0_spec, na.rm = TRUE),
      
      v_mean = mean(v[convergence == 0], na.rm = TRUE),
      v_sd = sd(v[convergence == 0], na.rm = TRUE),
      v_bias = mean(v[convergence == 0] - CV_true, na.rm = TRUE),
      
      t0_mean = mean(t0[convergence == 0], na.rm = TRUE),
      t0_sd = sd(t0[convergence == 0], na.rm = TRUE),
      t0_bias = mean(t0[convergence == 0] - t0_spec, na.rm = TRUE),
      
      c_value2_mean = mean(c_value2[convergence == 0], na.rm = TRUE),
      c_value2_sd = sd(c_value2[convergence == 0], na.rm = TRUE),
      c_value2_bias = mean(c_value2[convergence == 0] - c_value2_spec, na.rm = TRUE),
      
      # Hessian-based correlations (only positive definite Hessians)
      R0_v_corr_median = median(R0_v_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      R0_v_corr_mean = mean(R0_v_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      R0_v_corr_min = min(R0_v_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      R0_v_corr_max = max(R0_v_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      
      R0_t0_corr_median = median(R0_t0_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      R0_t0_corr_mean = mean(R0_t0_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      R0_t0_corr_min = min(R0_t0_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      R0_t0_corr_max = max(R0_t0_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      
      v_c_value2_corr_median = median(v_c_value2_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      v_c_value2_corr_mean = mean(v_c_value2_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      v_c_value2_corr_min = min(v_c_value2_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      v_c_value2_corr_max = max(v_c_value2_trans_cor[convergence == 0 & hess_pd == 1], na.rm = TRUE),
      
      .groups = 'drop'
    )
  
  # Save summary
 summary_file <- make_filename("summary", "csv")
  write.csv(result_summary, file.path(output_folder, summary_file), row.names = FALSE)
  
  print("Summary statistics by initial condition:")
  print(kable(result_summary, digits = 3))
}

# =============================================================================
# CREATE VISUALIZATIONS
# =============================================================================

# Filter valid results for plotting
valid_results <- results %>% 
  filter(convergence == 0 & hess_pd == 1 & !is.na(R0_v_trans_cor))

if (nrow(valid_results) > 0) {
  
  cat("\nCreating visualizations...\n")
  
  # 1. HESSIAN CORRELATION HEATMAP
  cat("Creating Hessian correlation heatmap...\n")
  
  # Calculate correlations by initial condition  
  hessian_corr_summary <- valid_results %>%
    group_by(I0) %>%
    summarize(
      n_valid = n(),
      R0_v_median = median(R0_v_trans_cor, na.rm = TRUE),
      R0_v_mean = mean(R0_v_trans_cor, na.rm = TRUE),
      R0_v_min = min(R0_v_trans_cor, na.rm = TRUE),
      R0_v_max = max(R0_v_trans_cor, na.rm = TRUE),
      
      R0_t0_median = median(R0_t0_trans_cor, na.rm = TRUE),
      R0_t0_mean = mean(R0_t0_trans_cor, na.rm = TRUE),
      R0_t0_min = min(R0_t0_trans_cor, na.rm = TRUE),
      R0_t0_max = max(R0_t0_trans_cor, na.rm = TRUE),
      
      v_c_median = median(v_c_value2_trans_cor, na.rm = TRUE),
      v_c_mean = mean(v_c_value2_trans_cor, na.rm = TRUE),
      v_c_min = min(v_c_value2_trans_cor, na.rm = TRUE),
      v_c_max = max(v_c_value2_trans_cor, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Create heatmap data using median values
  heatmap_data <- hessian_corr_summary %>%
    pivot_longer(cols = c(R0_v_median, R0_t0_median, v_c_median),
                 names_to = "correlation_type", 
                 values_to = "correlation_value") %>%
    mutate(
      correlation_type = case_when(
        correlation_type == "R0_v_median" ~ "R0-v",
        correlation_type == "R0_t0_median" ~ "R0-t0", 
        correlation_type == "v_c_median" ~ "v-c"
      )
    )
  
  # Create heatmap
  heatmap_plot <- ggplot(heatmap_data, 
                         aes(x = correlation_type, y = factor(I0), fill = correlation_value)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limit = c(-1, 1)) +
    geom_text(aes(label = sprintf("%.2f", correlation_value)), 
              color = "black", size = 4, fontface = "bold") +
    labs(
      title = "Hessian-based Parameter Correlations by Initial Condition",
      subtitle = paste("Single Epidemic Model (NPI strength =", c_value2_spec, ")"),
      x = "Parameter Pair",
      y = "Initial I0",
      fill = "Correlation"
    ) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold"))
  
  # Save heatmap
  heatmap_file <- make_filename("hessian_correlation_heatmap")
 ggsave(file.path(output_folder, heatmap_file), heatmap_plot, width = 10, height = 8)
  
  # 2. CORRELATION TRENDS PLOT
  cat("Creating correlation trends plot...\n")
  
  # Calculate absolute correlations for trends using median
  trends_data <- hessian_corr_summary %>%
    mutate(
      R0_v_abs = abs(R0_v_median),
      R0_t0_abs = abs(R0_t0_median),
      v_c_abs = abs(v_c_median)
    )
 
  # 3. CORRELATION BOXPLOT
  cat("Creating correlation boxplot...\n")
  
  # Prepare data for boxplot
  boxplot_data <- valid_results %>%
    dplyr::select(I0, R0_v_trans_cor, R0_t0_trans_cor, v_c_value2_trans_cor) %>%
    pivot_longer(cols = c(R0_v_trans_cor, R0_t0_trans_cor, v_c_value2_trans_cor),
                 names_to = "correlation_type",
                 values_to = "correlation_value") %>%
    mutate(
      correlation_type = case_when(
        correlation_type == "R0_v_trans_cor" ~ "R0-v",
        correlation_type == "R0_t0_trans_cor" ~ "R0-t0",
        correlation_type == "v_c_value2_trans_cor" ~ "v-c"
      )
    )
  
  boxplot <- ggplot(boxplot_data, 
                    aes(x = correlation_type, y = correlation_value, fill = factor(I0))) +
    geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    scale_fill_brewer(palette = "Set2", name = "Initial I0") +
    labs(
      title = "Distribution of Parameter Correlations",
      subtitle = paste("Single Epidemic Model (NPI strength =", c_value2_spec, ")"),
      x = "Parameter Pair",
      y = "Correlation Value"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12))
  
  ## Save boxplot
  boxplot_file <- make_filename("correlation_boxplot")
  ggsave(file.path(output_folder, boxplot_file), boxplot, width = 12, height = 8)
  
  print(heatmap_plot)
 
  print(boxplot)
  
} else {
  cat("No valid results with positive definite Hessian matrices found for visualization.\n")
}

# =============================================================================
# CREATE FINAL SUMMARY TABLE
# =============================================================================

# Create a comprehensive correlation summary
if (nrow(valid_results) > 0) {
  final_summary <- data.frame(
    Parameter_Pair = c("R0-v", "R0-t0", "v-c"),
    Median = c(
      median(valid_results$R0_v_trans_cor, na.rm = TRUE),
      median(valid_results$R0_t0_trans_cor, na.rm = TRUE),
      median(valid_results$v_c_value2_trans_cor, na.rm = TRUE)
    ),
    Mean = c(
      mean(valid_results$R0_v_trans_cor, na.rm = TRUE),
      mean(valid_results$R0_t0_trans_cor, na.rm = TRUE),
      mean(valid_results$v_c_value2_trans_cor, na.rm = TRUE)
    ),
    SD = c(
      sd(valid_results$R0_v_trans_cor, na.rm = TRUE),
      sd(valid_results$R0_t0_trans_cor, na.rm = TRUE),
      sd(valid_results$v_c_value2_trans_cor, na.rm = TRUE)
    ),
    Min = c(
      min(valid_results$R0_v_trans_cor, na.rm = TRUE),
      min(valid_results$R0_t0_trans_cor, na.rm = TRUE),
      min(valid_results$v_c_value2_trans_cor, na.rm = TRUE)
    ),
    Max = c(
      max(valid_results$R0_v_trans_cor, na.rm = TRUE),
      max(valid_results$R0_t0_trans_cor, na.rm = TRUE),
      max(valid_results$v_c_value2_trans_cor, na.rm = TRUE)
    ),
    Abs_Median = c(
      median(abs(valid_results$R0_v_trans_cor), na.rm = TRUE),
      median(abs(valid_results$R0_t0_trans_cor), na.rm = TRUE),
      median(abs(valid_results$v_c_value2_trans_cor), na.rm = TRUE)
    )
  )
  
  # Save final summary
 final_summary_file <- make_filename("correlation_summary", "csv")
  write.csv(final_summary, file.path(output_folder, final_summary_file), row.names = FALSE)
  
  cat("\nFinal Correlation Summary:\n")
  print(kable(final_summary, digits = 3, 
              caption = paste("Parameter Correlations -", model_name, npi_label)))
}

cat("\n=============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=============================================================================\n")
cat("Total results:", nrow(results), "\n")
cat("Valid results:", nrow(valid_results), "\n")
cat("Output folder:", output_folder, "\n")
cat("Key files created:\n")
cat("  -", results_file, "(full results)\n")
cat("  -", summary_file, "(summary by I0)\n")
cat("  -", final_summary_file, "(correlation summary)\n")
cat("  - Correlation heatmap, trends, and boxplot images\n")
cat("=============================================================================\n")


