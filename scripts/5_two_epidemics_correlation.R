# Simple Two Epidemic Correlation Investigation

# Improvements: Better naming, NaN handling, cleaner output organization

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
n_replicates <- 200         # Number of replicates per initial condition
N <- 100000                 # Population size
alpha_gamma_shape <- 0.5    # Shape and rate of gamma distribution

# True parameter values
R0_spec <- 3
delta_spec <- 1/5.5
rho_spec <- 0.5
gamma_spec <- 1/4
CV_true <- sqrt(1/alpha_gamma_shape)  # 

# Intervention parameters (same for both epidemics)
t0_spec <- 15        # When people start distancing
t1_spec <- 20        # When lockdown begins
t2_spec <- 99        # When lockdown ends
t3_spec <- t2_spec + 1
tfinal_spec <- t3_spec
c_value1_spec <- 1
c_value2_spec <- 0.2  # Main NPI strength (change this to test different strengths)
c_value3_spec <- 1

# Two-epidemic setup
I0_small <- 40              # Fixed initial infected for small epidemic
E0_small <- I0_small * 2.5  # Fixed initial exposed for small epidemic

# Initial condition values for the auxiliary epidemic
I0_values <- c( 20,40,80, 160, 320, 400)

# Create output folder with descriptive name
model_name <- "two_epidemic"
npi_label <- paste0("npi", c_value2_spec)
output_folder <- paste0(model_name, "_", npi_label)
dir.create(output_folder, showWarnings = FALSE)

cat("=============================================================================\n")
cat("TWO EPIDEMIC CORRELATION ANALYSIS\n")
cat("=============================================================================\n")
cat("Model:", model_name, "\n")
cat("NPI strength:", c_value2_spec, "\n")
cat("Small epidemic: I0 =", I0_small, ", E0 =", E0_small, "(fixed)\n")
cat("Large epidemic I0 values:", paste(I0_values, collapse = ", "), "\n")
cat("Replicates per condition:", n_replicates, "\n")
cat("Output folder:", output_folder, "\n")
cat("=============================================================================\n")

# =============================================================================
# HELPER FUNCTIONS 
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

# Loop through different initial conditions for the larger epidemic
for (i in 1:length(I0_values)) {
  
  I0_large <- I0_values[i]
  E0_large <- I0_large * 2.5
  
  # Create initial states for both epidemics
  initial_state_1 <- c(S = N - E0_large - I0_large, E = E0_large, I = I0_large, R = 0, C = 0)
  initial_state_2 <- c(S = N - E0_small - I0_small, E = E0_small, I = I0_small, R = 0, C = 0)
  
  # Set these in the global environment for the fitting function to access
  assign("initial_state_1", initial_state_1, envir = .GlobalEnv)
  assign("initial_state_2", initial_state_2, envir = .GlobalEnv)
  
  cat(sprintf("\nSimulating with larger epidemic: I0 = %d, E0 = %d\n", I0_large, E0_large))
  cat(sprintf("Smaller epidemic fixed at: I0 = %d, E0 = %d\n", I0_small, E0_small))
  
  # Run simulations for multiple replicates
  for (j in 1:n_replicates) {
    if (j %% 50 == 0) {
      cat(sprintf("Processing replicate %d of %d for I0_large = %d\n", j, n_replicates, I0_large))
    }
    
    # Set different seeds for each epidemic to ensure independence
    set.seed(154327 + i * 1000 + j)
    
    # Simulate the larger epidemic
    sim_large <- simulate_cases_reduced_model(
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
    
    # Set different seed for the smaller epidemic
    set.seed(20000 + i * 1000 + j)
    
    # Simulate the smaller epidemic
    sim_small <- simulate_cases_reduced_model(
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
    
    # Extract simulated data
    sim_data_large <- sim_large$sim_data
    sim_data_small <- sim_small$sim_data
    
    # Set up global times variable
    times <- sim_data_large$time
    assign("times", times, envir = .GlobalEnv)
    
    # Fit the model using both epidemics
    z_mle <- tryCatch({
      fit4_2epic_reduced.loglik.NPI(dat1 = sim_data_large, dat2 = sim_data_small)
    }, error = function(e) {
      cat("Error in fitting replicate", j, ":", e$message, "\n")
      return(NULL)
    })
    
    # Skip this dataset if fitting failed
    if (is.null(z_mle)) {
      cat("Skipping replicate", j, "due to fitting error\n")
      next
    }
    
    # Set up results dataframe for this iteration
    i_results <- as.data.frame(matrix(z_mle$parms, nrow = 1))
    colnames(i_results) <- names(z_mle$parms)
    
    # Add epidemic information
    i_results$I0_large <- I0_large
    i_results$E0_large <- E0_large
    i_results$I0_small <- I0_small
    i_results$E0_small <- E0_small
    
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
    group_by(I0_large) %>%
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
  
  print("Summary statistics by large epidemic initial condition:")
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
  
  hessian_corr_summary <- valid_results %>%
    group_by(I0_large) %>%
    summarize(
      n_valid = n(),
      R0_v_median = median(R0_v_trans_cor, na.rm = TRUE),
      R0_t0_median = median(R0_t0_trans_cor, na.rm = TRUE),
      v_c_median = median(v_c_value2_trans_cor, na.rm = TRUE),
      .groups = 'drop'
    )
  
  heatmap_data <- hessian_corr_summary %>%
    pivot_longer(cols = c(R0_v_median, R0_t0_median, v_c_median),
                 names_to = "correlation_type", 
                 values_to = "correlation_value") %>%
    mutate(
      correlation_type = case_when(
        correlation_type == "R0_v_median" ~ "R[0] - nu",
        correlation_type == "R0_t0_median" ~ "R[0] - t[0]", 
        correlation_type == "v_c_median" ~ "nu - c[1]"
      )
    )
  
  heatmap_plot <- ggplot(heatmap_data, 
                         aes(x = correlation_type, y = factor(I0_large), fill = correlation_value)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limit = c(-1, 1)) +
    geom_text(aes(label = sprintf("%.2f", correlation_value)), 
              color = "black", size = 4, fontface = "bold") +
    labs(
      title = "Hessian-based Parameter Correlations by Initial Condition",
      subtitle = paste("Two Epidemic Model (NPI strength =", c_value2_spec, ", small epidemic I0 =", I0_small, ")"),
      x = "Parameter Pair",
      y = "Large Epidemic Initial I0",
      fill = "Correlation"
    ) +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 14,face = "bold"),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  ggsave(file.path(output_folder, make_filename("hessian_correlation_heatmap")),
         heatmap_plot, width = 10, height = 8)
  
 
  # 2. CORRELATION BOXPLOT
  cat("Creating correlation boxplot...\n")
  
  boxplot_data <- valid_results %>%
    dplyr::select(I0_large, R0_v_trans_cor, R0_t0_trans_cor, v_c_value2_trans_cor) %>%
    pivot_longer(cols = c(R0_v_trans_cor, R0_t0_trans_cor, v_c_value2_trans_cor),
                 names_to = "correlation_type",
                 values_to = "correlation_value") %>%
    mutate(
      correlation_type = case_when(
        correlation_type == "R0_v_trans_cor" ~ "R[0] - nu",
        correlation_type == "R0_t0_trans_cor" ~ "R[0] - t[0]",
        correlation_type == "v_c_value2_trans_cor" ~ "nu - c[1]"
      )
    )
  
  boxplot <- ggplot(boxplot_data, 
                    aes(x = correlation_type, y = correlation_value, fill = factor(I0_large))) +
    geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    scale_fill_brewer(palette = "Set2", name = "Large Epidemic I0") +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    labs(
      title = "Distribution of Parameter Correlations",
      subtitle = paste("Two Epidemic Model (NPI strength =", c_value2_spec, ", small epidemic I0 =", I0_small, ")"),
      x = "Parameter Pair",
      y = "Correlation Value"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"))
  
  ggsave(file.path(output_folder, make_filename("correlation_boxplot")),
         boxplot, width = 12, height = 8)
  
  # Print plots for immediate inspection
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

#"\n=============================================================================\n")
#"ANALYSIS COMPLETE
#"\n=============================================================================\n")







# =============================================================================
#  CORRELATION COMPARISON PLOT ACROSS INTERVENTION STRENGTHS 
# =============================================================================
# --- 1. Load Libraries ---
library(tidyverse)
library(patchwork)
library(viridis)

# --- 2. Load Raw Data Files ---
# Note: Ensure these file paths are correct for your system.
single_npi02 <- read.csv("single_epidemic_npi0.2/single_epidemic_npi0.2_results.csv")
single_npi03 <- read.csv("single_epidemic_npi0.3/single_epidemic_npi0.3_results.csv")
single_npi04 <- read.csv("single_epidemic_npi0.4/single_epidemic_npi0.4_results.csv")

two_npi02 <- read.csv("two_epidemic_npi0.2/two_epidemic_npi0.2_results.csv")
two_npi03 <- read.csv("two_epidemic_npi0.3/two_epidemic_npi0.3_results.csv")
two_npi04 <- read.csv("two_epidemic_npi0.4/two_epidemic_npi0.4_results.csv")


single_npi02$NPI <- 0.2   #
single_npi03$NPI <- 0.3
single_npi04$NPI <- 0.4

two_npi02$NPI <- 0.2
two_npi03$NPI <- 0.3
two_npi04$NPI <- 0.4



# Function to process single epidemic data
process_single_epidemic <- function(data, npi_value) {
  data %>%
    filter(convergence == 0 & hess_pd == 1 ) %>%
    mutate(
      I0 = I0,
      NPI = npi_value,
      Model = "Single Epidemic"
    ) %>%
    dplyr::select(I0, NPI, Model,
                  R0_v_trans_cor, R0_t0_trans_cor, v_c_value2_trans_cor,
                  dataset_id)
}

# Function to process two epidemic data
process_two_epidemic <- function(data, npi_value) {
  data %>%
    filter(convergence == 0 & hess_pd == 1) %>%
    mutate(
      I0 = I0_large, 
      NPI = npi_value,
      Model = "Two Epidemics"
    ) %>%
    dplyr:: select(I0, NPI, Model,
                   R0_v_trans_cor, R0_t0_trans_cor, v_c_value2_trans_cor,
                   dataset_id)
}

# Apply the processing functions
single_02 <- process_single_epidemic(single_npi02, 0.2)
single_03 <- process_single_epidemic(single_npi03, 0.3)
single_04 <- process_single_epidemic(single_npi04, 0.4)

two_02 <- process_two_epidemic(two_npi02, 0.2)
two_03 <- process_two_epidemic(two_npi03, 0.3)
two_04 <- process_two_epidemic(two_npi04, 0.4)

# Combine all processed data
all_data <- bind_rows(
  single_02, single_03, single_04,
  two_02, two_03, two_04
)

# Filter for valid, non-NA results
valid_data <- all_data %>%
  filter(!is.na(R0_v_trans_cor) & !is.na(I0))


# 4. Prepare Data for Plotting ---
# This section reshapes the data processed above for the separate plots.

plot_data <- valid_data %>%
  pivot_longer(
    cols = c("R0_v_trans_cor", "R0_t0_trans_cor", "v_c_value2_trans_cor"),
    names_to = "correlation_type",
    values_to = "correlation_value"
  ) %>%
  mutate(
    # Use plotmath expressions for labels
    correlation_type_parsed = factor(
      correlation_type,
      levels = c("R0_v_trans_cor", "R0_t0_trans_cor", "v_c_value2_trans_cor"),
      labels = c("R[0]~-~nu", "R[0]~-~t[0]", "nu~-~c[1]")
    ),
    # Add a column to identify the baseline case for highlighting
    baseline_highlight = ifelse(I0 == 40, "Baseline (I0 = 40)", "Other Cases")
  )

# Create final data subsets for each model. specify NPI value to print NPI==0.2,0.3,0.4
data_single_epidemic <- plot_data %>% filter(Model == "Single Epidemic", NPI==0.3) 
data_two_epidemics <- plot_data %>% filter(Model == "Two Epidemics", NPI==0.3)


# --- 5. Generate Plots for "Single Epidemic" Model 

# Boxplot for Single Epidemic Model with Highlighting
p_boxplot_single <- ggplot(data_single_epidemic, aes(x = factor(I0), y = correlation_value, fill = baseline_highlight)) +
  geom_boxplot(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~correlation_type_parsed, ncol = 3, labeller = label_parsed) +
  scale_fill_manual(
    name = "Comparison Group",
    values = c("Baseline (I0 = 40)" = "steelblue", "Other Cases" = "#FF7F00")
  ) +
  labs(
    title = "Parameter Correlations for Single Epidemic Model",
    subtitle = "Baseline case I0=40 highlighted",
    x = expression(paste("Initial ", I(0))),
    y = "Correlation Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size=14,angle = 45, hjust = 1),
    legend.position = "none"
  )


print(p_boxplot_single)

# Boxplot for Two Epidemics Model with Highlighting
p_boxplot_two <- ggplot(data_two_epidemics, aes(x = factor(I0), y = correlation_value, fill = baseline_highlight)) +
  geom_boxplot(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~correlation_type_parsed, ncol = 3, labeller = label_parsed) +
  scale_fill_manual(
    name = "Comparison Group",
    values = c("Baseline (I0 = 40)" = "steelblue", "Other Cases" = "orange")
  ) +
  labs(
    title = "Parameter Correlations for Two Epidemics Model",
     subtitle = "Baseline case I0=40 highlighted",
    x = expression(paste("Initial ", I(0), " of the auxiliary epidemic")),
    y = "Correlation Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size=14,angle = 45, hjust = 1),
    legend.position = "none"
  )


print(p_boxplot_two)

# To save the plots, use ggsave()
# ggsave("single_epidemic_boxplots_highlighted.png", plot = p_boxplot_single, width = 10, height = 6)
# ggsave("two_epidemics_boxplots_highlighted.png", plot = p_boxplot_two, width = 10, height = 6)