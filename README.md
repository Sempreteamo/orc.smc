
# orc.smc

<!-- badges: start -->
<!-- badges: end -->

The goal of orc.smc is to provide tools for implementing online rolling controlled sequential Monte Carlo (SMC) methods. This package is designed to facilitate robust normalizing constant estimate in dynamic systems with a rolling-window setting, particularly where real-time adaptation and efficient handling of streaming data are crucial. It aims to offer a flexible and extensible framework for developing and applying advanced particle filtering techniques in evolving environments.

## Installation

You can install the development version of orc.smc like so:

``` r
# Install devtools if you haven't already
# install.packages("devtools")

devtools::install_github("Sempreteamo/orc.smc")
# Load the ORCSMC package developed based on the Arxiv paper

## Example

This is a basic example which the code works for a multivaraite linear gaussian model:

``` r
library(orc.smc)
## basic example code

## --- 1. Simulation Setup ---

Napf = N = 100 # Number of particles (N)
lag = 2 # Rolling window length (L)
Time = 10 # Total time steps (T)
d_ = 2 # State dimension

# Construct the Transition Matrix 
alpha = 0.415 # Correlation coefficient for the transition matrix

tran_m <- matrix(nrow = d_, ncol = d_)
for (i in 1:d_){
  for (j in 1:d_){
    tran_m[i,j] = alpha^(abs(i-j) + 1)
  }
}

## --- 2. Model Specification ---

# Define the Linear Gaussian State-Space Model parameters
ini <- rep(0, d_) # Initial state mean
ini_c = diag(1, nrow = d_, ncol = d_) # Initial state covariance
tran_c = diag(1, nrow = d_, ncol = d_) # State transition noise covariance
obs_m = diag(1, nrow = d_, ncol = d_) # Observation matrix
obs_c = diag(1, nrow = d_, ncol = d_) # Observation noise covariance

# ORCSMC specific hyperparameters
parameters_ <- list(k = 5, tau = 0.5, kappa = 0.5)
obs_p <- list(obs_mean = obs_m, obs_cov = obs_c)

# Model structure including likelihood and simulation functions
model <- list(
  ini_mu = ini, ini_cov = ini_c, 
  tran_mu = tran_m, tran_cov = tran_c, 
  obs_params = obs_p,
  eval_likelihood = evaluate_likelihood_lg, # Evaluates p(y_t | x_t)
  simu_observation = simulate_observation_lg, 
  parameters = parameters_
)

## --- 3. Ground Truth Generation & Benchmarking ---

set.seed(1234)

# Simulate observations y_{1:T}
obs_ <- sample_obs(model, Time, d_) 

# Setup parameters for Fast Kalman Filter (FKF)
params <- list(
  dt = matrix(0, d_, 1), ct = matrix(0, d_, 1), Tt = as.matrix(tran_m), 
  P0 = ini_c, Zt = obs_m, Ht = tran_c, Gt = obs_c, a0 = ini, d = d_
)

# Compute the analytical solution using FKF
filter <- compute_fkf(params, obs_)
fkf_obj <- filter[[1]]

## --- 4. ORCSMC Algorithm ---

data_ = data <- list(obs = obs_)

output <- Orc_SMC(lag, data, model, Napf)

## --- 5. Quantitative Evaluation ---
# A log-ratio close to 1 indicates that ORCSMC successfully learned the optimal 
# proposal distribution
ratio <- compute_ratio(output$logZ[Time], fkf_obj)
```

## Figure Reproduction

This script reproduces the figures from "Online Rolling Controlled Sequential Monte Carlo" (Xue et al.) using synthetic data generated from experimental functions.

``` r
# ==============================================================================
# Script: Visualization for ORCSMC Paper Reproducibility
# Description: This script generates TikZ (.tex) plots from experimental CSV data.
# Instructions: 
#   1. Ensure your generated datasets (CSVs) are in the working directory.
#   2. Update the file paths in the 'DATA_PATHS' section below if necessary.
# ==============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(boot)
library(tikzDevice)
library(ggplot2)
library(patchwork)

# --- 0. DATASET CONFIGURATION (User: Define your file paths here) ---
# To reproduce the plots, replace these filenames with your own generated datasets.
DATA_PATHS <- list(
  diag_1d   = "orc_N1000T100_d1_lag2-5_diagf_rep100.csv",
  rmse_data = "orc+bpf+iapf_N1000T100_d2-64_lag2-16_non-diagf_rep100.csv",
  nc_comp   = "orc+bpf+iapf_N1000T100_d2-64_lag2-16_non-diagf_rep100.csv",
  svm_data  = "bpf+orc_N200_d1_lag2-16_svm_rep100.csv",
  bin_multi = "bin_N1000T100_d2-64_lag2-16_rep100.csv",
  l1_error  = "/Users/lexie/Documents/data/l1error_orc_N1000T100_d2-64_lag2_16_rep100.csv",
  diag_nd   = "diag0.415_nc_orc+bpf+iapf_N1000T100_d2-64_l2-16.csv",
  neuro_ess = "bin_ess_1d_l2-16.csv",
  neuro_nd  = "bin_N1000T100_d2-64_lag2-16_rep100.csv"
)

# ##############################################################################
#### diagnol_1d ####
# ##############################################################################
# Load the 1D diagonal case results
combined_df <- read.csv(DATA_PATHS$diag_1d)

tikz(filename = "diag_1d.tex", width = 6, height = 4, pointsize = 10, sanitize = TRUE)
p <- ggplot(combined_df, aes(x = factor(lag), y = x)) +
  geom_boxplot(fill = "grey95", color = "black", outlier.size = 0.8, outlier.shape = 1, outlier.color = "black") +
  ylim(0, 2) + 
  labs(title = NULL, x = "Lag $L$", y = "Estimated Log-ratio Likelihood") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.1)
  )
print(p)
dev.off()

# ##############################################################################
#### rmse ####
# ##############################################################################
# Analysis of RMSE with Bootstrap Confidence Intervals
combined_df_all <- read.csv(DATA_PATHS$rmse_data)
orc_df <- combined_df_all %>% filter(method == "orc")

rmse_boot <- function(data, indices) {
  d_val <- data[indices] 
  sqrt(mean((d_val - 1)^2))
}

summary_df <- orc_df %>%
  group_by(d, lag) %>%
  summarise(
    rmse = sqrt(mean((x - 1)^2)),
    boot_obj = list(boot(x, rmse_boot, R = 1000)),
    ci_low = quantile(boot_obj[[1]]$t, 0.025, na.rm = TRUE),
    ci_high = quantile(boot_obj[[1]]$t, 0.975, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(x_label = factor(lag, levels = c(2, 3, 4, 5, 8, 16)))

tikz("rmse.tex", width = 6, height = 3.75, sanitize = TRUE)
p <- ggplot(summary_df, aes(x = x_label, y = rmse, color = as.factor(d))) +
  geom_point(aes(shape = as.factor(d)), position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(x = "Lag $L$", y = "RMSE of relative \n normalising constant estimate", color = "Dimension", shape = "Dimension") +
  theme_bw() +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), panel.grid.minor = element_blank()) +
  scale_color_viridis_d()
print(p)
dev.off()

# ##############################################################################
#### n_d NC comparison ####
# ##############################################################################
# Comparative analysis of Normalizing Constant (NC) estimates across dimensions
df <- read.csv(DATA_PATHS$nc_comp)
d_values <- c(2, 4, 8, 16, 32, 64); lag_values <- c(2, 4, 8, 16)

plot_data <- df %>%
  filter(d %in% d_values) %>%
  filter((method == "orc" & lag %in% lag_values) | method %in% c("iapf", "bpf")) %>%
  mutate(
    plot_group = factor(case_when(
      method == "orc" ~ paste0("ORCSMC(", lag, ")"),
      method == "iapf" ~ "CSMC",
      method == "bpf" ~ "BPF",
      TRUE ~ method), levels = c(paste0("ORCSMC(", lag_values, ")"), "CSMC", "BPF")),
    d_label = factor(paste0("$d=", d, "$"), levels = paste0("$d=", d_values, "$"))
  )

tikz("comparison_nc.tex", width = 5, height = 3.75, sanitize = FALSE)
p <- ggplot(plot_data, aes(x = plot_group, y = x)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.size = 0.5) +
  facet_wrap(~ d_label, ncol = 3, scales = "fixed") +
  ylim(0, 2) +
  labs(x = 'Method', y = "Relative normalising constant estimate") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(face = "bold", size = 12))
print(p)
dev.off()

# ##############################################################################
#### svm ####
# ##############################################################################
# Stochastic Volatility Model (SVM) Log-NC Estimation
df_to_plot <- read.csv(DATA_PATHS$svm_data)
plot_data <- df_to_plot %>%
  filter(method == "ORCSMC", lag %in% c(2, 4, 8, 16)) %>%
  mutate(lag = factor(lag, levels = c(2, 4, 8, 16)))

tikz("svm.tex", width = 5, height = 3.75, sanitize = TRUE)
p <- ggplot(plot_data, aes(x = lag, y = value)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.colour = "red") +
  labs(x = "Lag $L$", y = "Log-normalising constant estimate") +
  theme_bw() +
  theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank())
print(p)
dev.off()

# ##############################################################################
#### smoothing ####
# ##############################################################################
# Visualization of Smoothing Density Quantiles
d_list <- c(2, 4, 8, 16, 32, 64)
plot_list <- list()
for (d_ in d_list) {
  # Dynamically fetch simulation output and corresponding observation data
  particle_data_name <- paste0("output_lag4_d", d_)
  output <- get(particle_data_name)
  obs_file <- paste0("orc_N1000T100_d", d_, "_obs_smoothing.csv")
  obs_ <- read.csv(obs_file, row.names = 1)
  
  # Forward/Backward Matrix Calculations
  alpha = 0.415
  tran_m <- matrix(nrow = d_, ncol = d_)
  for (i in 1:d_) { for (j in 1:d_) { tran_m[i, j] = alpha^(abs(i - j) + 1) } }
  params <- list(dt = matrix(0, d_, 1), ct = matrix(0, d_, 1), Tt = as.matrix(tran_m), 
                 P0 = diag(1, d_), Zt = diag(1, d_), Ht = diag(1, d_), Gt = diag(1, d_), a0 = rep(0, d_), d = d_)
  
  filter <- compute_fkf(params, obs_)[[2]]
  fks_means <- filter$ahatt[1, ]; fks_sds <- sqrt(filter$Vt[1, 1,])
  
  list_of_dataframes <- list()
  for (t in 1:100) {
    particle_values <- output$H_forward[[t + 1]]$X[, 1]
    list_of_dataframes[[t]] <- data.frame(value = (particle_values - fks_means[t]) / fks_sds[t], time = t)
  }
  empirical_particles_df <- bind_rows(list_of_dataframes)
  
  p <- ggplot() +
    geom_density(data = empirical_particles_df, aes(x = value, group = time, color = "Standardized Particles"), linewidth = 0.1) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), aes(color = "Standard Normal"), linewidth = 0.5) +
    scale_color_manual(values = c("Standardized Particles" = "deepskyblue", "Standard Normal" = "red")) +
    labs(title = paste("d =", d_), x = NULL, y = NULL) +
    theme_bw() + theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  plot_list[[as.character(d_)]] <- p
}

tikz("smoothing_final.tex", width = 7, height = 5, sanitize = FALSE)
y_axis_label_plot <- ggplot() + labs(y = "Density") + theme_void() + theme(axis.title.y = element_text(size = 16, angle = 90, face = "bold"))
final_plot <- y_axis_label_plot + wrap_plots(plot_list, ncol = 3) + plot_layout(widths = c(0.05, 1)) +
  plot_annotation(caption = "Standardized Value", theme = theme(plot.caption = element_text(size = 16, hjust = 0.5, face = "bold")))
print(final_plot)
dev.off()

# ##############################################################################
#### multivaraiet binomial ####
# ##############################################################################
# Multivariate Binomial Log-NC Estimation Analysis
combined_df <- read.csv(DATA_PATHS$bin_multi, row.names = 1)
plot_data <- combined_df %>%
  filter(d %in% c(2, 4, 8, 16, 32, 64)) %>%
  mutate(lag_label = factor(lag, levels = c(2, 4, 8, 16)), 
         d_label = factor(paste("d =", d), levels = paste("d =", c(2, 4, 8, 16, 32, 64))))

p <- ggplot(plot_data, aes(x = lag_label, y = x)) +
  geom_boxplot(fill = "skyblue", color = "black", outlier.size = 0.5, outlier.shape = 1, outlier.color = "black") +
  facet_wrap(~ d_label, scales = "free_y", ncol = 3) +
  labs(x = "Lag $L$", y = "Log-normalising constant estimate") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), strip.text = element_text(face = "bold", size = 12))
print(p)

# ##############################################################################
#### l1-error ####
# ##############################################################################
# Wasserstein-1 Distance (L1-Error) across time steps
df_l1 <- read.csv(DATA_PATHS$l1_error) %>%
  pivot_longer(cols = c(X1, T.2, T), names_to = "Time_Raw", values_to = "Value") %>%
  mutate(d = factor(d, levels = c(2, 4, 8, 16, 32, 64)), lag = factor(lag, levels = c(2, 4, 8, 16)), Time_Factor = factor(Time_Raw, levels = c("X1", "T.2", "T")))

tikz("l1_error.tex", width = 6, height = 3.75, sanitize = FALSE)
p <- ggplot(df_l1, aes(x = lag, y = Value, fill = Time_Factor)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1.2, outlier.colour = "red", lwd = 0.5) +
  facet_wrap(~d, nrow = 2, ncol = 3, scales = "free_y", labeller = as_labeller(c("2"="$d=2$","4"="$d=4$","8"="$d=8$","16"="$d=16$","32"="$d=32$","64"="$d=64$"))) +
  scale_fill_manual(values = c("X1" = "#FF7F50", "T.2" = "#32CD32", "T" = "#4682B4"), labels = c("X1" = "$1$", "T.2" = "$T/2$", "T" = "$T$")) +
  labs(x = "Lag $L$", y = "Wasserstein-1 distance", fill = "Time step") +
  theme_bw() + theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

# ##############################################################################
#### diagnal_nd ####
# ##############################################################################
# N-dimensional diagonal case: NC Comparative Performance
df_diag_nd <- read.csv(DATA_PATHS$diag_nd)
plot_data <- df_diag_nd %>%
  filter(d %in% c(2, 4, 8, 16, 32, 64)) %>%
  filter((method == "orc" & lag %in% c(2, 4, 8, 16)) | method %in% c("iapf", "bpf")) %>%
  mutate(plot_group = factor(case_when(method == "orc" ~ paste0("ORCSMC(", lag, ")"), method == "iapf" ~ "CSMC", method == "bpf" ~ "BPF", TRUE ~ method), levels = c("ORCSMC(2)", "ORCSMC(4)", "ORCSMC(8)", "ORCSMC(16)", "CSMC", "BPF")),
         d_label = factor(paste0("$d=", d, "$"), levels = paste0("$d=", c(2, 4, 8, 16, 32, 64), "$")))

tikz("diag_nd.tex", width = 5, height = 3.75, sanitize = FALSE)
p <- ggplot(plot_data, aes(x = plot_group, y = x)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.size = 0.5) +
  facet_wrap(~ d_label, ncol = 3, scales = "fixed") + ylim(0, 2) +
  labs(x = 'Method', y = "Relative normalising constant estimate") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(face = "bold", size = 12))
print(p)
dev.off()

# ##############################################################################
#### neuro1d ####
# ##############################################################################
# Neuro-scientific dataset: ESS Evolution and NC Distribution
ess_df_long <- read.csv(DATA_PATHS$neuro_ess) %>% select(-X) %>% pivot_longer(cols = starts_with("l"), names_to = "Lag", values_to = "ESS") %>% mutate(Lag = gsub("l", "", Lag))

p_ess <- ggplot(ess_df_long, aes(x = Time, y = ESS, color = Lag)) +
  geom_line(linewidth = 0.7) + geom_hline(yintercept = 500, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("2" = "#7a0177", "4" = "#2c7fb8", "8" = "#41b6c4", "16" = "#fed976")) +
  labs(title = "(a) Evolution of ESS.", x = "Time step", y = "ESS") + theme_bw()

logz_df <- read.csv(DATA_PATHS$neuro_nd) %>% rename(LogZ = x) %>% select(-X) %>% mutate(Lag = factor(Lag, levels = c(2, 4, 8, 16)))
p_logz <- ggplot(logz_df, aes(x = Lag, y = LogZ)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = 1, outlier.color = "red", outlier.size = 2) +
  labs(title = "(b) Log-normalising constant.", x = "Lag $L$", y = "Log-normalising constant") + theme_bw()

tikz("neuro.tex", width = 6.75, height = 3, sanitize = FALSE)
print(p_ess + p_logz + plot_layout(ncol = 2))
dev.off()

# ##############################################################################
#### neuro_nd ####
# ##############################################################################
# Multivariate Neuro-scientific data analysis
df_neuro_nd <- read.csv(DATA_PATHS$neuro_nd) %>% rename(LogZ = x) %>% select(-X) %>%
  mutate(d = factor(d, levels = c(2, 4, 8, 16, 32, 64)), lag = factor(lag, levels = c(2, 4, 8, 16)))

p <- ggplot(df_neuro_nd, aes(x = lag, y = LogZ)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.shape = 1, outlier.size = 1.2, lwd = 0.8) +
  facet_wrap(~d, nrow = 2, ncol = 3, scales = "free_y", labeller = as_labeller(c("2"="$d=2$","4"="$d=4$","8"="$d=8$","16"="$d=16$","32"="$d=32$","64"="$d=64$"))) +
  labs(x = "Lag $L$", y = "Log-normalising constant estimate") +
  theme_bw() + theme(strip.text = element_text(size = 12, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

tikz("multivariate_logz.tex", width = 6, height = 3.5, sanitize = FALSE)
print(p)
dev.off()
```

