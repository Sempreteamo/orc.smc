
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

Napf = N = 1000 # Number of particles (N)
lag = 4 # Rolling window length (L)
Time = 100 # Total time steps (T)
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

## Example 2: Fully Adaptive ORC-SMC with Benchmarking

This section extends the basic model to showcase the **Fully Adaptive ORCSMC** algorithm. It dynamically adjusts the rolling window ($B_t$), the number of learning steps ($K_t$), and the particle swarm size ($N_t$) based on the computational budget $\gamma$. 

``` r
## (Assuming the Model Specification and FKF setup from Example 1 have been run)

## --- 1. Adaptive Algorithm Hyperparameters ---
B_max <- 5        # Maximum rolling window length
K_min <- 1        # Minimum learning iterations
K_max <- 5        # Maximum learning iterations
K1    <- 2        # Initial K at t=1
gamma <- 10000    # Computational budget parameter
N1    <- floor(gamma / (1 * (K1 + 1))) # Initial number of particles

output <- adaptative_Orc_SMC(
      B_max = B_max, 
      K_min = K_min, 
      K_max = K_max, 
      K1    = K1, 
      gamma = gamma, 
      data  = data_, 
      model = model, 
      N1    = N1
    )

ratio_vec_orc <- compute_ratio(output$logZ[Time], fkf_obj)
```

## Experiment for Figure 1
The experiment is repeated 50 times
``` r
d_values   <- c(2, 4, 8, 16, 32, 64)
lag_values <- c("2", "4", "8", "16") 
Napf       <- 1000  
N_bpf      <- 320000 
Time       <- 100
alpha      <- 0.415
K_iterations <- 5

results_list <- list()


for (d in d_values) {
  
  model <- list(
    ini_mu = rep(0, d), ini_cov = diag(1, d),
    tran_mu = diag(1, d), tran_cov = diag(1, d),
    obs_params = list(obs_mean = diag(1, d), obs_cov = diag(1, d)),
    eval_likelihood = evaluate_likelihood_lg,
    simu_observation = simulate_observation_lg,
    parameters = list(k = 5, tau = 0.5, kappa = 0.5)
  )
  
  set.seed(1234)
  obs_   <- sample_obs(model, Time, d)
  data_  <- list(obs = obs_)
  
  
  params_fkf <- list(dt=matrix(0,d,1), ct=matrix(0,d,1), Tt=as.matrix(model$tran_mu),
                     P0=diag(1,d), Zt=diag(1,d), Ht=diag(1,d), Gt=diag(1,d), a0=rep(0,d), d=d)
  fkf_logZ <- compute_fkf(params_fkf, obs_)[[1]]
  
  #run orcsmc#
  for (l_char in lag_values) {
    lag_val <- as.numeric(l_char)
    output  <- Orc_SMC(lag_val, data_, model, Napf)
    x_val   <- compute_ratio(output$logZ[Time], fkf_logZ)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      X = NA, x = x_val, d = d, lag = l_char, method = "orc"
    )
  }
  
  #run bpf#
  output_bpf <- run_bpf(data = data_, model, N = N_bpf)
  x_val_bpf  <- compute_ratio(output_bpf$logZ, fkf_logZ)
  
  results_list[[length(results_list) + 1]] <- data.frame(
    X = NA, x = x_val_bpf, d = d, lag = "none", method = "bpf"
  )
  
  #run csmc#
  output_iapf <- run_CSMC(data = data_, Napf = Napf, K = K_iterations, model = model)
  
  
  x_val_iapf <- compute_ratio(output_iapf$logZ_final, fkf_logZ)
  
  
  results_list[[length(results_list) + 1]] <- data.frame(
    X = NA, x = x_val_iapf, d = d, lag = "none", method = "csmc"
  )
}


final_df <- bind_rows(results_list)

final_df$X <- 1:nrow(final_df)

head(final_df)

write.csv(final_df, "diag0.415_nc_orc+bpf+iapf_N1000T100_d2-64_l2-16.csv", row.names = FALSE)


```

## Experiment for Figure 2 and 3
The experiment is repeated 50 times
``` r
d_values   <- c(2, 4, 8, 16, 32, 64)
lag_values <- c("2", "4", "8", "16") 
Napf       <- 1000  
N_bpf      <- 320000 
Time       <- 100
alpha      <- 0.415
K_iterations <- 5

results_list <- list()


for (d in d_values) {
  
  
  tran_m <- matrix(0, d, d)
  for (i in 1:d) for (j in 1:d) tran_m[i, j] <- alpha^(abs(i - j) + 1)
  
  model <- list(
    ini_mu = rep(0, d), ini_cov = diag(1, d),
    tran_mu = tran_m, tran_cov = diag(1, d),
    obs_params = list(obs_mean = diag(1, d), obs_cov = diag(1, d)),
    eval_likelihood = evaluate_likelihood_lg,
    simu_observation = simulate_observation_lg,
    parameters = list(k = 5, tau = 0.5, kappa = 0.5)
  )

  set.seed(1234)
  obs_   <- sample_obs(model, Time, d)
  data_  <- list(obs = obs_)
  
  
  params_fkf <- list(dt=matrix(0,d,1), ct=matrix(0,d,1), Tt=as.matrix(tran_m),
                     P0=diag(1,d), Zt=diag(1,d), Ht=diag(1,d), Gt=diag(1,d), a0=rep(0,d), d=d)
  fkf_logZ <- compute_fkf(params_fkf, obs_)[[1]]
  

  for (l_char in lag_values) {
    lag_val <- as.numeric(l_char)
    output  <- Orc_SMC(lag_val, data_, model, Napf)
    x_val   <- compute_ratio(output$logZ[Time], fkf_logZ)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      X = NA, x = x_val, d = d, lag = l_char, method = "orc"
    )
  }
  

  output_bpf <- run_bpf(data = data_, model, N = N_bpf)
  x_val_bpf  <- compute_ratio(output_bpf$logZ, fkf_logZ)
  
  results_list[[length(results_list) + 1]] <- data.frame(
    X = NA, x = x_val_bpf, d = d, lag = "none", method = "bpf"
  )
  
  
  output_iapf <- run_CSMC(data = data_, Napf = Napf, K = K_iterations, model = model)
  
  
  x_val_iapf <- compute_ratio(output_iapf$logZ_final, fkf_logZ)
  
  
  results_list[[length(results_list) + 1]] <- data.frame(
    X = NA, x = x_val_iapf, d = d, lag = "none", method = "csmc"
  )
}


final_df <- bind_rows(results_list)

final_df$X <- 1:nrow(final_df)

head(final_df)

write.csv(final_df, "orc+bpf+iapf_N1000T100_d2-64_lag2-16_non-diagf_rep100.csv", row.names = FALSE)


```

## Experiment for Figure 4

``` r
library(dplyr)
library(tidyr)

d_values     <- c(2, 4, 8, 16, 32, 64)
lag_values   <- c("2", "4", "8", "16") 
Napf         <- 1000  
Time         <- 100
alpha        <- 0.415
n_repeats    <- 50 

all_results <- list()

for (d in d_values) {
  
  tran_m <- matrix(0, d, d)
  for (i in 1:d) for (j in 1:d) tran_m[i, j] <- alpha^(abs(i - j) + 1)
  
  model <- list(
    ini_mu = rep(0, d), ini_cov = diag(1, d),
    tran_mu = tran_m, tran_cov = diag(1, d),
    obs_params = list(obs_mean = diag(1, d), obs_cov = diag(1, d)),
    eval_likelihood = evaluate_likelihood_lg,
    simu_observation = simulate_observation_lg,
    parameters = list(k = 5, tau = 0.5, kappa = 0.5)
  )
  
  set.seed(1234)
  obs_ <- sample_obs(model, Time, d)
  params_fkf <- list(dt=matrix(0,d,1), ct=matrix(0,d,1), Tt=as.matrix(tran_m),
                     P0=diag(1,d), Zt=diag(1,d), Ht=diag(1,d), Gt=diag(1,d), 
                     a0=rep(0,d), d=d)
  filter_res <- compute_fkf(params_fkf, obs_)[[2]]
  
  for (l_char in lag_values) {
    lag_val <- as.numeric(l_char)
    
    for (r in 1:n_repeats) {
      output <- Orc_SMC(lag_val, list(obs = obs_), model, Napf)
      
      # Calculate L1 error at t=1, t=T/2, t=T
      test_times <- c(1, floor(Time/2), Time)
      
      for (t in test_times) {
        p_vals <- output$H_forward[[t + 1]]$X[, 1] # First coordinate
        m_t <- filter_res$ahatt[1, t]
        s_t <- sqrt(filter_res$Vt[1, 1, t])
        
        calc_w1 <- function(particles, true_mean, true_sd) {
          n <- length(particles)
          sorted_p <- sort(particles)
          # Quantiles of the target Gaussian distribution
          target_q <- qnorm(seq(1/(n+1), n/(n+1), length.out = n), mean = true_mean, sd = true_sd)
          return(mean(abs(sorted_p - target_q)))
        }
        
        err <- calc_w1(p_vals, m_t, s_t)
        
        all_results <- rbind(all_results, data.frame(
          d = d,
          lag = lag_val,
          rep = r,
          Time_Point = ifelse(t==1, "1", ifelse(t==Time, "T", "T/2")),
          Value = err
        ))
      }
    }
  }
}

long_df <- bind_rows(all_results)

fig <- long_df %>%
  group_by(d, lag, rep, Time_Point) %>%
  summarise(Value = mean(Value), .groups = 'drop') %>%
  pivot_wider(names_from = Time_Point, values_from = Value) %>%
  rename(X1 = `1`, T.2 = `T/2`, T = `T`)

fig <- fig %>% select(X1, T.2, T, d, lag)

head(fig)

write.csv(fig,"l1error_orc_N1000T100_d2-64_lag2_16_rep100.csv", row.names = FALSE)
```

## Experiment for Figure 5
The experiment is repeated 50 times
``` r
alpha_svm  <- 0.986
sigma_svm  <- 0.13
beta_svm   <- 0.69
Time_svm   <- 945
N_particles <- 200 
n_repeats   <- 50  

ini_var_svm <- sigma_svm^2 / (1 - alpha_svm^2)

model_svm <- list(
  ini_mu  = 0, 
  ini_cov = as.matrix(ini_var_svm), 
  tran_mu = as.matrix(alpha_svm), 
  tran_cov = as.matrix(sigma_svm^2), 
  obs_params = list(beta = beta_svm, den_cov = beta_svm^2), 
  eval_likelihood = evaluate_likelihood_svm,
  simu_observation = simulate_observation_svm, 
  parameters = list(k = 5, tau = 0.5, kappa = 0.5) 
)

set.seed(123)
obs_svm <- sample_obs(model_svm, Time_svm, d = 1)
data_svm <- list(obs = obs_svm)

lag_values <- c(2, 4, 8, 16)
results_list <- list()

for (l_val in lag_values) {
  for (r in 1:n_repeats) {
    
    output <- Orc_SMC(l_val, data_svm, model_svm, N_particles)
    
    val <- output$logZ[Time_svm]
    
    
    results_list[[length(results_list) + 1]] <- data.frame(
      value = val,
      method = "ORCSMC",
      lag = l_val
    )
  }
}


fig <- do.call(rbind, results_list)

head(fig)
write.csv(fig,'bpf+orc_N200_d1_lag2-16_svm_rep100.csv', row.names = FALSE)

```

## Experiment for Figure 6

``` r
library(dplyr)
library(tidyr)

alpha_neuro  <- 0.99
sigma2_neuro <- 0.11
M_neurons    <- 50
Time         <- 100
Napf         <- 1000
lag_list     <- c(2, 4, 8, 16) 
d_values     <- c(2, 4, 8, 16, 32, 64)
n_repeats    <- 50  

#d = 1
model_1d <- list(
  ini_mu           = 0, 
  ini_cov          = as.matrix(1.0),
  tran_mu          = as.matrix(alpha_neuro), 
  tran_cov         = as.matrix(sigma2_neuro), 
  obs_params       = M_neurons,
  eval_likelihood  = evaluate_likelihood_bin,
  simu_observation = simulate_observation_bin,
  parameters       = list(k = 5, tau = 0.5, kappa = 0.5)
)


set.seed(1234)
obs_data <- sample_obs(model_1d, Time, d = 1)
ess_data<- list(
  X = 1:Time, 
  Time = 1:Time
)
logz_df <- data.frame()

for (L in lag_list) {
  #cat("Running Lag =", L, "\n")
  for (i in 1:n_repeats) {
    # Suppress internal output for cleaner execution
    
    res <- Orc_SMC(L, list(obs = obs_data), model_1d, Napf)
    
    
    # Check if ESS was returned correctly to prevent the "0 row" error
    if (!is.null(res$ess_history) && length(res$ess_history) == Time) {
      if (i == 1) { # Only store one trajectory per Lag for the ESS plot
        col_name <- paste0("l", L)
        ess_data[[col_name]] <- res$ess_history
      }
    }
    
    # Store LogZ for all repeats (for the boxplot)
    logz_df <- rbind(logz_df, data.frame(LogZ = res$logZ[Time], Lag = factor(L, levels = lag_list)))
  }
}

logz_df$X <- 1:nrow(logz_df)
neuro_1d <- logz_df %>%
  
  mutate(
    X = 1:n(),
    x = LogZ
  ) %>%
  
  select(X, x, Lag)

ess <- as.data.frame(ess_data)

write.csv(ess, "bin_ess_1d_l2-16.csv", row.names = FALSE)
write.csv(neuro_1d, "bin_N1000T100_d1_lag2-16_rep100.csv", row.names = FALSE)
```

## Experiment for Figure 7

``` r
library(dplyr)
library(tidyr)

alpha_neuro  <- 0.99
sigma2_neuro <- 0.11
M_neurons    <- 50
Time         <- 100
Napf         <- 1000
lag_list     <- c(2, 4, 8, 16) 
d_values     <- c(2, 4, 8, 16, 32, 64)
n_repeats    <- 50  

multivariate_results <- list()

for (d in d_values) {
  
  model_nd <- list(
    ini_mu           = rep(0, d), 
    ini_cov          = diag(1, d),
    tran_mu          = diag(alpha_neuro, d), 
    tran_cov         = diag(sigma2_neuro, d), 
    obs_params       = M_neurons,
    eval_likelihood  = evaluate_likelihood_bin,
    simu_observation = simulate_observation_bin,
    parameters       = list(k = 5, tau = 0.5, kappa = 0.5)
  )
  
  obs_nd  <- sample_obs(model_nd, Time, d)
  data_nd <- list(obs = obs_nd)
  
  for (l_val in lag_list) {
    for (r in 1:n_repeats) {
      output_nd <- Orc_SMC(l_val, data_nd, model_nd, Napf)
      
      multivariate_results[[length(multivariate_results) + 1]] <- data.frame(
        X = NA, x = output_nd$logZ[Time], d = d, lag = as.character(l_val), method = "orc"
      )
    }
  }
}

neuro_nd <- do.call(rbind, multivariate_results)

neuro_nd$X <- 1:nrow(neuro_nd) 

neuro_nd$d   <- as.numeric(as.character(neuro_nd$d))
neuro_nd$lag <- as.numeric(as.character(neuro_nd$lag))

head(neuro_nd)

write.csv(neuro_nd, "bin_N1000T100_d2-64_lag2-16_rep100.csv", row.names = FALSE)
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
  neuro_1d  = "bin_N1000T100_d1_lag2-16_rep100.csv",
  neuro_ess = "bin_ess_1d_l2-16.csv",
  neuro_nd  = "bin_N1000T100_d2-64_lag2-16_rep100.csv",
  svm_data  = "bpf+orc_N200_d1_lag2-16_svm_rep100.csv",
  l1_error  = "l1error_orc_N1000T100_d2-64_lag2_16_rep100.csv",
  diag_nd   = "diag0.415_nc_orc+bpf+iapf_N1000T100_d2-64_l2-16.csv",
  nc_comp   = rmse_data = "orc+bpf+iapf_N1000T100_d2-64_lag2-16_non-diagf_rep100.csv",
  df_diag_nd = "adaptive_orc_bpf_10000T100d2-64.csv"
)

# ##############################################################################
#### Figure 1 ####
# ##############################################################################
# Comparative analysis of Normalizing Constant (NC) estimates across dimensions for diagnal case.
df_diag_nd <- read.csv(DATA_PATHS$diag_nd)
plot_data <- df_diag_nd %>%
  filter(d %in% c(2, 4, 8, 16, 32, 64)) %>%
  filter((method == "orc" & lag %in% c(2, 4, 8, 16)) | method %in% c("csmc", "bpf")) %>%
  mutate(plot_group = factor(case_when(method == "orc" ~ paste0("ORCSMC(", lag, ")"), method == "csmc" ~ "CSMC", method == "bpf" ~ "BPF", TRUE ~ method), levels = c("ORCSMC(2)", "ORCSMC(4)", "ORCSMC(8)", "ORCSMC(16)", "CSMC", "BPF")),
         d_label = factor(paste0("$d=", d, "$"), levels = paste0("$d=", c(2, 4, 8, 16, 32, 64), "$")))

plot_data$x <- as.numeric(as.character(plot_data$x))
tikz("diag_nd.tex", width = 5, height = 3.75, sanitize = FALSE)
p <- ggplot(plot_data, aes(x = plot_group, y = x)) +
  
  geom_boxplot(color = "black", outlier.size = 0.5) +
  
  facet_wrap(~ d_label, ncol = 3, scales = "fixed") +
  ylim(0, 2) +
  labs(x = 'Method', y = "Relative normalising constant estimate") +
  #
  theme_bw() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p)
dev.off()

# ##############################################################################
#### Figure 2 ####
# ##############################################################################
# Comparative analysis of Normalizing Constant (NC) estimates across dimensions for non-diagnal case.
df <- read.csv(DATA_PATHS$nc_comp)
d_values <- c(2, 4, 8, 16, 32, 64); lag_values <- c(2, 4, 8, 16)

plot_data <- df %>%
  filter(d %in% d_values) %>%
  filter((method == "orc" & lag %in% lag_values) | method %in% c("csmc", "bpf")) %>%
  mutate(
    plot_group = factor(case_when(
      method == "orc" ~ paste0("ORCSMC(", lag, ")"),
      method == "csmc" ~ "CSMC",
      method == "bpf" ~ "BPF",
      TRUE ~ method), levels = c(paste0("ORCSMC(", lag_values, ")"), "CSMC", "BPF")),
    d_label = factor(paste0("$d=", d, "$"), levels = paste0("$d=", d_values, "$"))
  )

tikz("comparison_nc.tex", width = 5, height = 3.75, sanitize = FALSE)
p <- ggplot(plot_data, aes(x = plot_group, y = x)) +
  
  geom_boxplot(fill = "white", color = "black", outlier.size = 0.5) +
  facet_wrap(~ d_label, ncol = 3, scales = "fixed") +
  ylim(0, 2) +
  labs(x = 'Method', y = "Relative normalising constant estimate") +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p)
dev.off()

# ##############################################################################
#### Figure 3 ####
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

tikz("rmse.tex", width = 5, height = 2.75, sanitize = TRUE)
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
#### Figure 4 ####
# ##############################################################################
# Wasserstein-1 Distance (L1-Error) across time steps
df_l1 <- read.csv(DATA_PATHS$l1_error) %>%
  pivot_longer(cols = c(X1, T.2, T), names_to = "Time_Raw", values_to = "Value") %>%
  mutate(
    d = factor(d, levels = c(2, 4, 8, 16, 32, 64)), 
    lag = factor(lag, levels = c(2, 4, 8, 16)), 
    Time_Factor = factor(Time_Raw, levels = c("X1", "T.2", "T"))
  )

tikz("l1_error.tex", width = 6, height = 3.75, sanitize = FALSE)

p <- ggplot(df_l1, aes(x = lag, y = Value, fill = Time_Factor)) +
  geom_boxplot(
    outlier.shape = 19, 
    outlier.size = 0.5
  ) +
  scale_fill_grey(start = 0.7, end = 0.9, labels = c("X1" = "$1$", "T.2" = "$T/2$", "T" = "$T$")) +
  facet_wrap(~d, nrow = 2, ncol = 3, scales = "free_y", 
             labeller = as_labeller(c("2"="$d=2$","4"="$d=4$","8"="$d=8$","16"="$d=16$","32"="$d=32$","64"="$d=64$"))) +
  labs(x = "Lag $L$", y = "Wasserstein-1 distance", fill = "Time step")

print(p)
dev.off()

# ##############################################################################
#### Figure 5 ####
# ##############################################################################
# Stochastic Volatility Model (SVM) Log-NC Estimation
df_to_plot <- read.csv(DATA_PATHS$svm_data)
plot_data <- df_to_plot %>%
  filter(method == "ORCSMC", lag %in% c(2, 4, 8, 16)) %>%
  mutate(lag = factor(lag, levels = c(2, 4, 8, 16)))

plot_data$lag <- factor(plot_data$lag, levels = c(2, 4, 8, 16))

theme_set(theme_bw())

tikz("svm.tex", width = 4, height = 2.5, sanitize = FALSE)
p <- ggplot(plot_data, aes(x = lag, y = value)) +
  geom_boxplot(outlier.size = 0.5) + 
  labs(x = "Lag $L$", y = "Log-normalising constant \n estimate \n") 
print(p)
dev.off()

# ##############################################################################
#### Figure 6 ####
# ##############################################################################
# Neuro-scientific dataset: ESS Evolution and NC Distribution
theme_set(theme_bw())

ess_df_long <- read.csv(DATA_PATHS$neuro_ess) %>% 
  select(-X) %>% 
  pivot_longer(cols = starts_with("l"), names_to = "Lag", values_to = "ESS") %>% 
  mutate(Lag = factor(as.numeric(gsub("l", "", Lag)), levels = c(2, 4, 8, 16)))

logz_df <- read.csv(DATA_PATHS$neuro_1d) %>% 
  rename(LogZ = x) %>%          
  select(-X) %>%                
  mutate(Lag = factor(Lag, levels = c(2, 4, 8, 16)))

p_ess <- ggplot(ess_df_long, aes(x = Time, y = ESS, color = Lag)) +
  geom_line(linewidth = 0.7) + 
  geom_hline(yintercept = 500, linetype = "dashed") +
  scale_color_viridis_d(option = "viridis") + 
  labs(x = "Time step", y = "ESS")

p_logz <- ggplot(logz_df, aes(x = Lag, y = LogZ)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 0.5) + 
  labs(
    x = "Lag $L$", 
    y = "Log-normalising constant \n estimate \n"
  )


tikz("neuro_ess.tex", width = 3.3, height = 3, sanitize = FALSE)
print(p_ess)
dev.off()

tikz("neuro_logz.tex", width = 3.3, height = 3, sanitize = FALSE)
print(p_logz)
dev.off()


# ##############################################################################
#### Figure 7 ####
# ##############################################################################
# Multivariate Neuro-scientific data analysis
df_neuro_nd <- read.csv(DATA_PATHS$neuro_nd) %>% rename(LogZ = x) %>% select(-X) %>%
  mutate(d = factor(d, levels = c(2, 4, 8, 16, 32, 64)), lag = factor(lag, levels = c(2, 4, 8, 16)))

theme_set(theme_bw())

p <- ggplot(df_neuro_nd, aes(x = lag, y = LogZ)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 0.5) + 
  facet_wrap(~d, nrow = 2, ncol = 3, scales = "free_y", 
             labeller = as_labeller(c("2"="$d=2$","4"="$d=4$","8"="$d=8$","16"="$d=16$","32"="$d=32$","64"="$d=64$"))) +
  labs(x = "Lag $L$", y = "Log-normalising constant estimate") +
  theme_bw()

tikz("multivariate_logz.tex", width = 6, height = 3.5, sanitize = FALSE)
print(p)
dev.off()

# ##############################################################################
#### Figure 8 ####
# ##############################################################################

df_all_history <- bind_rows(history_list)

df_summary <- df_all_history %>%
  group_by(Time) %>%
  summarise(
    mean_K = mean(K, na.rm = TRUE),
    min_K  = min(K, na.rm = TRUE),
    max_K  = max(K, na.rm = TRUE),
    mean_B = mean(B, na.rm = TRUE),
    min_B  = min(B, na.rm = TRUE),
    max_B  = max(B, na.rm = TRUE),
    mean_N = mean(N, na.rm = TRUE),
    min_N  = min(N, na.rm = TRUE),
    max_N  = max(N, na.rm = TRUE)
  )

theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

p_K <- ggplot(df_summary, aes(x = Time)) +
  geom_ribbon(aes(ymin = min_K, ymax = max_K), fill = "gray90") +
  geom_line(aes(y = mean_K), color = "black", linewidth = 0.8) +
  labs(x = NULL, y = "$K_t$") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_B <- ggplot(df_summary, aes(x = Time)) +
  geom_ribbon(aes(ymin = min_B, ymax = max_B), fill = "gray90") +
  geom_line(aes(y = mean_B), color = "black", linewidth = 0.8) +
  labs(x = NULL, y = "$B_t$") + # 仅保留符号
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_N <- ggplot(df_summary, aes(x = Time)) +
  geom_ribbon(aes(ymin = min_N, ymax = max_N), fill = "gray90") +
  geom_line(aes(y = mean_N), color = "black", linewidth = 0.8) +
  labs(x = "Time step $t$ when $d = 2$", y = "$N_t$") # 仅保留符号

tikz("adaptive_convergence2.tex", width = 5.5, height = 5.5, sanitize = FALSE)

grid.arrange(p_K, p_B, p_N, ncol = 1, heights = c(1, 1, 1.1))

dev.off()

# ##############################################################################
#### Figure 9 ####
# ##############################################################################
df_diag_nd <- read.csv(DATA_PATHS$df_diag_nd)
plot_data <- df_diag_nd %>%
  filter(d %in% c(2, 4, 8, 16, 32, 64)) %>%
  filter(method %in% c("orc", "bpf")) %>%  
  mutate(
    plot_group = case_when(
      method == "orc" ~ "\\begin{tabular}{c}Adaptive\\\\ORCSMC\\end{tabular}",
      method == "bpf" ~ "BPF",
      TRUE ~ method
    ),
    plot_group = factor(plot_group, levels = c("\\begin{tabular}{c}Adaptive\\\\ORCSMC\\end{tabular}", "BPF")),
    
   
    d_label = factor(
      paste0("$d=", d, "$"),
      levels = paste0("$d=", c(2,4,8,16,32,64), "$")
    )
  )


plot_data$x <- as.numeric(plot_data$x)


tikz("diag_nd_only_orc_bpf.tex", width = 5, height = 3.75, sanitize = FALSE)

p <- ggplot(plot_data, aes(x = plot_group, y = x)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ d_label, ncol = 3, scales = "fixed") +
  coord_cartesian(ylim = c(0, 2)) + 
  labs(
    x = 'Method',
    y = "Relative normalising constant estimate"
  ) +
  theme_bw() +
  
  theme(
    axis.text.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))
  )

print(p)
dev.off()

