
# orc.smc

<!-- badges: start -->
<!-- badges: end -->

The goal of orc.smc is to provide **tools for implementing online rolling controlled sequential Monte Carlo (SMC) methods**. 
This package is designed to facilitate robust normalizing constant estimate in dynamic systems with a rolling-window setting, particularly where real-time adaptation and efficient handling of streaming data are crucial. It aims to offer a flexible and extensible framework for developing and applying advanced particle filtering techniques in evolving environments.


## Installation

You can install the development version of orc.smc like so:

```r
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
log_ratio <- compute_ratio(output$logZ[Time], fkf_obj)
```
