
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


## Example

This is a basic example which the code works for a multivaraite linear gaussian model:

``` r
library(orc.smc)
## basic example code


Napf = N = 100
lag = 2
Time = 10
d_ = 2

alpha = 0.415
tran_m <- matrix(nrow = d_, ncol = d_)
for (i in 1:d_){
  for (j in 1:d_){
    tran_m[i,j] = alpha^(abs(i-j) + 1)
  }
}
ini <- rep(0, d_)

tran_c = diag(1, nrow = d_, ncol = d_)
ini_c = diag(1, nrow = d_, ncol = d_)
obs_m = diag(1, nrow = d_, ncol = d_)
obs_c = diag(1, nrow = d_, ncol = d_)
parameters_ <- list(k = 5, tau = 0.5, kappa = 0.5)
obs_p <- list(obs_mean = obs_m, obs_cov = obs_c)


model <- list(ini_mu = ini, ini_cov = ini_c, tran_mu = tran_m, tran_cov = tran_c, obs_params = obs_p,
  eval_likelihood = evaluate_likelihood_lg, simu_observation = simulate_observation_lg,
  parameters = parameters_, dist = 'lg')

set.seed(1234)
obs_ <- sample_obs(model, Time, d_) #provided by users

a0_ = ini      # Initial state mean
P0_ = ini_c    # Initial state covariance
Zt_ = obs_m    # Observation matrix (C)
Ht_ = tran_c   # Observation noise covariance (R)
Gt_ = obs_c  # Process noise covariance (Q)

dt_ <- ct_ <- matrix(0, d_, 1)
Tt_ <- as.matrix(tran_m)
params <- list(dt = dt_, ct = ct_, Tt = Tt_, P0 = P0_ , Zt = Zt_,
                Ht = Ht_, Gt = Gt_, a0 = a0_, d = d_)

filter <- compute_fkf(params, obs_)
fkf_obj <- filter[[1]]

data_ = data <- list(obs = obs_)

output <- Orc_SMC(lag, data, model, Napf)

log_ratio <- compute_ratio(output$logZ[Time], fkf_obj)
```
