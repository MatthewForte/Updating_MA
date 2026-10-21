library(metafor)
library(brms)
library(bayesmeta)
library(invgamma)
library(cmdstanr)
library(posterior)

meta_mod <- cmdstan_model("bayes_meta_trunc_tau.stan")

fit_meta_trunc_tau <- function(y, sigma,
                               chains = 4,
                               iter_warmup = 1000,
                               iter_sampling = 1000,
                               seed = 1234) {
  stopifnot(length(y) == length(sigma))
  
  K <- length(y)
  
  data_list <- list(
    K = K,
    y = as.vector(y),
    sigma = as.vector(sigma)
  )
  
  fit <- meta_mod$sample(
    data = data_list,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = 0
  )
  
  # Get draws as a data frame
  draws_df <- as_draws_df(fit$draws(c("mu", "tau", "tau2", "theta_new")))
  
  # Posterior means
  mu_mean   <- mean(draws_df$mu)
  tau2_mean <- mean(draws_df$tau2)
  
  # 95% posterior predictive interval for a new study's true effect
  pred_int  <- quantile(draws_df$theta_new, probs = c(0.025, 0.975), names = FALSE)
  
  list(
    fit       = fit,              # full cmdstanr fit object if you want to dig deeper
    mu_mean   = mu_mean,
    tau2_mean = tau2_mean,
    pred_lower = pred_int[1],
    pred_upper = pred_int[2]
  )
}

# What I'm not varying
# sample size/vi -> all set to 0.04
# True mean effect -> set to 0.1
#############################################
# Simulation Setup
#############################################
# Progress bar and simulations
simulations <- 10

# Setting random seed
set.seed(30)

# Set up data frame
sim_data <- data.frame(
  sim = numeric(),
  method = character(),
  num_studies = numeric(),
  final_k = numeric(),
  true_ts = numeric(),
  true_mean = numeric(),
  est_ts = numeric(),
  est_ub = numeric(), 
  est_lb = numeric(), 
  est_mean = numeric()
)

# Setting based on I2 of 20%, 50%, 80%
varying_ts <- c(0.01, 0.04, 0.16)
num_total_studies <- 10

pb <- txtProgressBar(min = 0, max = (simulations*length(varying_ts)), style = 3)
progress <- 1

for (sim_ts in varying_ts) {
  for (j in 1:simulations) {
    
    # Number of studies
    k <- num_total_studies
    
    # Sample sizes (all equal for now)
    n <- rep(100, k)
    
    # Within study variances
    v <- 4 / round(n)
    
    # True mean of distribution of effect sizes
    theta <- 0.1
    
    # True variance of distribution of effect sizes
    tau_squared <- sim_ts
    
    # True effect sizes for each of k studies
    theta_i <- rnorm(k, mean=theta, sd=sqrt(tau_squared))
    
    # Observed effect sizes for each of k studies
    t <- c()
    for (i in 1:k) {
      ti <- rnorm(1, mean=theta_i[i], sd=sqrt(v[i]))
      t <- append(t, ti)
    }
    
    #############################################
    # RE Throughout
    #############################################
    re_ts_list <- c()
    re_ub_list <- c()
    re_lb_list <- c()
    re_mean_list <- c()
    for (i in 2:k) {
      re_model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
      re_ub <- predict(re_model)$pi.ub[[1]]
      re_lb <- predict(re_model)$pi.lb[[1]]
      re_ts <- re_model$tau2
      re_mean <- re_model$b
      re_ts_list <- append(re_ts_list, re_ts)
      re_ub_list <- append(re_ub_list, re_ub)
      re_lb_list <- append(re_lb_list, re_lb)
      re_mean_list <- append(re_mean_list, re_mean)
    }
    
    num_mas <- length(re_ts_list)
    sim <- rep(j, num_mas)
    method <- rep("RE", num_mas)
    num_studies <- c(2:k)
    final_k <- rep(k, num_mas)
    true_ts <- rep(tau_squared, num_mas)
    true_mean <- rep(theta, num_mas)
    est_ts <- re_ts_list
    est_ub <- re_ub_list
    est_lb <- re_lb_list
    est_mean <- re_mean_list
    
    temp <- as.data.frame(cbind(sim, method, num_studies, 
                                final_k, true_ts, true_mean, est_ts, est_ub, est_lb, est_mean))
    
    sim_data <- rbind(sim_data, temp)
    
    #############################################
    # Bayesian - Maybe can't do a constant prior
      # Trying just a very tight prior
    #############################################
    bayes_ts_list <- c()
    bayes_ub_list <- c()
    bayes_lb_list <- c()
    bayes_mean_list <- c()

    # Running MA's
    for (i in 2:k) {
      bayes_model <- fit_meta_trunc_tau(t[1:i], sqrt(v[1:i]))
      bayes_mean <- bayes_model$mu_mean
      bayes_ts <- bayes_model$tau2_mean
      bayes_lb <- bayes_model$pred_lower
      bayes_ub <- bayes_model$pred_upper
      bayes_ts_list <- append(bayes_ts_list, bayes_ts)
      bayes_ub_list <- append(bayes_ub_list, bayes_ub)
      bayes_lb_list <- append(bayes_lb_list, bayes_lb)
      bayes_mean_list <- append(bayes_mean_list, bayes_mean)
    }
    
    # Storing stuff
    num_mas <- length(bayes_ts_list)
    sim <- rep(j, num_mas)
    method <- rep("Bayesian", num_mas)
    num_studies <- c(2:k)
    final_k <- rep(k, num_mas)
    true_ts <- rep(tau_squared, num_mas)
    true_mean <- rep(theta, num_mas)
    est_ts <- bayes_ts_list
    est_ub <- bayes_ub_list
    est_lb <- bayes_lb_list
    est_mean <- bayes_mean_list
    
    temp <- as.data.frame(cbind(sim, method, num_studies, 
                                final_k, true_ts, true_mean, est_ts, est_ub, est_lb, est_mean))
    
    sim_data <- rbind(sim_data, temp)
    
    #############################################
    # Update Progress bar
    #############################################
    setTxtProgressBar(pb, progress)
    progress <- progress + 1
    
  }
}
close(pb)
write.csv(sim_data, file = "sim_data_bayesian_single_prior.csv", row.names = FALSE)