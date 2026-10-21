# Load in libraries
library(metafor)

#############################################
# Simulation Setup
#############################################

# Setting constant prior - Needs to be in terms of I^2
constant_prior <- c(0.1, 0.2, 0.3, 0.9)

# Setting weights
weight_list <- c(1, 3, 4)

# Number of simulations
simulations <- 1000

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
  est_mean = numeric(),
  prior = numeric(),
  isquared = numeric()
)

# I-squared to tau-squared mapping based on our expected sample size
# Expected sample size is 220
# 10: 0.0020
# 20: 0.0045
# 30: 0.00779
# 40: 0.012
# 50: 0.018
# 60: 0.027
# 70: 0.042
# 80: 0.072
# 90: 0.16


# Setting based on I2 of 0%, 30%, 60%, 90%
  # I'm rounding to make them cleaner
varying_ts <- c(0, 0.008, 0.03, 0.16)
num_total_studies <- 50

# Progress bar
pb <- txtProgressBar(min = 0, max = (simulations*length(varying_ts)), style = 3)
progress <- 1

# Sample sizes
  # Must be integers and must be even
n <- sample(seq(40, 400, by = 2), num_total_studies, replace = TRUE)

# enforce even

for (sim_ts in varying_ts) {
  for (j in 1:simulations) {
    
    # Number of studies
    k <- num_total_studies
    
    # True mean of distribution of effect sizes
    theta <- 0.1
    
    # True variance of distribution of effect sizes
    tau_squared <- sim_ts
    
    # True effect sizes for each of k studies
    theta_i <- rnorm(k, mean=theta, sd=sqrt(tau_squared))
    
    # Split sample sizes into two groups
    n_c <- n / 2
    n_t <- n / 2
    
    # Participant data
    t <- c()
    v <- c()
    for (i in 1:k) {
      # Generate group observations
      control_obs <- rnorm(n_c[i], mean=0, sd=1)
      treatment_obs <- rnorm(n_t[i], mean=theta_i[i], sd=1)
      
      # Calculate g and variance
      es <- escalc(measure="SMD", m1i=mean(treatment_obs), sd1i=sd(treatment_obs), n1i=n_t[i], 
                          m2i=mean(control_obs), sd2i=sd(control_obs), n2i=n_c[i])
      
      g <- es$yi[1]
      var_g <- es$vi
      t <- append(t, g)
      v <- append(v, var_g)
    }
    
    #############################################
    # RE Throughout
    #############################################
    re_ts_list <- c()
    re_ub_list <- c()
    re_lb_list <- c()
    re_mean_list <- c()
    re_isquared_list <- c()
    mas_to_do <- c(2:20, 30, 50)
    for (i in mas_to_do) {
      re_model <- rma(yi=t[1:i], vi=v[1:i], method="REML", control = list(stepadj = 0.25, maxiter = 10000))
      re_ub <- predict(re_model)$pi.ub[[1]]
      re_lb <- predict(re_model)$pi.lb[[1]]
      re_ts <- re_model$tau2
      re_mean <- re_model$b
      re_isquared <- re_model$I2
      re_ts_list <- append(re_ts_list, re_ts)
      re_ub_list <- append(re_ub_list, re_ub)
      re_lb_list <- append(re_lb_list, re_lb)
      re_mean_list <- append(re_mean_list, re_mean)
      re_isquared_list <- append(re_isquared_list, re_isquared)
    }
    num_mas <- length(re_ts_list)
    sim <- rep(j, num_mas)
    method <- rep("REML", num_mas)
    num_studies <- c(2:20, 30, 50)
    final_k <- rep(k, num_mas)
    true_ts <- rep(tau_squared, num_mas)
    true_mean <- rep(theta, num_mas)
    est_ts <- re_ts_list
    est_ub <- re_ub_list
    est_lb <- re_lb_list
    est_mean <- re_mean_list
    prior <- rep(NA, num_mas)
    isquared <- re_isquared_list
    
    temp <- as.data.frame(cbind(sim, method, num_studies, 
                                final_k, true_ts, true_mean, est_ts, est_ub, est_lb, est_mean, prior, isquared))
    
    sim_data <- rbind(sim_data, temp)
    
    #############################################
    # Weighted RE
    #############################################
    for (p in constant_prior) {
      for (weight in weight_list) {
        w_ts_list <- c()
        w_ub_list <- c()
        w_lb_list <- c()
        w_mean_list <- c()
        w_isquared_list <- c()
        mas_to_do <- c(2:20, 30, 50)
        for (i in mas_to_do) {
          # Set weight
          if (weight == 1) {
            w <- 1/i
          } else if (weight == 2) {
            w <- 1/(i-1)
          } else if (weight == 3) {
            w <- 1/sqrt(i)
          } else {
            w <- 1/(i^2)
          }
          # Set prior
          v_bar <- mean(v[1:i])
          c <- (p/(1-p)) * v_bar
          
          # Fit the model to get an initial ts
          init_model <- rma(yi=t[1:i], vi=v[1:i], method="REML", control = list(stepadj = 0.25, maxiter = 10000))
          init_ts <- init_model$tau2
          
          # Create weighted ts
          w_ts <- (w*c) + (1-w)*init_ts
          
          # Re-run the model with this as "fixed" ts
          w_model <- rma(yi=t[1:i], vi=v[1:i], tau2=w_ts)
          w_ub <- predict(w_model)$pi.ub[[1]]
          w_lb <- predict(w_model)$pi.lb[[1]]
          w_mean <- w_model$b
          w_isquared <- w_model$I2
          w_ts_list <- append(w_ts_list, w_ts)
          w_ub_list <- append(w_ub_list, w_ub)
          w_lb_list <- append(w_lb_list, w_lb)
          w_mean_list <- append(w_mean_list, w_mean)
          w_isquared_list <- append(w_isquared_list, w_isquared)
        }
        num_mas <- length(w_ts_list)
        sim <- rep(j, num_mas)
        method <- rep(paste("Weighted_", weight, sep=""), num_mas)
        num_studies <- c(2:20, 30, 50)
        final_k <- rep(k, num_mas)
        true_ts <- rep(tau_squared, num_mas)
        true_mean <- rep(theta, num_mas)
        est_ts <- w_ts_list
        est_ub <- w_ub_list
        est_lb <- w_lb_list
        est_mean <- w_mean_list
        prior <- rep(p, num_mas)
        isquared <- w_isquared_list
        
        temp <- as.data.frame(cbind(sim, method, num_studies, 
                                    final_k, true_ts, true_mean, est_ts, est_ub, est_lb, est_mean, prior, isquared))
        
        sim_data <- rbind(sim_data, temp)
      }
    }

    #############################################
    # Update Progress bar
    #############################################
    setTxtProgressBar(pb, progress)
    progress <- progress + 1
    
  }
}

write.csv(sim_data, file = "may_5_sim_data.csv", row.names = FALSE)
