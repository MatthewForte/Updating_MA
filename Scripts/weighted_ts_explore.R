# Idea: A weighted average for tau, basically using a single number prior
# A good option to start for this would be 0.1 - psych database of ts values
# Things to tune:
  # Value of c, our constant number prior
  # Value of w, what the weights will be
  # Other weights to try: 1/k^p with p between 0.5 and 1, 1/sqrt(k-1)

# Load in libraries
library(metafor)

#############################################
# Simulation Setup
#############################################

# Setting constant prior
#constant_prior <- c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25)
constant_prior <- c(0.01, 0.02, 0.05, 0.1, 0.2)

# Setting weights
weight_list <- c(1, 2, 3, 4)

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
  prior = numeric()
)

# Setting based on I2 of 0%, 20%, 50%, 80%
varying_ts <- c(0, 0.01, 0.04, 0.16)
num_total_studies <- 50

# Progress bar
pb <- txtProgressBar(min = 0, max = (simulations*length(varying_ts)), style = 3)
progress <- 1

# Sample sizes
  # The order of the sample sizes matters and we want the order to stay constant across simulations
  # That is, if some simulations start with a big study and others a small, that doesn't make sense
  # So sample sizes are the same in each simulation, just not all equal
  # Starting with 3 small studies, then the rest of small-medium

n1 <- c(40, 40, 40)
n2 <- runif(num_total_studies - 3, min=40, max=400)
n <- c(n1, n2)

for (sim_ts in varying_ts) {
  for (j in 1:simulations) {
    
    # Number of studies
    k <- num_total_studies
    
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
    mas_to_do <- c(2:20, 30, 50)
    for (i in mas_to_do) {
      re_model <- rma(yi=t[1:i], vi=v[1:i], method="REML", control = list(stepadj = 0.25, maxiter = 10000))
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
    num_studies <- c(2:20, 30, 50)
    final_k <- rep(k, num_mas)
    true_ts <- rep(tau_squared, num_mas)
    true_mean <- rep(theta, num_mas)
    est_ts <- re_ts_list
    est_ub <- re_ub_list
    est_lb <- re_lb_list
    est_mean <- re_mean_list
    prior <- rep(NA, num_mas)
    
    temp <- as.data.frame(cbind(sim, method, num_studies, 
                                final_k, true_ts, true_mean, est_ts, est_ub, est_lb, est_mean, prior))
    
    sim_data <- rbind(sim_data, temp)
    
    #############################################
    # Weighted RE
    #############################################
    for (c in constant_prior) {
      for (weight in weight_list) {
        w_ts_list <- c()
        w_ub_list <- c()
        w_lb_list <- c()
        w_mean_list <- c()
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
      
          # Fit the model to get an initial ts
          init_model <- rma(yi=t[1:i], vi=v[1:i], method="REML", control = list(stepadj = 0.25, maxiter = 10000))
          init_ts <- init_model$tau2
      
          # Create weighted ts
          w_ts <- (w*c) + (1-w)*init_ts
      
          # Re-run the model with this as "fixed" ts
          w_model <- rma(yi=t[1:i], vi=v[1:i], method="REML", tau2=w_ts)
          w_ub <- predict(w_model)$pi.ub[[1]]
          w_lb <- predict(w_model)$pi.lb[[1]]
          w_mean <- w_model$b
          w_ts_list <- append(w_ts_list, w_ts)
          w_ub_list <- append(w_ub_list, w_ub)
          w_lb_list <- append(w_lb_list, w_lb)
          w_mean_list <- append(w_mean_list, w_mean)
        }
        num_mas <- length(w_ts_list)
        sim <- rep(j, num_mas)
        method <- rep(paste("Weighted_", weight, sep=""), num_mas)
        num_studies <- num_studies <- c(2:20, 30, 50)
        final_k <- rep(k, num_mas)
        true_ts <- rep(tau_squared, num_mas)
        true_mean <- rep(theta, num_mas)
        est_ts <- w_ts_list
        est_ub <- w_ub_list
        est_lb <- w_lb_list
        est_mean <- w_mean_list
        prior <- rep(c, num_mas)
    
        temp <- as.data.frame(cbind(sim, method, num_studies, 
                                    final_k, true_ts, true_mean, est_ts, est_ub, est_lb, est_mean, prior))
    
        sim_data <- rbind(sim_data, temp)
      }
    }
    
    #############################################
    # Max Weighted RE
    #############################################
    #for (c in constant_prior) {
    #  for (weight in weight_list) {
    #    max_w_ts_list <- c()
    #    max_w_ub_list <- c()
    #    max_w_lb_list <- c()
    #    max_w_mean_list <- c()
    #    for (i in 2:k) {
    #      # Set weight
    #      if (weight == 1) {
    #        w <- 1/i
    #      } else if (weight == 2) {
    #        w <- 1/(i-1)
    #      } else if (weight == 3) {
    #        w <- 1/sqrt(i)
    #      } else {
    #        w <- 1/(i^2)
    #      }
    #  
    #      # Fit the model to get an initial ts
    #      init_model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
    #      init_ts <- init_model$tau2
    #  
    #      # Create max weighted ts
    #      w_ts <- (w*c) + (1-w)*init_ts
    #      max_ts <- max(w_ts, init_ts)
    #  
    #      # Re-run the model with this as "fixed" ts
    #      max_w_model <- rma(yi=t[1:i], vi=v[1:i], method="REML", tau2=max_ts)
    #      max_w_ub <- predict(max_w_model)$pi.ub[[1]]
    #      max_w_lb <- predict(max_w_model)$pi.lb[[1]]
    #      max_w_mean <- max_w_model$b
    #      max_w_ts_list <- append(max_w_ts_list, max_ts)
    #      max_w_ub_list <- append(max_w_ub_list, max_w_ub)
    #      max_w_lb_list <- append(max_w_lb_list, max_w_lb)
    #      max_w_mean_list <- append(max_w_mean_list, max_w_mean)
    #    }
    #    num_mas <- length(max_w_ts_list)
    #    sim <- rep(j, num_mas)
    #    method <- rep(paste("Max_Weighted_", weight, sep=""), num_mas)
    #    num_studies <- c(2:k)
    #    final_k <- rep(k, num_mas)
    #    true_ts <- rep(tau_squared, num_mas)
    #    true_mean <- rep(theta, num_mas)
    #    est_ts <- max_w_ts_list
    #    est_ub <- max_w_ub_list
    #    est_lb <- max_w_lb_list
    #    est_mean <- max_w_mean_list
    #    prior <- rep(c, num_mas)
    #
    #    temp <- as.data.frame(cbind(sim, method, num_studies, 
    #                                final_k, true_ts, true_mean, est_ts, est_ub, est_lb, est_mean, prior))
    #
    #    sim_data <- rbind(sim_data, temp)
    #  }
    #}
    
    #############################################
    # Update Progress bar
    #############################################
    setTxtProgressBar(pb, progress)
    progress <- progress + 1
    
  }
}

write.csv(sim_data, file = "weighted_re_v2.csv", row.names = FALSE)
