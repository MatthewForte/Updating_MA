library(metafor)
# LINEX 
# I think what is happening - LINEX is overestimating variance, but consistently throughout
  # That is, it overestimates at 2 studies, but also 10 studies
  # So it's not really forcing monotonicity
  # Need the inflation factor to shrink as the number of studies increase
# This simulation will look at only the LINEX part of it, and try a bunch of different inflation factors
  # Original sim set a to -2

# Setting c to 2? maybe try 3 after if things look weird
# c / k
# c / k-1
# c / sqrt(k)


######################################################################
# Simulation
######################################################################

# Progress bar and simulations
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
  est_mean = numeric()
)

# Setting based on I2 of 20%, 50%, 80%
varying_ts <- c(0.01, 0.04, 0.16)
num_total_studies <- 10

# Tuning parameter
alphas <- c(-3, -2.5, -2, -1.5, -1, 0.1, 0.01, 0.001)
c <- 2

# Progress bar
pb <- txtProgressBar(min = 0, max = (simulations*length(varying_ts)*length(alphas)), style = 3)
progress <- 1

for (alpha in alphas) {
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
      # LINEX Loss
      #############################################
      # Have to use a parametric bootstrap to calculate this
    
      # First
      linex_ts_list <- c()
      linex_ub_list <- c()
      linex_lb_list <- c()
      linex_mean_list <- c()
      for (i in 2:k) {
        # First, fit a RE model and get the mean
        init_model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
        re_theta <- init_model$beta[1]
        re_ts <- init_model$tau2
      
        # Then, get the estimated study variances (vi + ts)
        variances <- v[1:i] + re_ts
      
        # Now, bootstrap
        boot_ts_list <- c()
        for (b in 1:1000) {
          boot_Ti <- rnorm(i, mean=re_theta, sd=sqrt(variances))
          boot_model <- rma(yi=boot_Ti, vi=v[1:i], method="REML") # Should this still be vi?
          boot_ts <- boot_model$tau2
          boot_ts_list <- append(boot_ts_list, boot_ts)
        }
        a = alpha
        if (alpha == 0.1) {
          a = (-c/i)
        }
        if (alpha == 0.01) {
          a = -c/(i-1)
        }
        if (alpha == 0.001) {
          a = -c/sqrt(i)
        }
        expectation <- sum(exp(-a*boot_ts_list)) / 1000
        linex_ts <- max((-1/a)*log(expectation), 0) # This is the ts
      
        # Refit model with estimated ts set to linex (easier for calculations on UB)
        temp_model <- rma(yi = t[1:i], vi = v[1:i], tau2 = linex_ts)
        linex_ub <- predict(temp_model)$pi.ub[[1]]
        linex_lb <- predict(temp_model)$pi.lb[[1]]
        linex_mean <- temp_model$b
        linex_ts_list <- append(linex_ts_list, linex_ts)
        linex_ub_list <- append(linex_ub_list, linex_ub)
        linex_lb_list <- append(linex_lb_list, linex_lb)
        linex_mean_list <- append(linex_mean_list, linex_mean)
      
      }
      num_mas <- length(linex_ts_list)
      sim <- rep(j, num_mas)
      method <- rep("LINEX", num_mas)
      num_studies <- c(2:k)
      final_k <- rep(k, num_mas)
      true_ts <- rep(tau_squared, num_mas)
      true_mean <- rep(theta, num_mas)
      est_ts <- linex_ts_list
      est_ub <- linex_ub_list
      est_lb <- linex_lb_list
      est_mean <- linex_mean_list
    
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
}
close(pb)
write.csv(sim_data, file = "linex_tuning.csv", row.names = FALSE)

