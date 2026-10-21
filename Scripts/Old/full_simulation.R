library(metafor)
library(brms)
library(bayesmeta)
library(invgamma)

# What I'm not varying
  # sample size/vi -> all set to 0.04
  # True mean effect -> set to 0.1
#############################################
# Simulation Setup
#############################################
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
      a = -2 # This may need to be tuned
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
close(pb)
write.csv(sim_data, file = "sim_data_1.csv", row.names = FALSE)




# Calculate bias and variance at each of the time points for each method





















# Bayesian - what if you put in a prior as just a number, not a distribution
#############################################
# Bayesian - Uniform
#############################################
# b_unif_estimate <- c()
# b_unif_ub <- c()
# b_unif_lb <- c()
# tau_prior <- function(t){dunif(t, min=0, max=1)}
# for (i in 2:k) {
#   # This is a normal prior for mu/theta
#   tryCatch({model <- bayesmeta(y=t[1:i], 
#                                sigma = sqrt(v[1:i]),
#                                label = study_id[1:i],
#                                mu.prior.mean=0, mu.prior.sd=50,
#                                tau.prior=tau_prior)
#   estimate <- model$summary[3, 2]
#   lb <- model$summary[5, 2]
#   ub <- model$summary[6, 2]
#   },
#   error = function(e) {
#     estimate <- NA
#     lb <- NA
#     ub <- NA
#   }
#   )
#   b_unif_estimate <- append(b_unif_estimate, estimate)
#   b_unif_ub <- append(b_unif_ub, ub)
#   b_unif_lb <- append(b_unif_lb, lb)
# }
# b_unif_mean_cor <- cor(2:k, b_unif_estimate, method="spearman")
# b_unif_ub_cor <- cor(2:k, b_unif_ub, method="spearman")
# b_unif_lb_cor <- cor(2:k, b_unif_lb, method="spearman")
# b_unif_mean_cor_list <- append(b_unif_mean_cor_list, b_unif_mean_cor)
# b_unif_ub_cor_list <- append(b_unif_ub_cor_list, b_unif_ub_cor)
# b_unif_lb_cor_list <- append(b_unif_lb_cor_list, b_unif_lb_cor)