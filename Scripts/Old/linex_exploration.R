library(metafor)
library(brms)
library(bayesmeta)
library(invgamma)

# Try making n much longer
# If a lot of variation in bootstrap ts, you wanna penalize more
# When does tau square stabilize in REML?
  # 
# SE of the estimate of tau
# Tipton Shuster 2017 has the method of moments estimator of tau squared
# Might need to go back to 20 studies and figure out when things stabilize


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

# Something to save linex bootstrap list
total_boot_ts_df <- data.frame(sim=numeric(), num_studies=numeric(), 
                               true_ts=numeric(), boot_ts=numeric(), stringsAsFactors=FALSE) 
  
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
      
      # Save the bootstrap list
      sim_to_add <- rep(j, 1000)
      k_to_add <- rep(i, 1000)
      ts_to_add <- rep(sim_ts, 1000)
      temp_df <- data.frame(sim_to_add, k_to_add, ts_to_add, boot_ts_list)
      colnames(temp_df) <- colnames(total_boot_ts_df)
      total_boot_ts_df <- rbind(total_boot_ts_df, temp_df)
      
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



# Clean up some data types on sim_data
sim_data$sim <- as.numeric(sim_data$sim)
sim_data$num_studies <- as.numeric(sim_data$num_studies)
sim_data$final_k <- as.numeric(sim_data$final_k)
sim_data$true_ts <- as.numeric(sim_data$true_ts)
sim_data$true_mean <- as.numeric(sim_data$true_mean)
sim_data$est_ts <- as.numeric(sim_data$est_ts)
sim_data$est_ub <- as.numeric(sim_data$est_ub)
sim_data$est_lb <- as.numeric(sim_data$est_lb )
sim_data$est_mean <- as.numeric(sim_data$est_mean)


#################################################
# Here is where im messing with the linex stuff
#################################################
sim_to_check <- 7

# Basically, look at when we have high ts
  # What we see if LINEX using a=-2 vs. RE
  # Now, recalculate LINEX for say, sim=1 using a bunch of different a's
  # Plot them all on the same line
temp <- sim_data[sim_data$true_ts == 0.16, ]
plot_data <- temp[temp$sim == sim_to_check, ]
ggplot(plot_data, aes(x=num_studies, y=est_ts, group=method, color=method)) +
  geom_line() +
  geom_hline(yintercept = 0.16)

# Recalculate linex using different a
  # Need to do this for each k
temp_linex <- total_boot_ts_df[total_boot_ts_df$true_ts == 0.16, ]
temp_linex <- temp_linex[temp_linex$sim == sim_to_check, ]

new_linex_ts <- c()
for (i in unique(temp_linex$num_studies)) {
  temp_loop <- temp_linex[temp_linex$num_studies == i, ]
  # Recalculate linex here
  new_a <- 3
  expectation <- sum(exp(-new_a*temp_loop$boot_ts)) / 1000
  linex_ts <- max((-1/new_a)*log(expectation), 0) # This is the ts
  new_linex_ts <- append(new_linex_ts, linex_ts)
}

# Now, add that as a line to the plot
num_studies <- 2:10
method <- rep("LINEX_new", 9)
new_linex_plot_data <- as.data.frame(cbind(num_studies, method, new_linex_ts))
new_linex_plot_data$num_studies <- as.numeric(new_linex_plot_data$num_studies)
new_linex_plot_data$new_linex_ts <- as.numeric(new_linex_plot_data$new_linex_ts)

  
ggplot(plot_data, aes(x=num_studies, y=est_ts, group=method, color=method)) +
  geom_line() +
  geom_hline(yintercept = 0.16) +
  geom_line(data = new_linex_plot_data, aes(x = num_studies, y = new_linex_ts), color = "green")
  

