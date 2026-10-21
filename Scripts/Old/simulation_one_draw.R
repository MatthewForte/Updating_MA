library(metafor)
library(ggplot2)
library(truncnorm)
library(ReIns)
library(patchwork)
library(tidyverse)

#############################################
# Simulation Setup
#############################################

# Setting random seed - 103, 130, 191, 225, 232, 279, 287, 295, 318, 359, 379, 400, 444, 462, 489
# 359, 489
# TS MSE - 130
set.seed(3)

# Number of studies
k <- 20

# Sample sizes (all equal for now)
n <- rep(100, k)

# Within study variances
v <- 1 / round(n)

# True mean of distribution of effect sizes
theta <- 0

# True variance of distribution of effect sizes
tau_squared <- 0.05

# True effect sizes for each of k studies
theta_i <- rnorm(k, mean=theta, sd=sqrt(tau_squared))

# Observed effect sizes for each of k studies
t <- c()
for (i in 1:k) {
  ti <- rnorm(1, mean=theta_i[i], sd=sqrt(v[i]))
  t <- append(t, ti)
}

############################################################################################################################
# Just 1 draw of all this
############################################################################################################################

###################################################
# Method 1: FE until k=5 then RE
###################################################

# Keep track of per each MA
mean_estimate <- c()
lower_bound_list <- c()
upper_bound_list <- c()
ci_lb_list <- c()
ci_ub_list <- c()
contains_neg_pos_list <- c()
pi_width_list <- c()
width_difference_list <- c()
lb_diff_list <- c()
ub_diff_list <- c()
percent_above_zero_list <- c()
ts_mse_list <- c()
ts_est_list <- c()
percent_contained_list <- c()


for (i in 2:k) {
  
  if (i < 5) {
    model <- rma(yi=t[1:i], vi=v[1:i], method="FE")
    
    estimate <- model$beta[1]
    mean_estimate <- append(mean_estimate, estimate)
    
    lower_bound_list <- append(lower_bound_list, 0)
    upper_bound_list <- append(upper_bound_list, 0)
    
    ci_lb <- model$ci.lb
    ci_ub <- model$ci.ub
    ci_lb_list <- append(ci_lb_list, ci_lb)
    ci_ub_list <- append(ci_ub_list, ci_ub)
    
    pi_lb <- model$beta[1]
    pi_ub <- model$beta[1]
    distribution_true_lb <- -1.96*sqrt(tau_squared)
    distribution_true_ub <- 1.96*sqrt(tau_squared)
    ts_est <- 0
    ts_var <- 0
    ts_mse <- (ts_est - tau_squared)^2
    
    # TS
    ts_mse_list <- append(ts_mse_list, ts_mse)
    ts_est_list <- append(ts_est_list, ts_est)
    
    # Does PI contain negative and positive values?
    contains_neg_pos <- 0
    contains_neg_pos_list <- append(contains_neg_pos_list, contains_neg_pos)
    
    # Width of PI
    pi_width <- 0
    pi_width_list <- append(pi_width_list, pi_width)
    
    # Width of PI compared to true width of distribution
    width_difference <- abs(distribution_true_ub - distribution_true_lb)
    width_difference_list <- append(width_difference_list, width_difference)
    
    # Difference between PI lb and true lb
    lb_diff <- abs(distribution_true_lb)
    lb_diff_list <- append(lb_diff_list, lb_diff)
    
    ## Difference between PI ub and true ub
    ub_diff <- abs(distribution_true_ub)
    ub_diff_list <- append(ub_diff_list, ub_diff)
    
    ## % of interval above 0
    percent_above_zero <- ifelse(pi_lb > 0, 1, 0)
    percent_above_zero_list <- append(percent_above_zero_list, percent_above_zero)
    
    ## How much of true distribution is contained in PI
    percent_contained <- 0
    percent_contained_list <- append(percent_contained_list, percent_contained)
    
  }
  else {
    model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
    estimate <- model$beta[1]
    mean_estimate <- append(mean_estimate, estimate)
    pi_lb <- predict(model)$pi.lb[[1]]
    pi_ub <- predict(model)$pi.ub[[1]]
    lower_bound_list <- append(lower_bound_list, pi_lb)
    upper_bound_list <- append(upper_bound_list, pi_ub)
    ci_lb <- model$ci.lb
    ci_ub <- model$ci.ub
    ci_lb_list <- append(ci_lb_list, ci_lb)
    ci_ub_list <- append(ci_ub_list, ci_ub)

    distribution_true_lb <- -1.96*sqrt(tau_squared)
    distribution_true_ub <- 1.96*sqrt(tau_squared)
    ts_est <- model$tau2
    ts_var <- (model$se.tau2)^2
    ts_mse <- (ts_est - tau_squared)^2 + ts_var
    
    # TS
    ts_mse_list <- append(ts_mse_list, ts_mse)
    ts_est_list <- append(ts_est_list, ts_est)
    
    ## Does PI contain negative and positive values?
    if (0 < pi_ub && 0 > pi_lb) {
      contains_neg_pos <- 1
    } else {
      contains_neg_pos <- 0
    }
    contains_neg_pos_list <- append(contains_neg_pos_list, contains_neg_pos)
    
    ## Width of the PI
    pi_width <- pi_ub - pi_lb
    pi_width_list <- append(pi_width_list, pi_width)
    
    ## Width of the PI compared to true width of distribution
    true_width <- distribution_true_ub - distribution_true_lb
    width_difference <- abs(pi_width - true_width)
    width_difference_list <- append(width_difference_list, width_difference)
    
    ## Difference between PI lb and true lb
    lb_diff <- abs(pi_lb - distribution_true_lb)
    lb_diff_list <- append(lb_diff_list, lb_diff)
    
    ## Difference between PI ub and true ub
    ub_diff <- abs(pi_ub - distribution_true_ub)
    ub_diff_list <- append(ub_diff_list, ub_diff)
    
    ## % of interval above 0
    percent_above_zero <- abs(0 - pi_ub) / abs(pi_width)
    percent_above_zero_list <- append(percent_above_zero_list, percent_above_zero)
    
    ## How much of true distribution is contained in PI
    overlap <- max(0, min(distribution_true_ub, pi_ub) - max(distribution_true_lb, pi_lb))
    percent_contained <- overlap / true_width
    percent_contained_list <- append(percent_contained_list, percent_contained)
  }
}

p <- 2
x <- seq(-p*sqrt(tau_squared), p*sqrt(tau_squared), length.out = 50)
y <- dnorm(x, 0, sqrt(tau_squared))
x0 <- 10
path1 <- data.frame(x = y + x0, y = x)

# Gather data
method <- rep("FE until k=5", k-1)
simulation_data <- as.data.frame(cbind(method, 2:k, mean_estimate, lower_bound_list, 
                                       upper_bound_list, ci_lb_list, ci_ub_list, contains_neg_pos_list, 
                                       pi_width_list, width_difference_list, lb_diff_list, ub_diff_list, 
                                       percent_above_zero_list, ts_mse_list, ts_est_list, percent_contained_list))
names(simulation_data)[2] <- "studies"
simulation_data$studies <- as.numeric(simulation_data$studies)
simulation_data$mean_estimate <- as.numeric(simulation_data$mean_estimate)
simulation_data$lower_bound_list <- as.numeric(simulation_data$lower_bound_list)
simulation_data$upper_bound_list <- as.numeric(simulation_data$upper_bound_list)
simulation_data$ci_lb_list <- as.numeric(simulation_data$ci_lb_list)
simulation_data$ci_ub_list <- as.numeric(simulation_data$ci_ub_list)
simulation_data$contains_neg_pos_list <- as.numeric(simulation_data$contains_neg_pos_list)
simulation_data$pi_width_list <- as.numeric(simulation_data$pi_width_list)
simulation_data$width_difference_list <- as.numeric(simulation_data$width_difference_list)
simulation_data$lb_diff_list <- as.numeric(simulation_data$lb_diff_list)
simulation_data$ub_diff_list <- as.numeric(simulation_data$ub_diff_list)
simulation_data$percent_above_zero_list <- as.numeric(simulation_data$percent_above_zero_list)
simulation_data$ts_mse_list <- as.numeric(simulation_data$ts_mse_list)
simulation_data$ts_est_list <- as.numeric(simulation_data$ts_est_list)
simulation_data$percent_contained_list <- as.numeric(simulation_data$percent_contained_list)

method_one_data <- simulation_data
method_one_data <- method_one_data[method_one_data$studies < 11, ]

# Plot
plot_method_one <- ggplot(method_one_data, aes(x=studies, y=mean_estimate)) +
  geom_pointrange(aes(ymin=lower_bound_list, ymax=upper_bound_list), color="#07bbc1", size=0.1) +
  geom_pointrange(aes(ymin=ci_lb_list, ymax=ci_ub_list), color="black", size=0.1) +
  geom_path(aes(x,y), data = path1, color = "#F8766D") +
  geom_hline(yintercept=0) + 
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  labs(title="Switch from FE to RE Model at k=5", x=NULL, y=NULL) +
  theme_minimal()

###################################################
# Method 2: FE until Q test significant, then RE
###################################################
# Keep track of per each MA
mean_estimate <- c()
lower_bound_list <- c()
upper_bound_list <- c()
ci_lb_list <- c()
ci_ub_list <- c()
tracker <- c()
contains_neg_pos_list <- c()
pi_width_list <- c()
width_difference_list <- c()
lb_diff_list <- c()
ub_diff_list <- c()
percent_above_zero_list <- c()
ts_mse_list <- c()
ts_est_list <- c()
percent_contained_list <- c()

# Perform MA with 2 studies, then 3, etc
for (i in 2:k) {
  if (i == 2) {
    test <- rma(yi=t[1:i], vi=v[1:i], method="REML")
    q_pval <- test$QEp
    if (q_pval < 0.05) { # DO RE
      model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
      estimate <- model$beta[1]
      mean_estimate <- append(mean_estimate, estimate)
      pi_lb <- predict(model)$pi.lb[[1]]
      pi_ub <- predict(model)$pi.ub[[1]]
      lower_bound_list <- append(lower_bound_list, pi_lb)
      upper_bound_list <- append(upper_bound_list, pi_ub)
      ci_lb <- model$ci.lb
      ci_ub <- model$ci.ub
      ci_lb_list <- append(ci_lb_list, ci_lb)
      ci_ub_list <- append(ci_ub_list, ci_ub)
      distribution_true_lb <- -1.96*sqrt(tau_squared)
      distribution_true_ub <- 1.96*sqrt(tau_squared)
      ts_est <- model$tau2
      ts_var <- (model$se.tau2)^2
      ts_mse <- (ts_est - tau_squared)^2 + ts_var
      ts_mse_list <- append(ts_mse_list, ts_mse)
      ts_est_list <- append(ts_est_list, ts_est)
      if (0 < pi_ub && 0 > pi_lb) {
        contains_neg_pos <- 1
      } else {
        contains_neg_pos <- 0
      }
      contains_neg_pos_list <- append(contains_neg_pos_list, contains_neg_pos)
      pi_width <- pi_ub - pi_lb
      pi_width_list <- append(pi_width_list, pi_width)
      true_width <- distribution_true_ub - distribution_true_lb
      width_difference <- abs(pi_width - true_width)
      width_difference_list <- append(width_difference_list, width_difference)
      lb_diff <- abs(pi_lb - distribution_true_lb)
      lb_diff_list <- append(lb_diff_list, lb_diff)
      ub_diff <- abs(pi_ub - distribution_true_ub)
      ub_diff_list <- append(ub_diff_list, ub_diff)
      percent_above_zero <- abs(0 - pi_ub) / abs(pi_width)
      percent_above_zero_list <- append(percent_above_zero_list, percent_above_zero)
      overlap <- max(0, min(distribution_true_ub, pi_ub) - max(distribution_true_lb, pi_lb))
      percent_contained <- overlap / true_width
      percent_contained_list <- append(percent_contained_list, percent_contained)
      tracker <- append(tracker, "RE")
    } else { # DO FE
      model <- rma(yi=t[1:i], vi=v[1:i], method="FE")
      estimate <- model$beta[1]
      mean_estimate <- append(mean_estimate, estimate)
      pi_lb <- model$beta[1]
      pi_ub <- model$beta[1]
      lower_bound_list <- append(lower_bound_list, pi_lb)
      upper_bound_list <- append(upper_bound_list, pi_ub)
      ci_lb <- model$ci.lb
      ci_ub <- model$ci.ub
      ci_lb_list <- append(ci_lb_list, ci_lb)
      ci_ub_list <- append(ci_ub_list, ci_ub)
      distribution_true_lb <- -1.96*sqrt(tau_squared)
      distribution_true_ub <- 1.96*sqrt(tau_squared)
      ts_est <- 0
      ts_var <- 0
      ts_mse <- (ts_est - tau_squared)^2
      ts_mse_list <- append(ts_mse_list, ts_mse)
      ts_est_list <- append(ts_est_list, ts_est)
      contains_neg_pos <- 0
      contains_neg_pos_list <- append(contains_neg_pos_list, contains_neg_pos)
      pi_width <- 0
      pi_width_list <- append(pi_width_list, pi_width)
      width_difference <- abs(distribution_true_ub - distribution_true_lb)
      width_difference_list <- append(width_difference_list, width_difference)
      lb_diff <- abs(distribution_true_lb)
      lb_diff_list <- append(lb_diff_list, lb_diff)
      ub_diff <- abs(distribution_true_ub)
      ub_diff_list <- append(ub_diff_list, ub_diff)
      percent_above_zero <- ifelse(pi_lb > 0, 1, 0)
      percent_above_zero_list <- append(percent_above_zero_list, percent_above_zero)
      percent_contained <- 0
      percent_contained_list <- append(percent_contained_list, percent_contained)
      tracker <- append(tracker, "FE")
    }
  } else {
    if (tracker[i-2] == "RE") { # DO RE
      model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
      estimate <- model$beta[1]
      mean_estimate <- append(mean_estimate, estimate)
      pi_lb <- predict(model)$pi.lb[[1]]
      pi_ub <- predict(model)$pi.ub[[1]]
      lower_bound_list <- append(lower_bound_list, pi_lb)
      upper_bound_list <- append(upper_bound_list, pi_ub)
      ci_lb <- model$ci.lb
      ci_ub <- model$ci.ub
      ci_lb_list <- append(ci_lb_list, ci_lb)
      ci_ub_list <- append(ci_ub_list, ci_ub)
      distribution_true_lb <- -1.96*sqrt(tau_squared)
      distribution_true_ub <- 1.96*sqrt(tau_squared)
      ts_est <- model$tau2
      ts_var <- (model$se.tau2)^2
      ts_mse <- (ts_est - tau_squared)^2 + ts_var
      ts_mse_list <- append(ts_mse_list, ts_mse)
      ts_est_list <- append(ts_est_list, ts_est)
      if (0 < pi_ub && 0 > pi_lb) {
        contains_neg_pos <- 1
      } else {
        contains_neg_pos <- 0
      }
      contains_neg_pos_list <- append(contains_neg_pos_list, contains_neg_pos)
      pi_width <- pi_ub - pi_lb
      pi_width_list <- append(pi_width_list, pi_width)
      true_width <- distribution_true_ub - distribution_true_lb
      width_difference <- abs(pi_width - true_width)
      width_difference_list <- append(width_difference_list, width_difference)
      lb_diff <- abs(pi_lb - distribution_true_lb)
      lb_diff_list <- append(lb_diff_list, lb_diff)
      ub_diff <- abs(pi_ub - distribution_true_ub)
      ub_diff_list <- append(ub_diff_list, ub_diff)
      percent_above_zero <- abs(0 - pi_ub) / abs(pi_width)
      percent_above_zero_list <- append(percent_above_zero_list, percent_above_zero)
      overlap <- max(0, min(distribution_true_ub, pi_ub) - max(distribution_true_lb, pi_lb))
      percent_contained <- overlap / true_width
      percent_contained_list <- append(percent_contained_list, percent_contained)
      tracker <- append(tracker, "RE")
    } else { # Check test again
      test <- rma(yi=t[1:i], vi=v[1:i], method="REML")
      q_pval <- test$QEp
      if (q_pval < 0.05) { # DO RE
        model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
        estimate <- model$beta[1]
        mean_estimate <- append(mean_estimate, estimate)
        pi_lb <- predict(model)$pi.lb[[1]]
        pi_ub <- predict(model)$pi.ub[[1]]
        lower_bound_list <- append(lower_bound_list, pi_lb)
        upper_bound_list <- append(upper_bound_list, pi_ub)
        ci_lb <- model$ci.lb
        ci_ub <- model$ci.ub
        ci_lb_list <- append(ci_lb_list, ci_lb)
        ci_ub_list <- append(ci_ub_list, ci_ub)
        distribution_true_lb <- -1.96*sqrt(tau_squared)
        distribution_true_ub <- 1.96*sqrt(tau_squared)
        ts_est <- model$tau2
        ts_var <- (model$se.tau2)^2
        ts_mse <- (ts_est - tau_squared)^2 + ts_var
        ts_mse_list <- append(ts_mse_list, ts_mse)
        ts_est_list <- append(ts_est_list, ts_est)
        if (0 < pi_ub && 0 > pi_lb) {
          contains_neg_pos <- 1
        } else {
          contains_neg_pos <- 0
        }
        contains_neg_pos_list <- append(contains_neg_pos_list, contains_neg_pos)
        pi_width <- pi_ub - pi_lb
        pi_width_list <- append(pi_width_list, pi_width)
        true_width <- distribution_true_ub - distribution_true_lb
        width_difference <- abs(pi_width - true_width)
        width_difference_list <- append(width_difference_list, width_difference)
        lb_diff <- abs(pi_lb - distribution_true_lb)
        lb_diff_list <- append(lb_diff_list, lb_diff)
        ub_diff <- abs(pi_ub - distribution_true_ub)
        ub_diff_list <- append(ub_diff_list, ub_diff)
        percent_above_zero <- abs(0 - pi_ub) / abs(pi_width)
        percent_above_zero_list <- append(percent_above_zero_list, percent_above_zero)
        overlap <- max(0, min(distribution_true_ub, pi_ub) - max(distribution_true_lb, pi_lb))
        percent_contained <- overlap / true_width
        percent_contained_list <- append(percent_contained_list, percent_contained)
        tracker <- append(tracker, "RE")
      } else { # DO FE
        model <- rma(yi=t[1:i], vi=v[1:i], method="FE")
        estimate <- model$beta[1]
        mean_estimate <- append(mean_estimate, estimate)
        pi_lb <- model$beta[1]
        pi_ub <- model$beta[1]
        lower_bound_list <- append(lower_bound_list, pi_lb)
        upper_bound_list <- append(upper_bound_list, pi_ub)
        ci_lb <- model$ci.lb
        ci_ub <- model$ci.ub
        ci_lb_list <- append(ci_lb_list, ci_lb)
        ci_ub_list <- append(ci_ub_list, ci_ub)
        distribution_true_lb <- -1.96*sqrt(tau_squared)
        distribution_true_ub <- 1.96*sqrt(tau_squared)
        ts_est <- 0
        ts_var <- 0
        ts_mse <- (ts_est - tau_squared)^2
        ts_mse_list <- append(ts_mse_list, ts_mse)
        ts_est_list <- append(ts_est_list, ts_est)
        contains_neg_pos <- 0
        contains_neg_pos_list <- append(contains_neg_pos_list, contains_neg_pos)
        pi_width <- 0
        pi_width_list <- append(pi_width_list, pi_width)
        width_difference <- abs(distribution_true_ub - distribution_true_lb)
        width_difference_list <- append(width_difference_list, width_difference)
        lb_diff <- abs(distribution_true_lb)
        lb_diff_list <- append(lb_diff_list, lb_diff)
        ub_diff <- abs(distribution_true_ub)
        ub_diff_list <- append(ub_diff_list, ub_diff)
        percent_above_zero <- ifelse(pi_lb > 0, 1, 0)
        percent_above_zero_list <- append(percent_above_zero_list, percent_above_zero)
        percent_contained <- 0
        percent_contained_list <- append(percent_contained_list, percent_contained)
        tracker <- append(tracker, "FE")
      }
    }
  }
}



# Gather data
method <- rep("FE until Q test significant", k-1)
simulation_data <- as.data.frame(cbind(method, 2:k, mean_estimate, lower_bound_list, 
                                       upper_bound_list, ci_lb_list, ci_ub_list, contains_neg_pos_list, 
                                       pi_width_list, width_difference_list, lb_diff_list, ub_diff_list, 
                                       percent_above_zero_list, ts_mse_list, ts_est_list, percent_contained_list))
names(simulation_data)[2] <- "studies"
simulation_data$studies <- as.numeric(simulation_data$studies)
simulation_data$mean_estimate <- as.numeric(simulation_data$mean_estimate)
simulation_data$lower_bound_list <- as.numeric(simulation_data$lower_bound_list)
simulation_data$upper_bound_list <- as.numeric(simulation_data$upper_bound_list)
simulation_data$ci_lb_list <- as.numeric(simulation_data$ci_lb_list)
simulation_data$ci_ub_list <- as.numeric(simulation_data$ci_ub_list)
simulation_data$contains_neg_pos_list <- as.numeric(simulation_data$contains_neg_pos_list)
simulation_data$pi_width_list <- as.numeric(simulation_data$pi_width_list)
simulation_data$width_difference_list <- as.numeric(simulation_data$width_difference_list)
simulation_data$lb_diff_list <- as.numeric(simulation_data$lb_diff_list)
simulation_data$ub_diff_list <- as.numeric(simulation_data$ub_diff_list)
simulation_data$percent_above_zero_list <- as.numeric(simulation_data$percent_above_zero_list)
simulation_data$ts_mse_list <- as.numeric(simulation_data$ts_mse_list)
simulation_data$ts_est_list <- as.numeric(simulation_data$ts_est_list)
simulation_data$percent_contained_list <- as.numeric(simulation_data$percent_contained_list)

method_two_data <- simulation_data
method_two_data <- method_two_data[method_two_data$studies < 11, ]

p <- 2
x <- seq(-p*sqrt(tau_squared), p*sqrt(tau_squared), length.out = 50)
y <- dnorm(x, 0, sqrt(tau_squared))
x0 <- 10
path1 <- data.frame(x = y + x0, y = x)

# Plot
plot_method_two <- ggplot(method_two_data, aes(x=studies, y=mean_estimate)) +
  geom_pointrange(aes(ymin=lower_bound_list, ymax=upper_bound_list), color="#07bbc1", size=0.1) +
  geom_pointrange(aes(ymin=ci_lb_list, ymax=ci_ub_list), color="black", size=0.1) +
  geom_path(aes(x,y), data = path1, color = "#F8766D") +
  geom_hline(yintercept=0) + 
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  labs(y="95% Prediction Interval",
       x=NULL,
       title="Switch from FE to RE Model When Q Test Significant") +
  theme_minimal()

###################################################
# Method 3: RE all the time
###################################################
# Keep track of per each MA
mean_estimate <- c()
lower_bound_list <- c()
upper_bound_list <- c()
ci_lb_list <- c()
ci_ub_list <- c()
contains_neg_pos_list <- c()
pi_width_list <- c()
width_difference_list <- c()
lb_diff_list <- c()
ub_diff_list <- c()
percent_above_zero_list <- c()
ts_mse_list <- c()
ts_est_list <- c()
percent_contained_list <- c()

for (i in 2:k) {
  model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
  estimate <- model$beta[1]
  mean_estimate <- append(mean_estimate, estimate)
  pi_lb <- predict(model)$pi.lb[[1]]
  pi_ub <- predict(model)$pi.ub[[1]]
  distribution_true_lb <- -1.96*sqrt(tau_squared)
  distribution_true_ub <- 1.96*sqrt(tau_squared)
  lower_bound_list <- append(lower_bound_list, pi_lb)
  upper_bound_list <- append(upper_bound_list, pi_ub)
  ci_lb <- model$ci.lb
  ci_ub <- model$ci.ub
  ci_lb_list <- append(ci_lb_list, ci_lb)
  ci_ub_list <- append(ci_ub_list, ci_ub)
  ts_est <- model$tau2
  ts_var <- (model$se.tau2)^2
  ts_mse <- (ts_est - tau_squared)^2 + ts_var
  
  # TS
  ts_mse_list <- append(ts_mse_list, ts_mse)
  ts_est_list <- append(ts_est_list, ts_est)
  
  ## How many PI contain positive and negative values?
  if (0 < pi_ub && 0 > pi_lb) {
    contains_neg_pos <- 1
  } else {
    contains_neg_pos <- 0
  }
  contains_neg_pos_list <- append(contains_neg_pos_list, contains_neg_pos)
  
  ## Width of the PI
  pi_width <- pi_ub - pi_lb
  pi_width_list <- append(pi_width_list, pi_width)
  
  ## Width of the PI compared to true width of distribution
  true_width <- distribution_true_ub - distribution_true_lb
  width_difference <- abs(pi_width - true_width)
  width_difference_list <- append(width_difference_list, width_difference)
  
  ## Difference between PI lb and true lb
  lb_diff <- abs(pi_lb - distribution_true_lb)
  lb_diff_list <- append(lb_diff_list, lb_diff)
  
  ## Difference between PI ub and true ub
  ub_diff <- abs(pi_ub - distribution_true_ub)
  ub_diff_list <- append(ub_diff_list, ub_diff)
  
  ## % of interval above 0
  percent_above_zero <- abs(0 - pi_ub) / abs(pi_width)
  percent_above_zero_list <- append(percent_above_zero_list, percent_above_zero)
  
  ## How much of true distribution is contained in PI
  overlap <- max(0, min(distribution_true_ub, pi_ub) - max(distribution_true_lb, pi_lb))
  percent_contained <- overlap / true_width
  percent_contained_list <- append(percent_contained_list, percent_contained)
}


# Gather data
method <- rep("RE always", k-1)
simulation_data <- as.data.frame(cbind(method, 2:k, mean_estimate, lower_bound_list, 
                                       upper_bound_list, ci_lb_list, ci_ub_list, contains_neg_pos_list, 
                                       pi_width_list, width_difference_list, lb_diff_list, ub_diff_list, 
                                       percent_above_zero_list, ts_mse_list, ts_est_list, percent_contained_list))
names(simulation_data)[2] <- "studies"
simulation_data$studies <- as.numeric(simulation_data$studies)
simulation_data$mean_estimate <- as.numeric(simulation_data$mean_estimate)
simulation_data$lower_bound_list <- as.numeric(simulation_data$lower_bound_list)
simulation_data$upper_bound_list <- as.numeric(simulation_data$upper_bound_list)
simulation_data$ci_lb_list <- as.numeric(simulation_data$ci_lb_list)
simulation_data$ci_ub_list <- as.numeric(simulation_data$ci_ub_list)
simulation_data$contains_neg_pos_list <- as.numeric(simulation_data$contains_neg_pos_list)
simulation_data$pi_width_list <- as.numeric(simulation_data$pi_width_list)
simulation_data$width_difference_list <- as.numeric(simulation_data$width_difference_list)
simulation_data$lb_diff_list <- as.numeric(simulation_data$lb_diff_list)
simulation_data$ub_diff_list <- as.numeric(simulation_data$ub_diff_list)
simulation_data$percent_above_zero_list <- as.numeric(simulation_data$percent_above_zero_list)
simulation_data$ts_mse_list <- as.numeric(simulation_data$ts_mse_list)
simulation_data$ts_est_list <- as.numeric(simulation_data$ts_est_list)
simulation_data$percent_contained_list <- as.numeric(simulation_data$percent_contained_list)

method_three_data <- simulation_data
method_three_data <- method_three_data[method_three_data$studies < 11, ]

p <- 2
x <- seq(-p*sqrt(tau_squared), p*sqrt(tau_squared), length.out = 50)
y <- dnorm(x, 0, sqrt(tau_squared))
x0 <- 10
path1 <- data.frame(x = y + x0, y = x)

# Plot
plot_method_three <- ggplot(method_three_data, aes(x=studies, y=mean_estimate)) +
  geom_pointrange(aes(ymin=lower_bound_list, ymax=upper_bound_list), color="#07bbc1", size=0.1) +
  geom_pointrange(aes(ymin=ci_lb_list, ymax=ci_ub_list), color="black", size=0.1) +
  geom_path(aes(x,y), data = path1, color = "#F8766D") +
  geom_hline(yintercept=0) + 
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  labs(x="Number of Studies Included",
       y=NULL,
       title="RE Model Used Throughout",
       caption=expression("Red curve is density of true effect size distribution shown to 2 standard deviations. "*tau^2*" set to 0.05")) +
  theme_minimal()




###################################################
# Method 4: FE always
###################################################

# Keep track of per each MA
mean_estimate <- c()
lower_bound_list <- c()
upper_bound_list <- c()
ci_lb_list <- c()
ci_ub_list <- c()
contains_neg_pos_list <- c()
pi_width_list <- c()
width_difference_list <- c()
lb_diff_list <- c()
ub_diff_list <- c()
percent_above_zero_list <- c()
ts_mse_list <- c()
ts_est_list <- c()
percent_contained_list <- c()


for (i in 2:k) {
  
  model <- rma(yi=t[1:i], vi=v[1:i], method="FE")
    
  estimate <- model$beta[1]
  mean_estimate <- append(mean_estimate, estimate)
    
  lower_bound_list <- append(lower_bound_list, 0)
  upper_bound_list <- append(upper_bound_list, 0)
    
  ci_lb <- model$ci.lb
  ci_ub <- model$ci.ub
  ci_lb_list <- append(ci_lb_list, ci_lb)
  ci_ub_list <- append(ci_ub_list, ci_ub)
    
  pi_lb <- model$beta[1]
  pi_ub <- model$beta[1]
  distribution_true_lb <- -1.96*sqrt(tau_squared)
  distribution_true_ub <- 1.96*sqrt(tau_squared)
  ts_est <- 0
  ts_var <- 0
  ts_mse <- (ts_est - tau_squared)^2
    
  # TS
  ts_mse_list <- append(ts_mse_list, ts_mse)
  ts_est_list <- append(ts_est_list, ts_est)
  
  # Does PI contain negative and positive values?
  contains_neg_pos <- 0
  contains_neg_pos_list <- append(contains_neg_pos_list, contains_neg_pos)
    
  # Width of PI
  pi_width <- 0
  pi_width_list <- append(pi_width_list, pi_width)
    
  # Width of PI compared to true width of distribution
  width_difference <- abs(distribution_true_ub - distribution_true_lb)
  width_difference_list <- append(width_difference_list, width_difference)
    
  # Difference between PI lb and true lb
  lb_diff <- abs(distribution_true_lb)
  lb_diff_list <- append(lb_diff_list, lb_diff)
  
  ## Difference between PI ub and true ub
  ub_diff <- abs(distribution_true_ub)
  ub_diff_list <- append(ub_diff_list, ub_diff)
    
  ## % of interval above 0
  percent_above_zero <- ifelse(pi_lb > 0, 1, 0)
  percent_above_zero_list <- append(percent_above_zero_list, percent_above_zero)
    
  ## How much of true distribution is contained in PI
  percent_contained <- 0
  percent_contained_list <- append(percent_contained_list, percent_contained)
}

method <- rep("FE always", k-1)
simulation_data <- as.data.frame(cbind(method, 2:k, mean_estimate, lower_bound_list, 
                                       upper_bound_list, ci_lb_list, ci_ub_list, contains_neg_pos_list, 
                                       pi_width_list, width_difference_list, lb_diff_list, ub_diff_list, 
                                       percent_above_zero_list, ts_mse_list, ts_est_list, percent_contained_list))
names(simulation_data)[2] <- "studies"
simulation_data$studies <- as.numeric(simulation_data$studies)
simulation_data$mean_estimate <- as.numeric(simulation_data$mean_estimate)
simulation_data$lower_bound_list <- as.numeric(simulation_data$lower_bound_list)
simulation_data$upper_bound_list <- as.numeric(simulation_data$upper_bound_list)
simulation_data$ci_lb_list <- as.numeric(simulation_data$ci_lb_list)
simulation_data$ci_ub_list <- as.numeric(simulation_data$ci_ub_list)
simulation_data$contains_neg_pos_list <- as.numeric(simulation_data$contains_neg_pos_list)
simulation_data$pi_width_list <- as.numeric(simulation_data$pi_width_list)
simulation_data$width_difference_list <- as.numeric(simulation_data$width_difference_list)
simulation_data$lb_diff_list <- as.numeric(simulation_data$lb_diff_list)
simulation_data$ub_diff_list <- as.numeric(simulation_data$ub_diff_list)
simulation_data$percent_above_zero_list <- as.numeric(simulation_data$percent_above_zero_list)
simulation_data$ts_mse_list <- as.numeric(simulation_data$ts_mse_list)
simulation_data$ts_est_list <- as.numeric(simulation_data$ts_est_list)
simulation_data$percent_contained_list <- as.numeric(simulation_data$percent_contained_list)

method_four_data <- simulation_data




###################################################
# Final Plot
###################################################

final_plot <- plot_method_one / plot_method_two / plot_method_three

###################################################
# More plots
###################################################

plot_data <- rbind(method_one_data, method_two_data, method_three_data)

# Percent Above 0 - Probability of Positive Effect
#ggplot(plot_data, aes(x=studies, y=percent_above_zero_list, group=method, color=method)) +
#  geom_point() + 
#  geom_line() +
#  xlim(2, 10) +
#  labs(x="Number of Studies", 
#       y="Probability of Positive Effect", 
#       caption="True probability of positive effect = 50%",
#       color="Method") +
#  theme_minimal()

temp <- plot_data[plot_data$method=="RE always", ]
new_plot_data <- rbind(temp, method_four_data)

new_plot_data[new_plot_data$method == "RE always", ]$method <- "RE"
new_plot_data[new_plot_data$method == "FE always", ]$method <- "FE"

# TS MSE
#ggplot(new_plot_data, aes(x=studies, y=ts_mse_list, group=method, color=method)) +
#  geom_point() + 
#  geom_line() +
#  xlim(2, 10) +
#  labs(x="Number of Studies", 
#       y=expression("MSE of "*tau^2),
#       color="Method") +
#  theme_minimal() +
#  theme(legend.position="left")

final_plot



