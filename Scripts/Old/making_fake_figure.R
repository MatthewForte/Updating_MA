library(metafor)
library(ggplot2)
library(truncnorm)
library(ReIns)
library(patchwork)
library(tidyverse)

#############################################
# Simulation Setup
#############################################

# Setting random seed
set.seed(3)

# Number of studies
k <- 10

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


###################################################
# RE all the time
###################################################
# Keep track of per each MA
mean_estimate <- c()
lower_bound_list <- c()
upper_bound_list <- c()
se_list <- c()

for (i in 2:k) {
  model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
  estimate <- model$beta[1]
  mean_estimate <- append(mean_estimate, estimate)
  pi_lb <- predict(model)$pi.lb[[1]]
  pi_ub <- predict(model)$pi.ub[[1]]
  lower_bound_list <- append(lower_bound_list, pi_lb)
  upper_bound_list <- append(upper_bound_list, pi_ub)
  se <- model$se
  se_list <- append(se_list, se)
}


# Gather data
method <- rep("RE always", k-1)
simulation_data <- as.data.frame(cbind(method, 2:k, mean_estimate, lower_bound_list, 
                                       upper_bound_list))
names(simulation_data)[2] <- "studies"
simulation_data$studies <- as.numeric(simulation_data$studies)
simulation_data$mean_estimate <- as.numeric(simulation_data$mean_estimate)
simulation_data$lower_bound_list <- as.numeric(simulation_data$lower_bound_list)
simulation_data$upper_bound_list <- as.numeric(simulation_data$upper_bound_list)
p <- 2
x <- seq(-p*sqrt(tau_squared), p*sqrt(tau_squared), length.out = 50)
y <- dnorm(x, 0, sqrt(tau_squared))
x0 <- 10
path1 <- data.frame(x = y + x0, y = x)

# Plot
re_plot <- ggplot(simulation_data, aes(x=studies, y=mean_estimate)) +
  geom_pointrange(aes(ymin=lower_bound_list, ymax=upper_bound_list), color="#07bbc1", size=0.1) +
  geom_path(aes(x,y), data = path1, color = "#F8766D") +
  geom_hline(yintercept=0) + 
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x="Number of Studies Included",
       y=NULL) +
  theme_minimal()

###################################################
# RE all the time
###################################################
# Fake data
# True tau_squared is 0.05
# Final estimated tau_squared is 0.0246
init_ts <- 0.1
method <- rep("Fake", 9)
studies <- 2:k
mean_list <- simulation_data$mean_estimate
se <- se_list
ts_list <- c()
ts_list <- append(ts_list, init_ts)
for (i in 1:8) {
  new_ts <- ts_list[i] - ((init_ts - 0.0246) / 8)
  ts_list <- append(ts_list, new_ts)
}

lb_list <- mean_list - 2.306*sqrt(ts_list + se^2)
ub_list <- mean_list + 2.306*sqrt(ts_list + se^2)

fake_data <- as.data.frame(cbind(method, studies, mean_list, lb_list, ub_list))
fake_data$studies <- as.numeric(fake_data$studies)
fake_data$mean_list <- as.numeric(fake_data$mean_list)
fake_data$lb_list <- as.numeric(fake_data$lb_list)
fake_data$ub_list <- as.numeric(fake_data$ub_list)

fake_plot <- ggplot(fake_data, aes(x=studies, y=mean_list)) +
  geom_pointrange(aes(ymin=lb_list, ymax=ub_list), color="#07bbc1", size=0.1) +
  geom_path(aes(x,y), data = path1, color = "#F8766D") +
  geom_hline(yintercept=0) + 
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x="Number of Studies Included",
       y=NULL,
       caption=expression("Red curve is density of true effect size distribution shown to 2 standard deviations. "*tau^2*" set to 0.05")) +
  theme_minimal()


final_plot <- re_plot / fake_plot



