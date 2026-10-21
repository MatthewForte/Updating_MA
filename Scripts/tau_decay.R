library(metafor)
library(patchwork)
library(tidyverse)

# Number of studies
k <- 20

# Empty dtf
sim_data <- data.frame(
  sim = numeric(),
  num_studies = numeric(),
  ts_true = numeric(),
  ts_est = numeric(), 
  ts_variance = numeric()
)

# Number of simulations
simulations <- 10000

# Setting random seed
set.seed(40)

# Sample sizes
n1 <- c(40, 40, 40)
n2 <- runif(k - 3, min=40, max=400)
n2 <- n2 / 2
n2 <- round(n2) * 2
n <- c(n1, n2)

# TS
varying_ts <- c(0, 0.01, 0.04, 0.16)

# Progress bar
pb <- txtProgressBar(min = 0, max = (simulations*length(varying_ts)), style = 3)
progress <- 1

for (sim_ts in varying_ts) {
  for (j in 1:simulations) {
    
    # Tau squared
    tau_squared <- sim_ts
    
    # True mean of distribution of effect sizes
    theta <- 0.1
    
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
    
    # Perform MA for 2:k studies, capture ts estimate
    ts_list <- c()
    mean_list <- c()
    ts_var_list <- c()
    for (i in 2:k) {
      model <- rma(yi=t[1:i], vi=v[1:i], method="REML", control = list(stepadj = 0.25, maxiter = 10000))
      ts_est <- model$tau2
      ts_list <- append(ts_list, ts_est)
      mean_est <- model$b[1]
      pval <- model$pval
      tau_var <- (model$se.tau2)^2
      ts_var_list <- append(ts_var_list, tau_var)
    }
    
    # Save row
    sim <- rep(j, k-1)
    num_studies <- 2:k
    ts_true <- rep(sim_ts, k-1)
    temp <- as.data.frame(cbind(sim, num_studies, ts_true, ts_list, ts_var_list))
    sim_data <- rbind(sim_data, temp)
    
    # Update progress bar
    setTxtProgressBar(pb, progress)
    progress <- progress + 1
  }
}

# Quick plots
plot_dat <- sim_data %>%
  group_by(num_studies, ts_true) %>%
  summarize(avg_var = mean(ts_var_list))


ggplot(plot_dat, aes(x=num_studies, y=avg_var, color=as.factor(ts_true))) +
  geom_point(size=1) +
  geom_line(linewidth=0.5) +
  labs(
    x = "Number of studies",
    y = expression("Average Variance of " * hat(tau)[plain(REML)]^2),
    color = expression("True " * tau^2)
  ) +
  theme_minimal(base_size=14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    panel.spacing = unit(1.2, "lines")
  )

# Temp dat
index <- c(2:20)
weight_1 <- 1/index
weight_3 <- 1/sqrt(index)
weight_4 <- 1/(index^2)
weight <- c(weight_1, weight_3, weight_4)
method <- c(rep("Weighted_1", length(index)), rep("Weighted_3", length(index)),
            rep("Weighted_4", length(index)))
num_studies <- c(index, index, index)


temp <- as.data.frame(cbind(num_studies, weight, method))
temp$num_studies <- as.numeric(temp$num_studies)
temp$weight <- as.numeric(temp$weight)

ggplot(temp, aes(x=num_studies, y=weight, color=method)) +
  geom_point(size=1) +
  geom_line(linewidth=0.5) +
  labs(
    x = "Number of studies",
    y = "Weight Value",
    color = "Weight Function"
  ) +
  scale_color_discrete(
    labels = c(
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ) |> lapply(function(x) parse(text = x))
  ) +
  theme_minimal(base_size=14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    panel.spacing = unit(1.2, "lines")
  )

