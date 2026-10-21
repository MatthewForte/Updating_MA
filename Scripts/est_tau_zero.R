library(metafor)
library(patchwork)
library(tidyverse)
######################################################################
# Idea
######################################################################
# In the paper, I want to say that estimating tau to be zero is bad
# Especially early on, because then that is mathematically a "model switch"
# Just empirically show this is a problem
  # Like, at k=2, we estimate tau^2 to be 0 X% of the time when it really isn't
  # This is a bigger problem at say k=2 than k=10 bc estimates are less stable and bc more research is coming
  # Also gives us a sense of when estimates start to stabilize

# Number of studies
k <- 10

# Empty dtf
sim_data <- data.frame(
  sim = numeric(),
  num_studies = numeric(),
  ts_true = numeric(),
  ts_est = numeric(), 
  sig_and_pos = numeric()
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
    sig_and_pos_list <- c()
    for (i in 2:k) {
      model <- rma(yi=t[1:i], vi=v[1:i], method="REML", control = list(stepadj = 0.25, maxiter = 10000))
      ts_est <- model$tau2
      ts_list <- append(ts_list, ts_est)
      mean_est <- model$b[1]
      pval <- model$pval
  
      # Is the estimate significant and positive?
      if (pval < 0.05 && mean_est > 0) {
        sig_and_pos <- 1
      } else {
        sig_and_pos <- 0
      }
      sig_and_pos_list <- append(sig_and_pos_list, sig_and_pos)
      
    }
  
    # Save row
    sim <- rep(j, k-1)
    num_studies <- 2:k
    ts_true <- rep(sim_ts, k-1)
    temp <- as.data.frame(cbind(sim, num_studies, ts_true, ts_list, sig_and_pos_list))
    sim_data <- rbind(sim_data, temp)
    
    # Update progress bar
    setTxtProgressBar(pb, progress)
    progress <- progress + 1
  }
}

######################################################################
# Get 0 counts and plot them
######################################################################

# Get counts of ts estimated as 0
zero_counts <- sim_data %>%
  group_by(num_studies, ts_true) %>%
  summarise(
    n_zero = sum(ts_list == 0)
  )

zero_counts$per_zero <- (zero_counts$n_zero / simulations) * 100
zero_counts$num_studies <- as.numeric(zero_counts$num_studies)

ggplot(zero_counts, aes(x=num_studies, y=per_zero)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~ ts_true, nrow = 1, labeller = label_bquote(tau^2 == .(ts_true))) +  
  labs(
    x = "Number of Studies",
    y = "% Estimated as 0"
    #title = expression("% of " * tau^2 * " Estimates Estimated as 0")
  ) +
  theme_minimal() + 
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines")
  ) + 
  xlim(c(0, 10)) +
  scale_x_continuous(breaks = 1:10)
  

