library(metadat)
library(metafor)
library(ggplot2)
library(patchwork)
library(metaplus)

####################################################################
# Adding one study at a time
####################################################################

# Set up data frame
example_data <- data.frame(
  method = character(),
  prior = numeric(),
  num_studies = numeric(),
  est_mean = numeric(), 
  est_ts = numeric(),
  est_ci_lb = numeric(),
  est_ci_ub = numeric(),
  est_pi_lb = numeric(),
  est_pi_ub = numeric()
)

# Read in data and calculate effect sizes
data <- metadat::dat.egger2001
data <- data[order(data$year, data$id), ]

es <- escalc(measure = "OR",
             ai = ai, bi = n1i - ai,
             ci = ci, di = n2i - ci,
             data = data)

# Number of studies
num_studies <- c(2:nrow(es))

# REML model estimation
og_ts_list <- c()
og_pi_ub_list <- c()
og_pi_lb_list <- c()
og_mean_list <- c()
og_ci_lb_list <- c()
og_ci_ub_list <- c()
for (i in 2:nrow(es)) {
  model <- rma(yi=yi[1:i], vi=vi[1:i], method="REML", data=es)
  
  og_ts <- model$tau2
  og_pi_ub <- predict(model)$pi.ub[[1]]
  og_pi_lb <- predict(model)$pi.lb[[1]]
  og_mean <- model$b
  og_ci_lb <- model$ci.lb
  og_ci_ub <- model$ci.ub
  
  og_ts_list <- append(og_ts_list, og_ts)
  og_pi_ub_list <- append(og_pi_ub_list, og_pi_ub)
  og_pi_lb_list <- append(og_pi_lb_list, og_pi_lb)
  og_mean_list <- append(og_mean_list, og_mean)
  og_ci_lb_list <- append(og_ci_lb_list, og_ci_lb)
  og_ci_ub_list <- append(og_ci_ub_list, og_ci_ub)
}

temp <- as.data.frame(cbind(
  rep("REML", nrow(es)-1), 
  rep(NA, nrow(es)-1), 
  c(2:nrow(es)), 
  og_mean_list, 
  og_ts_list, 
  og_ci_lb_list, 
  og_ci_ub_list, 
  og_pi_ub_list,
  og_pi_lb_list
))

colnames(temp) <- colnames(example_data)
example_data <- rbind(example_data, temp)

# Weighted model estimation
  # Weights: 1/sqrt(k), 1/(k-1)
  # Priors: 0.10, 0.20, 0.30
prior_list <- c(0.1, 0.2, 0.3)
weight_list <- c("1/sqrt(k)", "1/(k-1)")

for (weight in weight_list) {
  for (prior in prior_list) {
    
    w_ts_list <- c()
    w_pi_ub_list <- c()
    w_pi_lb_list <- c()
    w_mean_list <- c()
    w_ci_ub_list <- c()
    w_ci_lb_list <- c()
    
    for (i in 2:nrow(es)) {
      
      # construct c
      v_bar <- mean(es$vi[1:i])
      c <- (prior / (1 - prior))*v_bar
      
      # Get weight
      if (weight == "1/sqrt(k)") {
        w <- 1/sqrt(i)
      } else {
        w <- 1 / (i-1)
      }
      
      # Do the MA
      init_model <- rma(yi=yi[1:i], vi=vi[1:i], method="REML", data=es)
      init_ts <- init_model$tau2
      w_ts <- (w*c) + (1-w)*init_ts
      w_model <- rma(yi=yi[1:i], vi=vi[1:i], method="REML", tau2=w_ts, data=es)
      
      # Save stuff
      w_pi_ub <- predict(w_model)$pi.ub[[1]]
      w_pi_lb <- predict(w_model)$pi.lb[[1]]
      w_mean <- w_model$b
      w_ci_lb <- w_model$ci.lb
      w_ci_ub <- w_model$ci.ub
      w_ts_list <- append(w_ts_list, w_ts)
      w_pi_ub_list <- append(w_pi_ub_list, w_pi_ub)
      w_pi_lb_list <- append(w_pi_lb_list, w_pi_lb)
      w_mean_list <- append(w_mean_list, w_mean)
      w_ci_lb_list <- append(w_ci_lb_list, w_ci_lb)
      w_ci_ub_list <- append(w_ci_ub_list, w_ci_ub)
    }
    
    temp <- as.data.frame(cbind(
      rep(weight, nrow(es)-1), 
      rep(prior, nrow(es)-1), 
      c(2:nrow(es)), 
      w_mean_list, 
      w_ts_list, 
      w_ci_lb_list, 
      w_ci_ub_list, 
      w_pi_ub_list,
      w_pi_lb_list
    ))
    
    colnames(temp) <- colnames(example_data)
    example_data <- rbind(example_data, temp)
  }
}

# Clean datatypes
example_data$prior <- as.numeric(example_data$prior)
example_data$num_studies <- as.numeric(example_data$num_studies)
example_data$est_mean <- as.numeric(example_data$est_mean)
example_data$est_ts <- as.numeric(example_data$est_ts)
example_data$est_ci_lb <- as.numeric(example_data$est_ci_lb)
example_data$est_ci_ub <- as.numeric(example_data$est_ci_ub)
example_data$est_pi_lb <- as.numeric(example_data$est_pi_lb)
example_data$est_pi_ub <- as.numeric(example_data$est_pi_ub)


################### Plot 1: Tau squared estimate sequence

ggplot(example_data, aes(x = num_studies, y = est_ts, 
                         color = interaction(method, prior),
                         group = interaction(method, prior)
                         )) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
    x = "Number of studies",
    y = "Estimate",
    color = "Method",
    title = "Tau Squared Estimates"
  ) +
  guides(color = guide_legend(nrow = 1)) +
  theme_minimal(base_size = 14) +
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

# Now just the "best"
  # prior = 0.2
example_data <- example_data[example_data$prior %in% c(NA, 0.2), ]

ggplot(example_data, aes(x = num_studies, y = est_ts, 
                         color = interaction(method, prior),
                         group = interaction(method, prior)
                        )) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
    x = "Number of studies",
    y = "Estimate",
    color = "Method",
    title = "Tau Squared Estimates"
  ) +
  guides(color = guide_legend(nrow = 1)) +
  theme_minimal(base_size = 14) +
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

################### Plot 2: Prediction interval plots

temp1 <- example_data[example_data$method == "REML", ]
temp2 <- example_data[example_data$method == "1/sqrt(k)", ]
temp3 <- example_data[example_data$method == "1/(k-1)", ]

# Prediction interval plots
reml_pi_plot <- ggplot(temp1, aes(x=num_studies, y=est_mean)) +
  geom_pointrange(aes(ymin=est_pi_lb, ymax=est_pi_ub), color="red", size=0.1) +
  geom_pointrange(aes(ymin=est_ci_lb, ymax=est_ci_ub), color="black", size=0.1) +
  geom_hline(yintercept=0, color="green") + 
  labs(x="Number of Studies Included",
       y="95% Prediction Interval") +
  theme_minimal()

w1_pi_plot <- ggplot(temp2, aes(x=num_studies, y=est_mean)) +
  geom_pointrange(aes(ymin=est_pi_lb, ymax=est_pi_ub), color="red", size=0.1) +
  geom_pointrange(aes(ymin=est_ci_lb, ymax=est_ci_ub), color="black", size=0.1) +
  geom_hline(yintercept=0, color="green") + 
  labs(x="Number of Studies Included",
       y="95% Prediction Interval") +
  theme_minimal()

w2_pi_plot <- ggplot(temp3, aes(x=num_studies, y=est_mean)) +
  geom_pointrange(aes(ymin=est_pi_lb, ymax=est_pi_ub), color="red", size=0.1) +
  geom_pointrange(aes(ymin=est_ci_lb, ymax=est_ci_ub), color="black", size=0.1) +
  geom_hline(yintercept=0, color="green") + 
  labs(x="Number of Studies Included",
       y="95% Prediction Interval") +
  theme_minimal()

combined_plot <- reml_pi_plot / w1_pi_plot / w2_pi_plot

####################################################################
# Adding studies by publication year
####################################################################

# Order: 3 studies, 4, 5, 8, 10, 13, 14, 16
num_studies <- c(3, 4, 5, 8, 10, 13, 14, 16)

# Set up data frame
example_data <- data.frame(
  method = character(),
  prior = numeric(),
  num_studies = numeric(),
  est_mean = numeric(), 
  est_ts = numeric(),
  est_ci_lb = numeric(),
  est_ci_ub = numeric(),
  est_pi_lb = numeric(),
  est_pi_ub = numeric()
)

# REML model estimation
og_ts_list <- c()
og_pi_ub_list <- c()
og_pi_lb_list <- c()
og_mean_list <- c()
og_ci_lb_list <- c()
og_ci_ub_list <- c()
for (i in num_studies) {
  model <- rma(yi=yi[1:i], vi=vi[1:i], method="REML", data=es)
  
  og_ts <- model$tau2
  og_pi_ub <- predict(model)$pi.ub[[1]]
  og_pi_lb <- predict(model)$pi.lb[[1]]
  og_mean <- model$b
  og_ci_lb <- model$ci.lb
  og_ci_ub <- model$ci.ub
  
  og_ts_list <- append(og_ts_list, og_ts)
  og_pi_ub_list <- append(og_pi_ub_list, og_pi_ub)
  og_pi_lb_list <- append(og_pi_lb_list, og_pi_lb)
  og_mean_list <- append(og_mean_list, og_mean)
  og_ci_lb_list <- append(og_ci_lb_list, og_ci_lb)
  og_ci_ub_list <- append(og_ci_ub_list, og_ci_ub)
}

temp <- as.data.frame(cbind(
  rep("REML", length(num_studies)), 
  rep(NA, length(num_studies)), 
  num_studies, 
  og_mean_list, 
  og_ts_list, 
  og_ci_lb_list, 
  og_ci_ub_list, 
  og_pi_ub_list,
  og_pi_lb_list
))

colnames(temp) <- colnames(example_data)
example_data <- rbind(example_data, temp)

# Weighted model estimation
# Weights: 1/sqrt(k), 1/(k-1)
# Priors: 0.10, 0.20, 0.30
prior_list <- c(0.1, 0.2, 0.3)
weight_list <- c("1/sqrt(k)", "1/(k-1)")

for (weight in weight_list) {
  for (prior in prior_list) {
    
    w_ts_list <- c()
    w_pi_ub_list <- c()
    w_pi_lb_list <- c()
    w_mean_list <- c()
    w_ci_ub_list <- c()
    w_ci_lb_list <- c()
    
    for (i in num_studies) {
      
      # construct c
      v_bar <- mean(es$vi[1:i])
      c <- (prior / (1 - prior))*v_bar
      
      # Get weight
      if (weight == "1/sqrt(k)") {
        w <- 1/sqrt(i)
      } else {
        w <- 1 / (i-1)
      }
      
      # Do the MA
      init_model <- rma(yi=yi[1:i], vi=vi[1:i], method="REML", data=es)
      init_ts <- init_model$tau2
      w_ts <- (w*c) + (1-w)*init_ts
      w_model <- rma(yi=yi[1:i], vi=vi[1:i], method="REML", tau2=w_ts, data=es)
      
      # Save stuff
      w_pi_ub <- predict(w_model)$pi.ub[[1]]
      w_pi_lb <- predict(w_model)$pi.lb[[1]]
      w_mean <- w_model$b
      w_ci_lb <- w_model$ci.lb
      w_ci_ub <- w_model$ci.ub
      w_ts_list <- append(w_ts_list, w_ts)
      w_pi_ub_list <- append(w_pi_ub_list, w_pi_ub)
      w_pi_lb_list <- append(w_pi_lb_list, w_pi_lb)
      w_mean_list <- append(w_mean_list, w_mean)
      w_ci_lb_list <- append(w_ci_lb_list, w_ci_lb)
      w_ci_ub_list <- append(w_ci_ub_list, w_ci_ub)
    }
    
    temp <- as.data.frame(cbind(
      rep(weight, length(num_studies)), 
      rep(prior, length(num_studies)), 
      num_studies, 
      w_mean_list, 
      w_ts_list, 
      w_ci_lb_list, 
      w_ci_ub_list, 
      w_pi_ub_list,
      w_pi_lb_list
    ))
    
    colnames(temp) <- colnames(example_data)
    example_data <- rbind(example_data, temp)
  }
}

# Clean datatypes
example_data$prior <- as.numeric(example_data$prior)
example_data$num_studies <- as.numeric(example_data$num_studies)
example_data$est_mean <- as.numeric(example_data$est_mean)
example_data$est_ts <- as.numeric(example_data$est_ts)
example_data$est_ci_lb <- as.numeric(example_data$est_ci_lb)
example_data$est_ci_ub <- as.numeric(example_data$est_ci_ub)
example_data$est_pi_lb <- as.numeric(example_data$est_pi_lb)
example_data$est_pi_ub <- as.numeric(example_data$est_pi_ub)


################### Plot 1: Tau squared estimate sequence

ggplot(example_data, aes(x = num_studies, y = est_ts, 
                         color = interaction(method, prior),
                         group = interaction(method, prior)
                        )) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
    x = "Number of studies",
    y = "Estimate",
    color = "Method",
    title = "Tau Squared Estimates"
  ) +
  guides(color = guide_legend(nrow = 1)) +
  theme_minimal(base_size = 14) +
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

# Now just the "best"
# prior = 0.2
example_data <- example_data[example_data$prior %in% c(NA, 0.2), ]

ggplot(example_data, aes(x = num_studies, y = est_ts, 
                         color = interaction(method, prior),
                         group = interaction(method, prior)
                        )) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
    x = "Number of studies",
    y = "Estimate",
    color = "Method",
    title = "Tau Squared Estimates"
  ) +
  guides(color = guide_legend(nrow = 1)) +
  theme_minimal(base_size = 14) +
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

################### Plot 2: Prediction interval plots

temp1 <- example_data[example_data$method == "REML", ]
temp2 <- example_data[example_data$method == "1/sqrt(k)", ]
temp3 <- example_data[example_data$method == "1/(k-1)", ]

# Prediction interval plots
reml_pi_plot <- ggplot(temp1, aes(x=num_studies, y=est_mean)) +
  geom_pointrange(aes(ymin=est_pi_lb, ymax=est_pi_ub), color="red", size=0.1) +
  geom_pointrange(aes(ymin=est_ci_lb, ymax=est_ci_ub), color="black", size=0.1) +
  geom_hline(yintercept=0, color="green") + 
  labs(x="Number of Studies Included",
       y="95% Prediction Interval") +
  theme_minimal()

w1_pi_plot <- ggplot(temp2, aes(x=num_studies, y=est_mean)) +
  geom_pointrange(aes(ymin=est_pi_lb, ymax=est_pi_ub), color="red", size=0.1) +
  geom_pointrange(aes(ymin=est_ci_lb, ymax=est_ci_ub), color="black", size=0.1) +
  geom_hline(yintercept=0, color="green") + 
  labs(x="Number of Studies Included",
       y="95% Prediction Interval") +
  theme_minimal()

w2_pi_plot <- ggplot(temp3, aes(x=num_studies, y=est_mean)) +
  geom_pointrange(aes(ymin=est_pi_lb, ymax=est_pi_ub), color="red", size=0.1) +
  geom_pointrange(aes(ymin=est_ci_lb, ymax=est_ci_ub), color="black", size=0.1) +
  geom_hline(yintercept=0, color="green") + 
  labs(x="Number of Studies Included",
       y="95% Prediction Interval") +
  theme_minimal()

combined_plot <- reml_pi_plot / w1_pi_plot / w2_pi_plot

