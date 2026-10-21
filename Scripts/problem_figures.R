library(metafor)
library(ggplot2)
library(tidyverse)

################################################
# Generate some data
################################################

# Setting random seed
set.seed(40)

# Number of studies
k <- 5

repeat {

  # True mean of distribution of effect sizes
  theta <- 0.05

  # True variance of distribution of effect sizes
  tau_squared <- 0.12

  # Study level standard errors
  se <- runif(k, 0.08, 0.15)

  # True effect sizes for each of k studies
  theta_i <- rnorm(k, mean=theta, sd=sqrt(tau_squared))

  # Observed effects
  yi <- rnorm(k, mean=theta_i, sd=se)

  results <- list()

  for (i in 2:k) {
    yi_k <- yi[1:i]
    se_k <- se[1:i]
  
    ce_model <- rma(yi=yi_k, vi=se_k, method="FE")
    re_model <- rma(yi=yi_k, vi=se_k, method="REML")
  
    results[[i]] <- data.frame(
      num_studies = c(i, i),
      method = c("CE", "RE"),
      mean = c(coef(ce_model), coef(re_model)),
      ci_lb = c(ce_model$ci.lb, re_model$ci.lb),
      ci_ub = c(ce_model$ci.ub, re_model$ci.ub),
      pi_lb = c(predict(ce_model)$pi.lb[[1]], predict(re_model)$pi.lb[[1]]),
      pi_ub = c(predict(ce_model)$pi.ub[[1]], predict(re_model)$pi.ub[[1]])
    )
  }
  plot_dat <- bind_rows(results)
  
  check5 <- plot_dat %>% filter(num_studies == 5)
  ce_sig <- with(check5 %>% filter(method == "CE"), ci_ub < 0 | ci_lb > 0)
  re_sig <- with(check5 %>% filter(method == "RE"), ci_ub < 0 | ci_lb > 0)
  
  
  
  check1 <- plot_dat %>% filter(num_studies == 1)
  ce_sig_1 <- with(check1 %>% filter(method == "CE"), ci_ub < 0 | ci_lb > 0)
  check2 <- plot_dat %>% filter(num_studies == 2)
  ce_sig_2 <- with(check2 %>% filter(method == "CE"), ci_ub < 0 | ci_lb > 0)
  check3 <- plot_dat %>% filter(num_studies == 3)
  ce_sig_3 <- with(check3 %>% filter(method == "CE"), ci_ub < 0 | ci_lb > 0)
  check4 <- plot_dat %>% filter(num_studies == 4)
  ce_sig_4 <- with(check4 %>% filter(method == "CE"), ci_ub < 0 | ci_lb > 0)
  
  if (isTRUE(ce_sig) && !isTRUE(re_sig) && isTRUE(ce_sig_4)) break
}

plot_dat <- plot_dat %>%
  filter(!(num_studies < 5 & method == "RE"))

plot_dat <- plot_dat %>%
  filter(num_studies < 6)

################################################
# Make Plot 1
################################################

theta_true <- 0.05
tau2_true  <- 0.12
tau_true   <- sqrt(tau2_true)

true_dist <- data.frame(
  y = seq(theta_true - 3*tau_true, theta_true + 3*tau_true, length.out = 300)
) |>
  mutate(
    density = dnorm(y, mean = theta_true, sd = tau_true),
    x = 5.35 + density / max(density) * 0.45
  )

re_pi_dat <- subset(plot_dat, method == "RE" & num_studies == 5)

# Plot
ggplot(plot_dat, aes(x=as.numeric(num_studies), y=mean, color=method, shape=method)) +
  geom_hline(yintercept=0) + 
  geom_errorbar(
    aes(ymin=ci_lb, ymax=ci_ub),
    position = position_dodge(width=0.35),
    width=0.15, 
    linewidth=0.7
  ) +
  geom_point(
    position=position_dodge(width=0.35),
    size=2.5
  ) +
  labs(
    x = "Number of Studies Included",
    y = "Estimated Effect Size",
    color = "Model",
    shape = "Model"
    #title = "Inference changes when switching from CE to RE at k = 5"
  ) +
  theme_minimal(base_size = 14) +
  geom_segment(
    data=re_pi_dat,
    aes(
      x = 5.0875, xend=5.0875,
      y = ci_ub, yend=pi_ub
    ),
    color = "#00BFC4",
    linetype = "dashed",
    linewidth=0.7,
    inherit.aes=FALSE
  ) +
  geom_segment(
    data = re_pi_dat,
    aes(
      x = 5.0875, xend = 5.0875,
      y = pi_lb, yend = ci_lb
    ),
    color = "#00BFC4",
    linetype = "dashed",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  scale_linetype_manual(
    values = c(
      "95% CI" = "solid",
      "95% PI" = "dashed"
    ),
    name = "Interval"
  ) +
  geom_line(
    data = true_dist,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = theta_true,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  coord_cartesian(xlim = c(1.8, 6.35))
  
################################################
# Make Plot 2
################################################

find_oscillating_example <- function(max_tries = 50000,
                                     k_max = 10,
                                     k_start = 2,
                                     mu = 0.12,
                                     tau2 = 0.02,
                                     min_oscillations = 4) {
  
  best <- NULL
  
  for (s in 1:max_tries) {
    set.seed(s)
    
    # More variable vi helps produce unstable PI width
    vi <- sample(c(0.005, 0.01, 0.02, 0.05, 0.08, 0.12), 
                 size = k_max, replace = TRUE)
    
    theta <- rnorm(k_max, mu, sqrt(tau2))
    yi <- rnorm(k_max, theta, sqrt(vi))
    
    out <- lapply(k_start:k_max, function(k) {
      fit <- rma(yi[1:k], vi[1:k], method = "REML", control = list(
        maxiter = 1000,
        stepadj = 0.5,
        threshold = 1e-8
      ))
      
      pred <- predict(fit)
      
      if (fit$tau2 == 0) {
        pilb <- as.numeric(fit$b)
        piub <- as.numeric(fit$b)
      } else {
        pilb <- pred$pi.lb[[1]]
        piub <- pred$pi.ub[[1]]
      }
      
      
      data.frame(
        k = k,
        est = as.numeric(fit$b),
        se = fit$se,
        tau2_hat = fit$tau2,
        pi_lb = pilb,
        pi_ub = piub
      )
    }) |> 
      bind_rows() |>
      mutate(
        # conclusion based on prediction interval
        pi_sig = pi_lb > 0 | pi_ub < 0
      )
    
    n_osc <- sum(diff(as.integer(out$pi_sig)) != 0)
    
    if (n_osc >= min_oscillations) {
      out$seed <- s
      out$n_oscillations <- n_osc
      out$pattern <- paste(as.integer(out$pi_sig), collapse = "")
      out$yi <- I(list(yi))
      out$vi <- I(list(vi))
      return(out)
    }
    
    # store best in case nothing reaches min_oscillations
    if (is.null(best) || n_osc > best$n_oscillations[1]) {
      out$seed <- s
      out$n_oscillations <- n_osc
      out$pattern <- paste(as.integer(out$pi_sig), collapse = "")
      out$yi <- I(list(yi))
      out$vi <- I(list(vi))
      best <- out
    }
  }
  
  best
}

mean_to_use = 0.12
ts_to_use <- 0.03

ex <- find_oscillating_example(
  max_tries = 20000,
  k_max = 10,
  k_start = 2,
  mu = mean_to_use,
  tau2 = ts_to_use,
  min_oscillations = 4
)

ex$x_plot <- factor(ex$k)
levels(ex$x_plot) <- c(levels(ex$x_plot), "True")

# Plot it
ggplot(ex, aes(x = x_plot, y = est)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(
    aes(ymin = pi_lb, ymax = pi_ub), 
    width = 0.15
  ) +
  geom_point(size = 2.5) +
  labs(
    x = "Number of Studies",
    y = "Estimated Effect Size"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

################################################
# Re-analyze using our method
################################################

t <- ex$yi[[1]]
v <- ex$vi[[1]]

est_list <- c()
ci_ub_list <- c()
ci_lb_list <- c()
for (i in 2:10) {
  
  reml_model <- rma(yi=t[1:i], vi=v[1:i], method="REML")
  reml_ts <- reml_model$tau2
  
  weight <- 1/sqrt(i)
  
  v_bar <- mean(v[1:i])
  
  prior <- (0.1/(1-0.1))*v_bar
  
  ts_weighted <- (weight*prior) + ((1-weight)*reml_ts)
  
  final_model <- rma(yi=t[1:i], vi=v[1:i], tau2=ts_weighted)
  
  
  # Need to save estimate, ci_ub, ci_lb
  est <- reml_model$b[1]
  ci_lb <- final_model$ci.lb
  ci_ub <- final_model$ci.ub
  est_list <- append(est_list, est)
  ci_ub_list <- append(ci_ub_list, ci_ub)
  ci_lb_list <- append(ci_lb_list, ci_lb)
}

reanalyzed_data <- as.data.frame(cbind(est_list, ci_ub_list, ci_lb_list, 2:10))
names(reanalyzed_data) <- c("est", "ci_ub", "ci_lb", "k")

reanalyzed_data$x_plot <- factor(reanalyzed_data$k)
levels(reanalyzed_data$x_plot) <- c(levels(reanalyzed_data$x_plot), "True")
vi_vec <- c(0.02, 0.02, 0.08, 0.02, 0.10, 0.03, 0.09, 0.02, 0.06, 0.03)
true_k <- 10
true_tau2 <- 0.12
true_mu <- 0.32
true_se <- sqrt(1 / sum(1 / (vi_vec[1:true_k] + true_tau2)))

true_dat <- data.frame(
  x_plot = factor("True", levels = levels(ex$x_plot)),
  mean = true_mu,
  ci_lb = true_mu - 1.96 * true_se,
  ci_ub = true_mu + 1.96 * true_se
)



ggplot(reanalyzed_data, aes(x = x_plot, y = est)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(
    aes(ymin = ci_lb, ymax = ci_ub, color = sig), 
    width = 0.15
  ) +
  geom_point(aes(color = sig), size = 2.5) +
  geom_point(
    data = true_dat,
    aes(x = x_plot, y = mean),
    inherit.aes = FALSE,
    size = 2.8,
    shape = 18,
    color = "black"
  ) +
  scale_color_manual(
    values = c("FALSE" = "red3", "TRUE" = "darkgreen")
  ) +
  labs(
    x = "Number of Studies",
    y = "Estimated Effect Size"
  ) +
  theme_minimal() +
  theme(legend.position = "none")








