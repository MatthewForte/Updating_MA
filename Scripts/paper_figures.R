library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)
library(grid)
library(gt)
library(ggpattern)

# Load in simulation data
sim_data <- read.csv("Data/june_26_sim_data.csv")
sim_data <- sim_data[sim_data$prior != 0.9 | is.na(sim_data$prior), ]

# Labels
method_labels <- c(
  REML = "REML",
  Weighted_3 = "1/sqrt(k)",
  Weighted_1 = "1/k",
  Weighted_4 = "1/k^2"
)

###################################################################
# Plot 1: MSE across k, by prior and sim tau_squared
###################################################################

# Get MSE by prior, sim tau_squared
mse_data <- sim_data
mse_data$MSE_ts <- ((mse_data$est_ts - mse_data$true_ts)^2)*1000
mse_data <- mse_data %>%
  group_by(method, num_studies, true_ts, prior) %>%
  summarize(avg_mse = mean(MSE_ts))

# Get the plot data together
all_priors <- sort(unique(na.omit(mse_data$prior)))

re_rows <- mse_data %>%
  ungroup() %>%
  filter(method == "REML") %>%
  select(-prior) %>%
  crossing(prior = all_priors)

non_re_rows <- mse_data %>%
  ungroup() %>%
  filter(method != "REML")

plot_dat <- bind_rows(non_re_rows, re_rows) %>%
  mutate(
    num_studies = as.numeric(num_studies),
    true_ts = factor(true_ts),
    true_ts_lab = paste0("tau^2 == ", true_ts),
    prior_lab = paste0("I[P]^2 == ", prior * 100, " * '%'"),
    method = factor(method, levels = names(method_labels))
  )

# Make it only up to k=10
plot_dat <- plot_dat[plot_dat$num_studies <= 10, ]

# Make plot
ggplot(plot_dat, aes(x = num_studies, y = avg_mse, 
                     color = method, group = method, shape=method)) +
  geom_line(linewidth = 0.5, show.legend=FALSE) +
  geom_point(size = 1) +
  facet_grid(
    true_ts_lab ~ prior_lab,
    scales = "free_y",
    labeller = labeller(
      prior_lab = label_parsed,
      true_ts_lab = as_labeller(c(
        "tau^2 == 0"     = "atop(paste(I^2, ' \u2248 0%'), tau^2 == 0)",
        "tau^2 == 0.008" = "atop(paste(I^2, ' \u2248 30%'), tau^2 == 0.008)",
        "tau^2 == 0.03"  = "atop(paste(I^2, ' \u2248 60%'), tau^2 == 0.03)",
        "tau^2 == 0.16"  = "atop(paste(I^2, ' \u2248 90%'), tau^2 == 0.16)"
      ), label_parsed)
    )
  ) + 
  scale_color_manual(
    name = "Method",
    values = c(
      REML = "black",
      Weighted_1 = "#d95f02",
      Weighted_3 = "#e7298a",
      Weighted_4 = "#66a61e"
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  scale_shape_manual(
    name = "Method",
    values = c(
      REML = 16,
      Weighted_3 = 17,
      Weighted_1 = 18,
      Weighted_4 = 8
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  labs(
    x = "Number of studies",
    y = "MSE(x1000)",
    color = "Method"
    #title = "MSE(x1000) Across Number of Studies"
  ) +
  guides(
    shape = "none",
    color = guide_legend(
      nrow=1,
      override.aes=list(
        shape=c(16, 17, 18, 8),
        linetype=0,
        size=3
      )
    )
  ) +
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

###################################################################
# Plot 2: MAE across k, by prior and sim tau_squared
###################################################################

# Get MAE by prior, sim tau_squared
mae_data <- sim_data
mae_data$MAE_ts <- abs(mae_data$est_ts - mae_data$true_ts)*1000
mae_data <- mae_data %>%
  group_by(method, num_studies, true_ts, prior) %>%
  summarize(avg_mse = mean(MAE_ts))

# Get the plot data together
all_priors <- sort(unique(na.omit(mae_data$prior)))

re_rows <- mae_data %>%
  ungroup() %>%
  filter(method == "REML") %>%
  select(-prior) %>%
  crossing(prior = all_priors)

non_re_rows <- mae_data %>%
  ungroup() %>%
  filter(method != "REML")

plot_dat <- bind_rows(non_re_rows, re_rows) %>%
  mutate(
    num_studies = as.numeric(num_studies),
    true_ts = factor(true_ts),
    true_ts_lab = paste0("tau^2 == ", true_ts),
    prior_lab = paste0("I[P]^2 == ", prior * 100, " * '%'"),
    method = factor(method, levels = names(method_labels))
  )

# Make it only up to k=10
plot_dat <- plot_dat[plot_dat$num_studies <= 10, ]

# Make plot
ggplot(plot_dat, aes(x = num_studies, y = avg_mse, 
                     color = method, group = method, shape=method)) +
  geom_line(linewidth = 0.5, show.legend=FALSE) +
  geom_point(size = 1) +
  facet_grid(
    true_ts_lab ~ prior_lab,
    scales = "free_y",
    labeller = labeller(
      prior_lab = label_parsed,
      true_ts_lab = as_labeller(c(
        "tau^2 == 0"     = "atop(paste(I^2, ' \u2248 0%'), tau^2 == 0)",
        "tau^2 == 0.008" = "atop(paste(I^2, ' \u2248 30%'), tau^2 == 0.008)",
        "tau^2 == 0.03"  = "atop(paste(I^2, ' \u2248 60%'), tau^2 == 0.03)",
        "tau^2 == 0.16"  = "atop(paste(I^2, ' \u2248 90%'), tau^2 == 0.16)"
      ), label_parsed)
    )
  ) + 
  scale_color_manual(
    name = "Method",
    values = c(
      REML = "black",
      Weighted_1 = "#d95f02",
      Weighted_3 = "#e7298a",
      Weighted_4 = "#66a61e"
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  scale_shape_manual(
    name = "Method",
    values = c(
      REML = 16,
      Weighted_3 = 17,
      Weighted_1 = 18,
      Weighted_4 = 8
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  labs(
    x = "Number of studies",
    y = "MAE(x1000)",
    color = "Method"
    #title = "MAE(x1000) Across Number of Studies"
  ) +
  guides(
    shape = "none",
    color = guide_legend(
      nrow=1,
      override.aes=list(
        shape=c(16, 17, 18, 8),
        linetype=0,
        size=3
      )
    )
  ) +
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

###################################################################
# Plot 3: Bias across k, by prior and sim tau_squared
###################################################################
# Get bias by prior, sim tau_squared
bias_data <- sim_data
bias_data <- bias_data %>%
  group_by(method, num_studies, true_ts, prior) %>%
  summarize(avg_ts_est = mean(est_ts))
bias_data$bias_ts <- (bias_data$avg_ts_est - bias_data$true_ts)*1000

# Data handling
all_priors <- sort(unique(na.omit(bias_data$prior)))

re_rows <- bias_data %>%
  ungroup() %>%
  filter(method == "REML") %>%
  select(-prior) %>%
  crossing(prior = all_priors)

non_re_rows <- bias_data %>%
  ungroup() %>%
  filter(method != "REML")

plot_dat <- bind_rows(non_re_rows, re_rows) %>%
  mutate(
    num_studies = as.numeric(num_studies),
    true_ts = factor(true_ts),
    true_ts_lab = paste0("tau^2 == ", true_ts),
    prior_lab = paste0("I[P]^2 == ", prior * 100, " * '%'"),
    method = factor(method, levels = names(method_labels)),
  )

# Make it only up to k=10
plot_dat <- plot_dat[plot_dat$num_studies <= 10, ]

# Make plot
ggplot(plot_dat, aes(x = num_studies, y = bias_ts, 
                     color = method, group = method, shape=method)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  facet_grid(
    true_ts_lab ~ prior_lab,
    scales = "free_y",
    labeller = labeller(
      prior_lab = label_parsed,
      true_ts_lab = as_labeller(c(
        "tau^2 == 0"     = "atop(paste(I^2, ' \u2248 0%'), tau^2 == 0)",
        "tau^2 == 0.008" = "atop(paste(I^2, ' \u2248 30%'), tau^2 == 0.008)",
        "tau^2 == 0.03"  = "atop(paste(I^2, ' \u2248 60%'), tau^2 == 0.03)",
        "tau^2 == 0.16"  = "atop(paste(I^2, ' \u2248 90%'), tau^2 == 0.16)"
      ), label_parsed)
    )
  ) + 
  scale_color_manual(
    name = "Method",
    values = c(
      REML = "black",
      Weighted_1 = "#d95f02",
      Weighted_3 = "#e7298a",
      Weighted_4 = "#66a61e"
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  scale_shape_manual(
    name = "Method",
    values = c(
      REML = 16,
      Weighted_3 = 17,
      Weighted_1 = 18,
      Weighted_4 = 8
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  labs(
    x = "Number of studies",
    y = "Bias(x1000)",
    color = "Method"
    #title = "Bias Across Number of Studies"
  ) +
  guides(
    shape = "none",
    color = guide_legend(
      nrow=1,
      override.aes=list(
        shape=c(16, 17, 18, 8),
        linetype=0,
        size=3
      )
    )
  ) +
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

###################################################################
# Plot 4: Average MSE gain over REML
###################################################################

# Get MSE by prior, sim tau_squared
mse_data <- sim_data
mse_data$MSE_ts <- ((mse_data$est_ts - mse_data$true_ts)^2)*1000
mse_data <- mse_data %>%
  group_by(method, num_studies, true_ts, prior) %>%
  summarize(avg_mse = mean(MSE_ts))

re_mse <- mse_data %>%
  filter(method == "REML") %>%
  select(true_ts, num_studies, re_mse = avg_mse)

mse_improvement_data <- mse_data %>%
  left_join(re_mse, by = c("true_ts", "num_studies"))

mse_improvement_data <- mse_improvement_data %>%
  mutate(
    pct_improvement = 100 * (re_mse - avg_mse) / re_mse
  )

mse_improvement_data <- mse_improvement_data %>%
  filter(method.x != "REML")

mse_improvement_data <- mse_improvement_data %>%
  arrange(method.x, prior, true_ts, num_studies)

avg_improve_data <- mse_improvement_data %>%
  ungroup() %>%
  arrange(method.x, prior, true_ts, num_studies) %>%
  group_by(method.x, prior, true_ts) %>%
  mutate(
    avg_improvement = cumsum(pct_improvement) / row_number(),
    cumulative_improvement = cumsum(pct_improvement)
  ) %>%
  ungroup()

plot_dat <- avg_improve_data %>%
  filter(num_studies %in% c(5, 10, 20))

plot_dat <- plot_dat %>%
  mutate(
    method = factor(method.x, levels = c("Weighted_1", "Weighted_3", "Weighted_4")),
    num_studies = factor(num_studies, levels = c(5, 10, 20)),
    true_ts_lab = paste0("tau^2 == ", true_ts),
    prior_lab = paste0("I[P]^2 == ", prior * 100, " * '%'")
  )

plot_dat$method.x <- factor(
  plot_dat$method.x,
  levels = c("Weighted_3", "Weighted_1", "Weighted_4")
)

# Make plot
ggplot(plot_dat, aes(x = num_studies, y = avg_improvement, 
                     fill = method.x, pattern=method.x)) +
  geom_col_pattern(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "grey30",              # softer than black
    pattern_fill = "grey30",       # softer pattern color
    pattern_colour = "grey30",
    pattern_density = 0.2,         # ↓ less busy
    pattern_spacing = 0.05,        # ↑ more airy
    pattern_key_scale_factor = 0.8
  ) +
  facet_grid(
    true_ts_lab ~ prior_lab,
    scales = "free_y",
    labeller = labeller(
      prior_lab = label_parsed,
      true_ts_lab = as_labeller(c(
        "tau^2 == 0"     = "atop(paste(I^2, ' \u2248 0%'), tau^2 == 0)",
        "tau^2 == 0.008" = "atop(paste(I^2, ' \u2248 30%'), tau^2 == 0.008)",
        "tau^2 == 0.03"  = "atop(paste(I^2, ' \u2248 60%'), tau^2 == 0.03)",
        "tau^2 == 0.16"  = "atop(paste(I^2, ' \u2248 90%'), tau^2 == 0.16)"
      ), label_parsed)
    )
  ) +
  scale_fill_manual(
    name = "Method",
    values = c(
      Weighted_3 = "#F4A3C4",  # light pink
      Weighted_1 = "#FDC086",  # pastel orange
      Weighted_4 = "#B2DF8A"   # light green
    ),
    breaks = c("Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      "1/sqrt(k)",
      "1/k",
      "1/k^2"
    ))
  ) +
  scale_pattern_manual(
    name = "Method",
    values = c(
      Weighted_3 = "stripe",
      Weighted_1 = "circle",
      Weighted_4 = "none"
    ),
    breaks = c("Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      "1/sqrt(k)",
      "1/k",
      "1/k^2"
    ))
  ) +
  labs(
    x = "Number of Studies",
    y = "Average of MSE(x1000) Percent Improvements",
    fill = "Method"
    #title = "MSE(x1000) Improvement of Shrinkage Estimator Over REML"
  ) +
  guides(
    fill = "none",
    pattern = guide_legend(
      nrow = 1,
      override.aes = list(
        fill = c("#F4A3C4", "#FDC086", "#B2DF8A"),
        pattern = c("stripe", "circle", "none"),
        colour = "grey30"
      )
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
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

###################################################################
# Plot 6: MSSD + MASD
###################################################################

sequence_data <- sim_data

# Can only go to 20 because we jump to 30, 50 so not a sequence
  # IE if we keep those it will treat 20 to 30 studies as 1 "consecutive" jump
  # Could be true if we add 10 studies in batch but not how we are framing things
sequence_data <- sequence_data[sequence_data$num_studies <= 20, ]
true_ts_values <- unique(sequence_data$true_ts)
method_values <- unique(sequence_data$method)
prior_values <- unique(sequence_data$prior)
prior_values <- prior_values[!is.na(prior_values)]
num_sims <- max(sequence_data$sim)

# Want to get the MSSD for each sequence
# Subset to a true_ts, method, prior
  # Then loop through each sim and get the metric
  # Will end up with a list of metrics?

metrics_data <- data.frame(
  true_ts = numeric(),
  method = character(),
  prior = numeric(),
  sim = numeric(),
  stop_points = numeric(),
  mssd = numeric(),
  masd = numeric()
)


sequence_stop_points <- c(3, 8)

for (ts in true_ts_values) {
  for (method in method_values) {
    for (prior in prior_values) {
      
      # If method is REML, and we already have a row with REML and ts, go to next iteration
      if(method == "REML" & any(metrics_data$method == "REML" & metrics_data$true_ts == ts)) {
        next
      }
      # Subset data
      temp <- sequence_data[sequence_data$true_ts == ts, ]
      temp <- temp[temp$method == method, ]
      if (method != "REML") {
        temp <- temp[temp$prior == prior, ]
      }
      # For each simulation run, get the MSSD of the sequence
      for (s in 1:num_sims) {
        temp1 <- temp[temp$sim == s, ]
        # Now, get the MSSD of the sequence, but actually cumulatively
        # So first get it with a sequence of 2, then 3, etc
        diffs <- diff(temp1$est_ts)
        mssd_list <- c()
        masd_list <- c()
        for (i in sequence_stop_points) {
          mssd <- mean(diffs[1:i]^2)
          mssd_list <- append(mssd_list, mssd)
          masd <- mean(abs(diffs[1:i]))
          masd_list <- append(masd_list, masd)
          
        }
        # Put them somewhere
        rows_to_add <- as.data.frame(cbind(
          rep(ts, 4),
          rep(method, 4),
          rep(prior, 4),
          rep(s, 4), 
          sequence_stop_points, 
          mssd_list,
          masd_list
        ))
        colnames(rows_to_add) <- colnames(metrics_data)
        metrics_data <- rbind(metrics_data, rows_to_add)
      }
    }
  }
}

# Clean data types
metrics_data$true_ts <- as.numeric(metrics_data$true_ts)
metrics_data$prior <- as.numeric(metrics_data$prior)
metrics_data$sim <- as.numeric(metrics_data$sim)
metrics_data$stop_points <- as.numeric(metrics_data$stop_points)
metrics_data$mssd <- as.numeric(metrics_data$mssd)
metrics_data$masd <- as.numeric(metrics_data$masd)

# Average across simulations
plot_dat <- metrics_data %>%
  group_by(true_ts, method, prior, stop_points) %>%
  summarize(avg_mssd = mean(mssd), avg_masd = mean(masd))
  
plot_dat$stop_points <- plot_dat$stop_points + 2

# Duplicate REML data for each prior - data in there now is for prior 0.1
reml_rows <- plot_dat[plot_dat$method == "REML", ]
reml1 <- reml_rows[reml_rows$true_ts == 0, ]
reml2 <- reml_rows[reml_rows$true_ts == 0.008, ]
reml3 <- reml_rows[reml_rows$true_ts == 0.03, ]
reml4 <- reml_rows[reml_rows$true_ts == 0.16, ]


# Format
plot_dat <- plot_dat %>% 
  mutate(
    true_ts = factor(true_ts),
    true_ts_lab = paste0("tau^2 == ", true_ts),
    prior_lab = paste0("I[P]^2 == ", prior * 100, " * '%'"),
    method = factor(method, levels = names(method_labels)),
    avg_mssd = avg_mssd*1000,
    avg_masd = avg_masd*1000
  )

# Make plot
ggplot(plot_dat, aes(x = stop_points, y = avg_mssd, 
                     color = method, group = method, shape=method)) +
  geom_line(linewidth = 0.5, show.legend=FALSE) +
  geom_point(size = 1) +
  facet_grid(
    true_ts_lab ~ prior_lab,
    scales = "free_y",
    labeller = labeller(
      prior_lab = label_parsed,
      true_ts_lab = as_labeller(c(
        "tau^2 == 0"     = "atop(paste(I^2, ' \u2248 0%'), tau^2 == 0)",
        "tau^2 == 0.008" = "atop(paste(I^2, ' \u2248 30%'), tau^2 == 0.008)",
        "tau^2 == 0.03"  = "atop(paste(I^2, ' \u2248 60%'), tau^2 == 0.03)",
        "tau^2 == 0.16"  = "atop(paste(I^2, ' \u2248 90%'), tau^2 == 0.16)"
      ), label_parsed)
    )
  ) + 
  scale_color_manual(
    name = "Method",
    values = c(
      REML = "black",
      Weighted_1 = "#d95f02",
      Weighted_3 = "#e7298a",
      Weighted_4 = "#66a61e"
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  scale_shape_manual(
    name = "Method",
    values = c(
      REML = 16,
      Weighted_3 = 17,
      Weighted_1 = 18,
      Weighted_4 = 8
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  labs(
    x = "Number of studies",
    y = "Average MSSD(x1000)",
    color = "Method"
    #title = "Average MSSD(x1000) Across Number of Studies"
  ) +
  guides(
    shape = "none",
    color = guide_legend(
      nrow=1,
      override.aes=list(
        shape=c(16, 17, 18, 8),
        linetype=0,
        size=3
      )
    )
  ) +
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
  
# Make plot 2
ggplot(plot_dat, aes(x = stop_points, y = avg_masd, 
                     color = method, group = method, shape=method)) +
  geom_line(linewidth = 0.5, show.legend=FALSE) +
  geom_point(size = 1) +
  facet_grid(
    true_ts_lab ~ prior_lab,
    scales = "free_y",
    labeller = labeller(
      prior_lab = label_parsed,
      true_ts_lab = as_labeller(c(
        "tau^2 == 0"     = "atop(paste(I^2, ' \u2248 0%'), tau^2 == 0)",
        "tau^2 == 0.008" = "atop(paste(I^2, ' \u2248 30%'), tau^2 == 0.008)",
        "tau^2 == 0.03"  = "atop(paste(I^2, ' \u2248 60%'), tau^2 == 0.03)",
        "tau^2 == 0.16"  = "atop(paste(I^2, ' \u2248 90%'), tau^2 == 0.16)"
      ), label_parsed)
    )
  ) + 
  scale_color_manual(
    name = "Method",
    values = c(
      REML = "black",
      Weighted_1 = "#d95f02",
      Weighted_3 = "#e7298a",
      Weighted_4 = "#66a61e"
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  scale_shape_manual(
    name = "Method",
    values = c(
      REML = 16,
      Weighted_3 = 17,
      Weighted_1 = 18,
      Weighted_4 = 8
    ),
    breaks = c("REML", "Weighted_3", "Weighted_1", "Weighted_4"),
    labels = parse(text = c(
      REML = "REML",
      Weighted_3 = "1/sqrt(k)",
      Weighted_1 = "1/k",
      Weighted_4 = "1/k^2"
    ))
  ) +
  labs(
    x = "Number of studies",
    y = "Average MASD(x1000)",
    color = "Method"
    #title = "Average MASD(x1000) Across Number of Studies"
  ) +
  guides(
    shape = "none",
    color = guide_legend(
      nrow=1,
      override.aes=list(
        shape=c(16, 17, 18, 8),
        linetype=0,
        size=3
      )
    )
  ) +
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





