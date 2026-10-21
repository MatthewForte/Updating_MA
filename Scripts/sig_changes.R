library(dplyr)

# Read in data
sim_data <- read.csv("Data/june_5_sim_data.csv")
sim_data <- sim_data[sim_data$method %in% c("REML", "Weighted_3"), ]
sim_data <- sim_data[sim_data$num_studies < 6, ]

# Set up dataframe
count_data <- data.frame(
  method = character(),
  prior = numeric(),
  true_ts = numeric(),
  count = numeric()
)

####################################################################################
# Percent of simulation repetitions with at least one significance reversal (mean)
####################################################################################

######### REML
reml_data <- sim_data[sim_data$method == "REML", ]

for (ts in unique(reml_data$true_ts)) {
  ts_data <- reml_data[reml_data$true_ts == ts, ]
  counter <- 0
  for (i in 1:4000) {
    loop_data <- ts_data[ts_data$sim == i, ]
    for (j in 1:3) {
      if ((loop_data$pval[j] < 0.05 && loop_data$pval[j+1] > 0.05) | (loop_data$pval[j] > 0.05 && loop_data$pval[j+1] < 0.05)) {
        counter <- counter + 1
        break
      }
    }
  }
  # Store method, prior (NA), true_ts, counter in a row and add to dataframe
  row_to_store <- c("REML", NA, ts, counter)
  names(row_to_store) <- names(count_data)
  count_data <- rbind(count_data, row_to_store)
}

######### Prior = 0.1
prior1 <- sim_data[sim_data$method != "REML", ]
prior1 <- prior1[prior1$prior == 0.1, ]
for (method in unique(prior1$method)) {
  method_data <- prior1[prior1$method == method, ]
  for (ts in unique(reml_data$true_ts)) {
    ts_data <- method_data[method_data$true_ts == ts, ]
    counter <- 0
    for (i in 1:4000) {
      loop_data <- ts_data[ts_data$sim == i, ]
      for (j in 1:3) {
        if ((loop_data$pval[j] < 0.05 && loop_data$pval[j+1] > 0.05) | (loop_data$pval[j] > 0.05 && loop_data$pval[j+1] < 0.05)) {
          counter <- counter + 1
          break
        }
      }
    }
    # Store method, prior (NA), true_ts, counter in a row and add to dataframe
    row_to_store <- c(method, 0.1, ts, counter)
    names(row_to_store) <- names(count_data)
    count_data <- rbind(count_data, row_to_store)
  }
}

######### Prior = 0.2
prior2 <- sim_data[sim_data$method != "REML", ]
prior2 <- prior2[prior2$prior == 0.2, ]
for (method in unique(prior1$method)) {
  method_data <- prior1[prior1$method == method, ]
  for (ts in unique(reml_data$true_ts)) {
    ts_data <- method_data[method_data$true_ts == ts, ]
    counter <- 0
    for (i in 1:4000) {
      loop_data <- ts_data[ts_data$sim == i, ]
      for (j in 1:18) {
        if ((loop_data$pval[j] < 0.05 && loop_data$pval[j+1] > 0.05) | (loop_data$pval[j] > 0.05 && loop_data$pval[j+1] < 0.05)) {
          counter <- counter + 1
          break
        }
      }
    }
    # Store method, prior (NA), true_ts, counter in a row and add to dataframe
    row_to_store <- c(method, 0.2, ts, counter)
    names(row_to_store) <- names(count_data)
    count_data <- rbind(count_data, row_to_store)
  }
}

######### Prior = 0.3
prior3 <- sim_data[sim_data$method != "REML", ]
prior3 <- prior3[prior3$prior == 0.3, ]
for (method in unique(prior1$method)) {
  method_data <- prior1[prior1$method == method, ]
  for (ts in unique(reml_data$true_ts)) {
    ts_data <- method_data[method_data$true_ts == ts, ]
    counter <- 0
    for (i in 1:4000) {
      loop_data <- ts_data[ts_data$sim == i, ]
      for (j in 1:18) {
        if ((loop_data$pval[j] < 0.05 && loop_data$pval[j+1] > 0.05) | (loop_data$pval[j] > 0.05 && loop_data$pval[j+1] < 0.05)) {
          counter <- counter + 1
          break
        }
      }
    }
    # Store method, prior (NA), true_ts, counter in a row and add to dataframe
    row_to_store <- c(method, 0.3, ts, counter)
    names(row_to_store) <- names(count_data)
    count_data <- rbind(count_data, row_to_store)
  }
}


names(count_data) <- c("method", "prior", "true_ts", "count")
count_data$percent <- (as.numeric(count_data$count) / 4000)*100

####################################################################################
# Percent of simulation repetitions with at least one yes to no sig reversal (mean)
####################################################################################

######### REML
reml_data <- sim_data[sim_data$method == "REML", ]

for (ts in unique(reml_data$true_ts)) {
  ts_data <- reml_data[reml_data$true_ts == ts, ]
  counter <- 0
  for (i in 1:4000) {
    loop_data <- ts_data[ts_data$sim == i, ]
    for (j in 1:3) {
      if ((loop_data$pval[j] > 0.05 && loop_data$pval[j+1] < 0.05)) {
        counter <- counter + 1
        break
      }
    }
  }
  # Store method, prior (NA), true_ts, counter in a row and add to dataframe
  row_to_store <- c("REML", NA, ts, counter)
  names(row_to_store) <- names(count_data)
  count_data <- rbind(count_data, row_to_store)
}

######### Prior = 0.1
prior1 <- sim_data[sim_data$method != "REML", ]
prior1 <- prior1[prior1$prior == 0.1, ]
for (method in unique(prior1$method)) {
  method_data <- prior1[prior1$method == method, ]
  for (ts in unique(reml_data$true_ts)) {
    ts_data <- method_data[method_data$true_ts == ts, ]
    counter <- 0
    for (i in 1:4000) {
      loop_data <- ts_data[ts_data$sim == i, ]
      for (j in 1:3) {
        if ((loop_data$pval[j] > 0.05 && loop_data$pval[j+1] < 0.05)) {
          counter <- counter + 1
          break
        }
      }
    }
    # Store method, prior (NA), true_ts, counter in a row and add to dataframe
    row_to_store <- c(method, 0.1, ts, counter)
    names(row_to_store) <- names(count_data)
    count_data <- rbind(count_data, row_to_store)
  }
}

######### Prior = 0.2
prior2 <- sim_data[sim_data$method != "REML", ]
prior2 <- prior2[prior2$prior == 0.2, ]
for (method in unique(prior1$method)) {
  method_data <- prior1[prior1$method == method, ]
  for (ts in unique(reml_data$true_ts)) {
    ts_data <- method_data[method_data$true_ts == ts, ]
    counter <- 0
    for (i in 1:4000) {
      loop_data <- ts_data[ts_data$sim == i, ]
      for (j in 1:18) {
        if ((loop_data$pval[j] > 0.05 && loop_data$pval[j+1] < 0.05)) {
          counter <- counter + 1
          break
        }
      }
    }
    # Store method, prior (NA), true_ts, counter in a row and add to dataframe
    row_to_store <- c(method, 0.2, ts, counter)
    names(row_to_store) <- names(count_data)
    count_data <- rbind(count_data, row_to_store)
  }
}

######### Prior = 0.3
prior3 <- sim_data[sim_data$method != "REML", ]
prior3 <- prior3[prior3$prior == 0.3, ]
for (method in unique(prior1$method)) {
  method_data <- prior1[prior1$method == method, ]
  for (ts in unique(reml_data$true_ts)) {
    ts_data <- method_data[method_data$true_ts == ts, ]
    counter <- 0
    for (i in 1:4000) {
      loop_data <- ts_data[ts_data$sim == i, ]
      for (j in 1:18) {
        if ((loop_data$pval[j] > 0.05 && loop_data$pval[j+1] < 0.05)) {
          counter <- counter + 1
          break
        }
      }
    }
    # Store method, prior (NA), true_ts, counter in a row and add to dataframe
    row_to_store <- c(method, 0.3, ts, counter)
    names(row_to_store) <- names(count_data)
    count_data <- rbind(count_data, row_to_store)
  }
}


names(count_data) <- c("method", "prior", "true_ts", "count")
count_data$percent <- (as.numeric(count_data$count) / 4000)*100

####################################################################################
# Percent of simulation repetitions with at least one no to yes sig reversal (mean)
####################################################################################

######### REML
reml_data <- sim_data[sim_data$method == "REML", ]

for (ts in unique(reml_data$true_ts)) {
  ts_data <- reml_data[reml_data$true_ts == ts, ]
  counter <- 0
  for (i in 1:4000) {
    loop_data <- ts_data[ts_data$sim == i, ]
    for (j in 1:3) {
      if ((loop_data$pval[j] < 0.05 && loop_data$pval[j+1] > 0.05)) {
        counter <- counter + 1
        break
      }
    }
  }
  # Store method, prior (NA), true_ts, counter in a row and add to dataframe
  row_to_store <- c("REML", NA, ts, counter)
  names(row_to_store) <- names(count_data)
  count_data <- rbind(count_data, row_to_store)
}

######### Prior = 0.1
prior1 <- sim_data[sim_data$method != "REML", ]
prior1 <- prior1[prior1$prior == 0.1, ]
for (method in unique(prior1$method)) {
  method_data <- prior1[prior1$method == method, ]
  for (ts in unique(reml_data$true_ts)) {
    ts_data <- method_data[method_data$true_ts == ts, ]
    counter <- 0
    for (i in 1:4000) {
      loop_data <- ts_data[ts_data$sim == i, ]
      for (j in 1:3) {
        if ((loop_data$pval[j] < 0.05 && loop_data$pval[j+1] > 0.05)) {
          counter <- counter + 1
          break
        }
      }
    }
    # Store method, prior (NA), true_ts, counter in a row and add to dataframe
    row_to_store <- c(method, 0.1, ts, counter)
    names(row_to_store) <- names(count_data)
    count_data <- rbind(count_data, row_to_store)
  }
}

names(count_data) <- c("method", "prior", "true_ts", "count")
count_data$percent <- (as.numeric(count_data$count) / 4000)*100



####################################################################################
# Average number of reversals per repitition
####################################################################################

sign_change_summary <- sim_data %>%
  mutate(
    sig = pval < 0.05,
    prior = ifelse(is.na(prior), "REML", as.character(prior))
  ) %>%
  arrange(method, prior, true_ts, sim, num_studies) %>%
  group_by(method, prior, true_ts, sim) %>%
  summarise(
    num_sign_changes = sum(sig != lag(sig), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(method, prior, true_ts) %>%
  summarise(
    avg_sign_changes = mean(num_sign_changes),
    sd_sign_changes = sd(num_sign_changes),
    prop_any_sign_change = mean(num_sign_changes > 0),
    .groups = "drop"
  )

sign_change_summary



sign_change_summary <- sim_data %>%
  mutate(
    sig = pval < 0.05,
    prior = ifelse(is.na(prior), "REML", as.character(prior))
  ) %>%
  arrange(method, prior, true_ts, sim, num_studies) %>%
  group_by(method, prior, true_ts, sim) %>%
  summarise(
    num_sig_to_nonsig = sum(lag(sig) == TRUE & sig == FALSE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(method, prior, true_ts) %>%
  summarise(
    avg_sig_to_nonsig = mean(num_sig_to_nonsig),
    sd_sig_to_nonsig = sd(num_sig_to_nonsig),
    prop_any_sig_to_nonsig = mean(num_sig_to_nonsig > 0),
    .groups = "drop"
  )


sign_change_summary <- sim_data %>%
  mutate(
    sig = pval < 0.05,
    prior = ifelse(is.na(prior), "REML", as.character(prior))
  ) %>%
  arrange(method, prior, true_ts, sim, num_studies) %>%
  group_by(method, prior, true_ts, sim) %>%
  summarise(
    num_nonsig_to_sig = sum(lag(sig) == FALSE & sig == TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(method, prior, true_ts) %>%
  summarise(
    avg_nonsig_to_sig = mean(num_nonsig_to_sig),
    sd_nonsig_to_sig = sd(num_nonsig_to_sig),
    prop_any_nonsig_to_sig = mean(num_nonsig_to_sig > 0),
    .groups = "drop"
  )






