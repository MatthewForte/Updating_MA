library(dplyr)
library(ggplot2)
library(gt)
library(glue)
library(purrr)
library(tidyr)
library(plotly)
library(crosstalk)
library(htmltools)
library(stringr)

# Notes:
  # For some runs, each of the 9 MAs all estimate ts as 0 and we cannot calculate a correlation
    # Thus, our final average correlation remove these runs, meaning they are not quite an average of 1000 runs
    # and maybe an average of, say, 990 runs, but only for spearman and kendall correlation
  # For isotonic R^2, it's the same thing
  # Sometimes it seems like people will set the R2 to 0 or 1 though (python defaults do this)

######################################################################
# Calculating Metrics Over Simulation Runs
######################################################################

# Source functions
source("metric_function.R")

# Load in simulation data
sim_data <- read.csv("Data/weighted_re_ind_studies_isquared.csv")

# Get rid of the max stuff
sim_data <- sim_data[sim_data$method != "Max_Weighted_1", ]
sim_data <- sim_data[sim_data$method != "Max_Weighted_2", ]
sim_data <- sim_data[sim_data$method != "Max_Weighted_3", ]
sim_data <- sim_data[sim_data$method != "Max_Weighted_4", ]

# Get values I varied in the simulation
true_ts_values <- unique(sim_data$true_ts)
final_k_values <- c(5, 10, 20, 30, 50)
method_values <- unique(sim_data$method)
prior_values <- unique(na.omit(sim_data$prior))

# Create empty dataframe
averaged_metrics <- data.frame(
  method = character(),
  true_ts = character(),
  final_k = character(),
  prior = numeric(),
  avg_ts_spearman = character(),
  avg_ts_kendall = character(),
  avg_ts_mono_ratio = character(),
  avg_ts_isotonic_r2 = character(),
  avg_ts_dcor = character(),
  avg_ts_var = character(),
  avg_ts_mssd = character(),
  avg_ts_masd = character(),
  avg_ts_maj_raw = character(),
  avg_ts_maj_norm = character(),
  avg_ub_spearman = character(),
  avg_ub_kendall = character(),
  avg_ub_mono_ratio = character(),
  avg_ub_isotonic_r2 = character(),
  avg_ub_dcor = character(),
  avg_ub_var = character(),
  avg_ub_mssd = character(),
  avg_ub_masd = character(),
  avg_ub_maj_raw = character(),
  avg_ub_maj_norm = character()
)

for (true_ts in true_ts_values) {
  for (final_k in final_k_values) {
    for (method in method_values) {
      for (prior in prior_values) {
        # For each combination of final_k, method, true_ts, get a dataframe
        temp <- sim_data[sim_data$method == method, ]
        temp <- temp[temp$num_studies <= final_k, ]
        temp <- temp[temp$true_ts == true_ts, ]
        if (method != "RE") {
          temp <- temp[temp$prior == prior, ]
        }
        # If method is RE and true_ts and final_k value are already in averaged metrics, go to the next iteration of the loop
      
        # For each simulation, you have a sequence of ts's and ub's
          # Get the metrics for each sequence
        ts_spearman_list <- c()
        ts_kendall_list <- c()
        ts_mono_ratio_list <- c()
        ts_isotonic_r2_list <- c()
        ts_dcor_list <- c()
        ts_var_list <- c()
        ts_mssd_list <- c()
        ts_masd_list <- c()
        ts_maj_raw_list <- c()
        ts_maj_norm_list <- c()
        ub_spearman_list <- c()
        ub_kendall_list <- c()
        ub_mono_ratio_list <- c()
        ub_isotonic_r2_list <- c()
        ub_dcor_list <- c()
        ub_var_list <- c()
        ub_mssd_list <- c()
        ub_masd_list <- c()
        ub_maj_raw_list <- c()
        ub_maj_norm_list <- c()
      
        for (i in 1:max(temp$sim)) {
          temp1 <- temp[temp$sim == i, ]
        
          metrics <- get_all_metrics(temp1$est_ts)
          spearman <- metrics$spearman
          kendall <- metrics$kendall
          mono_ratio <- metrics$mono_ratio
          isotonic_r2 <- metrics$isotonic_r2
          dcor <- metrics$dcor
          var <- metrics$var
          mssd <- metrics$mssd
          masd <- metrics$masd
          maj_raw <- metrics$maj_raw[[1]]
          maj_norm <- metrics$maj_norm[[1]]
        
          ts_spearman_list <- append(ts_spearman_list, spearman)
          ts_kendall_list <- append(ts_kendall_list, kendall)
          ts_mono_ratio_list <- append(ts_mono_ratio_list, mono_ratio)
          ts_isotonic_r2_list <- append(ts_isotonic_r2_list, isotonic_r2)
          ts_dcor_list <- append(ts_dcor_list, dcor)
          ts_var_list <- append(ts_var_list, var)
          ts_mssd_list <- append(ts_mssd_list, mssd)
          ts_masd_list <- append(ts_masd_list, masd)
          ts_maj_raw_list <- append(ts_maj_raw_list, maj_raw)
          ts_maj_norm_list <- append(ts_maj_norm_list, maj_norm)
        
          metrics <- get_all_metrics(temp1$est_ub)
          spearman <- metrics$spearman
          kendall <- metrics$kendall
          mono_ratio <- metrics$mono_ratio
          isotonic_r2 <- metrics$isotonic_r2
          dcor <- metrics$dcor
          var <- metrics$var
          mssd <- metrics$mssd
          masd <- metrics$masd
          maj_raw <- metrics$maj_raw[[1]]
          maj_norm <- metrics$maj_norm[[1]]
        
          ub_spearman_list <- append(ub_spearman_list, spearman)
          ub_kendall_list <- append(ub_kendall_list, kendall)
          ub_mono_ratio_list <- append(ub_mono_ratio_list, mono_ratio)
          ub_isotonic_r2_list <- append(ub_isotonic_r2_list, isotonic_r2)
          ub_dcor_list <- append(ub_dcor_list, dcor)
          ub_var_list <- append(ub_var_list, var)
          ub_mssd_list <- append(ub_mssd_list, mssd)
          ub_masd_list <- append(ub_masd_list, masd)
          ub_maj_raw_list <- append(ub_maj_raw_list, maj_raw)
          ub_maj_norm_list <- append(ub_maj_norm_list, maj_norm)
        }
        # Clean up isotonic R2, set NaN to NA
        ts_isotonic_r2_list[is.nan(ts_isotonic_r2_list)] <- NA
      
        # Now, I have a list of metrics, 1 for each simulation run
        # I want to average them and get an average metric
        avg_ts_spearman <- mean(ts_spearman_list, na.rm=TRUE)
        avg_ts_kendall <- mean(ts_kendall_list, na.rm=TRUE)
        avg_ts_mono_ratio <- mean(ts_mono_ratio_list)
        avg_ts_isotonic_r2 <- mean(ts_isotonic_r2_list, na.rm=TRUE)
        avg_ts_dcor <- mean(ts_dcor_list)
        avg_ts_var <- mean(ts_var_list)
        avg_ts_mssd <- mean(ts_mssd_list)
        avg_ts_masd <- mean(ts_masd_list)
        avg_ts_maj_raw <- mean(ts_maj_raw_list, na.rm=TRUE)
        avg_ts_maj_norm <- mean(ts_maj_norm_list, na.rm=TRUE)
      
        avg_ub_spearman <- mean(ub_spearman_list)
        avg_ub_kendall <- mean(ub_kendall_list)
        avg_ub_mono_ratio <- mean(ub_mono_ratio_list)
        avg_ub_isotonic_r2 <- mean(ub_isotonic_r2_list)
        avg_ub_dcor <- mean(ub_dcor_list)
        avg_ub_var <- mean(ub_var_list)
        avg_ub_mssd <- mean(ub_mssd_list)
        avg_ub_masd <- mean(ub_masd_list)
        avg_ub_maj_raw <- mean(ub_maj_raw_list)
        avg_ub_maj_norm <- mean(ub_maj_norm_list)
      
        # Put a row in a data frame
        row_to_insert <- c(method, true_ts, final_k, prior, avg_ts_spearman, avg_ts_kendall, avg_ts_mono_ratio,
                      avg_ts_isotonic_r2, avg_ts_dcor, avg_ts_var, avg_ts_mssd, avg_ts_masd, 
                      avg_ts_maj_raw, avg_ts_maj_norm, avg_ub_spearman, avg_ub_kendall, avg_ub_mono_ratio,
                      avg_ub_isotonic_r2, avg_ub_dcor, avg_ub_var, avg_ub_mssd, avg_ub_masd, avg_ub_maj_raw,
                      avg_ub_maj_norm)
      
        row_to_insert <- as.data.frame(t(row_to_insert), stringAsFactors=FALSE)
        colnames(row_to_insert) <- colnames(averaged_metrics)
        averaged_metrics <- rbind(averaged_metrics, row_to_insert)
      }
    }
  }
}

# Clean data types
cols_to_stay <- c("method")
averaged_metrics[, !names(averaged_metrics) %in% cols_to_stay] <- 
  lapply(averaged_metrics[, !names(averaged_metrics) %in% cols_to_stay], as.numeric)

######################################################################
# Make tables
######################################################################

############################## Metric tables (averaged_data)

method_labels <- c(
  RE         = "RE",
  Weighted_1 = "1/k",
  Weighted_2 = "1/(k-1)",
  Weighted_3 = "1/sqrt(k)",
  Weighted_4 = "1/(k^2)"
)

metric_cols_ts <- averaged_metrics %>%
  select(starts_with("avg_ts_")) %>%
  names()

metric_cols_ub <- averaged_metrics %>%
  select(starts_with("avg_ub_")) %>%
  names()

metrics_table_list <- list()
i <- 1
for (ts in true_ts_values) {
  for (k in final_k_values) {
    for (pr in prior_values) {
      
      # TS table
      temp <- averaged_metrics %>%
        filter(true_ts == ts, prior == pr, final_k == k) %>%
        select(true_ts, prior, final_k, method, all_of(metric_cols_ts))
      
      tbl_dat <- temp %>%
        select(method, all_of(metric_cols_ts)) %>%
        pivot_longer(-method, names_to = "metric", values_to = "value") %>%
        mutate(
          metric = str_remove(metric, "^avg_ts_"),
          metric = str_replace_all(metric, "_", " ")
        ) %>%
        pivot_wider(names_from = method, values_from = value)
      
      max_metrics <- c("mono ratio", "isotonic r2", "dcor")
      method_cols <- setdiff(names(tbl_dat), "metric")
      hi <- tbl_dat %>%
        mutate(.row = row_number()) %>%
        rowwise() %>%
        mutate(
          .winners = list({
            vals <- c_across(all_of(method_cols))
            
            # if all missing, highlight nothing
            if (all(is.na(vals))) character(0) else {
              want_max <- metric %in% max_metrics
              target <- if (want_max) max(vals, na.rm = TRUE) else min(vals, na.rm = TRUE)
              
              # return ALL tied winners
              method_cols[!is.na(vals) & vals == target]
            }
          })
        ) %>%
        ungroup() %>%
        select(.row, metric, .winners) %>%
        tidyr::unnest(.winners) %>%
        rename(method = .winners)
      
      table <- tbl_dat %>%
        gt(rowname_col = "metric") %>%
        cols_label(.list = method_labels) %>%
        fmt_number(columns = all_of(method_cols), decimals = 4) %>%
        tab_header(
          title = "TS Metrics",
          subtitle = glue("Number of studies: {k} | Tau-squared: {ts} | Prior: {pr}")
        )
      
      for (k in seq_len(nrow(hi))) {
        table <- table %>%
          tab_style(
            style = cell_fill(color = "#E6F4EA"),
            locations = cells_body(
              columns = hi$method[k],
              rows    = hi$.row[k]
            )
          )
      }
      
      key <- paste0("iter_", i)
      metrics_table_list[[key]] <- table
      
      # UB table
      temp <- averaged_metrics %>%
        filter(true_ts == ts, prior == pr, final_k == k) %>%
        select(true_ts, prior, final_k, method, all_of(metric_cols_ub))
      
      tbl_dat <- temp %>%
        select(method, all_of(metric_cols_ub)) %>%
        pivot_longer(-method, names_to = "metric", values_to = "value") %>%
        mutate(
          metric = str_remove(metric, "^avg_ub_"),
          metric = str_replace_all(metric, "_", " ")
        ) %>%
        pivot_wider(names_from = method, values_from = value)
      
      hi <- tbl_dat %>%
        mutate(.row = row_number()) %>%
        rowwise() %>%
        mutate(
          .winners = list({
            vals <- c_across(all_of(method_cols))
            
            # if all missing, highlight nothing
            if (all(is.na(vals))) character(0) else {
              want_max <- metric %in% max_metrics
              target <- if (want_max) max(vals, na.rm = TRUE) else min(vals, na.rm = TRUE)
              
              # return ALL tied winners
              method_cols[!is.na(vals) & vals == target]
            }
          })
        ) %>%
        ungroup() %>%
        select(.row, metric, .winners) %>%
        tidyr::unnest(.winners) %>%
        rename(method = .winners)
      
      table <- tbl_dat %>%
        gt(rowname_col = "metric") %>%
        cols_label(.list = method_labels) %>%
        fmt_number(columns = everything(), decimals = 4) %>%
        tab_header(
          title = "UB Metrics",
          subtitle = glue("Number of studies: {k} | Tau-squared: {ts} | Prior: {pr}")
        )
      
      for (k in seq_len(nrow(hi))) {
        table <- table %>%
          tab_style(
            style = cell_fill(color = "#E6F4EA"),
            locations = cells_body(
              columns = hi$method[k],
              rows    = hi$.row[k]
            )
          )
      }
      
      key <- paste0("iter_", i + 1)
      metrics_table_list[[key]] <- table
    
      i <- i + 2
    }
  }
}

############################## MSE Tables (sim_data)

mse_data <- sim_data

mse_data$MSE_ts <- ((mse_data$est_ts - mse_data$true_ts)^2)*1000

mse_data <- mse_data %>%
  group_by(method, num_studies, true_ts, prior) %>%
  summarize(avg_mse = mean(MSE_ts))

mse_table_list <- list()
i <- 1
for (ts in true_ts_values) {
  for (prior in prior_values) {
  
    # Make temp
    temp <- mse_data[mse_data$true_ts == ts, ]
    temp <- temp[temp$prior== prior | is.na(temp$prior), ]
  
    wide_temp <- temp %>%
      ungroup() %>%
      select(num_studies, method, avg_mse) %>%
      pivot_wider(names_from = method, values_from = avg_mse) %>%
      arrange(num_studies)
  
    wide_temp <- wide_temp %>%
      select(
        num_studies, RE, Weighted_1, Weighted_2, Weighted_3, Weighted_4
        #Max_Weighted_1, Max_Weighted_2, Max_Weighted_3, Max_Weighted_4
      )
  
    method_cols <- setdiff(names(wide_temp), "num_studies")
    other_cols  <- setdiff(method_cols, "RE")
    
    hi <- wide_temp %>%
      mutate(.row = row_number()) %>%
      pivot_longer(cols = all_of(method_cols), names_to = "method", values_to = "val") %>%
      group_by(.row) %>%
      mutate(re_val = val[method == "RE"][1]) %>%
      ungroup() %>%
      # highlight any non-RE method strictly less than RE
      filter(method != "RE", !is.na(val), !is.na(re_val), val < re_val) %>%
      select(.row, method)
    
    rows_where_anything_beats_RE <- unique(hi$.row)
    rows_where_RE_wins <- setdiff(seq_len(nrow(wide_temp)), rows_where_anything_beats_RE)
    hi <- bind_rows(
      hi,
      tibble(.row = rows_where_RE_wins, method = "RE")
    )
  
    table <- wide_temp %>%
      gt() %>%
      cols_label(num_studies="Studies", Weighted_1="1/k", Weighted_2="1/(k-1)",
                 Weighted_3="1/sqrt(k)", Weighted_4="1/(k^2)") %>%
      fmt_number(columns=all_of(method_cols), decimals=4) %>%
      tab_header(
        title = "MSE",
        subtitle = glue("Tau-squared: {ts} | Prior: {prior}")
      )
  
    for (k in seq_len(nrow(hi))) {
      table <- table %>%
        tab_style(
          style = cell_fill(color = "#E6F4EA"),
          locations = cells_body(
            columns = hi$method[k],
            rows    = hi$.row[k]
          )
        )
    }

    # Add to table list
    key <- paste0("iter_", i)
    mse_table_list[[key]] <- table
    i <- i + 1
  }
}

############################## Bias Tables (sim_data)

# Calculating bias
bias_data <- sim_data

bias_data <- bias_data %>%
  group_by(method, num_studies, true_ts, prior) %>%
  summarize(avg_ts_est = mean(est_ts))
  
bias_data$bias_ts <- bias_data$avg_ts_est - bias_data$true_ts

bias_table_list <- list()
i <- i
for (ts in true_ts_values) {
  for (prior in prior_values) {
    # Make temp
    temp <- bias_data[bias_data$true_ts == ts, ]
    temp <- temp[temp$prior== prior | is.na(temp$prior), ]
  
    wide_temp <- temp %>%
      ungroup() %>%
      select(num_studies, method, bias_ts) %>%
      pivot_wider(names_from = method, values_from = bias_ts) %>%
      arrange(num_studies)
  
    wide_temp <- wide_temp %>%
      select(
        num_studies, RE, Weighted_1, Weighted_2, Weighted_3, Weighted_4
        #Max_Weighted_1, Max_Weighted_2, Max_Weighted_3, Max_Weighted_4
      )

    method_cols <- setdiff(names(wide_temp), "num_studies")
    other_cols  <- setdiff(method_cols, "RE")
    
    hi <- wide_temp %>%
      mutate(.row = row_number()) %>%
      pivot_longer(cols = all_of(method_cols), names_to = "method", values_to = "val") %>%
      group_by(.row) %>%
      mutate(re_val = val[method == "RE"][1]) %>%
      ungroup() %>%
      # highlight any non-RE method strictly less than RE
      filter(method != "RE", !is.na(val), !is.na(re_val), abs(val) < abs(re_val)) %>%
      select(.row, method)
    
    rows_where_anything_beats_RE <- unique(hi$.row)
    rows_where_RE_wins <- setdiff(seq_len(nrow(wide_temp)), rows_where_anything_beats_RE)
    hi <- bind_rows(
      hi,
      tibble(.row = rows_where_RE_wins, method = "RE")
    )
  
    table <- wide_temp %>%
      gt() %>%
      cols_label(num_studies="Studies", Weighted_1="1/k", Weighted_2="1/(k-1)",
                 Weighted_3="1/sqrt(k)", Weighted_4="1/(k^2)") %>%
      fmt_number(columns=all_of(method_cols), decimals=4) %>%
      tab_header(
        title = "Bias",
        subtitle = glue("Tau-squared: {ts} | Prior: {prior}")
      )
  
    for (k in seq_len(nrow(hi))) {
      table <- table %>%
        tab_style(
          style = cell_fill(color = "#E6F4EA"),
          locations = cells_body(
            columns = hi$method[k],
            rows    = hi$.row[k]
          )
        )
    }
  
    # Add to table list
    key <- paste0("iter_", i)
    bias_table_list[[key]] <- table
    i <- i + 1
  
  }
}

######################################################################
# View and save tables
######################################################################
# 84 tables
#lapply(metrics_table_list, print)

# 21 tables
#lapply(mse_table_list, print)

# 21 tables
#lapply(bias_table_list, print)

#for (i in seq_along(mse_table_list)) {
#  gtsave(
#    data = mse_table_list[[i]],
#    filename = file.path("Tables", paste0("table_", i, ".png"))
#  )
#}

#for (i in seq_along(bias_table_list)) {
#  gtsave(
#    data = bias_table_list[[i]],
#    filename = file.path("Tables", paste0("table_", i, ".png"))
#  )
#}

#for (i in seq_along(metrics_table_list)) {
#  gtsave(
#    data = metrics_table_list[[i]],
#    filename = file.path("Tables", paste0("table_", i, ".png"))
#  )
#}


  
  
  
  
  
  
  