library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)

######################################################################
# README
######################################################################

# This file is an example of a shiny app used to display data from 
# simulation results. One of the papers in my disseratation involves
# a simulation study, comparing the performance of our proposed method
# to previous, common practice estimators. This file shows both data
# manipulation and visualization, as well as overall coding proficiency.
# I chose to include this file rather than the simulation file itself,
# as the simulation file is mostly just modeling fitting over a grid of
# different parameters, but I would be happy to share any other syntax files.

######################################################################
# App
######################################################################
sim_data <- read.csv("Data/apr8_sim.csv")
mse_data <- sim_data
mse_data$MSE_ts <- ((mse_data$est_ts - mse_data$true_ts)^2)*1000
mse_data <- mse_data %>%
  group_by(method, num_studies, true_ts, prior) %>%
  summarize(avg_mse = mean(MSE_ts))

bias_data <- sim_data
bias_data <- bias_data %>%
  group_by(method, num_studies, true_ts, prior) %>%
  summarize(avg_ts_est = mean(est_ts))
bias_data$bias_ts <- bias_data$avg_ts_est - bias_data$true_ts

ui <- fluidPage(
  titlePanel("MSE by Number of Studies"),
  sidebarLayout(
    sidebarPanel(
      selectInput("ts", "Tau-squared", choices = NULL, selected = NULL),
      selectInput("prior", "Prior", choices = NULL, selected = NULL),
      selectInput("kmax", "Max k (show up to this many studies)", choices = NULL),
      checkboxInput("log_y", "Log scale y-axis", FALSE),
      checkboxInput("exclude_ts0", "Exclude tau-squared = 0 (condition on heterogeneity)", TRUE),
      radioButtons(
        "outcome", "Plot",
        choices = c("MSE" = "mse", "Bias" = "bias"),
        selected = "mse",
        inline = TRUE
      ),
      radioButtons(
        "k_weight",
        "Weight k values",
        choices = c(
          "Uniform over k" = "uniform",
          "Emphasize small k (1/k)" = "invk"
        ),
        selected = "invk",
        inline = TRUE
      ),
      helpText("Hover points for values. Pan/zoom with mouse; double-click to reset.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Methods", plotlyOutput("mse_plot", height = "600px")),
        tabPanel("Recommendation",
                 selectInput("ts_reco", "Tau-squared used for recommendation", choices = NULL),
                 verbatimTextOutput("reco_text"),
                 plotlyOutput("reco_plot", height = "450px")
        ),
        tabPanel("Safe vs RE",
                 numericInput("eps", "Safety tolerance ε (treat losses smaller than ε as ties)", value = 0, min = 0, step = 0.0001),
                 checkboxInput("safe_only", "Only show combos that never lose to RE", TRUE),
                 plotlyOutput("safety_plot", height = "450px"),
                 plotlyOutput("safe_rank_plot", height = "450px")
        ),
        tabPanel("Winner map",
                 checkboxInput("facet_prior", "Show all priors (facet)", FALSE),
                 plotlyOutput("winner_map", height = "650px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  method_labels <- c(
    RE         = "RE",
    Weighted_1 = "1/k",
    Weighted_2 = "1/(k-1)",
    Weighted_3 = "1/sqrt(k)",
    Weighted_4 = "1/(k^2)"
  )
  pretty_levels <- unname(method_labels)
  weighted_methods <- c("Weighted_1", "Weighted_2", "Weighted_3", "Weighted_4")
  
  # Populate dropdowns once
  observe({
    ts_vals <- sort(unique(mse_data$true_ts))
    pr_vals <- sort(unique(mse_data$prior))
    k_vals  <- sort(unique(mse_data$num_studies))
    updateSelectInput(session, "ts", choices = ts_vals, selected = ts_vals[1])
    updateSelectInput(session, "prior", choices = pr_vals, selected = pr_vals[1])
    updateSelectInput(session, "kmax",  choices = k_vals,  selected = max(k_vals))
    updateSelectInput(session, "ts_reco",choices = ts_vals, selected = ts_vals[1])
  })
  
  filtered_dat <- reactive({
    req(input$ts, input$prior, input$kmax, input$outcome)
    
    base_df <- if (input$outcome == "mse") mse_data else bias_data
    
    # Choose y column name + labels
    y_col <- if (input$outcome == "mse") "avg_mse" else "bias_ts"
    y_lab <- if (input$outcome == "mse") "MSE" else "Bias"
    
    df <- base_df %>%
      mutate(
        num_studies = as.numeric(num_studies),
        method      = trimws(as.character(method))
      ) %>%
      filter(true_ts == input$ts) %>%
      filter(prior == input$prior | is.na(prior)) %>%
      filter(num_studies <= as.numeric(input$kmax)) %>%
      group_by(true_ts, prior, num_studies, method) %>%
      summarize(y = mean(.data[[y_col]]), .groups = "drop") %>%
      mutate(
        method_pretty = recode(method, !!!method_labels),
        method_pretty = factor(method_pretty, levels = pretty_levels),
        method = factor(method, levels = names(method_labels)),
        y_label = y_lab
      ) %>%
      arrange(method_pretty, num_studies)
    
    df
  })
  
  winner_dat <- reactive({
    req(input$outcome)
    
    weighted_methods <- c("Weighted_1","Weighted_2","Weighted_3","Weighted_4")
    
    base_df <- if (input$outcome == "mse") mse_data else bias_data
    loss_col <- if (input$outcome == "mse") "avg_mse" else "bias_ts"
    
    df <- base_df %>%
      mutate(
        num_studies = as.numeric(num_studies),
        method = trimws(as.character(method))
      ) %>%
      filter(method %in% weighted_methods)
    
    # Prior handling: either facet all priors, or use selected prior (+ NA priors if they exist)
    if (isTRUE(input$facet_prior)) {
      # keep all priors
    } else {
      req(input$prior)
      df <- df %>% filter(prior == input$prior | is.na(prior))
    }
    
    df <- df %>%
      group_by(true_ts, prior, num_studies, method) %>%
      summarize(val = mean(.data[[loss_col]]), .groups = "drop") %>%
      mutate(loss = if (input$outcome == "mse") val else abs(val))
    
    df <- df %>% filter(num_studies <= as.numeric(input$kmax))
    
    # Winner per (tau2, prior, k)
    winners <- df %>%
      group_by(true_ts, prior, num_studies) %>%
      slice_min(loss, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(
        method_pretty = recode(method, !!!method_labels),
        method_pretty = factor(method_pretty, levels = pretty_levels),
        true_ts = factor(true_ts) # discrete tiles look cleaner
      )
    
    winners
  })
  
  reco <- reactive({
    req(input$k_weight)
    
    # Use the pre-aggregated MSE data (avg_mse), and compare ONLY weighted methods
    wdat <- mse_data %>%
      mutate(
        num_studies = as.numeric(num_studies),
        avg_mse     = as.numeric(avg_mse),
        method      = trimws(as.character(method))
      ) %>%
      filter(method %in% weighted_methods)
    
    # weights over k: uniform vs 1/k
    wdat <- wdat %>%
      mutate(k_wt = if (input$k_weight == "invk") 1/num_studies else 1)
    
    if (isTRUE(input$exclude_ts0)) {
      wdat <- wdat %>% filter(true_ts > 0)
    }
    
    # ORACLE best at each (tau2, k) over (method, prior)
    oracle <- wdat %>%
      group_by(true_ts, num_studies) %>%
      summarize(oracle_mse = min(avg_mse, na.rm = TRUE), .groups = "drop")
    
    # regret and average over k for each (method, prior, tau2)
    reg_tau <- wdat %>%
      left_join(oracle, by = c("true_ts","num_studies")) %>%
      mutate(regret = avg_mse - oracle_mse) %>%
      group_by(method, prior, true_ts) %>%
      summarize(
        mean_regret = weighted.mean(regret, w = k_wt, na.rm = TRUE),
        mean_mse    = weighted.mean(avg_mse, w = k_wt, na.rm = TRUE),
        .groups = "drop"
      )
    
    # minimax over tau2 for each (method, prior)
    reg_worst <- reg_tau %>%
      group_by(method, prior) %>%
      summarize(
        worst_regret = max(mean_regret, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(worst_regret)
    
    list(reg_tau = reg_tau, reg_worst = reg_worst)
  })
  
  reco_by_ts <- reactive({
    req(input$ts_reco, input$k_weight)
    
    weighted_methods <- c("Weighted_1","Weighted_2","Weighted_3","Weighted_4")
    
    wdat <- mse_data %>%
      mutate(
        num_studies = as.numeric(num_studies),
        avg_mse     = as.numeric(avg_mse),
        method      = trimws(as.character(method))
      ) %>%
      filter(method %in% weighted_methods) %>%
      filter(true_ts == input$ts_reco) %>%
      mutate(k_wt = if (input$k_weight == "invk") 1/num_studies else 1)
    
    # Expected MSE across k for each (method, prior) at this tau^2
    score <- wdat %>%
      group_by(method, prior) %>%
      summarize(
        exp_mse = weighted.mean(avg_mse, w = k_wt, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(exp_mse)
    
    score
  })
  
  safe_vs_re <- reactive({
    req(input$eps, input$k_weight)
    
    eps <- as.numeric(input$eps)
    
    weighted_methods <- c("Weighted_1","Weighted_2","Weighted_3","Weighted_4")
    
    # RE benchmark by (tau2, k) — RE may have prior = NA, so ignore prior here
    re <- mse_data %>%
      mutate(num_studies = as.numeric(num_studies),
             avg_mse = as.numeric(avg_mse),
             method = trimws(as.character(method))) %>%
      filter(method == "RE") %>%
      group_by(true_ts, num_studies) %>%
      summarize(re_mse = mean(avg_mse, na.rm = TRUE), .groups = "drop")
    
    if (isTRUE(input$exclude_ts0)) {
      re <- re %>% filter(true_ts > 0)
    }
    
    # Weighted combos by (tau2, prior, k)
    w <- mse_data %>%
      mutate(num_studies = as.numeric(num_studies),
             avg_mse = as.numeric(avg_mse),
             method = trimws(as.character(method))) %>%
      filter(method %in% weighted_methods) %>%
      filter(num_studies <= as.numeric(input$kmax)) %>%
      group_by(true_ts, prior, num_studies, method) %>%
      summarize(combo_mse = mean(avg_mse, na.rm = TRUE), .groups = "drop") %>%
      inner_join(re, by = c("true_ts","num_studies")) %>%
      mutate(
        delta = re_mse - combo_mse,   # positive means combo beats RE
        k_wt  = if (input$k_weight == "invk") 1/num_studies else 1
      )
    
    if (isTRUE(input$exclude_ts0)) {
      w <- w %>% filter(true_ts > 0)
    }
    
    # Safety summary: "never lose to RE" => min(delta) >= -eps
    safety <- w %>%
      group_by(method, prior) %>%
      summarize(
        min_delta  = min(delta, na.rm = TRUE),
        mean_delta = mean(delta, na.rm = TRUE),
        frac_win   = mean(delta >  eps, na.rm = TRUE),
        frac_tie   = mean(abs(delta) <= eps, na.rm = TRUE),
        frac_loss  = mean(delta < -eps, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        safe = (min_delta >= -eps),
        method_pretty = recode(method, !!!method_labels),
        combo = paste0(method_pretty, " | prior=", prior)
      )
    
    # Now: among SAFE combos, compute expected MSE across k and tau2
    # (1) average over k within each tau2
    by_tau <- w %>%
      group_by(method, prior, true_ts) %>%
      summarize(
        exp_mse_over_k = weighted.mean(combo_mse, w = k_wt, na.rm = TRUE),
        .groups = "drop"
      )
    
    # (2) average across tau2 grid (uniform across your simulated tau2 values)
    overall <- by_tau %>%
      group_by(method, prior) %>%
      summarize(
        exp_mse = mean(exp_mse_over_k, na.rm = TRUE),
        worst_tau_mse = max(exp_mse_over_k, na.rm = TRUE),  # optional robustness
        .groups = "drop"
      ) %>%
      left_join(safety %>% select(method, prior, safe, min_delta, mean_delta, frac_win, combo), by = c("method","prior")) %>%
      mutate(method_pretty = recode(method, !!!method_labels)) %>%
      arrange(desc(safe), exp_mse)
    
    list(w = w, safety = safety, overall = overall)
  })
  
  weighted_methods <- c("Weighted_1","Weighted_2","Weighted_3","Weighted_4")
  comp_dat <- reactive({
    df <- filtered_dat() %>%
      mutate(method = as.character(method)) %>%
      group_by(num_studies, method) %>%
      summarize(avg_mse = mean(avg_mse), .groups = "drop")
    
    # Get RE per k
    re <- df %>%
      filter(method == "RE") %>%
      select(num_studies, re_mse = avg_mse)
    
    # Get best weighted per k
    bestw <- df %>%
      filter(method %in% weighted_methods) %>%
      group_by(num_studies) %>%
      slice_min(avg_mse, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        num_studies,
        best_w_method = method,
        best_w_mse = avg_mse
      )
    
    # Join and compute comparison metrics
    comp <- re %>%
      inner_join(bestw, by = "num_studies") %>%
      mutate(
        diff = re_mse - best_w_mse,                         # + means weighted better
        ratio = best_w_mse / re_mse,                        # <1 means weighted better
        pct_improve = 100 * (re_mse - best_w_mse) / re_mse  # + means weighted better
      )
    
    list(df = df, comp = comp)
  })
  
  output$mse_plot <- renderPlotly({
    df <- filtered_dat()
    
    # guardrail: if a method has only 1 point, line won't appear (still shows the point)
    g <- ggplot(df, aes(x = num_studies, y = y, color = method_pretty, group = method_pretty)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1) +
      scale_color_manual(
        values = c(
          "RE"        = "#1b9e77",
          "1/k"       = "#d95f02",
          "1/(k-1)"   = "#7570b3",
          "1/sqrt(k)" = "#e7298a",
          "1/(k^2)"   = "#66a61e"
        )
      ) +
      theme_minimal(base_size = 14) +     
      labs(
        title = "MSE by Number of Studies",
        subtitle = paste0("Tau-squared: ", input$ts, " | Prior: ", input$prior),
        x = "Number of studies",
        y = unique(df$y_label),
        color = "Method"
      )
    
    if (isTRUE(input$log_y)) {
      g <- g + scale_y_log10()
    }
    
    title_txt <- if (input$outcome == "mse") "MSE by Number of Studies" else "Bias by Number of Studies"
  
    ggplotly(g, tooltip = c("x", "y", "colour")) %>%
      layout(
        title = list(
          text = paste0(
            title_txt,
            "<br><sup>Tau-squared: ", input$ts, " | Prior: ", input$prior, "</sup>"
          )
        ),
        legend = list(orientation = "v")
      )

  })
  
  output$reco_text <- renderText({
    best <- reco_by_ts() %>% slice(1)
    best_method_pretty <- method_labels[[best$method]]
    
    paste0(
      "Recommendation for tau-squared = ", input$ts_reco, "\n",
      "Criterion: minimum expected MSE over k\n",
      "k weighting: ", if (input$k_weight == "invk") "1/k (emphasize small k)" else "uniform", "\n\n",
      "Recommended weights: ", best_method_pretty, "\n",
      "Recommended prior: ", best$prior, "\n",
      "Expected MSE (over k): ", round(best$exp_mse, 6)
    )
  })
  
  output$reco_plot <- renderPlotly({
    df <- reco_by_ts() %>%
      mutate(
        method_pretty = recode(method, !!!method_labels),
        combo = paste0(method_pretty, " | prior=", prior)
      )
    
    g <- ggplot(df, aes(x = reorder(combo, exp_mse), y = exp_mse)) +
      geom_col() +
      coord_flip() +
      labs(
        title = paste0("Best (weights, prior) at tau-squared = ", input$ts_reco),
        subtitle = paste0("Expected MSE over k (", if (input$k_weight == "invk") "1/k weighting" else "uniform weighting", ")"),
        x = NULL,
        y = "Expected MSE over k"
      ) +
      theme_minimal(base_size = 13)
    
    ggplotly(g, tooltip = c("x","y"))
  })
  
  output$winner_map <- renderPlotly({
    w <- winner_dat()
    
    title_txt <- if (input$outcome == "mse") "Winner map (min MSE)" else "Winner map (min |Bias|)"
    
    g <- ggplot(w, aes(x = true_ts, y = num_studies, fill = method_pretty,
                       text = paste0("tau2: ", true_ts,
                                     "<br>k: ", num_studies,
                                     "<br>winner: ", method_pretty))) +
      geom_tile() +
      labs(
        title = title_txt,
        subtitle = if (isTRUE(input$facet_prior)) "Faceted by prior" else paste0("Prior: ", input$prior),
        x = "Tau-Squared",
        y = "Number of studies (k)",
        fill = "Winner"
      ) +
      theme_minimal(base_size = 13)
    
    if (isTRUE(input$facet_prior)) {
      g <- g + facet_wrap(~prior)
    }
    
    ggplotly(g, tooltip = "text")
  })
  
  output$safety_plot <- renderPlotly({
    df <- safe_vs_re()$safety %>%
      arrange(desc(safe), min_delta)
    
    if (isTRUE(input$safe_only)) df <- df %>% filter(safe)
    
    g <- ggplot(df, aes(x = reorder(combo, min_delta), y = min_delta)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
      geom_col() +
      coord_flip() +
      labs(
        title = "Safety vs RE (worst-case across τ² and k)",
        subtitle = paste0(
          "min(RE MSE − combo MSE); ≥ 0 means never worse than RE",
          if (isTRUE(input$exclude_ts0)) " (computed over τ²>0 only)" else ""
        ),
        x = NULL,
        y = "Worst-case improvement margin"
      ) +
      theme_minimal(base_size = 13)
    
    ggplotly(g, tooltip = c("x","y"))
  })
  
  output$safe_rank_plot <- renderPlotly({
    df <- safe_vs_re()$overall
    
    if (isTRUE(input$safe_only)) df <- df %>% filter(safe)
    
    g <- ggplot(df, aes(x = reorder(paste0(method_pretty, " | prior=", prior), exp_mse), y = exp_mse)) +
      geom_col() +
      coord_flip() +
      labs(
        title = "Rank within safe set",
        subtitle = paste0("Expected MSE (avg over τ²; k-weighting: ", if (input$k_weight=="invk") "1/k" else "uniform", ")"),
        x = NULL,
        y = "Expected MSE"
      ) +
      theme_minimal(base_size = 13)
    
    ggplotly(g, tooltip = c("x","y"))
  })
  
}

shinyApp(ui, server)
