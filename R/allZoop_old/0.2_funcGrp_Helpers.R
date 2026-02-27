# Helper functions for temperature sensitivity data analysis for funcGrps
# 09/02/2026
# Josh Hill



# Packages ----
library(tidyverse)
library(glmmTMB)



############# FUNCTIONS FOR DATA ANALYSIS ############# 



# BOOTSTRAP MODEL FUNCTION ----
# Using glmmTMB...
# Calculate Q10 with CI and prediction ribbons for plotting using bootstrapping for confidence intervals
boot_Q10 <- function(df) {
  require(glmmTMB, quietly = TRUE)
  
  df_id <- sample(1:nrow(df), nrow(df), replace = TRUE) # A random sample of rows of df (with replacement) of same size as df
  d <- df[df_id,] # Get those random rows
  out <- glmmTMB(y ~ x * funcGrp + (1|primRef) + (1|taxa),  data = d) # fit the model
  
  return(out)
}
# END OF BOOTSTRAP MODEL FUNCTION



# GET Q10 FUNCTION ----
# Get Q10s from bootstrap models
get_Q10s <- function(m) {
  require(dplyr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  
  # Define groups and only progress for 1 group at a time
  data <- m$frame %>% 
    group_by(funcGrp) %>% 
    group_split()
  
  Q10_by_group <- function(d) {
    grp <- unique(d$funcGrp)
    stopifnot(length(grp) == 1)
    # Define temperatures
    d <- d %>%  
      mutate(temp_col = 1/x)
    T1 <- min(d$temp_col, na.rm = TRUE) # Take the minimum temperature per group
    T2 <- T1 + 10 # Take minimum temp and add specified temperature
    
    # Create prediction data for just two temps
    newdat <- data.frame(x = c(1 / T1, 1 / T2), # A dataframe with JUST the predictor comprising the reciprocals of two temps
                         funcGrp = grp) 
    
    # Get predictions with std err
    fit <- predict(m, newdata = newdat, # Predict the responses (in ln space) for the two temps
                   se.fit = FALSE, 
                   re.form = NA) # Maintain population-level predictions (ie., ignore random effects)
    # Calculate Q10 value
    Q10 <- exp(fit[2] - fit[1]) # The ratio of rates at the warmer vs cooler temps
    return(tibble(Group = grp, Q10 = Q10))
  }
  map_dfr(data, Q10_by_group)
}
# END OF GETQ10 FUNCTION



# EXTRACTING RESULTS FUNCTION
# Extract predictions and Q10 values for each group to prepare for Arrhenius plots
get_results <- function(df, groups, coefs_list, Q10pdat) {
  
  results <- map(groups, function(g) {
    # Extract x values for the current zooplankton group
    x_vals <- df %>% filter(funcGrp == g) %>% pull(x)
    
    # Skip if no data for a group
    if (length(x_vals) == 0 || all(is.na(x_vals))) return(NULL)
    
    # Create a sequence of 100 evenly-spaced x values spanning the data range to generate smooth pred curves
    x_seq <- seq(min(x_vals, na.rm = TRUE), 
                 max(x_vals, na.rm = TRUE), 
                 length.out = 100)
    
    # Generate predictions for each bootstrap iteration
    preds <- imap_dfr(coefs_list, function(fx, i) {
      # Calculate group-specific intercept (base intercept + group adjustment)
      intercept <- fx["(Intercept)"] + ifelse(paste0("funcGrp", g) %in% names(fx), fx[paste0("funcGrp", g)], 0)
      
      # Calculate group-specific slope (base slope + group interaction effect)
      slope <- fx["x"] + ifelse(paste0("x:funcGrp", g) %in% names(fx), fx[paste0("x:funcGrp", g)], 0)
      
      # Generate predicted values using the linear model... back transform from ln space using exp(intercept + slope * x)
      tibble(
        x = x_seq,
        predicted = exp(intercept + slope * x_seq),
        boot_id = i  # Track the bootstrap iteration
      )
    })
    
    # Summarise predictions across all bootstrap iterations
    # Give me the median prediction and confidence intervals
    pred_summary <- preds %>%
      filter(!is.na(predicted)) %>%
      group_by(x) %>%
      group_split() %>%  # Split into separate dataframes by x
      map_dfr(function(group_data) {
        tibble(
          x = unique(group_data$x),
          predicted = median(group_data$predicted), # median prediction
          conf.low = quantile(group_data$predicted, 0.025),
          conf.high = quantile(group_data$predicted, 0.975)
        )
      })
    
    # Extract the Q10 values
    Q10_row <- Q10pdat %>% filter(Group == g)
    Q10_list <- list(Q10 = Q10_row$median, 
                     CI_lower = Q10_row$lower_CI, 
                     CI_upper = Q10_row$upper_CI)
    
    # Return both the prediction curve and Q10 statistics for this group
    list(predictions = pred_summary, Q10 = Q10_list)
  })
  
  # Name the results list by group
  names(results) <- groups
  return(results)
}
# END OF EXTRACTING RESULTS FUNCTION



######################### FUNCTIONS FOR PLOTTING #########################


# ARRHENIUS PLOT FUNCTION ----
arrhenius_plot <- function(mdat, rate_col, results, 
                           group_order = NULL, x_limits = NULL, ncol = 3) {
  
  # Prepare data for plotting
  plot_data <- map(names(results), function(g) {
    if (is.null(results[[g]])) return(NULL)
    
    df_data <- mdat %>% filter(funcGrp == g)
    df_pred <- results[[g]]$predictions
    df_Q10 <- results[[g]]$Q10
    
    #Calculate sample size and number of taxa
    n_obs <- nrow(df_data)
    
    # Create the label for annotation
    list(
      data = df_data %>% mutate(funcGrp = g),
      pred = df_pred %>% mutate(funcGrp = g),
      Q10_label = paste0(
        "Q[10] == ", round(df_Q10$Q10, 2),
        " * ' (95% CI: ", round(df_Q10$CI_lower, 2), " - ",
        round(df_Q10$CI_upper, 2), "); n = ", n_obs, "'")
      )
  })
  
  # Remove NULLs
  plot_data <- compact(plot_data)
  
  # Combine all data into single dataframes
  all_data <- bind_rows(map(plot_data, \(x) x$data))
  all_pred <- bind_rows(map(plot_data, \(x) x$pred))
  all_labels <- tibble(
    funcGrp = map_chr(plot_data, \(x) as.character(unique(x$data$funcGrp))),
    label = map_chr(plot_data, \(x) x$Q10_label),
    x = 0.0037,
    y = Inf
  )
  
  # Re-order funcGrps if specified
  if (!is.null(group_order)) {
    all_data <- all_data %>%
      mutate(funcGrp = fct_relevel(funcGrp, group_order))
    
    all_pred <- all_pred %>%
      mutate(funcGrp = fct_relevel(funcGrp, group_order))
    
    all_labels <- all_labels %>%
      mutate(funcGrp = fct_relevel(funcGrp, group_order))
  }
  
  # Set x-axis limits... unless specified
  if (is.null(x_limits)) {
    x_limits <- c(max(all_data$x, na.rm = TRUE), 
                  min(all_data$x, na.rm = TRUE))
  }
  
  # Create the plot
  p <- ggplot() +
    # Add the raw data
    geom_point(data = all_data, aes(x = x, y = .data[[rate_col]]), 
               alpha = 0.3, size = 2) +
    # Add the confidence ribbons
    geom_ribbon(data = all_pred, aes(x = x, ymin = conf.low, ymax = conf.high),
                fill = "grey", alpha = 0.3) +
    # Add the fit
    geom_line(data = all_pred, aes(x = x, y = predicted), 
              colour = "darkblue", linewidth = 1, linetype = "dashed") +
    # Add the annotation
    geom_text(data = all_labels, aes(x = x, y = y, label = label), 
              parse = TRUE, hjust = 0, vjust = 1.5, size = 4, colour = "darkblue") +
    scale_x_continuous(name = expression(bold("1 / Temp (K"^-1*")")),
                       trans = "reverse",
                       sec.axis = sec_axis(~1/. - 273.15, name = "Temp (°C)")) +
    facet_wrap(~ funcGrp, scales = "free_y", ncol = ncol) +
    coord_cartesian(xlim = x_limits) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.title.x.top = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 12, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(t = 5, r = 15, b = 5, l = 5))
  
  return(p)
}
# END OF ARRHENIUS PLOT FUNCTION



# Q10 PLOT FUNCTION ----
plotQ10 <- function(Q10pdat, data_type, colours){
  
  # Get the rate-specific colour
  rateCol <- colours [[data_type]]

  p <- ggplot() +
    geom_errorbar(data = Q10pdat, aes(x = Group, ymin = lower_CI, ymax = upper_CI),
                  width = .2,
                  linewidth = 1,
                  colour = rateCol) +
    geom_point(data = Q10pdat, aes(x = Group, y = median),
               size = 3,
               colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.title.x.top = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12)) +  
    labs(
      x = expression(bold("Zooplankton group")),
      y = bquote(bold(.(data_type)~"rate Q"[10])))
  
  return(p)
}
# END OF Q10 PLOT FUNCTION


