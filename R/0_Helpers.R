# Helper functions for temperature sensitivity data analysis
# 16/09/25
# Josh Hill



# Packages ----
library(tidyverse)
library(glmmTMB)



############# Functions ############# 

# Carbon weight converter ----
# Convert carbon weights to mg
convert_CW <- function(weight, unit) {
  if(is.na(weight) || is.na(unit)) return(NA_real_)
  if(unit == "mg") return(weight)              # maintain mg
  if(unit == "ug") return(weight / 1000)       # µg to mg
  if(unit == "ng") return(weight / 1e6)        # ng to mg
  if(unit == "g")  return(weight * 1000)       # g to mg
  if(unit == "kg") return(weight * 1e6)        # kg to mg
  
  return(NA_real_)
}
# END OF CARBON WEIGHT CONVERTER



# CLEARANCE CONVERTER ----
# Convert clearance rates to ml/ind/hr
convert_clearance <- function(rate, unit) {
  if (is.na(rate) || is.na(unit)) return(list(rate = NA_real_, unit = NA_character_))
  
  if (unit == "ml/ind/hr")  return(list(rate = rate,                unit = "ml/ind/hr")) # maintain ml/ind/hr
  if (unit == "ml/ind/day") return(list(rate = rate / 24,           unit = "ml/ind/hr")) # day to hr
  if (unit == "l/ind/day")  return(list(rate = (rate * 1000) / 24,  unit = "ml/ind/hr")) # l to ml, day to hr
  if (unit == "l/ind/hr")   return(list(rate = rate * 1000,         unit = "ml/ind/hr")) # l to ml
  if (unit == "ml/mgC/day") return(list(rate = rate / 24,           unit = "ml/mgC/hr")) # day to hr
  
  return(list(rate = NA_real_, unit = NA_character_))
}
# END OF CLEARANCE CONVERTER



# INGESTION CONVERTER ----
# Convert ingestion rates to mgC/ind/hr
convert_ingestion <- function(rate, unit) {
  if (is.na(rate) || is.na(unit)) return(list(rate = NA_real_, unit = NA_character_))
  
  if (unit == "mgC/ind/hr") return(list(rate = rate,                unit = "mgC/ind/hr")) # maintain mgC/ind/hr
  if (unit == "ugC/ind/hr") return(list(rate = rate / 1000,         unit = "mgC/ind/hr")) # µg to mg
  if (unit == "ugC/ind/day") return(list(rate = (rate / 1000) / 24, unit = "mgC/ind/hr")) # µg to mg, day to hr
  if (unit == "ngC/ind/hr") return(list(rate = rate / 1e6,          unit = "mgC/ind/hr")) # ng to mg
  if (unit == "ngC/ind/day") return(list(rate = (rate / 1e6) / 24,  unit = "mgC/ind/hr")) # ng to mg, day to hr
  if (unit == "ugC/mgC/hr") return(list(rate = rate / 1000,         unit = "mgC/mgC/hr")) # ngC to mgC
  if (unit == "ugC/ugC/day") return(list(rate = rate / 24,          unit = "mgC/mgC/hr")) # mgC to mgC ratio cancels out, day to hr
  
  return(list(rate = NA_real_, unit = NA_character_))
}
# END OF INGESTION CONVERTER



# RESPIRATION CONVERTER ----
# Convert respiration rates to ulO2/ind/hr
convert_respiration <- function(rate, unit, genus) {
  if (is.na(rate) || is.na(unit)) return(list(rate = NA_real_, unit = NA_character_))
  
  if (unit == "ulO2/ind/hr") return(list(rate = rate,                  unit = "ulO2/ind/hr")) # maintain ulO2/ind/hr
  if (unit == "ulO2/ind/day") return(list(rate = rate / 24,            unit = "ulO2/ind/hr")) # day to hr
  if (unit == "ulO2/mgC/hr") return(list(rate = rate,                  unit = "ulO2/mgC/hr")) # maintain ulO2/mgC/hr
  if (unit == "ulO2/ugC/day") return(list(rate = (rate / 1000) / 24,   unit = "ulO2/mgC/hr")) # ugC to mgC
  if (unit == "umol/ind/hr") return(list(rate = rate * 22.4,           unit = "ulO2/ind/hr")) # umol to ulO2, as per the ideal gas law (i.e., 1 mol of gas = 22.4 L)
  if (unit == "nlO2/ind/hr") return(list(rate = rate / 1000,           unit = "ulO2/ind/hr")) # nlO2 to ulO2
  
  # Stoichiometry
  if (genus == "Euphausia" & unit == "ugC/ind/day") 
    return(list(rate = ((rate / 24) / (12.011 * 0.9)) * 22.4,         unit = "ulO2/ind/hr")) # ugC to ulO2, where 0.9 = RQ (Ross1982), day to hr
  
  return(list(rate = NA_real_, unit = NA_character_))
}
# END OF RESPIRATION CONVERTER



# TAXONOMIC COVERAGE SUMMARY FUNCTION ----
summarise_taxonomic_coverage <- function(data, dataset_name) {
  coverage <- data %>%
    summarise(
      Dataset = dataset_name,
      n_species = n_distinct(taxa), # change taxa to species
      n_genera = n_distinct(genus),
      n_families = n_distinct(family, na.rm = TRUE),
      n_orders = n_distinct(order, na.rm = TRUE),
      n_classes = n_distinct(class, na.rm = TRUE),
      n_phyla = n_distinct(phylum, na.rm = TRUE),
      n_observations = n()
    )
  return(coverage)
}
# END OF TAXONOMIC SUMMARY FUNCTION



# THERMAL SENSITIVITY FUNCTIONS ----
# Two functions that are used to estimate thermal sensitivities per group

# Calculate Q10 with CI and prediction ribbons for plotting using bootstrapping for confidence intervals
boot_Q10 <- function(df) {
  require(glmmTMB, quietly = TRUE)
  
  df_id <- sample(1:nrow(df), nrow(df), replace = TRUE) # A random sample of rows of df (with replacement) of same size as df
  d <- df[df_id,] # Get those random rows
  out <- glmmTMB(y ~ x * zoopGrp + (1|primRef) + (1|taxa),  data = d) # fit the model
}

# Get Q10s from bootstrap models
get_Q10s <- function(m) {
  require(dplyr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  
  # Define groups and only progress for 1 group at a time
  data <- m$frame %>% 
    group_by(zoopGrp) %>% 
    group_split()
  
  Q10_by_group <- function(d) {
    grp <- unique(d$zoopGrp)
    stopifnot(length(grp) == 1)
    # Define temperatures
    d <- d %>%  
      mutate(temp_col = 1/x)
    T1 <- min(d$temp_col, na.rm = TRUE) # Take the minimum temperature per group
    T2 <- T1 + 10 # Take minimum temp and add specified temperature
    
    # Create prediction data for just two temps
    newdat <- data.frame(x = c(1 / T1, 1 / T2), # A dataframe with JUST the predictor comprising the reciprocals of two temps
                         zoopGrp = grp) 
    
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
# END OF THERMAL SENSITIVITY ANALYSERS



# EXTRACTING RESULTS FUNCTION
# Extract predictions and Q10 values for each group to prepare for Arrhenius plots
get_results <- function(df, groups, coefs_list, Q10pdat) {
  
  results <- map(groups, function(g) {
    # Extract x values for the current zooplankton group
    x_vals <- df %>% filter(zoopGrp == g) %>% pull(x)
    
    # Skip if no data for a group
    if (length(x_vals) == 0 || all(is.na(x_vals))) return(NULL)
    
    # Create a sequence of 100 evenly-spaced x values spanning the data range to generate smooth pred curves
    x_seq <- seq(min(x_vals, na.rm = TRUE), 
                 max(x_vals, na.rm = TRUE), 
                 length.out = 100)
    
    # Generate predictions for each bootstrap iteration
    preds <- imap_dfr(coefs_list, function(fx, i) {
      # Calculate group-specific intercept (base intercept + group adjustment)
      intercept <- fx["(Intercept)"] + ifelse(paste0("zoopGrp", g) %in% names(fx), fx[paste0("zoopGrp", g)], 0)
      
      # Calculate group-specific slope (base slope + group interaction effect)
      slope <- fx["x"] + ifelse(paste0("x:zoopGrp", g) %in% names(fx), fx[paste0("x:zoopGrp", g)], 0)
      
      # Generate predicted values using the linear model... back transform from ln space using exp(intercept + slope * x)
      tibble(
        x = x_seq,
        predicted = exp(intercept + slope * x_seq),
        boot_id = i  # Track the bootstrap iteration
      )
    })
    
    # Summarise predictions across all bootstrap iterations
    # Gives us the median prediction and confidence intervals
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
# END OF EXTRACTING RESULTS Fn



# ARRHENIUS PLOT FUNCTION
arrhenius_plot <- function(mdat, rate_col, boot_models, Q10pdat, 
                           group_order = NULL, x_limits = NULL, ncol = 3) {
  
  # Prepare data for plotting
  plot_data <- map(names(results), function(g) {
    if (is.null(results[[g]])) return(NULL)
    
    df_data <- mdat %>% filter(zoopGrp == g)
    df_pred <- results[[g]]$predictions
    df_Q10 <- results[[g]]$Q10
    
    #Calculate sample size and number of taxa
    n_obs <- nrow(df_data)

    list(
      data = df_data %>% mutate(zoopGrp = g),
      pred = df_pred %>% mutate(zoopGrp = g),
      Q10_label = paste0("Q10 = ", round(df_Q10$Q10, 2),
                         " (95% CI: ", round(df_Q10$CI_lower, 2), " - ",
                         round(df_Q10$CI_upper, 2), "); n = ", n_obs)
    )
  })
  
  # Remove NULLs
  plot_data <- compact(plot_data)
  
  # Combine all data into single dataframes
  all_data <- bind_rows(map(plot_data, \(x) x$data))
  all_pred <- bind_rows(map(plot_data, \(x) x$pred))
  all_labels <- tibble(
    zoopGrp = map_chr(plot_data, \(x) as.character(unique(x$data$zoopGrp))),
    label = map_chr(plot_data, \(x) x$Q10_label),
    x = 0.0037,
    y = Inf
  )
  
  # Re-order zoopGrps if specified
  if (!is.null(group_order)) {
    all_data <- all_data %>%
      mutate(zoopGrp = fct_relevel(zoopGrp, group_order))
    
    all_pred <- all_pred %>%
      mutate(zoopGrp = fct_relevel(zoopGrp, group_order))
    
    all_labels <- all_labels %>%
      mutate(zoopGrp = fct_relevel(zoopGrp, group_order))
  }
  
  # Set x-axis limits... unless specified
  if (is.null(x_limits)) {
    x_limits <- c(max(all_data$x, na.rm = TRUE), 
                  min(all_data$x, na.rm = TRUE))
  }
  
  # Create the plot
  p <- ggplot() +
    geom_point(data = all_data, aes(x = x, y = .data[[rate_col]]), 
               alpha = 0.3, size = 2) +
    geom_ribbon(data = all_pred, aes(x = x, ymin = conf.low, ymax = conf.high),
                fill = "skyblue", alpha = 0.3) +
    geom_line(data = all_pred, aes(x = x, y = predicted), 
              colour = "darkblue", linewidth = 1, linetype = "dashed") +
    geom_text(data = all_labels, aes(x = x, y = y, label = label), 
              hjust = 0, vjust = 1.5, size = 4, colour = "darkblue") +
    scale_x_continuous(name = expression(bold("1 / Temp (K"^-1*")")),
                       trans = "reverse",
                       sec.axis = sec_axis(~1/. - 273.15, name = "Temp (°C)")) +
    facet_wrap(~ zoopGrp, scales = "free_y", ncol = ncol) +
    coord_cartesian(xlim = x_limits) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.title.x.top = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 12, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(p)
}
# END OF ARRHENIUS PLOT FUNCTION

