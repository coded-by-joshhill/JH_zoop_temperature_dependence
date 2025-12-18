# Old thermal sensitivity analyser
# This used to estimate Q10 but was scrapped because it would cause anomolies in the confidence intervals due to calculating the confidence intervals from pooling std errors.


# OLD THERMAL SENSITIVITY ANALYSER ----
# Calculate Q10 with CI AND prediction ribbons for plotting
# This function does both calculations at once for efficiency
analyse_thermal_sensitivity <- function(model, data,
                                        temp_col = "temp_K", # for Q10 calculation
                                        x_col = "x", # for ribbon plotting
                                        delta_T = 10, # For Q10 calculation
                                        alpha = 0.05) # For 95% CI
{ 
  #### Calculate Q10 ----
  # Define groups and only progress for 1 group at a time
  grp <- unique(data$zoopGrp)
  stopifnot(length(grp) == 1)
  
  # Number of observations
  n_obs <- nrow(data)
  
  # Define temperatures
  T1 <- min(data[[temp_col]], na.rm = TRUE) # Take the minimum temperature per group
  T2 <- T1 + delta_T # Take minimum temp and add specified temperature (10 by default)
  
  # Create prediction data for just two temps
  newdat <- data.frame(x = c(1 / T1, 1 / T2), # A dataframe with JUST the predictor comprising the reciprocals of two temps
                       zoopGrp = grp) 
  
  # Get predictions with std err
  fit <- predict(model, newdata = newdat, # Predict the responses (in ln space) for the two temps
                 se.fit = TRUE, 
                 re.form = NA) # Maintain population-level predictions (ie., ignore random effects)
  
  # Extract preds and SEs
  pred1 <- fit$fit[1]
  pred2 <- fit$fit[2]
  se1 <- fit$se.fit[1]
  se2 <- fit$se.fit[2]
  
  # Calculate Q10 value
  Q10 <- exp(pred2 - pred1) # The ratio of rates at the warmer vs cooler temps
  
  # Calculate SE of the difference
  se_diff <- sqrt(se1^2 + se2^2)
  
  # CI on log scale then transform
  log_Q10 <- pred2 - pred1
  crit <- qnorm(1 - alpha / 2)
  CI_lower_log <- log_Q10 - crit * se_diff
  CI_upper_log <- log_Q10 + crit * se_diff
  
  # Store Q10 results
  Q10_results <- list(
    Q10 = as.numeric(Q10),
    CI_lower = as.numeric(exp(CI_lower_log)),
    CI_upper = as.numeric(exp(CI_upper_log)),
    T1_C = T1 - 273.15,
    T2_C = T2 - 273.15,
    se_Q10 = as.numeric(Q10 * se_diff),
    n_obs = n_obs
  )
  
  #### Get prediction ribbon for plotting ----
  # Get predictions with standard errors for all temps in dataset
  preds_raw <- predict(model,
                       newdata = data.frame(x = data[[x_col]],
                                            zoopGrp = grp),
                       se.fit = TRUE,
                       re.form = NA # Maintain population-level predictions (ie., ignore random effects)
  )
  
  # Create dataframe with predictions and CIs
  preds_df <- data.frame(
    x = data[[x_col]],
    predicted_log = preds_raw$fit,
    conf.low_log = preds_raw$fit - crit * preds_raw$se.fit,
    conf.high_log = preds_raw$fit + crit * preds_raw$se.fit
  ) %>%
    mutate(
      predicted = exp(predicted_log),
      conf.low = exp(conf.low_log),
      conf.high = exp(conf.high_log)
    )
  
  # Return both Q10 and predictions
  return(list(
    Q10 = Q10_results,
    predictions = preds_df
  ))
}
# END OF THERMAL SENSITIVITY ANALYSER



# MASS-SPECIFIC RATE VIOLIN PLOT FUNCTION
plotRates <- function(mdat, data_type, colours){
  
  # Get the rate-specific colour and unit
  rateCol <- colours [[data_type]]
  
  p <- mdat %>%
    ggplot(aes(x = zoopGrp, y = log(Cspecific_rate))) +
    geom_violin(fill = rateCol, 
                alpha = 0.5) +
    geom_jitter(width = 0.2, 
                alpha = 0.3, 
                size = 1) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.title.x.top = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12)) +
    labs(
      x = expression(bold("Zooplankton group"))) 
  
  return(p)
  
}
# END OF RATE VIOLIN PLOT FUNCTION