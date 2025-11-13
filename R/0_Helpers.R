# Helper functions for clearance and ingestion data analysis
# 16/09/25
# Josh Hill



# Packages ----
library(tidyverse)



# Functions ----

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
convert_clearance <- function(rate, unit, BMC_mg) {
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
convert_ingestion <- function(rate, unit, BMC_mg) {
  if (is.na(rate) || is.na(unit)) return(list(rate = NA_real_, unit = NA_character_))
  
  if (unit == "mgC/ind/hr") return(list(rate = rate,                unit = "mgC/ind/hr")) # maintain mgC/ind/hr
  if (unit == "ugC/ind/hr") return(list(rate = rate / 1000,         unit = "mgC/ind/hr")) # µg to mg
  if (unit == "ugC/ind/day") return(list(rate = (rate / 1000) / 24, unit = "mgC/ind/hr")) # µg to mg, day to hr
  if (unit == "ngC/ind/hr") return(list(rate = rate / 1e6,          unit = "mgC/ind/hr")) # ng to mg
  if (unit == "ngC/ind/day") return(list(rate = (rate / 1e6) / 24,  unit = "mgC/ind/hr")) # ng to mg, day to hr
  if (unit == "ngC/mgC/hr") return(list(rate = rate / 1e6,          unit = "mgC/mgC/hr")) # ngC to mgC
  if (unit == "ugC/ugC/day") return(list(rate = rate / 24,          unit = "mgC/mgC/hr")) # mgC to mgC ratio cancels out, day to hr
  
  return(list(rate = NA_real_, unit = NA_character_))
}
# END OF INGESTION CONVERTER



# RESPIRATION CONVERTER ----
# Convert respiration rates to ulO2/ind/hr
convert_respiration <- function(rate, unit, genus) {
  if (is.na(rate) || is.na(unit)) return(list(rate = NA_real_, unit = NA_character_))
  
  if (unit == "ulO2/ind/hr") return(list(rate = rate,                 unit = "ulO2/ind/hr")) # maintain ulO2/ind/hr
  if (unit == "ulO2/ind/day") return(list(rate = rate / 24,           unit = "ulO2/ind/hr")) # day to hr
  if (unit == "ulO2/mgC/hr") return(list(rate = rate,                 unit = "ulO2/mgC/hr")) # maintain ulO2/mgC/hr
  if (unit == "ulO2/ugC/day") return(list(rate = rate / 1000,         unit = "ulO2/mgC/hr")) # ugC to mgC
  if (unit == "umol/ind/hr") return(list(rate = rate * 22.4,          unit = "ulO2/ind/hr")) # umol to ulO2, as per the ideal gas law (i.e., 1 mol of gas = 22.4 L)
  if (unit == "pmol/ind/hr") return(list(rate = rate * 2.24e-17,      unit = "ulO2/ind/hr")) # pmol to ulO2, as above, but scaled down
  if (unit == "nlO2/ind/hr") return(list(rate = rate / 1000,          unit = "ulO2/ind/hr")) # nlO2 to ulO2
  
  # Stoichiometry
  if (genus == "Acartia" || unit == "ugC/ugC/hr") 
    return(list(rate = ((rate / (12.011 * 0.66)) * (22.4 / 1000)),      unit = "ulO2/mgC/hr")) # ugC to ulO2, where 0.66 = RQ (Mayzaud2005), ugC to mgC for individuals mass
  if (genus == "Euphausia" || unit == "ugC/ind/day") 
    return(list(rate = ((rate / (12.011 * 1.45)) * (22.4 / 24)), unit = "ulO2/ind/hr")) # ugC to ulO2, where 1.45 = RQ (Mayzaud2005), day to hr
  
  return(list(rate = NA_real_, unit = NA_character_))
}
# END OF RESPIRATION CONVERTER


# THERMAL SENSITIVITY ANALYSER ----
# Calculate Q10 with CI AND prediction ribbons for plotting
# This function does both calculations at once for efficiency
analyse_thermal_sensitivity <- function(model, data,
                                        temp_col = "temp_K", # for Q10 calculation
                                        x_col = "x", # for ribbon plotting
                                        delta_T = 10, # For Q10 calculation
                                        alpha = 0.05) {
  # Calculate Q10 ----
    # Define temperatures
    T1 <- median(data[[temp_col]], na.rm = TRUE) # use median because it is the mid point of temperatures
    T2 <- T1 + delta_T # Take temp1 and add delta_T (i.e., specified temp difference)
  
    # Create prediction data for just two temps
    newdat <- data.frame(x = c(1 / T1, 1 / T2)) # A dataframe with JUST the predictor comprising the reciprocals of two temps
  
    # Get predictions with std err
    fit <- predict(model, newdata = newdat, se.fit = TRUE, re.form = NA) # Predict the responses (in ln space) for the two temps
  
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
      se_Q10 = as.numeric(Q10 * se_diff)
    )

  # Get prediction ribbon for plotting ----
    # Get predictions with standard errors for all temps in dataset
    preds_raw <- predict(model,
      newdata = data.frame(x = data[[x_col]]),
      se.fit = TRUE,
      re.form = NA
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

