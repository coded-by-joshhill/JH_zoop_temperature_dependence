# Helper functions for temperature sensitivity data analysis for all zooplankton
# Josh Hill
# 09/02/2026



# Packages ----
library(tidyverse)



############# FUNCTIONS FOR DATA CLEANING AND ANALYSIS ############# 



# BODY MASS AS CARBON CALCULATOR
# Estimate carbon weight using an allometric equation
calc_BMC <- function(bodyLength_mm) {
  # Following Jaspers et al. 2009, use the slope (2.455 μgC um), intercept (-6.96 μgC), and trunk length (TL) to estimate carbon mass
  # where, logC (μgC) = 2.455 log TL(μm) -6.96.... located in Table 1 and Figure 2 - DOI: 10.1093/plankt/fbp002 
  # therefore... C (μgC) = 10^-6.96 * (bodyLength (mm) * 1000) ^ 2.455
  carbonMass_ugC = 10^(-6.96) * (bodyLength_mm * 1000) ^ 2.455
  
  return(carbonMass_ugC)
}
# END OF BMC CALCULATOR



# CARBON WEIGHT CONVERTER ----
# Convert carbon weights to mg
convert_CW <- function(weight, unit) {
  if(is.na(weight) || is.na(unit)) return(NA_real_)
  if(unit == "mg") return(weight)         # maintain mg
  if(unit == "ug") return(weight / 1000)  # µg to mg
  if(unit == "ng") return(weight / 1e6)   # ng to mg
  if(unit == "g")  return(weight * 1000)  # g to mg
  if(unit == "kg") return(weight * 1e6)   # kg to mg
  
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


# GROWTH CALCULATOR ----
# Estimate growth
calcGrowthRate <- function(mass1, mass0, time) {
  # Following McConville and colleagues (2017), estimate growth using:
  # G = (M1 - M0) / t
  # Where M1 is mass at a time, M0 is mass at the previous time point, and t is the time period between the two measurements
  G = (mass1 - mass0) / time
  
  return(G)
}
# END OF GROWTH CALCULATOR



# RESPIRATION CONVERTER ----
# Convert respiration rates to μlO2/ind/hr
convert_respiration <- function(rate, unit, genus) {
  if (is.na(rate) || is.na(unit)) return(list(rate = NA_real_, unit = NA_character_))
  
  if (unit == "ulO2/ind/hr") return(list(rate = rate,                  unit = "ulO2/ind/hr")) # maintain μlO2/ind/hr
  if (unit == "ulO2/ind/day") return(list(rate = rate / 24,            unit = "ulO2/ind/hr")) # day to hr
  if (unit == "ulO2/mgC/hr") return(list(rate = rate,                  unit = "ulO2/mgC/hr")) # maintain μlO2/mgC/hr
  if (unit == "ulO2/ugC/day") return(list(rate = (rate / 1000) / 24,   unit = "ulO2/mgC/hr")) # μgC to mgC
  if (unit == "umol/ind/hr") return(list(rate = rate * 22.4,           unit = "ulO2/ind/hr")) # μmol to μlO2, as per the ideal gas law (i.e., 1 mol of gas = 22.4 L)
  if (unit == "nlO2/ind/hr") return(list(rate = rate / 1000,           unit = "ulO2/ind/hr")) # nlO2 to μlO2
  
  # Stoichiometry
  if (genus == "Euphausia" & unit == "ugC/ind/day") 
    return(list(rate = ((rate / 24) / (12.011 * 0.9)) * 22.4,         unit = "ulO2/ind/hr")) # μgC to μlO2, where 0.9 = RQ (Ross1982), day to hr
  
  return(list(rate = NA_real_, unit = NA_character_))
}
# END OF RESPIRATION CONVERTER



# TAXONOMIC COVERAGE SUMMARY FUNCTION ----
summarise_taxonomic_coverage <- function(data, dataset_name) {
  coverage <- data %>%
    summarise(
      Dataset    = dataset_name,
      n_species  = n_distinct(species),
      n_genera   = n_distinct(genus),
      n_families = n_distinct(family, na.rm = TRUE),
      n_orders   = n_distinct(order, na.rm = TRUE),
      n_classes  = n_distinct(class, na.rm = TRUE),
      n_phylas   = n_distinct(phylum, na.rm = TRUE),
      n_zoopGrps = n_distinct(zoopGrp, na.rm = TRUE),
      n_observations = n()
    )
  return(coverage)
}
# END OF TAXONOMIC SUMMARY FUNCTION



######################### FUNCTIONS FOR ALLZOOP ANALYSIS AND PLOTTING #########################



# THERMAL SENSITIVITY ANALYSER ----
# Calculate Q10 with CI AND prediction ribbons for plotting
analyse_thermal_sensitivity <- function(model, data,
                                        temp_col = "temp_K", # for Q10 calculation
                                        x_col = "x", # for ribbon plotting
                                        alpha = 0.05) # For 95% CI
{ 
  #### Calculate Q10 ----
  # Define groups and only progress for 1 group at a time
  grp <- unique(data$group)
  stopifnot(length(grp) == 1)
  
  # Number of observations
  n_obs <- nrow(data)
  
  # Get degrees of freedom from model
  df_residual <- df.residual(model)
  
  # Define temperatures
  T1 <- min(data[[temp_col]], na.rm = TRUE) # Take the minimum temperature per group
  T2 <- T1 + 10 # Take minimum temp and add specified temperature (hardwired to 10 degC)
  
  # Create prediction data for just two temps
  newdat <- data.frame(x = c(1 / T1, 1 / T2), # A dataframe with JUST the predictor comprising the reciprocals of two temps
                       group = grp) 
  
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
  crit <- qt(1 - alpha / 2, df = df_residual) 
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
  
  #### Get confidence ribbon for plotting ----
  # Get predictions with standard errors for all temps in dataset
  preds_raw <- predict(model,
                       newdata = data.frame(
                         x = data[[x_col]],
                         group = grp,
                         primRef = NA,
                         taxa = NA),
                       se.fit = TRUE,
                       re.form = NA # ignore random effects)
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



# ALL ZOOP ARRHENIUS PLOT FUNCTION ----
arrhenius_plot <- function(arrhenius_obj,
                           x_limits = NULL) {

  # Create Arrhenius plot (in log space because our rates vary in orders of magnitude)
  p <- ggplot() +
    # raw points 
    geom_point(data = arrhenius_obj$all_data, 
               aes(x = x, y = log(Cspecific_rate)), 
               alpha = 0.3, size = 1) +
    # prediction ribbon
    geom_ribbon(data = arrhenius_obj$all_pred,
                aes(x = x, ymin = conf.low_log, ymax = conf.high_log),
                fill = "lightgrey", alpha = 0.3) +
    # predicted line
    geom_line(data = arrhenius_obj$all_pred, 
              aes(x = x, y = predicted_log), 
              colour = "darkblue", linewidth = 1, linetype = "dashed") +
    # reverse x-axis for Arrhenius
    scale_x_continuous(name = expression(bold("1 / Temp (K"^-1*")")),
                       trans = "reverse",
                       sec.axis = sec_axis(~1/. - 273.15, name = "Temp (°C)")) +
    coord_cartesian(xlim = c(0.0037, 0.0033)) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 11),
          axis.title.y = element_text(size = 12),
          axis.title.x.top = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(t = 5, r = 15, b = 5, l = 5))
  
  return(p)

}
# END OF ARRHENIUS PLOT FUNCTION
