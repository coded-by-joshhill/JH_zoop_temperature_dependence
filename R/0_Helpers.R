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
  if (unit == "ml/mgC/hr")  return(list(rate = rate,                unit = "ml/mgC/hr")) # maintain ml/mgC/hr
  if (unit == "ml/mgC/day") return(list(rate = rate / 24,           unit = "ml/mgC/hr")) # day to hr
  
  return(list(rate = NA_real_, unit = NA_character_))
}
# END OF CLEARANCE CONVERTER



# INGESTION CONVERTER ----
# Convert ingestion rates to mgC/ind/hr
convert_ingestion <- function(rate, unit) {
  if (is.na(rate) || is.na(unit)) return(list(rate = NA_real_, unit = NA_character_))
  
  if (unit == "mgC/ind/hr") return(list(rate = rate,                unit = "mgC/ind/hr")) # maintain mgC/ind/hr
  if (unit == "mgC/mgC/hr") return(list(rate = rate,                unit = "mgC/mgC/hr")) # maintain mgC/mgC/hr
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
  # Relative Growth Rate (RGR) = (log(M1) - log(M0)) / t
  # Where M1 is mass at a time, M0 is mass at the previous time point, and t is the time period (delta time) between the two measurements
  RGR = (log(mass1) - log(mass0)) / time
  
  return(RGR)
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
    return(list(rate = ((rate / 24) / (12.011 * 0.9)) * 22.4,         unit = "ulO2/ind/hr")) # μgC to μlO2, where 0.9 = RQ (Ross1982), day to hr, 22.4 = uL/umol as molar volume of oxygen and 12.011 is ug/umol of molar mass carbon
  
  return(list(rate = NA_real_, unit = NA_character_))
}
# END OF RESPIRATION CONVERTER



# TAXONOMIC COVERAGE SUMMARY FUNCTION ----
summarise_taxonomic_coverage <- function(data, dataset_name) {
  coverage <- data %>%
    summarise(
      Dataset    = dataset_name,
      n_phylas   = n_distinct(phylum, na.rm = TRUE),
      n_classes  = n_distinct(class, na.rm = TRUE),
      n_orders   = n_distinct(order, na.rm = TRUE),
      n_families = n_distinct(family, na.rm = TRUE),
      n_genera   = n_distinct(genus),
      n_species  = n_distinct(species),
      n_observations = n(),
      n_records = n_distinct(primRef)
    )
  return(coverage)
}
# END OF TAXONOMIC SUMMARY FUNCTION



# GROSS GROWTH EFFICIENCY CALCULATOR ----
calcGGE <- function(params) {
  # Solving y = mx + b for each physiological rate
  # where, m = slope, b = intercept
  # x = 15degC (i.e., mid point of our temp plots)
  # and GGE = Growth : Ingestion
  # Note, our estimates are log-transformed rate data, so we will need to back-transform
  
  # Get growth parameters
  G <- params %>% 
    filter(rate_name == "Growth")
  
  # Get ingestion paramaters
  I <- params %>% 
    filter(rate_name == "Ingestion")
  
  # Back-transform with exp() before calculating the ratio because my response was log-transformed
  GGE <- exp((G$slope * 15) + (G$intercept)) / exp((I$slope * 15) + (I$intercept))
  
  return(GGE)
  
}
# END OF GGE CALCULATOR




