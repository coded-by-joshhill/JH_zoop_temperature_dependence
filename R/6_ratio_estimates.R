# Estimating ratios of the physiological terms
# Josh Hill
# 11/05/2026



# Packages and helpers ----
library(tidyverse)
source("R/0_Helpers.R")



# Read in the data ----
allZ_params <- readRDS(file = "Data/modelParameters/allZestimates.rds")
allZ_params


# Calculate ratios ----

# Estimate Gross Growth Efficiency based on GGE = G:I
GGE <- data.frame(meanGGE = calcGGE(allZ_params))
GGE




# GVR CALCULATOR ----
calcGvR <- function(params) {
  # Solving y = mx + b for each physiological rate
  # where, m = slope, b = intercept
  # x = 15degC (i.e., mid point of our temp plots)
  # and GvR = Growth : Respiration
  # Note, our estimates are log-transformed rate data, so we will need to back-transform
  
  # Get growth parameters
  G <- params %>% 
    filter(rate_name == "Growth")
  
  # Get respiration paramaters
  R <- params %>% 
    filter(rate_name == "Respiration")
  
  I <- params %>% 
    filter(rate_name == "Ingestion")
  
  RQ <- 0.9 # based on Ross 1982...could try 0.97 based on McKinnon 
  
  RQconv <- (1 / RQ) * (12.011 / 22.4) * (1 / 1000)
  
  G_mgC <- exp((G$slope * 15) + (G$intercept))
  R_uLO2 <- exp((R$slope * 15) + (R$intercept))
  R_mgC <- R_uLO2 * RQconv
  I_mgC <- exp((I$slope * 15) + (I$intercept))
  
  # Back-transform with exp() before calculating the ratio because my response was log-transformed
  GvR <- G_mgC / R_mgC
  RvI <- R_mgC / I_mgC
  
  return(GvR)
  
}
# END OF GvE CALCULATOR

# Estimate <add ratio name> G:R
GvR <- data.frame(meanGvR = calcGvR(allZ_params))
GvR






