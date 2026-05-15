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

# RATIO CALCULATOR ----
calcRatios <- function(params) {
  # Solving y = mx + b for each physiological rate
  # where, m = slope, b = intercept
  # x = 15degC (i.e., mid point of our temp plots)
  # Note, our estimates are log-transformed rate data, so we will need to back-transform
  
  # Get rate-specific parameters
  G <- params %>% 
    filter(rate_name == "Growth")
  R <- params %>% 
    filter(rate_name == "Respiration") # this is in ulO2 mgC hr...will need to convert to carbon
  I <- params %>% 
    filter(rate_name == "Ingestion")
  
  refT = 15 # Reference temp
  
  RQ <- 0.8 # Respiratory quotient based on Hirst and Sheader and Ikeda and Motoda
  
  RQconv <- RQ * (12.011 / 22.4) * (1 / 1000) # create the RQ conversion
  
  # Back-transform with exp() before calculating the ratio because my response was log-transformed
  G_mgC <- exp((G$slope * refT) + G$intercept)
  R_uLO2 <- exp((R$slope * refT) + R$intercept)
  R_mgC <- R_uLO2 * RQconv # take uLO2 and convert to mgC
  I_mgC <- exp((I$slope * refT) + I$intercept)
  
  # Build a dataframe with the ratios 
  data.frame(
    ratio = c("G:I", # Gross (or ecological) growth efficiency - what is the proportion of energy ingested relative to used for growth
              "G:R", # Metabolic efficiency - what is the proportion of energy respired relative to energy for growth?
              "R:I"), # Metabolic expense - what is the proportion of ingested energy relative to energy respired?
    value = c(round(G_mgC / I_mgC, digits = 2), # GGE = G:I
              round(G_mgC / R_mgC, digits = 2), # G:R
              round(R_mgC / I_mgC, digits = 2)), # R:I
    description = c("Gross growth efficiency",
                    "Metabolic efficiency",
                    "Metabolic expense"),
    referenceTemp = refT
  )
}
# END OF RATIO CALCULATOR

# Estimate ratios ----
ratios <- calcRatios(allZ_params)
ratios








