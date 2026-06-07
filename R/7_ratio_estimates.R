# Estimating ratios of the physiological terms
# Josh Hill
# 11/05/2026



# Packages and helpers ----
library(tidyverse)
source("R/0_Helpers.R")



# Read in the paramater data ----
allZ_params <- readRDS(file = "Data/modelParameters/allZestimates.rds")
sizeGrp_params <- readRDS(file = "Data/modelParameters/sizeGrpestimates.rds")
funcGrp_params <- readRDS(file = "Data/modelParameters/funcGrpestimates.rds")
zoopGrp_clearance_params <- readRDS(file = "Data/modelParameters/ZGrpClearanceEstimates.rds")
zoopGrp_ingestion_params <- readRDS(file = "Data/modelParameters/ZGrpIngestionEstimates.rds")
zoopGrp_growth_params <- readRDS(file = "Data/modelParameters/ZGrpGrowthEstimates.rds")
zoopGrp_respiration_params <- readRDS(file = "Data/modelParameters/ZGrpRespirationEstimates.rds")
zoopGrp_excretion_params <- readRDS(file = "Data/modelParameters/ZGrpExcretionEstimates.rds")

# Combine the zoopGrp params into one object
zoopGrp_params <- bind_rows(
  readRDS("Data/modelParameters/ZGrpClearanceEstimates.rds")   %>% mutate(rate_name = "Clearance"),
  readRDS("Data/modelParameters/ZGrpIngestionEstimates.rds")   %>% mutate(rate_name = "Ingestion"),
  readRDS("Data/modelParameters/ZGrpGrowthEstimates.rds")      %>% mutate(rate_name = "Growth"),
  readRDS("Data/modelParameters/ZGrpRespirationEstimates.rds") %>% mutate(rate_name = "Respiration"),
  readRDS("Data/modelParameters/ZGrpExcretionEstimates.rds")   %>% mutate(rate_name = "Excretion"))

# Read in the mdat objects (we'll use these to get the number of taxa the efficienices represent)
allZ_mdat <- readRDS(file = "Data/modelParameters/allZoop_mdat.rds")
sizeGrp_mdat <- readRDS(file = "Data/modelParameters/sizeGrp_mdat.rds")
funcGrp_mdat <- readRDS(file = "Data/modelParameters/funcGrp_mdat.rds")
zoopGrp_mdat <- readRDS(file = "Data/modelParameters/zoopGrp_mdat.rds")



# Calculate ratios ----

# RATIO CALCULATOR ----
calcRatios <- function(params, groupVar = NULL) {
  # Solving y = mx + b for each physiological rate
  # where, m = slope, b = intercept
  # x = 15degC (i.e., mid point of our temp plots)
  # Note, our estimates are log-transformed rate data, so we will need to back-transform
  refT = 15 # Reference temp
  
  RQ <- 0.8 # Respiratory quotient based on Hirst and Sheader and Ikeda and Motoda
  
  RQconv <- RQ * (12.011 / 22.4) * (1 / 1000) # create the RQ conversion and convert from ug to mg
  
  calcSingle <- function(df, grpLabel = NULL) { # Get rate-specific parameters
  G <- df %>% filter(rate_name == "Growth")
  R <- df %>% filter(rate_name == "Respiration") # this is in ulO2 mgC hr...will need to convert to carbon
  I <- df %>% filter(rate_name == "Ingestion")

  if(nrow(G) == 0 | nrow(R) == 0 | nrow(I) == 0) {
    warning(paste("NO RATES FOR GROUPING VARIABLE:", grpLabel))
    return(data.frame())  # return empty data frame instead of NULL when no data exists
    return(NULL)
  }
  
  # Back-transform with exp() before calculating the ratio because my response was log-transformed
  G_mgC <- exp((G$slope * refT) + G$intercept) 
  R_uLO2 <- exp((R$slope * refT) + R$intercept)
  R_mgC <- R_uLO2 * RQconv
  I_mgC <- exp((I$slope * refT) + I$intercept)
  
  # Build a dataframe with the ratios 
  data.frame(
    ratio = c("G+R:I", # Assimilation efficiency - what is the proportion of ingested energy that is absorbed across the gut wall?
              "G:I", # Gross (or ecological) growth efficiency - what is the proportion of ingested energy used for growth
              "G:G+R" # Net growth efficiency - what is the proportion of energy that is absorbed used ?
    ),
    value = c(round(((G_mgC + R_mgC) / I_mgC) * 100, digits = 0), # AE -> assuming A ~= to G + R
              round((G_mgC / I_mgC) * 100, digits = 0), # GGE = G:I
              round((G_mgC / (G_mgC + R_mgC)) * 100, digits = 0) # NGE -> assuming A ~= to G + R
    ),
    description = c("AE (G+R:I)",
                    "GGE (G:I)",
                    "NGE (G:G+R)"
    ),
    referenceTemp = refT
    )
  }
  
if(is.null(groupVar)) {
  return(calcSingle(params))
}

params %>% 
  group_by(across(all_of(groupVar))) %>% 
  group_modify(~ calcSingle(.x, grpLabel = as.character(.y[[groupVar]]))) %>% 
  ungroup()
}
# END OF RATIO CALCULATOR



# Estimate ratios ----
ratios_allZ    <- calcRatios(allZ_params)
ratios_allZ

ratios_sizeGrp <- calcRatios(sizeGrp_params, groupVar = "sizeGrp")
ratios_sizeGrp

ratios_funcGrp <- calcRatios(funcGrp_params, groupVar = "funcGrp")
ratios_funcGrp

ratios_zoopGrp <- calcRatios(zoopGrp_params, groupVar = "zoopGrp")
ratios_zoopGrp

# Count taxa per level from mdats
n_taxa_all <- bind_rows(
  allZ_mdat %>%
    filter(rate_name %in% c("Growth", "Respiration", "Ingestion")) %>%
    summarise(n_taxa = n_distinct(taxa)) %>%
    mutate(group = "All zooplankton"),
  sizeGrp_mdat %>%
    filter(rate_name %in% c("Growth", "Respiration", "Ingestion")) %>%
    group_by(sizeGrp) %>%
    summarise(n_taxa = n_distinct(taxa), .groups = "drop") %>%
    rename(group = sizeGrp),
  funcGrp_mdat %>%
    filter(rate_name %in% c("Growth", "Respiration", "Ingestion")) %>%
    group_by(funcGrp) %>%
    summarise(n_taxa = n_distinct(taxa), .groups = "drop") %>%
    rename(group = funcGrp),
  zoopGrp_mdat %>%
    filter(rate_name %in% c("Growth", "Respiration", "Ingestion")) %>%
    group_by(zoopGrp) %>%
    summarise(n_taxa = n_distinct(taxa), .groups = "drop") %>%
    rename(group = zoopGrp)
)

# Combine them together 
ratios_combined <- bind_rows(
  ratios_allZ    %>% mutate(level = "All zooplankton", group = "All zooplankton"),
  ratios_sizeGrp %>% rename(group = sizeGrp) %>% mutate(level = "Size group"),
  ratios_funcGrp %>% rename(group = funcGrp) %>% mutate(level = "Functional group"),
  ratios_zoopGrp %>% rename(group = zoopGrp) %>% mutate(level = "Zooplankton group")) %>%
  select(group, description, value) %>%
  pivot_wider(names_from = description, values_from = value) %>%
  left_join(n_taxa_all, by = "group") %>%
  relocate(n_taxa, .after = everything())

ratios_combined

