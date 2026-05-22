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
  
  RQconv <- RQ * (12.011 / 22.4) * (1 / 1000) # create the RQ conversion and convert from ug to mg
  
  # Back-transform with exp() before calculating the ratio because my response was log-transformed
  G_mgC <- exp((G$slope * refT) + G$intercept) 
  R_uLO2 <- exp((R$slope * refT) + R$intercept)
  R_mgC <- R_uLO2 * RQconv
  I_mgC <- exp((I$slope * refT) + I$intercept)
  
  # Build a dataframe with the ratios 
  data.frame(
    ratio = c("G:I", # Gross (or ecological) growth efficiency - what is the proportion of ingested energy used for growth
              "G:R", # Metabolic efficiency - what is the proportion of energy respired relative to energy for growth?
              "R:I", # Metabolic expense - what is the proportion of ingested energy relative to energy respired?
              "G+R:I"# Assimilation efficiency - what is the proportion of ingested energy that is digested?
    ),
    value = c(round(G_mgC / I_mgC, digits = 2), # GGE = G:I
              round(G_mgC / R_mgC, digits = 2), # G:R
              round(R_mgC / I_mgC, digits = 2), # R:I
              round((G_mgC + R_mgC) / I_mgC, digits = 2) # AE
    ),
    description = c("Gross growth efficiency",
                    "Metabolic efficiency",
                    "Metabolic expense",
                    "Assimilation efficiency"
    ),
    referenceTemp = refT
  )
}
# END OF RATIO CALCULATOR

# Estimate ratios ----
ratios <- calcRatios(allZ_params)
ratios


temp_seq <- seq(-2, 32, by = 0.5)

RQ     <- 0.8
RQconv <- RQ * (12.011 / 22.4) * (1 / 1000)

rateCols <- c(
  "Ingestion" = "green2",
  "Respiration" = "red2",
  "Growth" = "slateblue"
)

pred_df <- allZ_params %>%
  filter(rate_name != "Clearance",
         rate_name != "Excretion") %>%
  rowwise() %>%
  mutate(temp = list(temp_seq)) %>%
  unnest(temp) %>%
  mutate(rate_pred = case_when(
    rate_name == "Respiration" ~ exp(slope * temp + intercept) * RQconv,
    TRUE                       ~ exp(slope * temp + intercept)
  ))

# Build annotation dataframe
slope_labels <- allZ_params %>%
  filter(rate_name %in% c("Ingestion", "Growth", "Respiration")) %>%
  mutate(
    label = paste0("slope = ", round(slope, 3)),
    temp  = 25,  # x position of label
    rate_pred = case_when(
      rate_name == "Respiration" ~ exp(slope * 25 + intercept) * RQconv,
      TRUE                       ~ exp(slope * 25 + intercept)
    )
  )

ggplot(pred_df, aes(x = temp, y = log10(rate_pred), colour = rate_name, linetype = rate_name)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 15, linetype = "dashed", colour = "grey60") +
  geom_text(data = slope_labels,
            aes(x = temp, y = log10(rate_pred), label = label, colour = rate_name),
            nudge_y = .2, # adjust as needed
            hjust = .5,
            size = 3.5,
            show.legend = FALSE) +
  scale_colour_manual(values = rateCols) +
  labs(
    x        = "Temperature (°C)",
    y        = expression("log10 Rate (mgC mgC"^-1~"hr"^-1*")"),
    colour   = NULL,
    linetype = NULL
  ) +
  theme_bw()

