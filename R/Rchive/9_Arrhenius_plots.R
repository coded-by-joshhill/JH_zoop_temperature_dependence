# Generating Arrhenius plots
# Josh Hill
# 27/01/2026


  # Here I read in the saved Arrhenius plot object and
  # Generate Arrhenius plots showing the relationship between C-specific rates and temperature (for supp materials)



# Packages and helpers ----
library(tidyverse)
source("R/0_Helpers.R")


##### CLEARANCE RATE ----
# Read in the plotting object
obj <- readRDS("Data/clearance_arrhenius_plot_data.rds")

# Generate the plot
clearance_plot <- arrhenius_plot(
  mdat = obj$mdat,
  rate_col = "Cspecific_rate",
  results = obj$results,
  group_order = obj$group_order,
  x_limits = c(0.0037, 0.0033)) +
  labs(y = expression(bold("Clearance rate (ml mgC"^-1*" h"^-1*")")))

clearance_plot

# Save it
ggsave("Output/Figure_4.png", plot = clearance_plot,
       width = 11.2, height = 6.5,  dpi = 300)
ggsave("Output/Figure_4.pdf", plot = clearance_plot,
       width = 11.2, height = 6.5,  dpi = 300)



##### INGESTION RATE ----
# Read in the plotting object
obj <- readRDS("Data/ingestion_arrhenius_plot_data.rds")

# Generate the plot
ingestion_plot <- arrhenius_plot(
  mdat = obj$mdat,
  rate_col = "Cspecific_rate",
  results = obj$results,
  group_order = obj$group_order,
  x_limits = c(0.0037, 0.0033)) +
  labs(y = expression(bold("Ingestion rate (mgC mgC"^-1*" h"^-1*")")))

ingestion_plot

# Save it
ggsave("Output/Figure_5.png", plot = ingestion_plot,
       width = 11.2, height = 4,  dpi = 300)
ggsave("Output/Figure_5.pdf", plot = ingestion_plot,
       width = 11.2, height = 4,  dpi = 300)



##### GROWTH RATE ----
# Read in the plotting object
obj <- readRDS("Data/growth_arrhenius_plot_data.rds")

# Generate the plot
growth_plot <- arrhenius_plot(
  mdat = obj$mdat,
  rate_col = "Cspecific_rate",
  results = obj$results,
  group_order = obj$group_order,
  x_limits = c(0.0037, 0.0033)) +
  labs(y = expression(bold("Growth rate (mgC mgC"^-1*" h"^-1*")")))

growth_plot

# Save it
ggsave("Output/Figure_6.png", plot = growth_plot,
       width = 11.2, height = 9,  dpi = 300)
ggsave("Output/Figure_6.pdf", plot = growth_plot,
       width = 11.2, height = 9,  dpi = 300)



##### RESPIRATION RATE ----
# Read in the plotting object
obj <- readRDS("Data/respiration_arrhenius_plot_data.rds")

# Generate the plot
respiration_plot <- arrhenius_plot(
  mdat = obj$mdat,
  rate_col = "Cspecific_rate",
  results = obj$results,
  group_order = obj$group_order,
  x_limits = c(0.0037, 0.0033)) +
  labs(y = expression(bold("Respiration rate (" * mu * "LO"[2] * " mgC"^-1 * " h"^-1 * ")")))

respiration_plot

# Save it
ggsave("Output/Figure_7.png", plot = respiration_plot,
       width = 11.2, height = 6.5,  dpi = 300)
ggsave("Output/Figure_7.pdf", plot = respiration_plot,
       width = 11.2, height = 6.5,  dpi = 300)

