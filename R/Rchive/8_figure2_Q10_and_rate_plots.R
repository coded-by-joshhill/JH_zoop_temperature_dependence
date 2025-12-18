# Calculating respiration rate Q10s and plotting
# Josh Hill
# 8/12/25



  # Here I read in all Q10 estimations
  # Plot Q10 and CMass-Specific rates
  # Plot Arrhenius plots



# Packages and helpers ----
library(tidyverse)
library(patchwork)
library(ggpubr)
source("R/0_Helpers.R")



# Define color palette ----
mycols <- c("Clearance"   = "#66c2a5", 
            "Ingestion"   = "#fc8d62", 
            "Growth"      = "#8da0cb", 
            "Respiration" = "#e78ac3")



# Read in clearance rate data ----
Q10pdat <- readRDS("Data/Q10_summary_clearance.rds")
mdat <- readRDS("Data/clearance_mdat.rds")

  # Plot it up
  # Dot plot of Q10
  Q10plot <- plotQ10(Q10pdat, data_type = "Clearance", colours = mycols) +
    theme(axis.title.x = element_blank())
  Q10plot
  
  # Violin plot of log(Cspecific_rate)
  ratePlot <- plotRates(mdat, data_type = "Clearance", colours = mycols) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    labs(y = expression(bold("Clearance rate ln(mL mgC"^-1*" h"^-1*")")))
  ratePlot
  
(clearancePlots <- ratePlot / Q10plot +
    plot_layout(axis_titles = "collect_x"))


  
# Read in ingestion rate data ----
Q10pdat <- readRDS("Data/Q10_summary_ingestion.rds")
mdat <- readRDS("Data/ingestion_mdat.rds")

  # Plot it up
  # Dot plot of Q10
  Q10plot <- plotQ10(Q10pdat, data_type = "Ingestion", colours = mycols) +
    theme(axis.title.x = element_blank())
  Q10plot
  
  # Violin plot of log(Cspecific_rate)
  ratePlot <- plotRates(mdat, data_type = "Ingestion", colours = mycols) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    labs(y = expression(bold("Ingestion rate ln(mgC mgC"^-1*" h"^-1*")")))
  ratePlot
  
(ingestionPlots <- ratePlot / Q10plot +
    plot_layout(axis_titles = "collect_x"))



# Read in growth rate data ----
Q10pdat <- readRDS("Data/Q10_summary_growth.rds")
mdat <- readRDS("Data/growth_mdat.rds")
  
  # Plot it up
  # Dot plot of Q10
  Q10plot <- plotQ10(Q10pdat, data_type = "Growth", colours = mycols) +
    theme(axis.title.x = element_blank())
  Q10plot
  
  # Violin plot of log(Cspecific_rate)
  ratePlot <- plotRates(mdat, data_type = "Growth", colours = mycols) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    labs(y = expression(bold("Growth rate ln(mgC mgC"^-1*" h"^-1*")")))
  ratePlot
  
(growthPlots <- ratePlot / Q10plot +
    plot_layout(axis_titles = "collect_x"))



# Read in respiration rate data ----
Q10pdat <- readRDS("Data/Q10_summary_respiration.rds")
mdat <- readRDS("Data/respiration_mdat.rds")
  
  # Plot it up
  # Dot plot of Q10
  Q10plot <- plotQ10(Q10pdat, data_type = "Respiration", colours = mycols) +
    theme(axis.title.x = element_blank())
  Q10plot
  
  # Violin plot of log(Cspecific_rate)
  ratePlot <- plotRates(mdat, data_type = "Respiration", colours = mycols) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    labs(y = expression(bold("Respiration rate ln(" * mu * "LO"[2] * " mgC"^-1 * " h"^-1 * ")")))
  ratePlot
  
(respirationPlots <- ratePlot / Q10plot +
    plot_layout(axis_titles = "collect_x"))


# Combine plots ----
jointPlots <-ggarrange(clearancePlots, ingestionPlots, growthPlots, respirationPlots,
                     ncol = 4, nrow = 1, 
                     labels = c("A", "B", "C", "D"), 
                     common.legend = TRUE, 
                     legend = "bottom")
  
jointPlots <- annotate_figure(jointPlots,
                          bottom = text_grob("Zooplankton group", face = "bold", size = 14))
jointPlots

# Save it
# ggsave("Output/Figure_2.png", plot = jointPlots, width = 15, height = 8, background = "white")

