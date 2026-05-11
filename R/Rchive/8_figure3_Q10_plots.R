# Calculating respiration rate Q10s and plotting
# Josh Hill
# 8/12/25



  # Here I read in all Q10 estimates
  # Plot Q10



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

  # Plot it up
  Q10plot <- plotQ10(Q10pdat, data_type = "Clearance", colours = mycols) +
    theme(axis.title.x = element_blank())
  Q10plot
  
(clearancePlot <- Q10plot)


  
# Read in ingestion rate data ----
Q10pdat <- readRDS("Data/Q10_summary_ingestion.rds")

  # Plot it up
  Q10plot <- plotQ10(Q10pdat, data_type = "Ingestion", colours = mycols) +
    theme(axis.title.x = element_blank())
  Q10plot
  
(ingestionPlot <- Q10plot)



# Read in growth rate data ----
Q10pdat <- readRDS("Data/Q10_summary_growth.rds")

  # Plot it up
  Q10plot <- plotQ10(Q10pdat, data_type = "Growth", colours = mycols) +
    theme(axis.title.x = element_blank())
  Q10plot
  
(growthPlot <- Q10plot)



# Read in respiration rate data ----
Q10pdat <- readRDS("Data/Q10_summary_respiration.rds")

  # Plot it up
  Q10plot <- plotQ10(Q10pdat, data_type = "Respiration", colours = mycols) +
    theme(axis.title.x = element_blank())
  Q10plot
  
(respirationPlot <- Q10plot)
  
  

# Combine plots ----
jointPlots <-ggarrange(clearancePlot, ingestionPlot, growthPlot, respirationPlot,
                     ncol = 2, nrow = 2, 
                     labels = c("a", "b", "c", "d"), 
                     common.legend = TRUE, 
                     legend = "bottom")
  
jointPlots <- annotate_figure(jointPlots,
                          bottom = text_grob("Zooplankton group", face = "bold", size = 14))
jointPlots

# Save it
# ggsave("Output/Figure_3.png", plot = jointPlots, width = 8, height = 8, background = "white")
# ggsave("Output/Figure_3.pdf", plot = jointPlots, width = 8, height = 8)


