# Pulling all Arrhenius and Q10 plots together for allZoop estimates
# Josh Hill
# 10/02/2026


  # Here I read in the saved Arrhenius plot object and
  # Generate Arrhenius plots showing the relationship between C-specific rates and temperature (for supp materials)



# Packages and helpers ----
library(tidyverse)
library(ggpubr)
library(ggbreak)
source("R/0_Helpers.R")




##### CLEARANCE RATE ----
# Read in the plotting object
arrhenius_obj_clear <- readRDS("Data/allZoop_analysis/allZoop_clearance_arrhenius_data.rds")


# Generate the plot
arrhenius_plot_clearance <- arrhenius_plot(
  arrhenius_obj_clear,
  x_limits = c(0.0037, 0.0033)) +
  # add Q10 annotation
  geom_text(aes(x = 0.0037, y = 12, 
                label = arrhenius_obj_clear$Q10_label),
            parse = TRUE,
            hjust = 0, vjust = 1.5,
            size = 3, colour = "darkblue") +
  # add y-lab
  labs(y = expression(bold(atop("ln(Clearance rate)", "(ml mgC"^-1*" h"^-1*")"))))

arrhenius_plot_clearance



##### INGESTION RATE ----
# Read in the plotting object
arrhenius_obj_ingest <- readRDS("Data/allZoop_analysis/allZoop_ingestion_arrhenius_data.rds")


# Generate the plot
arrhenius_plot_ingestion <- arrhenius_plot(
  arrhenius_obj_ingest,
  x_limits = c(0.0037, 0.0033)) +
  # add Q10 annotation
  geom_text(aes(x = 0.0037, y = 0, 
                label = arrhenius_obj_ingest$Q10_label),
            parse = TRUE,
            hjust = 0, vjust = 1.5,
            size = 3, colour = "darkblue") +
  # add y-lab
  labs(y = expression(bold(atop("ln(Ingestion rate)", "(mgC mgC"^-1*" h"^-1*")"))))

arrhenius_plot_ingestion



##### GROWTH RATE ----
# Read in the plotting object
arrhenius_obj_grow <- readRDS("Data/allZoop_analysis/allZoop_growth_arrhenius_data.rds")


# Generate the plot
arrhenius_plot_growth <- arrhenius_plot(
  arrhenius_obj_grow,
  x_limits = c(0.0037, 0.0033)) +
  # add Q10 annotation
  geom_text(aes(x = 0.0037, y = 0, 
                label = arrhenius_obj_grow$Q10_label),
            parse = TRUE,
            hjust = 0, vjust = 1.5,
            size = 3, colour = "darkblue") +
  # add y-lab
  labs(y = expression(bold(atop("ln(Growth rate)", "(mgC mgC"^-1*" h"^-1*")"))))
  
arrhenius_plot_growth



##### RESPIRATION RATE ----
# Read in the plotting object
arrhenius_obj_resp <- readRDS("Data/allZoop_analysis/allZoop_respiration_arrhenius_data.rds")


# Generate the plot
arrhenius_plot_respiration <- arrhenius_plot(
  arrhenius_obj_resp,
  x_limits = c(0.0037, 0.0033)) +
  # add Q10 annotation
  geom_text(aes(x = 0.0037, y = 5, 
                label = arrhenius_obj_resp$Q10_label),
            parse = TRUE,
            hjust = 0, vjust = 1.5,
            size = 3, colour = "darkblue") +
  # add y-lab
  labs(y = expression(bold(atop("ln(Respiration rate)", "(" * mu * "LO"[2] * " mgC"^-1 * " h"^-1 * ")"))))
  
arrhenius_plot_respiration



# Combine plots ----
jointArrhenius <-ggarrange(arrhenius_plot_clearance, arrhenius_plot_ingestion, arrhenius_plot_growth, arrhenius_plot_respiration,
                       ncol = 2, nrow = 2, 
                       labels = c("a", "b", "c", "d"), 
                       legend = "none")

jointArrhenius

# Save it
ggsave("Output/Figure_2_allZoopArrhenius.png", plot = jointArrhenius, width = 180, height = 140, unit = "mm",  dpi = 300)
ggsave("Output/Figure_2_allZoopArrhenius.pdf", plot = jointArrhenius, width = 180, height = 140, unit = "mm",  dpi = 300)



# Q10 plots ----

# Read in Q10 data and assign data names ----
Q10pdat_clear <- readRDS("Data/allZoop_analysis/allZoop_clearance_Q10_data.rds") %>% 
  mutate(data = "Clearance")
Q10pdat_ing <- readRDS("Data/allZoop_analysis/allZoop_ingestion_Q10_data.rds") %>% 
  mutate(data = "Ingestion")
Q10pdat_grow <- readRDS("Data/allZoop_analysis/allZoop_growth_Q10_data.rds") %>% 
  mutate(data = "Growth")
Q10pdat_resp <- readRDS("Data/allZoop_analysis/allZoop_respiration_Q10_data.rds") %>% 
  mutate(data = "Respiration")

# Combine the datasets
Q10pdat <- bind_rows(Q10pdat_clear, Q10pdat_ing, Q10pdat_grow, Q10pdat_resp)

# Plot it up

# Define color palette ----
mycols <- c("Clearance"   = "#66c2a5",
            "Ingestion"   = "#fc8d62",
            "Growth"      = "#8da0cb",
            "Respiration" = "#e78ac3")

# Generate Q10 plot
Q10plot <- ggplot() +
  geom_errorbar(data = Q10pdat, aes(x = data, ymin = Q10lwr, ymax = Q10upr),
                width = .10,
                linewidth = .75,
                colour = mycols) +
  geom_point(data = Q10pdat, aes(x = data, y = Q10),
             size = 2.5,
             colour = "black") +
  scale_x_discrete(limits = c("Clearance", "Ingestion", "Growth", "Respiration")) +
  scale_y_break(c(5.3, 4.5), scales = 0.3, ticklabels = c(4.5, 8)) +  # break y-axis
  scale_y_continuous(breaks = c(0, 2, 4, 6), # set-up y-axis explicitly
  limits = c(0, 11.5)) +  # set limits for y-axis
  theme_bw() +
  theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x.top = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        axis.text.y.right = element_blank(),      # Remove right y-axis text
        axis.ticks.y.right = element_blank()) +   # Remove right y-axis ticks
  labs(
    x = expression(bold("Zooplankton rate process")),
    y = bquote(bold("Mass-specific Q"[10])))
Q10plot


# Save it
ggsave("Output/Figure_2_allZoopQ10.png", plot = Q10plot, width = 175, height = 70, unit = "mm",  dpi = 300)
ggsave("Output/Figure_2_allZoopQ10.pdf", plot = Q10plot, width = 175, height = 70, unit = "mm",  dpi = 300)



