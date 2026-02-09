# Calculating allZoop clearance rate Q10s
# Josh Hill
# 09/02/2026



  # Here I read in clean feeding rate data
  # Subset the data by target groups
  # Subset by rate type (clearance)
  # Calculate interspecific Q10s
  # Save those Q10s and variance into a dataframe
  # Save Arrhenius plot object for later plotting



# Packages and helpers ----
library(tidyverse)
library(parallel)
library(DHARMa)
source("R/0_Helpers.R")



# Read in the data ----
dat <- readRDS("Data/clear_ingest_data.rds")


  # Check the temperature range for each sizeGrp
  dat %>% 
    select(temp_C, rate_name) %>% 
    filter(rate_name == "ClearanceRate") %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
  


# Exclude data and prep for analysis ----
usedat <- dat %>% 
    filter(rate_name == "ClearanceRate") %>%
    ungroup() %>% 
    select(primRef, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
    drop_na(Cspecific_rate)
  

  
  # Quickly view distribution of Cspecific_rates by sizeGrp
  usedat %>%
    ggplot(aes(x = Cspecific_rate)) +
    geom_histogram(bins = 50, fill = "pink", color = "grey") +
    scale_x_log10() + # use log 10 to see extreme outliers
    theme_bw() +
    labs(
      x = "Clearance rate (Cspecific_rate, log10 scale)",
      y = "Count",
      title = "Distribution of Cspecific_rate across all zooplankton")
      # I will use this to remove extreme outliers, particularly those that don't make biological sense

  
  
# Convert temp to Kelvin and prep data for modelling ----
mdat <- usedat %>% 
    filter( # exclude values that are not biologically reasonable or are extreme outliers
      (Cspecific_rate < 15000)) %>% 
    filter(!is.na(temp_C)) %>% 
    mutate(temp_K = temp_C + 273.15,
           x = 1/temp_K,
           y = log(Cspecific_rate)) %>% 
    filter(is.finite(y))

  
mdat %>% 
  ggplot() +
  geom_point(aes(x = x, y = y)) +
  theme_bw()



# Fit the model ----
m1 <- glmmTMB(y ~ x + (1|primRef) + (1|taxa),  data = mdat) # fit the model

  
  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim)
  testOutliers(sim) 
  summary(m1)



# Get Q10 estimates and variance ----
results <- analyse_thermal_sensitivity(m1, mdat)


# Build Q10 summary table with diagnostics
Q10dat <- data.frame(
    n_obs    = nrow(mdat),                     # number of observations
    temp_min = min(mdat$temp_C, na.rm = TRUE), # min temperature in group
    temp_max = max(mdat$temp_C, na.rm = TRUE), # max temperature in group
    Q10      = results$Q10$Q10,                # Q10 estimate
    Q10lwr   = results$Q10$CI_lower,           # Upper Q10
    Q10upr   = results$Q10$CI_upper)           # Lower Q10

# Show me the results
Q10dat

# saveRDS(Q10dat, "Data/allZoop_analysis/allZoop_clearance_Q10_data.rds")



# Prepare data for plotting ----
all_data <- mdat  # raw observations
all_pred <- results$predictions  # model predictions


# Q10 annotation
Q10_label <- paste0("Q10 = ", round(results$Q10$Q10, 2),
                    " (95% CI: ", round(results$Q10$CI_lower, 2), " - ",
                    round(results$Q10$CI_upper, 2), ")")

# Arrhenius object 
arrhenius_obj <- list(
  all_data    = all_data,
  all_pred    = all_pred,
  Q10_label   = Q10_label)

# Save object for plotting
saveRDS(arrhenius_obj, "Data/allZoop_analysis/allZoop_clearance_arrhenius_data.rds")

# Create Arrhenius plot (in log space because our rates vary in orders of magnitude)
arrhenius_plot <- ggplot() +
  # raw points 
  geom_point(data = arrhenius_obj$all_data, 
             aes(x = x, y = log(Cspecific_rate)), 
             alpha = 0.3, size = 2) +
  # prediction ribbon
  geom_ribbon(data = arrhenius_obj$all_pred,
              aes(x = x, ymin = conf.low_log, ymax = conf.high_log),
              fill = "lightgrey", alpha = 0.3) +
  # predicted line
  geom_line(data = arrhenius_obj$all_pred, 
            aes(x = x, y = predicted_log), 
            colour = "darkblue", linewidth = 1, linetype = "dashed") +
  # Q10 annotation (log scale, optional: you can annotate raw Q10)
  geom_text(aes(x = 0.0037, y = 12, 
                label = arrhenius_obj$Q10_label),
            hjust = 0, vjust = 1.5,
            size = 4, colour = "darkblue") +
  # reverse x-axis for Arrhenius
  scale_x_continuous(name = expression(bold("1 / Temp (K"^-1*")")),
                     trans = "reverse",
                     sec.axis = sec_axis(~1/. - 273.15, name = "Temp (°C)")) +
  labs(y = expression(bold("ln(Clearance rate) (ml mgC"^-1*" h"^-1*")"))) +
  coord_cartesian(xlim = c(0.0037, 0.0033)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x.top = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 5, r = 15, b = 5, l = 5))

arrhenius_plot



