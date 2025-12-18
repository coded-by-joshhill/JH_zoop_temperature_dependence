# Calculating ingestion rate Q10s and plotting
# Josh Hill
# 1/12/25



  # Here I read in clean feeding rate data
  # Subset the data by target groups
  # Subset by rate type (ingestion)
  # Calculate interspecific Q10s
  # Save those Q10s and variance into a dataframe
  # Plot the Q10s and C-specific rates for available groups
  # Generate Arrhenius plots showing the relationship between C-specific rates and temperature (for supp materials)



# Packages and helpers ----
library(tidyverse)
library(parallel)
library(DHARMa) 
library(ggeffects)
library(patchwork)
source("R/0_Helpers.R")



# Read in the data ----
dat <- readRDS("Data/clear_ingest_data.rds") %>% 
  mutate(zoopGrp = as.factor(zoopGrp)) # set to factor


  # Check the temperature range for each zoopGrp
  dat %>% 
    select(zoopGrp, temp_C, rate_name) %>% 
    group_by(zoopGrp) %>% 
    filter(rate_name == "IngestionRate") %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
  


# Exclude data and prep for analysis ----
usedat <- dat %>% 
    filter(rate_name == "IngestionRate") %>%
    group_by(zoopGrp) %>% 
    filter(n() >= 20, # Exclude zoopGrps that don't have suitable data or temp ranges
           max(temp_C) - min(temp_C) >= 5) %>% 
    ungroup() %>% 
    select(primRef, zoopGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
    drop_na(Cspecific_rate)
  

  
  # Quickly view distribution of Cspecific_rates by zoopGrp
  usedat %>%
    ggplot(aes(x = Cspecific_rate)) +
    geom_histogram(bins = 50, fill = "pink", color = "grey") +
    scale_x_log10() + # use log 10 to see extreme outliers
    facet_wrap(~zoopGrp, scales = "free") + 
    theme_bw() +
    labs(
      x = "Clearance rate (Cspecific_rate, log10 scale)",
      y = "Count",
      title = "Distribution of Cspecific_rate for each zooplankton group")
      # I will use this to remove extreme outliers, particularly those that don't make biological sense

  
# Custom grouping order 
group_order <- c("Euphausiacea", "Copepoda", "Ctenophora", 
                 "Cnidaria", "Thaliacea", "Appendicularia")

  
  
# Convert temp to Kelvin and prep data for modelling ----
mdat <- usedat %>% 
    filter( # exclude values that are not biologically reasonable
      (zoopGrp == "Copepoda"       & Cspecific_rate < 0.2) |
      (zoopGrp == "Euphausiacea"   & Cspecific_rate < 0.06) |
      (zoopGrp == "Thaliacea"      & Cspecific_rate < 0.4)) %>%
    filter(!is.na(temp_C)) %>% 
    mutate(temp_K = temp_C + 273.15,
           x = 1/temp_K,
           y = log(Cspecific_rate),
           zoopGrp = fct_relevel(zoopGrp, group_order)) %>% 
    filter(is.finite(y))

mdat %>% 
  ggplot() +
  geom_point(aes(x = x, y = y, colour = zoopGrp)) +
  theme_bw()



# Fit the bootstrap models using parallel processing ----

# Check no. cores
detectCores()
ncores <- detectCores() - 2

# Specify cores for parallel processing
mirai::daemons(ncores)

# Using glmmTMB fit (y ~ x * zoopGrp + (1|primRef) + (1|taxa)).... this takes a while
boot_models <- map(1:9999, in_parallel(\(x) boot_Q10(mdat),
                                      # Specify locally-defined functions
                                      boot_Q10 = boot_Q10,
                                      glmmTMB = glmmTMB,
                                      mdat = mdat)) 

# Get Q10 estimates from bootstrap models... this takes even longer :)
Q10_estimates <- boot_models |>
  map_dfr(in_parallel(\(df) get_Q10s(df), 
                      # Specify locally-defined function
                      get_Q10s = get_Q10s))


  # Get confidence intervals
  Q10pdat <- Q10_estimates %>% 
    group_by(Group) %>% 
    summarise(median = median(Q10, na.rm = TRUE),
              lower_CI = quantile(Q10, 0.025),
              upper_CI = quantile(Q10, 0.975)) %>% 
    ungroup()

# Tear down existing daemons
mirai::daemons(0)

Q10pdat

# Save as RDS
# saveRDS(Q10_estimates, "Data/Q10_estimates_ingestion.rds") # estimates
# saveRDS(mdat, "Data/ingestion_mdat.rds") # save modelling dataframe for plotting
# saveRDS(Q10pdat, "Data/Q10_summary_ingestion.rds") # median and confidence intervals



# Create Arrhenius plots ----
# Extract coefficients from all bootstrap models
coefs_list <- map(boot_models, ~fixef(.x)$cond)

# Get all zooplankton groups
groups <- levels(mdat$zoopGrp)

# Extract results using get_results function
results <- get_results(mdat, groups, coefs_list, Q10pdat)

# Plot it up
ingestion_plot <- arrhenius_plot(
  mdat = mdat,
  rate_col = "Cspecific_rate",
  boot_models = boot_models,
  Q10pdat = Q10pdat,
  group_order = group_order,
  x_limits = c(0.0037, 0.0033)) +
  labs(y = expression(bold("Ingestion rate (mgC mgC"^-1*" h"^-1*")")))

ingestion_plot

# ggsave("Output/Figure_Supp4.png", plot = ingestion_plot, width = 12, height = 4.5)
