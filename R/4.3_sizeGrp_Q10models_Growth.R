# Calculating sizeGrp growth rate Q10s and plotting
# Josh Hill
# 09/02/2026



  # Here I read in clean growth rate data
  # Subset the data by target groups
  # Calculate interspecific Q10s
  # Save those Q10s and variance into a dataframe
  # Save Arrhenius plot object for later plotting



# Packages and helpers ----
library(tidyverse)
library(parallel)
source("R/0.1_sizeGrp_Helpers.R")



# Read in the data ----
dat <- readRDS("Data/grwth_dat.rds")


  # Check the temperature range for each sizeGrp
  dat %>% 
    select(sizeGrp, temp_C, rate_name) %>% 
    group_by(sizeGrp) %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
  


# Exclude data and prep for analysis ----
usedat <- dat %>% 
    group_by(sizeGrp) %>% 
    filter(n() >= 15, # Exclude sizeGrps that don't have suitable data or temp ranges
           max(temp_C) - min(temp_C) >= 5) %>% 
    ungroup() %>% 
    select(primRef, sizeGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
    drop_na(Cspecific_rate)
  

  
  # Quickly view distribution of Cspecific_rates by sizeGrp
  usedat %>%
    ggplot(aes(x = Cspecific_rate)) +
    geom_histogram(bins = 50, fill = "pink", color = "grey") +
    scale_x_log10() + # use log 10 to see extreme outliers
    facet_wrap(~sizeGrp, scales = "free") + 
    theme_bw() +
    labs(
      x = "Clearance rate (Cspecific_rate, log10 scale)",
      y = "Count",
      title = "Distribution of Cspecific_rate for each zooplankton group")
      # I will use this to remove extreme outliers, particularly those that don't make biological sense

  
# Custom grouping order 
group_order <- c("Mesoplankton", "Macroplankton")
  
  
# Convert temp to Kelvin and prep data for modelling ----
mdat <- usedat %>% 
    # filter( # exclude values that are not biologically reasonable and remove extreme outliers
    #   (sizeGrp == "Mesoplankton"    & Cspecific_rate < 0.05) |
    #   (sizeGrp == "Macroplankton"   & Cspecific_rate < 0.05)) %>%
    filter(!is.na(temp_C)) %>% 
    mutate(temp_K = temp_C + 273.15,
           x = 1/temp_K,
           y = log(Cspecific_rate),
           sizeGrp = fct_relevel(sizeGrp, group_order)) %>% 
    filter(is.finite(y))

mdat %>% 
  ggplot() +
  geom_point(aes(x = x, y = y, colour = sizeGrp)) +
  theme_bw()



# Fit the bootstrap models using parallel processing ----

# Check no. cores
detectCores()
ncores <- detectCores() - 2

# Specify cores for parallel processing
mirai::daemons(ncores)

# Using glmmTMB fit (y ~ x * sizeGrp + (1|primRef) + (1|taxa)).... this takes a while
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
# saveRDS(Q10_estimates, "Data/sizeGrp_Q10_estimates_growth.rds") # estimates
# saveRDS(mdat, "Data/sizeGrp_growth_mdat.rds") # save modelling dataframe for plotting
# saveRDS(Q10pdat, "Data/sizeGrp_Q10_summary_growth.rds") # median and confidence intervals



# Create Arrhenius plots data ----
# Extract coefficients from all bootstrap models
coefs_list <- map(boot_models, ~fixef(.x)$cond)

# Get all zooplankton groups
groups <- levels(mdat$sizeGrp)

# Extract results using get_results function
results <- get_results(mdat, groups, coefs_list, Q10pdat)

# Arrhenius object 
arrhenius_obj <- list(
  mdat    = mdat[, c("sizeGrp", "x", "Cspecific_rate")],
  results = results,
  group_order = group_order)

# Save object for plotting
# saveRDS(arrhenius_obj, "Data/sizeGrp_ingestion_arrhenius_plot_data.rds")


# Generate the plot
growth_plot <- arrhenius_plot(
  mdat = arrhenius_obj$mdat,
  rate_col = "Cspecific_rate",
  results = arrhenius_obj$results,
  group_order = arrhenius_obj$group_order,
  x_limits = c(0.0037, 0.0033)) +
  labs(y = expression(bold("Growth rate (mgC mgC"^-1*" h"^-1*")")))

growth_plot
