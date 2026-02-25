# Calculating sizeGrp clearance rate Q10s
# Josh Hill
# 09/02/2026



  # Here I read in clean feeding rate data
  # Subset the data by target groups
  # Subset by rate type (clearance)
  # Calculate interspecific Q10s
  # Save the Q10 and variance into a dataframe
  # Save Arrhenius plot object for later plotting



# Packages and helpers ----
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
source("R/0_Helpers.R") # use the sizeGrp helpers functions



# Read in the data ----
dat <- readRDS("Data/clear_ingest_data.rds")


  # Check the temperature range for each sizeGrp
  dat %>% 
    select(sizeGrp, temp_C, rate_name) %>% 
    group_by(sizeGrp) %>% 
    filter(rate_name == "ClearanceRate") %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
  


# Exclude data and prep for analysis ----
usedat <- dat %>% 
  filter(rate_name == "ClearanceRate") %>%
  group_by(sizeGrp) %>% 
  filter(n() >= 15, # Exclude sizeGrps that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  select(primRef, sizeGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
  drop_na
  

  
  # Quickly view distribution of Cspecific_rates by sizeGrp
  # I will use these to remove extreme outliers, particularly those that don't make biological sense
  usedat %>%
    ggplot(aes(x = log(Cspecific_rate))) + # log because we will transform the data for analysis
    geom_histogram(bins = 50, fill = "pink", colour = "grey") +
    facet_wrap(~sizeGrp, scales = "free") + 
    theme_bw() +
    labs(
      x = "Clearance rate (Cspecific_rate, log10 scale)",
      y = "Count",
      title = "Distribution of Cspecific_rate for each zooplankton group")
      # distribution looks pretty normal for each group

  usedat %>%
    filter( # exclude values that are not biologically reasonable or are extreme outliers
      (Cspecific_rate < 10000)) %>% 
    ggplot() +
    geom_point(aes(x = temp_C, y = Cspecific_rate)) +
    facet_wrap(~sizeGrp, scales = "free") + 
    theme_bw() +
    labs(
      x = "Temp C",
      y = "Clearance rate")
  
  
# Custom grouping order 
group_order <- c("Mesoplankton", "Macroplankton")

  
# Convert temp to Kelvin and prep data for modelling ----
mdat <- usedat %>% 
  filter( # exclude values that are not biologically reasonable or are extreme outliers
    (sizeGrp == "Mesoplankton"    & Cspecific_rate < 10000) |
      (sizeGrp == "Macroplankton"   & Cspecific_rate < 10000)) %>%
  rename(group = sizeGrp) %>% # rename sizeGrps for Q10 function later...
  filter(!is.na(temp_C)) %>% 
  mutate(temp_K = temp_C + 273.15, 
         x = 1/temp_K,
         y = log(Cspecific_rate)) %>% 
  filter(is.finite(y))

mdat %>% 
  ggplot() +
  geom_point(aes(x = x, y = y, colour = group)) +
  facet_wrap(~group) + 
  theme_bw()



# Fit the model ----
m1 <- glmmTMB(y ~ x + group + (1|primRef) + (1|taxa),  data = mdat) # with primRef and taxa as random intercepts

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim)
  summary(m1)
  # Family: gaussian  ( identity )
  # Formula:          y ~ x + group + (1 | primRef) + (1 | taxa)
  # Data: mdat
  # 
  # AIC       BIC    logLik -2*log(L)  df.resid 
  # 2869.2    2899.1   -1428.6    2857.2      1058 
  # 
  # Random effects:
  #   
  #   Conditional model:
  #   Groups   Name        Variance Std.Dev.
  # primRef  (Intercept) 1.6510   1.2849  
  # taxa     (Intercept) 0.8140   0.9022  
  # Residual             0.6607   0.8128  
  # Number of obs: 1064, groups:  primRef, 72; taxa, 76
  # 
  # Dispersion estimate for gaussian family (sigma^2): 0.661 
  # 
  # Conditional model:
  #   Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)          30.8749     4.2992   7.182 6.89e-13 ***
  #   x                 -7461.3220  1238.9605  -6.022 1.72e-09 ***
  #   groupMesoplankton     0.6848     0.4160   1.646   0.0997 .  
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  r.squaredGLMM(m1)
  # R2m       R2c
  # 0.1846485 0.8276591



# Get Q10 estimates and variance ----
# Run analysis ----
# First split the data by group
groups <- mdat %>%
    group_by(group) %>% # group by sizeGroups
    group_split() # split the data by sizeGrps

# Get results for each group
results <- imap(groups, ~analyse_thermal_sensitivity(m1, .x))

# Get Q10 estimates
Q10dat <- imap_dfr(results, ~{
    df <- groups[[.y]]
    g  <- unique(df$group) # distinct groupings depending on the analysis
  
    tibble(
      Group    = g,
      n_obs    = nrow(df),                         # number of observations
      temp_min = min(df$temp_C, na.rm = TRUE),     # min temperature
      temp_max = max(df$temp_C, na.rm = TRUE),     # max temperature
      Q10      = .x$Q10$Q10,                       # Q10 estimate
      Q10lwr   = .x$Q10$CI_lower,                  # Lower limit
      Q10upr   = .x$Q10$CI_upper                   # Upper limit
  )
}) %>%
  mutate(Group = factor(Group))

# Show me the results
Q10dat
# Group         n_obs temp_min temp_max   Q10 Q10lwr Q10upr
# Macroplankton   805     -1         29  2.64  0.967   7.22
# Mesoplankton    259     -1.7       30  2.66  1.01    6.99

# Save it
# saveRDS(Q10dat, "Data/sizeGrp_analysis/sizeGrp_clearance_Q10_data.rds")



# Prepare data for plotting ----
# Create arrhenius objects for each group
# Combine all group data
all_data_combined <- imap_dfr(groups, ~{
  .x %>% mutate(group = unique(.x$group))
})

all_pred_combined <- imap_dfr(results, ~{
  df <- groups[[.y]]
  g  <- unique(df$group)
  .x$predictions %>% mutate(group = g)
})

# Create Q10 annotations dataframe
Q10_label <- imap_dfr(results, ~{
  df <- groups[[.y]]
  g  <- unique(df$group)
  
  tibble(
    group = g,
    Q10_label = paste0(
      "Q[10] == ", round(.x$Q10$Q10, 2),
      " * ' (95% CI: ", round(.x$Q10$CI_lower, 2), " - ",
      round(.x$Q10$CI_upper, 2), "); n = ", .x$Q10$n_obs, "'"
    ),
    x = 0.0037,
    y = 12
  )
})

# Create the arrhenius object (standard format for your plotting function)
arrhenius_obj <- list(
  all_data = all_data_combined,
  all_pred = all_pred_combined,
  Q10_label = Q10_label,
  Q10dat = Q10dat  # Include summary table
)


# Save object for plotting
# saveRDS(arrhenius_obj, "Data/sizeGrp_analysis/sizeGrp_clearance_arrhenius_data.rds")


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
            parse = TRUE,
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



