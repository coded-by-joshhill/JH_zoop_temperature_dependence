# Calculating clearance rate Q10s and plotting
# Josh Hill
# 19/11/25



  # Here I read in clean feeding rate data
  # Subset the data by target groups
  # Subset by rate type (clearance)
  # Calculate interspecific Q10s
  # Save those Q10s and variance into a dataframe
  # Plot the C-specific rates and Q10s
  # Plot the relationship between C-specific rates and temperature - these will be used as supplementary materials figures



# Packages and helpers ----
library(tidyverse)
library(DHARMa) 
library(glmmTMB)
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
    filter(rate_name == "ClearanceRate") %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
  


# Exclude data and prep for analysis ----
usedat <- dat %>% 
    filter(rate_name == "ClearanceRate") %>%
    group_by(zoopGrp) %>% 
    filter(n() >= 20, # Exclude and zoopGrps that don't have suitable data or temp ranges
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

  
  
# Convert temp to Kelvin and prep data for model ----
mdat <- usedat %>% 
    filter(
      (zoopGrp == "Appendicularia" & Cspecific_rate < 15000) |
      (zoopGrp == "Cnidaria"       & Cspecific_rate < 5000) |
      (zoopGrp == "Copepoda"       & Cspecific_rate < 3000) |
      (zoopGrp == "Ctenophora"     & Cspecific_rate < 10000) |
      (zoopGrp == "Euphausiacea"   & Cspecific_rate < 50) |
      (zoopGrp == "Thaliacea"      & Cspecific_rate < 5000)) %>%
    filter(!is.na(temp_C)) %>% 
    mutate(temp_K = temp_C + 273.15,
           x = 1/temp_K,
           y = log(Cspecific_rate)) %>% 
    filter(is.finite(y))

mdat %>% 
  ggplot() +
  geom_point(aes(x = x, y = y, colour = zoopGrp)) +
  theme_bw()




# Fit the model ----
m1 <- glmmTMB(y ~ x * zoopGrp + (1|primRef) + (1|taxa),  data = mdat) 
  summary(m1)
  
  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim)
  testOutliers(sim) 


  
# Run analysis ----
results <- mdat %>%
  group_by(zoopGrp) %>%
  group_split() %>%
  lapply(function(df) analyse_thermal_sensitivity(m1, df))

# Name the list by group
names(results) <- unique(mdat$zoopGrp)

# Build Q10 summary table with diagnostics
clearanceQ10s <- lapply(names(results), function(g) {
  df <- mdat %>% filter(zoopGrp == g)
  
  data.frame(
    Group    = g,
    n_obs    = nrow(df),                     # number of observations
    temp_min = min(df$temp_C, na.rm = TRUE), # min temperature in group
    temp_max = max(df$temp_C, na.rm = TRUE), # max temperature in group
    Q10      = results[[g]]$Q10$Q10,
    Q10lwr   = results[[g]]$Q10$CI_lower,
    Q10upr   = results[[g]]$Q10$CI_upper
  )
}) %>% bind_rows() %>% 
  mutate(Group = factor(Group))

  
# Show me the results ----
clearanceQ10s



# Plot it up ---- 
# Violin plot of log(Cspecific_rate)
ratePlot <- mdat %>%
  ggplot(aes(x = zoopGrp, y = log(Cspecific_rate))) +
  geom_violin(fill = "skyblue", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +  # optional: show individual points
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x.top = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12)) +
  labs(
    x = expression(bold("Zooplankton group")),
    y = expression(bold("log10 Clearance rate (ml mgC"^-1*" h"^-1*")"))) 


# Dot plot of Q10
Q10plot <- clearanceQ10s %>%
  ggplot(aes(x = Group, y = Q10)) +
  geom_errorbar(aes(ymin = Q10lwr, ymax = Q10upr), width = 0.2, color = "skyblue") +
  geom_point(size = 3, color = "darkblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x.top = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12)) +  
  labs(
    x = expression(bold("Zooplankton group")),
    y = expression(bold("Q10")))

(figure2 <- ratePlot + Q10plot +
  plot_layout(axis_titles = "collect_x"))

# ggsave("Output/Figure_2.png", plot = figure2, width = 12, height = 8)


# Produce Arrhenius plots for supplementary material ----
# Prepare combined data for faceting
plot_data <- lapply(names(results), function(g) {
  df_data <- mdat %>% filter(zoopGrp == g)
  df_pred <- results[[g]]$predictions
  df_Q10 <- results[[g]]$Q10
  
  list(
    data = df_data %>% mutate(zoopGrp = g),
    pred = df_pred %>% mutate(zoopGrp = g),
    Q10_label = paste0("Q10 = ", round(df_Q10$Q10, 2),
                       " (95% CI: ", round(df_Q10$CI_lower, 2), " - ",
                       round(df_Q10$CI_upper, 2), ")")
  )
})

# Combine all data
all_data <- bind_rows(lapply(plot_data, function(x) x$data))
all_pred <- bind_rows(lapply(plot_data, function(x) x$pred))
all_labels <- data.frame(
  zoopGrp = names(results),
  label = sapply(plot_data, function(x) x$Q10_label),
  x = 0.0037,
  y = Inf  # This will place it at the top
)


# Create plot
arrhenius_plots <- ggplot() +
  # raw points
  geom_point(data = all_data, 
             aes(x = x, y = Cspecific_rate), 
             alpha = 0.3, size = 2) +
  # prediction ribbon
  geom_ribbon(data = all_pred,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "lightgrey", alpha = 0.3) +
  # predicted line
  geom_line(data = all_pred, 
            aes(x = x, y = predicted), 
            colour = "darkblue", linewidth = 1, linetype = "dashed") +
  # Q10 annotations
  geom_text(data = all_labels, 
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1.5,
            size = 4, colour = "darkblue") +
  # reverse x-axis to account for Arrhenius styling
  scale_x_continuous(name = expression(bold("1 / Temp (K"^-1*")")),
                     trans = "reverse",
                     sec.axis = sec_axis(~1/. - 273.15, name = "Temp (°C)")) +
  labs(y = expression(bold("Clearance rate (ml mgC"^-1*" h"^-1*")"))) +
  facet_wrap(~ zoopGrp, scales = "free_y", ncol = 3) +
  coord_cartesian(xlim = c(0.0037, 0.0033)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x.top = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

arrhenius_plots

# ggsave("Output/Figure_Supp1.png", plot = arrhenius_plots, width = 12, height = 8)
