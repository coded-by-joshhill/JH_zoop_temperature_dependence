# Calculating funcGrps Q10s
# Josh Hill
# 03/11/2026



  # Here I read in all cleaned data
  # Combine data and subset by funcGrp
  # Fit glmm with random effects to test if the effect of temp on biological rate process
  # Calculate Q10s from the model slopes



# Packages and helpers ----
library(tidyverse)
library(glmmTMB) # for modelling
library(DHARMa) # for diagnostics
library(emmeans) # Estimated marginal means
library(performance)
library(patchwork)
theme_set(new = theme_bw())



# Read in the data and prep for analysis ----

# Feeding data
dat <- readRDS("Data/clear_ingest_data.rds")


# Clearance data
cleardat <- dat %>% 
  filter(rate_name == "ClearanceRate") %>% # Filter for clearance rate
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Clearance"))


# Ingestion data
ingdat <- dat %>% 
  filter(rate_name == "IngestionRate") %>% # Filter for clearance rate
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Ingestion"))


# Growth data
grwdat <- readRDS("Data/grwth_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Growth"))


# Respiration data
respdat <- readRDS("Data/resp_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Respiration"))


# Custom grouping order
group_order <- c("Crustaceans", 
                 "GelPreds", 
                 "GelFilter")


# Combine them into one dataframe
usedat <- rbind(cleardat, ingdat, grwdat, respdat) %>% 
  filter_out(funcGrp == "OTHER" | is.na(temp_C)) %>%  # filter out the functional group "OTHER" and remove any NAs in temp_C
  mutate(funcGrp = fct_relevel(funcGrp, group_order)) # reorder funcGrp
  

# Check the temperature range
usedat %>% 
  group_by(funcGrp, rate_name) %>% 
  select(temp_C) %>% 
  summarise( 
    temp_range = paste0(min(temp_C), "-", max(temp_C))) %>% 
  arrange(rate_name, funcGrp)


# Quickly view distribution of raw Cspecific_rates
# I will use these to remove extreme outliers, particularly those that don't make biological sense
usedat %>%
  ggplot() +
  geom_point(aes(x = temp_C, y = Cspecific_rate)) +
  theme_bw() +
  facet_wrap(~ rate_name, scales = "free") + 
  labs(
    x = "Temp C",
    y = "Clearance rate")


usedat %>%
  ggplot(aes(x = log(Cspecific_rate))) + # log because we will transform the data for analysis
  geom_histogram(bins = 50, fill = "pink", colour = "grey") +
  theme_bw() +
  facet_wrap(~rate_name, scales = "free") +
  labs(
    x = "Mass-specific rate (Cspecific_rate, log scale)",
    y = "Count",
    title = "Distribution of Cspecific_rate across all zooplankton")
    # distribution looks pretty normal but ingestion and clearance are slightly left skewed



# Tidy up data and prep data for modelling ----
mdat <- usedat %>% 
  group_by(funcGrp) %>% 
  filter(n() >= 15, # Exclude funcGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  mutate(ln_Cspecific_rate = log(Cspecific_rate)) # log transform mass-specific rate


# Quick look
mdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = funcGrp)) +
  facet_wrap(~ rate_name, scale = "free") +
  theme_bw()
  # Looks pretty tidy. Some clear relationships here

summary(mdat)


# My main question here is...
  # How does temperature dependence vary across zooplankton func groups for each rate? AND
  # How does temperature dependence vary across rate processes?

# I am also mainly interested in the interactions between at least temp:rate and temp:group



# Fit the models ----

# A complex model with 3 way interactions for temp, funcGrp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * funcGrp + # 3-way interaction
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim) # looks pretty good
  summary(m1)
  # Kind of difficult to interpret as before...

  
# A simpler model without the 3-way interactions, just 2-way interactions
m2 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                temp_C * funcGrp +  # between temp and different funcGrps for all rates
                rate_name * funcGrp + # between rates and funcGrp 
                (temp_C | primRef) + (1 | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m2)
  plot(sim) # Looks fine but seems less scattered compared to m1
  summary(m2)

  
# Compare models
performance::compare_performance(m1, m2)
# Models are not all mutually nested...we'll just treat AIC/BIC as descriptive...


# Likelihood ratios test of the models
# Refit with ML for valid test on fixed/dispersion and random effect structures
m1_m1 <- update(m1, REML = FALSE)
m2_m2 <- update(m2, REML = FALSE)

anova(m1_m1, m2_m2) # test if the three-way interaction is better than the two-way interaction
  # m1 with 3-way interaction is slightly better

# we will progress with m1
summary(m1)


# Extract slopes using and calculate Q10 for each funcGrp ----
funcGrp_slopes <- emtrends(m1, ~ rate_name * funcGrp, var = "temp_C")
summary(funcGrp_slopes, infer = TRUE) # test whether each slope is different from zero for each zoopGrp
pairs(funcGrp_slopes, by = "rate_name") # pairwise test whether slopes differ significantly across rate types and grps
# yes, some slightly significant differences for:
  # Clearance - crust - gelPreds
  # Clearance - crust - gelFilter
  # Growth - crust - gelPreds

# Get n
n_obs <- mdat %>%
  count(rate_name, funcGrp)

# Get Q10
funcGrp_slopes_Q10 <- as.data.frame(funcGrp_slopes) %>% 
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  left_join(n_obs, by = c("rate_name", "funcGrp")) %>% 
  arrange(rate_name, funcGrp)
funcGrp_slopes_Q10



# Prep for plotting ----

# Define group colours
grp_cols <- c("Crustaceans" = "#0077BB",
              "GelPreds" = "#EE7733",
              "GelFilter" = "#009988")


# Set min and max temp
minTempC <- min(mdat$temp_C)
maxTempC <- max(mdat$temp_C)


# FUNCTIOIN TO PLOT MODEL ----
PlotLMM = function(model){
  temp_seq <- seq(minTempC, maxTempC, length.out = 100)
  
  # Build newdata grid manually
  newdat <- expand.grid(
    temp_C    = temp_seq,
    funcGrp   = unique(mdat$funcGrp), # levels of funcGrp
    rate_name = unique(mdat$rate_name)) # levels of rates
  

  # Zooplankton-level (population) predictions
  pop_preds <- newdat
  pred <- predict(model, newdata = newdat, se.fit = TRUE,
                  type = "response",
                  re.form = NA)
  pop_preds$estimate <- pred$fit
  pop_preds$conf.low  <- pred$fit - 1.96 * pred$se.fit
  pop_preds$conf.high <- pred$fit + 1.96 * pred$se.fit
  

    # Plot it up...
  ggplot() +
    geom_ribbon(data = pop_preds,
                aes(x = temp_C, ymin = conf.low, ymax = conf.high, fill = funcGrp),
                alpha = 0.25) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate, colour = funcGrp),
              linewidth = 1) +
    geom_point(data = mdat,
               aes(x = temp_C, y = ln_Cspecific_rate, colour = funcGrp),
               alpha = 0.2) +
    facet_wrap(~rate_name, scales = "free",
               labeller = as_labeller(c(
                 "Clearance"   = "bold(Clearance~rate~(ml~mgC^-1~h^-1))",
                 "Ingestion"   = "bold(Ingestion~rate~(mgC~mgC^-1~h^-1))",
                 "Growth"      = "bold(Growth~rate~(mgC~mgC^-1~h^-1))",
                 "Respiration" = "bold(Respiration~rate~(µlO[2]~mgC^-1~h^-1))"
               ), 
               label_parsed))+ 
    coord_cartesian(xlim = c(-2, 32)) +
    scale_fill_manual(values = grp_cols, labels = c("Crustaceans" , "Gelatinious predators", "Gelatinious filter-feeders")) +
    scale_colour_manual(values = grp_cols, labels = c("Crustaceans" , "Gelatinious predators", "Gelatinious filter-feeders")) +
    labs(x = "Temp (°C)",
         y = "ln(Carbon-mass specific rate)",
         fill = "Functional group",
         colour = "Functional group") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "whitesmoke", colour = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.position = "top"
    )
}


# Plot it
tempPlot <- PlotLMM(m1)
tempPlot


# Plot funcGrp Q10s
funcGQ10plot <- ggplot() +
  geom_errorbar(data = funcGrp_slopes_Q10, 
                aes(x = funcGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = funcGrp),
                width = .09,
                linewidth = 1) +
  geom_point(data = funcGrp_slopes_Q10, aes(x = funcGrp, y = Q10),
             size = 3,
             colour = "black") +
  geom_text(data = funcGrp_slopes_Q10,
            aes(x = funcGrp, y = Q10, label = sprintf("%.2f", Q10)),
            nudge_x = 0.3,
            nudge_y = -0.0) +
  geom_text(data = funcGrp_slopes_Q10,
            aes(x = funcGrp, y = Q10_lwr, label = paste0("n = ", n)),
            nudge_y = -0.3,   
            size = 3.5, colour = "grey40") +
  facet_wrap(~rate_name, scales = "free",
             labeller = as_labeller(c(
               "Clearance"   = "bold(Clearance~rate~(ml~mgC^-1~h^-1))",
               "Ingestion"   = "bold(Ingestion~rate~(mgC~mgC^-1~h^-1))",
               "Growth"      = "bold(Growth~rate~(mgC~mgC^-1~h^-1))",
               "Respiration" = "bold(Respiration~rate~(µlO[2]~mgC^-1~h^-1))"
             ), 
             label_parsed)) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "Functional group",
       y = bquote(bold("Carbon-mass specific Q"[10])),
       colour = "Functional group") +
  theme(
    strip.background = element_rect(fill = "whitesmoke", colour = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )
funcGQ10plot

library(ggbreak) # to break y-axis on massive Q10 variance
funcG10plot_break <- funcGQ10plot +
  geom_text(data = filter(funcGrp_slopes_Q10, 
                          (rate_name == "Ingestion" & Q10_upr > 6) |
                          (rate_name == "Respiration" & Q10_upr > 6) |
                          (rate_name == "Growth" & Q10_upr > 6) |
                          (rate_name == "Clearance" & Q10_upr > 6)),
            aes(x = funcGrp, y = 5.85, label = paste0("↑ ", round(Q10_upr, 0))),
            colour = "grey40", size = 4,
            nudge_x = 0.3) +
  coord_cartesian(ylim = c(-1, 6))
funcG10plot_break

tempPlot/funcG10plot_break

# Save it
ggsave("Output/Figure4/Figure4_tempPlot.pdf", tempPlot, width = 160, height = 150, units = "mm", dpi = 300)
ggsave("Output/Figure4/Figure4_Q10Plot.pdf", funcG10plot_break, width = 160, height = 120, units = "mm", dpi = 300)

