# Calculating overall rate Q10s
# Josh Hill
# 25/02/2026



  # Here I read in all cleaned data
  # Combine data
  # Fit glmm with random effects to test if the effect of temp on biological rate process
  # Calculate Q10s from the model slopes



# Packages and helpers ----
library(tidyverse)
library(glmmTMB) # for modelling
library(DHARMa) # for diagnostics
library(MuMIn) # for Rsqr
library(emmeans) # Estimated marginal means
library(performance)
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


# Combine them into one dataframe
usedat <- rbind(cleardat, ingdat, grwdat, respdat)


# Check the temperature range
usedat %>% 
  select(temp_C) %>% 
  summarise( 
    temp_range = paste0(min(temp_C), "-", max(temp_C)))
    # -1.8-31 degC


# Quickly view distribution of raw Cspecific_rates
# I will use these to remove extreme outliers, particularly those that don't make biological sense
usedat %>%
  # exclude values that are not biologically reasonable or are extreme outliers
  filter(
    (rate_name == "Clearance" & Cspecific_rate < 15000) |
    (rate_name == "Ingestion" & Cspecific_rate < 0.15) |
    (rate_name == "Growth" & Cspecific_rate < 0.075) |
    (rate_name == "Respiration" & Cspecific_rate < 60)
  ) %>%
  ggplot() +
  geom_point(aes(x = temp_C, y = Cspecific_rate)) +
  theme_bw() +
  facet_wrap(~ rate_name, scales = "free") + 
  labs(
    x = "Temp C",
    y = "Clearance rate")


usedat %>%
  # exclude values that are not biologically reasonable or are extreme outliers
  filter(
    (rate_name == "Clearance" & Cspecific_rate < 15000) |
      (rate_name == "Ingestion" & Cspecific_rate < 0.15) |
      (rate_name == "Growth" & Cspecific_rate < 0.075) |
      (rate_name == "Respiration" & Cspecific_rate < 60)
  ) %>%
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
  # exclude values that are not biologically reasonable or are extreme outliers
  filter(
    (rate_name == "Clearance" & Cspecific_rate < 15000) |
      (rate_name == "Ingestion" & Cspecific_rate < 0.15) |
      (rate_name == "Growth" & Cspecific_rate < 0.075) |
      (rate_name == "Respiration" & Cspecific_rate < 60)
  ) %>%
  mutate(ln_Cspecific_rate = log(Cspecific_rate)) # log transform mass-specific rate


# Quick look
mdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate)) +
  facet_wrap(~ rate_name, scale = "free") +
  theme_bw()
  # Looks pretty tidy. Some clear relationships here

summary(mdat)


# My main question here is...
  # How does temperature dependence vary across each rate for all zooplankton?



# Fit the models ----

# A complex model with 2 way interactions for temp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ temp_C * rate_name + 
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim) # looks pretty good
  summary(m1)
  r.squaredGLMM(m1)
  # Kind of difficult to interpret I'll fit a simpler model to tease this apart

  
# A simpler model without interactions
m2 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C + rate_name + 
                (temp_C | primRef) + (1 | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m2)
  plot(sim) # Looks fine but seems less scattered compared to m1
  summary(m2)
  r.squaredGLMM(m2)


# m1 but without random slopes, just intercepts
m3 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                (1 | primRef) + (1 | taxa), # with primRef and taxa as random intercepts
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m3)
  plot(sim) # Looks fine but seems less scattered compared to m1
  summary(m3)
  r.squaredGLMM(m3)


performance::compare_performance(m1, m2, m3)


# Likelihood ratios test of the models
anova(m1, m2, m3)
# likelihood ratios test shows m1 has significantly more explanatory power 
# BIC is also slightly better
# I will use m1 on the basis of the chisqr test and AIC and BIC

summary(m1)



# Prep for plotting ----

# Set min and max temp
minTempC <- min(mdat$temp_C)
maxTempC <- max(mdat$temp_C)


# FUNCTIOIN TO PLOT MODEL ----
PlotLMM = function(model){
  temp_seq <- seq(minTempC, maxTempC, length.out = 100)
  
  # Build newdata grid manually
  newdat <- expand.grid(
    temp_C    = temp_seq,
    rate_name = unique(mdat$rate_name) # levels of rates
  )
  
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
    geom_point(data = mdat,
               aes(x = temp_C, y = ln_Cspecific_rate), colour = "grey27",
               alpha = 0.2) +
    geom_ribbon(data = pop_preds,
                aes(x = temp_C, ymin = conf.low, ymax = conf.high), fill = "midnightblue",
                alpha = 0.25) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate), colour = "midnightblue",
              linewidth = 1) +

    facet_wrap(~rate_name, scales = "free",
               labeller = as_labeller(c(
                 "Clearance"   = "bold(Maximum~clearance~rate~(ml~mgC^-1~h^-1))",
                 "Ingestion"   = "bold(Maximum~ingestion~rate~(mgC~mgC^-1~h^-1))",
                 "Growth"      = "bold(Growth~rate~(mgC~mgC^-1~h^-1))",
                 "Respiration" = "bold(Respiration~rate~(µlO[2]~mgC^-1~h^-1))"
               ), 
               label_parsed))+ 
    coord_cartesian(xlim = c(-2, 32)) +
    labs(x = "Temp (°C)",
         y = "ln(Carbon-mass specific rate)") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "whitesmoke", colour = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold")
    )
}



# Plot it
tempPlot <- PlotLMM(m1)
tempPlot


# Generate Q10s for zooplankton in general ----
# Extract slopes and calculate Q10 for overall zooplankton
slopes <- emtrends(m1, ~ rate_name, var = "temp_C")

slopes_Q10 <- as.data.frame(slopes) %>% 
  mutate(Q10 = exp(10* temp_C.trend),
         Q10_lwr = exp(10 * asymp.LCL),
         Q10_upr = exp(10 * asymp.UCL))
slopes_Q10


# Define color palette ----
rate_cols <- c("Clearance"   = "#66c2a5",
               "Ingestion"   = "#fc8d62",
               "Growth"      = "#8da0cb",
               "Respiration" = "#e78ac3")


# Plot allZoop Q10s
allZoopQ10plot <- ggplot() +
  geom_errorbar(data = slopes_Q10, 
                aes(x = rate_name, ymin = Q10_lwr, ymax = Q10_upr), 
                colour = "darkgray",
                width = 0.05, linewidth = 1) +
  geom_point(data = slopes_Q10, 
             aes(x = rate_name, y = Q10, fill = rate_name),
             size = 4, colour = "black", shape = 21) +
  scale_fill_manual(values = rate_cols, guide = "none") +
  labs(x = "Biological rate process",
       y = bquote(bold("Carbon-mass specific Q"[10]))) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  )

allZoopQ10plot


library(patchwork)

tempPlot/allZoopQ10plot +
  plot_layout(guides = "collect")
