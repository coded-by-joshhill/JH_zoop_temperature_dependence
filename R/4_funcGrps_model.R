# Calculating funcGrps Q10s
# Josh Hill
# 26/02/2026



  # Here I read in all cleaned data
  # Combine data and subset by funcGrps
  # Fit glmm with random effects to test if the effect of temp on biological rate process
  # Calculate Q10s from the model slopes
  # Save the Q10 and variance into a dataframe
  # Save Arrhenius plot object for later plotting



# Packages and helpers ----
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(marginaleffects)
library(emmeans) # Estimated marginal means
theme_set(new = theme_bw())



# Read in the data and prep for analysis ----
# Feeding data
dat <- readRDS("Data/clear_ingest_data.rds")

cleardat <- dat %>% 
  filter(rate_name == "ClearanceRate") %>% # Filter for clearance rate
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = "Clearance")


ingdat <- dat %>% 
  filter(rate_name == "IngestionRate") %>% # Filter for clearance rate
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = "Ingestion")


# Growth data
grwdat <- readRDS("Data/grwth_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = "Growth")


# Respiration data
respdat <- readRDS("Data/resp_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = "Respiration")


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
  filter_out(funcGrp == "OTHER" | is.na(temp_C)) %>% 
  mutate(ln_Cspecific_rate = log(Cspecific_rate)) # log transform mass-specific rate


# Quick look
mdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = funcGrp)) +
  facet_wrap(~ rate_name) +
  theme_bw()



# Fit the models ----


# A model with 3 way interactions for temp, funcGrp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ temp_C * rate_name * funcGrp + 
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim) # looks pretty good, there is a slight deviation in the observed residuals at the lower end of the QQplot
  summary(m1)
  r.squaredGLMM(m1)
  
  
# A simpler model without the 3-way interaction
m2 <- glmmTMB(ln_Cspecific_rate ~ temp_C * rate_name + temp_C * funcGrp + rate_name * funcGrp +
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m2)
  plot(sim) # looks pretty schmick
  summary(m2)
  r.squaredGLMM(m2)

# Likelihood ratios test of the two models
anova(m1, m2) # anova of complex vs simpler models
# likelihood ratios test shows m1 has significantly more explanatory power
# AIC is also slightly better despite the BIC being slightly higher
# I will use m1 on the basis of the chisqr test and AIC...



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
    rate_name = unique(mdat$rate_name) # levels of rates
  )
  
  # Population-level predictions (fixed effects only)
  pop_preds <- newdat
  pred <- predict(model, newdata = newdat, se.fit = TRUE,
                  re.form = NA, allow.new.levels = TRUE)
  pop_preds$estimate <- pred$fit
  pop_preds$conf.low  <- pred$fit - 1.96 * pred$se.fit
  pop_preds$conf.high <- pred$fit + 1.96 * pred$se.fit
  
  # Plot it up...
  ggplot() +
    geom_ribbon(data = pop_preds,
                aes(x = temp_C, ymin = conf.low, ymax = conf.high, fill = funcGrp),
                alpha = 0.3) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate, colour = funcGrp),
              linewidth = 1.5) +
    geom_point(data = mdat,
               aes(x = temp_C, y = ln_Cspecific_rate, colour = funcGrp),
               alpha = 0.4) +
    facet_wrap(~rate_name, scales = "free")
}

A <- PlotLMM(m1)
A



# Extract slopes and calculate Q10 for each funcGrp
funcGrp_slopes <- emtrends(m1, ~ rate_name * funcGrp, var = "temp_C")

funcGrp_slopes_Q10 <- as.data.frame(funcGrp_slopes) |>
  mutate(Q10 = exp(10 * temp_C.trend),
         Q10_lwr = exp(10 * asymp.LCL),
         Q10_upr = exp(10 * asymp.UCL)) %>% 
  arrange(rate_name, funcGrp)
funcGrp_slopes_Q10



# Plot funcGrp Q10s
sizeGQ <- ggplot() +
  geom_errorbar(data = funcGrp_slopes_Q10, 
                aes(x = funcGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = funcGrp),
                width = .05,
                linewidth = 1) +
  geom_point(data = funcGrp_slopes_Q10, aes(x = funcGrp, y = Q10),
             size = 3,
             colour = "black") +
  facet_wrap(~rate_name, scale = "free")
sizeGQ

library(patchwork)

A/sizeGQ
