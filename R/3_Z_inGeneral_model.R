# Calculating overall rate Q10s
# Josh Hill
# 03/05/2026



  # Here I read in all cleaned data
  # Combine data
  # Fit glmm with random effects to test if the effect of temp on biological rate process
  # Calculate Q10s from the model slopes



# Packages and helpers ----
library(tidyverse)
library(glmmTMB) # for modelling
library(DHARMa) # for diagnostics
library(emmeans) # estimated marginal means
library(MuMIn) # for R2s
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
  filter(rate_name == "IngestionRate") %>% # Filter for ingestion rate
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


# Excretion data
excredat <- readRDS("Data/excrete_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Excretion"))


# Combine them into one...
usedat <- rbind(cleardat, ingdat, grwdat, respdat, excredat)


# Check the temperature range
usedat %>% 
  select(temp_C) %>% 
  summarise( 
    temp_range = paste0(min(temp_C), "-", max(temp_C)))
    # -1.9-31.4 degC


# Quickly view the distribution of raw Cspecific_rates
# I will use these to remove extreme outliers, particularly those that don't make biological sense
usedat %>%
  ggplot() +
  geom_point(aes(x = temp_C, y = Cspecific_rate)) +
  theme_bw() +
  facet_wrap(~ rate_name, scales = "free") + 
  labs(
    x = "Temp C",
    y = "C-specific rate")
  # all seems reasonable and the data looks standard

# Check distribution of data
usedat %>%
  ggplot(aes(x = log(Cspecific_rate))) + # log because we will transform the data for analysis
  geom_histogram(bins = 50, fill = "pink", colour = "grey") +
  theme_bw() +
  facet_wrap(~rate_name, scales = "free") +
  labs(
    x = "Mass-specific rate (Cspecific rates, log scale)",
    y = "Count",
    title = "Distribution of Cspecific rates across all zooplankton")
    # distribution looks pretty normal but ingestion and clearance are slightly left skewed



# Tidy up data and prep data for modelling ----
mdat <- usedat %>% 
  mutate(ln_Cspecific_rate = log(Cspecific_rate)) # log transform mass-specific rate and save as new variable


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

# A complex model with 2 way interactions for temp and rate with random slopes and intercepts
m1 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  diagnose(m1) # I suspect the model does not converge well due to the random effects structure
  summary(m1) # Check random effects
  # I'll drop the random slope for primary reference because the variance is almost 0

  
# A simpler model with only random intercept for primary reference but intercept and slope for taxa
m2 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                (1 | primRef) + (temp_C | taxa),
              data = mdat)
  # Great, model converged...

  # Check diagnostics
  sim <- simulateResiduals(m2)
  plot(sim) # Looks fine, bit of deviation on the tail
  summary(m2)
  # we will swap the random intercepts over and check the model
  
  
# A simpler model with only random intercepts but swapping single intercept structure to taxa...
m3 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                (temp_C | primRef) + (1 | taxa),
              data = mdat)

  # Check diagnostics
  sim <- simulateResiduals(m3)
  plot(sim) 
  summary(m3)

  
# Compare models 2 and 3...m1 did not converge
performance::compare_performance(m2, m3)
# Models are not all mutually nested...we'll just treat AIC/BIC as descriptive...
# Seems like m2 is the best fit so far...we'll refit with ML and check AIC

# Likelihood ratios test of the models
# Refit with ML for valid test on fixed/dispersion and random effect structures
m2_m2 <- update(m2, REML = FALSE)
m3_m3 <- update(m3, REML = FALSE)

anova(m2_m2, m3_m3) # test which random effects structure is a better fit
  # m2 is better

# We will proceed with m2 on the basis of AIC and likelihood ratio tests
summary(m2)
r.squaredGLMM(m2)
# R2m       R2c
# 0.9184957 0.9825001

# Extract slopes and calculate Q10 for all zooplankton ----
slopes <- emmeans::emtrends(m2, ~ rate_name, var = "temp_C")
summary(slopes, infer = TRUE) # test whether each slope is different from zero for each zoopGrp
  # all significantly different to zero
emmeans::contrast(slopes, method = "pairwise", adjust = "mvt") # pairwise test whether slopes differ significantly across rate types
  # clearance - ingest
  # clearance - growth
  # clearance - respiration
# Makes sense... clearance is technically not an energy budget term compared to the others

# Extract intercepts ----
intercepts <- data.frame(emmeans(m2, ~ rate_name, var = "temp_C"))


# Build a parameter dataframe to estimate ratios
params <- data.frame(slopes) %>% 
  rename(slope = temp_C.trend) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL)) %>% 
  left_join(intercepts) %>% 
  rename(intercept = emmean) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL))
  

# Save the parameter data for later, we'll use this to estimate ratios of the terms
# saveRDS(params, file = "Data/modelParameters/allZestimates.rds")


# Get n_obs
n_obs <- mdat %>%
  count(rate_name)

# Get Q10
slopes_Q10 <- as.data.frame(slopes) %>% 
  # Calculate Q10s
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  # Calculate equivalent activation energies (Ea)
  mutate(refT = 15 + 273.15, # reference temp (15degC) in Kelvin
         k = 8.617e-5, # Boltzmann's constant as eV/K-1
         R = 8.314 / 1000, # Ideal gas constant as kJ mol-1
         Ea_eV = k * log(Q10) * (refT * (refT + 10)) / 10, # Ea as eV
         Ea_kJ.mol =  (R) * log(Q10) * (refT * (refT + 10)) / 10 # Ea as kJ/mol-1
         ) %>%
  select(- c(k, R, df, refT, asymp.LCL, asymp.UCL)) %>% 
  left_join(n_obs, by = "rate_name")
slopes_Q10



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
                 "Clearance"   = "bold(Clearance~rate~(ml~mgC^-1~h^-1))",
                 "Ingestion"   = "bold(Ingestion~rate~(mgC~mgC^-1~h^-1))",
                 "Growth"      = "bold(Growth~rate~(mgC~mgC^-1~h^-1))",
                 "Respiration" = "bold(Respiration~rate~(µlO[2]~mgC^-1~h^-1))",
                 "Excretion"   = "bold(Excretion~rate~('mgN-NH'[4]^'+'~mgC^-1~h^-1))"
               ), 
               label_parsed))+ 
    coord_cartesian(xlim = c(-2, 32)) +
    labs(x = "Temp (°C)",
         y = "ln (Carbon-mass specific rate)") +
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


# Plot allZoop Q10s
allZoopQ10plot <- ggplot() +
  geom_errorbar(data = slopes_Q10, 
                aes(x = rate_name, ymin = Q10_lwr, ymax = Q10_upr), 
                colour = "grey",
                width = 0.05, linewidth = 1) +
  geom_point(data = slopes_Q10, 
             aes(x = rate_name, y = Q10),
             size = 4, colour = "midnightblue") +
  geom_text(data = slopes_Q10,
            aes(x = rate_name, y = Q10, label = sprintf("%.2f", Q10)),
            nudge_x = 0.22,
            nudge_y = 0.01) +
  geom_text(data = slopes_Q10,
            aes(x = rate_name, y = Q10_lwr, label = paste0("n = ", n)),
            nudge_y = -0.3,   
            size = 3.5, colour = "grey40") +
  labs(x = "Biological rate process",
       y = bquote(bold("Carbon-mass specific Q"[10]))) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  )

allZoopQ10plot


fig2 <- tempPlot/allZoopQ10plot +
  plot_layout(guides = "collect")
fig2


# Save the plots ----
# ggsave("Output/Figure2/Figure2_tempPlot.pdf", tempPlot, width = 160, height = 150, units = "mm", dpi = 300)
# ggsave("Output/Figure2/Figure2_Q10Plot.pdf", allZoopQ10plot, width = 160, height = 70, units = "mm", dpi = 300)
