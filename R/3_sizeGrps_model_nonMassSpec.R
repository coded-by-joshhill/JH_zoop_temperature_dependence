# Calculating sizeGrps Q10s
# Josh Hill
# 25/02/2026



  # Here I read in all cleaned data
  # Combine data and subset by sizeGrp
  # Fit glmm with random effects to test if the effect of temp on biological rate process
  # Calculate Q10s from the model slopes
  # Save the Q10 and variance into a dataframe
  # Save Arrhenius plot object for later plotting



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
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, rate_value_clean, rate_unit_clean, temp_C, BMC_mg) %>% 
  filter(rate_unit_clean == "ml/ind/hr") %>% 
  drop_na(rate_value_clean) %>% 
  mutate(rate_name = factor("Clearance"))


# Ingestion data
ingdat <- dat %>% 
  filter(rate_name == "IngestionRate") %>% # Filter for clearance rate
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, rate_value_clean, rate_unit_clean, temp_C, BMC_mg) %>% 
  filter(rate_unit_clean == "mgC/ind/hr") %>% 
  drop_na(rate_value_clean) %>% 
  mutate(rate_name = factor("Ingestion"))


# Growth data
grwdat <- readRDS("Data/grwth_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, rate_value_clean, rate_unit_clean, temp_C, BMC_mg) %>% 
  drop_na(rate_value_clean) %>% 
  mutate(rate_name = factor("Growth"))


# Respiration data
respdat <- readRDS("Data/resp_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, rate_value_clean, rate_unit_clean, temp_C, BMC_mg) %>% 
  filter(rate_unit_clean == "ulO2/ind/hr") %>% 
  drop_na(rate_value_clean) %>% 
  mutate(rate_name = factor("Respiration"))


# Custom grouping order # Custom grofactor()uping order 
group_order <- c("Mesoplankton", "Macroplankton") # meso before macro because they're smaller


# Combine them into one dataframe
usedat <- rbind(cleardat, ingdat, grwdat, respdat) %>% 
  filter_out(sizeGrp == "OTHER" | is.na(temp_C)) %>%  # filter out the size group "OTHER" and remove any NAs in temp_C
  mutate(sizeGrp = fct_relevel(sizeGrp, group_order)) # reorder sizeGrp
  

# Check the temperature range
usedat %>% 
  select(temp_C) %>% 
  summarise( 
    temp_range = paste0(min(temp_C), "-", max(temp_C)))
    # -1.8-31 degC


# Quickly view distribution of raw rate_value_cleans
# I will use these to remove extreme outliers, particularly those that don't make biological sense
usedat %>%
  # exclude values that are not biologically reasonable or are extreme outliers
  filter(
    (rate_name == "Clearance" & rate_value_clean < 500000) |
    (rate_name == "Ingestion" & rate_value_clean < 0.1) |
    (rate_name == "Growth" & rate_value_clean < 0.4) |
    (rate_name == "Respiration" & rate_value_clean < 100000)
  ) %>%
  ggplot() +
  geom_point(aes(x = temp_C, y = rate_value_clean)) +
  theme_bw() +
  facet_wrap(~ rate_name, scales = "free") + 
  labs(
    x = "Temp C",
    y = "Clearance rate")


usedat %>%
  # exclude values that are not biologically reasonable or are extreme outliers
  # filter(
  #   (rate_name == "Clearance" & rate_value_clean < 15000) |
  #     (rate_name == "Ingestion" & rate_value_clean < 0.15) |
  #     (rate_name == "Growth" & rate_value_clean < 0.075) |
  #     (rate_name == "Respiration" & rate_value_clean < 60)
  # ) %>%
  ggplot(aes(x = log(rate_value_clean))) + # log because we will transform the data for analysis
  geom_histogram(bins = 50, fill = "pink", colour = "grey") +
  theme_bw() +
  facet_wrap(~rate_name, scales = "free") +
  labs(
    x = "Absolute rates (rate_value_clean, log scale)",
    y = "Count",
    title = "Distribution of raw rates across all zooplankton")
    # distribution looks pretty normal but ingestion and clearance are slightly left skewed



# Tidy up data and prep data for modelling ----
mdat <- usedat %>% 
  # exclude values that are not biologically reasonable or are extreme outliers
  filter(
    (rate_value_clean > 0)) %>%
  mutate(ln_rate_value_clean = log(rate_value_clean)) # log transform absolute rates


# Quick look
mdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_rate_value_clean, , colour = log(BMC_mg), shape = sizeGrp), 
             size = 3) +
  facet_wrap(~ rate_name, scale = "free") +
  ggtitle("Raw data: log(absoluteRates) by temp and sizeGrps, coloured by carbonBodyMass")
  # Looks pretty tidy. Some clear relationships here

mdat %>% 
  ggplot() +
  geom_point(aes(x = log(BMC_mg), y = ln_rate_value_clean, , colour = temp_C, shape = sizeGrp), 
             size = 3) +
  facet_wrap(~ rate_name, scale = "free") +
  ggtitle("Raw data: log(absoluteRates) by carbonBodyMass and sizeGrps, coloured by temp of experiment")

summary(mdat)


# My main question here is...
  # How does temperature dependence vary across zooplankton groups for each rate? AND
  # How does temperature dependence vary across rate processes?

  # So... I will need to model the logAbsoluteRates as a function of 
    # temperature, the rate and the groups...
    # I will need a complex model with interactions to test the difference between the effect of temp on rates and temp on groups.

# My response is log transformed, continuous data and 2 rates are normally distributed and two are slightly left skewed... 
  # This should be fine to Gaussian but I could also check a Gamma with link = log family on the normal absoluteRate data

# I am also mainly interested in the interactions between at least temp:rate and temp:group



# Fit the models ----

# A complex model with 3 way interactions for temp, sizeGrp and rate with random effects
m1 <- glmmTMB(ln_rate_value_clean ~ temp_C * rate_name * sizeGrp + 
                (temp_C | primRef) + (1 | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim) # looks pretty good
  summary(m1)
  r.squaredGLMM(m1)
  # Kind of difficult to interpret I'll fit a simpler model to tease this apart

  
# A simpler model without the 3-way interactions, just 2-way interactions
m2 <- glmmTMB(ln_rate_value_clean ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                temp_C * sizeGrp +  # between temp and different sizeGrps for all rates
                rate_name * sizeGrp + # between rates and sizeGrp 
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m2)
  plot(sim) # Looks fine but seems less scattered compared to m1
  summary(m2)
  r.squaredGLMM(m2)
  # temp:rate - looks like there is significantly different temp dependence among rates for all zoops
  # rate:grp - appears to be no signif differences among rates and sizeGrp but this doesn't include temperature in the interaction...
  # my random effects seem to be soaking up a fairly decent amount of variance, though there is less for the slope compared to intercept


# m1 but without random slopes, just intercepts
m3 <- glmmTMB(ln_rate_value_clean ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                temp_C * sizeGrp +  # between temp and different sizeGrps for all rates
                rate_name * sizeGrp + # between rates and sizeGrp 
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
# likelihood ratios test shows m1 has significantly more explanatory power than both m2 and m3
# AIC is also slightly better despite the BIC being slightly higher
# I will use m1 on the basis of the chisqr test and AIC...

summary(m2)
  


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
    sizeGrp   = unique(mdat$sizeGrp), # levels of sizeGrp
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
    geom_ribbon(data = pop_preds,
                aes(x = temp_C, ymin = conf.low, ymax = conf.high, fill = sizeGrp),
                alpha = 0.25) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate, colour = sizeGrp),
              linewidth = 1) +
    geom_point(data = mdat,
               aes(x = temp_C, y = ln_rate_value_clean, colour = sizeGrp),
               alpha = 0.2) +
    facet_wrap(~rate_name, scales = "free",
               labeller = as_labeller(c(
                 "Clearance"   = "bold(Clearance~rate~(ml~ind^-1~h^-1))",
                 "Ingestion"   = "bold(Ingestion~rate~(mgC~ind^-1~h^-1))",
                 "Growth"      = "bold(Growth~rate~(mgC~ind^-1~h^-1))",
                 "Respiration" = "bold(Respiration~rate~(µlO[2]~ind^-1~h^-1))"
               ), 
               label_parsed))+ 
    coord_cartesian(xlim = c(-2, 32)) +
    labs(x = "Temp (°C)",
         y = "ln(Absolute rate)",
         fill = "Size group",
         colour = "Size group") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "whitesmoke", colour = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold")
    )
}



# Plot it
tempPlot <- PlotLMM(m2)
tempPlot +
  ggtitle("GLMM: log(absoluteRates) ~ temp * sizeGrp + temp * rate + rate * sizeGrp")

  
  
# Extract slopes using and calculate Q10 for each sizeGrp
sizeGrp_slopes <- emtrends(m2, ~ rate_name * sizeGrp, var = "temp_C")

sizeGrp_slopes_Q10 <- as.data.frame(sizeGrp_slopes) %>% 
  mutate(Q10 = exp(10 * temp_C.trend),
         Q10_lwr = exp(10 * asymp.LCL),
         Q10_upr = exp(10 * asymp.UCL))
sizeGrp_slopes_Q10



# Plot sizeGrp Q10s
sizeGQ10plot <- ggplot() +
  geom_errorbar(data = sizeGrp_slopes_Q10, 
                aes(x = sizeGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = sizeGrp),
                width = .05,
                linewidth = 1) +
  geom_point(data = sizeGrp_slopes_Q10, aes(x = sizeGrp, y = Q10),
             size = 3,
             colour = "black") +
  facet_wrap(~rate_name, scales = "free",
             labeller = as_labeller(c(
               "Clearance"   = "bold(Clearance~rate~(ml~ind^-1~h^-1))",
               "Ingestion"   = "bold(Ingestion~rate~(mgC~ind^-1~h^-1))",
               "Growth"      = "bold(Growth~rate~(mgC~ind^-1~h^-1))",
               "Respiration" = "bold(Respiration~rate~(µlO[2]~ind^-1~h^-1))"
             ), 
             label_parsed))+  
  labs(x = "Size group",
       y = bquote(bold("Absolute Q"[10])),
       colour = "Size group") +
  theme(
    strip.background = element_rect(fill = "whitesmoke", colour = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )
sizeGQ10plot


library(patchwork)

tempPlot/sizeGQ10plot +
  plot_layout(guides = "collect") 



# Lets generate Q10s for zooplankton in general ----
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
       y = bquote(bold("Absolute specific Q"[10]))) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  )

allZoopQ10plot

