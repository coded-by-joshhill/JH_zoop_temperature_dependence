# Calculating zoopGrps Q10s
# Josh Hill
# 25/02/2026

# As it stands this analyses isn't going in the manuscript because there is not enough data across all zoopGrps for each rate...so the model has to drop several terms



  # Here I read in all cleaned data
  # Combine data and subset by zoopGrp
  # Fit glmm with random effects to test if the effect of temp on biological rate process
  # Calculate Q10s from the model slopes



# Packages and helpers ----
library(tidyverse)
library(glmmTMB) # for modelling
library(DHARMa) # for diagnostics
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


# Custom grouping order
group_order <- c("Ctenophores",
                 "Cnidarians",
                 "Chaetognaths",
                 "Amphipods",
                 "Copepods",
                 "Decapods",
                 "Euphausiids",
                 "Mysids",
                 "Appendicularians",
                 "Thaliaceans")


# Combine them into one dataframe
usedat <- rbind(cleardat, ingdat, grwdat, respdat) %>% 
  filter_out(zoopGrp == "OTHER" | is.na(temp_C)) %>%  # filter out the size group "OTHER" and remove any NAs in temp_C
  mutate(zoopGrp = fct_relevel(zoopGrp, group_order)) # reorder zoopGrp
  

# Check the temperature range
usedat %>% 
  group_by(zoopGrp, rate_name) %>% 
  select(temp_C) %>% 
  summarise( 
    temp_range = paste0(min(temp_C), "-", max(temp_C))) %>% 
  arrange(rate_name, zoopGrp) %>% print(n = "Inf")
  # Clearance:
    # won't be able to get estimate for chaetognaths
  # Ingestion:
    # Won't be able to get estimate for chaetognaths and appendicularians
  # Growth:
    # won't be able to get estimates for mysids
  # Respiration:
    # wont be able to get estimates for amphipods, decapods


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
  # exclude groups lacking data
  filter_out(
    zoopGrp == "Chaetognaths" |
      zoopGrp == "Annelids" |
      zoopGrp == "Mysids" |
      zoopGrp == "Decapods" |
      zoopGrp == "Amphipods" |
      zoopGrp == "Appendicularians" |
      zoopGrp == "Thaliaceans"
  ) %>%
  mutate(ln_Cspecific_rate = log(Cspecific_rate)) # log transform mass-specific rate


# Quick look
mdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = zoopGrp)) +
  facet_wrap(~ rate_name, scale = "free") +
  theme_bw()
  # Looks pretty tidy. Some clear relationships here

summary(mdat)


# My main question here is...
  # How does temperature dependence vary across zooplankton groups for each rate? AND
  # How does temperature dependence vary across rate processes?

  # So... I will need to model the logMassSpecificRates as a function of 
    # temperature, the rate and the groups...
    # I will need a complex model with interactions to test the difference between the effect of temp on rates and temp on groups.

# My response is log transformed, continuous data and 2 rates are normally distributed and two are slightly left skewed... 
  # This should be fine to Gaussian but I could also check a Gamma with link = log family on the normal massSpecRate data

# I am also mainly interested in the interactions between at least temp:rate and temp:group



# Fit the models ----

# A complex model with 3 way interactions for temp, zoopGrp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ temp_C * rate_name * zoopGrp + 
                (temp_C | primRef) + (1 | taxa), # with primRef and taxa as random intercepts and slopes
              dispformula = ~rate_name,
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim) # looks pretty good
  summary(m1)
  # Kind of difficult to interpret I'll fit a simpler model to tease this apart

  
# A simpler model without the 3-way interactions, just 2-way interactions
m2 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                temp_C * zoopGrp +  # between temp and different zoopGrps for all rates
                rate_name * zoopGrp + # between rates and zoopGrp 
                (temp_C | primRef) + (1 | taxa), # with primRef and taxa as random intercepts and slopes
              dispformula = ~rate_name,
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m2)
  plot(sim) # Looks fine but seems less scattered compared to m1
  summary(m2)
  # temp:rate - looks like there is significantly different temp dependence among rates for all zoops
  # rate:grp - appears to be no signif differences among rates and zoopGrp but this doesn't include temperature in the interaction...
  # my random effects seem to be soaking up a fairly decent amount of variance, though there is less for the slope compared to intercept


# m1 but without random slopes, just intercepts
m3 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                temp_C * zoopGrp +  # between temp and different zoopGrps for all rates
                rate_name * zoopGrp + # between rates and zoopGrp 
                (1 | primRef) + (1 | taxa), # with primRef and taxa as random intercepts
              dispformula = ~rate_name,
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m3)
  plot(sim) # Looks fine but seems less scattered compared to m1
  summary(m3)
  r.squaredGLMM(m3)


performance::compare_performance(m1, m2, m3)


# Likelihood ratios test of the models
anova(m1, m2, m3)
# likelihood ratios test shows m2 has significantly more explanatory power than both m1 and m3
# BIC is also slightly better
# I will use m2 on the basis of the chisqr test and AIC...

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
    zoopGrp   = unique(mdat$zoopGrp), # levels of zoopGrp
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
                aes(x = temp_C, ymin = conf.low, ymax = conf.high, fill = zoopGrp),
                alpha = 0.25) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate, colour = zoopGrp),
              linewidth = 1) +
    geom_point(data = mdat,
               aes(x = temp_C, y = ln_Cspecific_rate, colour = zoopGrp),
               alpha = 0.2) +
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
         y = "ln(Carbon-mass specific rate)",
         fill = "Taxonomic group",
         colour = "Taxonomic group") +
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
tempPlot

  
# Extract slopes using and calculate Q10 for each zoopGrp
zoopGrp_slopes <- emtrends(m1, ~ rate_name * zoopGrp, var = "temp_C")

zoopGrp_slopes_Q10 <- as.data.frame(zoopGrp_slopes) %>% 
  mutate(Q10 = exp(10 * temp_C.trend),
         Q10_lwr = exp(10 * asymp.LCL),
         Q10_upr = exp(10 * asymp.UCL)) %>% 
  arrange(rate_name)
zoopGrp_slopes_Q10



# Plot zoopGrp Q10s
zoopGQ10plot <- ggplot() +
  geom_errorbar(data = zoopGrp_slopes_Q10, 
                aes(x = zoopGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = zoopGrp),
                width = .05,
                linewidth = 1) +
  geom_point(data = zoopGrp_slopes_Q10, aes(x = zoopGrp, y = Q10),
             size = 3,
             colour = "black") +
  facet_wrap(~rate_name, scales = "free",
             labeller = as_labeller(c(
               "Clearance"   = "bold(Maximum~clearance~rate~(ml~mgC^-1~h^-1))",
               "Ingestion"   = "bold(Maximum~ingestion~rate~(mgC~mgC^-1~h^-1))",
               "Growth"      = "bold(Growth~rate~(mgC~mgC^-1~h^-1))",
               "Respiration" = "bold(Respiration~rate~(µlO[2]~mgC^-1~h^-1))"
             ), 
             label_parsed))+  
  labs(x = "Size group",
       y = bquote(bold("Carbon-mass specific Q"[10])),
       colour = "Size group") +
  theme(
    strip.background = element_rect(fill = "whitesmoke", colour = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )
zoopGQ10plot


library(patchwork)

tempPlot/funcGQ10plot +
  plot_layout(guides = "collect") 


