# Calculating sizeGrps Q10s
# Josh Hill
# 11/03/2026



  # Here I read in all cleaned data
  # Combine data and subset by sizeGrp
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
group_order <- c("Mesoplankton", "Macroplankton") # meso before macro because they're smaller


# Combine them into one dataframe
usedat <- rbind(cleardat, ingdat, grwdat, respdat) %>% 
  filter_out(sizeGrp == "OTHER" | is.na(temp_C)) %>%  # filter out the size group "OTHER" and remove any NAs in temp_C
  mutate(sizeGrp = fct_relevel(sizeGrp, group_order)) # reorder sizeGrp
  

# Check the temperature range
usedat %>% 
  group_by(sizeGrp, rate_name) %>% 
  select(temp_C) %>% 
  summarise( 
    temp_range = paste0(min(temp_C), "-", max(temp_C))) %>% 
  arrange(rate_name, sizeGrp)
  # looks good


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
  # exclude values that are not biologically reasonable or are extreme outliers
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
  group_by(sizeGrp) %>% 
  filter(n() >= 15, # Exclude sizeGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  mutate(ln_Cspecific_rate = log(Cspecific_rate)) # log transform mass-specific rate


# Quick look
mdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = sizeGrp)) +
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

# A complex model with 3 way interactions for temp, sizeGrp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * sizeGrp + # three-way interaction
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              dispformula = ~rate_name,
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim) # looks pretty good
  summary(m1)
  # Kind of difficult to interpret I'll use emtrends to extract rate-specific slopes later
  
  
# A simpler model without the 3-way interactions, just 2-way interactions, same random effect structure
m2 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                temp_C * sizeGrp +  # between temp and different sizeGrps for all rates
                rate_name * sizeGrp + # between rates and sizeGrp 
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              dispformula = ~rate_name,
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m2)
  plot(sim) # Looks fine but seems less scattered compared to m1
  summary(m2)
  

# m1 but without random slopes, just intercepts
m3 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * sizeGrp + # three-way interaction
                (1 | primRef) + (1 | taxa), # with primRef and taxa as random intercepts
              dispformula = ~rate_name,
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m3)
  plot(sim) # Looks fine but seems less scattered compared to m1
  summary(m3)
  

# m1 but without dispersion model, including
m4 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * sizeGrp + # three-way interaction
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              dispformula = ~1,
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m4)
  plot(sim) # Looks fine but seems less scattered compared to m1
  summary(m4)
  

# Compare models
performance::compare_performance(m1, m2, m3, m4)
# Models are not all mutually nested...we'll just treat AIC/BIC as descriptive...


# Likelihood ratios test of the models
# Refit with REML for valid test on fixed/dispersion and random effect structures
m1_m1 <- update(m1, REML = FALSE)
m2_m2 <- update(m2, REML = FALSE)
m3_m3 <- update(m3, REML = FALSE)
m4_m4 <- update(m4, REML = FALSE)

anova(m1_m1, m2_m2) # test if the three-way interaction is better than the two-way interaction
  # m1 with 3-way interaction is better

anova(m1_m1, m4_m4) # test if dispersion submodel is justified
  # yes, m1 with dispersion structure is better

anova(m1_m1, m3_m3) # test if 3 way interaction with simpler RE structure is better
  # m1 RE structure is better

anova(m1_m1, m2_m2, m3_m3,m4_m4) # although models are not mutally nested, lets look at all models for descriptive purposes
# likelihood ratios test across the board shows m1 has significantly more explanatory power than other models
# AIC and BIC is also slightly better
# I will use m1 on the basis of the chisqr test and AIC...
summary(m1)



# Extract slopes using and calculate Q10 for each sizeGrp ----
sizeGrp_slopes <- emtrends(m1, ~ rate_name * sizeGrp, var = "temp_C")
summary(sizeGrp_slopes, infer = TRUE) 
test(sizeGrp_slopes) # test whether each slope is different from zero
pairs(sizeGrp_slopes, by = "rate_name") # pairwise test whether slopes differ significantly across rate types and grps
  # yes, some slightly significant differences

# Get Q10
sizeGrp_slopes_Q10 <- as.data.frame(sizeGrp_slopes) %>% 
  mutate(Q10 = exp(10 * temp_C.trend),
         Q10_lwr = exp(10 * asymp.LCL),
         Q10_upr = exp(10 * asymp.UCL)) %>% 
  arrange(rate_name, sizeGrp)
sizeGrp_slopes_Q10




# Prep for plotting ----

# Define group colours
grp_cols <- c("Mesoplankton" = "#fc8d59",
              "Macroplankton" = "#91bfdb")

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
               aes(x = temp_C, y = ln_Cspecific_rate, colour = sizeGrp),
               alpha = 0.2) +
    facet_wrap(~rate_name, scales = "free",
               labeller = as_labeller(c(
                 "Clearance"   = "bold(Clearance~rate~(ml~mgC^-1~h^-1))",
                 "Ingestion"   = "bold(Ingestion~rate~(mgC~mgC^-1~h^-1))",
                 "Growth"      = "bold(Growth~rate~(mgC~mgC^-1~h^-1))",
                 "Respiration" = "bold(Respiration~rate~(µlO[2]~mgC^-1~h^-1))"
               ), 
               label_parsed))+ 
    scale_fill_manual(values = grp_cols, labels = c("Mesoplankton", "Macroplankton")) +
    scale_colour_manual(values = grp_cols, labels = c("Mesoplankton", "Macroplankton")) +
    coord_cartesian(xlim = c(-2, 32)) +
    labs(x = "Temp (°C)",
         y = "ln(Carbon-mass specific rate)",
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
tempPlot <- PlotLMM(m1)
tempPlot



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
               "Clearance"   = "bold(Clearance~rate~(ml~mgC^-1~h^-1))",
               "Ingestion"   = "bold(Ingestion~rate~(mgC~mgC^-1~h^-1))",
               "Growth"      = "bold(Growth~rate~(mgC~mgC^-1~h^-1))",
               "Respiration" = "bold(Respiration~rate~(µlO[2]~mgC^-1~h^-1))"
             ), 
             label_parsed)) + 
  geom_text(data = sizeGrp_slopes_Q10,
            aes(x = sizeGrp, y = Q10, label = sprintf("%.2f", Q10)),
            nudge_x = 0.21,
            nudge_y = 0.0) +
  scale_colour_manual(values = grp_cols, labels = c("Mesoplankton", "Macroplankton")) +
  labs(x = "Size group",
       y = bquote(bold("Carbon-mass specific Q"[10])),
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

# Save it
ggsave("Output/Figure3/Figure3_tempPlot.pdf", tempPlot, width = 160, height = 150, units = "mm", dpi = 300)
ggsave("Output/Figure3/Figure3_Q10Plot.pdf", sizeGQ10plot, width = 160, height = 70, units = "mm", dpi = 300)

