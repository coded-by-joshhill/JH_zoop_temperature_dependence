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
theme_set(new = theme_bw())



# Read in the data and prep for analysis ----

# Feeding data
dat <- readRDS("Data/clear_ingest_data.rds")


# Clearance data
cleardat <- dat %>% 
  filter(rate_name == "ClearanceRate") %>% # Filter for clearance rate
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Clearance"))


# Ingestion data
ingdat <- dat %>% 
  filter(rate_name == "IngestionRate") %>% # Filter for clearance rate
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Ingestion"))


# Growth data
grwdat <- readRDS("Data/grwth_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Growth"))


# Respiration data
respdat <- readRDS("Data/resp_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, final_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
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
m1 <- glmmTMB(ln_Cspecific_rate ~ temp_C * rate_name * sizeGrp + 
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim) # looks pretty good
  summary(m1)
  r.squaredGLMM(m1)
  # Kind of difficult to interpret I'll fit a simpler model to tease this apart

  
# A simpler model without the 3-way interactions, just 2-way interactions
m2 <- glmmTMB(ln_Cspecific_rate ~ 
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
m3 <- glmmTMB(ln_Cspecific_rate ~ 
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



# Likelihood ratios test of the models
anova(m1, m2, m3)
# likelihood ratios test shows m1 has significantly more explanatory power than both m2 and m3
# AIC is also slightly better despite the BIC being slightly higher
# I will use m1 on the basis of the chisqr test and AIC...

summary(m1)
  # Family: gaussian  ( identity )
  # Formula:          ln_Cspecific_rate ~ temp_C * rate_name * sizeGrp + (temp_C |      primRef) + (temp_C | taxa)
  # Data: mdat
  # 
  # AIC       BIC    logLik -2*log(L)  df.resid 
  # 7881.0    8020.3   -3917.5    7835.0      3131 
  # 
  # Random effects:
  #   
  # Conditional model:
  #   Groups   Name        Variance  Std.Dev. Corr  
  #   primRef  (Intercept) 1.2713938 1.12756        
  #   temp_C      0.0017736 0.04211  -0.47 
  #   taxa     (Intercept) 0.5433531 0.73712        
  #   temp_C      0.0005054 0.02248  -0.86 
  #   Residual             0.5699156 0.75493        
  #  Number of obs: 3154, groups:  primRef, 151; taxa, 220
  # 
  # Dispersion estimate for gaussian family (sigma^2): 0.57 
  # 
  # Conditional model:
  #   Estimate Std. Error z value Pr(>|z|)    
  #   (Intercept)                                       3.905063   0.459163   8.505  < 2e-16 ***
  #   temp_C                                            0.107480   0.026223   4.099 4.16e-05 ***
  #   rate_nameIngestion                               -9.558626   0.482632 -19.805  < 2e-16 ***
  #   rate_nameGrowth                                  -9.535298   0.583282 -16.348  < 2e-16 ***
  #   rate_nameRespiration                             -3.358547   0.526216  -6.382 1.74e-10 ***
  #   sizeGrpMacroplankton                             -1.280160   0.550483  -2.326  0.02004 *  
  #   temp_C:rate_nameIngestion                         0.009251   0.027928   0.331  0.74046    
  #   temp_C:rate_nameGrowth                           -0.139961   0.033256  -4.209 2.57e-05 ***
  #   temp_C:rate_nameRespiration                      -0.006323   0.030730  -0.206  0.83696    
  #   temp_C:sizeGrpMacroplankton                       0.065112   0.031244   2.084  0.03717 *  
  #   rate_nameIngestion:sizeGrpMacroplankton           2.260457   0.542997   4.163 3.14e-05 ***
  #   rate_nameGrowth:sizeGrpMacroplankton              0.627717   0.692364   0.907  0.36460    
  #   rate_nameRespiration:sizeGrpMacroplankton         1.287128   0.564210   2.281  0.02253 *  
  #   temp_C:rate_nameIngestion:sizeGrpMacroplankton   -0.130810   0.031807  -4.113 3.91e-05 ***
  #   temp_C:rate_nameGrowth:sizeGrpMacroplankton      -0.038390   0.038863  -0.988  0.32324    
  #   temp_C:rate_nameRespiration:sizeGrpMacroplankton -0.088665   0.032963  -2.690  0.00715 ** 
  #   ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

  # 
  # temp:rate:sizegrp - there is significantly different temp dependence among all rates and grps for: ingestion and respiration
    # growth is not significantly different between meso and macroplankton but growth is significantly different to other rates...



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
  
  
# Extract slopes using and calculate Q10 for each sizeGrp
sizeGrp_slopes <- emtrends(m1, ~ rate_name * sizeGrp, var = "temp_C")

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
               "Clearance"   = "bold(Clearance~rate~(ml~mgC^-1~h^-1))",
               "Ingestion"   = "bold(Ingestion~rate~(mgC~mgC^-1~h^-1))",
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
       y = bquote(bold("Carbon-mass specific Q"[10]))) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  )

allZoopQ10plot

