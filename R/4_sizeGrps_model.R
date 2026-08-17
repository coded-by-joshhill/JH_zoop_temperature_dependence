# Calculating sizeGrps Q10s
# Josh Hill
# 03/05/2026



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

# Excretion data
excredat <- readRDS("Data/excrete_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Excretion"))

# Custom grouping order
group_order <- c("Mesoplankton", "Macroplankton") # meso before macro because they're smaller


# Combine them into one dataframe
usedat <- rbind(cleardat, ingdat, grwdat, respdat, excredat) %>% 
  filter_out(sizeGrp == "OTHER" | is.na(temp_C)) %>%  # filter out the size group "OTHER" and remove any NAs in temp_C
  mutate(sizeGrp = fct_relevel(sizeGrp, group_order)) # reorder sizeGrp
  

# Check the temperature range
usedat %>% 
  group_by(sizeGrp, rate_name) %>% 
  select(temp_C) %>% 
  summarise( 
    temp_range = paste0(min(temp_C), "-", max(temp_C))) %>% 
  arrange(rate_name, sizeGrp)
  # looks good for all groups



# Tidy up data and prep data for modelling ----
mdat <- usedat %>% 
  group_by(sizeGrp) %>% 
  filter(n() >= 15, # Exclude sizeGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  mutate(ln_Cspecific_rate = log(Cspecific_rate)) # log transform mass-specific rate

# Save mdat so we can estimate number of taxa for ratios later...
saveRDS(mdat, file = "Data/modelParameters/sizeGrp_mdat.rds")

# Quick look
mdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = sizeGrp)) +
  facet_wrap(~ rate_name, scale = "free")
  # Looks pretty tidy. Some clear relationships here
  # Note: 
  # clearance is ml/mgC/hr
  # ingestion is mgC/mgC/hr 
  # growth is mgC/mgC/hr
  # respiration is uLO2/mgC/hr
  # excretion is mgN-NH4+/mgC/hr
summary(mdat)


# My main questions here are...
  # How does temperature dependence vary across zooplankton groups for each rate? AND
  # How does temperature dependence vary across rate processes among groups?

# I am also mainly interested in the interactions between at least temp:rate and temp:group



# Fit the models ----

# A complex model with 3 way interactions for temp, sizeGrp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * sizeGrp + # three-way interaction
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 
  # model doesn't converge...will check the random effects variance

  # Check diagnostics
  summary(m1)
  # primRef random slope is close to zero so I will drop that and re-fit
  
  
# Refit m1 as m2 without primRef as random slope
m2 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * sizeGrp + # three-way interaction
                (1 | primRef) + (temp_C | taxa), # with primRef as random intercept and taxa as random intercept and slope
              data = mdat) 
  # Great, model converged

  # Check diagnostics
  sim <- simulateResiduals(m2)
  plot(sim) # Looks fine, could be a couple of outliers based on the residuals plot...it is probably the appendicularians in the excretion dataset
  # Plot without tests
  par(mfrow = c(1, 2))
  plotQQunif(sim, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
  plotResiduals(sim)
  par(mfrow = c(1, 1)) # Reset back to normal
  
  summary(m2)
  
  
# Fit a simpler model without the 3-way interactions, just 2-way interactions, same random effect structure
m3 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                temp_C * sizeGrp +  # between temp and different sizeGrps for all rates
                rate_name * sizeGrp + # between rates and sizeGrp 
                (1 | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m3)
  plot(sim) # Looks fine but seems less scattered compared to m1, QQ plot looks slightly better though
  summary(m3)
  

# Compare models with different fixed effects
performance::compare_performance(m2, m3)
# Models are not all mutually nested...we'll just treat AIC as descriptive...
  # 3-way interaction appears better so far (m2)


# Likelihood ratios test of the models
# Refit with ML to test the different fixed effect structures
m2_m2 <- update(m2, REML = FALSE)
m3_m3 <- update(m3, REML = FALSE)


anova(m2_m2, m3_m3) # test if the three-way interaction is better than the two-way interaction
  # m2 with 3-way interaction is better

# I will use m2 on the basis of the chisqr test and AIC...
summary(m2)



# Extract slopes using and calculate Q10 for each sizeGrp ----
sizeGrp_slopes <- emtrends(m2, ~ rate_name * sizeGrp, var = "temp_C")
summary(sizeGrp_slopes, infer = TRUE) %>% arrange(sizeGrp, rate_name) # test whether each slope is different from zero for each sizeGrp
# yes, all different to zero
emmeans::contrast(sizeGrp_slopes, by = "rate_name", method = "pairwise", adjust = "none") # pairwise test whether slopes differ significantly across rate types and grps
  # yes, clearance: meso-macro



# Extract intercepts ----
sizeGrp_intercepts <- data.frame(emmeans(m2, ~ rate_name * sizeGrp, var = "temp_C"))


# Build a parameter dataframe and to estimate ratios
params <- data.frame(sizeGrp_slopes) %>% 
  rename(slope = temp_C.trend) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL)) %>% 
  left_join(sizeGrp_intercepts) %>% 
  rename(intercept = emmean) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL))


# Save the parameter data for later, we'll use this to estimate ratios of the terms
saveRDS(params, file = "Data/modelParameters/sizeGrpestimates.rds")


# Get n
n_obs <- mdat %>%
  count(rate_name, sizeGrp)

n_obs %>% arrange(sizeGrp, rate_name)

# Get Q10
sizeGrp_slopes_Q10 <- as.data.frame(sizeGrp_slopes) %>% 
  # Calculate Q10s
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  # Calculate equivalent activation energies (Ea)
  mutate(refT = 15 + 273.15, # reference temp in Kelvin
         k = 8.617e-5, # Boltzmann's constant as eV/K-1
         R = 8.314 / 1000, # Ideal gas constant as kJ mol-1
         Ea_eV = round(k * log(Q10) * (refT * (refT + 10)) / 10, digits = 2), # Ea as eV
         Ea_kJ.mol = round((R) * log(Q10) * (refT * (refT + 10)) / 10) # Ea as kJ/mol-1
         ) %>%
  select(- c(k, R, df, refT, asymp.LCL, asymp.UCL)) %>%
  left_join(n_obs, by = c("rate_name", "sizeGrp")) %>% 
  arrange(rate_name, sizeGrp)
sizeGrp_slopes_Q10



# Prep for plotting ----

# Define group colours
grp_cols <- c("Mesoplankton" = "#0077BB",
              "Macroplankton" = "#EE7733")

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
    rate_name = unique(mdat$rate_name)) # levels of rates
  

  # Zooplankton-level (population) predictions
  pop_preds <- newdat
  pred <- predict(model, newdata = newdat, se.fit = TRUE,
                  type = "response",
                  re.form = NA)
  pop_preds$estimate <- pred$fit
  pop_preds$conf.low <- pred$fit - 1.96 * pred$se.fit
  pop_preds$conf.high <- pred$fit + 1.96 * pred$se.fit
  

    # Plot it up...
  ggplot() +
    geom_point(data = mdat,
               aes(x = temp_C, y = ln_Cspecific_rate, colour = sizeGrp),
               alpha = 0.2) +
    geom_ribbon(data = pop_preds,
                aes(x = temp_C, ymin = conf.low, ymax = conf.high, fill = sizeGrp),
                alpha = 0.25) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate, colour = sizeGrp),
              linewidth = 1) +

    facet_wrap(~rate_name, scales = "free",
               labeller = as_labeller(c(
                 "Clearance"   = "bold(Clearance~(ml~mgC^-1~h^-1))",
                 "Ingestion"   = "bold(Ingestion~(mgC~mgC^-1~h^-1))",
                 "Growth"      = "bold(Growth~(mgC~mgC^-1~h^-1))",
                 "Respiration" = "bold(Respiration~(µlO[2]~mgC^-1~h^-1))",
                 "Excretion"   = "bold(Excretion~('mgN-NH'[4]^'+'~mgC^-1~h^-1))"
               ), 
               label_parsed),
               ncol = 1)+ 
    scale_fill_manual(values = grp_cols, labels = c("Mesozooplankton", "Macrozooplankton")) +
    scale_colour_manual(values = grp_cols, labels = c("Mesozooplankton", "Macrozooplankton")) +
    coord_cartesian(xlim = c(-2, 32)) +
    labs(x = "Temp (°C)",
         y = "ln(Carbon-mass specific rate)",
         fill = "Size group",
         colour = "Size group") +
    theme(
      strip.background = element_rect(fill = "whitesmoke", colour = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.position = "top",
    )
}


# Plot it
tempPlot <- PlotLMM(m2)
tempPlot


# Plot sizeGrp Q10s
sizeGQ10plot <- ggplot() +
  geom_errorbar(data = sizeGrp_slopes_Q10, 
                aes(x = sizeGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = sizeGrp),
                width = .1,
                linewidth = 1) +
  geom_point(data = sizeGrp_slopes_Q10, aes(x = sizeGrp, y = Q10),
             size = 3,
             colour = "black") +
  geom_text(data = sizeGrp_slopes_Q10,
            aes(x = sizeGrp, y = Q10, label = sprintf("%.2f", Q10)),
            nudge_x = 0.21,
            nudge_y = 0.0) + 
  geom_text(data = sizeGrp_slopes_Q10,
            aes(x = sizeGrp, y = Q10_lwr, label = paste0("n = ", n)),
            nudge_y = -0.3,   
            size = 3.5, colour = "grey40") +
  facet_wrap(~rate_name, scales = "free",
             labeller = as_labeller(c(
               "Clearance"   = "bold(Clearance~(ml~mgC^-1~h^-1))",
               "Ingestion"   = "bold(Ingestion~(mgC~mgC^-1~h^-1))",
               "Growth"      = "bold(Growth~(mgC~mgC^-1~h^-1))",
               "Respiration" = "bold(Respiration~(µlO[2]~mgC^-1~h^-1))",
               "Excretion"   = "bold(Excretion~('mgN-NH'[4]^'+'~mgC^-1~h^-1))"
             ), 
             label_parsed), 
             ncol = 1) + 
  scale_colour_manual(values = grp_cols, labels = c("Mesozooplankton", "Macrozooplankton")) +
  labs(x = "Size group",
       y = bquote(bold("Carbon-mass specific Q"[10])),
       colour = "Size group") +
  scale_x_discrete(labels = c("Mesozooplankton", "Macrozooplankton")) +
  coord_cartesian(ylim = c(-1, 8)) +
  theme(
    strip.background = element_rect(fill = "whitesmoke", colour = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "top")
sizeGQ10plot

tempPlot+sizeGQ10plot

# Save it
ggsave("Output/Figure3_raw/Figure3_tempPlot.pdf", tempPlot, width = 75, height = 240, units = "mm", dpi = 300)
ggsave("Output/Figure3_raw/Figure3_Q10Plot.pdf", sizeGQ10plot, width = 75, height = 240, units = "mm", dpi = 300)
