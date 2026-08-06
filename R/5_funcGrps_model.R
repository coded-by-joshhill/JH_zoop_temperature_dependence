# Calculating funcGrps Q10s
# Josh Hill
# 05/05/2026



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


# Excretion data
excredat <- readRDS("Data/excrete_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Excretion"))


# Custom grouping order
group_order <- c("Crustaceans", 
                 "GelPreds", 
                 "GelFilter")


# Combine them into one dataframe
usedat <- rbind(cleardat, ingdat, grwdat, respdat, excredat) %>% 
  filter_out(funcGrp == "OTHER" | is.na(temp_C)) %>%  # filter out the functional group "OTHER" and remove any NAs in temp_C
  mutate(funcGrp = fct_relevel(funcGrp, group_order)) # reorder funcGrp
  

# Check the temperature range
usedat %>% 
  group_by(funcGrp, rate_name) %>% 
  select(temp_C) %>% 
  summarise( 
    temp_range = paste0(min(temp_C), "-", max(temp_C))) %>% 
  arrange(rate_name, funcGrp)
# GelFilter may have issues due to narrow temp range for clearance, ingestion, growth...
# GelPreds may have issues with ingestion



# Tidy up data and prep data for modelling ----
mdat <- usedat %>% 
  group_by(funcGrp) %>% 
  filter(n() >= 15, # Exclude funcGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  mutate(ln_Cspecific_rate = log(Cspecific_rate)) # log transform mass-specific rate

# Save mdat so we can estimate number of taxa for ratios later...
# saveRDS(mdat, file = "Data/modelParameters/funcGrp_mdat.rds")

# Quick look
mdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = funcGrp)) +
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
  # How does temperature dependence vary across zooplankton func groups for each rate? AND
  # How does temperature dependence vary across rate processes?

# I am also mainly interested in the interactions between at least temp:rate and temp:group



# Fit the models ----

# A complex model with 3 way interactions for temp, funcGrp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * funcGrp + # 3-way interaction
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 
  # Model did not converge....as with the previous size-group analysis, we will simplify the random effects structure first
  summary(m1)
  # primRef slope is closer to zero
  

# Refit m1 as m2 without primRef as random slope
m2 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * funcGrp + # three-way interaction
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

  
# Fit a simpler model without the 3-way interactions, just 2-way interactions and same random effect structure
  # This model does not initially converge with random slope for primRef, so we remove it
m3 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name + # interactions between temp and different rates across all zooplankton
                temp_C * funcGrp +  # between temp and different funcGrps for all rates
                rate_name * funcGrp + # between rates and funcGrp 
                (1 | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes (only slope for taxa, because model did not converge)
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m3)
  plot(sim) # Looks OK, virtually the same as m2
  summary(m3)

  
# Compare models
performance::compare_performance(m2, m3)
# Models are not all mutually nested...we'll just treat AIC as descriptive...
  # 3-way interaction is still best so far (m2)


# Likelihood ratios test of the models
# Refit with ML to test fixed effects and random effect structures
m2_m2 <- update(m2, REML = FALSE)
m3_m3 <- update(m3, REML = FALSE)

anova(m2_m2, m3_m3) # test if the three-way interaction is better than the two-way interaction without random slope for primRef
  # Yep, m2 with 3-way interaction is better

# we will progress with m2
summary(m2)


# Extract slopes using and calculate Q10 for each funcGrp ----
funcGrp_slopes <- emtrends(m2, ~ rate_name * funcGrp, var = "temp_C")
summary(funcGrp_slopes, infer = TRUE) %>% arrange(funcGrp, rate_name) # test whether each slope is different from zero for each zoopGrp
# No significant temp relationships for: (may due to limited data across the temp range?)
  # GelPreds, clearance
  # GelPreds, ingestion
  # GelPreds, growth
  # GelFilter, clearance
  # GelFilter, ingestion
  # GelFilter, growth
  # GelFilter, respiration
emmeans::contrast(funcGrp_slopes, by = "rate_name", method = "pairwise", adjust = "mvt") # pairwise test whether slopes differ significantly across rate types and grps
# Only sig differences for:
  # crustaceans - gelPreds, clearance
  # crustaceans - gelPreds, growth
  # crustaceans - gelPreds, excretion
  # gelPreds - gelFilters, excretion



# Extract intercepts ----
funcGrp_intercepts <- data.frame(emmeans(m2, ~ rate_name * funcGrp, var = "temp_C"))

# Build a parameter dataframe and to estimate ratios
params <- data.frame(funcGrp_slopes) %>% 
  rename(slope = temp_C.trend) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL)) %>% 
  left_join(funcGrp_intercepts) %>% 
  rename(intercept = emmean) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL))


# Save the parameter data for later, we'll use this to estimate ratios of the terms
# saveRDS(params, file = "Data/modelParameters/funcGrpestimates.rds")


# Get n
n_obs <- mdat %>%
  count(rate_name, funcGrp)
n_obs %>% arrange(funcGrp, rate_name)

# Get Q10
funcGrp_slopes_Q10 <- as.data.frame(funcGrp_slopes) %>% 
  # Calculate Q10s
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  # Calculate equivalent activation energies (Ea)
  mutate(refT = 15 + 273.15, # reference temp in Kelvin
         k = 8.617e-5, # Boltzmann's constant as eV/K-1
         R = 8.314 / 1000, # Ideal gas constant as kJ mol-1
         Ea_eV = round(k * log(Q10) * (refT * (refT + 10)) / 10, digits = 2), # Ea as eV
         Ea_kJ.mol =  round((R) * log(Q10) * (refT * (refT + 10)) / 10, digits = 2) # Ea as kJ/mol-1
         ) %>%
  select(- c(k, R, df, refT, asymp.LCL, asymp.UCL)) %>%
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
    geom_point(data = mdat,
               aes(x = temp_C, y = ln_Cspecific_rate, colour = funcGrp),
               alpha = 0.2) +
    geom_ribbon(data = pop_preds,
                aes(x = temp_C, ymin = conf.low, ymax = conf.high, fill = funcGrp),
                alpha = 0.25) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate, colour = funcGrp),
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
    coord_cartesian(xlim = c(-2, 32)) +
    scale_fill_manual(values = grp_cols, labels = c("Crustaceans" , "Gelatinous predators", "Gelatinous filter-feeders")) +
    scale_colour_manual(values = grp_cols, labels = c("Crustaceans" , "Gelatinous predators", "Gelatinous filter-feeders")) +
    labs(x = "Temp (°C)",
         y = "ln(Carbon-mass specific rate)",
         fill = "Functional group",
         colour = "Functional group") +
    theme(
      strip.background = element_rect(fill = "whitesmoke", colour = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.position = "top"
    )
}


# Plot it
tempPlot <- PlotLMM(m2)
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
               "Clearance"   = "bold(Clearance~(ml~mgC^-1~h^-1))",
               "Ingestion"   = "bold(Ingestion~(mgC~mgC^-1~h^-1))",
               "Growth"      = "bold(Growth~(mgC~mgC^-1~h^-1))",
               "Respiration" = "bold(Respiration~(µlO[2]~mgC^-1~h^-1))",
               "Excretion"   = "bold(Excretion~('mgN-NH'[4]^'+'~mgC^-1~h^-1))"
             ), 
             label_parsed),
             ncol = 1) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "Functional group",
       y = bquote(bold("Carbon-mass specific Q"[10])),
       colour = "Functional group") +
  theme(
    strip.background = element_rect(fill = "whitesmoke", colour = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )
funcGQ10plot
# confidence intervals are not good to visualise, we will break the y-axes to better compared estimates

library(ggbreak) # to break y-axis on massive Q10 variance
funcG10plot_break <- funcGQ10plot +
  geom_text(data = filter(funcGrp_slopes_Q10, 
                          (rate_name == "Ingestion" & Q10_upr > 5) |
                          (rate_name == "Respiration" & Q10_upr > 5) |
                          (rate_name == "Growth" & Q10_upr > 5) |
                          (rate_name == "Clearance" & Q10_upr > 5)),
            aes(x = funcGrp, y = 4.85, label = paste0("↑ ", round(Q10_upr, 0))),
            colour = "grey40", size = 4,
            nudge_x = 0.3) +
  coord_cartesian(ylim = c(-1, 5))
funcG10plot_break

tempPlot + funcG10plot_break

# Save it
ggsave("Output/Figure4/Figure4_tempPlot.pdf", tempPlot, width = 75, height = 240, units = "mm", dpi = 300)
ggsave("Output/Figure4/Figure4_Q10Plot.pdf", funcG10plot_break, width = 75, height = 240, units = "mm", dpi = 300)

