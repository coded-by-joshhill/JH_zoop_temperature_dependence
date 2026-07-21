# Calculating zoopGrps Q10s
# Josh Hill
# 5/05/2026



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
usedat <- rbind(cleardat, ingdat, grwdat, respdat, excredat) %>% 
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
    # wont be able to get estimates for amphipods, mysids
  # Excretion:
    # seems fine




# Prep data for modelling ----
mdat <- usedat %>%
  group_by(zoopGrp, rate_name) %>% 
  filter(n() >= 15, # Exclude zoopGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  mutate(ln_Cspecific_rate = log(Cspecific_rate)) # log transform mass-specific rate

# Save mdat so we can estimate number of taxa for ratios later...
# saveRDS(mdat, file = "Data/modelParameters/zoopGrp_mdat.rds")


# Quick look
mdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = zoopGrp)) +
  facet_wrap(~ rate_name, scale = "free")
  # Looks pretty tidy. Some clear relationships here
  # Note: 
  # clearance is ml/mgC/hr
  # ingestion is mgC/mgC/hr 
  # growth is mgC/mgC/hr
  # respiration is uLO2/mgC/hr
  # excretion is mgN-NH4+/mgC/hr
summary(mdat)


# My main question here is...
  # How does temperature dependence vary across zooplankton groups for each rate? AND
  # How does temperature dependence vary across rate processes?

# I am also mainly interested in the interactions between at least temp:rate and temp:group



# Fit the models ----

# A complex model with 3 way interactions for temp, zoopGrp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * zoopGrp + # 3-way interaction
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat)
# Columns are being dropped from rank-deficient conditional model - probably not enough data for a few groups across all rates
# Lets check the counts for each fixed effect before progressing
  
# Are all zoopGrp levels present across all rate types?
mdat %>%
  group_by(zoopGrp) %>%
  summarise(rates_n = n_distinct(rate_name),
            n = n())
# nope, which I suspected was the issue (some are rates_n = 1 etc..)

# Two-way counts to see where we are missing data
mdat %>%
  count(rate_name, zoopGrp) %>%
  pivot_wider(names_from = zoopGrp, values_from = n, values_fill = 0)
  # rate to group combinations are clearly lacking data to estimate coefficients in a single model

# Because there are many zeros across the groups and rates...I will fit per-rate models to estimate temperature dependence for available groups

  

# Per-rate models ----
m2 <- mdat %>% # using my model data
  group_by(rate_name) %>% # group by each rate 
  nest() %>% # nest each rate into their own tibbles
  mutate(
    fit = map(data, \(data) glmmTMB( # then fit my mixed effect model
      ln_Cspecific_rate ~ # with mass-specific rates as a function off...
        temp_C * zoopGrp + # a two-way interaction between temp and zoop taxonomic group
        (1 | primRef) + (1 | taxa), # with random intercepts only - there is insufficient data for random slopes
      data = data))) # using the nested rate-specific data


# Extract each model using the name
m_clearance <- m2$fit[m2$rate_name == "Clearance"][[1]]
m_ingestion <- m2$fit[m2$rate_name == "Ingestion"][[1]]
m_growth <- m2$fit[m2$rate_name == "Growth"][[1]]
m_respiration <- m2$fit[m2$rate_name == "Respiration"][[1]]
m_excretion <- m2$fit[m2$rate_name == "Excretion"][[1]]


# Get the data too
d_clearance <- m2$data[m2$rate_name == "Clearance"][[1]]
d_ingestion <- m2$data[m2$rate_name == "Ingestion"][[1]]
d_growth <- m2$data[m2$rate_name == "Growth"][[1]]
d_respiration <- m2$data[m2$rate_name == "Respiration"][[1]]
d_excretion <- m2$data[m2$rate_name == "Excretion"][[1]]


# Check diagnostics ----
# Clearance
sim <- simulateResiduals(m_clearance)
plot(sim) # doesn't look great
# Plot without tests
par(mfrow = c(1, 2))
plotQQunif(sim, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
plotResiduals(sim)
par(mfrow = c(1, 1)) # Reset back to normal

# we'll accept this as a limitation instead of over fitting this data
summary(m_clearance)


# Ingestion
sim <- simulateResiduals(m_ingestion)
plot(sim) # basically have the same issue here with ingestion, though QQ looks better
# Plot without tests
par(mfrow = c(1, 2))
plotQQunif(sim, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
plotResiduals(sim)
par(mfrow = c(1, 1)) # Reset back to normal
summary(m_ingestion)


# Growth
sim <- simulateResiduals(m_growth)
plot(sim) # looks fine
# Plot without tests
par(mfrow = c(1, 2))
plotQQunif(sim, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
plotResiduals(sim)
par(mfrow = c(1, 1)) # Reset back to normal

summary(m_growth)


# Respiration
sim <- simulateResiduals(m_respiration)
plot(sim) # looks pretty good overall
# Plot without tests
par(mfrow = c(1, 2))
plotQQunif(sim, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
plotResiduals(sim)
par(mfrow = c(1, 1)) # Reset back to normal
summary(m_respiration)


# Excretion
sim <- simulateResiduals(m_excretion)
plot(sim) # looks pretty good overall
# Plot without tests
par(mfrow = c(1, 2))
plotQQunif(sim, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
plotResiduals(sim)
par(mfrow = c(1, 1)) # Reset back to normal
summary(m_excretion)



# Extract coefficients and calculate Q10 for each zoopGrp ----
## Clearance ----
clear_slopes <- emtrends(m_clearance, ~zoopGrp, var = "temp_C")
clear_intercepts <- data.frame(emmeans(m_clearance, ~zoopGrp, var = "temp_C"))
summary(clear_slopes, infer = TRUE) # test whether each slope is different from zero for each zoopGrp
# Only copepod clearance rates are sig diff to zero
emmeans::contrast(clear_slopes, method = "pairwise", adjust = "mvt") # pairwise test whether slopes differ significantly across groups
  # No sig differences between any groups for clearance rate

# Build a parameter dataframe to estimate ratios
clearance_params <- data.frame(clear_slopes) %>% 
  rename(slope = temp_C.trend) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL)) %>% 
  left_join(clear_intercepts) %>% 
  rename(intercept = emmean) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL))

# Save the parameter data for later
saveRDS(clearance_params, file = "Data/modelParameters/ZGrpClearanceEstimates.rds")

# Get n
n_obs_clearance <- d_clearance %>%
  count(zoopGrp)

# Get Q10s
clearanceQ10 <- as.data.frame(clear_slopes) %>% 
  # Calculate Q10s
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  # Calculate equivalent activation energies (Ea)
  mutate(refT = 15 + 273.15, # reference temp in Kelvin
         k = 8.617e-5, # Boltzmann's constant as eV/K-1
         R = 8.314 / 1000, # Ideal gas constant as kJ mol-1
         Ea_eV = round(k * log(Q10) * (refT * (refT + 10)) / 10, digits = 2), # Ea as eV
         Ea_kJ.mol =  round((R) * log(Q10) * (refT * (refT + 10)) / 10), digits = 2) %>% # Ea as kJ/mol-1
  select(- c(k, R, df, refT, asymp.LCL, asymp.UCL)) %>%
  left_join(n_obs_clearance, by = "zoopGrp")
clearanceQ10
# Euphausiids probably lacking data... Q10=0.58 doesn't make sense



## Ingestion ----
ingest_slopes <- emtrends(m_ingestion, ~zoopGrp, var = "temp_C")
ingest_intercepts <- data.frame(emmeans(m_ingestion, ~zoopGrp, var = "temp_C"))
summary(ingest_slopes, infer = TRUE) # test whether each slope is different from zero for each zoopGrp
# Copepods and euphausiids are sig diff...thaliaceans marginally...
emmeans::contrast(ingest_slopes, method = "pairwise", adjust = "mvt") # pairwise test whether slopes differ significantly across groups
  # No sig differences between any groups for ingestion rate

# Build a parameter dataframe to estimate ratios
ingest_params <- data.frame(ingest_slopes) %>% 
  rename(slope = temp_C.trend) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL)) %>% 
  left_join(ingest_intercepts) %>% 
  rename(intercept = emmean) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL))

# Save the parameter data for later
saveRDS(ingest_params, file = "Data/modelParameters/ZGrpIngestionEstimates.rds")

# Get n
n_obs_ingestion <- d_ingestion %>%
  count(zoopGrp)

# Get Q10s
ingestionQ10 <- as.data.frame(ingest_slopes) %>% 
  # Calculate Q10s
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  # Calculate equivalent activation energies (Ea)
  mutate(refT = 15 + 273.15, # reference temp in Kelvin
         k = 8.617e-5, # Boltzmann's constant as eV/K-1
         R = 8.314 / 1000, # Ideal gas constant as kJ mol-1
         Ea_eV = round(k * log(Q10) * (refT * (refT + 10)) / 10, digits = 2), # Ea as eV
         Ea_kJ.mol =  round((R) * log(Q10) * (refT * (refT + 10)) / 10), digits = 2) %>% # Ea as kJ/mol-1
  select(- c(k, R, df, refT, asymp.LCL, asymp.UCL)) %>%
  left_join(n_obs_ingestion, by = "zoopGrp")
ingestionQ10
# Again, probably lacking data for euphausiids...Q10=4.16 seems somewhat extreme



## Growth ----
growth_slopes <- emtrends(m_growth, ~zoopGrp, var = "temp_C")
growth_intercepts <- data.frame(emmeans(m_growth, ~zoopGrp, var = "temp_C"))
summary(growth_slopes, infer = TRUE)# test whether each slope is different from zero for each zoopGrp
# Sig differences for amphipods, copepods, euphausiids and thaliaceans
emmeans::contrast(growth_slopes, method = "pairwise", adjust = "mvt") # pairwise test whether slopes differ significantly across groups
# Sig differences for:  
  # Cnidarians - Amphipods    -0.12759 0.0296 Inf  -4.310  0.0003
  # Cnidarians - Copepods     -0.10940 0.0261 Inf  -4.198  0.0005

# Build a parameter dataframe to estimate ratios
growth_params <- data.frame(growth_slopes) %>% 
  rename(slope = temp_C.trend) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL)) %>% 
  left_join(growth_intercepts) %>% 
  rename(intercept = emmean) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL))

# Save the parameter data for later
saveRDS(growth_params, file = "Data/modelParameters/ZGrpGrowthEstimates.rds")

# Get n
n_obs_growth <- d_growth %>%
  count(zoopGrp)

# Get Q10s
growthQ10 <- as.data.frame(growth_slopes) %>% 
  # Calculate Q10s
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  # Calculate equivalent activation energies (Ea)
  mutate(refT = 15 + 273.15, # reference temp in Kelvin
         k = 8.617e-5, # Boltzmann's constant as eV/K-1
         R = 8.314 / 1000, # Ideal gas constant as kJ mol-1
         Ea_eV = round(k * log(Q10) * (refT * (refT + 10)) / 10, digits = 2), # Ea as eV
         Ea_kJ.mol =  round((R) * log(Q10) * (refT * (refT + 10)) / 10), digits = 2) %>% # Ea as kJ/mol-1
  select(- c(k, R, df, refT, asymp.LCL, asymp.UCL)) %>%
  left_join(n_obs_growth, by = "zoopGrp")
growthQ10
# These results look pretty good! Interesting that cnidarians growth rate scales negatively with temperature



## Respiration ----
resp_slopes <- emtrends(m_respiration, ~zoopGrp, var = "temp_C")
resp_intercepts <- data.frame(emmeans(m_respiration, ~zoopGrp, var = "temp_C"))
summary(resp_slopes, infer = TRUE) # test whether each slope is different from zero for each zoopGrp
# All grps sig diff to zero
emmeans::contrast(resp_slopes, method = "pairwise", adjust = "mvt") # pairwise test whether slopes differ significantly across groups
# Sig differences for:
# Ctenophores - Cnidarians    0.1112 0.0232 Inf   4.801 <0.0001
# Ctenophores - Euphausiids   0.0899 0.0246 Inf   3.658  0.0020
# Copepods - Euphausiids      0.0420 0.0134 Inf   3.136  0.0115

# Build a paramater dataframe to estimate ratios
respiration_params <- data.frame(resp_slopes) %>% 
  rename(slope = temp_C.trend) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL)) %>% 
  left_join(resp_intercepts) %>% 
  rename(intercept = emmean) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL))

# Save the parameter data for later
saveRDS(respiration_params, file = "Data/modelParameters/ZGrpRespirationEstimates.rds")

# Get n
n_obs_respiration <- d_respiration %>%
  count(zoopGrp)

# Get Q10s
respQ10 <- as.data.frame(resp_slopes) %>% 
  # Calculate Q10s
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  # Calculate equivalent activation energies (Ea)
  mutate(refT = 15 + 273.15, # reference temp in Kelvin
         k = 8.617e-5, # Boltzmann's constant as eV/K-1
         R = 8.314 / 1000, # Ideal gas constant as kJ mol-1
         Ea_eV = round(k * log(Q10) * (refT * (refT + 10)) / 10, digits = 2), # Ea as eV
         Ea_kJ.mol =  round((R) * log(Q10) * (refT * (refT + 10)) / 10), digits = 2) %>% # Ea as kJ/mol-1
  select(- c(k, R, df, refT, asymp.LCL, asymp.UCL)) %>%
  left_join(n_obs_respiration, by = "zoopGrp")
respQ10
# These results look great too



## Excretion ----
excr_slopes <- emtrends(m_excretion, ~zoopGrp, var = "temp_C")
excr_intercepts <- data.frame(emmeans(m_excretion, ~zoopGrp, var = "temp_C"))
summary(excr_slopes, infer = TRUE) # test whether each slope is different from zero for each zoopGrp
# All sig diff except for mysiids and appendicularians
emmeans::contrast(excr_slopes, method = "pairwise", adjust = "mvt") # pairwise test whether slopes differ significantly across groups
# Sig differences for:
# Ctenophores - Mysids             0.141232 0.0455 Inf   3.103  0.0413
# Ctenophores - Thaliaceans        0.094533 0.0294 Inf   3.217  0.0287

# Build a paramater dataframe to estimate ratios
excr_params <- data.frame(excr_slopes) %>% 
  rename(slope = temp_C.trend) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL)) %>% 
  left_join(excr_intercepts) %>% 
  rename(intercept = emmean) %>% 
  select(-c(SE, df, asymp.LCL, asymp.UCL))

# Save the parameter data for later
saveRDS(excr_params, file = "Data/modelParameters/ZGrpExcretionEstimates.rds")

# Get n
n_obs_excretion <- d_excretion %>%
  count(zoopGrp)

# Get Q10s
excretionQ10 <- as.data.frame(excr_slopes) %>% 
  # Calculate Q10s
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  # Calculate equivalent activation energies (Ea)
  mutate(refT = 15 + 273.15, # reference temp in Kelvin
         k = 8.617e-5, # Boltzmann's constant as eV/K-1
         R = 8.314 / 1000, # Ideal gas constant as kJ mol-1
         Ea_eV = round(k * log(Q10) * (refT * (refT + 10)) / 10, digits = 2), # Ea as eV
         Ea_kJ.mol =  round((R) * log(Q10) * (refT * (refT + 10)) / 10), digits = 2) %>% # Ea as kJ/mol-1
  select(- c(k, R, df, refT, asymp.LCL, asymp.UCL)) %>%
  left_join(n_obs_excretion, by = "zoopGrp")
excretionQ10


# Reprint the results
clearanceQ10
ingestionQ10
growthQ10
respQ10
excretionQ10


# Prep for plotting ----

# Define group colours
grp_cols <- c(
  "Ctenophores"      = "#f46d43",
  "Cnidarians"       = "#fdae61",
  "Chaetognaths"     = "#fee090",
  "Amphipods"        = "#bdd7e7",
  "Copepods"         = "#6baed6",  
  "Decapods"         = "#3182bd",
  "Euphausiids"      = "#08519c",
  "Mysids"           = "#08306b",
  "Appendicularians" = "#a1d99b",
  "Thaliaceans"      = "#31a354")



# Set min and max temp
minTempC <- min(mdat$temp_C)
maxTempC <- max(mdat$temp_C)



# FUNCTIOIN TO PLOT MODEL ----
PlotLMM = function(model, data){
  temp_seq <- seq(minTempC, maxTempC, length.out = 100)
  
  # Build newdata grid manually
  newdat <- expand.grid(
    temp_C    = temp_seq,
    zoopGrp   = unique(data$zoopGrp)) # levels of zoopGrp
  

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
    geom_point(data = data,
               aes(x = temp_C, y = ln_Cspecific_rate, colour = zoopGrp),
               alpha = 0.2) +
    geom_ribbon(data = pop_preds,
                aes(x = temp_C, ymin = conf.low, ymax = conf.high, fill = zoopGrp),
                alpha = 0.25) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate, colour = zoopGrp),
              linewidth = 1) +
    coord_cartesian(xlim = c(-2, 32)) +
    labs(x = "Temp (°C)",
         fill = "Taxonomic group",
         colour = "Taxonomic group") +
    scale_fill_manual(values = grp_cols) +
    scale_colour_manual(values = grp_cols) +
    scale_x_continuous(breaks = c(0, 15, 30)) +
    theme(
      strip.background = element_rect(fill = "whitesmoke", colour = "black"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold")
    )
}



# Plot it
# Clearance
clearPlot <- PlotLMM(m_clearance, d_clearance) + 
  labs(y = expression(atop(bold("ln (Clearance)"),
                           bold("(ml mgC"^-1*" h"^-1*")")))) +
  theme(legend.position = "none") +
  facet_wrap(~zoopGrp, nrow = 1)
clearPlot


# Ingestion
ingPlot <- PlotLMM(m_ingestion, d_ingestion) + 
  labs(y = expression(atop(bold("ln (Ingestion)"),
                           bold("(mgC mgC"^-1*" h"^-1*")")))) +
  theme(legend.position = "none") +
  facet_wrap(~zoopGrp, nrow = 1)
ingPlot


# Growth
growPlot <- PlotLMM(m_growth, d_growth) + 
  labs(y = expression(atop(bold("ln (Growth)"),
                           bold("(mgC mgC"^-1*" h"^-1*")")))) +
  theme(legend.position = "none") +
  facet_wrap(~zoopGrp, nrow = 2)
growPlot


# Respiration
respPlot <- PlotLMM(m_respiration, d_respiration) + 
  labs(y = expression(atop(bold("ln (Respiration)"),
                           bold("(" * mu * "lO"[2] * " mgC"^-1 * " h"^-1 * ")")))) +
  theme(legend.position = "none") +
  facet_wrap(~zoopGrp, nrow = 1)
respPlot 


# Excretion
excrPlot <- PlotLMM(m_excretion, d_excretion) + 
  labs(y = expression(atop(bold("ln (Excretion)"),
                           bold('(mgN-NH'[4]^"+" * " mgC"^-1 * " h"^-1 * ")")))) +
  theme(legend.position = "none") +
  facet_wrap(~zoopGrp, nrow = 2)
excrPlot



tempPlots <- clearPlot / ingPlot / growPlot / respPlot / excrPlot + 
  plot_layout(
    axis_titles = "collect_x",
    heights = c(1, 1, 2, 1, 2, 0.5), # adjust relative heights of each subplot
  )
tempPlots
# Looks fine...not going to include a legend, there is no space and these colours are just hues of higher level group analyses and it's pretty obvious what they represent across the processes.



# Plot zoopGrp Q10s ----
# Clearance
clearQ10plot <- ggplot() +
  geom_errorbar(data = clearanceQ10, 
                aes(x = zoopGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = zoopGrp),
                width = .05,
                linewidth = 1) +
  geom_point(data = clearanceQ10, aes(x = zoopGrp, y = Q10),
             size = 3,
             colour = "black") +
  geom_text(data = clearanceQ10,
            aes(x = zoopGrp, y = Q10, label = sprintf("%.2f", Q10)),
            nudge_x = 0.25,
            nudge_y = 0.0) +
  scale_fill_manual(values = grp_cols) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "Taxonomic group",
       y = bquote(bold("Carbon-mass specific clearance Q"[10])),
       colour = "Taxonomic group") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none") +
  coord_cartesian(ylim = c(0, 5.5))
clearQ10plot


# Ingestion
ingQ10plot <- ggplot() +
  geom_errorbar(data = ingestionQ10, 
                aes(x = zoopGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = zoopGrp),
                width = .05,
                linewidth = 1) +
  geom_point(data = ingestionQ10, aes(x = zoopGrp, y = Q10),
             size = 3,
             colour = "black") +
  geom_text(data = ingestionQ10,
            aes(x = zoopGrp, y = Q10, label = sprintf("%.2f", Q10)),
            nudge_x = 0.20,
            nudge_y = 0.0) +
  scale_fill_manual(values = grp_cols) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "Taxonomic group",
       y = bquote(bold("Carbon-mass specific ingestion Q"[10])),
       colour = "Taxonomic group") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none") +
  coord_cartesian(ylim = c(0, 5.5))
ingQ10plot

library(ggbreak) # to break y-axis on massive Q10 CI
ingQ10plot_break <- ingQ10plot +
  geom_text(data = filter(ingestionQ10, 
                          (zoopGrp == "Copepods" & Q10_upr > 2) |
                            (zoopGrp == "Euphausiids" & Q10_upr > 2)),
            aes(x = zoopGrp, y = 5.5, label = paste0("↑ ", round(Q10_upr, 1))),
            colour = "grey40", size = 4,
            nudge_x = 0.3) +
  coord_cartesian(ylim = c(0, 5.5))
ingQ10plot_break


# Growth
growQ10plot <- ggplot() +
  geom_errorbar(data = growthQ10, 
                aes(x = zoopGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = zoopGrp),
                width = .09,
                linewidth = 1) +
  geom_point(data = growthQ10, aes(x = zoopGrp, y = Q10),
             size = 3,
             colour = "black") +
  geom_text(data = growthQ10,
            aes(x = zoopGrp, y = Q10, label = sprintf("%.2f", Q10)),
            nudge_x = 0.35,
            nudge_y = 0.0) +
  scale_fill_manual(values = grp_cols) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "Taxonomic group",
       y = bquote(bold("Carbon-mass specific growth Q"[10])),
       colour = "Taxonomic group") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none",
    legend.title = element_text(size = 10, face = "bold")) +
  scale_x_discrete(expand = expansion(add = 0.8))  # adds 0.8 units of padding on each side
growQ10plot

growQ10plot_break <- growQ10plot +
  geom_text(data = filter(growthQ10, 
                          (zoopGrp == "Decapods" & Q10_upr > 2) |
                          (zoopGrp == "Euphausiids" & Q10_upr > 2) |
                          (zoopGrp == "Thaliaceans" & Q10_upr > 2)),
            aes(x = zoopGrp, y = 5.5, label = paste0("↑ ", round(Q10_upr, 1))),
            colour = "grey40", size = 4,
            nudge_x = 0.3) +
  coord_cartesian(ylim = c(0, 5.5))
growQ10plot_break


# Respiration
respQ10plot <- ggplot() +
  geom_errorbar(data = respQ10, 
                aes(x = zoopGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = zoopGrp),
                width = .07,
                linewidth = 1) +
  geom_point(data = respQ10, aes(x = zoopGrp, y = Q10),
             size = 3,
             colour = "black") +
  geom_text(data = respQ10,
            aes(x = zoopGrp, y = Q10, label = sprintf("%.2f", Q10)),
            nudge_x = 0.3,
            nudge_y = 0.0) +
  scale_fill_manual(values = grp_cols) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "Taxonomic group",
       y = bquote(bold("Carbon-mass specific respiration Q"[10])),
       colour = "Taxonomic group") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none")
respQ10plot

# Break the CIs
respQ10plot_break <- respQ10plot +
  geom_text(data = filter(respQ10, 
                          (zoopGrp == "Ctenophores" & Q10_upr > 5) |
                            (zoopGrp == "Thaliaceans" & Q10_upr > 5)),
            aes(x = zoopGrp, y = 5, label = paste0("↑ ", round(Q10_upr, 1))),
            colour = "grey40", size = 4,
            nudge_x = 0.3) +
  coord_cartesian(ylim = c(0, 5.5))
respQ10plot_break


# Excretion
excrQ10plot <- ggplot() +
  geom_errorbar(data = excretionQ10, 
                aes(x = zoopGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = zoopGrp),
                width = .07,
                linewidth = 1) +
  geom_point(data = excretionQ10, aes(x = zoopGrp, y = Q10),
             size = 3,
             colour = "black") +
  geom_text(data = excretionQ10,
            aes(x = zoopGrp, y = Q10, label = sprintf("%.2f", Q10)),
            nudge_x = 0.35,
            nudge_y = 0.0) +
  scale_fill_manual(values = grp_cols) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "Taxonomic group",
       y = bquote(bold("Carbon-mass specific excretion Q"[10])),
       colour = "Taxonomic group") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none")
excrQ10plot

# Break the CIs
excrQ10plot_break <- excrQ10plot +
  geom_text(data = filter(excretionQ10, 
                          (zoopGrp == "Ctenophores" & Q10_upr > 5)),
            aes(x = zoopGrp, y = 5.5, label = paste0("↑ ", round(Q10_upr, 1))),
            colour = "grey40", size = 4,
            nudge_x = 0.3) +
  coord_cartesian(ylim = c(0, 5.5))
excrQ10plot_break


tempPlots
# Combine with patchwork
Q10Plots <- (clearQ10plot + ingQ10plot_break) / (growQ10plot_break + respQ10plot_break) + excrQ10plot_break +
  plot_layout(guides = "collect", 
              axis_titles = "collect_x") & 
  theme(legend.position = "none") 
  
Q10Plots
# Nice to see them together but will add post hoc results on with post digitising

# Can't save together properly so will do it separately
clearQ10plot
ingQ10plot_break
growQ10plot_break
respQ10plot_break
excrQ10plot_break

ggsave("Output/Figure5/Figure5_clearQ10.pdf", clearQ10plot, width = 90, height = 90, units = "mm", dpi = 300)
ggsave("Output/Figure5/Figure5_ingQ10.pdf", ingQ10plot_break, width = 90, height = 90, units = "mm", dpi = 300)
ggsave("Output/Figure5/Figure5_growQ10.pdf", growQ10plot_break, width = 90, height = 90, units = "mm", dpi = 300)
ggsave("Output/Figure5/Figure5_respQ10.pdf", respQ10plot_break, width = 90, height = 90, units = "mm", dpi = 300)
ggsave("Output/Figure5/Figure5_excrQ10.pdf", excrQ10plot_break, width = 175, height = 90, units = "mm", dpi = 300)


# # Save the legend from this plot for post digitising
# legend <- growQ10plot_break + excrQ10plot_break +
#   plot_layout(guides = "collect") & theme(legend.position = "right")
# legend
# # ggsave("Output/Figure5/Figure5_legend.pdf", legend, width = 180, height = 180, units = "mm", dpi = 300)


# Save temperature plots
ggsave("Output/Figure5/FigureS1_tempPlot.pdf", tempPlots, width = 170, height = 220, units = "mm", dpi = 300)
ggsave("Output/Figure5/FigureS1_tempPlot.png", tempPlots, width = 170, height = 220, units = "mm", dpi = 300)


