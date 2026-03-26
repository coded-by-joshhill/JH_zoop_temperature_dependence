# Calculating zoopGrps Q10s
# Josh Hill
# 13/03/2026



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
    # wont be able to get estimates for amphipods, decapods, mysids, annelids, molluscs


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



# Prep data for modelling ----
mdat <- usedat %>%
  group_by(zoopGrp, rate_name) %>% 
  filter(n() >= 15, # Exclude zoopGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
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

# I am also mainly interested in the interactions between at least temp:rate and temp:group



# Fit the models ----

# A complex model with 3 way interactions for temp, zoopGrp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ 
                temp_C * rate_name * zoopGrp + # 3-way interaction
                (temp_C | primRef) + (1 | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat)

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim) # looks pretty good
  summary(m1)

# Columns are being dropped from rank-deficient conditional model - probably not enough data for a few groups across all rates
# Lets check the counts for each fixed effect before progressing
  
  # Are all zoopGrp levels present across all rate types?
  mdat %>%
    group_by(zoopGrp) %>%
    summarise(rates_n = n_distinct(rate_name),
              n = n())
  # nope, which I suspected was the issue


# Two-way counts
mdat %>%
  count(rate_name, zoopGrp) %>%
  pivot_wider(names_from = zoopGrp, values_from = n, values_fill = 0)
  # 12 rate to group combinations lacking data to estimate coefficients

# Because there are many zeros across the groups and rates...I will fit per-rate models to estimate temperature dependence for available groups

  

# Per-rate models ----
m2 <- mdat %>% # using my model data
  group_by(rate_name) %>% # group by each rate 
  nest() %>% # nest each rate into their own tibbles
  mutate(
    fit = map(data, \(data) glmmTMB(
      ln_Cspecific_rate ~ 
        temp_C * zoopGrp + # two-way interaction between temp and zooplankton group
        (1 | primRef) + (1 | taxa), # random intercepts only - insufficient data for random slopes
      dispformula = ~1,
      data = data # use the nested rate-specific data slice
    )))


# Extract each model using the name
m_clearance   <- m2$fit[m2$rate_name == "Clearance"][[1]]
m_ingestion   <- m2$fit[m2$rate_name == "Ingestion"][[1]]
m_growth      <- m2$fit[m2$rate_name == "Growth"][[1]]
m_respiration <- m2$fit[m2$rate_name == "Respiration"][[1]]

# Get the data too
d_clearance   <- m2$data[m2$rate_name == "Clearance"][[1]]
d_ingestion   <- m2$data[m2$rate_name == "Ingestion"][[1]]
d_growth      <- m2$data[m2$rate_name == "Growth"][[1]]
d_respiration <- m2$data[m2$rate_name == "Respiration"][[1]]


# Check diagnostics
# Clearance
sim <- simulateResiduals(m_clearance)
plot(sim) # doesn't look great
# Let's update the model with a dispersion formula to try improve homoscedasticity
m_clearance2 <- update(m_clearance, dispformula = ~zoopGrp, data = d_clearance)
sim <- simulateResiduals(m_clearance2)
plot(sim) 
# doesn't seem to have improved anything, we will stick with the first model and accept this as a limitation instead of over fitting this data
summary(m_clearance)


# Ingestion
sim <- simulateResiduals(m_ingestion)
plot(sim) # basically have the same issue here with ingestion
m_ingestion2 <- update(m_ingestion, dispformula = ~zoopGrp, data = d_ingestion) # update with dispersion formula across each group
sim <- simulateResiduals(m_ingestion2)
plot(sim) # also no difference...stick with original model
summary(m_ingestion)


# Growth
sim <- simulateResiduals(m_growth)
plot(sim) # residuals are sitting above the line on the QQ plot, probably an artifact of slightly left-skewed data
  # the residuals vs predicted plot looks OK
summary(m_growth)


# Respiration
sim <- simulateResiduals(m_respiration)
plot(sim) # looks pretty good overall
summary(m_respiration)


# Extract slopes using and calculate Q10 for each zoopGrp ----
# Clearance
clear_slopes <- emtrends(m_clearance, ~zoopGrp, var = "temp_C")
summary(clear_slopes, infer = TRUE) # test whether each slope is different from zero for each zoopGrp
pairs(clear_slopes)
  # No sig differences between any groups for clearance rate

# Get n
n_obs_clearance <- d_clearance %>%
  count(zoopGrp)

clearanceQ10 <- as.data.frame(clear_slopes) %>% 
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  left_join(n_obs_clearance, by = "zoopGrp")
clearanceQ10


# Ingestion
ingest_slopes <- emtrends(m_ingestion, ~zoopGrp, var = "temp_C")
summary(ingest_slopes, infer = TRUE) # test whether each slope is different from zero for each zoopGrp
pairs(ingest_slopes)
  # No sig differences between any groups for ingestion rate

# Get n
n_obs_ingestion <- d_ingestion %>%
  count(zoopGrp)

ingestionQ10 <- as.data.frame(ingest_slopes) %>% 
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  left_join(n_obs_ingestion, by = "zoopGrp")
ingestionQ10


# Growth
growth_slopes <- emtrends(m_growth, ~zoopGrp, var = "temp_C")
summary(growth_slopes, infer = TRUE)# test whether each slope is different from zero for each zoopGrp
pairs(growth_slopes)
  # ctenophores - copepods
  # cnidarians - copepods
  # chaetognaths - copepods

# Get n
n_obs_growth <- d_growth %>%
  count(zoopGrp)

growthQ10 <- as.data.frame(growth_slopes) %>% 
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  left_join(n_obs_growth, by = "zoopGrp")
growthQ10


# Respiration
resp_slopes <- emtrends(m_respiration, ~zoopGrp, var = "temp_C")
summary(resp_slopes, infer = TRUE) # test whether each slope is different from zero for each zoopGrp
pairs(resp_slopes)
  # ctenophores - cnidarians
  # cnidarians- copepods
  # cnidarians - euphausiids

# Get n
n_obs_respiration <- d_respiration %>%
  count(zoopGrp)

respQ10 <- as.data.frame(resp_slopes) %>% 
  mutate(Q10 = round(exp(10* temp_C.trend), digits = 2),
         Q10_lwr = round(exp(10 * asymp.LCL), digits = 2),
         Q10_upr = round(exp(10 * asymp.UCL), digits = 2)) %>% 
  left_join(n_obs_respiration, by = "zoopGrp")
respQ10

clearanceQ10
ingestionQ10
growthQ10
respQ10



# Prep for plotting ----

# Define group colours
grp_cols <- c(
  "Ctenophores"  = "#E69F00",  # orange
  "Cnidarians"   = "#56B4E9",  # sky blue
  "Chaetognaths" = "#009E73",  # green
  "Amphipods"    = "#F0E442",  # yellow
  "Copepods"     = "#0072B2",  # dark blue
  "Decapods"     = "#D55E00",  # vermillion
  "Euphausiids"  = "#CC79A7",  # pink
  "Thaliaceans"  = "#000000")   # black


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
    geom_ribbon(data = pop_preds,
                aes(x = temp_C, ymin = conf.low, ymax = conf.high, fill = zoopGrp),
                alpha = 0.25) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate, colour = zoopGrp),
              linewidth = 1) +
    geom_point(data = data,
               aes(x = temp_C, y = ln_Cspecific_rate, colour = zoopGrp),
               alpha = 0.2) +
    coord_cartesian(xlim = c(-2, 32)) +
    labs(x = "Temp (°C)",
         fill = "Taxonomic group",
         colour = "Taxonomic group") +
    scale_fill_manual(values = grp_cols) +
    scale_colour_manual(values = grp_cols) +
    theme_bw() +
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
  labs(y = expression(atop(bold("ln (Clearance rate)"),
                           bold("(ml mgC"^-1*" h"^-1*")")))) +
  theme(legend.position = "none") +
  facet_wrap(~zoopGrp, nrow = 1)
clearPlot


# Ingestion
ingPlot <- PlotLMM(m_ingestion, d_ingestion) + 
  labs(y = expression(atop(bold("ln (Ingestion rate)"),
                           bold("(mgC mgC"^-1*" h"^-1*")")))) +
  theme(legend.position = "none") +
  facet_wrap(~zoopGrp, nrow = 1)
ingPlot


# Growth
growPlot <- PlotLMM(m_growth, d_growth) + 
  labs(y = expression(atop(bold("ln (Growth rate)"),
                           bold("(mgC mgC"^-1*" h"^-1*")")))) +
  facet_wrap(~zoopGrp, nrow = 2)
growPlot


# Respiration
respPlot <- PlotLMM(m_respiration, d_respiration) + 
  labs(y = expression(atop(bold("Respiration rate"),
                           bold("(" * mu * "lO"[2] * " mgC"^-1 * " h"^-1 * ")")))) +
  theme(legend.position = "none") +
  facet_wrap(~zoopGrp, nrow = 1)
respPlot 

tempPlots <- clearPlot / ingPlot / growPlot / respPlot + 
  plot_layout(
    axis_titles = "collect_x",
    guides = "collect", heights = c(1, 1, 2, 1)) # adjust relative heights of each subplot
tempPlots




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
    legend.position = "none")
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
    legend.position = "none")
ingQ10plot


# Growth
growQ10plot <- ggplot() +
  geom_errorbar(data = growthQ10, 
                aes(x = zoopGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = zoopGrp),
                width = .05,
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
    axis.text.x = element_text(angle = 25, hjust = 1))
growQ10plot


# Respiration
respQ10plot <- ggplot() +
  geom_errorbar(data = respQ10, 
                aes(x = zoopGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = zoopGrp),
                width = .05,
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


Q10Plots <- 
  clearQ10plot +
  ingQ10plot +
  growQ10plot +
  respQ10plot +
  plot_layout(guides = "collect")
Q10Plots


tempPlots  
Q10Plots 



