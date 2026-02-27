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
library(glmmTMB)
library(DHARMa) # for diagnostics
library(MuMIn) # for Rsqr
library(marginaleffects)
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



# Fit the models ----


# A model with 3 way interactions for temp, sizeGrp and rate with random effects
m1 <- glmmTMB(ln_Cspecific_rate ~ temp_C * rate_name * sizeGrp + 
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts and slopes
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m1)
  plot(sim) # looks pretty schmick
  summary(m1)
  r.squaredGLMM(m1)
  

# A simpler model without the 3-way interaction
m2 <- glmmTMB(ln_Cspecific_rate ~ temp_C * rate_name + temp_C * sizeGrp + rate_name * sizeGrp +
                (temp_C | primRef) + (temp_C | taxa), # with primRef and taxa as random intercepts
              data = mdat) 

  # Check diagnostics
  sim <- simulateResiduals(m2)
  plot(sim) # Looks fine but is less scattered compared to m1
  summary(m2)
  r.squaredGLMM(m2)

# Likelihood ratios test of the two models
anova(m1, m2) # anova of complex vs simpler models
# likelihood ratios test shows m1 has significantly more explanatory power
# AIC is also slightly better despite the BIC being slightly higher
# I will use m1 on the basis of the chisqr test and AIC...


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
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}



# Plot it
tempPlot <- PlotLMM(m1)
tempPlot
  
  


# Extract slopes using and calculate Q10 for each sizeGrp
sizeGrp_slopes <- emtrends(m1, ~ rate_name * sizeGrp, var = "temp_C")

sizeGrp_slopes_Q10 <- as.data.frame(sizeGrp_slopes) |>
  mutate(Q10 = exp(10 * temp_C.trend),
         Q10_lwr = exp(10 * asymp.LCL),
         Q10_upr = exp(10 * asymp.UCL))
sizeGrp_slopes_Q10



# Plot sizeGrp Q10s
sizeGQ <- ggplot() +
  geom_errorbar(data = sizeGrp_slopes_Q10, 
                aes(x = sizeGrp, ymin = Q10_lwr, ymax = Q10_upr, colour = sizeGrp),
                width = .05,
                linewidth = 1) +
  geom_point(data = sizeGrp_slopes_Q10, aes(x = sizeGrp, y = Q10),
             size = 3,
             colour = "black") +
  facet_wrap(~rate_name, scale = "free")



# Extract slopes and calculate Q10 for overall zooplankton
slopes <- emtrends(m1, ~ rate_name, var = "temp_C")

slopes_Q10 <- as.data.frame(slopes) |> 
  mutate(Q10 = exp(10* temp_C.trend),
         Q10_lwr = exp(10 * asymp.LCL),
         Q10_upr = exp(10 * asymp.UCL))
slopes_Q10


# Plot allZoop Q10s
Q <- ggplot() +
  geom_errorbar(data = slopes_Q10, 
                aes(x = rate_name, ymin = Q10_lwr, ymax = Q10_upr), 
                colour = "grey",
                width = .05,
                linewidth = 1) +
  geom_point(data = slopes_Q10, aes(x = rate_name, y = Q10),
             size = 2,
             colour = "black")

library(patchwork)

tempPlot/sizeGQ







# PLot working:

# FUNCTION TO PLOT MODEL ----
PlotLMM = function(model){
  temp_seq <- seq(minTempC, maxTempC, length.out = 100)
  
  # Build newdata grid manually
  newdat <- expand.grid(
    temp_C    = temp_seq,
    sizeGrp   = unique(mdat$sizeGrp),
    rate_name = unique(mdat$rate_name)
  )
  
  # Population-level predictions (fixed effects only)
  pop_preds <- newdat
  pred <- predict(model, newdata = newdat, se.fit = TRUE,
                  re.form = NA, allow.new.levels = TRUE)
  pop_preds$estimate <- pred$fit
  pop_preds$conf.low  <- pred$fit - 1.96 * pred$se.fit
  pop_preds$conf.high <- pred$fit + 1.96 * pred$se.fit
  
  # Rate-specific settings
  rate_settings <- list(
    "Clearance"   = list(col = "#66c2a5", ylab = expression(bold("ln(Clearance rate) (ml mgC"^-1~"h"^-1*")"))),
    "Ingestion"   = list(col = "#fc8d62", ylab = expression(bold("ln(Ingestion rate) (mgC mgC"^-1~"h"^-1*")"))),
    "Growth"      = list(col = "#8da0cb", ylab = expression(bold("ln(Growth rate) (mgC mgC"^-1~"h"^-1*")"))),
    "Respiration" = list(col = "#e78ac3", ylab = expression(bold("ln(Respiration rate) (µlO"[2]~"mgC"^-1~"h"^-1*")")))
  )
  
  # Define sizeGrp colours
  sizecols <- c("Mesoplankton" = "black",  # or whatever your levels are called
                "Macroplankton" = "snow3")
  
  # Function to build each panel
  plot_rate <- function(rate) {
    settings <- rate_settings[[rate]]
    
    ggplot() +
      geom_point(data = mdat %>% filter(rate_name == rate),
                 aes(x = temp_C, y = ln_Cspecific_rate, colour = sizeGrp),
                 alpha = 0.2) +
      geom_ribbon(data = pop_preds %>% filter(rate_name == rate),
                  aes(x = temp_C, ymin = conf.low, ymax = conf.high, fill = sizeGrp),
                  alpha = 0.25) +
      geom_line(data = pop_preds %>% filter(rate_name == rate),
                aes(x = temp_C, y = estimate, colour = sizeGrp),
                linewidth = 1) +
      scale_colour_manual(values = sizecols) +
      scale_fill_manual(values = sizecols) +
      coord_cartesian(xlim = c(-2, 32)) +
      labs(x = "Temp (°C)", y = settings$ylab,
           fill = "Size group", colour = "Size group",
           title = rate) +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5,
                                  colour = "white", size = 11),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        title = element_text(size = 10),
        strip.background = element_rect(fill = settings$col, colour = "black"),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  }
  
  # Build all 4 panels
  p1 <- plot_rate("Clearance")
  p2 <- plot_rate("Ingestion")
  p3 <- plot_rate("Growth")
  p4 <- plot_rate("Respiration")
  
  # Combine with patchwork
  (p1 + p2 + p3 + p4) + 
    plot_layout(guides = "collect", ncol = 2) &
    theme(legend.position = "right")
}

# Plot it
tempPlot <- PlotLMM(m1)
tempPlot

