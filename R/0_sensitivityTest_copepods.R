# Sensitivity test of growth Q10
# Josh Hill
# 5/3/2026



# Packages and helpers ----
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(emmeans) # Estimated marginal means

theme_set(new = theme_bw())




# Growth data
grwdat <- readRDS("Data/grwth_dat.rds") %>% 
  drop_na(Cspecific_rate) %>% 
  select(ref_no, life_stage, primRef, secRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, rate_value_clean, rate_unit_clean, temp_C, BMC_mg) %>%
  mutate(rate_name = factor("Growth"),
         ln_Cspecific_rate = log(Cspecific_rate))


grwdat %>% distinct(zoopGrp)
grwdat %>% distinct(primRef) %>% arrange(primRef) %>% print(n = "Inf")
  # The with growth data worth look# The paper with growth data worth looking at is: Hirst and Bunker which has a heap of data from different authors

grwdat %>% distinct(zoopGrp)



# Quick look
grwdat %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = sizeGrp)) +
  facet_wrap(~ rate_name, scale = "free")

grwdat %>% 
  filter(zoopGrp == "Copepods") %>%
  filter(primRef == "Kiorboe2014") %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = primRef)) +
  facet_wrap(~ rate_name, scale = "free")

mdat <- grwdat %>% 
  filter(zoopGrp == "Copepods") %>% 
  filter(primRef == "Kiorboe2014")
  

# Histogram of data
mdat %>% 
  ggplot(aes(x = ln_Cspecific_rate)) +
  geom_histogram()


# Lets fit a model with just copepods then with allZ
m1 <- glmmTMB(ln_Cspecific_rate ~ temp_C + (1 | taxa),
              data = mdat)

# Check diagnostics
sim <- simulateResiduals(m1)
plot(sim)
summary(m1)

# Set min and max temp
minTempC <- min(mdat$temp_C)
maxTempC <- max(mdat$temp_C)


# FUNCTIOIN TO PLOT MODEL ----
PlotLMM = function(model){
  temp_seq <- seq(minTempC, maxTempC, length.out = 100)
  
  # Build newdata grid manually
  newdat <- expand.grid(
    temp_C    = temp_seq
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
                aes(x = temp_C, ymin = conf.low, ymax = conf.high),
                alpha = 0.25) +
    geom_line(data = pop_preds,
              aes(x = temp_C, y = estimate),
              linewidth = 1) +
    geom_point(data = mdat,
               aes(x = temp_C, y = ln_Cspecific_rate),
               alpha = 0.2) +
    coord_cartesian(xlim = c(-2, 32)) +
    labs(x = "Temp (°C)",
         y = "ln(Carbon-mass specific growth rate)",
         title = "Copepods mass-specific growth rate")
 }

PlotLMM(m1)

# Extract slopes using and calculate Q10 for each sizeGrp
slopes <- emtrends(m1, var = "temp_C")

slopes_Q10 <- as.data.frame(slopes) %>% 
  mutate(Q10 = exp(10 * temp_C.trend))
slopes_Q10

