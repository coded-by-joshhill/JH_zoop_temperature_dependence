# Sensitivity test of growth Q10
# Josh Hill
# 5/3/2026



# Packages and helpers ----
library(tidyverse)
library(glmmTMB)
library(DHARMa)
theme_set(new = theme_bw())




# Growth data
grwdat <- readRDS("Data/grwth_dat.rds") %>% 
  drop_na(Cspecific_rate) %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, rate_value_clean, rate_unit_clean, temp_C, BMC_mg) %>%
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
  ggplot() +
  geom_point(aes(x = temp_C, y = ln_Cspecific_rate, , colour = primRef)) +
  facet_wrap(~ rate_name, scale = "free")

mdat <- grwdat %>% 
  filter(zoopGrp == "Copepods") 


# Histogram of data
mdat %>% 
  ggplot(aes(x = ln_Cspecific_rate)) +
  geom_histogram()


# Lets fit a model with just copepods then with allZ
m1 <- glmmTMB(ln_Cspecific_rate ~ temp_C +  
              (1 | primRef) + (1 | taxa),
              data = mdat)

# Check diagnostics
sim <- simulateResiduals(m1)
plot(sim)
summary(m1)




