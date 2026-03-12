datFinal

# Filter for clearance and select columns
datFinal2 <- datFinal %>% 
  filter(rate_name == "ClearanceRate") %>% 
  select(ref_no, primRef, taxa, temp_C, Cspecific_rate, BM_C, food_type, food_conc, food_conc_unit) %>% 
  drop_na(food_conc) # drop any values that have no food concentration values

# Tidy briefly
datFinal3 <- datFinal2 %>% 
  filter(Cspecific_rate < 2000) %>% 
  group_by(primRef, taxa, temp_C, BM_C) %>% 
  slice_min(Cspecific_rate, n = 1, with_ties = FALSE) %>% 
  ungroup()


datFinal3 %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = log(Cspecific_rate)))


datFinal4 <- datFinal2 %>% 
  filter(Cspecific_rate < 2000) %>% 
  group_by(primRef, taxa, temp_C, BM_C) %>% 
  slice_min(food_conc, n = 1, with_ties = FALSE) %>% 
  ungroup()

datFinal4 %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = log(Cspecific_rate)))


library(glmmTMB)

minRateMod <- glmmTMB(log(Cspecific_rate) ~ temp_C, data = datFinal3)
summary(minRateMod)

minFoodConcMod <- glmmTMB(log(Cspecific_rate) ~ temp_C, data = datFinal4)
summary(minFoodConcMod)

library(emmeans)

minRateSlopes <- emtrends(minRateMod, var = "temp_C")
minRateSlopes


minFoodSlopes <- emtrends(minFoodConcMod, var = "temp_C")
minFoodSlopes

# Check Q10s
minRateslopes_Q10 <- as.data.frame(minRateSlopes) %>% 
  mutate(Q10 = exp(10* temp_C.trend))
minRateslopes_Q10

minFoodlopes_Q10 <- as.data.frame(minFoodSlopes) %>% 
  mutate(Q10 = exp(10* temp_C.trend))
minFoodlopes_Q10


                            