
# Packages and helpers ----
library(tidyverse)
library(gt)
library(patchwork)
source("R/0_Helpers.R")



# Read in and filter the cleaned data ----
clearance <- readRDS("Data/clear_ingest_data.rds") %>% 
  filter(rate_name == "ClearanceRate") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)


ingestion <- readRDS("Data/clear_ingest_data.rds") %>% 
  filter(rate_name == "IngestionRate") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)


growth <- readRDS("Data/grwth_dat.rds") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)


respiration <- readRDS("Data/resp_dat.rds") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)



glimpse(respiration)
  


# BS figures looking at effect of carbon mass on Cspecific rate / temp
clearance %>% 
  select(Cspecific_rate, temp_C, BMC_mg, zoopGrp) %>% 
  filter(Cspecific_rate < 2500, zoopGrp == "Copepoda") %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = Cspecific_rate, colour = log(BMC_mg)), size = 3) +
  scale_color_viridis_c() +
  labs(title = "Copepod respiration",
       x = "Temperature °C",
       y = "Mass-specific respiration") +
  theme_classic()
  
  

# Fit a simple model
m1 <- glmmTMB(Cspecific_rate ~ temp_C * zoopGrp + (1|primRef) + (1|taxa),  data = respiration)

summary(m1)
performance::r2(m1)
# Model is performing well
# Can see there is extremely large variance in the model due to the random effect of taxon... 
# primaryRef seems to have very little effect on variance.

growth%>% 
  select(Cspecific_rate, temp_C, BMC_mg) %>% 
  filter(Cspecific_rate < 0.075) %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = Cspecific_rate, colour = log(BMC_mg)))


ingestion%>% 
  select(Cspecific_rate, temp_C, BMC_mg) %>% 
  filter(Cspecific_rate < 0.2) %>% 
  ggplot() +
  geom_point(aes(x = temp_C, y = Cspecific_rate, colour = log(BMC_mg)))

clearance%>% 
  select(Cspecific_rate, temp_C, BMC_mg) %>% 
  filter(Cspecific_rate < 10000) %>%
  ggplot() +
  geom_point(aes(x = temp_C, y = Cspecific_rate, colour = log(BMC_mg)))


#I also think I found a way to explain how we still have positive trends despite -1/4 scaling of mass to mass-specific metabolic rate.... it seems temperature has a more powerful effect on mass-specific metabolic rate compared to mass alone. We see greater carbon mass clustered clustered lowere than smaller masses despite the associated mass-specific metabolic rate being higher. This generally increases as temperature rises. 


