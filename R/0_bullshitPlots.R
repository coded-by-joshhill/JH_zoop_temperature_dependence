library(tidyverse)

# Read in the data and prep for analysis ----

# Feeding data
dat <- readRDS("Data/clear_ingest_data.rds")


# Clearance data
cleardat <- dat %>% 
  filter(rate_name == "ClearanceRate") %>% # Filter for clearance rate
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, rate_value_clean, rate_unit_clean, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Clearance"),
         BMC_15 = BMC_mg * 2.35^((15-temp_C) / 10))


cleardat %>% ggplot() +
  geom_point(aes(x = log10(BMC_15), y = log10(rate_value_clean), colour = primRef,
                 shape = zoopGrp),
             size = 3)

cleardat %>% ggplot() +
  geom_point(aes(x = log10(BMC_15), y = log10(Cspecific_rate), colour = zoopGrp),
             size = 3)

cleardat %>% ggplot() +
  geom_point(aes(x = temp_C, y = log10(Cspecific_rate), colour = zoopGrp),
             size = 3)

# clearance looks good!

# Ingestion data
ingdat <- dat %>% 
  filter(rate_name == "IngestionRate") %>% # Filter for ingestion rate
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, rate_value_clean, rate_unit_clean, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Ingestion"),
         BMC_15 = BMC_mg * 2.35^((15-temp_C) / 10))


ingdat %>% ggplot() +
  geom_point(aes(x = log10(BMC_15), y = log10(rate_value_clean), colour = zoopGrp), size = 3) +
  facet_wrap(~ rate_unit_clean)

ingdat %>% ggplot() +
  geom_point(aes(x = log10(BMC_15), y = log10(Cspecific_rate), colour = zoopGrp), size = 3)

ingdat %>% ggplot() +
  geom_point(aes(x = temp_C, y = log10(Cspecific_rate), colour = zoopGrp))

# Ingestion looks good too...


# Growth data
grwdat <- readRDS("Data/grwth_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Growth"),
         BMC_15 = BMC_mg * 2.35^((15-temp_C) / 10))

grwdat %>% ggplot() +
  geom_point(aes(x = log10(BMC_15), y = log10(Cspecific_rate * BMC_15), colour = zoopGrp), size = 3)

grwdat %>% ggplot() +
  geom_point(aes(x = log10(BMC_15), y = log10(Cspecific_rate), colour = zoopGrp), size = 3)

grwdat %>% ggplot() +
  geom_point(aes(x = temp_C, y = log10(Cspecific_rate), colour = zoopGrp))


# Growth looks right....


# Respiration data
respdat <- readRDS("Data/resp_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, rate_value_clean, rate_unit_clean, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  # mutate(BMC_mg = case_when(primRef == "Ross1982" ~ 10^(BMC_mg),
  #                          TRUE ~ BMC_mg)) %>%
  mutate(rate_name = factor("Respiration"),
         BMC_15 = BMC_mg * 2.35^((15-temp_C) / 10)) 

respdat %>% distinct(rate_unit_clean)

respdat %>%  
  filter(rate_unit_clean == "ulO2/ind/hr") %>% 
  # filter(primRef == "Ross1982") %>% 
  ggplot() +
  geom_point(aes(x = log(BMC_15), y = log(rate_value_clean), colour = zoopGrp))

respdat %>% ggplot() +
  geom_point(aes(x = log10(BMC_15), y = log10(Cspecific_rate), colour = primRef))
# Clearly an issue with ross data... sad but will just have to get rid of it for now

respdat %>% ggplot() +
  geom_point(aes(x = temp_C, y = log10(Cspecific_rate), colour = primRef))


# Excretion data
excredat <- readRDS("Data/excrete_dat.rds") %>% 
  select(primRef, sizeGrp, funcGrp, zoopGrp, taxa, rate_value, rate_unit, Cspecific_rate, Cspecific_unit, temp_C, BMC_mg) %>% 
  drop_na(Cspecific_rate) %>% 
  mutate(rate_name = factor("Excretion"))


# Combine them into one...
