# Supplementary table of primary authors
# Josh H
# 27/07/2026



# Packages ----
library(tidyverse)



# Read in the data ----
# Feeding data
dat <- readRDS("Data/clear_ingest_data.rds")


# Clearance data
cleardat <- dat %>% 
  filter(rate_name == "ClearanceRate") %>% # Filter for clearance rate
  select(primRef, primRef_URL) %>% 
  mutate(rate_name = factor("Clearance"))

# Ingestion data
ingdat <- dat %>% 
  filter(rate_name == "IngestionRate") %>% # Filter for ingestion rate
  select(primRef, primRef_URL) %>% 
  mutate(rate_name = factor("Ingestion"))

# Growth data
grwdat <- readRDS("Data/grwth_dat.rds") %>% 
  select(primRef, primRef_URL) %>% 
  mutate(rate_name = factor("Growth"))

# Respiration data
respdat <- readRDS("Data/resp_dat.rds") %>% 
  select(primRef, primRef_URL) %>% 
  mutate(rate_name = factor("Respiration"))

# Excretion data
excredat <- readRDS("Data/excrete_dat.rds") %>% 
  select(primRef, primRef_URL) %>% 
  mutate(rate_name = factor("Excretion"))


# Combine them into one...
usedat <- rbind(cleardat, ingdat, grwdat, respdat, excredat)


# Give me a table of unique primary authors

usedat %>% 
  distinct(primRef, primRef_URL) %>% 
  arrange(primRef) %>%
  drop_na(primRef) %>% 
  print(n = Inf)


