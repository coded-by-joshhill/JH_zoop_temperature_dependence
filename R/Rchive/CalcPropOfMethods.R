
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


# What proportion of the data had food type and concentration reported? ----
# Quick function to select data across each dataset
food_dat_summary <- function(data, dataset_name) {
  data %>%
    select(ref_no, food_type, food_conc, method) %>%
    mutate(Dataset = dataset_name)
}


# Bind all the data into a tibble
foodDat_binded <- bind_rows(
  food_dat_summary(respiration, "respiration"),
  food_dat_summary(growth, "growth"),
  food_dat_summary(clearance, "clearance"),
  food_dat_summary(ingestion, "ingestion"))


# Summarise the data and calculate the proportion of data that has food info reported
foodDat_binded %>%
  group_by(Dataset) %>%
  summarise(
    prop_food_type_reported = mean(!is.na(food_type)),
    prop_food_conc_reported = mean(!is.na(food_conc)))


# Summarise the experiment type data and calculate the proportions reported.
foodDat_binded %>% 
  mutate(method = if_else(is.na(method), "Not reported", method)) %>%
  group_by(Dataset, method) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Dataset) %>%
  mutate(prop_method = n / sum(n) * 100)


