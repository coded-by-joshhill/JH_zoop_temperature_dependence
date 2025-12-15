# Temperature sensitivity ~ carbon body mass
# Josh Hill
# 15/12/25



  # Read in Q10 data and cleaned rate and weight data
  # Estimate the mean BMC_mg weight for each taxa recorded
  # Bind the Q10 data to the weight data based on zoopGrp. 



# Packages and helpers ----
library(tidyverse)
library(glmmTMB)
source("R/0_Helpers.R")



# Read in the data ----
clearance_Q10pdat <- readRDS("Data/Q10_summary_clearance.rds") %>% 
  rename(zoopGrp = Group)

dat <- readRDS("Data/clear_ingest_data.rds") %>% 
  filter(rate_name == "ClearanceRate", !BMC_mg == "NA") %>% 
  mutate(zoopGrp = as.factor(zoopGrp)) %>%  # set to factor
  select(zoopGrp, taxa, BMC_mg)

# Estimate mean BMC
meanMass <- dat %>% 
  group_by(zoopGrp) %>% 
  summarise(log_mean_BMC_mg = log(mean(BMC_mg))) #%>%
  # drop_na(mean_BMC_mg)

datMean <- meanMass %>% 
  left_join(clearance_Q10pdat, by = "zoopGrp") %>% 
  rename(median_Q10 = median) %>% 
  drop_na()

datMean %>% 
  ggplot()+
  geom_point(aes(x = log_mean_BMC_mg, y = median_Q10))
  


