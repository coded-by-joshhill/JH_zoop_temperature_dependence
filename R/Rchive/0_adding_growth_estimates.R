# Cleaning growth data
# Josh Hill
# 14/11/25



# Here I read in the data
# Use worrms package to get AphiaID's and taxon classifications
# Harmonise weight data to mg
# Convert absolute rates to mass-specific
# Save as an RDS file



# Packages and helpers ----
library(tidyverse)
library(worrms)
library(janitor)
source("R/0_Helpers.R")



# Estimate growth rates first ---

# Goldstein2023
GoldsteinDat <- read_csv("https://www.dropbox.com/scl/fi/6o7ovfttvc7roz0kc27n1/Goldstein2023_growth_estimates.csv?rlkey=lfgibkpio5b21kw1gvu52avt7&st=ee82n8vz&dl=1",
                         skip = 1) %>% 
  mutate(BM_C = C_M0_mg, # carbon mass gets the initial weight measurement
         weight_unit = "mg", # update the weight unit
         rate_value = calcGrowthRate(C_M1_mg, C_M0_mg, time_days), # estimate growth rate with function
         rate_value = (rate_value /BM_C) / 24, # convert growth rate to mass specific and to hour
         rate_unit = "mgC/mgC/hr") %>% # update the rate unit
  select(-c(C_M1_mg, C_M0_mg, time_days)) # tidy the dataframe and prep for binding
glimpse(GoldsteinDat)


# Luskow2016
LuskowDat <- read_csv("https://www.dropbox.com/scl/fi/uyqkquws5jvw7gfbd95oa/Luskow2016_growth_estimates.csv?rlkey=4jsx0e9c1hqzwaebei1o0dyj7&st=4ash4wxl&dl=1",
                      skip = 1) %>% 
  mutate(C_M1_mg = DW_M1_mg * (13.2/100), # estimate carbon mass using Kiorboe 2013 DW concentration (13.2%)
         C_M0_mg = DW_M0_mg * (13.2/100),
         BM_dry = DW_M0_mg, # update dataset dry weight (mg)
         BM_C = C_M0_mg, # carbon mass gets the initial weight measurement (in carbon form)
         weight_unit = "mg", # update the weight unit
         rate_value = calcGrowthRate(C_M1_mg, C_M0_mg, time_days), # estimate growth rate with function
         rate_value = (rate_value /BM_C) / 24, # convert growth rate to mass specific and to hour
         rate_unit = "mgC/mgC/hr" # update the rate unit
  ) %>% 
  select(-c(DW_M1_mg, DW_M0_mg, C_M1_mg, C_M0_mg, time_days))
glimpse(LuskowDat)


# Rey2001
ReyDat <- read_csv("https://www.dropbox.com/scl/fi/896t9633il9zjkmev423p/Rey2001_growth_estimates.csv?rlkey=0zahgx268iu0zb38hzfvi7zrp&st=k0zencye&dl=1",
                   skip = 1) %>% 
  mutate(C_M1_mg = C_M1_ug / 1000, # Convert ugC to mgC
         C_M0_mg = C_M0_ug / 1000,
         BM_C = C_M0_mg, # carbon mass gets the initial weight measurement
         weight_unit = "mg", # update the weight unit
         rate_value = calcGrowthRate(C_M1_mg, C_M0_mg, time_days), # estimate growth rate with function
         rate_value = (rate_value /BM_C) / 24, # convert growth rate to mass specific and to hour
         rate_unit = "mgC/mgC/hr" # update the rate unit
  ) %>% 
  select(-c(C_M1_ug, C_M0_ug, C_M1_mg, C_M0_mg, time_days))
glimpse(ReyDat)


# Read in the main data ----
dat <- read_csv("https://www.dropbox.com/scl/fi/gdllcg9d1dx1dzf38pckd/Grwth_dat.csv?rlkey=xqayol7mkakxn2fdek5yvkkn8&st=49ooa3to&dl=1",
                skip = 1) %>%
  bind_rows(GoldsteinDat, LuskowDat, ReyDat) %>% # add in our estimated growth data
  mutate(ref_no = if_else(is.na(ref_no) | ref_no == "", # if the value is NA or empty...
                          paste0("Hill_", row_number()), # apply a unique reference number
                          ref_no), # otherwise keep what is there
         taxa = str_squish(taxa)) %>% # remove extra spaces from taxon names
  relocate(ref_no, .before = everything())  # move it before all columns
glimpse(dat)






