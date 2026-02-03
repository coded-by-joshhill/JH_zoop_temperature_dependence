# Cleaning growth data
# Josh Hill
# 03/02/26



  # Here I read in the data to estimate growth, read in the compiled Pata and Hunt data,
  # Use worrms package to get AphiaID's and taxon classifications
  # Harmonise weight data to mg
  # Convert absolute rates to mass-specific
  # Save as an RDS file



# Packages and helpers ----
library(tidyverse)
library(worrms)
source("R/0_Helpers.R")



# Read in the data ----
# First, estimate growth rates

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

  
  # Count number of initial pre-cleaned records
  dat %>% group_by(rate_name) %>% 
    summarise(count = n())
  # GrowthRate  686


  # Look at all unique taxon
  dat %>% 
    filter(rate_name == "GrowthRate") %>% 
    distinct(taxa) %>% 
    arrange(taxa) %>% 
    print(n = "Inf")
 
  
  # Look at all unique primary references for each rate type
  dat %>% 
    group_by(rate_name) %>% 
    distinct(primRef, rate_name) %>% 
    summarise(count = n())
      # 55 records

  

# Subset taxa data to get AphiaIDs and classifications ----
taxaDat <- dat %>% 
  select(taxa) %>% 
  distinct(taxa) %>% 
  arrange(taxa)
  
  
  # Get AphiaIDs
  taxaDatID <- taxaDat %>% 
    mutate(AphiaID = map_int(taxa, wm_name2id))  # extract AphiaIDs for taxa using the worrms package

    taxaDatID %>% 
      summary() # ensure there are no NAs


  # Add classifications using AphiaIDs
  taxaDatClass <- taxaDatID %>%
    mutate(classification = map(AphiaID, ~ { # get taxonomic classifications for each AphiaID 
      wm_classification(.x)}))
  

  # Unnest classification info and update taxaDat
  taxaDat <- taxaDatClass %>%
    select(taxa, classification) %>%
    unnest(classification) %>% # Unnest classification information into separate columns
    select(taxa, rank, scientificname) %>%
    pivot_wider(names_from = rank,
                values_from = scientificname) %>% 
    select(taxa, Phylum, Class, Order, Family, Genus, Species) %>% # only keep necessary columns
    clean_names(case = "snake")

  

# Join taxa info and harmonise weight data ----
datClean <- dat %>% 
    left_join(taxaDat, by = "taxa") %>% 
    rowwise() %>% 
    mutate(BMC_mg = convert_CW(BM_C, weight_unit)) %>% # harmonize C weight data to mg
    ungroup() %>% 
    relocate(c(phylum, class, order, family, genus, species), .before = taxa) %>% 
    relocate(BMC_mg, .after = BM_C) %>% 
    mutate(
      # Create custom size groupings following Grigoratou et al. 2025 Figure 1
      sizeGrp = case_when(
        # Mesoplankton: 0.2 um - 20 mm
        class == "Copepoda"       ~ "Mesoplankton",
        class == "Appendicularia"  ~ "Mesoplankton", # grouped here because we only have Oikopleura dioica
        # Macroplankton: 20 mm - 200 mm
        class == "Malacostraca"   ~ "Macroplankton",
        class == "Hydrozoa"       ~ "Macroplankton",
        class == "Scyphozoa"      ~ "Macroplankton",
        class == "Tentaculata"    ~ "Macroplankton",
        class == "Thaliacea"      ~ "Macroplankton",
        class == "Sagittoidea"    ~ "Macroplankton",
        .default = "OTHER"),
      
      # Create custom functional groups based on feeding modes
      funcGrp = case_when(
        # Crustaceans
        class == "Malacostraca"   ~ "Crustaceans",
        class == "Copepoda"       ~ "Crustaceans",
        # Gelatinous filter-feeders
        class == "Thaliacea"      ~ "GelFilter",
        class == "Appendicularia" ~ "GelFilter",
        # Gelatinous predators
        phylum == "Cnidaria"      ~ "GelPreds",
        phylum == "Ctenophora"    ~ "GelPreds",
        phylum == "Chaetognatha"  ~ "GelPreds", # grouped with GelPreds because  although not full-gelatinous, this makes more sense than crustaceans
        .default = "OTHER"),
      
      # Create custom groupings for general zoop groups following Ikeda 2014
      zoopGrp = case_when( 
        order == "Euphausiacea"   ~ "Euphausiacea",
        order == "Amphipoda"      ~ "Amphipoda",
        order == "Decapoda"       ~ "Decapoda",
        order == "Mysida"         ~ "Mysida",
        class == "Copepoda"       ~ "Copepoda",
        phylum == "Mollusca"      ~ "Mollusca",
        phylum == "Chaetognatha"  ~ "Chaetognatha",
        phylum == "Cnidaria"      ~ "Cnidaria",
        phylum == "Ctenophora"    ~ "Ctenophora",
        class == "Thaliacea"      ~ "Thaliacea",
        class == "Appendicularia" ~ "Appendicularia",
        .default = "OTHER"),
    ) %>% 
    relocate(zoopGrp, .before = phylum) %>% 
    relocate(sizeGrp, .before = zoopGrp) %>% 
    relocate(funcGrp, .before = sizeGrp)
    
    
  # Count unique ZoopGrps rates
  datClean %>% 
    group_by(zoopGrp) %>% 
    mutate(countZGrp = sum(zoopGrp > 1, na.rm = TRUE)) %>% 
    distinct(zoopGrp, countZGrp) %>% 
    arrange(countZGrp)
      # Mysida                14
      # Decapoda              28
      # Chaetognatha          28
      # Euphausiacea          34
      # Thaliacea             34
      # Amphipoda             63
      # Cnidaria              68
      # Ctenophora           132
      # Copepoda             285
  
  
  # Check the temperature range for each zoopGrp
  datClean %>% 
    select(zoopGrp, temp_C, rate_name) %>% 
    group_by(zoopGrp) %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
      # All have sensible ranges for estimating Q10, except for Mysids
  
# End data cleaning ----

     
    
# Prep data for analysis ----
datFinal <- datClean %>%
  # No need to harmoise this data (unless more is added)... so I'll use a simple mutate to maintain naming convention consistency
  mutate(Cspecific_rate = rate_value,
         final_unit = rate_unit,
         zoopGrp = as.factor(zoopGrp)) %>%
  relocate(Cspecific_rate, .after = rate_name) %>% 
  relocate(final_unit, .after = Cspecific_rate) %>% 
  filter(ref_no != "Pata_excl_1460") # remove this huge outlier by Kasuya2002
glimpse(datFinal)

# End conversion

# Final checks
# sizeGrp
datFinal %>% 
  ggplot() + 
  geom_point(aes(x = temp_C, 
                 y = log(Cspecific_rate), 
                 colour = sizeGrp))


# funcGrp
datFinal %>% 
  ggplot() + 
  geom_point(aes(x = temp_C, 
                 y = log(Cspecific_rate), 
                 colour = funcGrp))


# zoopGrp
datFinal %>% 
  ggplot() + 
  geom_point(aes(x = temp_C, 
                 y = log(Cspecific_rate), 
                 colour = zoopGrp))
  
# Save it as an RDS for later use
# saveRDS(datFinal, "Data/grwth_dat.rds")
