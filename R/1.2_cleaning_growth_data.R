# Cleaning growth data
# Josh Hill
# 03/02/26



  # Here I read in the data
  # Use worrms package to get AphiaID's and taxon classifications
  # Harmonise mass data to mg
  # Assign unique groupings
  # Convert absolute rates to mass-specific
  # Save as an RDS file



# Packages and helpers ----
library(tidyverse)
library(janitor)
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
         rate_value = rate_value / 24, # convert growth rate to hour
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
         rate_value = rate_value / 24, # convert growth rate to hour
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
         rate_value = rate_value  / 24, # convert growth rate to hour
         rate_unit = "mgC/mgC/hr" # update the rate unit
  ) %>% 
  select(-c(C_M1_ug, C_M0_ug, C_M1_mg, C_M0_mg, time_days))
glimpse(ReyDat)


# Read in the main data ----
dat <- read_csv("https://www.dropbox.com/scl/fi/gdllcg9d1dx1dzf38pckd/Grwth_dat.csv?rlkey=xqayol7mkakxn2fdek5yvkkn8&st=az686ihp&dl=1") %>%
  bind_rows(GoldsteinDat, LuskowDat, ReyDat) %>% # add in our estimated growth data
  mutate(ref_no = if_else(is.na(ref_no) | ref_no == "", # if the value is NA or empty...
                          paste0("Hill_", row_number()), # apply a unique reference number
                          ref_no), # otherwise keep what is there
         taxa = str_squish(taxa)) %>% # remove extra spaces from taxon names
  relocate(ref_no, .before = everything()) %>%   # move it before all columns
  filter_out(data_type == "Mean") %>%  # exclude any mean data
  mutate(taxa = case_when(taxa == "Acartia (Acartiura) clausi" ~ "Acartia (Acartiura) clausii", # fix name...
                          TRUE ~ taxa)) %>% # leave the rest as is...
  select(-class)

glimpse(dat)

  
  # Count number of initial pre-cleaned records
  dat %>% group_by(rate_name) %>% 
    summarise(count = n())
  # GrowthRate  488


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
      # 54 records

  

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
    mutate(BMC_mg = convert_CW(BM_C, weight_unit)) %>% # harmonise C weight data to mg
    ungroup() %>% 
    relocate(c(phylum, class, order, family, genus, species), .before = taxa) %>% 
    relocate(BMC_mg, .after = BM_C) %>% 
    mutate(
      # Create custom size groupings following Grigoratou et al. 2025 Figure 1
      sizeGrp = case_when(
        # Mesoplankton: 0.2 mm - 20 mm
        phylum == "Chaetognatha"  ~ "Mesoplankton",
        class == "Appendicularia" ~ "Mesoplankton", # grouped here because we only have Oikopleura dioica
        class == "Copepoda"       ~ "Mesoplankton",
        order == "Pteropoda"      ~ "Mesoplankton", # grouped here because they are Pteropods
        # Macroplankton: 20 mm - 200 mm
        phylum == "Annelida"      ~ "Macroplankton", # grouped here because Tomopteris carpenteri is a larger sp.
        phylum == "Cnidaria"      ~ "Macroplankton",
        phylum == "Ctenophora"    ~ "Macroplankton",
        class == "Hydrozoa"       ~ "Macroplankton",
        class == "Malacostraca"   ~ "Macroplankton",
        class == "Thaliacea"      ~ "Macroplankton",
        # Others - not classified and will be excluded from analyses
        order == "Oegopsida"      ~ "OTHER", # an order of Cephalopod (squid), excluded because likely too large and rare as zoops
        .default = "OTHER"),
      
      # Create custom functional groups based on feeding modes
      funcGrp = case_when(
        # Crustaceans
        class == "Copepoda"       ~ "Crustaceans",
        class == "Malacostraca"   ~ "Crustaceans",
        # Gelatinous filter-feeders
        class == "Appendicularia" ~ "GelFilter",
        class == "Thaliacea"      ~ "GelFilter",
        # Gelatinous predators
        phylum == "Chaetognatha"  ~ "GelPreds", # grouped here because they can be quite gelatinous and are highly predatory
        phylum == "Cnidaria"      ~ "GelPreds",
        phylum == "Ctenophora"    ~ "GelPreds",
        # Others - not classified and will be excluded from analyses
        phylum == "Annelida"      ~ "OTHER", # excluded because Tomopteris are quite gelatinious and are generally pretty rare
        phylum == "Mollusca"      ~ "OTHER", # excluded because most are pteropods and don't fit into feeding classification
        .default = "OTHER"),
      
      # Create custom groupings for general zoop groups following Ikeda 2014
      zoopGrp = case_when( 
        phylum == "Annelida"      ~ "Annelids",
        phylum == "Chaetognatha"  ~ "Chaetognaths",
        phylum == "Cnidaria"      ~ "Cnidarians",
        phylum == "Ctenophora"    ~ "Ctenophores",
        phylum == "Mollusca"      ~ "Molluscs",
        phylum == "Rotifera"      ~ "Rotifers",
        class == "Appendicularia" ~ "Appendicularians",
        class == "Copepoda"       ~ "Copepods",
        class == "Thaliacea"      ~ "Thaliaceans",
        order == "Euphausiacea"   ~ "Euphausiids",
        order == "Amphipoda"      ~ "Amphipods",
        order == "Decapoda"       ~ "Decapods",
        order == "Mysidacea"      ~ "Mysids",
        order == "Mysida"         ~ "Mysids",
        .default = "OTHER"),
    ) %>% 
    relocate(zoopGrp, .before = phylum) %>% 
    relocate(funcGrp, .before = zoopGrp) %>% 
    relocate(sizeGrp, .before = funcGrp)

glimpse(datClean)  
# End data cleaning ----

     
    
# Prep data for analysis ----
datFinal <- datClean %>%
  # Estimate dry mass based on the % from Kiorboe's 2013 Table 1 
  mutate(BM_dry_mg = case_when(zoopGrp == "Appendicularians" ~ BMC_mg / (10.3/100),
                               zoopGrp == "Chaetognaths" ~ BMC_mg / (36.7/100),
                               zoopGrp == "Cnidarians"   ~ BMC_mg / (13.2/100),
                               zoopGrp == "Copepods"     ~ BMC_mg / (48/100),
                               zoopGrp == "Ctenophores"  ~ BMC_mg / (5.1/100),
                               zoopGrp == "Euphausiids"  ~ BMC_mg / (41.9/100),
                               zoopGrp == "Thaliaceans"  ~ BMC_mg / (10.3/100))) %>% 
  # No need to harmonise this data (unless more is added)... so I'll use a simple mutate to maintain naming convention consistency
  # Carbon mass-specific rates (i.e., relative growth rates...)
  mutate(Cspecific_rate = rate_value, # rate_value is already mass-specific so merge to consistent naming convention
         Cspecific_unit = "mgC/mgC/hr", # update the units

         # Convert to dry-mass specific rates
         DrySpecific_rate = case_when(
           rate_name == "GrowthRate" & Cspecific_unit == "mgC/mgC/hr" ~ Cspecific_rate, # value stays the same because it is relative growth rate
           TRUE ~ NA_real_),
         DrySpecific_unit = case_when( # update the mass-specific units to match the mass-specific rates...
           rate_name == "GrowthRate" & !is.na(DrySpecific_rate) ~ "mgDry/mgDry/hr",
           TRUE ~ Cspecific_unit),
         sizeGrp = as.factor(sizeGrp),
         funcGrp = as.factor(funcGrp),  
         zoopGrp = as.factor(zoopGrp)) %>%
  select(-c(BM_wet, BM_dry, BM_C, weight_unit, weight_calc)) %>%  # tidy up the dataframe
  relocate(Cspecific_rate, .after = rate_name) %>% 
  relocate(Cspecific_unit, .after = Cspecific_rate)
  
glimpse(datFinal)

# End conversion

# Final checks
  # sizeGrp
  datFinal %>% 
    ggplot() + 
    geom_point(aes(x = temp_C, 
                   y = log(Cspecific_rate), # mass-specific
                   colour = primRef)) #+
    #facet_wrap(~ zoopGrp, scales = "free")
  
  # sizeGrp
  datFinal %>% 
    ggplot() + 
    geom_point(aes(x = temp_C, 
                   y = log(Cspecific_rate), # absolute
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
  
  # Count unique ZoopGrps rates
  datClean %>% 
    group_by(zoopGrp) %>% 
    mutate(countZGrp = sum(zoopGrp > 1, na.rm = TRUE)) %>% 
    distinct(zoopGrp, countZGrp) %>% 
    arrange(countZGrp)
    # 1 Mysids              10
    # 2 Chaetognaths        14
    # 3 Decapods            27
    # 4 Euphausiids         33
    # 5 Thaliaceans         34
    # 6 Amphipods           63
    # 7 Cnidarians          65
    # 8 Ctenophores        102
    # 9 Copepods           140
  
  
  # Count unique functional groups rates
  datClean %>% 
    group_by(funcGrp) %>% 
    mutate(countFuncGrp = sum(funcGrp > 1, na.rm = TRUE)) %>% 
    distinct(funcGrp, countFuncGrp) %>% 
    arrange(countFuncGrp)
    # 1 GelFilter             34
    # 2 GelPreds             181
    # 3 Crustaceans          273
  
  
  # Count unique size groups rates
  datClean %>% 
    group_by(sizeGrp) %>% 
    mutate(countSizeGrp = sum(sizeGrp > 1, na.rm = TRUE)) %>% 
    distinct(sizeGrp, countSizeGrp) %>% 
    arrange(countSizeGrp)
    # 1 Mesoplankton           154
    # 2 Macroplankton          334
  
  
  # Check the temperature range for each zoopGrp
  datClean %>% 
    select(zoopGrp, temp_C, rate_name) %>% 
    group_by(zoopGrp) %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
    # All have sensible ranges for estimating Q10, except for Mysids
  
# Save it as an RDS for later use
# saveRDS(datFinal, "Data/grwth_dat.rds")

  