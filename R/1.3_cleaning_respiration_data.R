# Cleaning respiration data
# Josh Hill
# 03/02/26



  # Here I read in the data
  # Use worrms package to get AphiaID's and taxon classifications
  # Harmonise weight data to mg
  # Assign unique groupings
  # Convert absolute rates to mass-specific
  # Save as an RDS file



# Packages and helpers ----
library(tidyverse)
library(janitor)
library(worrms)
source("R/0_Helpers.R")



# Read in the data ----
dat <- read_csv("https://www.dropbox.com/scl/fi/itv9vpnu8twxz2fyhmp1u/Resp_dat.csv?rlkey=9fig4vcw2cog4rc4qa6mtj7lj&st=6rhjjioa&dl=1", 
                skip = 1) %>%
  mutate(ref_no = if_else(is.na(ref_no) | ref_no == "", # if the value is NA or empty...
                          paste0("Hill_", row_number()), # apply a unique reference number
                          ref_no), # otherwise keep what is there
         taxa = str_squish(taxa)) %>% # remove extra spaces from taxon names
  relocate(ref_no, .before = everything())  # move it before all columns
glimpse(dat)

  
  # Count number of initial pre-cleaned records
  dat %>% group_by(rate_name) %>% 
    summarise(count = n())
  # RespirationRate  1036


  # Look at all unique taxon
  dat %>% 
    filter(rate_name == "RespirationRate") %>% 
    distinct(taxa) %>% 
    arrange(taxa) %>% 
    print(n = "Inf")
 
  
  # Look at all unique primary references for each rate type
  dat %>% 
    group_by(rate_name) %>% 
    distinct(primRef, rate_name) %>% 
    summarise(count = n())
      # 23 records

  

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
        class == "Appendicularia"  ~ "Mesoplankton", # grouped here because we only have Oikopleura dioica
        class == "Copepoda"        ~ "Mesoplankton",
        order == "Pteropoda"       ~ "Mesoplankton", # grouped here because they are Pteropods
        # Macroplankton: 20 mm - 200 mm
        phylum == "Annelida"      ~ "Macroplankton", # grouped here because Tomopteris carpenteri is a larger sp.
        phylum == "Chaetognatha"  ~ "Macroplankton",
        phylum == "Cnidaria"      ~ "Macroplankton",
        phylum == "Ctenophora"    ~ "Macroplankton",
        class == "Hydrozoa"       ~ "Macroplankton",
        class == "Malacostraca"   ~ "Macroplankton",
        class == "Thaliacea"      ~ "Macroplankton",
        order == "Oegopsida"      ~ "Macroplankton", # an order of Cephalopod, not explicitly reported as larvae so grouped here
        .default = "OTHER"),
      
      # Create custom functional groups based on feeding modes
      funcGrp = case_when(
        # Crustaceans and others
        phylum == "Annelida"      ~ "CrustOthers",
        phylum == "Chaetognatha"  ~ "CrustOthers", # grouped here because more functionally/taxonomically closer to a crustacean than a gelatinous predator
        phylum == "Mollusca"      ~ "CrustOthers", # grouped here because more functionally/taxonomically closer to a crustacean than a gelatinous predator
        class == "Copepoda"       ~ "CrustOthers",
        class == "Malacostraca"   ~ "CrustOthers",
        # Gelatinous filter-feeders
        class == "Appendicularia" ~ "GelFilter",
        class == "Thaliacea"      ~ "GelFilter",
        # Gelatinous predators
        phylum == "Cnidaria"      ~ "GelPreds",
        phylum == "Ctenophora"    ~ "GelPreds",
        
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
    relocate(sizeGrp, .before = zoopGrp) %>% 
    relocate(funcGrp, .before = sizeGrp)
    
  
  # Count unique ZoopGrps rates
  datClean %>% 
    group_by(zoopGrp) %>% 
    mutate(countZGrp = sum(zoopGrp > 1, na.rm = TRUE)) %>% 
    distinct(zoopGrp, countZGrp) %>% 
    arrange(countZGrp)
    # Mysids              1
    # Decapods            1
    # Annelids            2
    # Molluscs            4
    # Amphipods          29
    # Ctenophores        55
    # Thaliaceans        88
    # Euphausiids       236
    # Cnidarians        241
    # Copepods          379
  
  
  # Count unique functional groups rates
  datClean %>% 
    group_by(funcGrp) %>% 
    mutate(countFuncGrp = sum(funcGrp > 1, na.rm = TRUE)) %>% 
    distinct(funcGrp, countFuncGrp) %>% 
    arrange(countFuncGrp)
    # GelFilter             88
    # GelPreds             296
    # CrustOthers          652
  
  
  # Count unique size groups rates
  datClean %>% 
    group_by(sizeGrp) %>% 
    mutate(countSizeGrp = sum(sizeGrp > 1, na.rm = TRUE)) %>% 
    distinct(sizeGrp, countSizeGrp) %>% 
    arrange(countSizeGrp)
    # Mesoplankton           382
    # Macroplankton          654
  
  
  # Check the temperature range for each zoopGrp
  datClean %>% 
    select(zoopGrp, temp_C, rate_name) %>% 
    group_by(zoopGrp) %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
  # All have sensible ranges for estimating Q10, except for Amphipods, Annelids, Decapods, Molluscs, Mysids, 

# End data cleaning ----

     
    
# Harmonise and prep data for analysis ----
datFinal <- datClean %>%
      
  # Harmonise data with conversion function
  rowwise() %>%
  mutate(.conv = list(
    if(rate_name == "RespirationRate") {
      convert_respiration(rate_value, rate_unit, genus)
      } 
    else {
      list(rate = NA_real_, unit = NA_character_)
      }),
    rate_value_fin = .conv$rate, 
    rate_unit_fin  = .conv$unit) %>%
  ungroup() %>%
      
      
  # Convert to mass-specific rates
  mutate(Cspecific_rate = case_when(
    rate_name == "RespirationRate" & rate_unit_fin == "ulO2/mgC/hr"                   ~ rate_value_fin,
    rate_name == "RespirationRate" & rate_unit_fin == "ulO2/ind/hr" & !is.na(BMC_mg)  ~ rate_value_fin / BMC_mg,
    TRUE ~ NA_real_),
    final_unit = case_when(
      rate_name == "RespirationRate" & !is.na(Cspecific_rate) ~ "ulO2/mgC/hr",
      TRUE ~ rate_unit_fin),
    zoopGrp = as.factor(zoopGrp)) %>%
  select(-.conv, -rate_value_fin, -rate_unit_fin) %>%
  relocate(Cspecific_rate, .after = rate_name) %>% 
  relocate(final_unit, .after = Cspecific_rate)
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

# Save as RDS for later use
# saveRDS(datFinal, "Data/resp_dat.rds")
