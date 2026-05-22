# Cleaning respiration data
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
dat <- read_csv("https://www.dropbox.com/scl/fi/itv9vpnu8twxz2fyhmp1u/Resp_dat.csv?rlkey=9fig4vcw2cog4rc4qa6mtj7lj&st=zkw9mqxu&dl=1", 
                skip = 1) %>%
  mutate(ref_no = if_else(is.na(ref_no) | ref_no == "", # if the value is NA or empty...
                          paste0("Hill_", row_number()), # apply a unique reference number
                          ref_no), # otherwise keep what is there
         taxa = str_squish(taxa)) %>% # remove extra spaces from taxon names
  relocate(ref_no, .before = everything()) %>% # move it before all columns
  filter_out(data_type == "Mean") # exclude any mean data

glimpse(dat)

  
  # Count number of initial pre-cleaned records
  dat %>% group_by(rate_name) %>% 
    summarise(count = n())
  # RespirationRate  981


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
      # 20 records

  

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
    rate_value_clean = .conv$rate, 
    rate_unit_clean  = .conv$unit) %>%
  ungroup() %>%
    
  # Estimate dry mass based on the % from Kiorboe's 2013 Table 1 
  mutate(BM_dry_mg = case_when(zoopGrp == "Appendicularians" ~ BMC_mg / (10.3/100),
                               zoopGrp == "Chaetognaths" ~ BMC_mg / (36.7/100),
                               zoopGrp == "Cnidarians"   ~ BMC_mg / (13.2/100),
                               zoopGrp == "Copepods"     ~ BMC_mg / (48/100),
                               zoopGrp == "Ctenophores"  ~ BMC_mg / (5.1/100),
                               zoopGrp == "Euphausiids"  ~ BMC_mg / (41.9/100),
                               zoopGrp == "Thaliaceans"  ~ BMC_mg / (10.3/100))) %>% 
      
  # Convert to mass-specific rates
  mutate(Cspecific_rate = case_when(
    rate_name == "RespirationRate" & rate_unit_clean == "ulO2/mgC/hr" ~ rate_value_clean, # keep as is, already mass-specific
    rate_name == "RespirationRate" & rate_unit_clean == "ulO2/ind/hr" & !is.na(BMC_mg) ~ rate_value_clean / BMC_mg, # convert to mass-specific
    TRUE ~ NA_real_),
    Cspecific_unit = case_when(
      rate_name == "RespirationRate" & !is.na(Cspecific_rate) ~ "ulO2/mgC/hr", # update mass-specific units to match the rates...
      TRUE ~ rate_unit_clean),
    zoopGrp = as.factor(zoopGrp)) %>%
    
  # Convert to dry-mass specific rates  
  mutate(DrySpecific_rate = case_when(
    rate_name == "RespirationRate" & rate_unit_clean == "ulO2/ind/hr"  & !is.na(BM_dry_mg) ~ rate_value_clean / BM_dry_mg, # convert to dryspecific
    rate_name == "RespirationRate" & rate_unit_clean == "ulO2/mgC/hr"  & !is.na(BM_dry_mg) & !is.na(BMC_mg) ~ (rate_value_clean * BMC_mg) / BM_dry_mg, # C to dry
    TRUE ~ NA_real_),
    DrySpecific_unit = case_when(
      rate_name == "RespirationRate" & !is.na(DrySpecific_rate) ~ "ulO2/mgDry/hr",
      TRUE ~ rate_unit_clean)) %>%    
    
  select(-c(.conv, BM_wet, BM_dry, BM_C, weight_unit, weight_calc)) %>% # tidy up the dataframe
  relocate(Cspecific_rate, .after = rate_name) %>%
  relocate(Cspecific_unit, .after = Cspecific_rate) %>% 
  relocate(rate_value_clean, .after = Cspecific_unit) %>% 
  relocate(rate_unit_clean, .after = rate_value_clean)
glimpse(datFinal)

# End harmonisation and conversion ----


# Final checks
  # Check corrected bodymass (at 15DegC) against mass-specific rates
  datFinal %>% filter(rate_unit_clean == "ulO2/ind/hr") %>% 
    ggplot() + 
    geom_point(aes(x = log10(BMC_mg * 2.35^((15-temp_C)/10)), 
                   y = log10(rate_value_clean), # mass-specific
                   colour = zoopGrp),size = 1.5)
  
  # Check corrected bodymass (at 15DegC) against mass-specific rates
  datFinal %>% #filter(rate_unit_clean == "ulO2/ind/hr") %>% 
    ggplot() + 
    geom_point(aes(x = log10(BMC_mg * 2.35^((15-temp_C)/10)), 
                   y = log10(Cspecific_rate), # mass-specific
                   colour = zoopGrp),size = 1)
  
  
  # temp
  datFinal %>% 
    ggplot() + 
    geom_point(aes(x = temp_C, 
                   y = log(Cspecific_rate), # mass-specific
                   colour = zoopGrp))
  
  # Count unique ZoopGrps rates
  datClean %>% 
    group_by(zoopGrp) %>% 
    mutate(countZGrp = sum(zoopGrp > 1, na.rm = TRUE)) %>% 
    distinct(zoopGrp, countZGrp) %>% 
    arrange(countZGrp)
    # Amphipods          25
    # Ctenophores        51
    # Thaliaceans        79
    # Euphausiids       231
    # Cnidarians        240
    # Copepods          355
  
  
  # Count unique functional groups rates
  datClean %>% 
    group_by(funcGrp) %>% 
    mutate(countFuncGrp = sum(funcGrp > 1, na.rm = TRUE)) %>% 
    distinct(funcGrp, countFuncGrp) %>% 
    arrange(countFuncGrp)
    # GelFilter             79
    # GelPreds             291
    # Crustaceans          611
  
  
  # Count unique size groups rates
  datClean %>% 
    group_by(sizeGrp) %>% 
    mutate(countSizeGrp = sum(sizeGrp > 1, na.rm = TRUE)) %>% 
    distinct(sizeGrp, countSizeGrp) %>% 
    arrange(countSizeGrp)
    # Mesoplankton           355
    # Macroplankton          626
  
  
  # Check the temperature range for each zoopGrp
  datClean %>% 
    select(zoopGrp, temp_C, rate_name) %>% 
    group_by(zoopGrp) %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
  # All have sensible ranges for estimating Q10, except for amphipods

# Save as RDS for later use
# saveRDS(datFinal, "Data/resp_dat.rds")
