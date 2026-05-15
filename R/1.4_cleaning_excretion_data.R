# Cleaning excretion data
# Josh Hill
# 14/05/26



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
dat <- read_csv("https://www.dropbox.com/scl/fi/wqxedgj0omi8byr9c3nkp/Excret_dat.csv?rlkey=r37upcwvf6qtxt38f6m54l9tv&st=u6v5kek8&dl=1") %>%
  mutate(ref_no = if_else(is.na(ref_no) | ref_no == "", # if the value is NA or empty...
                          paste0("Hill_", row_number()), # apply a unique reference number
                          ref_no) # otherwise keep what is there
         ) %>% 
  relocate(ref_no, .before = everything()) %>%  # move it before all columns
  rename(
    taxa = scientificName,
    primRef = primaryReference,
    primRef_URL = primaryReferenceDOI,
    secRef = secondaryReference,
    secRef_URL = secondaryReferenceDOI,
    rate_value = traitValue,
    rate_name = traitName,
    rate_unit = traitUnit,
    temp_C = assocTemperature
    ) %>% 
  mutate(
    taxa = str_squish(taxa), # remove extra spaces from taxon names
    taxa = case_when(taxa == "Acartia (Acartiura) clausi" ~ "Acartia (Acartiura) clausii", # fix name...
                          TRUE ~ taxa)) %>%  # leave the rest as is...
  # Drop these...Worrms package cannot extact info
  filter(taxa != "Beroe",
         taxa != "Beroe cucumis",
         taxa != "Lucifer")

glimpse(dat)

  
  # Count number of initial pre-cleaned records
  dat %>% group_by(rate_name) %>% 
    summarise(count = n())
  #  excretionRateN  2042


  # Look at all unique taxon
  dat %>% 
    distinct(taxa) %>% 
    arrange(taxa) %>% 
    print(n = "Inf")
 
  
  # Look at all unique primary references for each rate type
  dat %>% 
    distinct(primRef) %>% 
    summarise(count = n())
      # 27 records

  

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
    mutate(BM_dry_mg = sizeAssocValue, # simplify dry mass variable...
           # Note, excretion rates are currently ug ammonia...
           rate_value = rate_value / 1000, # convert them to mg like all other rates...
           rate_unit = "mgN-NH4+/ind/hr", # update the unit to match updated rate values
           ) %>% 
    group_by(primRef, taxa, temp_C, BM_dry_mg) %>% # group the data by study, taxon, and associated temp and carbon mass
    slice_max(rate_value, n = 1, with_ties = FALSE) %>% # slice the maximum rate value from that data
    ungroup() %>% 
    relocate(c(phylum, class, order, family, genus, species), .before = taxa) %>% 
    mutate(
      # Create custom size groupings following Grigoratou et al. 2025 Figure 1
      sizeGrp = case_when(
        # Mesoplankton: 0.2 mm - 20 mm
        phylum == "Chaetognatha"  ~ "Mesoplankton",
        class == "Appendicularia" ~ "Mesoplankton", # grouped here because we only have Oikopleura dioica
        class == "Copepoda"       ~ "Mesoplankton",
        order == "Pteropoda"      ~ "Mesoplankton", # grouped here because they are Pteropods
        # Macroplankton: 20 mm - 200 mm
        phylum == "Annelida"      ~ "Macroplankton", 
        phylum == "Cnidaria"      ~ "Macroplankton",
        phylum == "Ctenophora"    ~ "Macroplankton",
        class == "Hydrozoa"       ~ "Macroplankton",
        class == "Malacostraca"   ~ "Macroplankton",
        class == "Thaliacea"      ~ "Macroplankton",
        # Others - not classified and will be excluded from analyses
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
    relocate(sizeGrp, .before = funcGrp) %>% 
    drop_na(temp_C) %>% 
    filter(funcGrp != "OTHER") # remove other groups...

# End data cleaning ----
glimpse(datClean)
     
    
# Harmonise and prep data for analysis ----
datFinal <- datClean %>%

  # Estimate dry mass based on the % from Kiorboe's 2013 Table 1 
  mutate(BMC_mg = case_when(zoopGrp == "Appendicularians" ~ BM_dry_mg * (10.3/100),
                            zoopGrp == "Chaetognaths" ~ BM_dry_mg * (36.7/100),
                            zoopGrp == "Cnidarians"   ~ BM_dry_mg * (13.2/100),
                            zoopGrp == "Copepods"     ~ BM_dry_mg * (48/100),
                            zoopGrp == "Ctenophores"  ~ BM_dry_mg * (5.1/100),
                            zoopGrp == "Euphausiids"  ~ BM_dry_mg * (41.9/100),
                            zoopGrp == "Thaliaceans"  ~ BM_dry_mg * (10.3/100),
                            zoopGrp == "Amphipods"    ~ BM_dry_mg * (34.5/100),
                            zoopGrp == "Molluscs"     ~ BM_dry_mg * (28.9/100),
                            zoopGrp == "Mysids"       ~ BM_dry_mg * (43.5/100)
                            )
         ) %>% 
      
  # Convert to Cmass-specific rates
  mutate(Cspecific_rate = case_when(
    rate_name == "excretionRateN" & rate_unit == "mgN-NH4+/ind/hr" & !is.na(BMC_mg) ~ rate_value / BMC_mg, # convert to mass-specific
    TRUE ~ NA_real_),
    Cspecific_unit = case_when(
      rate_name == "excretionRateN" & !is.na(Cspecific_rate) ~ "mgN-NH4+/mgC/hr", # update mass-specific units to match the rates...
      TRUE ~ rate_unit),
    zoopGrp = as.factor(zoopGrp)) %>%
    
  # Convert to dry-mass specific rates  
  mutate(DrySpecific_rate = case_when(
    rate_name == "excretionRateN" & rate_unit == "mgN-NH4+/ind/hr"  & !is.na(BM_dry_mg) ~ rate_value / BM_dry_mg, # convert to dryspecific
    TRUE ~ NA_real_),
    DrySpecific_unit = case_when(
      rate_name == "excretionRateN" & !is.na(DrySpecific_rate) ~ "ulO2/mgDry/hr",
      TRUE ~ rate_unit)) %>%    
    
  select(-c(sizeAssocName, sizeAssocValue, sizeAssocUnit)) %>% # tidy up the dataframe
  relocate(Cspecific_rate, .after = rate_name) %>%
  relocate(Cspecific_unit, .after = Cspecific_rate)
glimpse(datFinal)

# End harmonisation and conversion ----

# Final checks
  # sizeGrp
  datFinal %>% 
    ggplot() + 
    geom_point(aes(x = temp_C, 
                   y = log(Cspecific_rate), # mass-specific
                   colour = sizeGrp))
  
  # sizeGrp
  datFinal %>% 
    ggplot() + 
    geom_point(aes(x = temp_C, 
                   y = log(DrySpecific_rate), # absolute
                   colour = sizeGrp))
  
  
  # funcGrp
  datFinal %>% 
    ggplot() + 
    geom_point(aes(x = temp_C, 
                   y = log(Cspecific_rate), # mass-specific
                   colour = funcGrp))
  
 
  
  # zoopGrp
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
  # 1 Appendicularians        21
  # 2 Ctenophores             32
  # 3 Decapods                35
  # 4 Mysids                  67
  # 5 Amphipods               70
  # 6 Cnidarians              81
  # 7 Thaliaceans             98
  # 8 Chaetognaths           132
  # 9 Euphausiids            333
  # 10 Copepods               914
  
  
  # Count unique functional groups rates
  datClean %>% 
    group_by(funcGrp) %>% 
    mutate(countFuncGrp = sum(funcGrp > 1, na.rm = TRUE)) %>% 
    distinct(funcGrp, countFuncGrp) %>% 
    arrange(countFuncGrp)
  # 1 GelFilter            119
  # 2 GelPreds             245
  # 3 Crustaceans         1419
  
  
  # Count unique size groups rates
  datClean %>% 
    group_by(sizeGrp) %>% 
    mutate(countSizeGrp = sum(sizeGrp > 1, na.rm = TRUE)) %>% 
    distinct(sizeGrp, countSizeGrp) %>% 
    arrange(countSizeGrp)
  # 1 Macroplankton          716
  # 2 Mesoplankton          1067
  
  
  # Check the temperature range for each zoopGrp
  datClean %>% 
    select(zoopGrp, temp_C) %>% 
    group_by(zoopGrp) %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
  # All seem fine...

# Save as RDS for later use
# saveRDS(datFinal, "Data/excrete_dat.rds")
