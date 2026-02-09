# Cleaning clearance and ingestion data
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
dat <- read_csv("https://www.dropbox.com/scl/fi/b1pl3iys1kqtuiwbfd99h/IngClear_dat.csv?rlkey=4ce0im2s1latych74hos8uors&st=7pwjkk5p&dl=1", 
                skip = 1) %>%
  mutate(ref_no = paste0("Hill_", row_number()),
         taxa = str_squish(taxa)) %>% # create a unique identifier (e.g., Hill_row#)
  relocate(ref_no, .before = everything()) # move it before all columns
glimpse(dat)



# Estimate carbon body mass on the basis of body length ----
cMassDat <- dat %>% 
  select(ref_no, primRef, body_length_mm, BM_C, weight_unit) %>% 
  # Filter for records that need carbon mass estimated - larvaceans here
  filter(primRef %in% c("Lombard2009", "Broms2003", "Aguirre2006", "DadonPilosof2023")) %>% 
  arrange(ref_no) %>% 
  mutate(BM_C = calc_BMC(body_length_mm),
         weight_unit = "ug") %>% 
  select(-body_length_mm) 


  # Update dat
  dat <- dat %>%
    left_join(cMassDat, 
              by = c("ref_no", "primRef"), 
              suffix = c("", "_new")) %>%
    mutate(BM_C = coalesce(BM_C_new, BM_C), 
           weight_unit = coalesce(weight_unit_new, weight_unit)) %>%
    select(-ends_with("_new")) %>%
    arrange(as.numeric(str_extract(ref_no, "\\d+")))
  
  
  # Count number of initial pre-cleaned records
  dat %>% group_by(rate_name) %>% 
    summarise(count = n())
    # ClearanceRate  1580
    # IngestionRate   471


  # Look at all unique ClearanceRate taxon
  dat %>% 
    filter(rate_name == "ClearanceRate") %>% 
    distinct(taxa) %>% 
    arrange(taxa) %>% 
    print(n = "Inf")
      # Need to exclude unknown
      # These taxa were recorded as "zooplankton", likely crustaceans, but I think I have collected those data from raw sources anyway
  
  
  # Look at all unique IngestionRate taxon
  dat %>% 
    filter(rate_name == "IngestionRate") %>% 
    distinct(taxa) %>% 
    arrange(taxa) %>% 
    print(n = "Inf")
      # No unknown values to exclude
 
  
  # Look at all unique primary references for each rate type
  dat %>% 
    filter(!taxa == "Unknown") %>% # Exclude unknown zoops
    group_by(rate_name) %>% 
    distinct(primRef, rate_name) %>% 
    summarise(count = n())
      # ClearanceRate    81 records
      # IngestionRate    58 records

  

# Subset taxa data to get AphiaIDs and classifications ----
taxaDat <- dat %>% 
  select(taxa) %>% 
  filter(!taxa %in% "Unknown") %>% # remove the "unknown" species 
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
    filter(!taxa %in% "Unknown") %>% # remove the "unknown" species 
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

  
  # Count unique ZoopGrps rates for ClearanceRate
  datClean %>% 
    filter(rate_name == "ClearanceRate") %>% 
    group_by(zoopGrp) %>% 
    mutate(countZGrp = sum(zoopGrp > 1, na.rm = TRUE)) %>% 
    distinct(zoopGrp, countZGrp) %>% 
    arrange(countZGrp)
    # Chaetognaths             4
    # Ctenophores             58
    # Appendicularians        60
    # Euphausiids            201
    # Copepods               243
    # Thaliaceans            286
    # Cnidarians             334
  
  
  # Count unique functional groups rates for ClearanceRate
  datClean %>% 
    filter(rate_name == "ClearanceRate") %>% 
    group_by(funcGrp) %>% 
    mutate(countFuncGrp = sum(funcGrp > 1, na.rm = TRUE)) %>% 
    distinct(funcGrp, countFuncGrp) %>% 
    arrange(countFuncGrp)
    # GelFilter            346
    # GelPreds             396
    # Crustaceans          444
  
  
  # Count unique size groups rates for ClearanceRate
  datClean %>% 
    filter(rate_name == "ClearanceRate") %>% 
    group_by(sizeGrp) %>% 
    mutate(countSizeGrp = sum(sizeGrp > 1, na.rm = TRUE)) %>% 
    distinct(sizeGrp, countSizeGrp) %>% 
    arrange(countSizeGrp)
    # Mesoplankton           307
    # Macroplankton          879
  
  
  # Check the temperature range for each zoopGrp
  datClean %>% 
    select(zoopGrp, temp_C, rate_name) %>% 
    group_by(zoopGrp) %>% 
    filter(rate_name == "ClearanceRate") %>% 
    summarise(temp_range = paste0(min(temp_C), "-", max(temp_C)))
      # All have sensible ranges for estimating Q10, except for Chaetognaths
  
  
  # Count unique ZoopGrps rates for IngestionRate
  datClean %>% 
    filter(rate_name == "IngestionRate") %>% 
    group_by(zoopGrp) %>% 
    mutate(countZGrp = sum(zoopGrp > 1, na.rm = TRUE)) %>% 
    distinct(zoopGrp, countZGrp) %>% 
    arrange(countZGrp)
    # Appendicularians         6
    # Cnidarians               6
    # Ctenophores             35
    # Chaetognaths            79
    # Thaliaceans             81
    # Euphausiids             91
    # Copepods               173
  
  
  # Count unique functional groups rates for IngestionRate
  datClean %>% 
    filter(rate_name == "IngestionRate") %>% 
    group_by(funcGrp) %>% 
    mutate(countFuncGrp = sum(funcGrp > 1, na.rm = TRUE)) %>% 
    distinct(funcGrp, countFuncGrp) %>% 
    arrange(countFuncGrp)
    # GelFilter             87
    # GelPreds             120
    # Crustaceans          264
  
  
  # Count unique size groups rates for IngestionRate
  datClean %>% 
    filter(rate_name == "IngestionRate") %>% 
    group_by(sizeGrp) %>% 
    mutate(countSizeGrp = sum(sizeGrp > 1, na.rm = TRUE)) %>% 
    distinct(sizeGrp, countSizeGrp) %>% 
    arrange(countSizeGrp)
    # Macroplankton          213
    # Mesoplankton           258
  
  
  # Check the temperature range for each zoopGrp
  datClean %>% 
    select(zoopGrp, temp_C, rate_name) %>% 
    group_by(zoopGrp) %>% 
    filter(rate_name == "IngestionRate") %>% 
    summarise( 
      temp_range = paste0(min(temp_C), "-", max(temp_C)))
      # All have sensible ranges for estimating Q10, except for Cnidarians and Ctenophores

# End data cleaning ----
    
  
    
# Harmonise and prep data for analysis ----
datFinal <- datClean %>%
      
  # Harmonise data with conversion functions
  rowwise() %>%
  mutate(.conv = list(
    if(rate_name == "ClearanceRate") {
      convert_clearance(rate_value, rate_unit)
      } 
    else if(rate_name == "IngestionRate") {
      convert_ingestion(rate_value, rate_unit)
      } 
    else {
      list(rate = NA_real_, unit = NA_character_)
      }),
    rate_value_fin = .conv$rate, 
    rate_unit_fin  = .conv$unit) %>%
  ungroup() %>%
      
  # Convert to mass-specific rates
  mutate(Cspecific_rate = case_when(
    rate_name == "ClearanceRate" & rate_unit_fin == "ml/mgC/hr"                   ~ rate_value_fin,
    rate_name == "ClearanceRate" & rate_unit_fin == "ml/ind/hr" & !is.na(BMC_mg)  ~ rate_value_fin / BMC_mg,
    rate_name == "IngestionRate" & rate_unit_fin == "mgC/mgC/hr"                  ~ rate_value_fin,
    rate_name == "IngestionRate" & rate_unit_fin == "mgC/ind/hr" & !is.na(BMC_mg) ~ rate_value_fin / BMC_mg,
    TRUE ~ NA_real_),
    final_unit = case_when(
      rate_name == "ClearanceRate" & !is.na(Cspecific_rate) ~ "ml/mgC/hr",
      rate_name == "IngestionRate" & !is.na(Cspecific_rate) ~ "mgC/mgC/hr",
      TRUE ~ rate_unit_fin),
    sizeGrp = as.factor(sizeGrp),
    funcGrp = as.factor(funcGrp),  
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
                 colour = sizeGrp)) +
  facet_wrap(~ rate_name, scales ="free")


# funcGrp
datFinal %>% 
  ggplot() + 
  geom_point(aes(x = temp_C, 
                 y = log(Cspecific_rate), 
                 colour = funcGrp)) +
  facet_wrap(~ rate_name, scales ="free")


# zoopGrp
datFinal %>% 
  ggplot() + 
  geom_point(aes(x = temp_C, 
                 y = log(Cspecific_rate), 
                 colour = zoopGrp)) +
  facet_wrap(~ rate_name, scales ="free")

# Save it as an RDS for later use
# saveRDS(datFinal, "Data/clear_ingest_data.rds")
