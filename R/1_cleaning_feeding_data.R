# Cleaning clearance and ingestion data
# Josh Hill
# 16/09/25



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



# Read in the data ----
dat <- read_csv("https://www.dropbox.com/scl/fi/b1pl3iys1kqtuiwbfd99h/IngClear_dat.csv?rlkey=4ce0im2s1latych74hos8uors&st=8c3k443n&dl=1", 
                skip = 1) %>%
  mutate(ref_no = paste0("Hill_", row_number()),
         taxa = str_squish(taxa)) %>% # create a unique identifier (e.g., Hill_row#)
  relocate(ref_no, .before = everything()) # move it before all columns
glimpse(dat)


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
    select(taxa, Phylum, Class, Order, Family, Genus) %>% # only keep necessary columns
    clean_names(case = "snake")

  

# Join taxa info and harmonize weight data ----
datClean <- dat %>% 
    filter(!taxa %in% "Unknown") %>% # remove the "unknown" species 
    left_join(taxaDat, by = "taxa") %>% 
    rowwise() %>% 
    mutate(BMC_mg = convert_CW(BM_C, weight_unit)) %>% # harmonize C weight data to mg
    ungroup() %>% 
    relocate(c(phylum, class, order, family, genus), .before = taxa) %>% 
    relocate(BMC_mg, .after = BM_C)

  
  # Count unique classes for ClearanceRate
  datClean %>% 
    filter(rate_name == "ClearanceRate") %>% 
    group_by(class) %>% 
    mutate(countClass = sum(class > 1, na.rm = TRUE)) %>% 
    distinct(class, countClass) %>% 
    arrange(countClass)
      # Good data for:
        # Scyphozoa
        # Thaliacea
        # Copepoda
        # Malacostraca and,
        # maybe Appendicularia and Tentaculata...
  
  
    # Check the temperature range for each class
    datClean %>% 
      filter(rate_name == "ClearanceRate", class == "Scyphozoa") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(rate_name == "ClearanceRate", class == "Thaliacea") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(rate_name == "ClearanceRate", class == "Copepoda") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(rate_name == "ClearanceRate", class == "Malacostraca") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(rate_name == "ClearanceRate", class == "Appendicularia") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(rate_name == "ClearanceRate", class == "Tentaculata") %>% 
      select(temp_C) %>% 
      summary()
      # All have sensible ranges for estimating Q10
      
  
  
  # Count unique classes for IngestionRate
  datClean %>% 
    filter(rate_name == "IngestionRate") %>% 
    group_by(class) %>% 
    mutate(countClass = sum(class > 1, na.rm = TRUE)) %>% 
    distinct(class, countClass) %>% 
    arrange(countClass)
      # Good data for:
        # Copepoda
        # Malacostraca
        # Thaliacea
        # Sagittoidea and,
        # Maybe Tentaculata
  
  
  # Check the temperature range for each class
    datClean %>% 
      filter(rate_name == "IngestionRate", class == "Copepoda") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(rate_name == "IngestionRate", class == "Malacostraca") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(rate_name == "IngestionRate", class == "Thaliacea") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(rate_name == "IngestionRate", class == "Sagittoidea") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(rate_name == "IngestionRate", class == "Tentaculata") %>% 
      select(temp_C) %>% 
      summary()
        # All except Sagittoidea (17-21 degC) and Tentaculata (10 degC only) have sensible ranges for estimating Q10

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
      TRUE ~ rate_unit_fin)) %>%
  select(-.conv, -rate_value_fin, -rate_unit_fin) %>% 
  relocate(Cspecific_rate, .after = rate_name) %>% 
  relocate(final_unit, .after = Cspecific_rate)
glimpse(datFinal)

# End conversion

  # Final check
  datFinal %>% 
    ggplot() + 
    geom_point(aes(x = temp_C, 
                   y = log(Cspecific_rate), 
                   colour = order)) +
    facet_wrap(~ rate_name, scales ="free")

# saveRDS(datFinal, "Data/clear_ingest_data.rds")
