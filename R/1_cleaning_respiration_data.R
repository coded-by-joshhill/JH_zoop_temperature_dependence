# Cleaning respiration data
# Josh Hill
# 12/11/25



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
dat <- read_csv("https://www.dropbox.com/scl/fi/itv9vpnu8twxz2fyhmp1u/Resp_dat.csv?rlkey=9fig4vcw2cog4rc4qa6mtj7lj&st=p6gcxo9y&dl=1", 
                skip = 1) %>%
  mutate(ref_no = if_else(is.na(ref_no) | ref_no == "", # if the value is NA or empty...
                          paste0("Hill_", row_number()), # apply a unique reference number
                          ref_no), # otherwise keep what is there
         taxa = str_squish(taxa)) %>% # remove extra spaces from taxon names
  relocate(ref_no, .before = everything())  # move it before all columns
glimpse(dat)


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
      # 22 records

  

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
    select(taxa, Phylum, Class, Order, Family, Genus) %>% # only keep necessary columns
    clean_names(case = "snake")

  

# Join taxa info and harmonize weight data ----
datClean <- dat %>% 
  left_join(taxaDat, by = "taxa") %>% 
  rowwise() %>% 
  mutate(BMC_mg = convert_CW(BM_C, weight_unit)) %>% # harmonize C weight data to mg
  ungroup() %>% 
  relocate(c(phylum, class, order, family, genus), .before = taxa) %>% 
  relocate(BMC_mg, .after = BM_C) %>% 
  mutate(zoopGrp = case_when( # Create custom groupings following Ikeda2014
    order == "Euphausiacea"   ~ "Euphausiacea",
    order == "Amphipoda"      ~ "Amphipoda",
    order == "Decapoda"       ~ "Decapoda",
    order == "Mysida"         ~ "Mysida",
    class == "Copepoda"       ~ "Copepoda",
    phylum == "Mollusca"      ~ "Mollusca",
    phylum == "Chaetognatha"  ~ "Chaetognatha",
    phylum == "Cnidaria"      ~ "Cnidaria",
    phylum == "Ctenophora"    ~ "Ctenophora",
    phylum == "Annelida"      ~ "Annelida",
    class == "Thaliacea"      ~ "Thaliacea",
    class == "Appendicularia" ~ "Appendicularia",
    .default = "OTHER")) %>% 
  relocate(zoopGrp, .before = phylum)

  
  # Count unique ZoopGrps rates
  datClean %>% 
    group_by(zoopGrp) %>% 
    mutate(countZGrp = sum(zoopGrp > 1, na.rm = TRUE)) %>% 
    distinct(zoopGrp, countZGrp) %>% 
    arrange(countZGrp)
      # Mysida               1
      # Decapoda             1
      # Annelida             2
      # Mollusca             4
      # Amphipoda           29
      # Ctenophora          55
      # Thaliacea           88
      # Euphausiacea       236
      # Cnidaria           241
      # Copepoda           290
  
  
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
      TRUE ~ rate_unit_fin)) %>%
  #select(-.conv, -rate_value_fin, -rate_unit_fin) %>% 
  relocate(Cspecific_rate, .after = rate_name) %>% 
  relocate(final_unit, .after = Cspecific_rate)
glimpse(datFinal)

# End conversion

  # Final check
  datFinal %>% 
    ggplot() + 
    geom_point(aes(x = temp_C, 
                   y = log(Cspecific_rate), 
                   colour = zoopGrp))

# Save as RDS for later use
# saveRDS(datFinal, "Data/resp_dat.rds")
