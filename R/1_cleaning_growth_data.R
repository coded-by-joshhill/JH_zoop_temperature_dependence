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



# Read in the data ----
dat <- read_csv("https://www.dropbox.com/scl/fi/gdllcg9d1dx1dzf38pckd/Grwth_dat.csv?rlkey=xqayol7mkakxn2fdek5yvkkn8&st=ju2ms1t2&dl=1",
                skip = 1) %>%
  mutate(ref_no = if_else(is.na(ref_no) | ref_no == "", # if the value is NA or empty...
                          paste0("Hill_", row_number()), # apply a unique reference number
                          ref_no), # otherwise keep what is there
         taxa = str_squish(taxa)) %>% # remove extra spaces from taxon names
  relocate(ref_no, .before = everything())  # move it before all columns
glimpse(dat)


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
      # 56 records

  

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
  relocate(BMC_mg, .after = BM_C)

  
  # Count unique classes
  datClean %>% 
    group_by(class) %>% 
    mutate(countClass = sum(class > 1, na.rm = TRUE)) %>% 
    distinct(class, countClass) %>% 
    arrange(countClass)
      # Good data for:
        # Copepoda
        # Malacostraca
        # Tentaculata
        # Scyphozoa
        # Thaliacea and
        # maybe Sagittoidea and Hydrozoa
  
  
    # Check the temperature range for each class
    datClean %>% 
      filter(class == "Copepoda") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(class == "Malacostraca") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(class == "Tentaculata") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(class == "Scyphozoa") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(class == "Thaliacea") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(class == "Sagittoidea") %>% 
      select(temp_C) %>% 
      summary()
    datClean %>% 
      filter(class == "Hydrozoa") %>% 
      select(temp_C) %>% 
      summary()
        # All have sensible ranges for estimating Q10

# End data cleaning ----

     
    
# Prep data for analysis ----
datFinal <- datClean %>%
  # No need to harmoise this data (unless more is added)... so I'll use a simple mutate to maintain naming convention consistency
  mutate(Cspecific_rate = rate_value,
         final_unit = rate_unit) %>% 
  relocate(Cspecific_rate, .after = rate_name) %>% 
  relocate(final_unit, .after = Cspecific_rate)
glimpse(datFinal)

# End conversion

  # Final check
  datFinal %>% 
    ggplot() + 
    geom_point(aes(x = temp_C, 
                   y = log(Cspecific_rate), 
                   colour = primRef))

# saveRDS(datFinal, "Data/grwth_dat.rds")
