# Cleaning respiration data
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
# dat <- read_csv("https://www.dropbox.com/scl/fi/qsach648tlqobvimd34ra/Historic_Q10_dat.csv?rlkey=pfe18t6oo5plxbxzwcpeh42xx&st=57a1nii7&dl=1") %>%
dat <- read_csv("/Users/jth025/Documents/PhD Local/Local data/Historic_Q10_dat.csv") %>% 
  mutate(ref_no = paste0("Hill_", row_number()),
         taxa = str_squish(taxa)) %>% # create a unique identifier (e.g., Hill_row#)
  relocate(ref_no, .before = everything()) # move it before all columns
glimpse(dat)


  # Look at all unique taxon
  dat %>% 
    distinct(taxa) %>% 
    arrange(taxa) %>% 
    print(n = "Inf")
 
  
  # Look at all unique primary references for each rate type
  dat %>% 
    group_by(Q10Type) %>% 
    distinct(primRef) %>% 
    summarise(count = n())
      # 9 interspecific records
      # 36 intraspecific records

  

# Subset taxa data to get AphiaIDs and classifications ----
taxaDat <- dat %>% 
  select(taxa) %>% 
  distinct(taxa) %>% 
  arrange(taxa)
  
  # Get AphiaIDs
  taxaDatID <- taxaDat %>% 
    mutate(AphiaID = map_int(
      taxa, # extract AphiaIDs for taxa using the worrms package
      possibly(~ { Sys.sleep(0.3); wm_name2id(.x)}, # use Sys.sleep to prevent overloading API
               otherwise = NA_integer_)))  # Set to NA if unable to get AphiaID

    taxaDatID %>% 
      summary() # ensure there are no NAs
    
    # Need to add aphia ids for the broken ones.
    # taxaDatID <- taxaDatID %>%
    #   mutate(AphiaID = if_else(row_number() == 19, "1248", AphiaID))


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
