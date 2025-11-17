# Cleaning historic Q10 data
# Josh Hill
# 14/11/25



  # Here I read in the data
  # Use worrms package to get AphiaID's and taxon classifications
  # 
  # Convert absolute rates to mass-specific
  # Save as an RDS file



# Packages and helpers ----
library(tidyverse)
library(worrms)
library(janitor)
library(gtsummary)
source("R/0_Helpers.R")



# Read in the data ----
dat <- read_csv("https://www.dropbox.com/scl/fi/qsach648tlqobvimd34ra/Historic_Q10_dat.csv?rlkey=pfe18t6oo5plxbxzwcpeh42xx&st=dws86kwd&dl=1") %>%
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
    
    # Manually add Aphia ids for Ctenophora on row 19
    taxaDatID <- taxaDatID %>%
      mutate(AphiaID = if_else(row_number() == 19, 1248, AphiaID))


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

  

# Join taxa info back into Q10 data ----
datClean <- dat %>% 
  left_join(taxaDat, by = "taxa") %>% 
  relocate(c(phylum, class, order, family, genus), .before = taxa)
  glimpse(datClean)
  
  
  # Count number of unique phylums
  datClean %>% 
    group_by(phylum) %>% 
    mutate(countPhylum = sum(phylum > 1, na.rm = TRUE)) %>% 
    distinct(phylum, countPhylum) %>% 
    arrange(countPhylum)
      # Good data for:
        # Rotifera               1
        # Chaetognatha           1
        # Ctenophora             1
        # Cnidaria              10
        # Chordata              11
        # Mollusca              15
        # Arthropoda           162
  
  # Count number of unique class
  datClean %>% 
    group_by(class) %>% 
    mutate(countClass = sum(class > 1, na.rm = TRUE)) %>% 
    distinct(class, countClass) %>% 
    arrange(countClass)
  # Data for:
    # Thaliacea               1
    # Bivalvia                7
    # Gastropoda              7
    # Scyphozoa               9
    # Appendicularia         10
    # Malacostraca           64
    # Copepoda               98
  

  # Summarise the data by class and phylum
  datClean %>% 
    select(class, Q10, Q10Type) %>% 
    tbl_summary(by = class,
                statistic = list(all_continuous() ~ "{mean} ({sd})"))
  
  datClean %>% 
    select(phylum, Q10, Q10Type) %>% 
    tbl_summary(by = phylum,
                statistic = list(all_continuous() ~ "{mean} ({sd})"))

# End data cleaning ----

  

# saveRDS(datClean, "Data/historicQ10_dat.rds")
