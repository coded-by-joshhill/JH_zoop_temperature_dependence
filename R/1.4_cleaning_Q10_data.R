# Cleaning historic Q10 data
# Josh Hill
# 03/02/26



  # Here I read in the data
  # Use worrms package to get AphiaID's and taxon classifications
  # 
  # Convert absolute rates to mass-specific
  # Save as an RDS file



# Packages and helpers ----
library(tidyverse)
library(janitor)
library(worrms)
source("R/0_Helpers.R")



# Read in the data ----
dat <- read_csv("https://www.dropbox.com/scl/fi/qsach648tlqobvimd34ra/Historic_Q10_dat.csv?rlkey=pfe18t6oo5plxbxzwcpeh42xx&st=odca7bem&dl=1") %>%
  mutate(ref_no = paste0("Hill_", row_number()), # create a unique identifier (e.g., Hill_row#)
         taxa = str_squish(taxa),
         rate = recode(rate,
                        "AmmoniaExcretion" = "ExcretionAmmonia",
                       .default = rate)) %>% 
  relocate(ref_no, .before = everything()) %>%  # move it before all columns
  select(-temp_range_C)
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
  
  
  # Look at the Q10 types and rates
  dat %>% 
    group_by(rate) %>% 
    distinct(rate)

  

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
      summary() # ensure there are no NAs...
      # There's one...
    
      # Manually add Aphia ids for missing data
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
    relocate(c(phylum, class, order, family, genus), .before = taxa) %>% 
    mutate(
      # Create custom size groupings following Grigoratou et al. 2025 Figure 1
      sizeGrp = case_when(
        # Mesoplankton: 0.2 um - 20 mm
        phylum == "Mollusca"       ~ "Mesoplankton", # grouped here because half are bivalvia larvae and half pteropods and 1 unknown mollusc (likely larvae)
        phylum == "Rotifera"       ~ "Mesoplankton",
        class == "Appendicularia"  ~ "Mesoplankton", # grouped here because we only have Oikopleura dioica
        class == "Copepoda"        ~ "Mesoplankton",
        # Macroplankton: 20 mm - 200 mm
        phylum == "Annelida"      ~ "Macroplankton", # grouped here because Tomopteris carpenteri is a larger sp.
        phylum == "Chaetognatha"  ~ "Macroplankton",
        phylum == "Cnidaria"      ~ "Macroplankton",
        phylum == "Ctenophora"    ~ "Macroplankton",
        class == "Hydrozoa"       ~ "Macroplankton",
        class == "Malacostraca"   ~ "Macroplankton",
        class == "Thaliacea"      ~ "Macroplankton",
        .default = "OTHER"),
      
      # Create custom functional groups based on feeding modes
      funcGrp = case_when(
        # Crustaceans and others
        phylum == "Annelida"      ~ "CrustOthers",
        phylum == "Chaetognatha"  ~ "CrustOthers", # grouped here because more functionally/taxonomically closer to a crustacean than a gelatinous predator
        phylum == "Mollusca"      ~ "CrustOthers", # grouped here because more functionally/taxonomically closer to a crustacean than a gelatinous predator
        phylum == "Rotifera"      ~ "CrustOthers",
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
  glimpse(datClean)


  # Count number of unique zoopGrps Q10s
  datClean %>% 
    group_by(zoopGrp) %>% 
    mutate(countZoopGrp = sum(zoopGrp > 1, na.rm = TRUE)) %>% 
    distinct(zoopGrp, countZoopGrp) %>% 
    arrange(countZoopGrp)
        # Rotifers                    1
        # Mysids                      1
        # Chaetognaths                1
        # Ctenophores                 1
        # Thaliaceans                 1
        # Amphipods                   5
        # Cnidarians                 10
        # Appendicularians           10
        # Molluscs                   15
        # Euphausiids                29
        # Decapods                   29
        # Copepods                   98
  
  # Count number of unique "functional" groups Q10s
  datClean %>% 
    group_by(funcGrp) %>% 
    mutate(countFuncGrp = sum(funcGrp > 1, na.rm = TRUE)) %>% 
    distinct(funcGrp, countFuncGrp) %>% 
    arrange(countFuncGrp)
      # GelPreds              11
      # GelFilter             11
      # CrustOthers          179
  
  # Count number of unique size groups Q10s
  datClean %>% 
    group_by(sizeGrp) %>% 
    mutate(countSizeGrp = sum(sizeGrp > 1, na.rm = TRUE)) %>% 
    distinct(sizeGrp, countSizeGrp) %>% 
    arrange(countSizeGrp)
      # Macroplankton           77
      # Mesoplankton           124

# End data cleaning ----

  
# Save the data
# saveRDS(datClean, "Data/historicQ10_dat.rds")
