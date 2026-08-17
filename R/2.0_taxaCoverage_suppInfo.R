# Taxonomic coverage assessment
# Josh Hill
# 04/02/26


  # Here I assess taxonomic coverage and composition of zooplankton rate processes
  # Create visualisations of taxonomic coverage and composition
  # Generate additional supplementary information



# Packages and helpers ----
library(tidyverse)
library(patchwork)
source("R/0_Helpers.R")



# Read in and filter the cleaned data ----
clearance <- readRDS("Data/clear_ingest_data.rds") %>% 
  filter(rate_name == "ClearanceRate") 

ingestion <- readRDS("Data/clear_ingest_data.rds") %>% 
  filter(rate_name == "IngestionRate") 

growth <- readRDS("Data/grwth_dat.rds") 

respiration <- readRDS("Data/resp_dat.rds") 

excretion <- readRDS("Data/excrete_dat.rds") %>% 
  drop_na(Cspecific_rate) %>% 
  # Create matching columns, which we don't have data for in the excret dataset
  mutate(food_type = NA,
         food_conc = NA,
         method = NA) # Unknown values



# Specify order of groupings ----
# Size groups
Sgroup_order <- c("Mesoplankton", "Macroplankton")

# Functional groups
Fgroup_order <- c("Crustaceans", "GelPreds", "GelFilter")

# Zooplankton groups - simple to complex taxonomy
Zgroup_order <- c("Ctenophores", "Cnidarians", "Chaetognaths", "Molluscs", 
                 "Annelids", "Amphipods", "Copepods", "Decapods", "Euphausiids", "Mysids", 
                 "Appendicularians", "Thaliaceans")



# Define my colour palette for rates ----
mycols <- c("Clearance" = "#66c2a5", "Ingestion" = "#fc8d62", "Growth" = "#8da0cb", "Respiration" = "#e78ac3", "Excretion" = "#ebbb54")



# What proportion of the data are between 0 and 30 degC? ----
# Quick function to select data across each dataset
temp_dat_summary <- function(data, dataset_name) {
  data %>%
    select(temp_C, Cspecific_rate, primRef, primRef_URL) %>%
    mutate(Dataset = dataset_name)
}


# Bind all the data into a tibble and summarise temperature data using our summry fcn
tempDat_binded <- bind_rows(
  temp_dat_summary(respiration, "respiration"),
  temp_dat_summary(growth, "growth"),
  temp_dat_summary(clearance, "clearance"),
  temp_dat_summary(ingestion, "ingestion"),
  temp_dat_summary(excretion, "excretion"))


# Summarise the data and calculate the proportion of data between 0 and 30 degC
tempDat_binded %>%
  distinct(Cspecific_rate, .keep_all = TRUE) %>% 
  summarise(
    n_total = n(),
    n_temp_range = sum(temp_C >= -1 & temp_C <= 30, na.rm = TRUE),
    prop_temp_range = n_temp_range / n_total * 100)
# Data between 0 and 30 degC
# A tibble: 1 × 3
# n_total n_temp_range prop_temp_range
#    3905         3584            91.8



# What is the number of distinct records per rate? ----
tempDat_binded %>% group_by(Dataset) %>% distinct(primRef) %>% count()
# Dataset         n
# 1 clearance      64
# 2 excretion      25
# 3 growth         54
# 4 ingestion      43
# 5 respiration    20
  # Note, some of these across rate records may be identical

# What is the number of distinct records?
tempDat_binded %>% 
  mutate(primRef = case_when(primRef == "Bimstedt1985" ~ "Bamstedt1985", # Fix this spelling...
                             .default = primRef)) %>%
  distinct(primRef, primRef_URL) %>% # take the distinct record based on primary author and URL
  arrange(primRef) %>%
  drop_na(primRef) %>% 
  select(-primRef_URL) %>% 
  print(n = Inf) 
    # 172 distinct records across all datasets

# Build a primary reference table for the supplementaries
primRef_table <- tempDat_binded %>% 
  drop_na(primRef) %>% # drop any NAs
  mutate(primRef = case_when(primRef == "Bimstedt1985" ~ "Bamstedt1985", # Fix this spelling
                             .default = primRef)) %>% 
  distinct(Dataset, primRef, primRef_URL) %>% 
  arrange(Dataset, primRef) %>% 
  group_by(Dataset) %>% # group by the dataset...
  mutate(No = row_number()) %>% # create a number within each dataset
  ungroup() %>% 
  select(No, primRef, Dataset, -primRef_URL) %>% # select the number, the primaryRef and the dataset
  pivot_wider(names_from = Dataset, values_from = primRef) %>% # pivot the data wider by Dataset
  # Adjust the data to match other figures
  relocate(excretion, .after = respiration) %>% 
  relocate(growth, .after = ingestion)
primRef_table

# A function that cleans old suffix and applies tidy ones
cleanSuffix <- function(x) {
  # Strip any existing suffix that follows a 4 digit year
  base <- str_remove(x, "(?<=[0-9]{4})[a-z]$")
  # Count occurrences of each base value
  dup_counts <- table(base[!is.na(base)]) # build frequency table of the stripped base primRefs
  out <- base # copy the base primRefs
  counter <- list() # create a counter list to track iterations for new suffix
  # A for loop to map the stripping of old suffix and adding of new suffix across each dataset/column
  for (i in seq_along(base)) {
    b <- base[i] # Take base primRef
    if (is.na(b)) next # skip if NA
    if (dup_counts[[b]] > 1) { # if the base primRef occurs more than once...
      counter[[b]] <- (counter[[b]] %||% 0) + 1 # increase counter by 1 and...
      out[i] <- paste0(b, letters[counter[[b]]]) # paste a suffix to the end of the base primRef
    } # end if statement
  } # end the loop
  out # give me the final output
} # END OF SUFFIXCLEANER

# Clean primRef table and deal with old suffix from databases
primRef_table_clean <- primRef_table %>% # make a clean table with the existing one...
  mutate(across(-No, cleanSuffix)) # map across each column and apply the suffix cleaner
primRef_table_clean

# Save Table S2: primRef table as a csv ----
write_csv(primRef_table_clean, "Output/primRef_table.csv")



# What proportion of each dataset had the method report? ----
# Quick function to select data across each dataset
method_dat_summary <- function(data, dataset_name) {
  data %>%
    select(ref_no, food_type, food_conc, method) %>%
    mutate(Dataset = dataset_name)
}


# Bind all the data into a tibble and summarise methodology info
methodDat_binded <- bind_rows(
  method_dat_summary(respiration, "respiration"),
  method_dat_summary(growth, "growth"),
  method_dat_summary(clearance, "clearance"),
  method_dat_summary(ingestion, "ingestion"),
  method_dat_summary(excretion, "excretion"))


# Summarise the experimental method data and calculate the proportions reported and not reported.
methodDat_binded %>% 
  mutate(
    method = if_else(is.na(method), "Not reported", method),
    method = recode(method,
                         "NA" = "Not reported",
                         "Bottle incubation" = "Controlled experiment",
                         "In situ bottle incubation" = "In situ experiment",
                         "Meta-analysis" = "Not reported",
                         "Meta analysis" = "Not reported", # has no hypen so is different to above...
                         "Oxygen respirometry" = "Controlled experiment",
                         "Water bottle method" = "Controlled experiment",
                         "Water bottle" = "Controlled experiment",
                         "Observation chamber" = "Controlled experiment",
                         "Plankton kreisel tank" = "Controlled experiment",
                         "Mixed diet feeding experiments" = "Controlled experiment",
                         "Sealed-chamber" = "Controlled experiment",
                         "Respirometer chamber" = "Controlled experiment",
                         "Through-flow system" = "Controlled experiment",
                         "Feeding tank" = "Controlled experiment")) %>% 
  group_by(method) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop_method = round((n / sum(n)) * 100, digits = 2)) %>% 
  arrange(-prop_method)
# Proportion of methods used for original data collection
# method                    n         prop_method
# 1 Not reported           3318        80.3
# 2 Controlled experiment   760        18.4
# 3 In situ experiment       56         1.4
  # There were many without the method reported


# How many total observations?
3318 + 760 + 56
  # 4134



# Taxonomic coverage summary ----
# Calculate coverage for each dataset
  # The custom function takes the dataset and assigns a it name, then generates simple taxonomic summaries
clearance_coverage <- summarise_taxonomic_coverage(clearance, "Clearance")
ingestion_coverage <- summarise_taxonomic_coverage(ingestion, "Ingestion")
growth_coverage <- summarise_taxonomic_coverage(growth, "Growth")
respiration_coverage <- summarise_taxonomic_coverage(respiration, "Respiration")
excretion_coverage <- summarise_taxonomic_coverage(excretion, "Excretion")


# Combine outputs into a summary table
taxonomic_coverage <- bind_rows(clearance_coverage, 
                                ingestion_coverage, 
                                growth_coverage, 
                                respiration_coverage, 
                                excretion_coverage)
# Flip the dataframe structure so rates are columns and taxonomic info are rows...
taxonomic_coverage_long <- taxonomic_coverage %>%
  pivot_longer(cols = -Dataset, names_to = "Groups", values_to = "value") %>%
  mutate(Groups = case_when(Groups == "n_phylas"       ~ "Phyla",
                            Groups == "n_classes"      ~ "Classes",
                            Groups == "n_orders"       ~ "Orders",
                            Groups == "n_families"     ~ "Families",
                            Groups == "n_genera"       ~ "Genera",
                            Groups == "n_species"      ~ "Species",
                            Groups == "n_observations" ~ "Observations",
                            Groups == "n_records"      ~ "Total records")) %>%
  pivot_wider(names_from = Dataset,
              values_from = value) %>%
  mutate(Groups = factor(Groups, levels = c("Phyla", "Classes", "Orders", "Families", "Genera", "Species", "Observations", "Total records"))) %>%
  arrange(Groups)
taxonomic_coverage_long # looks good
# Table S1
# Groups             Clearance Ingestion Growth Respiration Excretion
# 1 Phyla                 5         3      5           4         5
# 2 Classes               8         5      7           7         8
# 3 Orders               12         7     15          15        15
# 4 Families             27        16     27          39        41
# 5 Genera               35        21     30          55        71
# 6 Species              65        44     44         103       110
# 7 Observations        618       307    488         978      1743
# 8 Total records        64        43     54          20        25

# Total nobs
618 + 307 + 488 + 978 + 1743
  # 4134, perfect. matches previous count

# number of observations pre-cleaning were... 
# Clearance   1580
# Ingestion    471
# Growth       686
# Respiration 1036
# Excretion - taken from separate database



# How many total unique groups/obs/records are there? ----
allDat <- bind_rows(clearance, ingestion, growth, respiration, excretion) %>% 
  mutate(primRef = factor(primRef))

allDat %>% 
  distinct(phylum, class, order, family, genus, species, Cspecific_rate, primRef) %>% 
  summarise(
    totalPhylum = n_distinct(phylum),
    totalClass = n_distinct(class),
    totalOrder = n_distinct(order),
    totalFamily = n_distinct(family),
    totalGenera = n_distinct(genus),
    totalSpecies = n_distinct(species))
#   totalPhylum totalClass totalOrder totalFamily totalGenera totalSpecies
#        5          9         22          65         110          216



# Taxonomic composition by zoopGrp ----
# Here we show the distribution of data across sizeGrps, funcGrps and zoopGrps for each rate process...
# and show the number of species per group

## SizeGroups ----
clearance_by_Sgroup <- clearance %>%
  group_by(sizeGrp) %>% 
  summarise(dataset = "Clearance",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(sizeGrp)

ingestion_by_Sgroup <- ingestion %>%
  group_by(sizeGrp) %>% 
  summarise(dataset = "Ingestion",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(sizeGrp)

growth_by_Sgroup <- growth %>%
  group_by(sizeGrp) %>% 
  summarise(dataset = "Growth",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(sizeGrp)

respiration_by_Sgroup <- respiration %>%
  group_by(sizeGrp) %>% 
  summarise(dataset = "Respiration",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  filter(!sizeGrp == "OTHER") %>% # exclude unclassified zooplankton
  arrange(sizeGrp)

excretion_by_Sgroup <- excretion %>%
  group_by(sizeGrp) %>% 
  summarise(dataset = "Excretion",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  filter(!sizeGrp == "OTHER") %>% # exclude unclassified zooplankton
  arrange(sizeGrp)


# Print summaries
clearance_by_Sgroup
ingestion_by_Sgroup
growth_by_Sgroup
respiration_by_Sgroup
excretion_by_Sgroup


# Combine the summaries and plot it up
all_Sgroups <- bind_rows(clearance_by_Sgroup, ingestion_by_Sgroup, growth_by_Sgroup, respiration_by_Sgroup, excretion_by_Sgroup) %>%
  mutate(dataset = factor(dataset, levels = c("Clearance",  "Ingestion", 
                                              "Growth", "Respiration", "Excretion")),
         sizeGrp = factor(sizeGrp, levels = Sgroup_order), # use group order
         sizeGrp = recode(sizeGrp,
                          "Mesoplankton" = "Mesozooplankton",
                          "Macroplankton" = "Macrozooplankton"))  


SgrpFreqPlot <- ggplot(all_Sgroups, 
                       aes(x = fct_rev(sizeGrp), y = n_observations, fill = dataset)) +
  geom_col(position = position_dodge(width = .9, reverse = TRUE)) +
  # Add the number of species annotation
  geom_text(aes(label = n_species),
            position = position_dodge(width = 0.9, reverse = TRUE),
            hjust = -0.5, size = 3.5) +
  scale_fill_manual(values = mycols) +
  scale_y_log10(expand = expansion(mult = c(0, 0.10))) +
  coord_flip() +
  labs(x = expression(bold("Size group")),
       y = expression(bold("Number of observations")),
       fill = "Rate process") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(size = 12, face = "bold"))

SgrpFreqPlot

# Save SizeGrp for Figure S1
ggsave("Output/FigureS1_raw/Figure_S1_sizeGrps.pdf", SgrpFreqPlot, width = 160, height = 80, units = "mm", dpi = 300)



## FuncGroups ----
clearance_by_Fgroup <- clearance %>%
  group_by(funcGrp) %>% 
  filter(n() >= 15, # Exclude funcGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Clearance",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(funcGrp)


ingestion_by_Fgroup <- ingestion %>%
  group_by(funcGrp) %>% 
  filter(n() >= 15, # Exclude funcGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Ingestion",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(funcGrp)


growth_by_Fgroup <- growth %>%
  group_by(funcGrp) %>% 
  filter(n() >= 15, # Exclude funcGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Growth",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(funcGrp)


respiration_by_Fgroup <- respiration %>%
  group_by(funcGrp) %>% 
  filter(n() >= 15, # Exclude funcGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Respiration",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  filter(!funcGrp == "OTHER") %>% # exclude unclassified zooplankton
  arrange(funcGrp)


excretion_by_Fgroup <- excretion %>%
  group_by(funcGrp) %>% 
  filter(n() >= 15, # Exclude funcGrp that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Excretion",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  filter(!funcGrp == "OTHER") %>% # exclude unclassified zooplankton
  arrange(funcGrp)


# Print summaries
clearance_by_Fgroup
ingestion_by_Fgroup
growth_by_Fgroup
respiration_by_Fgroup
excretion_by_Fgroup


# Combine the summaries and plot it up
all_Fgroups <- bind_rows(clearance_by_Fgroup, ingestion_by_Fgroup, growth_by_Fgroup, respiration_by_Fgroup, excretion_by_Fgroup) %>%
  mutate(dataset = factor(dataset, levels = c("Clearance",  "Ingestion", 
                                              "Growth", "Respiration", "Excretion")),
         funcGrp = factor(funcGrp, levels = Fgroup_order))  # use group order


# Plot it up
FgrpFreqPlot <- ggplot(all_Fgroups,
                       aes(x = fct_rev(funcGrp), y = n_observations, fill = dataset)) +
  geom_col(position = position_dodge(width = .9, reverse = TRUE)) +
  geom_text(
    aes(label = n_species),
    position = position_dodge(width = 0.9, reverse = TRUE),
    hjust = -0.5, size = 3.5) +
  scale_fill_manual(values = mycols) +
  scale_y_log10(expand = expansion(mult = c(0, 0.10))) +
  scale_x_discrete(labels = c("CrustOthers" = str_wrap("Crustaceans and others", width = 10),
                              "GelPreds" = str_wrap("Gelatinous predators", width = 10),
                              "GelFilter" = str_wrap("Gelatinous filter-feeders", width = 10))) +
  coord_flip() +
  labs(
    x = expression(bold("Functional group")),
    y = expression(bold("Number of observations")),
    fill = "Rate process") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 12, face = "bold"))
FgrpFreqPlot

# Save FuncGrps for Figure S1
ggsave("Output/FigureS1_raw/Figure_S1_funcGrps.pdf", FgrpFreqPlot, width = 160, height = 80, units = "mm", dpi = 300)



## ZoopGroups (major taxonomic groups) ----
clearance_by_Zgroup <- clearance %>%
  group_by(zoopGrp) %>% 
  filter(n() >= 15, # Exclude zoopGrps that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Clearance",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(zoopGrp)


ingestion_by_Zgroup <- ingestion %>%
  group_by(zoopGrp) %>% 
  filter(n() >= 15, # Exclude zoopGrps that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Ingestion",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(zoopGrp)


growth_by_Zgroup <- growth %>%
  group_by(zoopGrp) %>% 
  filter(n() >= 15, # Exclude zoopGrps that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Growth",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(zoopGrp)


respiration_by_Zgroup <- respiration %>%
  group_by(zoopGrp) %>% 
  filter(n() >= 15, # Exclude zoopGrps that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Respiration",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(zoopGrp)


excretion_by_Zgroup <- excretion %>%
  group_by(zoopGrp) %>% 
  filter(n() >= 15, # Exclude zoopGrps that don't have suitable data or temp ranges
         max(temp_C) - min(temp_C) >= 5) %>%
  summarise(dataset = "Excretion",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(zoopGrp)


# Print summaries
clearance_by_Zgroup
ingestion_by_Zgroup
growth_by_Zgroup
respiration_by_Zgroup
excretion_by_Zgroup

  
# Combine the summaries and plot it up
all_Zgroups <- bind_rows(clearance_by_Zgroup, ingestion_by_Zgroup, growth_by_Zgroup, respiration_by_Zgroup, excretion_by_Zgroup) %>%
  mutate(dataset = factor(dataset, levels = c("Clearance", "Ingestion", 
                                              "Growth", "Respiration", 
                                              "Excretion")),
         zoopGrp = factor(zoopGrp, levels = Zgroup_order))  # use group order

  
# ZoopGrps
ZgrpFreqPlot <- ggplot(all_Zgroups,
                       aes(x = fct_rev(zoopGrp), y = n_observations, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.9, preserve = "single", reverse = TRUE), width = 0.9) +
  geom_text(
    aes(label = n_species),
    position = position_dodge(width = 0.9, reverse = TRUE),
    hjust = -0.3, size = 3) +
  scale_fill_manual(values = mycols) +
  scale_y_log10(expand = expansion(mult = c(0, 0.2))) +
  coord_flip() +
  labs(
    x = expression(bold("Taxonomic group")),
    y = expression(bold("Number of observations")),
    fill = "Rate process") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank())
ZgrpFreqPlot

# Save it
ggsave("Output/FigureS1_raw/Figure_S1_zoopGrps.pdf", ZgrpFreqPlot, width = 160, height = 100, units = "mm", dpi = 300)


# Remove y-axis labels from individual plots
SgrpFreqPlot_no_y <- SgrpFreqPlot + labs(y = NULL)
FgrpFreqPlot_no_y <- FgrpFreqPlot + labs(y = NULL)
