# Taxonomic coverage assessment
# Josh Hill
# 04/02/26


  # Here I assess taxonomic coverage and composition of zooplankton rate processes
  # Create visualisations of taxonomic coverage and composition



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

# Size groups
Sgroup_order <- c("Mesoplankton", "Macroplankton")

# Functional groups
Fgroup_order <- c("Crustaceans", "GelPreds", "GelFilter")

# Zooplankton groups - simple to complex taxonomy
Zgroup_order <- c("Ctenophores", "Cnidarians", "Chaetognaths", "Molluscs", 
                 "Annelids", "Amphipods", "Copepods", "Decapods", "Euphausiids", "Mysids", 
                 "Appendicularians", "Thaliaceans")


# Define my colour palette
mycols <- c("Clearance" = "#66c2a5", "Ingestion" = "#fc8d62", "Growth" = "#8da0cb", "Respiration" = "#e78ac3")



# What proportion of the data is between 0 and 30 degC? ----
# Quick function to select data across each dataset
temp_dat_summary <- function(data, dataset_name) {
  data %>%
    select(temp_C, Cspecific_rate, primRef) %>%
    mutate(Dataset = dataset_name)
}


# Bind all the data into a tibble and summarise temperature data
tempDat_binded <- bind_rows(
  temp_dat_summary(respiration, "respiration"),
  temp_dat_summary(growth, "growth"),
  temp_dat_summary(clearance, "clearance"),
  temp_dat_summary(ingestion, "ingestion"))


# Summarise the data and calculate the proportion of data between 0 and 30 degC
tempDat_binded %>%
  distinct(Cspecific_rate, .keep_all = TRUE) %>% 
  summarise(
    n_total = n(),
    n_temp_range = sum(temp_C >= 0 & temp_C <= 30, na.rm = TRUE),
    prop_temp_range = n_temp_range / n_total * 100)
# Data between 0 and 30 degC
# n_total n_temp_range prop_temp_range
#    2327         2201            94.6



# What proportion of each dataset had the method report? ----
# Quick function to select data across each dataset
method_dat_summary <- function(data, dataset_name) {
  data %>%
    select(ref_no, food_type, food_conc, method, locality) %>%
    mutate(Dataset = dataset_name)
}


# Bind all the data into a tibble and summarise methodology info
methodDat_binded <- bind_rows(
  method_dat_summary(respiration, "respiration"),
  method_dat_summary(growth, "growth"),
  method_dat_summary(clearance, "clearance"),
  method_dat_summary(ingestion, "ingestion"))


# Table S1 ----
# Summarise the experiment method data and calculate the proportions reported and not reported.
methodDat_binded %>% 
  mutate(
    method = if_else(is.na(method), "Not reported", method),
    method = recode(method,
                         "NA"                             = "Not reported",
                         "Bottle incubation"              = "Controlled experiment",
                         "In situ bottle incubation"      = "In situ experiment",
                         "Meta-analysis"                  = "Not reported",
                         "Meta analysis"                  = "Not reported", # has no hypen so is different to above...
                         "Oxygen respirometry"            = "Controlled experiment",
                         "Water bottle method"            = "Controlled experiment",
                         "Water bottle"                   = "Controlled experiment",
                         "Observation chamber"            = "Controlled experiment",
                         "Plankton kreisel tank"          = "Controlled experiment",
                         "Mixed diet feeding experiments" = "Controlled experiment",
                         "Sealed-chamber"                 = "Controlled experiment",
                         "Respirometer chamber"           = "Controlled experiment",
                         "Through-flow system"            = "Controlled experiment",
                         "Feeding tank"                   = "Controlled experiment"
                         )
         ) %>% 
  group_by(method) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop_method = round(n / sum(n) * 100, digit = 0)) %>% 
  arrange(-prop_method)
# Table S2 - Proportion of methods used for original data collection
# method                    n         prop_method
# Not reported           1581       66 
# Controlled experiment   758       32 
# In situ experiment       56        2


# How many total observations?
1581 + 758 + 56
  # 2395

# What is the number of distinct records per rate?
tempDat_binded %>% group_by(Dataset) %>% distinct(primRef) %>% count()
# Dataset         n
# 1 clearance      64
# 2 growth         54
# 3 ingestion      43
# 4 respiration    19

# What is the number of distinct records?
tempDat_binded %>% distinct(primRef) %>% nrow()
# 140


# Table 1 ---
# Taxonomic coverage summary 
# Calculate coverage for each dataset
  # The custom function takes the dataset and assigns a it name, then generates simple taxonomic summaries
clearance_coverage <- summarise_taxonomic_coverage(clearance, "Clearance")
ingestion_coverage <- summarise_taxonomic_coverage(ingestion, "Ingestion")
growth_coverage <- summarise_taxonomic_coverage(growth, "Growth")
respiration_coverage <- summarise_taxonomic_coverage(respiration, "Respiration")


# Combine outputs into a summary table
taxonomic_coverage <- bind_rows(clearance_coverage, ingestion_coverage, growth_coverage, respiration_coverage)
taxonomic_coverage

# Table S2
# Dataset     n_phylas n_classes n_orders n_families n_genera n_species n_observations n_records
# 1 Clearance          5         8       12         27       35        65            618        64
# 2 Ingestion          3         5        7         16       21        44            307        43
# 3 Growth             5         7       15         27       30        44            488        54
# 4 Respiration        4         7       15         39       55       103            982        19

# Total nobs
618 + 307 + 488 + 982 
  # 2395, perfect. matches previous count

# number of observations pre-cleaning were... 
# Clearance   1580
# Ingestion    471
# Growth       686
# Respiration 1036

allDat <- bind_rows(clearance, ingestion, growth, respiration) %>% 
  mutate(primRef = factor(primRef))

allDat %>% # How many total unique groups/obs/records are there?
  distinct(phylum, class, order, family, genus, species, Cspecific_rate, primRef) %>% 
  summarise(
    totalPhylum = n_distinct(phylum),
    totalClass = n_distinct(class),
    totalOrder = n_distinct(order),
    totalFamily = n_distinct(family),
    totalGenera = n_distinct(genus),
    totalSpecies = n_distinct(species))

#   totalPhylum totalClass totalOrder totalFamily totalGenera totalSpecies
#   1           5          9         21          54          82          167



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


# Print summaries
clearance_by_Sgroup
ingestion_by_Sgroup
growth_by_Sgroup
respiration_by_Sgroup


# Combine the summaries and plot it up
all_Sgroups <- bind_rows(clearance_by_Sgroup, ingestion_by_Sgroup, growth_by_Sgroup, respiration_by_Sgroup) %>%
  mutate(dataset = factor(dataset, levels = c("Clearance",  "Ingestion", "Growth", "Respiration")),
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  coord_flip() +
  labs(x = expression(bold("Zooplankton size")),
       y = expression(bold("Number of observations")),
       fill = "Dataset") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(size = 12, face = "bold"))

SgrpFreqPlot

# Save it
ggsave("Output/Figure_1_sizeGrps.pdf", SgrpFreqPlot, width = 160, height = 80, units = "mm", dpi = 300)



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


  # Print summaries
  clearance_by_Fgroup
  ingestion_by_Fgroup
  growth_by_Fgroup
  respiration_by_Fgroup


# Combine the summaries and plot it up
all_Fgroups <- bind_rows(clearance_by_Fgroup, ingestion_by_Fgroup, growth_by_Fgroup, respiration_by_Fgroup) %>%
  mutate(dataset = factor(dataset, levels = c("Clearance",  "Ingestion", "Growth", "Respiration")),
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  scale_x_discrete(labels = c("CrustOthers" = str_wrap("Crustaceans and others", width = 10),
                              "GelPreds" = str_wrap("Gelatinous predators", width = 10),
                              "GelFilter" = str_wrap("Gelatinous filter-feeders", width = 10))) +
  coord_flip() +
  labs(
    x = expression(bold("Functional group")),
    y = expression(bold("Number of observations")),
    fill = "Dataset") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 12, face = "bold"))

FgrpFreqPlot

# Save it
ggsave("Output/Figure_1_funcGrps.pdf", FgrpFreqPlot, width = 160, height = 80, units = "mm", dpi = 300)



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


  # Print summaries
  clearance_by_Zgroup
  ingestion_by_Zgroup
  growth_by_Zgroup
  respiration_by_Zgroup

  
# Combine the summaries and plot it up
all_Zgroups <- bind_rows(clearance_by_Zgroup, ingestion_by_Zgroup, growth_by_Zgroup, respiration_by_Zgroup) %>%
  mutate(dataset = factor(dataset, levels = c("Clearance", "Ingestion", "Growth", "Respiration")),
         zoopGrp = factor(zoopGrp, levels = Zgroup_order))  # use group order

  
# ZoopGrps
ZgrpFreqPlot <- ggplot(all_Zgroups,
                   aes(x = fct_rev(zoopGrp), y = n_observations, fill = dataset)) +
  geom_col(position = position_dodge(width = .9, reverse = TRUE)) +
  geom_text(
    aes(label = n_species),
    position = position_dodge(width = 0.9, reverse = TRUE),
    hjust = -0.5, size = 3.5) +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  coord_flip() +
  labs(
    x = expression(bold("Zooplankton group")),
    y = expression(bold("Number of observations")),
    fill = "Dataset") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "whitesmoke", colour = "black")) +
  facet_wrap(~ dataset, nrow = 1)

ZgrpFreqPlot

# Save it
ggsave("Output/Figure_1_zoopGrps.pdf", ZgrpFreqPlot, width = 160, height = 80, units = "mm", dpi = 300)




# Remove y-axis labels from individual plots
SgrpFreqPlot_no_y <- SgrpFreqPlot + labs(y = NULL)
FgrpFreqPlot_no_y <- FgrpFreqPlot + labs(y = NULL)


# Combine plots
Fig1_combined <- SgrpFreqPlot_no_y + FgrpFreqPlot_no_y + ZgrpFreqPlot +
  plot_layout(
    ncol = 1,
    guides = "collect") &
  theme(legend.position = "top")

Fig1_combined

ggsave("Output/Figure_1_FreqPlots.pdf", Fig1_combined, width = 180, height = 165, units = "mm", dpi = 300)
ggsave("Output/Figure_1_FreqPlots.png", Fig1_combined, width = 180, height = 165, units = "mm", dpi = 300)



