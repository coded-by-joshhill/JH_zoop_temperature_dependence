# Taxonomic coverage assessment
# Josh Hill
# 04/02/26


  # Here I assess taxonomic coverage and composition of zooplankton rate processes
  # Create visualisations of taxonomic coverage and composition



# Packages and helpers ----
library(tidyverse)
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



# Taxonomic coverage summary ----
# Calculate coverage for each dataset
  # The custom function takes the dataset and assigns a it name, then generates simple summaries
clearance_coverage <- summarise_taxonomic_coverage(clearance, "Clearance")
ingestion_coverage <- summarise_taxonomic_coverage(ingestion, "Ingestion")
growth_coverage <- summarise_taxonomic_coverage(growth, "Growth")
respiration_coverage <- summarise_taxonomic_coverage(respiration, "Respiration")


# Combine outputs into a summary table
taxonomic_coverage <- bind_rows(clearance_coverage, ingestion_coverage, growth_coverage, respiration_coverage)
taxonomic_coverage

# number of observations pre-cleaning were... 
# Clearance   1580
# Ingestion    471
# Growth       686
# Respiration 1036

# Table S1
# Dataset     n_species n_genera n_families n_orders n_classes n_phylas n_zoopGrps n_observations
# Clearance          74       38         28       12         8        5          7           1186
# Ingestion          58       29         22       10         7        5          7            471
# Growth             48       32         28       16         7        5          9            685
# Respiration       117       65         47       20        10        6         10           1036



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
         sizeGrp = factor(sizeGrp, levels = Sgroup_order))  # use group order


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
    strip.text = element_text(size = 10, face = "bold")) +
  facet_wrap(~ dataset, nrow = 1)

ZgrpFreqPlot

# Save it
ggsave("Output/Figure_1_zoopGrps.pdf", ZgrpFreqPlot, width = 160, height = 80, units = "mm", dpi = 300)



library(patchwork)

# Remove y-axis labels from individual plots
SgrpFreqPlot_no_y <- SgrpFreqPlot + labs(y = NULL)
FgrpFreqPlot_no_y <- FgrpFreqPlot + labs(y = NULL)


# Combine plots
Fig1_combined <- SgrpFreqPlot_no_y + FgrpFreqPlot_no_y + ZgrpFreqPlot +
  plot_layout(
    ncol = 1,
    guides = "collect"
  ) &
  theme(legend.position = "none")

Fig1_combined

ggsave("Output/Figure_1_FreqPlots.pdf", Fig1_combined, width = 180, height = 160, units = "mm", dpi = 300)




# What proportion of the data is between 0 and 30 degC? ----
# Quick function to select data across each dataset
temp_dat_summary <- function(data, dataset_name) {
  data %>%
    select(temp_C) %>%
    mutate(Dataset = dataset_name)
}


# Bind all the data into a tibble
tempDat_binded <- bind_rows(
  temp_dat_summary(respiration, "respiration"),
  temp_dat_summary(growth, "growth"),
  temp_dat_summary(clearance, "clearance"),
  temp_dat_summary(ingestion, "ingestion"))


# Summarise the data and calculate the proportion of data between 0 and 30 degC
tempDat_binded %>%
  summarise(
    n_total = n(),
    n_temp_range = sum(temp_C >= 0 & temp_C <= 30, na.rm = TRUE),
    prop_temp_range = n_temp_range / n_total * 100)
    # 94.3% of the data is between 0 and 30



# What proportion of each dataset had the method report? ----
# Quick function to select data across each dataset
food_dat_summary <- function(data, dataset_name) {
  data %>%
    select(ref_no, food_type, food_conc, method, locality) %>%
    mutate(Dataset = dataset_name)
}


# Bind all the data into a tibble
foodDat_binded <- bind_rows(
  food_dat_summary(respiration, "respiration"),
  food_dat_summary(growth, "growth"),
  food_dat_summary(clearance, "clearance"),
  food_dat_summary(ingestion, "ingestion"))


# Table S1 ----
# Summarise the experiment method data and calculate the proportions reported and not reported.
foodDat_binded %>% 
  mutate(method = if_else(is.na(method), "Not reported", method)) %>%
  group_by(method) %>%
  mutate(method = recode(method,
                    "NA"                   = "Not reported",
                    "Bottle incubation"              = "Controlled experiment",
                    "In situ bottle incubation"      = "In situ experiment",
                    "Meta-analysis"                  = "Not reported",
                    "Meta analysis"                  = "Not reported",
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
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop_method = n / sum(n) * 100) %>% 
  arrange(-prop_method) %>% 
  print(n = "Inf")
      # Table S1 - Proportion of methods used for original data collection
      # method                    n         prop_method
      # 1 Not reported         1915       56.7 
      # 2 Controlled experiment  1311       38.8 
      # 3 In situ experiment     152        4.50




