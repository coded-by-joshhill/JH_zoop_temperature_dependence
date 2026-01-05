# Taxonomic autocorrelation assessment
# Josh Hill
# 21/11/25


  # Here I assess taxonomic coverage and composition of zooplankton rate processes
  # Create visualisations of taxonomic coverage and composition



# Packages and helpers ----
library(tidyverse)
library(gt)
library(patchwork)
source("R/0_Helpers.R")



# Read in and filter the cleaned data ----
clearance <- readRDS("Data/clear_ingest_data.rds") %>% 
  filter(rate_name == "ClearanceRate") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)


ingestion <- readRDS("Data/clear_ingest_data.rds") %>% 
  filter(rate_name == "IngestionRate") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)


growth <- readRDS("Data/grwth_dat.rds") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)


respiration <- readRDS("Data/resp_dat.rds") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)



# Taxonomic coverage summary ----
# Calculate coverage for each dataset using custom function
  # The function takes the dataset and assigns a it name, then generates simple summaries
clearance_coverage <- summarise_taxonomic_coverage(clearance, "Clearance")
ingestion_coverage <- summarise_taxonomic_coverage(ingestion, "Ingestion")
growth_coverage <- summarise_taxonomic_coverage(growth, "Growth")
respiration_coverage <- summarise_taxonomic_coverage(respiration, "Respiration")


# Combine outputs into a summary table
taxonomic_coverage <- bind_rows(clearance_coverage, ingestion_coverage, growth_coverage, respiration_coverage)
taxonomic_coverage



# Create a table for supplementary materials ----
coverage_table <- taxonomic_coverage %>%
  gt() %>%
  cols_label(
    Dataset = "Rate Process",
    n_species = "Species",
    n_genera = "Genera",
    n_families = "Families",
    n_orders = "Orders",
    n_classes = "Classes",
    n_phyla = "Phyla",
    n_observations = "Observations") %>%
  cols_align(align = "center", columns = -Dataset) %>%
  tab_options(
    table.font.size = px(17),
    table.width = pct(100),
    table.border.top.style = "solid",
    table.border.top.width = px(2),
    table.border.top.color = "black",
    table.border.bottom.style = "none",
    table.border.bottom.width = px(2),
    table.border.bottom.color = "black",
    table_body.border.bottom.style = "solid", 
    table_body.border.bottom.color = "black", 
    table_body.hlines.style = "none",
    
    heading.border.bottom.style = "solid",
    heading.border.bottom.width = px(2),
    heading.border.bottom.color = "black",
    
    column_labels.font.weight = "bold",
    column_labels.border.top.style = "solid",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.width = px(2),
    column_labels.border.bottom.color = "black")

coverage_table

# Save the table
# gtsave(coverage_table, "Output/Table_S1.png", vwidth = 650, vheight = 500)



# Taxonomic composition by zoopGrp ----
# Show the distribution of species across major taxonomic groups relative to each rate process

# Calculate for each dataset
clearance_by_group <- clearance %>%
  group_by(zoopGrp) %>% 
  summarise(dataset = "Clearance",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(-n_species)

ingestion_by_group <- ingestion %>%
  group_by(zoopGrp) %>% 
  summarise(dataset = "Ingestion",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(-n_species)

growth_by_group <- growth %>% 
  group_by(zoopGrp) %>% 
  summarise(dataset = "Growth",
            n_species = n_distinct(taxa),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(-n_species)
  
respiration_by_group <- respiration %>% 
  group_by(zoopGrp) %>% 
  summarise(dataset = "Respiration",
            n_species = n_distinct(species),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(-n_species)


  # Print summaries
  clearance_by_group
  ingestion_by_group
  respiration_by_group
  growth_by_group

  

# Combine the summaries and plot it up ----
all_groups <- bind_rows(clearance_by_group, ingestion_by_group, respiration_by_group, growth_by_group) %>%
  mutate(Dataset = factor(dataset, levels = c("Clearance", "Ingestion", "Growth", "Respiration")))

  
# Define color palette
mycols <- c("Clearance" = "#66c2a5", "Ingestion" = "#fc8d62", "Growth" = "#8da0cb", "Respiration" = "#e78ac3")
  

# Create plot
p_composition <- ggplot(all_groups, 
                        aes(x = reorder(zoopGrp, n_species), y = n_species, fill = Dataset)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = n_species), 
            position = position_dodge(width = 0.7), 
            hjust = -0.5, size = 3.5) +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip() +
  labs(x = (expression(bold("Zooplankton group"))),
       y = (expression(bold("Number of species")))) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        legend.title = element_text(size = 14, face = "bold"),
        panel.grid.major.y = element_blank())

p_composition

# Save it
ggsave("Output/Figure_Supp1.png", p_composition, width = 10, height = 8, dpi = 300, bg = "white")
