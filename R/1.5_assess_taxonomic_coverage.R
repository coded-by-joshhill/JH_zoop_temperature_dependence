# Taxonomic autocorrelation assessment
# Josh Hill
# 21/11/25


  # Here I assess taxonomic coverage and composition of zooplankton rate processes
  # Create visualisations of taxonomic coverage and composition



# Packages and helpers ----
library(tidyverse)
library(glmmTMB)
library(gt)
library(patchwork)
source("R/0_Helpers.R")



# Read in the data ----
feeding <- readRDS("Data/clear_ingest_data.rds")
respiration <- readRDS("Data/resp_dat.rds")
growth <- readRDS("Data/grwth_dat.rds")



# Taxonomic coverage summary ----
# Assess the taxonomic breadth of each dataset to demonstrate coverage of major zooplankton lineages

# Calculate coverage for each dataset
feeding_coverage <- summarise_taxonomic_coverage(feeding, "Feeding")
respiration_coverage <- summarise_taxonomic_coverage(respiration, "Respiration")
growth_coverage <- summarise_taxonomic_coverage(growth, "Growth")

# Combine into summary table
taxonomic_coverage <- bind_rows(feeding_coverage, respiration_coverage, growth_coverage)

taxonomic_coverage

# Create a table for supplementary materials
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
    n_observations = "Observations"
  ) %>%
  tab_header(
    title = "Taxonomic coverage across rate processes",
    subtitle = "Number of unique taxonomic units represented in each dataset"
  ) %>%
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
    column_labels.border.bottom.color = "black"
  )

coverage_table

# Save the table
gtsave(coverage_table, "Output/taxonomic_coverage_table.png", vwidth = 1000, vheight = 400)



# Taxonomic composition by zoopGrp ----
# Show the distribution of species across major taxonomic groups relative to each rate process


# Calculate for each dataset
feeding_by_group <- feeding %>%
  group_by(zoopGrp) %>% 
  summarise(dataset = "Feeding",
            n_species = n_distinct(taxa),
            n_observations = n(),
            .groups = "drop") %>% 
  arrange(-n_species)
  
respiration_by_group <- respiration %>% 
  group_by(zoopGrp) %>% 
  summarise(dataset = "Respiration",
            n_species = n_distinct(taxa),
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

  # Print summaries
  feeding_by_group
  respiration_by_group
  growth_by_group

  
# Combine the summaries and create visualisation
all_groups <- bind_rows(feeding_by_group, respiration_by_group, growth_by_group) %>%
  mutate(Dataset = factor(dataset, levels = c("Feeding", "Growth", "Respiration")))

# Create visualisation
p_composition <- ggplot(all_groups, aes(x = reorder(zoopGrp, n_species), y = n_species, fill = Dataset)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = n_species), 
            position = position_dodge(width = 0.7), 
            hjust = -0.2, size = 3) +
  scale_fill_manual(values = c("Feeding" = "#66c2a5", "Growth" = "#8da0cb", "Respiration" = "#e78ac3")) +
  coord_flip() +
  labs(
    title = "Taxonomic composition across rate processes",
    subtitle = "Number of species per zooplankton group",
    x = "Zooplankton group",
    y = "Number of species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, colour = "grey30"),
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

p_composition

# Save it
ggsave("Output/taxonomic_composition.png", p_composition, width = 10, height = 8, dpi = 300, bg = "white")





