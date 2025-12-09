# Taxonomic Bias Assessment
# Josh Hill
# 2025-12-09

# This script compares the taxonomic composition of species in the rate process datasets
# against those represented in the 18S genetic database. It uses existing, WoRMS-standardised
# taxonomy from your cleaned RDS datasets to identify where genetic data coverage is lacking.

# Load packages
library(tidyverse)

# ==== 1. READ AND PREPARE DATASETS ====

# Load pre-classified RDS datasets (already includes taxa, genus, family, order, class, phylum)
clearance <- readRDS("Data/clear_ingest_data.rds") %>%
  filter(rate_name == "ClearanceRate", taxa != "Unknown")
ingestion <- readRDS("Data/clear_ingest_data.rds") %>%
  filter(rate_name == "IngestionRate", taxa != "Unknown")
growth <- readRDS("Data/grwth_dat.rds")
respiration <- readRDS("Data/resp_dat.rds")

# Extract unique species and taxonomic classifications
extract_taxonomic_frame <- function(df, dataset_name) {
  df %>%
    select(taxa, species, genus, family, order, class, phylum) %>%
    distinct() %>%
    mutate(dataset = dataset_name)
}

tax_clearance <- extract_taxonomic_frame(clearance, "Clearance")
tax_ingestion <- extract_taxonomic_frame(ingestion, "Ingestion")
tax_growth <- extract_taxonomic_frame(growth, "Growth")
tax_respiration <- extract_taxonomic_frame(respiration, "Respiration")


# Combine them all
tax_datasets <- bind_rows(tax_clearance, tax_ingestion, tax_growth, tax_respiration)

# ==== 2. PARSE GENETIC DATABASE ====

fasta_path <- "/Users/jth025/Documents/PhD Local/Local data/Genetic/MZGfasta-18s__T4000000__o00__C.fasta"
lines <- readLines(fasta_path)
headers <- lines[grepl("^>", lines)]
species_raw <- sapply(headers, function(h) {
  parts <- strsplit(sub("^>", "", h), "\t")[[1]][1]
  name <- strsplit(parts, "__")[[1]][2]
  gsub("_", " ", name)
})
species_genetic <- unique(species_raw)
genetic_df <- tibble(species = species_genetic) %>%
  mutate(genus = word(species, 1))

# ==== 3. MATCH TAXONOMY ====

# Match at multiple taxonomic levels by comparing species and genus
tax_datasets <- tax_datasets %>%
  mutate(
    in_genetic_species = taxa %in% genetic_df$species,
    in_genetic_genus = genus %in% genetic_df$genus
  )

# Summarise genetic representation for each dataset
summary_table <- tax_datasets %>%
  group_by(dataset) %>%
  summarise(
    total_species = n_distinct(taxa),
    species_in_genetic = sum(in_genetic_species),
    genus_in_genetic = sum(in_genetic_genus),
    species_coverage_pct = round(100 * species_in_genetic / total_species, 1),
    genus_coverage_pct = round(100 * genus_in_genetic / total_species, 1)
  )

print(summary_table)

# ==== 4. VISUALISE TAXONOMIC COVERAGE ====

coverage_long <- summary_table %>%
  pivot_longer(cols = c(species_coverage_pct, genus_coverage_pct),
               names_to = "level", values_to = "coverage")

ggplot(coverage_long, aes(x = dataset, y = coverage, fill = level)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("species_coverage_pct" = "#1f77b4",
                               "genus_coverage_pct" = "#e377c2"),
                    labels = c("Species-level", "Genus-level")) +
  labs(
    title = "Taxonomic coverage of genetic database across rate process datasets",
    x = "Rate process dataset",
    y = "Proportion of taxa represented in genetic database (%)",
    fill = "Taxonomic level"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

# ==== 5. EXPORTS ====

write_csv(summary_table, "Output/taxonomic_bias_summary.csv")
ggsave("Output/taxonomic_bias_coverage.png", width = 8, height = 6, dpi = 300, bg = "white")

cat("\nTaxonomic bias assessment complete.\nSummary table and figure saved to Output folder.\n")