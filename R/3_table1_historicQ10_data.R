# Table 1 - Historic Q10 values
# Josh Hill
# 18/11/2025



  # Here I read in the historic Q10 data
  # Separate the data into appropriate groupings
  # Synthesise the current state of Q10 values into one concise table
    # Firstly, split the data wider by the types of Q10 (rates) making rates the columns
    # zoopGrp will be the rows



# Packages and helpers ----
library(tidyverse)
library(gt)



# Read in the data ----
dat <- readRDS("Data/historicQ10_dat.rds") %>% 
  filter(!phylum == "Rotifera") %>% 
  mutate(zoopGrp = case_when( # Create custom groupings following Ikeda2014
    order == "Euphausiacea"   ~ "Euphausiacea",
    order == "Amphipoda"      ~ "Amphipoda",
    order == "Decapoda"       ~ "Decapoda",
    order == "Mysidacea"      ~ "Mysida",
    class == "Copepoda"       ~ "Copepoda",
    phylum == "Mollusca"      ~ "Mollusca (larvae)",
    phylum == "Chaetognatha"  ~ "Chaetognatha",
    phylum == "Cnidaria"      ~ "Cnidaria",
    phylum == "Ctenophora"    ~ "Ctenophora",
    class == "Thaliacea"      ~ "Thaliacea",
    class == "Appendicularia" ~ "Appendicularia",
    .default = "OTHER")) %>% 
  relocate(zoopGrp, .before = phylum)
glimpse(dat)



# Table 1: Summary of Q10s ----

# Prepare dataframe for table
# Prepare dataframe for table
tableDat <- dat %>%
  filter(rate %in% c(
    "Clearance", "Ingestion", "Growth", "Respiration",
    "HouseProduction", "ExcretionAmmonia", "ExcretionPhosphate")) %>%
  mutate(taxon_id = coalesce(taxa, genus, family, order, class, phylum)) %>%  # Get most specific taxon
  group_by(zoopGrp, rate) %>%
  summarise(
    Q10_summary = if (all(is.na(Q10))) {
      NA_character_
    } else {
      n_count <- sum(!is.na(Q10))
      n_taxa  <- n_distinct(taxon_id[!is.na(Q10)])  # Count unique taxa with Q10 values
      
      paste0(
        round(mean(Q10, na.rm = TRUE), 2), # mean Q10
        if (n_taxa == 1) "<sup>*</sup>" else "", # add * if only 1 unique taxon
        "<br>(", # open parentheses always
        if (n_count > 1) {
          paste0(
            # Q10 variance... if any
            round(min(Q10, na.rm = TRUE), 2), "–", 
            round(max(Q10, na.rm = TRUE), 2), ", "
          )
        } else {
          ""
        },
        "n = ", n_count, # number of Q10s per rate
        ")")
    },
    .groups = "drop") %>%
  pivot_wider(
    names_from = rate,
    values_from = Q10_summary) %>% 
  select(
    zoopGrp, Clearance, Ingestion, Growth, Respiration, 
    HouseProduction, ExcretionAmmonia, ExcretionPhosphate) %>% 
  arrange(factor(
    zoopGrp,
    levels = c(
      "Amphipoda", "Euphausiacea", "Copepoda", 
      "Decapoda", "Mysida", "Chaetognatha", 
      "Ctenophora", "Cnidaria", "Mollusca (larvae)", 
      "Thaliacea", "Appendicularia", "OTHER")))



# Create the table
table1 <- tableDat %>%
  gt() %>%
  fmt_markdown(columns = everything()) %>%
  cols_align(align = "center", columns = -zoopGrp) %>%
  cols_label(
    zoopGrp = "Zooplankton group",
    Clearance = "Clearance",
    Ingestion = "Ingestion",
    Growth = "Growth",
    Respiration = "Respiration",
    HouseProduction = html("House<br>Production"),
    ExcretionAmmonia = html("Ammonia<br>Excretion"),
    ExcretionPhosphate = html("Phosphate<br>Excretion")) %>%
  sub_missing(
    columns = everything(),
    missing_text = "—") %>%
  
  tab_source_note(
    source_note = md("Values are mean interspecific Q<sub>10</sub> with range in parentheses."))%>% 
  tab_source_note(
    source_note = md("\\* Single species observation only (i.e., intraspecific Q<sub>10</sub>), *n* = number of Q<sub>10</sub> values, – = data deficient.")) %>%
  
  tab_options(
    table.font.size = px(12),
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
table1

# Save table ----
# gtsave(table1, "Output/Table_1.png", vwidth = 1200, vheight = 700)

  

# Plot it up ---- 
  # to visualise the distribution of historic temp dependence estimates
# Custom grouping orders 
group_order <- c("Amphipoda", "Euphausiacea", "Copepoda", "Decapoda", "Mysida", "Chaetognatha",  
                 "Ctenophora", "Cnidaria", "Mollusca (larvae)", "Thaliacea", "Appendicularia")

rate_order <- c("Clearance", "Ingestion", "Growth", "Respiration", 
                "HouseProduction", "ExcretionAmmonia", "ExcretionPhosphate")

# Create plotting dataframe
pdat <- dat %>% 
  filter(rate %in% c("Clearance", "Ingestion", "Growth", "Respiration", 
                     "HouseProduction", "ExcretionAmmonia", "ExcretionPhosphate")) %>% 
  select(rate, zoopGrp, Q10) %>% 
  mutate(rate = fct_relevel(rate, rate_order),
         rate = fct_recode(rate,
                           "House Production"  = "HouseProduction",
                           "Ammonia Excretion" = "ExcretionAmmonia",
                           "Phosphate Excretion" = "ExcretionPhosphate"),
    zoopGrp = fct_relevel(zoopGrp, group_order))



# Prepare summary (mean Q10 per zoopGrp by rate)
summary_data <- pdat %>%
  group_by(zoopGrp, rate) %>%
  summarise(mean_Q10 = mean(Q10, na.rm = TRUE), .groups = "drop")

meanQ10s<- ggplot() +
  geom_jitter(data = pdat,
              aes(x = zoopGrp, y = Q10),
              color = "darkgrey",
              width = 0.15, size = 1.5, alpha = 0.6) +
  geom_point(data = summary_data,
             aes(x = zoopGrp, y = mean_Q10),
             color = "black", 
             size = 2) +
  facet_wrap(~rate, scales = "fixed", ncol = 2) +
  theme_bw() +
  labs(x = expression(bold("Zooplankton group")),
       y = expression(bold("Temperature senstivitiy (Q"[10] *")"))) +
  scale_y_continuous(breaks = seq(1, 8, by = 1)) +
  theme(axis.text = element_text(size = 10),
        legend.position = "right") +
  coord_flip()
meanQ10s

# ggsave("Output/Figure_Supp2.png", plot = meanQ10s, width = 180, height = 220, units = "mm")

  
  