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
    phylum == "Mollusca"      ~ "Mollusca",
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
tableDat <- dat %>%
  filter(rate %in% c(
    "Ingestion", "Clearance", "Respiration",
    "Growth", "HouseProduction", "ExcretionAmmonia", "ExcretionPhosphate"
  )) %>%
  mutate(taxon_id = coalesce(taxa, genus, family, order, class, phylum)) %>%  # Get most specific taxon
  group_by(zoopGrp, rate) %>%
  summarise(
    Q10_summary = if (all(is.na(Q10))) {
      NA_character_
    } else {
      n_count <- sum(!is.na(Q10))
      n_taxa <- n_distinct(taxon_id[!is.na(Q10)])  # Count unique taxa with Q10 values
      paste0(
        round(mean(Q10, na.rm = TRUE), 2), # mean Q10
        if (n_taxa == 1) "<sup>*</sup>" else "", # add asterisk if only 1 unique taxon
        "<br>(", round(min(Q10, na.rm = TRUE), 2), "–", round(max(Q10, na.rm = TRUE), 2), # Q10 variance
        ", n = ", n_count, ")" # number of Q10s per rate
      ) 
    },
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = rate,
    values_from = Q10_summary
  ) %>% 
  select(zoopGrp, Ingestion, Clearance, Respiration, Growth, 
         HouseProduction, ExcretionAmmonia, ExcretionPhosphate) %>% 
  arrange(factor(zoopGrp, levels = c("Copepoda", "Euphausiacea", "Amphipoda", 
                                     "Decapoda", "Mysida", "Chaetognatha", 
                                     "Cnidaria", "Ctenophora", "Mollusca", 
                                     "Thaliacea", "Appendicularia", "OTHER")))


# Create the table
table1 <- tableDat %>%
  gt() %>%
  fmt_markdown(columns = everything()) %>%
  cols_align(align = "center", columns = -zoopGrp) %>%
  cols_label(
    zoopGrp = "Zooplankton Group",
    Ingestion = "Ingestion",
    Clearance = "Clearance",
    Respiration = "Respiration",
    Growth = "Growth",
    HouseProduction = html("House<br>Production"),
    ExcretionAmmonia = html("Ammonia<br>Excretion"),
    ExcretionPhosphate = html("Phosphate<br>Excretion")
  ) %>%
  sub_missing(
    columns = everything(),
    missing_text = "—"
  ) %>%
  tab_footnote(
    footnote = "* Single species observation only (i.e., intraspecific Q10)",
  ) %>%
  tab_source_note(
    source_note = md("Values are mean interspecific Q10 with range in parentheses. *n* = number of Q10 values.")
  ) %>%
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
table1

# Save table ----
gtsave(table1, "Output/Table_1.png", vwidth = 1200, vheight = 700)

  



  ####################### Below is working and dump code for now ########################
  
  # Create plotting dataframe
  pdat <- dat %>% 
    filter(rate %in% c("Ingestion", "Clearance", 
                       "Respiration", 
                       "Growth", "HouseProduction", 
                       "ExcretionAmmonia", "ExcretionPhosphate")) 
  # Box plot
  ggplot(pdat, aes(x = zoopGrp, y = Q10, color = zoopGrp)) +
    geom_point(position = position_jitter(width = 0.15), size = 1.5, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") + 
    facet_wrap(~rate, scales = "free_y") +
    theme_bw() +
    labs(
      title = "Q10 Values by Zooplankton Group and Rate",
      x = "Rate",
      y = "Q10",
      color = "Zooplankton Group"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  
  # Prepare summary (mean Q10 per zoopGrp × rate)
  summary_data <- pdat %>%
    group_by(zoopGrp, rate) %>%
    summarise(
      mean_Q10 = mean(Q10, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot() +
    # Small dots: individual Q10s
    geom_jitter(
      data = pdat,
      aes(x = zoopGrp, y = Q10, color = zoopGrp),
      width = 0.15, size = 1.5, alpha = 0.6
    ) +
    # Large dots: mean Q10 (black)
    geom_point(
      data = summary_data,
      aes(x = zoopGrp, y = mean_Q10),
      color = "black",   # black large dots
      size = 2
    ) +
    facet_wrap(~rate, scales = "fixed") +
    theme_bw() +
    labs(
      title = "Mean Q10 with Individual Measurements by Zooplankton Group and Rate",
      x = "Rate",
      y = "Q10",
      color = "Zooplankton Group"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  