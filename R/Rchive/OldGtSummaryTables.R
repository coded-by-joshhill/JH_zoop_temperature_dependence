# GT tables 
# Create a table for supplementary materials ----
coverage_table <- taxonomic_coverage %>%
  gt() %>%
  cols_label(
    Dataset = "Rate process",
    n_species = "Species",
    n_genera = "Genera",
    n_families = "Families",
    n_orders = "Orders",
    n_classes = "Classes",
    n_phylas = "Phyla",
    n_zoopGrps = "Zooplankton groups",
    n_preObservations = "Observations pre-data cleaning",
    n_observations = "Observations post-data cleaning") %>%
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