# Old plots

# crossbar plot for synthesises historic Q10s
ggplot(pdat, 
       aes(x = zoopGrp, y = Q10, color = zoopGrp)) +
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




#### Frequency plot examples

respiration %>% 
  mutate(zoopGrp = factor(zoopGrp, levels = Zgroup_order)) %>%   # use group order
  ggplot(aes(x = log(Cspecific_rate), y = fct_rev(zoopGrp))) +
  geom_density_ridges_gradient(
    aes(fill = after_stat(x)), scale = 2) +
  scale_fill_gradientn(
    colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
    name = "Mass-specific\nrespiration rate") +  # Change this
  theme_bw() +
  labs(title = "Respiration",
       x = "ln mass-specific respiration rate",
       y = "Zooplankton group")


# Ridge plot examples
respiration %>% 
  mutate(zoopGrp = factor(zoopGrp, levels = Zgroup_order)) %>%   # use group order
  group_by(zoopGrp) %>%
  mutate(n_species = n_distinct(taxa),
         species_category = cut(n_species, 
                                breaks = c(0, 5, 15, 25, Inf), 
                                labels = c("Low (1-5)", "Medium (6-15)", "High (16-25)", "Very high (25+)"))) %>%
  ggplot(aes(x = log(Cspecific_rate), y = fct_rev(zoopGrp), fill = species_category)) +
  geom_density_ridges(scale = 2, alpha = 0.7) +
  scale_fill_manual(
    values = c("#2166AC", "#92C5DE", "#F4A582", "#B2182B"),
    name = "Number of\nspecies") +
  theme_bw() +
  labs(title = "Carbon-specific respiration rate across all temperatures",
       x = "ln C-specific respiration rate",
       y = "Zooplankton group")