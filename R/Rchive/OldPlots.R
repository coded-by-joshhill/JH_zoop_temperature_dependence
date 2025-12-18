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