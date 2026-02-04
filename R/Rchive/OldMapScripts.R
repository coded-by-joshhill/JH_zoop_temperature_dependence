# Old maps scripts


# Submap 1: Europe ----  
Europe <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040", alpha = 0.5,
          size = 2) +
  geom_sf(data = Q10pts,
          aes(fill = Q10),
          shape = 21,
          colour = "#404040", alpha = 0.5,
          size = 1.5) +
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(-10, 30), ylim = c(34, 63), crs = 4326) +  # Europe bounds
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank(),
        plot.tag = element_text(size = 14, face = "bold"))
Europe
# ggsave("Output/Figure_1.2.png", plot = Europe, height = 3.5)



# Submap 2: Northern Oceans
NorthOceans <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040", alpha = 0.5,
          size = 2) +
  geom_sf(data = Q10pts,
          aes(fill = Q10),
          shape = 21,
          colour = "#404040", alpha = 0.5,
          size = 1.5) +
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(-130, -50), ylim = c(20, 70), crs = 4326) +  # Europe bounds
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank(),        
        plot.tag = element_text(size = 14, face = "bold"))

NorthOceans
# ggsave("Output/Figure_1.3.png", plot = NorthOceans, height = 3.5)


# Submap 3: NWPacific
NWPacific <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040", alpha = 0.5,
          size = 2) +
  geom_sf(data = Q10pts,
          aes(fill = Q10),
          shape = 21,
          colour = "#404040", alpha = 0.5,
          size = 1.5) +
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(110, 155), ylim = c(20, 46), crs = 4326) +  # Europe bounds
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank(),
        plot.tag = element_text(size = 14, face = "bold"))
NWPacific
# ggsave("Output/Figure_1.4.png", plot = NWPacific, height = 4)


# Submap 4: Southern Ocean
SOcean <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040", alpha = 0.5,
          size = 2) +
  geom_sf(data = Q10pts,
          aes(fill = Q10),
          shape = 21,
          colour = "#404040", alpha = 0.5,
          size = 1.5) +
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(-90, 160), ylim = c(-80, -30), crs = 4326) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        plot.tag = element_text(size = 14, face = "bold")) +
  guides(fill = guide_legend(title = "Data type"),
         color = guide_legend(title = "Data type"))
SOcean
# ggsave("Output/Figure_1.5.png", plot = SOcean, width = 8)



# Combine plots for visual ----
(Figure_1 <-
   (mainmap | Europe +
      plot_layout(widths = c(.5, 1))) /
   (NorthOceans | NWPacific + 
      plot_layout(widths = c(.5, .5))) /
   (SOcean + 
      plot_layout(widths = c(.8, .1))))

# This looks OK... but for now I'll do the illustration in Inkscape...




# Regional representation of mapped data ----
# Take my two datasets and assess the proportion of each rate for each region
region_props <- bind_rows(
  mdat %>% 
    drop_na(lat, lon) %>%
    mutate(type = rate_name),
  Q10Dat %>%
    drop_na(lat, lon) %>%
    mutate(type = "Q10")
) %>%
  mutate(
    region = case_when(
      lon >= -10 & lon <= 30 & lat >= 34 & lat <= 63 ~ "Europe",
      lon >= -100 & lon <= -60 & lat >= 20 & lat <= 70 ~ "NW Atlantic",
      lon >= 110 & lon <= 155 & lat >= 20 & lat <= 46 ~ "NW Pacific",
      lon >= -160 & lon <= -100 & lat >= 20 & lat <= 60 ~ "NE Pacific",
      lon >= -90 & lon <= 160 & lat >= -80 & lat <= -50 ~ "Southern Ocean",
      lon >= 110 & lon <= 180 & lat >= -50 & lat <= 0 ~ "Australia & Oceania",
      TRUE ~ "Other"
    )
  ) %>%
  count(type, region) %>%
  group_by(type) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(type, desc(prop))

region_props %>% print(n = "Inf")
