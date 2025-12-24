# Mapping experiments
# Josh Hill
# 17/11/25
  
  # Here I read in the data and shapefiles
  # Plot the feeding rates globally
  # Plot the decadal observations of experiments
  # Plot the number of observations of taxa



# Packages and helpers ----
library(tidyverse)
library(sf)
library(rnaturalearth)
library(patchwork)



# Read in the data ----
feedDat <- readRDS("Data/clear_ingest_data.rds") %>% 
    group_by(zoopGrp) %>% 
    filter(n() >= 20, # Exclude zoopGrps that don't have suitable data or temp ranges
           max(temp_C) - min(temp_C) >= 5) %>% 
    ungroup() %>% 
    drop_na(Cspecific_rate)


growDat <- readRDS("Data/grwth_dat.rds") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)


respDat <- readRDS("Data/resp_dat.rds") %>% 
  group_by(zoopGrp) %>% 
  filter(n() >= 20,
         max(temp_C) - min(temp_C) >= 5) %>% 
  ungroup() %>% 
  drop_na(Cspecific_rate)
  

Q10Dat <- readRDS("Data/historicQ10_dat.rds") %>% 
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
  

  # Assign rate dataframes to an object
  dfs <- list(feedDat, growDat, respDat)

  
  
# Join and tidy rate data and select variables of interest ----
mdat <- reduce(dfs, full_join) %>% 
    select(primRef, phylum, class, order, family, genus, taxa, rate_name, locality, lat, lon) %>% 
    drop_na(lat, lon) %>% 
    mutate(rate_name = as_factor(rate_name))
  
  
  # Assign appropriate names for feeding type
  levels(mdat$rate_name) <- c("Clearance", "Ingestion", "Growth", "Respiration") 
  
  

# Read in shapefiles ----
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_transform(crs = "+proj=robin")

  
# Convert rate data to sf for plotting
pts <- mdat %>% 
  drop_na(lat, lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(crs = "+proj=robin") # transform point data to Robinson projection


# Convert Q10 data to sf for plotting
Q10pts <- Q10Dat %>% 
  select(rate, lat, lon) %>%
  mutate(Q10 = "Q10") %>% 
  drop_na(lat, lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(crs = "+proj=robin") # transform point data to Robinson projection


# Define color palette
mycols <- c("Clearance" = "#66c2a5", 
            "Ingestion" = "#fc8d62", 
            "Growth" = "#8da0cb", 
            "Respiration" = "#e78ac3",
            "Q10" = "#a6d854")


# Create graticules
graticule <- st_graticule(
  lat = seq(-90, 90, by = 30),
  lon = seq(-180, 180, by = 60),
  crs = 4326) %>%
  st_transform(crs = "ESRI:54030")


# Create an outline
outline <- st_graticule(
  lat = c(-89.9, 89.9),
  lon = c(-180, 180),
  crs = 4326
) %>%
  st_transform(crs = "ESRI:54030")



# Main map ----
mainmap <- ggplot() +
  geom_sf(data = graticule, colour = "lightgrey", linewidth = 0.2) + 
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
  theme_void() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        plot.tag = element_text(size = 14, face = "bold")) +
  guides(fill = guide_legend(title = "Data type"),
         color = guide_legend(title = "Data type"))
mainmap
# ggsave("Output/Figure_1.1.png", plot = mainmap, height = 10)



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
# ggsave("Output/Figure_1.2.png", plot = Europe, height = 8)



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
# ggsave("Output/Figure_1.3.png", plot = NorthOceans, height = 8)


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
# ggsave("Output/Figure_1.4.png", plot = NWPacific, height = 6)


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



# Combine plots ----
(Figure_1 <-
    (mainmap | Europe +
       plot_layout(widths = c(.5, 1))) /
    (NorthOceans | NWPacific + 
       plot_layout(widths = c(.5, .5))) /
    (SOcean + 
       plot_layout(widths = c(.8, .1))))

# ggsave("Output/Figure_1.png", plot = Figure_1, height = 8)



# Regional representation of mapped data ----
region_props <- bind_rows(
  mdat %>%
    drop_na(lat, lon) %>%
    mutate(type = rate_name),
  Q10Dat %>%
    drop_na(lat, lon) %>%
    mutate(type = "Q10")) %>%
  mutate(
    region = case_when(
      lon >= -10  & lon <= 30  & lat >= 34  & lat <= 63  ~ "Europe",
      lon >= -100 & lon <= -60 & lat >= 20  & lat <= 70  ~ "NW Atlantic",
      lon >= 110  & lon <= 155 & lat >= 20  & lat <= 46  ~ "NW Pacific",
      lon >= -160  & lon <= -100 & lat >= 20  & lat <= 60  ~ "NE Pacific",
      lon >= -90  & lon <= 160 & lat >= -80 & lat <= -50 ~ "Southern Ocean",
      lon >= 110  & lon <= 180 & lat >= -50 & lat <= 0 ~ "Australia & Oceania",
      TRUE ~ "Other")) %>%
  count(type, region) %>%
  group_by(type) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(type, desc(prop) )

region_props %>% print(n = "Inf")




