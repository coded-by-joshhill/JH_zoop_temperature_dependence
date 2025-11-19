# Mapping experiments
# Josh Hill
# 17/11/25
  
  # Here I read in the data and shapefiles
  # Plot the feeding rates globally
  # Plot the decadal observations of experiments
  # Plot the number of observations of taxa



# Packages and helpers ----
library(tidyverse)
library(ggh4x)
library(terra)
library(sf)
library(rnaturalearth)
library(patchwork)



# Read in the data ----
feedDat <- readRDS("Data/clear_ingest_data.rds")
growDat <- readRDS("Data/grwth_dat.rds")
respDat <- readRDS("Data/resp_dat.rds")
Q10Dat <- readRDS("Data/historicQ10_dat.rds")
  
  # Assign rate dataframes to an object
  dfs <- list(feedDat, growDat, respDat)

  
  
# Join and tidy rate data and select variables of interest ----
mdat <- reduce(dfs, full_join) %>% 
    select(primRef, phylum, class, order, family, genus, taxa, rate_name, locality, lat, lon) %>% 
    drop_na(lat, lon) %>% 
    mutate(rate_name = as_factor(rate_name))
  
  # Assign appropriate names for feeding type
  levels(mdat$rate_name) <- c("Clearance", "Ingestion", "Growth", "Respiration") 
  
  
# Set study area and read in shapefiles ----
  # Available shp files:
    # Add url here
  
  

# Read in data/shapefiles
# oceans <- st_read("/Users/jth025/Documents/PhD Local/Local data/Shapefiles/ne_10m_ocean/ne_10m_ocean.shp") %>%
#   st_transform(crs = "+proj=robin") %>% # transform to Robinson projection
#   st_make_valid()
  

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



# MAINPLOT ----
mainmap <- ggplot() +
  geom_sf(data = graticule, colour = "lightgrey", linewidth = 0.2) + 
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040",
          size = 2) +
  geom_sf(data = Q10pts,
          aes(fill = Q10),
          shape = 25,
          colour = "#404040",
          size = 2) +
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  theme_void() +
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Data type"),
         color = guide_legend(title = "Data type")) 
  
mainmap


# SUBPLOT 1: Europe ----  
Europe <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040",
          size = 2) +
  geom_sf(data = Q10pts,
          aes(fill = Q10),
          shape = 25,
          colour = "#404040",
          size = 2) +
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(-20, 40), ylim = c(35, 80), crs = 4326) +  # Europe bounds
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank())

Europe


# SUBPLOT 2: Northern Oceans
NorthOceans <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040",
          size = 2) +
  geom_sf(data = Q10pts,
          aes(fill = Q10),
          shape = 25,
          colour = "#404040",
          size = 2) +
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(-160, -50), ylim = c(20, 70), crs = 4326) +  # Europe bounds
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank())

NorthOceans


# SUBPLOT 3: NWPacific
NWPacific <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040",
          size = 2) +
  geom_sf(data = Q10pts,
          aes(fill = Q10),
          shape = 25,
          colour = "#404040",
          size = 2) +
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(100, 160), ylim = c(20, 50), crs = 4326) +  # Europe bounds
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank())

NWPacific


# SUBPLOT 4: Southern Ocean
SOcean <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040",
          size = 2) +
  geom_sf(data = Q10pts,
          aes(fill = Q10),
          shape = 25,
          colour = "#404040",
          size = 2) +
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(-90, 160), ylim = c(-80, -30), crs = 4326) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank())

SOcean



# Combine plots
(Figure_1 <- (mainmap | Europe )/
  (NorthOceans | NWPacific) /
  (SOcean) +
  plot_layout(height = c(2, 1.5, 1)) +
  # plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12, face = "bold")))

# ggsave("Output/Figure_1.png", plot = Figure_1, width = 10, height = 10)


# Using patchwork is like spreading frozen butter on fresh bread... leaving this as is for now...

