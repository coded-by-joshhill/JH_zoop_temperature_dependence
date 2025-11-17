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


  
# Convert to sf
pts <- mdat %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(crs = "+proj=robin") # transform point data to Robinson projection

# Define color palette
mycols <- c("Clearance" = "#d64545", 
            "Ingestion" = "#6b8e23", 
            "Growth" = "#2e7d9a", 
            "Respiration" = "#9370db")

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
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) +
  theme_void() +
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Feeding type"),
         color = guide_legend(title = "Feeding type")) 
  
mainmap


# SUBPLOT 1: Europe ----  
Europe <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040",
          size = 2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(-20, 40), ylim = c(35, 80), crs = 4326) +  # Europe bounds
  theme_bw() +
  theme(legend.position = "null",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(title = "Feeding type"))

Europe


# SUBPLOT 2: Northern Oceans
NorthOceans <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040",
          size = 2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(-160, -50), ylim = c(20, 70), crs = 4326) +  # Europe bounds
  theme_bw() +
  theme(legend.position = "null",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(title = "Feeding type"))

NorthOceans


# SUBPLOT 3: NWPacific
NWPacific <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040",
          size = 2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(100, 160), ylim = c(20, 50), crs = 4326) +  # Europe bounds
  theme_bw() +
  theme(legend.position = "null",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(title = "Feeding type"))

NWPacific


# SUBPLOT 4: Southern Ocean
SOcean <- ggplot() +
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040",
          size = 2) +
  scale_fill_manual(values = mycols) +
  coord_sf(xlim = c(-120, 180), ylim = c(-80, -30), crs = 4326) +
  theme_bw() +
  theme(legend.position = "null",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(title = "Feeding type"))

SOcean



# Combine plots

mainmap / (NorthOceans | Europe | NWPacific) / (SOcean) +
  plot_layout(heights = c(4, 2, 2), widths = c(1, 1, 1))

