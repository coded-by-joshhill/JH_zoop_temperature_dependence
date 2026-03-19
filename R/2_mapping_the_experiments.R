# Mapping experiments
# Josh Hill
# 04/02/26
  
  # Here I read in the data and shapefiles
  # Plot the biological rates globally
  # Plot the number of observations of taxa



# Packages and helpers ----
library(tidyverse)
library(sf)
library(rnaturalearth)



# Read in the data ----
# Feeding data
feedDat <- readRDS("Data/clear_ingest_data.rds") 


# Growth data
growDat <- readRDS("Data/grwth_dat.rds") 


# Respiration data
respDat <- readRDS("Data/resp_dat.rds") 
  

# Assign rate data to a single list
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


# Define color palette
mycols <- c("Clearance" = "#66c2a5", 
            "Ingestion" = "#fc8d62", 
            "Growth" = "#8da0cb", 
            "Respiration" = "#e78ac3")


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
Fig1_map <- ggplot() +
  # Add graticlues
  geom_sf(data = graticule, colour = "lightgrey", linewidth = 0.2) + 
  # Add the world shapefile
  geom_sf(data = world, fill = "darkgrey", colour = "white", linewidth = 0.1) +
  # Add the data points
  geom_sf(data = pts,
          aes(fill = rate_name),
          shape = 21,
          colour = "#404040", stroke = .15, alpha = .75,
          size = 1.5) +
  # Add the outline for the shapefile
  geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = mycols) + # apply my colours
  theme_void() +
  theme(legend.position = "none", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        plot.tag = element_text(size = 14, face = "bold")) +
  guides(fill = guide_legend(title = "Data type"),
         color = guide_legend(title = "Data type"))
Fig1_map
# ggsave("Output/Figure_1_map.pdf", Fig1_map, width = 140, units = "mm", dpi = 300)

