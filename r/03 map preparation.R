### Required packages ###

# install.packages("pacman")
pacman::p_load(ggplot2, spdep, RColorBrewer, graphics, sp, lattice, sf, 
               install = FALSE)

### Loading data ###

estaciones_coords <- read.csv(file = file.path("data", "estaciones_coords.csv"))

# Convert to an sf object
locations <- st_as_sf(data.frame("lng" = estaciones_coords$lng, 
                                 "lat" = estaciones_coords$lat),
                      coords = c("lng", "lat"),
                      crs = 4326)

# Cartography of the Valencian Community
load(file.path("data", "CartoCV.Rdata"))
carto_muni_sf <- st_as_sf(carto_muni)

# Not all stations in estaciones_coords have data in data,
# nor in df_coefs:

load(file.path("outputs", "data.rda"))
load(file.path("outputs", "df_coefs.rda"))

# estaciones_data <- match(unique(data$COD_ESTACION), unique(estaciones_coords$codigo))
estaciones_data <- match(unique(df_coefs$code[df_coefs$pollutant == "O3"]), unique(estaciones_coords$codigo))

# Convert to an sf object
locations <- st_as_sf(data.frame("lng" = estaciones_coords$lng[estaciones_data], 
                                 "lat" = estaciones_coords$lat[estaciones_data]),
                      coords = c("lng", "lat"),
                      crs = 4326)

# Plot
ggplot() +
  geom_sf(data = carto_muni_sf, fill = "#FFF9F0", color = "grey50") +
  geom_sf(data = locations, color = "black", size = 1) +
  theme_minimal()

# Project to a metric CRS
locations_utm <- st_transform(locations, 25830)  # UTM zone 30N

# Euclidean distances
dist <- st_distance(locations_utm)
dist_km <- dist/1000
# save(dist_km, file = file.path("outputs", "dist_km.rda"))
