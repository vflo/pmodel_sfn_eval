library(readr)
library(sf)
library(sp)
library(rgdal)
EC_sites_coord <- read_csv("DATA/EC_sites_coord.csv")

# Define coordinate reference system
prj4string <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
my.projection <- st_crs(prj4string)

# Create sf object
lat_long_sf <- st_as_sf(EC_sites_coord, coords = c("lon", "lat"), crs = my.projection)
st_crs(lat_long_sf)

plot(lat_long_sf)

# Export shapefile
st_write(lat_long_sf, "DATA/EC_sites_coord/EC_sites_coord.shp", driver="ESRI Shapefile")

