library(readr)
library(sf)
library(sp)
library(rgdal)
sfn_sites_coord <- read_csv("DATA/sfn_sites_coord.csv")
sfn_sites_coord <- sfn_sites_coord %>% 
  rename(lon = si_long,
         lat = si_lat) %>% 
  dplyr::select(-si_code)

# Define coordinate reference system
prj4string <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
my.projection <- st_crs(prj4string)

# Create sf object
lat_long_sf <- st_as_sf(sfn_sites_coord, coords = c("lon", "lat"), crs = my.projection)
st_crs(lat_long_sf)

plot(lat_long_sf)

# Export shapefile
st_write(lat_long_sf, "DATA/sfn_sites_coord/sfn_sites_coord.shp", driver="ESRI Shapefile")

