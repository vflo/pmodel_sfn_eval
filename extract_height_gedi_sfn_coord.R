library(raster)
library(sf)
library(stars)


sfn_points <- st_read("DATA/sfn_sites_coord")
r_height <- read_stars("DATA/GEDI/GEDI03_rh100_mean_2019108_2020287_002_02.tif") 
r_height_crs <- r_height %>% st_transform(crs = st_crs(4326))

height_points <- st_extract(r_height, sfn_points, bilinear = TRUE)
