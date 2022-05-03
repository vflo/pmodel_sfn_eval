library(tidyverse)
library(data.table)
library(R.utils)
source("dist_merge.R")
load(file = "DATA/sfn_meta.RData")
#sfn_meta$site_md %>% filter(si_code %in% c('AUS_WOM', 'FRA_FON', 'FRA_PUE', 'NLD_LOO', 'RUS_FYO', 'USA_UMB_CON','AUS_RIC_EUC_ELE')) %>% dplyr::select(si_code,si_lat,si_long) %>% write_csv(file='DATA/sfn_sites_coord.csv')
sfn_meta_site <- sfn_meta$site_md %>% 
  # filter(si_code %in% c('AUS_WOM', 'FRA_FON', 'FRA_PUE', 'NLD_LOO', 'RUS_FYO', 'USA_UMB_CON','AUS_RIC_EUC_ELE'))%>%
  dplyr::select(si_code,si_lat,si_long)

coord <-read_csv("DATA/sfn_sites_coord.csv")

#### MODIS ####
fapar <-read_csv("DATA/FAPAR.csv")

fapar %>% 
  mutate(coord = sub("{\"type\":\"Point\",\"coordinates\":[","",`.geo`,fixed = TRUE),
         coord = sub("]}","",coord)) %>% 
  separate(coord,c("si_long",'si_lat'),",") %>% 
  dplyr::select(-`.geo`) %>% 
  distinct() %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat))->foo

FAPAR<-dist_merge(foo, sfn_meta_site, 'si_lat', 'si_long', 'si_lat', 'si_long')

foo <- FAPAR %>% dplyr::select(-c(si_long,si_lat)) %>% left_join(sfn_meta_site)

FAPAR <- foo %>% dplyr::select(-si_code) %>% left_join(sfn_meta_site)

write_csv(FAPAR, file="DATA/FAPAR_sites.csv")


#### NOAA ####
fapar <-read_csv("DATA/FAPAR_noaa.csv")

fapar %>% 
  mutate(coord = sub("{\"type\":\"Point\",\"coordinates\":[","",`.geo`,fixed = TRUE),
         coord = sub("]}","",coord)) %>% 
  separate(coord,c("si_long",'si_lat'),",") %>% 
  dplyr::select(-`.geo`) %>% 
  distinct() %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat),
         QC = intToBin(FparLai_QC),
         FparLai_QC = str_sub(QC,-2,-1),
         FparLai_QC = strtoi(FparLai_QC, base = 2))->foo

FAPAR<-dist_merge(foo, sfn_meta_site, 'si_lat', 'si_long', 'si_lat', 'si_long')

foo <- FAPAR %>% dplyr::select(-c(si_long,si_lat)) %>% left_join(sfn_meta_site)

FAPAR <- foo %>% dplyr::select(-si_code,-QC) %>% left_join(sfn_meta_site)

write_csv(FAPAR, file="DATA/FAPAR_sites_noaa.csv")
