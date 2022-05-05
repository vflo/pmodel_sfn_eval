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

FAPAR_reverse<-dist_merge_mod(sfn_meta_site, foo, 'si_long', 'si_lat', 'si_long', 'si_lat')

FAPAR <- foo %>% 
  left_join(FAPAR_reverse%>% 
              dplyr::select(si_code, y_lat, y_long) %>% 
              rename(si_lat = y_lat,
                     si_long = y_long))

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

FAPAR_reverse<-dist_merge_mod(sfn_meta_site, foo, 'si_long', 'si_lat', 'si_long', 'si_lat')

FAPAR <- foo %>% 
  dplyr::select(-QC) %>% 
  left_join(FAPAR_reverse%>% 
              dplyr::select(si_code, y_lat, y_long) %>% 
              rename(si_lat = y_lat,
                     si_long = y_long))



write_csv(FAPAR, file="DATA/FAPAR_sites_noaa.csv")
