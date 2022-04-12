library(tidyverse)
library(data.table)
source("dist_merge.R")
# load(file = "DATA/sfn_meta.RData")
# sfn_meta$site_md %>% filter(si_code %in% c('AUS_WOM', 'FRA_FON', 'FRA_PUE', 'NLD_LOO', 'RUS_FYO', 'USA_UMB_CON','AUS_RIC_EUC_ELE')) %>% dplyr::select(si_code,si_lat,si_long) %>% write_csv(file='DATA/sfn_sites_coord.csv')
sfn_meta_site <- sfn_meta$site_md %>% filter(si_code %in% c('AUS_WOM', 'FRA_FON', 'FRA_PUE', 'NLD_LOO', 'RUS_FYO', 'USA_UMB_CON','AUS_RIC_EUC_ELE'))%>%dplyr::select(si_code,si_lat,si_long)

fapar <-read_csv("DATA/FAPAR_terra.csv")
coord <-read_csv("DATA/sfn_sites_coord.csv")


fapar %>% 
  mutate(coord = sub(".*[[[*]","",`.geo`),
         coord = sub("*]]]}*.","",coord)) %>% 
  separate(coord,c("si_long",'si_lat'),",") %>% 
  dplyr::select(-`.geo`) %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat))->foo


FAPAR<-dist_merge(foo, sfn_meta_site, 'si_lat', 'si_long', 'si_lat', 'si_long')

write_csv(FAPAR, file="DATA/FAPAR.csv")
