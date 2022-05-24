library(tidyverse)
library(data.table)
library(R.utils)
source("dist_merge.R")
load(file = "DATA/sfn_meta.RData")
sfn_meta_site <- sfn_meta$site_md %>% 
  dplyr::select(si_code,si_lat,si_long)

coord <-read_csv("DATA/sfn_sites_coord.csv")

#### GEE data ####
# df <-read_csv("DATA/Precipitation_monthly_SFN_1990_2018.csv")
# df <-read_csv("DATA/sw_in_monthly.csv")
df <-read_csv("DATA/temp_monthly.csv")

df %>% 
  distinct() %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat))->foo

df_reverse<-dist_merge_mod(sfn_meta_site, foo, 'si_long', 'si_lat', 'si_long', 'si_lat')

temp <- foo %>% 
  left_join(df_reverse%>% 
              dplyr::select(si_code, y_lat, y_long) %>% 
              rename(si_lat = y_lat,
                     si_long = y_long))

write_csv(temp, file="DATA/temp_monthly_average_sites.csv")
