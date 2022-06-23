library(tidyverse)
library(data.table)
library(R.utils)
source("dist_merge.R")
load(file = "DATA/sfn_meta.RData")
sfn_meta_site <- sfn_meta$site_md %>% 
  dplyr::select(si_code,si_lat,si_long)

coord <-read_csv("DATA/sfn_sites_coord.csv")

#### SW_IN data ####
path <- list.files("DATA/SW_IN_ERA5_HOURLY",full.names = TRUE)

path %>% 
  # split(seq(nrow(.))) %>%
  purrr::map_df(function(x){
  df <-read_csv(x)

  df %>% 
    distinct() %>% 
    mutate(si_long = as.numeric(si_long),
          si_lat = as.numeric(si_lat))->foo

  df_reverse<-dist_merge_mod(sfn_meta_site, foo[c(1:1000),], 'si_long', 'si_lat', 'si_long', 'si_lat')

  temp <- foo %>%  
    left_join(df_reverse%>% 
                dplyr::select(si_code, y_lat, y_long) %>% 
                rename(si_lat = y_lat,
                      si_long = y_long)) %>% 
    dplyr::select(-c(Hour, si_lat,si_long))
    
  return(temp)

  })->df


split(df, df$si_code) %>%
  purrr::map_df(function(x){
    
    x %>%
      dplyr::rename(TIMESTAMP = "Date") %>% 
      mutate(Date = lubridate::date(TIMESTAMP),
             Hour = lubridate::hour(TIMESTAMP),
             Date = case_when(Hour == 0~lag(Date),
                              TRUE~Date)) %>% 
      dplyr::arrange('TIMESTAMP') %>% 
      group_by(si_code,Date) %>% 
      mutate(SW_in_ERA5 = (mean - lag(mean))/(60*60),
             SW_in_ERA5 = case_when(Hour == 1~0,
                                    TRUE ~ SW_in_ERA5)
             ) %>% #the value is moved to the previous hour as it represents the start of each time step
      ungroup() %>%
      mutate(SW_in_ERA5 = lead(SW_in_ERA5)) %>% 
      dplyr::select(TIMESTAMP,si_code,SW_in_ERA5)->sw_ERA5
    
    save(sw_ERA5, file=paste("DATA/SW_IN_ERA5_HOURLY_SITE/", unique(x$si_code),".RData", sep =""))
  })


