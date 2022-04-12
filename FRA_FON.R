library(sapfluxnetr)
library(rpmodel)
library(tidyverse)
library(lubridate)
library(tsibble)
library(rgdal)

sapply(list('calc_sfn_aggr_E.R','calc_sfn_aggr_env.R','calc_sfn_scale.R'), source, .GlobalEnv)

#### LOAD SAPFLUXNET DATA ####
# path <- "../0.1.5/RData/plant"
# sfn_meta <- sapfluxnetr::read_sfn_metadata(folder = path, .write_cache = TRUE)
# save(sfn_meta, file = "sfn_meta.RData")
# load(file = "sfn_meta.RData")

# 
# FRA_FON <- read_sfn_data("FRA_FON", folder = path) %>%
#   sfn_metrics(period = "1 hour", .funs = list(~ mean(., na.rm = TRUE)),
#               solar = TRUE, interval = "general") %>%
#   sapfluxnetr::metrics_tidyfier(metadata = sfn_meta,interval = "general")
  
# 
# save(FRA_FON, file = "FRA_FON.RData")
load(file = "FRA_FON.RData")

#### SITE SWC ####
# WGScoor <-  FRA_FON %>% select(long = si_long, lat = si_lat)
# coordinates(WGScoor)=~long+lat
# proj4string(WGScoor)<- CRS("+proj=longlat +datum=WGS84")
# LLcoor<-spTransform(WGScoor,CRS("+proj=longlat"))
# raster::shapefile(LLcoor, "FRA_FON.shp")

load(file = "df_hourly.RData")
df_swc <- df %>% 
  pivot_wider(names_from = variable, values_from = value, values_fn = mean) %>% 
  rename(TIMESTAMP = timestamp)

FRA_FON <- FRA_FON %>% left_join(df_swc, by="TIMESTAMP")

#### SPATIAL AND TEMPORAL DATA AGGREGATION ####
sfn <- calc_sfn_scale(FRA_FON)
env_hour <- calc_sfn_aggr_env(FRA_FON,"1 hour")
sfn_aggr <- calc_sfn_aggr_E(sfn,"1 week")
env <- calc_sfn_aggr_env(FRA_FON,"1 week")
sfn_aggr_day <- calc_sfn_aggr_E(sfn,"1 day")
env_day <- calc_sfn_aggr_env(FRA_FON,"1 day")


#### HOURLY ####
df_hour <- sfn %>% 
  left_join(env_hour)%>% 
  mutate(ppfd = ppfd_in/1e6,
         patm = rpmodel::calc_patm(FRA_FON$si_elev %>% unique()),
         swvl = swvl1*7/289 + swvl2*(28-7)/289 + swvl3*(100-28)/289 + swvl4*(289-100)/289,
         REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE))) %>% 
  filter(!is.na(REW))

pm <- df_hour %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=400, fapar=1, ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = FALSE, meanalpha = 1, elv=FRA_FON$si_elev %>% unique() )
  })

pm_swc <- df_hour %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=400, fapar=1, ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = TRUE,meanalpha = 1, elv=FRA_FON$si_elev %>% unique() )
  })

R =  8.31446261815324 #universal gas constant m3 Pa K−1 mol−1
pm_E <- 1.6*pm$gs*(df_hour$vpd*1e3)*3600 #mol m-2 h-1
pm_E_swc <- 1.6*pm_swc$gs*(df_hour$vpd*1e3)*3600 #mol m-2 h-1
df_res <- df_hour %>% bind_cols(pm_E = pm_E) %>% bind_cols(pm_E_swc = pm_E_swc)%>% 
  mutate(E_stand_mol= E_stand/18.2) #mol m-2 h-1

cor(df_res$E_stand_mol, df_res$pm_E, use = "complete.obs")
df_res %>%
  ggplot(aes(E_stand_mol,pm_E))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_smooth(method= "lm", color = "grey30")+
  theme_classic()

cor(df_res$E_stand_mol, df_res$pm_E_swc, use = "complete.obs")
df_res %>%
  ggplot(aes(E_stand_mol,pm_E_swc))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_smooth(method= "lm", color = "grey30")+
  theme_classic()

df_res[4000:4100,]%>% 
  ggplot()+
  geom_line(aes(TIMESTAMP,E_stand_mol))+
  geom_line(aes(TIMESTAMP,pm_E), color = "red2")+
  geom_line(aes(TIMESTAMP,pm_E_swc), color = "blue")+
  theme_classic()


#### DAILY ####
df_day <- sfn_aggr_day  %>% 
  left_join(env_day)%>% 
  mutate(ppfd = ppfd_in/1e6,
         patm = rpmodel::calc_patm(FRA_FON$si_elev %>% unique()),
         swvl = swvl1*7/289 + swvl2*(28-7)/289 + swvl3*(100-28)/289 + swvl4*(289-100)/289,
         REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE))) %>%  
  filter(!is.na(REW))

pm_day <- df_day %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=400, fapar=1, ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = FALSE,meanalpha = 1, elv=FRA_FON$si_elev %>% unique() )
  })

pm_day_swc <- df_day %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=400, fapar=1, ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = TRUE,meanalpha = 1, elv=FRA_FON$si_elev %>% unique() )
  })

R =  8.31446261815324 #universal gas constant m3 Pa K−1 mol−1
pm_E <- 1.6*pm_day$gs*(df_day$vpd*1e3)*3600 #mol m-2 h-1
pm_E_swc <- 1.6*pm_day_swc$gs*(df_day$vpd*1e3)*3600 #mol m-2 h-1
df_res <- df_day %>% bind_cols(pm_E = pm_E) %>% bind_cols(pm_E_swc = pm_E_swc)%>% 
  mutate(E_stand_mol= E_stand/18.2) #mol m-2 h-1

cor(df_res$E_stand_mol, df_res$pm_E, use = "complete.obs")
df_res %>%
  ggplot(aes(E_stand_mol,pm_E))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_smooth(method= "lm", color = "grey30")+
  theme_classic()

cor(df_res$E_stand_mol, df_res$pm_E_swc, use = "complete.obs")
df_res %>%
  ggplot(aes(E_stand_mol,pm_E_swc))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_smooth(method= "lm", color = "grey30")+
  theme_classic()


df_res%>% 
  ggplot()+
  geom_line(aes(timestamp_aggr,E_stand_mol))+
  geom_line(aes(timestamp_aggr,pm_E), color = "red2")+
  geom_line(aes(timestamp_aggr,pm_E_swc), color = "blue")+
  theme_classic()


#### WEEKLY ####
df <- sfn_aggr %>% 
  left_join(env) %>% 
  mutate(ppfd = ppfd_in/1e6,
         patm = rpmodel::calc_patm(FRA_FON$si_elev %>% unique()),
         swvl = swvl1*7/289 + swvl2*(28-7)/289 + swvl3*(100-28)/289 + swvl4*(289-100)/289,
         REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE))) %>% 
  filter(!is.na(REW))

pm_weekly <- df %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
  rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=400, fapar=1, ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = FALSE, elv=FRA_FON$si_elev %>% unique() )
  })

pm_weekly_swc <- df %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=400, fapar=1, ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = TRUE, elv=FRA_FON$si_elev %>% unique() )
  })

R =  8.31446261815324 #universal gas constant m3 Pa K−1 mol−1
pm_E <- 1.6*pm_weekly$gs*(df$vpd*1e3)*3600 #mol m-2 h-1
pm_E_swc <- 1.6*pm_weekly_swc$gs*(df$vpd*1e3)*3600 #mol m-2 h-1
df_res <- df %>% bind_cols(pm_E = pm_E) %>% bind_cols(pm_E_swc = pm_E_swc)%>% 
  mutate(E_stand_mol= E_stand/18.2) #mol m-2 h-1

cor(df_res$E_stand_mol, df_res$pm_E, use = "complete.obs")
df_res %>%
  ggplot(aes(E_stand_mol,pm_E))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_smooth(method= "lm", color = "grey30")+
  theme_classic()

cor(df_res$E_stand_mol, df_res$pm_E_swc, use = "complete.obs")
df_res %>%
  ggplot(aes(E_stand_mol,pm_E_swc))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_smooth(method= "lm", color = "grey30")+
  theme_classic()


df_res%>% 
  ggplot()+
  geom_line(aes(timestamp_aggr,E_stand_mol))+
  geom_line(aes(timestamp_aggr,pm_E), color = "red2")+
  geom_line(aes(timestamp_aggr,pm_E_swc), color = "blue")+
  geom_line(aes(timestamp_aggr,REW), color = "green")+
  theme_classic()

#He de comprobar unidades!!!! mirar bien si es gas o 




