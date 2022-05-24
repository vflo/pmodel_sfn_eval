library(sapfluxnetr)
library(rpmodel)
library(tidyverse)
library(lubridate)
library(tsibble)
library(rgdal)
library(Metrics)
library(ghibli)
library(splashTools)
library(xts)
library(medfate)
library(mgcv)
library(qgam)
library(GauPro)
source("dist_merge.R")

sapply(list('calc_sfn_aggr_E.R','calc_sfn_aggr_env.R','calc_sfn_scale.R',
            'swvl_aggregation.R','calc_optim_swc.R','do_phydro_params.R',
            'data_processing_functions.R'), source, .GlobalEnv)

#### LOAD SAPFLUXNET DATA ####
path <- "../0.1.5/RData/plant"
# sfn_meta <- sapfluxnetr::read_sfn_metadata(folder = path, .write_cache = TRUE)
# save(sfn_meta, file = "DATA/sfn_meta.RData")

#load SAPFLUXNET metadata
load(file = "DATA/sfn_meta.RData")

#create table of SFN coordinates
sfn_meta$site_md %>% dplyr::select(si_code,si_lat,si_long) %>% write_csv(file='DATA/sfn_sites_coord.csv')

#load ERA5-land swc
# load(file = "DATA/df_hourly.RData")
# df_swc <- df %>% 
#   pivot_wider(names_from = variable, values_from = value, values_fn = mean) %>% 
#   mutate(TIMESTAMP = lubridate::date(timestamp))
# rm(df)
# gc()

load(file = "DATA/swc_ERA5_land.RData")
df_swc <- swc_ERA5_land %>% 
  mutate(TIMESTAMP = lubridate::date(TIMESTAMP_daylight),
         swvl1 = swc_ERA5_land,
         swvl2 = swc_ERA5_land,
         swvl3 = swc_ERA5_land,
         swvl4 = swc_ERA5_land)
rm(swc_ERA5_land)
gc()

#load FLUXNET sites names and file names
flx_files <- read_csv("DATA/flx_files.csv")

#load PSI
load("DATA/PSI_DF.RData")
PSI_DF <- as_tibble(PSI_DF)

#load SW ERA5
load("DATA/sw_ERA5_18.RData")
sw_ERA5 <- as_tibble(sw_ERA5_18)

#load FAPAR
fapar <-read_csv("DATA/FAPAR_sites.csv")
fapar_noaa <-read_csv("DATA/FAPAR_sites_noaa.csv")

#species data load
# data("SpParamsUS")
# data("SpParamsMED")
# 
# SpParams <- SpParamsUS %>% bind_rows(SpParamsMED)
SpParams <- read_csv("DATA/sp_traits.csv")

#soil data
# load('DATA/clay.RData')
# load('DATA/sand.RData')



