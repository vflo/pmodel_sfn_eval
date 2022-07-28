# library(sapfluxnetr)
library(rpmodel)
library(tidyverse)
library(lubridate)
library(tsibble)
# library(rgdal)
# library(Metrics)
# library(ghibli)
library(splashTools)
library(xts)
library(medfate)
# library(mgcv)
# library(furrr)
# library(future)
# library(qgam)
# library(GauPro)
# source("dist_merge.R")
library(furrr)
plan('multisession', workers = 4)
options('future.global.maxsize'=2*1024*1024^2)

sapply(list('calc_sfn_aggr_E.R','calc_sfn_aggr_env.R',
            'calc_sfn_scale.R','swvl_aggregation.R',
            'calc_optim_swc.R','do_phydro_params.R',
            'data_processing_functions.R','calc_models.R',
            'stomatal_optimization_functions.R', 
            'pmodel_ecrit.R','photosynthetic_functions.R',
            'hydraulic_functions.R'), source, .GlobalEnv)

#### LOAD DATA LIST ####
path <- "DATA/PREPARED_DATA_LIST"

list_files <- list.files(path, full.names = TRUE)


#precipitation data
precip_monthly_average <- read_csv("DATA/precipitation_monthly_average_sites.csv")

precip_monthly_average <- precip_monthly_average %>% 
  mutate(
    Date = paste0(Date,"01"),
    timestamp = as.Date(Date, format = "%Y%m%d"),
    precip_month = mean*lubridate::days_in_month(timestamp))

#sw_in data
sw_in_monthly_average <- read_csv("DATA/sw_in_monthly_average_sites.csv")

sw_in_monthly_average <- sw_in_monthly_average %>% 
  mutate(
    Date = paste0(Date,"01"),
    timestamp = as.Date(Date, format = "%Y%m%d"),
    sw_in_month = mean/(60*60*24))

#temp data
temp_monthly_average <- read_csv("DATA/temp_monthly_average_sites.csv")

temp_monthly_average <- temp_monthly_average %>% 
  mutate(
    Date = paste0(Date,"01"),
    timestamp = as.Date(Date, format = "%Y%m%d"),
    temp_month = mean-273.15)