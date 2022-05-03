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
source("dist_merge.R")

sapply(list('calc_sfn_aggr_E.R','calc_sfn_aggr_env.R','calc_sfn_scale.R','swvl_aggregation.R','calc_optim_swc.R','do_phydro_params.R'), source, .GlobalEnv)

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
data("SpParamsUS")
data("SpParamsMED")

SpParams <- SpParamsUS %>% bind_rows(SpParamsMED)

#soil data
load('DATA/clay.RData')
load('DATA/sand.RData')

#### MODEL CALCULATION ####
as.list(flx_files$sfn_sites)[c(-1)] %>%
  purrr::map(function(x){
  #load SAPFLUXNET site and aggregation at daily level
  sfn <- read_sfn_data(x, folder = path) %>%
          sfn_metrics(period = "1 day", .funs = list(~ mean(., na.rm = TRUE)),
                      solar = TRUE, interval = "general") %>%
          sapfluxnetr::metrics_tidyfier(metadata = sfn_meta, interval = "general")

  if (!is.na(sfn$st_soil_depth %>% unique())){sensors_d = sfn$st_soil_depth %>% unique()}else {sensors_d = 0.5}
  #load FLUXNET and extract simplified variables if there is EC and SFN overlap
  meanalpha <- NA
  flx_file <- flx_files[which(flx_files$sfn_sites == x),"flx_sites"]
  if(!is.na(flx_file)){
    flx <- splashTools::readFluxdata(paste0("DATA/EC/", flx_file$flx_sites), sfn$si_elev %>% unique(), sensors_d = sensors_d)
    
    #calculate alpha value to use in the PMODEL swc limitation
    meanalpha <- sum(flx$aet,na.rm = TRUE)/sum(flx$pet,na.rm = TRUE)
    
    #join SAPFLUXNET and FLUXNET datasets
    sfn <- sfn %>% left_join(flx %>% fortify.zoo %>% as_tibble %>% rename(TIMESTAMP = Index),
                             by = "TIMESTAMP")
    
    #remove FLUXNET dataset
    rm(flx_file)
    rm(flx)
    gc()
  }
  
  #join soil water content from ERA5-land
  # df_swc <- dist_merge(df_swc, sfn %>% select(si_code,si_lat,si_long) %>% unique(), 'si_lat', 'si_long', 'si_lat', 'si_long')
  sfn <- sfn %>% 
    left_join(df_swc, by=c("TIMESTAMP","si_code"))
  
  #if there is no column add it filled with NA
  if(!any(colnames(sfn) == "ws")){sfn <- sfn %>% mutate(ws = c(NA))}
  if(!any(colnames(sfn) == "netrad")){sfn <- sfn %>% mutate(netrad = c(NA))}
  if(!any(colnames(sfn) == "sw_in")){sfn <- sfn %>% mutate(sw_in = c(NA))}
  if(!any(colnames(sfn) == "ppfd_in")){sfn <- sfn %>% mutate(ppfd_in = c(NA))}
  if(!any(colnames(sfn) == "vpd")){sfn <- sfn %>% mutate(vpd = c(NA))}
  if(!any(colnames(sfn) == "ta")){sfn <- sfn %>% mutate(ta = c(NA))}
  if(!any(colnames(sfn) == "swc_shallow")){sfn <- sfn %>% mutate(swc_shallow = c(NA))}
  if(!any(colnames(sfn) == "CO2")){sfn <- sfn %>% mutate(CO2 = 400)}
  if(!any(colnames(sfn) == "aet")){sfn <- sfn %>% mutate(aet = NA)}
  if(!any(colnames(sfn) == "cond")){sfn <- sfn %>% mutate(cond = NA)}
  if(!any(colnames(sfn) == "netr")){sfn <- sfn %>% mutate(netr = NA)}
  if(!any(colnames(sfn) == "PPFD")){sfn <- sfn %>% mutate(PPFD = NA)}
  if(!any(colnames(sfn) == "tc")){sfn <- sfn %>% mutate(tc = c(NA))}
  if(!any(colnames(sfn) == "VPD")){sfn <- sfn %>% mutate(VPD = c(NA))}
  
  #use in site swc if present
  if(all(is.na(sfn$swc_shallow))){
    sfn <- sfn %>% 
      mutate(swc_shallow = swvl2)
  }else{
    swc_model <- lm(swc_shallow ~ swvl2,data=sfn)
    sfn <- sfn %>% 
      mutate(swc_shallow = case_when(is.na(swc_shallow) ~ swc_model$coefficients[1]+swc_model$coefficients[2]*swvl2,
                                     !is.na(swc_shallow) ~ swc_shallow))
  }
  
  #join soil properties
  sfn <- sfn %>% 
    left_join(clay %>% dplyr::select(-c(lat,lon)),by='si_code') %>% 
    left_join(sand %>% dplyr::select(-c(lat,lon)),by='si_code') %>% 
    mutate(st_clay_perc = case_when(is.na(st_clay_perc)~ clay,
                                    !is.na(st_clay_perc)~ st_clay_perc),
           st_sand_perc = case_when(is.na(st_sand_perc)~ sand,
                                    !is.na(st_sand_perc)~ st_sand_perc))
  
  #join SW ERA5
  
  sfn <- sfn %>% left_join(sw_ERA5) %>% mutate(sw_ERA5 = sw_ERA5_18/2) # transform into mean daily value (divide by two to account for night value)
  
  ##CALIBRATE SOIL POTENTIAL by adding SWC
  ## THE OPTIMUM SWC value to add
  is_psi_data <- any(unique(sfn$si_code) %in% unique(PSI_DF$id_sfn) &
                       any(lubridate::date(PSI_DF$TIMESTAMP) %in% lubridate::date(sfn$TIMESTAMP)) & 
                       any(PSI_DF$time_psi == "pre-dawn"))
  if(is_psi_data){
    opt_swc <- optim_swc(sfn, PSI_DF, type ="swc")
  }else{
    opt_swc <- c(1,0)
  }
  
  #### FAPAR ####
  fapar_noaa_filled <- fapar_noaa %>% 
    filter(si_code == unique(sfn$si_code),
           FparLai_QC <= 1
    )%>%
    group_by(date,si_code) %>% 
    summarise_all(mean, na.rm = TRUE) %>% 
    dplyr::select(TIMESTAMP = date, Fapar = Fpar,Lai = Lai) %>%
    mutate(Fapar = Fapar*0.001, Lai = Lai*0.001, source = "NOAA")
  
  fapar_filled <- fapar %>% 
    filter(si_code == unique(sfn$si_code),
           FparLai_QC <= 63
           ) %>%
    group_by(date,si_code) %>% 
    summarise_all(mean, na.rm = TRUE) %>% 
    dplyr::select(TIMESTAMP = date, Fapar = Fpar,Lai = Lai) %>%
    mutate(Fapar = Fapar*0.01,Lai = Lai*0.1, source = "MODIS")
  
  fapar_filled <- fapar_noaa_filled %>% bind_rows(fapar_filled)
  
  min_data = min(fapar_filled$TIMESTAMP)
  
  # Build the model
  model_fapar <- qgam(Fapar ~ s(TIMESTAMP, bs = 'tp', k = 180), data = fapar_filled %>% 
                 dplyr::select(TIMESTAMP,Fapar) %>% 
                 mutate(TIMESTAMP = TIMESTAMP %>% 
                          as.numeric()), qu = 0.5)
  model_lai <- qgam(Lai ~ s(TIMESTAMP, bs = 'tp', k = 180), data = fapar_filled %>% 
                        dplyr::select(TIMESTAMP,Lai) %>% 
                        mutate(TIMESTAMP = TIMESTAMP %>% 
                                 as.numeric()), qu = 0.5)

  pred_FAPAR <- predict(model_fapar, sfn %>% 
                          dplyr::select(TIMESTAMP) %>%  
                          mutate(TIMESTAMP = as.numeric(lubridate::ymd(TIMESTAMP)))%>% 
                          unique())
  pred_LAI <- predict(model_lai, sfn %>% 
                          dplyr::select(TIMESTAMP) %>%  
                          mutate(TIMESTAMP = as.numeric(lubridate::ymd(TIMESTAMP)))%>% 
                          unique())

  
  sfn <- sfn %>% 
    left_join(tibble(TIMESTAMP = sfn$TIMESTAMP %>% 
                       unique(), 
                     FAPAR = pred_FAPAR) %>% 
                filter(TIMESTAMP > min_data), 
              by = "TIMESTAMP") %>% 
    left_join(tibble(TIMESTAMP = sfn$TIMESTAMP %>% 
                       unique(),  LAI = pred_LAI) %>% 
                filter(TIMESTAMP > min_data) %>% 
                mutate(LAI = case_when(LAI<0~0,
                                       LAI>=0~LAI)),
              by = "TIMESTAMP") %>% 
    left_join(fapar_filled, by = "TIMESTAMP")%>% 
    mutate(FAPAR = case_when(FAPAR<0~0,
                             FAPAR>=0~FAPAR))
  
  #### PHYDRO PREPARATION ####
  PHYDRO_TRUE <- FALSE
  density <- rpmodel::calc_density_h2o(25, rpmodel::calc_patm(sfn$si_elev %>% unique()))
  visco <- rpmodel::calc_viscosity_h2o(25, rpmodel::calc_patm(sfn$si_elev %>% unique()))
  
  sp_data <- SpParams %>% filter(Name %in% (sfn$pl_species %>% unique()))
  
  #If there is species data -> create par_plant_list to calculate PHYDRO
  if(nrow(sp_data)>0){
    sp_data %>% 
      split(seq(nrow(.))) %>%
      purrr::map(do_phydro_params,sfn = sfn, visco = visco, density = density) %>% bind_rows() -> par_plant_std
    
    par_plant_std <- par_plant_std %>% 
      summarise(v_huber = weighted.mean(v_huber, sp_basal_area),
                height = weighted.mean(height, sp_basal_area),
                LAI = weighted.mean(LAI, sp_basal_area),
                conductivity = weighted.mean(conductivity, sp_basal_area),
                psi50 = weighted.mean(psi50, sp_basal_area),
                b = weighted.mean(b, sp_basal_area))
    
    if(is.na(unique(par_plant_std$height))){par_plant_std$height <- unique(sfn$st_height)}
      
    if(!is.na(par_plant_std$psi50)){PHYDRO_TRUE <- TRUE}
  }
 
  #### SPATIAL AND TEMPORAL DATA AGGREGATION ####
  sfn_agr <- calc_sfn_scale(sfn)
  sfn_aggr_week <- calc_sfn_aggr_E(sfn_agr,"1 week")
  env_week <- calc_sfn_aggr_env(sfn,"1 week")
  env_week <- env_week %>% mutate(CO2 = oce::fillGap(CO2))
  sfn_aggr_month <- calc_sfn_aggr_E(sfn_agr,"1 month")
  env_month <- calc_sfn_aggr_env(sfn,"1 month")
  env_month <- env_month %>% mutate(CO2 = oce::fillGap(CO2))
  sfn_aggr_day <- calc_sfn_aggr_E(sfn_agr,"1 day")
  env_day <- calc_sfn_aggr_env(sfn,"1 day")
  env_day <- env_day %>% mutate(CO2 = oce::fillGap(CO2))
  
  
  
  #### PMODEL WEEKLY ####
  #DATA preparation
  df <- sfn_aggr_week %>%
    left_join(env_week, by = c("timestamp_aggr", "TIMESTAMP")) %>% 
    bind_cols(tibble(opt_swc_slope = opt_swc[1],opt_swc_int =  opt_swc[2])) %>% 
    # filter(!is.na(E_stand)) %>% 
    # filter(lubridate::date(TIMESTAMP) >= lubridate::ymd(20140101), lubridate::date(TIMESTAMP) < lubridate::ymd(20150101)) %>%
    mutate(PPFD = PPFD*1e6/86400,
           ppfd_ERA5 = sw_ERA5 * 2.114, #transform sw ERA5 to ppfd
           netr = netr*1e6/86400,
           patm = rpmodel::calc_patm(sfn$si_elev %>% unique()),
           swvl = swvl_aggregation(st_soil_depth,swvl1,swvl2,swvl3,swvl4),
           # swvl = swc_shallow,
           swvl =  opt_swc_slope*swvl+opt_swc_int,
           REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE)),
           vpd = case_when(is.na(vpd)~ VPD/1000,
                           !is.na(vpd)~ vpd),
           ta = case_when(is.na(ta)~ tc,
                          !is.na(ta)~ ta),
           ppfd_in = case_when(is.na(ppfd_in)~ PPFD, #if there is ppfd in SFN, use if not use PPFD from EC and if not use ppfd calculated from SW from ERA5
                               is.na(PPFD)~ ppfd_ERA5,
                              !is.na(ppfd_in)~ ppfd_in),
           netrad = case_when(is.na(netrad)~ netr,
                              !is.na(netrad)~ netrad),
           CO2 = case_when(is.na(CO2) ~ 400,
                              !is.na(CO2) ~ CO2)) %>%
    filter(!is.na(REW))
  
  # fit_ws <- lm(ws~0+ws_ERA, data =df) %>% summary()
  if(any(names(df) == "ws")){
  df <- df %>% mutate(ws = case_when(is.na(ws)~2,
                                     !is.na(ws)~ws))
  }else{df <- df %>% mutate(ws =2)}

  #MODEL CALCULATION
  pm_weekly <- df %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, 
                       ppfd = x$ppfd_in/1e6, soilm =x$REW, do_soilmstress = FALSE, elv=sfn$si_elev %>% unique() )
      return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
    })
  colnames(pm_weekly)[-1] <- paste(colnames(pm_weekly[,-1]), "pmodel", sep = "_")
  
  pm_weekly_smith <- df %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, ppfd = x$ppfd_in/1e6, soilm =x$REW,
                       do_soilmstress = FALSE, elv=sfn$si_elev %>% unique(), method_jmaxlim = "smith19" )
      return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
    })
  colnames(pm_weekly_smith)[-1] <- paste(colnames(pm_weekly_smith[,-1]), "pmodel_smith", sep = "_")

  df_res <- df %>% 
    left_join(pm_weekly, by = "timestamp_aggr") %>% 
    left_join(pm_weekly_smith, by = "timestamp_aggr") %>% 
    mutate(`E Sap flow`= E_stand/18.2, #mol m-2 h-1
           PMODEL = 1.6*gs_pmodel*(vpd*1000)*3600,
           `PMODEL Smith` = 1.6*gs_pmodel_smith*(vpd*1000)*3600,#mol m-2 h-1
    ) %>% 
    bind_cols(si_code = x)
    
  # if there is AET/PET then calculate PMODEL with SWC limitation
  if(!is.na(meanalpha)){
    pm_weekly_swc <- df %>% 
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, ppfd = x$ppfd_in/1e6, soilm =x$REW,
                                do_soilmstress = TRUE, elv=sfn$si_elev %>% unique(),meanalpha = meanalpha)
        return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
      })
    colnames(pm_weekly_swc)[-1] <- paste(colnames(pm_weekly_swc[,-1]), "pmodel_swc", sep = "_")
    df_res <- df_res %>% 
      left_join(pm_weekly_swc, by = "timestamp_aggr") %>%  
      mutate(`PMODEL swc limitation`= 1.6*gs_pmodel_swc*(vpd*1000)*3600, #mol m-2 h-1
             aet = aet * 55.5/24
      )
  }
  
if(PHYDRO_TRUE & opt_swc[1] != 1 & opt_swc[2] != 0){
  df %>%  
    filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),!is.na(netrad), !is.na(vpd),FAPAR > 0, ppfd_in > 0) %>%
    split(seq(nrow(.)))%>%
    purrr::map(function(x){
      print(x$timestamp_aggr)
      psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, sand = x$st_sand_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "VG")
      if(is.na(psi_soil)){
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, sand = x$st_sand_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
        }
      res <- pmodel_hydraulics_numerical(tc = x$ta, ppfd =x$ppfd_in/x$LAI, vpd = x$vpd, u = x$ws, ustar=NA, nR=x$netrad, co2=x$CO2, LAI = x$LAI,
                                         elv = x$si_elev, fapar = (x$FAPAR), 
                                         kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, par_cost = NULL, 
                                         opt_hypothesis = "PM", gs_approximation = "Ohm")
      
      return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
    })%>% bind_rows()->df_Ohm
  colnames(df_Ohm)[-1] <- paste(colnames(df_Ohm[,-1]), "Ohm", sep = "_")
  
  df %>%  
    filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),!is.na(netrad), !is.na(vpd),FAPAR > 0, ppfd_in > 0) %>%
    split(seq(nrow(.)))%>%
    purrr::map(function(x){
      print(paste("PM:", x$timestamp_aggr))
      psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, sand = x$st_sand_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "VG")
      if(is.na(psi_soil)){
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, sand = x$st_sand_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
      }
      res <- pmodel_hydraulics_numerical(tc = x$ta, ppfd =x$ppfd_in/x$LAI, vpd = x$vpd, u = x$ws, ustar=NA, nR=x$netrad, co2=x$CO2, LAI = x$LAI,
                                         elv = x$si_elev, fapar = (x$FAPAR), 
                                         kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, 
                                         par_cost = list(
                                           alpha = 0.1,       # cost of Jmax
                                           gamma = 1          # cost of hydraulic repair
                                         ), 
                                         opt_hypothesis = "PM", gs_approximation = "PM")
      return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
    })%>% bind_rows()->df_PM
  colnames(df_PM)[-1] <- paste(colnames(df_PM[,-1]), "PM", sep = "_")

  
  df_res <- df_res %>% 
    left_join(df_PM, by = "timestamp_aggr") %>%
    left_join(df_Ohm, by = "timestamp_aggr") %>%
    mutate(PHydro = E_Ohm/1e6*3600*LAI , #transform to mol m-2soil h-1
           `PHydro Penman-Monteith` = E_PM*3600*LAI#
           # `PHydro Penman-Monteith` = E_PM*55.5*3600 *LAI#,
           # `PHydro Penman-Monteith simple` = E_simple_PM/1e6*3600*LAI
    )
}
  #### SAVE DATA ####

  save(df_res, file = paste0("DATA/OUTCOME_WEEKLY/",x,".RData"))
  
  
  #### PMODEL MONTHLY ####
  df <- sfn_aggr_month %>%
    left_join(env_month, by = c("timestamp_aggr", "TIMESTAMP")) %>% 
    bind_cols(tibble(opt_swc_slope = opt_swc[1],opt_swc_int =  opt_swc[2])) %>% 
    # filter(!is.na(E_stand)) %>% 
    # filter(lubridate::date(TIMESTAMP) >= lubridate::ymd(20140101), lubridate::date(TIMESTAMP) < lubridate::ymd(20150101)) %>%
    mutate(PPFD = PPFD*1e6/86400,
           ppfd_ERA5 = sw_ERA5 * 2.114, #transform sw ERA5 to ppfd
           netr = netr*1e6/86400,
           patm = rpmodel::calc_patm(sfn$si_elev %>% unique()),
           swvl = swvl_aggregation(st_soil_depth,swvl1,swvl2,swvl3,swvl4),
           # swvl = swc_shallow,
           swvl =  opt_swc_slope*swvl+opt_swc_int,
           REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE)),
           vpd = case_when(is.na(vpd)~ VPD/1000,
                           !is.na(vpd)~ vpd),
           ta = case_when(is.na(ta)~ tc,
                          !is.na(ta)~ ta),
           ppfd_in = case_when(is.na(ppfd_in)~ PPFD, #if there is ppfd in SFN, use if not use PPFD from EC and if not use ppfd calculated from SW from ERA5
                               is.na(PPFD)~ ppfd_ERA5,
                               !is.na(ppfd_in)~ ppfd_in),
           netrad = case_when(is.na(netrad)~ netr,
                              !is.na(netrad)~ netrad),
           CO2 = case_when(is.na(CO2) ~ 400,
                           !is.na(CO2) ~ CO2)) %>%
    filter(!is.na(REW))
  # fit_ws <- lm(ws~0+ws_ERA, data =df) %>% summary()
  if(any(names(df) == "ws")){
    df <- df %>% mutate(ws = case_when(is.na(ws)~2,
                                       !is.na(ws)~ws))
  }else{df <- df %>% mutate(ws =2)}
  
  
  pm_monthly <- df %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, 
                              ppfd = x$ppfd_in/1e6, soilm =x$REW, do_soilmstress = FALSE, elv=sfn$si_elev %>% unique() )
      return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
    })
  colnames(pm_monthly)[-1] <- paste(colnames(pm_monthly[,-1]), "pmodel", sep = "_")

  pm_monthly_smith <- df %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, ppfd = x$ppfd_in/1e6, soilm =x$REW,
                              do_soilmstress = FALSE, elv=sfn$si_elev %>% unique(), method_jmaxlim = "smith19" )
      return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
    })
  colnames(pm_monthly_smith)[-1] <- paste(colnames(pm_monthly_smith[,-1]), "pmodel_smith", sep = "_")
  
  df_res <- df %>% 
    left_join(pm_monthly, by = "timestamp_aggr") %>% 
    left_join(pm_monthly_smith, by = "timestamp_aggr") %>% 
    mutate(`E Sap flow`= E_stand/18.2, #mol m-2 h-1
           PMODEL = 1.6*gs_pmodel*(vpd*1000)*3600,
           `PMODEL Smith` = 1.6*gs_pmodel_smith*(vpd*1000)*3600,#mol m-2 h-1
          ) %>% 
    bind_cols(si_code = x)
  
  # if there is AET/PET then calculate PMODEL with SWC limitation
  if(!is.na(meanalpha)){
    pm_monthly_swc <- df %>% 
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, ppfd = x$ppfd_in/1e6, soilm =x$REW,
                                do_soilmstress = TRUE, elv=sfn$si_elev %>% unique(),meanalpha = meanalpha)
        return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
      })
    colnames(pm_monthly_swc)[-1] <- paste(colnames(pm_monthly_swc[,-1]), "pmodel_swc", sep = "_")
    df_res <- df_res %>% 
      left_join(pm_monthly_swc, by = "timestamp_aggr") %>%  
      mutate(`PMODEL swc limitation`= 1.6*gs_pmodel_swc*(vpd*1000)*3600, #mol m-2 h-1
             aet = aet * 55.5/24
      )
  }  
  
  if(PHYDRO_TRUE & opt_swc[1] != 1 & opt_swc[2] != 0){
    df %>%  
      filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),!is.na(netrad), !is.na(vpd), FAPAR > 0, ppfd_in > 0) %>%
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        print(x$timestamp_aggr)
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, sand = x$st_sand_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "VG")
        if(is.na(psi_soil)){
          psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, sand = x$st_sand_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
        }
        res <- pmodel_hydraulics_numerical(tc = x$ta, ppfd =x$ppfd_in/x$LAI, vpd = x$vpd, u = x$ws, ustar=NA, nR=x$netrad, co2=x$CO2, LAI = x$LAI,
                                           elv = x$si_elev, fapar = (x$FAPAR), 
                                           kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, par_cost = NULL, 
                                           opt_hypothesis = "PM", gs_approximation = "Ohm")
        
        return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
      })%>% bind_rows()->df_Ohm
    colnames(df_Ohm)[-1] <- paste(colnames(df_Ohm[,-1]), "Ohm", sep = "_")
    
    df %>%  
      filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),!is.na(netrad), !is.na(vpd), ppfd_in>0, FAPAR>0) %>%
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        print(paste("PM:", x$timestamp_aggr))
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, sand = x$st_sand_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "VG")
        if(is.na(psi_soil)){
          psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, sand = x$st_sand_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
        }
        res <- pmodel_hydraulics_numerical(tc = x$ta, ppfd =x$ppfd_in/x$LAI, vpd = x$vpd, u = x$ws, ustar=NA, nR=x$netrad, co2=x$CO2, LAI = x$LAI,
                                           elv = x$si_elev, fapar = (x$FAPAR), 
                                           kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, 
                                           par_cost = list(
                                             alpha = 0.1,       # cost of Jmax
                                             gamma = 1          # cost of hydraulic repair
                                           ), 
                                           opt_hypothesis = "PM", gs_approximation = "PM")
        return(tibble(timestamp_aggr = x$timestamp_aggr) %>% bind_cols(res))
      })%>% bind_rows()->df_PM
    colnames(df_PM)[-1] <- paste(colnames(df_PM[,-1]), "PM", sep = "_")
    
    
    df_res <- df_res %>% 
      left_join(df_PM, by = "timestamp_aggr") %>%
      left_join(df_Ohm, by = "timestamp_aggr") %>%
      mutate(PHydro = E_Ohm/1e6*3600*LAI , #transform to mol m-2soil h-1
             `PHydro Penman-Monteith` = E_PM*3600*LAI#
             # `PHydro Penman-Monteith` = E_PM*55.5*3600 *LAI#,
             # `PHydro Penman-Monteith simple` = E_simple_PM/1e6*3600*LAI
      )
  }
  #### REMOVE DATA ####
  rm(sfn)
  gc()
  
  save(df_res, file = paste0("DATA/OUTCOME_MONTHLY/",x,".RData"))
  #return(df_res)
  
})#->DF


value <- 1
DF[[value]] %>% 
  # filter(!is.na(PMODEL)) %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
  )) %>%
  ggplot(aes(x=name, y = (value-`E Sap flow`),fill =name))+
  ggbeeswarm::geom_quasirandom(shape =21,size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  coord_flip()+
  ghibli::scale_fill_ghibli_d("MononokeMedium",direction = -1)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow])))+
  xlab("MODEL")+
  theme_classic()-> P1

DF[[value]] %>% 
  # filter(!is.na(PMODEL)) %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
  )) %>% 
  as_tibble() %>%
  dplyr::select(name,value,`E Sap flow`) %>% 
  na.omit() %>% 
  group_by(name) %>% 
  summarise(cor_value = cor(`E Sap flow`,value)) %>% 
  ggplot(aes(x=name, y = cor_value, fill =name))+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  coord_flip()+
  ghibli::scale_fill_ghibli_d("MononokeMedium",direction = -1)+
  ylab("r Pearson's Correlation")+
  xlab("MODEL")+
  theme_classic()->P2

DF[[value]] %>% 
  # filter(TIMESTAMP > as.Date(lubridate::ymd("2012/01/01"))) %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
  )) %>% 
  ggplot()+
  geom_line(aes(timestamp_aggr,value,color=name),size = 0.7,alpha=0.9, show.legend = FALSE)+
  geom_line(aes(timestamp_aggr,`E Sap flow`),color="grey10", size = 1)+
  ghibli::scale_color_ghibli_d("MononokeMedium",direction = -1)+
  ylab(expression(paste(mol[H2O]," ", m[soil]^{-2}," ", h^{-1})))+
  theme_classic()->P3

gridExtra::grid.arrange(P3,gridExtra::arrangeGrob(P1,P2,ncol = 2),ncol=1)
