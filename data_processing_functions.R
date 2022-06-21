###################################################
# DATA PROCESSING FUNCTIONS
###################################################
#Created on 12 MAY 2022
#author: vflo

add_miss_var <- function(sfn=sfn){
  if(!any(colnames(sfn) == "ws")){sfn <- sfn %>% mutate(ws = c(NA))}
  if(!any(colnames(sfn) == "rh")){sfn <- sfn %>% mutate(rh = c(NA))}
  if(!any(colnames(sfn) == "netrad")){sfn <- sfn %>% mutate(netrad = c(NA))}
  if(!any(colnames(sfn) == "sw_in")){sfn <- sfn %>% mutate(sw_in = c(NA))}
  if(!any(colnames(sfn) == "ppfd_in")){sfn <- sfn %>% mutate(ppfd_in = c(NA))}
  if(!any(colnames(sfn) == "vpd")){sfn <- sfn %>% mutate(vpd = c(NA))}
  if(!any(colnames(sfn) == "ta")){sfn <- sfn %>% mutate(ta = c(NA))}
  if(!any(colnames(sfn) == "swc_shallow")){sfn <- sfn %>% mutate(swc_shallow = c(NA))}
  if(!any(colnames(sfn) == "st_soil_depth")){sfn <- sfn %>% mutate(st_soil_depth = c(NA))}
  if(!any(colnames(sfn) == "CO2")){sfn <- sfn %>% mutate(CO2 = 400)}
  if(!any(colnames(sfn) == "aet")){sfn <- sfn %>% mutate(aet = NA)}
  if(!any(colnames(sfn) == "cond")){sfn <- sfn %>% mutate(cond = NA)}
  if(!any(colnames(sfn) == "netr")){sfn <- sfn %>% mutate(netr = NA)}
  if(!any(colnames(sfn) == "PPFD")){sfn <- sfn %>% mutate(PPFD = NA)}
  if(!any(colnames(sfn) == "tc")){sfn <- sfn %>% mutate(tc = c(NA))}
  if(!any(colnames(sfn) == "VPD")){sfn <- sfn %>% mutate(VPD = c(NA))}
  sfn
}


include_fapar_lai <- function(sfn, fapar_noaa, fapar){
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

  fapar_filled <- fapar_noaa_filled %>% bind_rows(fapar_filled) %>% 
    filter(TIMESTAMP>= min(sfn$TIMESTAMP, na.rm = TRUE),TIMESTAMP< max(sfn$TIMESTAMP, na.rm = TRUE))

  min_data = min(fapar_filled$TIMESTAMP)

# Build FAPAR model
  k_ <- as.integer(nrow(fapar_filled)*0.1)
  if(k_ >= 3){
    
    model_fapar <- qgam(Fapar ~ s(TIMESTAMP, bs = 'tp', k = k_), data = fapar_filled %>%
                        dplyr::select(TIMESTAMP,Fapar) %>%
                        mutate(TIMESTAMP = TIMESTAMP %>%
                                 as.numeric()), qu = 0.6)
    model_lai <- qgam(Lai ~ s(TIMESTAMP, bs = 'tp', k = k_), data = fapar_filled %>%
                      dplyr::select(TIMESTAMP,Lai) %>%
                      mutate(TIMESTAMP = TIMESTAMP %>%
                               as.numeric()), qu = 0.6)
  # x <- fapar_filled$TIMESTAMP %>% as.numeric()
  # y_fapar <- fapar_filled$Fapar
  # y_lai <- fapar_filled$Lai
  # model_fapar <- GauPro(x, y_fapar, parallel=TRUE)
  # model_lai <- GauPro(x, y_lai, parallel=TRUE)
  # 
  # pred_FAPAR <- model_fapar$pred(sfn[,"TIMESTAMP"] %>% 
  #                         mutate(TIMESTAMP = as.numeric(lubridate::ymd(TIMESTAMP)))%>% 
  #                         unique() %>% 
  #                         pull())
  # pred_LAI <- model_lai$pred(sfn[,"TIMESTAMP"] %>%
  #                              mutate(TIMESTAMP = as.numeric(lubridate::ymd(TIMESTAMP)))%>% 
  #                              unique() %>% 
  #                              pull())
    pred_FAPAR <- predict(model_fapar, tibble(TIMESTAMP = sfn[,"TIMESTAMP"] %>%
                                                mutate(TIMESTAMP = as.numeric(lubridate::ymd(TIMESTAMP)))%>% 
                                                unique() %>% pull()))
    pred_LAI <- predict(model_lai, tibble(TIMESTAMP = sfn[,"TIMESTAMP"] %>%
                                          mutate(TIMESTAMP = as.numeric(lubridate::ymd(TIMESTAMP)))%>% 
                                          unique() %>% pull()))
  
    sfn %>% 
      left_join(tibble(TIMESTAMP = sfn$TIMESTAMP %>% 
                       unique(), 
                       FAPAR = pred_FAPAR) %>% 
                  filter(TIMESTAMP > min_data)%>%
                  mutate(FAPAR = case_when(FAPAR<0~0,
                                           FAPAR>=0~FAPAR)),
                by = "TIMESTAMP") %>% 
      left_join(tibble(TIMESTAMP = sfn$TIMESTAMP %>% 
                       unique(),  LAI = pred_LAI) %>%
                  filter(TIMESTAMP > min_data) %>% 
                  mutate(LAI = case_when(LAI<0~0,
                                        LAI>=0~LAI)),
                by = "TIMESTAMP") %>% 
    left_join(fapar_filled, by = "TIMESTAMP")
    
  }else{
    sfn %>%
      left_join(fapar_filled, by = "TIMESTAMP") %>% 
      mutate(FAPAR = mean(Fapar,na.rm = TRUE),
             LAI = mean(Lai,na.rm = TRUE))
  }
}


sp_params <- function(sfn, SpParams){
  density <- rpmodel::calc_density_h2o(25, rpmodel::calc_patm(sfn$si_elev %>% unique()))
  visco <- rpmodel::calc_viscosity_h2o(25, rpmodel::calc_patm(sfn$si_elev %>% unique()))
  sp_data <- SpParams %>% filter(pl_species %in% (sfn$pl_species %>% unique()))
  par_plant_std <- NULL
  
  #If there is species data -> create par_plant_list to calculate PHYDRO
  if(nrow(sp_data)>0){
    sp_data %>% 
      split(seq(nrow(.))) %>%
      purrr::map(do_phydro_params, sfn = sfn, visco = visco, density = density) %>% bind_rows() -> par_plant_std
  
    par_plant_std <- par_plant_std %>% 
      summarise(st_basal_area = weighted.mean(st_basal_area, sp_basal_area),
                pl_sapw_area = weighted.mean(pl_sapw_area, sp_basal_area, na.rm = TRUE),
                pl_ba = weighted.mean(plant_ba, sp_basal_area, na.rm = TRUE),
                # v_huber = weighted.mean(v_huber, sp_basal_area),
                height = weighted.mean(height, sp_basal_area),
                LAI = weighted.mean(LAI, sp_basal_area),
                conductivity = weighted.mean(conductivity, sp_basal_area), #Kg m-1 s-1 MPa
                conductivity_Kl = weighted.mean(conductivity_Kl, sp_basal_area),
                psi50 = weighted.mean(psi50, sp_basal_area),
                b = weighted.mean(b, sp_basal_area),
                c = weighted.mean(c, sp_basal_area),
                d = weighted.mean(d, sp_basal_area),
                v25 = weighted.mean(v25, sp_basal_area),
                j25 = weighted.mean(j25, sp_basal_area),
                sp_basal_area = sum(sp_basal_area)) %>% 
      mutate(K = conductivity/height*55.5*pl_sapw_area/pl_ba*st_basal_area*1e-4)#mol m-2(ground) s-1 MPa-1
  
    # if(is.na(unique(par_plant_std$height))){par_plant_std$height <- unique(sfn$st_height)}
  
  }
  
  return(par_plant_std)
  
}



mean_alpha_calc <- function(sfn, sw_in_monthly_average, temp_monthly_average, precip_monthly_average, env_month, soil){
  sw_in_month <- sw_in_monthly_average %>% 
    filter(si_code == sfn$si_code %>% unique) %>% 
    as_tsibble(index = timestamp) %>% 
    dplyr::select(sw_in_month) %>%
    index_by(timestamp_aggr = ~ yearmonth(.)) %>% 
    as_tibble() %>% 
    dplyr::select(-timestamp, sw_in_month) %>% 
    as_tsibble(index = timestamp_aggr) %>% 
    as.ts() %>%  
    as.xts()

  temp_month <- temp_monthly_average %>% 
    filter(si_code == sfn$si_code %>% unique) %>% 
    as_tsibble(index = timestamp) %>% 
    dplyr::select(temp_month) %>%
    index_by(timestamp_aggr = ~ yearmonth(.)) %>% 
    as_tibble() %>% 
    dplyr::select(-timestamp, temp_month) %>% 
    as_tsibble(index = timestamp_aggr) %>% 
    as.ts() %>%  
    as.xts()

  precip_month <- precip_monthly_average %>% 
    filter(si_code == sfn$si_code %>% unique) %>% 
    as_tsibble(index = timestamp) %>%
    dplyr::select(precip_month) %>%
    index_by(timestamp_aggr = ~ yearmonth(.)) %>% 
    as_tibble() %>% 
    dplyr::select(-timestamp, precip_month) %>% 
    as_tsibble(index = timestamp_aggr) %>% 
    as.ts() %>%  
    as.xts()

  si_elev <- unique(env_month$si_elev)
  lat <- unique(sfn$si_lat)

  timeind<-as.Date(format(time(temp_month),'%Y-%m-%d'))

  time(sw_in_month)<-timeind;
  time(temp_month)<-timeind;
  time(precip_month)<-timeind

  # soil <- soil %>% mutate(depth = 7)
  
  splash_out <- rsplash::splash.point(sw_in_month,
                                    temp_month,
                                    precip_month,
                                    lat,
                                    si_elev,
                                    0,
                                    0,
                                    soil_data = as.numeric(soil))
  
  # rsplash::unSWC(as.numeric(soil),7,splash_out$wn,units="v/v")
  
  mean(splash_out$aet,na.rm = TRUE)/mean(splash_out$pet,na.rm = TRUE)
}


data_prep <- function(df_sfn, env, opt_swc){
  df_sfn %>% 
  left_join(env, by = c("timestamp_aggr", "TIMESTAMP")) %>%
  bind_cols(tibble(opt_swc_slope = opt_swc[1],opt_swc_int =  opt_swc[2])) %>% 
  # filter(!is.na(E_stand)) %>% 
  # filter(lubridate::date(TIMESTAMP) >= lubridate::ymd(20140101), lubridate::date(TIMESTAMP) < lubridate::ymd(20150101)) %>%
  mutate(E_sapflow = E_stand/18.2, #mol m-2 h-1
         E_sapflow_sd = E_stand_sd/18.2, #mol m-2 h-1
         PPFD = PPFD*1e6/86400,
         ppfd_ERA5 = sw_ERA5 * 2.114, #transform sw ERA5 to ppfd
         netr = netr*1e6/86400,
         patm = rpmodel::calc_patm(si_elev %>% unique()),
         swvl = swvl_aggregation(st_soil_depth,swvl1,swvl2,swvl3,swvl4),
         REW = swvl/max_swvl,
         # swvl = swc_shallow,
         swvl =  opt_swc_slope*swvl+opt_swc_int,
         # REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE)),
         vpd = case_when(is.na(vpd)~ VPD/1000,
                         !is.na(vpd)~ vpd),
         ta = case_when(is.na(ta)~ tc,
                        !is.na(ta)~ ta),
         ppfd_in = case_when(is.na(ppfd_in) & !is.na(PPFD) ~ PPFD, #if there is ppfd in SFN, use if not use PPFD from EC and if not use ppfd calculated from SW from ERA5
                             is.na(ppfd_in) & is.na(PPFD) ~ ppfd_ERA5,
                             !is.na(ppfd_in) ~ ppfd_in),
         netrad = case_when(is.na(netrad)~ netr,
                            !is.na(netrad)~ netrad),
         CO2 = case_when(is.na(CO2) ~ 400,
                         !is.na(CO2) ~ CO2),
         patm = rpmodel::calc_patm(si_elev),
         Gs_sapflow = E_sapflow / (vpd*1000) / 1.6 / 3600*1e6, #umol m-2(ground) s-1 Pa-1
         Gs_sapflow_sd = E_sapflow_sd / (vpd*1000) / 1.6 / 3600*1e6) %>% 
  filter(!is.na(REW))->df

# fit_ws <- lm(ws~0+ws_ERA, data =df) %>% summary()
  if(any(names(df) == "ws")){
    df <- df %>% mutate(ws = case_when(is.na(ws)~2,
                                     !is.na(ws)~ws))
  }else{df <- df %>% mutate(ws =2)}
  
  return(df)
}


data_prep_hourly <- function(df_sfn, env, opt_swc){
  df_sfn %>% 
    left_join(env) %>%
    bind_cols(tibble(opt_swc_slope = opt_swc[1],opt_swc_int =  opt_swc[2])) %>% 
    # filter(!is.na(E_stand)) %>% 
    # filter(lubridate::date(TIMESTAMP) >= lubridate::ymd(20140101), lubridate::date(TIMESTAMP) < lubridate::ymd(20150101)) %>%
    mutate(E_sapflow = E_stand/18.2, #mol m-2 h-1
           E_sapflow_sd = E_stand_sd/18.2, #mol m-2 h-1
           PPFD = PPFD*1e6/86400,
           ppfd_ERA5 = sw_ERA5 * 2.114, #transform sw ERA5 to ppfd
           netr = netr*1e6/86400,
           patm = rpmodel::calc_patm(si_elev %>% unique()),
           swvl = swvl_aggregation(st_soil_depth,swvl1,swvl2,swvl3,swvl4),
           REW = swvl/max_swvl,
           # swvl = swc_shallow,
           swvl =  opt_swc_slope*swvl+opt_swc_int,
           # REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE)),
           vpd = case_when(is.na(vpd)~ VPD/1000,
                           !is.na(vpd)~ vpd),
           ta = case_when(is.na(ta)~ tc,
                          !is.na(ta)~ ta),
           ppfd_in = case_when(is.na(ppfd_in) & !is.na(PPFD) ~ PPFD, #if there is ppfd in SFN, use if not use PPFD from EC and if not use ppfd calculated from SW from ERA5
                               is.na(ppfd_in) & is.na(PPFD) ~ ppfd_ERA5,
                               !is.na(ppfd_in) ~ ppfd_in),
           netrad = case_when(is.na(netrad)~ netr,
                              !is.na(netrad)~ netrad),
           CO2 = case_when(is.na(CO2) ~ 400,
                           !is.na(CO2) ~ CO2),
           patm = rpmodel::calc_patm(si_elev),
           Gs_sapflow = E_sapflow / (vpd*1000) / 1.6 / 3600*1e6, #umol m-2(ground) s-1 Pa-1
           Gs_sapflow_sd = E_sapflow_sd / (vpd*1000) / 1.6 / 3600*1e6) %>% 
    filter(!is.na(REW))->df
  
  # fit_ws <- lm(ws~0+ws_ERA, data =df) %>% summary()
  if(any(names(df) == "ws")){
    df <- df %>% mutate(ws = case_when(is.na(ws)~2,
                                       !is.na(ws)~ws))
  }else{df <- df %>% mutate(ws =2)}
  
  return(df)
}



LIST <- function(...) {
  nms <- sapply(as.list(substitute(list(...))), deparse)[-1]
  setNames(list(...), nms)
}
