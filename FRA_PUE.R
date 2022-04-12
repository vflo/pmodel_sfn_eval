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


sapply(list('calc_sfn_aggr_E.R','calc_sfn_aggr_env.R','calc_sfn_scale.R','swvl_aggregation.R','calc_optim_swc.R'), source, .GlobalEnv)

#### LOAD SAPFLUXNET DATA ####
# path <- "../0.1.5/RData/plant"
# sfn_meta <- sapfluxnetr::read_sfn_metadata(folder = path, .write_cache = TRUE)
# save(sfn_meta, file = "DATA/sfn_meta.RData")
# load(file = "DATA/sfn_meta.RData")
# sfn_meta$site_md %>% filter(si_code %in% c('AUS_WOM', 'FRA_FON', 'FRA_PUE', 'NLD_LOO', 'RUS_FYO', 'USA_UMB_CON')) %>% dplyr::select(si_code,si_lat,si_long) %>% write_csv(file='sfn_sites_coord.csv')
# 
# 
# FRA_PUE <- read_sfn_data("FRA_PUE", folder = path) %>%
#   sfn_metrics(period = "1 hour", .funs = list(~ mean(., na.rm = TRUE)),
#               solar = TRUE, interval = "general") %>%
#   sapfluxnetr::metrics_tidyfier(metadata = sfn_meta,interval = "general")
# 
# save(FRA_PUE, file = "DATA/FRA_PUE.RData")
# load(file = "DATA/FRA_PUE.RData")
# 
# FRA_PUE_2010 <- FRA_PUE %>% filter(lubridate::date(TIMESTAMP) >= lubridate::ymd(20100101))
# rm(FRA_PUE)
# gc()
# save(FRA_PUE_2010, file = "DATA/FRA_PUE_2010.RData")
load(file = "DATA/FRA_PUE_2014.RData")
FRA_PUE <- FRA_PUE_2014
rm(FRA_PUE_2014)
gc()


Fr_Pue <-splashTools::readFluxdata("DATA/FLX_FR-Pue_FLUXNET2015_FULLSET_HH_2000-2014_2-4.csv",FRA_PUE$si_elev %>% unique())

meanalpha <- sum(Fr_Pue$aet,na.rm = TRUE)/sum(Fr_Pue$pet,na.rm = TRUE)


FRA_PUE <- FRA_PUE %>% left_join(Fr_Pue %>% fortify.zoo %>% as_tibble %>% rename(TIMESTAMP = Index))
#### SITE SWC ####
# WGScoor <-  FRA_PUE %>% select(long = si_long, lat = si_lat)
# coordinates(WGScoor)=~long+lat
# proj4string(WGScoor)<- CRS("+proj=longlat +datum=WGS84")
# LLcoor<-spTransform(WGScoor,CRS("+proj=longlat"))
# raster::shapefile(LLcoor, "FRA_PUE.shp")

load(file = "DATA/df_hourly.RData")
df_swc <- df %>% 
  pivot_wider(names_from = variable, values_from = value, values_fn = mean) %>% 
  rename(TIMESTAMP = timestamp)# %>% 
  # mutate(ws_ERA = sqrt(v10^2 + u10^2))
rm(df)
gc()

FRA_PUE <- FRA_PUE %>% left_join(df_swc %>% filter(abs(si_lat - 43.74) == min(abs(si_lat - 43.74))), by=c("TIMESTAMP"#,'si_lat','si_long'
                                              ))

##CALIBRATE SOIL POTENTIAL by adding SWC
psi_FRA_PUE <-  read_csv("DATA/psi_FRA_PUE.csv")

## THE OPTIMUM SWC value to add
opt_swc <- optim_swc(FRA_PUE, psi_FRA_PUE, type ="swc")

#### FAPAR ####
# fapar <-read_csv("DATA/FAR_PUE_FAPAR_RAW_2015_MODIS.csv")
fapar <-read_csv("DATA/FRA_PUE_FAPAR.csv")
fapar_filled <-df_swc %>% 
              filter(lubridate::date(TIMESTAMP) >= lubridate::ymd(20140101)) %>% 
              dplyr::select(TIMESTAMP) %>% 
              mutate(TIMESTAMP = date(TIMESTAMP)) %>% 
              unique() %>% 
  left_join(fapar %>% 
              dplyr::select(TIMESTAMP = date, Fapar = Fpar,Lai = Lai) %>%
              mutate(Fapar = Fapar*0.01,
                     Lai = Lai*0.1)) %>%
  mutate(Fapar = oce::fillGap(Fapar),
         Lai = oce::fillGap(Lai)) %>% 
  filter(!is.na(Fapar),!is.na(Lai))
FAPAR <- smooth.spline(fapar_filled$TIMESTAMP, fapar_filled$Fapar, df = 8)
LAI <- smooth.spline(fapar_filled$TIMESTAMP, fapar_filled$Lai, df = 8)

FRA_PUE <- FRA_PUE %>% left_join(fapar_filled%>% cbind(FAPAR = broom::augment(FAPAR)$.fitted) %>% cbind(LAI = broom::augment(LAI)$.fitted)) 

FRA_PUE <- FRA_PUE %>% 
  filter(lubridate::date(TIMESTAMP) <= lubridate::ymd(20150101))

#### SPATIAL AND TEMPORAL DATA AGGREGATION ####
sfn_hour <- calc_sfn_scale(FRA_PUE)
env_hour <- calc_sfn_aggr_env(FRA_PUE,"1 hour")
env_hour <- env_hour %>% mutate(CO2 = oce::fillGap(CO2))
sfn_aggr_week <- calc_sfn_aggr_E(sfn_hour,"1 week")
env_week <- calc_sfn_aggr_env(FRA_PUE,"1 week")
env_week <- env_week %>% mutate(CO2 = oce::fillGap(CO2))
sfn_aggr_month <- calc_sfn_aggr_E(sfn_hour,"1 month")
env_month <- calc_sfn_aggr_env(FRA_PUE,"1 month")
env_month <- env_month %>% mutate(CO2 = oce::fillGap(CO2))
sfn_aggr_day <- calc_sfn_aggr_E(sfn_hour,"1 day")
env_day <- calc_sfn_aggr_env(FRA_PUE,"1 day")
env_day <- env_day %>% mutate(CO2 = oce::fillGap(CO2))


#### PHYDRO PREPARATION ####
density <- rpmodel::calc_density_h2o(25, rpmodel::calc_patm(FRA_PUE$si_elev %>% unique()))
visco <- rpmodel::calc_viscosity_h2o(25, rpmodel::calc_patm(FRA_PUE$si_elev %>% unique()))
par_plant_std = list(
  # Ks0=1e-12, # m2
  v_huber=FRA_PUE$pl_sapw_area[1]/FRA_PUE$pl_leaf_area[1]/10^4, # m2sapwood m-2leaf
  height=FRA_PUE$pl_height[1], # m
  LAI = FRA_PUE$st_lai[1], # m2leaf m-2soil
  conductivity = 2.63*visco/(1e9*density*55.5), #kmax leaf = 2.63 mmol s-1 m-2 MPa-1; it should be provided as m3/m2/s/Pa: to back transform -> kmax*viscosity_water/(1e3*1e6*density*55.5)
  # conductivity_scalar=1,  # Scales (Ks0*v_huber/H)
  psi50 = -4.2, # MPa
  b=2
)


#### PMODEL DAILY ####
df_day <- sfn_aggr_day %>% 
  left_join(env_day) %>% 
  mutate(ppfd = ppfd_in/1e6,
         patm = rpmodel::calc_patm(FRA_PUE$si_elev %>% unique()),
         swvl = swvl_aggregation(st_soil_depth,swvl1,swvl2,swvl3,swvl4),
         swvl = opt_swc[1]*swvl+opt_swc[2],
         # swvl = swc_shallow,
         REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE))) %>% 
  filter(!is.na(REW))
fit_ws <- lm(ws~0+ws_ERA, data =df_day) %>% summary()
df_day <- df_day %>% mutate(ws_ERA = case_when(is.na(ws)~ws_ERA*coefficients(fit_ws)[,'Estimate'],
                                       !is.na(ws)~ws)) %>% 
  dplyr::select(-ws)

pm_day <- df_day %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=x[["CO2"]], fapar=x[["FAPAR"]], ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = FALSE,meanalpha = 1, elv=FRA_PUE$si_elev %>% unique() )
  })

pm_day_swc <- df_day %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000,co2=x[["CO2"]], fapar=x[["FAPAR"]], ppfd = x[["ppfd"]], soilm =x[["REW"]], 
                     do_soilmstress = TRUE,meanalpha = meanalpha, elv=FRA_PUE$si_elev %>% unique() )
  })

pm_day_smith <- df_day %>% 
  split(seq(nrow(.))) %>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=x[["CO2"]], fapar=x[["FAPAR"]], ppfd = x[["ppfd"]], soilm =x[["REW"]],
                     do_soilmstress = FALSE, elv=FRA_PUE$si_elev %>% unique(), method_jmaxlim = "smith19" )
  })

pm_day_none  <- df_day %>%
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=x[["CO2"]], fapar=x[["FAPAR"]], 
                     ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = FALSE,
                     elv=FRA_PUE$si_elev %>% unique(), method_jmaxlim = "none")
  })

pm_E <- 1.6*pm_day$gs*(df_day$vpd*1e3)*3600 #mol m-2 h-1
pm_E_swc <- 1.6*pm_day_swc$gs*(df_day$vpd*1e3)*3600 #mol m-2 h-1
pm_E_smith <- 1.6*pm_day_smith$gs*(df_day$vpd*1e3)*3600 #mol m-2 h-1
pm_E_none <- 1.6*pm_day_none$gs*(df_day$vpd*1e3)*3600 #mol m-2 h-1


#### PHYDRO  DAILY ####
#Phydro numerical PM
# df_day %>% na.omit() %>%
#   filter(vpd >0) %>%
#   split(seq(nrow(.)))%>%
#   purrr::map(function(x){
#     psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_sand_perc, sand = x$st_clay_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
#     res <- pmodel_hydraulics_numerical(tc = x$ta, ppfd =x$ppfd_in, vpd = x$vpd*1000, u = x$ws_ERA, ustar=NA, nR=x$netrad, co2=400, LAI = x$LAI,
#                                        elv = x$si_elev, fapar = x[["FAPAR"]], kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, par_cost = NULL,
#                                        opt_hypothesis = "PM", gs_approximation = "PM")
#     return(res)
#   })%>% bind_rows()->df_PM
# colnames(df_PM) <- paste(colnames(df_PM), "PM", sep = "_")
# df_PM <- df_PM %>% bind_cols(df_day %>% na.omit() %>% filter(vpd >0) %>% dplyr::select(TIMESTAMP) )
# 
# df_PM_m <- df_PM %>%
#   mutate(jmax_PM_m = pracma::movavg(jmax_PM, n=50, type = 'e'),
#          vcmax_PM_m = pracma::movavg(vcmax_PM, n=50, type = 'e'))
# 
# #Phydro instantaneous PM
# df_day %>% 
#   left_join(df_PM_m %>% dplyr::select(TIMESTAMP,jmax_PM_m,vcmax_PM_m)) %>% 
#   na.omit() %>%
#   filter(vpd >0) %>% 
#   split(seq(nrow(.)))%>%  
#   purrr::map(function(x){
#     psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_sand_perc, sand = x$st_clay_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
#     res <- pmodel_hydraulics_instantaneous(tc = x$ta, ppfd =x$ppfd_in, vpd = x$vpd*1000, u = x$ws_ERA, ustar=NA, nR=x$netrad, co2=400, vcmax = x$vcmax_PM_m, jmax = x$jmax_PM_m, LAI = x$LAI,
#                                        elv = x$si_elev, fapar = x[["FAPAR"]], kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, par_cost = NULL, 
#                                        opt_hypothesis = "PM", gs_approximation = "PM")
#     return(res)
#   })%>% bind_rows()->df_PM
# colnames(df_PM) <- paste(colnames(df_PM), "PM", sep = "_")
# df_PM <- df_PM %>% bind_cols(df_day %>% na.omit() %>% filter(vpd >0) %>% dplyr::select(TIMESTAMP) )

#Phydro numerical Ohm
df_day %>% na.omit() %>%
  filter(vpd >0) %>% 
  split(seq(nrow(.)))%>% 
  purrr::map(function(x){
    # print(x$TIMESTAMP)
    psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_sand_perc, sand = x$st_clay_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
    res <- pmodel_hydraulics_numerical(tc = x$ta, ppfd = (x$ppfd_in/x$LAI), vpd = x$vpd*1000, u = x$ws_ERA, ustar=NA, nR=x$netrad, co2=x[["CO2"]], LAI = x$LAI,
                                       elv = x$si_elev, fapar =(x[["FAPAR"]]/x$LAI), kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, par_cost = NULL, 
                                       opt_hypothesis = "PM", gs_approximation = "Ohm")
    return(res)
  })%>% bind_rows()->df_Ohm
colnames(df_Ohm) <- paste(colnames(df_Ohm), "Ohm", sep = "_")
df_Ohm <- df_Ohm %>% bind_cols(df_day %>% na.omit() %>% filter(vpd >0) %>% dplyr::select(TIMESTAMP) )

df_Ohm_m <- df_Ohm %>%
  mutate(jmax_Ohm_m = pracma::movavg(jmax_Ohm, n=30, type = 'e'),
         vcmax_Ohm_m = pracma::movavg(vcmax_Ohm, n=30, type = 'e'))

#Phydro instantaneous Ohm
df_day %>% 
  left_join(df_Ohm_m %>% dplyr::select(TIMESTAMP,jmax_Ohm_m,vcmax_Ohm_m)) %>% 
  na.omit() %>%
  filter(vpd >0) %>% 
  split(seq(nrow(.)))%>%  
  purrr::map(function(x){
    psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_sand_perc, sand = x$st_clay_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
    res <- pmodel_hydraulics_instantaneous(tc = x$ta, ppfd = (x$ppfd_in/x$LAI) , vpd = x$vpd*1000, u = x$ws_ERA, ustar=NA, nR=x$netrad, co2=x[["CO2"]], vcmax = x$vcmax_Ohm_m, jmax = x$jmax_Ohm_m, LAI = x$LAI,
                                           elv = x$si_elev, fapar = (x[["FAPAR"]]/x$LAI), kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, par_cost = NULL, 
                                           opt_hypothesis = "PM", gs_approximation = "Ohm")
    return(res)
  })%>% bind_rows()->df_Ohm
colnames(df_Ohm) <- paste(colnames(df_Ohm), "Ohm", sep = "_")
df_Ohm <- df_Ohm %>% bind_cols(df_day %>% na.omit() %>% filter(vpd >0) %>% dplyr::select(TIMESTAMP) )


## Plotting results
df_res <- df_day %>% 
  bind_cols(PMODEL = pm_E) %>% 
  bind_cols(`PMODEL swc limitation` = pm_E_swc) %>%  
  bind_cols(`PMODEL Smith` = pm_E_smith) %>%
  bind_cols(`PMODEL none` = pm_E_none) %>%
  left_join(df_Ohm)%>% 
  # left_join(df_PM) %>% 
  mutate(`E Sap flow`= E_stand/18.2, #mol m-2 h-1
         PHydro = E_Ohm/10^6*3600 / LAI#, #transform to mol m-2soil h-1
         # `PHydro Penman-Monteith` = E_PM/10^6*3600 / LAI
         ) #transform to mol m-2soil h-1


df_res %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith','PHydro'#,'PHydro Penman-Monteith'
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

df_res %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith','PHydro'#,'PHydro Penman-Monteith'
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

df_res %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith','PHydro'#,'PHydro Penman-Monteith'
                        )) %>% 
  ggplot()+
  geom_line(aes(timestamp_aggr,value,color=name),size = 0.7,alpha=0.9, show.legend = FALSE)+
  geom_line(aes(timestamp_aggr,`E Sap flow`),color="grey10", size = 1)+
  ghibli::scale_color_ghibli_d("MononokeMedium",direction = -1)+
  ylab(expression(paste(mol[H2O]," ", m[soil]^{-2}," ", h^{-1})))+
  theme_classic()->P3

gridExtra::grid.arrange(P3,gridExtra::arrangeGrob(P1,P2,ncol = 2),ncol=1)


#### PMODEL WEEKLY ####
df <- sfn_aggr_week %>% 
  left_join(env_week) %>% 
  mutate(ppfd = ppfd_in/1e6,
         patm = rpmodel::calc_patm(FRA_PUE$si_elev %>% unique()),
         swvl = swvl_aggregation(st_soil_depth,swvl1,swvl2,swvl3,swvl4),
         swvl = opt_swc[1]*swvl+opt_swc[2],
         # swvl = swc_shallow,
         REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE))) %>% 
  filter(!is.na(REW))
fit_ws <- lm(ws~0+ws_ERA, data =df) %>% summary()
df <- df %>% mutate(ws_ERA = case_when(is.na(ws)~ws_ERA*coefficients(fit_ws)[,'Estimate'],
                                               !is.na(ws)~ws)) %>% 
  dplyr::select(-ws)

pm_weekly <- df %>% 
  split(seq(nrow(.))) %>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=x[["CO2"]], fapar=x[["FAPAR"]], 
                     ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = FALSE, elv=FRA_PUE$si_elev %>% unique() )
  })

pm_weekly_swc <- df %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=x[["CO2"]], fapar=x[["FAPAR"]], ppfd = x[["ppfd"]], soilm =x[["REW"]],
                     do_soilmstress = TRUE, elv=FRA_PUE$si_elev %>% unique(),meanalpha = meanalpha)
  })

pm_weekly_smith <- df %>% 
  split(seq(nrow(.))) %>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=x[["CO2"]], fapar=x[["FAPAR"]], ppfd = x[["ppfd"]], soilm =x[["REW"]],
                     do_soilmstress = FALSE, elv=FRA_PUE$si_elev %>% unique(), method_jmaxlim = "smith19" )
  })

pm_E <- 1.6*pm_weekly$gs*(df$vpd*1e3)*3600 #mol m-2 h-1
pm_E_swc <- 1.6*pm_weekly_swc$gs*(df$vpd*1e3)*3600 #mol m-2 h-1
pm_E_smith <- 1.6*pm_weekly_smith$gs*(df$vpd*1e3)*3600 #mol m-2 h-1

#### PHYDRO  WEEKLY ####
df %>% split(seq(nrow(.)))%>% 
  purrr::map(function(x){
    psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_sand_perc, sand = x$st_clay_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
    res <- pmodel_hydraulics_numerical(tc = x$ta, ppfd =(x$ppfd_in/ x$LAI), vpd = x$vpd*1000, u = x$ws_ERA, ustar=NA, nR=x$netrad, co2=x[["CO2"]], LAI = x$LAI,
                                       elv = x$si_elev, fapar = (x[["FAPAR"]]/ x$LAI), kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, par_cost = NULL, 
                                       opt_hypothesis = "PM", gs_approximation = "PM")
    return(res)
  })%>% bind_rows()->df_PM
colnames(df_PM) <- paste(colnames(df_PM), "PM", sep = "_")

df %>% split(seq(nrow(.)))%>% 
  purrr::map(function(x){
    psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_sand_perc, sand = x$st_clay_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
    res <- pmodel_hydraulics_numerical(tc = x$ta, ppfd =(x$ppfd_in/ x$LAI), vpd = x$vpd*1000, u = x$ws_ERA, ustar=NA, nR=x$netrad, co2=x[["CO2"]], LAI = x$LAI,
                                       elv = x$si_elev, fapar = (x[["FAPAR"]]/ x$LAI), kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, par_cost = NULL, 
                                       opt_hypothesis = "PM", gs_approximation = "Ohm")
    return(res)
  })%>% bind_rows()->df_Ohm
colnames(df_Ohm) <- paste(colnames(df_Ohm), "Ohm", sep = "_")

df_res <- df %>% 
  bind_cols(PMODEL = pm_E) %>% 
  bind_cols(`PMODEL swc limitation` = pm_E_swc) %>%  
  bind_cols(`PMODEL Smith` = pm_E_smith) %>% 
  # bind_cols(pm_E_swc_smith = pm_E_swc_smith)%>%
  bind_cols(df_Ohm) %>% 
  bind_cols(df_PM) %>%
  mutate(`E Sap flow`= E_stand/18.2, #mol m-2 h-1
         PHydro = E_Ohm/10^6*3600 / LAI, #transform to mol m-2soil h-1
         `PHydro Penman-Monteith` = E_PM/10^6*3600 / LAI,
         aet = aet * 55.5/24
         ) #transform to mol m-2soil h-1


df_res %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith','PHydro','aet','PHydro Penman-Monteith'
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

df_res %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith','PHydro','aet','PHydro Penman-Monteith'
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

df_res %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith','PHydro','aet','PHydro Penman-Monteith'
                        )) %>% 
  ggplot()+
  geom_line(aes(timestamp_aggr,value,color=name),size = 0.7,alpha=0.9, show.legend = FALSE)+
  geom_line(aes(timestamp_aggr,`E Sap flow`),color="grey10", size = 1)+
  ghibli::scale_color_ghibli_d("MononokeMedium",direction = -1)+
  ylab(expression(paste(mol[H2O]," ", m[soil]^{-2}," ", h^{-1})))+
  theme_classic()->P3

gridExtra::grid.arrange(P3,gridExtra::arrangeGrob(P1,P2,ncol = 2),ncol=1)



#### PMODEL MONTHLY ####
df <- sfn_aggr_month %>% 
  left_join(env_month) %>% 
  mutate(ppfd = ppfd_in/1e6,
         patm = rpmodel::calc_patm(FRA_PUE$si_elev %>% unique()),
         swvl = swvl_aggregation(st_soil_depth,swvl1,swvl2,swvl3,swvl4),
         swvl = opt_swc[1]*swvl+opt_swc[2],
         # swvl = swc_shallow,
         REW = (swvl-min(swvl,na.rm = TRUE))/(max(swvl,na.rm = TRUE)-min(swvl,na.rm = TRUE))) %>% 
  filter(!is.na(REW))
fit_ws <- lm(ws~0+ws_ERA, data =df) %>% summary()
df <- df %>% mutate(ws_ERA = case_when(is.na(ws)~ws_ERA*coefficients(fit_ws)[,'Estimate'],
                                       !is.na(ws)~ws)) %>% 
  dplyr::select(-ws)

pm_month <- df %>% 
  split(seq(nrow(.))) %>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=x[["CO2"]], fapar=x[["FAPAR"]], ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = FALSE, elv=FRA_PUE$si_elev %>% unique() )
  })

pm_month_smith <- df %>% 
  split(seq(nrow(.))) %>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=x[["CO2"]], fapar=x[["FAPAR"]], ppfd = x[["ppfd"]], soilm =x[["REW"]],
                     do_soilmstress = FALSE, elv=FRA_PUE$si_elev %>% unique(), method_jmaxlim = "smith19"  )
  })

pm_month_swc <- df %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
    rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=x[["CO2"]], fapar=x[["FAPAR"]], ppfd = x[["ppfd"]], soilm =x[["REW"]], do_soilmstress = TRUE, elv=FRA_PUE$si_elev %>% unique() )
  })

# R =  8.31446261815324 #universal gas constant m3 Pa K−1 mol−1
pm_E <- 1.6*pm_month$gs*(df$vpd*1e3)*3600 #mol m-2 h-1
pm_E_swc <- 1.6*pm_month_swc$gs*(df$vpd*1e3)*3600 #mol m-2 h-1
pm_E_smith <- 1.6*pm_month_smith$gs*(df$vpd*1e3)*3600 #mol m-2 h-1
# pm_E_swc_smith <- 1.6*pm_weekly_swc_smith$gs*(df$vpd*1e3)*3600 #mol m-2 h-1


#### PHYDRO  MONTHLY ####
df %>% split(seq(nrow(.)))%>% 
  purrr::map(function(x){
    psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_sand_perc, sand = x$st_clay_perc, om = NA, bd = 1.6, rfc = 70), W =x$swvl), model = "SX")
    res <- pmodel_hydraulics_numerical(tc = x$ta, ppfd =(x$ppfd_in/x$LAI), vpd = x$vpd*1000, u = x$ws_ERA, ustar=NA, nR=x$netrad, co2=400, LAI = x$LAI,
                                       elv = x$si_elev, fapar = (x[["FAPAR"]]/x$LAI), kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, rdark = 0, par_plant = par_plant_std, par_cost = NULL, 
                                       opt_hypothesis = "PM", gs_approximation = "Ohm")
    return(res)
  })%>% bind_rows()->df_Ohm
colnames(df_Ohm) <- paste(colnames(df_Ohm), "Ohm", sep = "_")

df_res <- df %>% 
  bind_cols(PMODEL = pm_E) %>% 
  bind_cols(`PMODEL swc limitation` = pm_E_swc) %>%  
  bind_cols(`PMODEL Smith` = pm_E_smith) %>% 
  # bind_cols(pm_E_swc_smith = pm_E_swc_smith)%>%
  bind_cols(df_Ohm) %>% 
  # bind_cols(df_PM) %>% 
  mutate(`E Sap flow`= E_stand/18.2, #mol m-2 h-1
         PHydro = E_Ohm/10^6*3600 / LAI#, #transform to mol m-2soil h-1
         # `PHydro Penman-Monteith` = E_PM*55.5*3600 / LAI
         ) #transform to mol m-2soil h-1


df_res %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith','PHydro'#,'PHydro Penman-Monteith'
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

df_res %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith','PHydro'#,'PHydro Penman-Monteith'
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

df_res %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith','PHydro'#,'PHydro Penman-Monteith'
                        )) %>% 
  ggplot()+
  geom_line(aes(timestamp_aggr,value,color=name),size = 0.7,alpha=0.9, show.legend = FALSE)+
  geom_line(aes(timestamp_aggr,`E Sap flow`),color="grey10", size = 1)+
  ghibli::scale_color_ghibli_d("MononokeMedium",direction = -1)+
  ylab(expression(paste(mol[H2O]," ", m[soil]^{-2}," ", h^{-1})))+
  theme_classic()->P3

gridExtra::grid.arrange(P3,gridExtra::arrangeGrob(P1,P2,ncol = 2),ncol=1)
