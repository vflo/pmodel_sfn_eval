source("init_SFN_MODEL_PREDICTION.R")

#### MODEL CALCULATION ####

  
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

  # CALCULATE MEANALPHA
  meanalpha <- mean_alpha_calc(sfn, sw_in_monthly_average, temp_monthly_average, precip_monthly_average, env_month, soil)
  
  #### PMODEL WEEKLY ####
  #DATA preparation
  df <- data_prep(sfn_aggr_week, env_week)


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
  
if(PHYDRO_TRUE){
  df %>%  
    filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd),FAPAR > 0, ppfd_in > 0) %>%
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
  
  df_res <- df_res %>% 
    left_join(df_Ohm, by = "timestamp_aggr") %>%
    mutate(PHydro = E_Ohm/1e6*3600*LAI) #transform to mol m-2soil h-1

  #if there is net radiation then calculate phydro with penman monteith
  if(any(!is.na(df$netrad))){
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
      mutate(`PHydro Penman-Monteith` = E_PM*3600*LAI) #transform to mol m-2soil h-1

  }

}
  #### SAVE DATA ####
  if(nrow(df_res)>0){
    save(df_res, file = paste0("DATA/OUTCOME_WEEKLY/",x,".RData"))
  }
  
  
  
  #### PMODEL MONTHLY ####
  df <- data_prep(sfn_aggr_month, env_month)
  
  
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
  
  if(PHYDRO_TRUE){
    df %>%  
      filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd), FAPAR > 0, ppfd_in > 0) %>%
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
    
    df_res <- df_res %>% 
      left_join(df_Ohm, by = "timestamp_aggr") %>%
      mutate(PHydro = E_Ohm/1e6*3600*LAI) #transform to mol m-2soil h-1
    
    #if there is net radiation then calculate phydro with penman monteith
    if(any(!is.na(df$netrad))){
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
        mutate(`PHydro Penman-Monteith` = E_PM*3600*LAI) #transform to mol m-2soil h-1
      
    }
    
  }
  #### REMOVE DATA ####
  rm(sfn)
  gc()
  
  if(nrow(df_res)>0){
    save(df_res, file = paste0("DATA/OUTCOME_MONTHLY/",x,".RData"))
  }
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

