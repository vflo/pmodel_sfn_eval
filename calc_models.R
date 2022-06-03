#### PMODEL MODEL ####
calc_pmodel <- function(df_temp){
  df_temp %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, 
                              ppfd = x$ppfd_in/1e6, do_soilmstress = FALSE, elv=x$si_elev %>% unique()) %>% 
        as_tibble()

      res_temp <- x %>% 
        bind_cols(low_swp = FALSE) %>% 
        bind_cols(res) %>% 
        bind_cols(model_type = "pmodel") %>% 
        mutate(E = 1.6*gs*(vpd*1000)*3600) #mol m-2 h-1
      return(res_temp)
    })
}

#### PMODEL SWC LIMITATION MODEL ####
calc_pmodel_swc <- function(df_temp, meanalpha, soil){
  if(!is.na(meanalpha)){
    df_temp %>% 
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        FC <- medfate::soil_thetaFC(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                   sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70)))
        res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, ppfd = x$ppfd_in/1e6, soilm =x$swvl/FC,
                              do_soilmstress = TRUE, elv=x$si_elev %>% unique(), meanalpha = meanalpha) %>% 
          as_tibble()
        res_temp <- x %>% 
          bind_cols(low_swp = FALSE) %>% 
          bind_cols(res) %>% 
          bind_cols(model_type = "pmodel_swc") %>% 
          mutate(E = 1.6*gs*(vpd*1000)*3600, #mol m-2 h-1
                             aet = aet * 55.5/24) #transform to mol m-2soil h-1
        return(res_temp)
    })
  }else{df_temp; print("Pmodel using soil water limitation was not computed")}
}

#### PHYDRO MODEL ####
calc_phydro <- function(df_temp, PHYDRO_TRUE, par_plant_std, soil, sensi){
  #filter
  temp <- df_temp %>%  
    filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd),FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
  
    #sensitivity analysis
    par_plant_std$psi50 <- par_plant_std$psi50 * sensi$sensitivity_psi
    par_plant_std$conductivity <- par_plant_std$conductivity * sensi$sensitivity_K
    
  if(PHYDRO_TRUE && nrow(temp)>0){
     temp%>%
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        # print(x$timestamp_aggr)
        if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                           sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
                                                    W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
                                      model = "SX")
        res <- model_numerical(tc = x$ta, 
                               ppfd =x$ppfd_in/x$LAI, 
                               vpd = x$vpd*1000, 
                               nR =x$netrad, 
                               co2=x$CO2, 
                               LAI = x$LAI,
                               elv = x$si_elev, 
                               fapar = (x$FAPAR),
                               kphio = calc_ftemp_kphio(x$ta), 
                               psi_soil = psi_soil, 
                               par_plant = par_plant_std,
                               stomatal_model = "phydro") %>% 
          as_tibble()

        res_temp <- x %>% 
          bind_cols(res) %>% 
          bind_cols(sensi) %>%
          bind_cols(model_type = "phydro") %>% 
          mutate(E = E*3600*LAI) #transform to mol m-2soil h-1
          
        return(res_temp)
      })
  }else{
    df_temp; print("PHYDRO model was not computed")
    return(temp)
  }
} 

#### SPERRY MODEL ####
calc_sperry <- function(df_temp, PHYDRO_TRUE, par_plant_std, soil, sensi){
  #filter
  temp <- df_temp %>%  
    filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd),FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
  
  #sensitivity analysis
  par_plant_std$d <- par_plant_std$d * sensi$sensitivity_psi
  par_plant_std$conductivity <- par_plant_std$conductivity * sensi$sensitivity_K
  
  if(PHYDRO_TRUE && nrow(temp)>0){
    temp%>%
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        # print(x$timestamp_aggr)
        if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                           sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
                                                    W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
                                      model = "SX")
        res <- model_numerical(tc = x$ta, 
                               ppfd =x$ppfd_in/x$LAI, 
                               vpd = x$vpd*1000,
                               nR=x$netrad, 
                               co2=x$CO2, 
                               LAI = x$LAI,
                               elv = x$si_elev, 
                               fapar = (x$FAPAR),
                               kphio = calc_ftemp_kphio(x$ta), 
                               psi_soil = psi_soil, 
                               # rdark = 0, 
                               par_plant = par_plant_std, 
                               stomatal_model = "sperry") %>% 
          as_tibble()

        res_temp <- x %>% 
          bind_cols(res) %>% 
          bind_cols(sensi) %>% 
          bind_cols(model_type = "sperry") %>% 
          mutate(E = E/1e3*3600*LAI) #transform to mol m-2soil h-1
        return(res_temp)
      })
  }else{
    df_temp; print("SPERRY model was not computed")
    return(temp)
  }
} 

#### WANG MODEL ####
calc_wang <- function(df_temp, PHYDRO_TRUE, par_plant_std, soil, sensi){
  #filter
  temp <- df_temp %>%  
    filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd),FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
  
  #sensitivity analysis
  par_plant_std$d <- par_plant_std$d * sensi$sensitivity_psi
  par_plant_std$conductivity <- par_plant_std$conductivity * sensi$sensitivity_K
  
  if(PHYDRO_TRUE && nrow(temp)>0){
    temp%>%
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        # print(x$timestamp_aggr)
        if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                           sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
                                                    W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
                                      model = "SX")
        res <- model_numerical(tc = x$ta, 
                               ppfd =x$ppfd_in/x$LAI, 
                               vpd = x$vpd*1000,
                               nR=x$netrad, 
                               co2=x$CO2, 
                               LAI = x$LAI,
                               elv = x$si_elev, 
                               fapar = (x$FAPAR),
                               kphio = calc_ftemp_kphio(x$ta), 
                               psi_soil = psi_soil, 
                               # rdark = 0, 
                               par_plant = par_plant_std, 
                               stomatal_model = "wang") %>% 
          as_tibble()
        
        res_temp <- x %>% 
          bind_cols(res) %>% 
          bind_cols(sensi) %>%
          bind_cols(model_type = "wang") %>% 
          mutate(E = E/1e3*3600*LAI)#transform to mol m-2soil h-1
        
        return(res_temp)
      })
  }else{
    df_temp; print("WANG model was not computed")
    return(temp)
  }
} 


#### WAP MODEL ####
calc_wap <- function(df_temp, PHYDRO_TRUE, par_plant_std, soil, sensi){
  #filter
  temp <- df_temp %>%  
    filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd),FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
  
  #sensitivity analysis
  par_plant_std$d <- par_plant_std$d * sensi$sensitivity_psi
  par_plant_std$conductivity <- par_plant_std$conductivity * sensi$sensitivity_K
  
  if(PHYDRO_TRUE && nrow(temp)>0){
    temp%>%
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        # print(x$timestamp_aggr)
        if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                           sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
                                                    W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
                                      model = "SX")
        res <- model_numerical(tc = x$ta, 
                               ppfd =x$ppfd_in/x$LAI, 
                               vpd = x$vpd*1000,
                               nR=x$netrad, 
                               co2=x$CO2, 
                               LAI = x$LAI,
                               elv = x$si_elev, 
                               fapar = (x$FAPAR),
                               kphio = calc_ftemp_kphio(x$ta), 
                               psi_soil = psi_soil, 
                               # rdark = 0, 
                               par_plant = par_plant_std, 
                               stomatal_model = "wap") %>% 
          as_tibble()

        res_temp <- x %>% 
          bind_cols(res) %>% 
          bind_cols(sensi) %>% 
          bind_cols(model_type = "wap") %>% 
          mutate(E = E /1e3*3600*LAI)  #transform to mol m-2soil h-1
          
        return(res_temp)
      })
  }else{
    df_temp; print("WAP model was not computed")
    return(temp)
    }
  
} 




#### PMODEL ECRIT ####
calc_pmodel_ecrit <- function(df_temp, PHYDRO_TRUE, par_plant_std, soil, sensi){
  #filter
  temp <- df_temp %>%  
    filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd),FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
  
  #sensitivity analysis
  par_plant_std$psi50 <- par_plant_std$psi50 * sensi$sensitivity_psi
  par_plant_std$conductivity <- par_plant_std$conductivity * sensi$sensitivity_K
  
  if(PHYDRO_TRUE && nrow(temp)>0){
    temp%>%
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        # print(x$timestamp_aggr)
        if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                           sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
                                                    W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
                                      model = "SX")
        res <- rpmodel_ecrit(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, psi_soil = psi_soil,
                             K = par_plant_std$conductivity, d = par_plant_std$d, c = par_plant_std$c, h = par_plant_std$height,
                             rs = 0.1, LAI = x$LAI, ppfd = x$ppfd_in, do_soilmstress = FALSE, elv=x$si_elev %>% unique()) %>% 
          as_tibble()
        
        res_temp <- x %>% 
          bind_cols(res) %>% 
          bind_cols(sensi) %>%
          bind_cols(model_type = "pmodel_ecrit") %>% 
          mutate(E = E*3600*1e-6) #transform to mol m-2soil h-1
        
        return(res_temp)
      })
  }else{
    df_temp; print("PHYDRO model was not computed")
    return(temp)
  }
} 
