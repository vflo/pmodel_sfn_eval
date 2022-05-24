
calc_pmodel <- function(df_temp){
  df_temp %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, 
                              ppfd = x$ppfd_in/1e6, soilm =x$REW, do_soilmstress = FALSE, elv=sfn$si_elev %>% unique()) %>% 
        as_tibble()
      colnames(res) <- paste(colnames(res), "pmodel", sep = "_")
      res_temp <- x %>% bind_cols(res) %>% mutate(E_pmodel = 1.6*gs_pmodel*(vpd*1000)*3600) #mol m-2 h-1
      return(res_temp)
    })
}

calc_pmodel_swc <- function(df_temp, meanalpha){
  if(!is.na(meanalpha)){
    df_temp %>% 
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, ppfd = x$ppfd_in/1e6, soilm =x$REW,
                              do_soilmstress = TRUE, elv=sfn$si_elev %>% unique(), meanalpha = meanalpha) %>% 
          as_tibble()
        colnames(res) <- paste(colnames(res), "pmodel_swc", sep = "_")
        res_temp <- x %>% bind_cols(res) %>% mutate(E_pmodel_swc= 1.6*gs_pmodel_swc*(vpd*1000)*3600, #mol m-2 h-1
                                                    aet = aet * 55.5/24)
        return(res_temp)
    })
  }else{df_temp; print("Pmodel using soil water limitation was not computed")}
}


calc_phydro <- function(df_temp, PHYDRO_TRUE, par_plant_std, soil){
  if(PHYDRO_TRUE){
    df_temp %>%  
      filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd),FAPAR > 0, ppfd_in > 0) %>%
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        # print(x$timestamp_aggr)
        if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                           sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
                                                    W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
                                      model = "SX")
        res <- model_numerical(tc = x$ta, ppfd =x$ppfd_in/x$LAI, vpd = x$vpd*1000, nR=x$netrad, co2=x$CO2, LAI = x$LAI,elv = x$si_elev, 
                               fapar = (x$FAPAR),kphio = calc_ftemp_kphio(x$ta), psi_soil = psi_soil, par_plant = par_plant_std,
                               stomatal_model = "phydro") %>% 
          as_tibble()
        colnames(res) <- paste(colnames(res), "phydro", sep = "_")
        res_temp <- x %>% bind_cols(res) %>% mutate(E_phydro = E_phydro*3600*LAI) %>%  #transform to mol m-2soil h-1
          dplyr::select(1,2,42:51)

        return(res_temp)
      })
    }else{df_temp; print("Phydro was not computed")}
} 

calc_sperry <- function(df_temp, PHYDRO_TRUE, par_plant_std, soil){
  if(PHYDRO_TRUE){
    df_temp %>%  
      filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd),FAPAR > 0, ppfd_in > 0) %>%
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
        colnames(res) <- paste(colnames(res), "sperry", sep = "_")
        res_temp <- x %>% 
          bind_cols(res) %>% 
          mutate(E_sperry = E_sperry/1e3*3600*LAI) %>%  #transform to mol m-2soil h-1
          dplyr::select(1,2,42:47)
        return(res_temp)
      })
  }else{df_temp; print("Sperry model was not computed")}
} 
