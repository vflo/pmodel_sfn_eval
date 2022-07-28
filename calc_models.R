#### PMODEL MODEL ####
calc_pmodel <- function(df_temp, soil, par_plant){
  df_temp %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
      psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                         sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
                                                  W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
                                    model = "SX")
      res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, 
                              ppfd = x$ppfd_in, do_soilmstress = FALSE, elv=x$si_elev %>% unique()) %>% 
        as_tibble()
      
      
      res_temp <- x %>% 
        bind_cols(low_swp = FALSE) %>% 
        bind_cols(res) %>%
        # bind_cols(p_leaf = psi_l) %>%
        bind_cols(psi_soil = psi_soil) %>%
        bind_cols(model_type = "pmodel") %>% 
        mutate(E = 1.6*gs*(vpd*1000)*3600*1e-6,
               gs = gs) #umol m-2(ground) s-1 Pa-1
      return(res_temp)
    })
}

#### PMODEL SWC LIMITATION MODEL ####
calc_pmodel_swc <- function(df_temp, meanalpha, soil, par_plant){
  if(!is.na(meanalpha)){
    df_temp %>% 
      split(seq(nrow(.)))%>%
      purrr::map_df(function(x){
        if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
        psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                           sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
                                                    W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
                                      model = "SX")
        FC <- medfate::soil_thetaFC(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
                                                         sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70)))
        res <- rpmodel::rpmodel(tc=x$ta, vpd=x$vpd*1000, co2=x$CO2, fapar=x$FAPAR, ppfd = x$ppfd_in, soilm =x$REW,
                                do_soilmstress = TRUE, elv=x$si_elev %>% unique(), meanalpha = meanalpha) %>% 
          as_tibble()
        
        res_temp <- x %>% 
          bind_cols(low_swp = FALSE) %>% 
          bind_cols(res) %>% 
          # bind_cols(p_leaf = psi_l) %>%
          bind_cols(psi_soil = psi_soil) %>%
          bind_cols(model_type = "pmodel_swc") %>% 
          mutate(E = 1.6*gs*(vpd*1000)*3600*1e-6, #mol m-2 h-1
                 gs = gs) # umol m-2(ground) s-1 Pa-1
        return(res_temp)
      })
  }else{df_temp; print("Pmodel using soil water limitation was not computed")}
}


#### PMODEL SUBDAILY MODEL ####
calc_pmodel_subdaily <- function(df_temp, soil, par_plant){
  df_temp %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
      medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, 
                                             clay = x$st_clay_perc,
                                             sand = x$st_sand_perc, 
                                             om = soil$OM, 
                                             bd = soil$bd, 
                                             rfc = 70),
                                      W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, 
                                                                sand = x$st_sand_perc, 
                                                                om = soil$OM)), 
                        model = "SX")->psi_soil
      return(tibble(psi_soil = psi_soil))
      })->psi_soil
      
      res <- rpmodel_subdaily(TIMESTAMP = df_temp$TIMESTAMP, tc=df_temp$ta, vpd=df_temp$vpd*1000, co2=df_temp$CO2, 
                              fapar= df_temp$FAPAR, ppfd = df_temp$ppfd_in, sw_in = df_temp$sw_in,
                              do_soilmstress = FALSE, elv=df_temp$si_elev %>% unique()) %>% 
        as_tibble()
      
      
      res_temp <- df_temp %>% 
        bind_cols(low_swp = FALSE) %>% 
        bind_cols(res) %>%
        # bind_cols(p_leaf = psi_l) %>%
        bind_cols(psi_soil = psi_soil) %>%
        bind_cols(model_type = "pmodel") %>% 
        mutate(E = 1.6*gs*(vpd*1000)*3600*1e-6,
               gs = gs) #umol m-2(ground) s-1 Pa-1
      return(res_temp)
}

#### PHYDRO SCHEMES MODELS ####
calc_phydro_schemes <- function(df_temp, PHYDRO_TRUE, par_plant, soil, sensi, stomatal_model){
  #filter
  temp <- df_temp %>%  
    filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),  !is.na(vpd), vpd>0,FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
  
  #sensitivity analysis
  par_plant$psi50 <- par_plant$psi50 * sensi$sensitivity_psi
  par_plant$conductivity <- par_plant$conductivity * sensi$sensitivity_K
  
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
        res <- model_numerical(tc             = x$ta, 
                               ppfd           = x$ppfd_in, 
                               vpd            = x$vpd*1000, 
                               co2            = x$CO2, 
                               LAI            = x$LAI,
                               elv            = x$si_elev, 
                               fapar          = (x$FAPAR),
                               kphio          = calc_ftemp_kphio(x$ta), 
                               psi_soil       = psi_soil, 
                               par_plant      = par_plant,
                               stomatal_model = stomatal_model) %>% 
          as_tibble()
        
        
        res_temp <- x %>% 
          bind_cols(res) %>% 
          bind_cols(sensi) %>% 
          bind_cols(psi_soil = psi_soil) %>%
          bind_cols(model_type = stomatal_model) %>% 
          mutate(E = E*3600,
                 gs_alt = E/(1.6*(vpd*1000)*3600*1e-6)) #transform to mol m-2 (ground) h-1
        
        return(res_temp)
      })
  }else{
    df_temp; print("PHYDRO model was not computed")
    return(temp)
  }
} 
# 
# #### PHYDRO ANALYTICAL MODEL ####
# calc_phydro_a <- function(df_temp, PHYDRO_TRUE, par_plant, soil, sensi){
#   #filter
#   temp <- df_temp %>%  
#     filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),  !is.na(vpd), vpd>0,FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
#   
#   #sensitivity analysis
#   par_plant$psi50 <- par_plant$psi50 * sensi$sensitivity_psi
#   par_plant$conductivity <- par_plant$conductivity * sensi$sensitivity_K
#   
#   if(PHYDRO_TRUE && nrow(temp)>0){
#     temp%>%
#       split(seq(nrow(.)))%>%
#       purrr::map_df(function(x){
#         # print(x$timestamp_aggr)
#         if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
#         psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
#                                                            sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
#                                                     W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
#                                       model = "SX")
#         res <- pmodel_hydraulics_analytical(
#                                tc             = x$ta, 
#                                ppfd           = x$ppfd_in, 
#                                vpd            = x$vpd*1000, 
#                                co2            = x$CO2, 
#                                LAI            = x$LAI,
#                                elv            = x$si_elev, 
#                                fapar          = (x$FAPAR),
#                                kphio          = calc_ftemp_kphio(x$ta), 
#                                rdark = 0,
#                                psi_soil       = psi_soil, 
#                                par_plant      = par_plant,
#                                par_cost       = NULL) %>% 
#           as_tibble()
#         
#         
#         res_temp <- x %>% 
#           bind_cols(res) %>% 
#           bind_cols(sensi) %>% 
#           bind_cols(psi_soil = psi_soil) %>%
#           bind_cols(model_type = "phydro_a") %>% 
#           mutate(E = E*3600,
#                  gs_alt = E/(1.6*(vpd*1000)*3600*1e-6)) #transform to mol m-2 (ground) h-1
#         
#         return(res_temp)
#       })
#   }else{
#     df_temp; print("PHYDRO model was not computed")
#     return(temp)
#   }
# } 


#### PHYDRO WANG MODEL ####
# calc_phydro_wang <- function(df_temp, PHYDRO_TRUE, par_plant, soil, sensi){
#   #filter
#   temp <- df_temp %>%  
#     filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in), !is.na(vpd), vpd>0,
#            FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
#   
#   #sensitivity analysis
#   par_plant$psi50 <- par_plant$psi50 * sensi$sensitivity_psi
#   par_plant$conductivity <- par_plant$conductivity * sensi$sensitivity_K
#   
#   if(PHYDRO_TRUE && nrow(temp)>0){
#     temp%>%
#       split(seq(nrow(.)))%>%
#       purrr::map_df(function(x){
#         # print(x$timestamp_aggr)
#         if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
#         psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
#                                                            sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
#                                                     W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
#                                       model = "SX")
#         res <- model_numerical(tc             = x$ta, 
#                                ppfd           = x$ppfd_in, 
#                                vpd            = x$vpd*1000, 
#                                co2            = x$CO2, 
#                                LAI            = x$LAI,
#                                elv            = x$si_elev, 
#                                fapar          = (x$FAPAR),
#                                kphio          = calc_ftemp_kphio(x$ta), 
#                                psi_soil       = psi_soil, 
#                                par_plant      = par_plant,
#                                stomatal_model = "phydro_wang") %>% 
#           as_tibble()
#         
#         
#         res_temp <- x %>% 
#           bind_cols(res) %>% 
#           bind_cols(sensi) %>% 
#           bind_cols(psi_soil = psi_soil) %>%
#           bind_cols(model_type = "phydro_wang") %>% 
#           mutate(E = E*3600,
#                  gs_alt = E/(1.6*(vpd*1000)*3600*1e-6)) #transform to mol m-2 (ground) h-1
#         
#         return(res_temp)
#       })
#   }else{
#     df_temp; print("PHYDRO - WANG model was not computed")
#     return(temp)
#   }
# } 
# 
# 
# #### PHYDRO WAP MODEL ####
# calc_phydro_wap <- function(df_temp, PHYDRO_TRUE, par_plant, soil, sensi){
#   #filter
#   temp <- df_temp %>%  
#     filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),  !is.na(vpd), vpd>0,FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
#   
#   #sensitivity analysis
#   par_plant$psi50 <- par_plant$psi50 * sensi$sensitivity_psi
#   par_plant$conductivity <- par_plant$conductivity * sensi$sensitivity_K
#   
#   if(PHYDRO_TRUE && nrow(temp)>0){
#     temp%>%
#       split(seq(nrow(.)))%>%
#       purrr::map_df(function(x){
#         # print(x$timestamp_aggr)
#         if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
#         psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
#                                                            sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
#                                                     W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
#                                       model = "SX")
#         res <- model_numerical(tc             = x$ta, 
#                                ppfd           = x$ppfd_in, 
#                                vpd            = x$vpd*1000, 
#                                co2            = x$CO2, 
#                                LAI            = x$LAI,
#                                elv            = x$si_elev, 
#                                fapar          = (x$FAPAR),
#                                kphio          = calc_ftemp_kphio(x$ta), 
#                                psi_soil       = psi_soil, 
#                                par_plant      = par_plant,
#                                stomatal_model = "phydro_wap") %>% 
#           as_tibble()
#         
#         
#         res_temp <- x %>% 
#           bind_cols(res) %>% 
#           bind_cols(sensi) %>% 
#           bind_cols(psi_soil = psi_soil) %>%
#           bind_cols(model_type = "phydro_wap") %>% 
#           mutate(E = E*3600,
#                  gs_alt = E/(1.6*(vpd*1000)*3600*1e-6)) #transform to mol m-2 (ground) h-1
#         
#         return(res_temp)
#       })
#   }else{
#     df_temp; print("PHYDRO - WANG model was not computed")
#     return(temp)
#   }
# } 

# #### SPERRY MODEL ####
# calc_sperry <- function(df_temp, PHYDRO_TRUE, par_plant, soil, sensi){
#   #filter
#   temp <- df_temp %>%  
#     filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),  !is.na(vpd), vpd>0,FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
#   
#   #sensitivity analysis
#   par_plant$d <- par_plant$d * sensi$sensitivity_psi
#   par_plant$conductivity <- par_plant$conductivity * sensi$sensitivity_K
#   
#   if(PHYDRO_TRUE && nrow(temp)>0){
#     temp%>%
#       split(seq(nrow(.)))%>%
#       purrr::map_df(function(x){
#         # print(x$timestamp_aggr)
#         if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
#         psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
#                                                            sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
#                                                     W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
#                                       model = "SX")
#         res <- model_numerical(tc             = x$ta, 
#                                ppfd           = x$ppfd_in, 
#                                vpd            = x$vpd*1000,
#                                co2            = x$CO2, 
#                                LAI            = x$LAI,
#                                elv            = x$si_elev, 
#                                fapar          = (x$FAPAR),
#                                kphio          = calc_ftemp_kphio(x$ta), 
#                                psi_soil       = psi_soil, 
#                                par_plant      = par_plant, 
#                                stomatal_model = "sperry") %>% 
#           as_tibble()
#         
#         res_temp <- x %>% 
#           bind_cols(res) %>% 
#           bind_cols(sensi) %>%
#           bind_cols(psi_soil = psi_soil) %>%
#           bind_cols(model_type = "sperry") %>% 
#           mutate(E = E*3600, #transform to mol m-2 (ground) h-1
#                  # gs_org = gs,
#                  gs = gs/patm*1e6 ,
#                  gs_alt = E/(1.6*(vpd*1000)*3600*1e-6)
#                  ) #transform to umol (Co2) m-2 (ground) s-1 Pa-1
#         return(res_temp)
#       })
#   }else{
#     df_temp; print("SPERRY model was not computed")
#     return(temp)
#   }
# } 
# 
# #### WANG MODEL ####
# calc_wang <- function(df_temp, PHYDRO_TRUE, par_plant, soil, sensi){
#   #filter
#   temp <- df_temp %>%  
#     filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),  !is.na(vpd), vpd>0,FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
#   
#   #sensitivity analysis
#   par_plant$d <- par_plant$d * sensi$sensitivity_psi
#   par_plant$conductivity <- par_plant$conductivity * sensi$sensitivity_K
#   
#   if(PHYDRO_TRUE && nrow(temp)>0){
#     temp%>%
#       split(seq(nrow(.)))%>%
#       purrr::map_df(function(x){
#         # print(x$timestamp_aggr)
#         if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
#         psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
#                                                            sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
#                                                     W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
#                                       model = "SX")
#         res <- model_numerical(tc = x$ta, 
#                                ppfd =x$ppfd_in, 
#                                vpd = x$vpd*1000,
#                                co2=x$CO2, 
#                                LAI = x$LAI,
#                                elv = x$si_elev, 
#                                fapar = (x$FAPAR),
#                                kphio = calc_ftemp_kphio(x$ta), 
#                                psi_soil = psi_soil, 
#                                # rdark = 0, 
#                                par_plant = par_plant, 
#                                stomatal_model = "wang") %>% 
#           as_tibble()
# 
#         res_temp <- x %>% 
#           bind_cols(res) %>% 
#           bind_cols(sensi) %>%
#           bind_cols(psi_soil = psi_soil) %>%
#           bind_cols(model_type = "wang") %>% 
#           mutate(E = E*3600, #transform to mol m-2 (ground) h-1
#                  gs = gs/patm*1e6,
#                  gs_alt = E/(1.6*(vpd*1000)*3600*1e-6))
#         
#         return(res_temp)
#       })
#   }else{
#     df_temp; print("WANG model was not computed")
#     return(temp)
#   }
# } 
# 
# 
# #### WAP MODEL ####
# calc_wap <- function(df_temp, PHYDRO_TRUE, par_plant, soil, sensi){
#   #filter
#   temp <- df_temp %>%  
#     filter(!is.na(ta),!is.na(FAPAR),!is.na(ppfd_in),  !is.na(vpd), vpd>0,FAPAR > 0, ppfd_in > 0, LAI>0, swvl > 0.01)
#   
#   #sensitivity analysis
#   par_plant$d <- par_plant$d * sensi$sensitivity_psi
#   par_plant$conductivity <- par_plant$conductivity * sensi$sensitivity_K
#   
#   if(PHYDRO_TRUE && nrow(temp)>0){
#     temp%>%
#       split(seq(nrow(.)))%>%
#       purrr::map_df(function(x){
#         # print(x$timestamp_aggr)
#         if( is.na(x$st_soil_depth)){ x$st_soil_depth <- soil$depth*100}
#         psi_soil <- medfate::soil_psi(medfate::soil(tibble(widths = x$st_soil_depth*10, clay = x$st_clay_perc, 
#                                                            sand = x$st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
#                                                     W =x$swvl/soil_thetaSATSX(clay = x$st_clay_perc, sand = x$st_sand_perc, om = soil$OM)), 
#                                       model = "SX")
#         res <- model_numerical(tc = x$ta, 
#                                ppfd =x$ppfd_in, 
#                                vpd = x$vpd*1000,
#                                co2=x$CO2, 
#                                LAI = x$LAI,
#                                elv = x$si_elev, 
#                                fapar = (x$FAPAR),
#                                kphio = calc_ftemp_kphio(x$ta), 
#                                psi_soil = psi_soil, 
#                                # rdark = 0, 
#                                par_plant = par_plant, 
#                                stomatal_model = "wap") %>% 
#           as_tibble()
# 
#         
#         res_temp <- x %>% 
#           bind_cols(res) %>% 
#           bind_cols(sensi) %>% 
#           bind_cols(psi_soil = psi_soil) %>%
#           bind_cols(model_type = "wap") %>% 
#           mutate(E = E*3600, #transform to mol m-2 (ground) h-1
#                  gs = gs/patm*1e6,
#                  gs_alt = E/(1.6*(vpd*1000)*3600*1e-6))
#         
#         return(res_temp)
#       })
#   }else{
#     df_temp; print("WAP model was not computed")
#     return(temp)
#   }
#   
# } 