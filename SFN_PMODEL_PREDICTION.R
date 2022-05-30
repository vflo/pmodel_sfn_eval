source("init_SFN_MODEL_PREDICTION.R")

#### MODEL CALCULATION ####
as.list(list_files) %>%
  purrr::map(function(x){
    load(x)
    print(gsub(".*/","",x = x))
    
    sfn <- sfn_list$sfn
    soil <- sfn_list$soil
    si_code <- sfn_list$sfn$si_code %>% unique()
    opt_swc <- sfn_list$opt_swc
    PHYDRO_TRUE <- sfn_list$PHYDRO_TRUE
    par_plant_std <- sfn_list$par_plant_std
    
  #### SPATIAL AND TEMPORAL DATA AGGREGATION ####
    sfn_agr <- calc_sfn_scale(sfn)
    sfn_aggr_week <- calc_sfn_aggr_E(sfn_agr,"1 week")
    env_week <- calc_sfn_aggr_env(sfn,"1 week")
    env_week <- env_week %>% mutate(CO2 = oce::fillGap(CO2))
    sfn_aggr_month <- calc_sfn_aggr_E(sfn_agr,"1 month")
    env_month <- calc_sfn_aggr_env(sfn,"1 month")
    env_month <- env_month %>% mutate(CO2 = oce::fillGap(CO2))
    # sfn_aggr_day <- calc_sfn_aggr_E(sfn_agr,"1 day")
    # env_day <- calc_sfn_aggr_env(sfn,"1 day")
    # env_day <- env_day %>% mutate(CO2 = oce::fillGap(CO2))

  # CALCULATE MEANALPHA
  
    meanalpha <- mean_alpha_calc(sfn, sw_in_monthly_average, temp_monthly_average, precip_monthly_average, env_month, soil)
  
  #### PMODEL WEEKLY ####
  #DATA preparation
    df <- data_prep(sfn_aggr_week, env_week, opt_swc) %>% bind_cols(si_code = si_code)

    # par_plant_std$conductivity <- par_plant_std$conductivity*0.15
    # par_plant_std$d <- par_plant_std$d*1.2
    sensitivity <- tibble(sensitivity_K = c(seq(0.5,1.5,0.1), rep(1,11)),
                          sensitivity_psi = c(rep(1,11), seq(0.5,1.5,0.1))) %>% distinct()

  #### MODEL CALCULATION ####
  ## PMODEL ##
    print("Calculating PMODEL")
    pmodel <- calc_pmodel(df)%>% as_tibble()

  # PMODEL SWC limitation if there is AET/PET then calculate it
    print("Calculating PMODEL with swc limitation")
    pmodel_swc <- calc_pmodel_swc(df, meanalpha, soil)%>% as_tibble()

  # PHydro if there is species information 
    print("Calculating PHYDRO")
    sensitivity %>% 
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        calc_phydro(df, PHYDRO_TRUE, par_plant_std, soil, x) %>% as_tibble()
      })  %>%  bind_rows() -> phydro 
    
  # Sperry model
    print("Calculating SPERRY model")
    sensitivity %>% 
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        calc_sperry(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
      }) %>%  bind_rows() -> sperry
    
  # Wang model
    print("Calculating WANG model")
    sensitivity %>% 
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        calc_wang(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
      })  %>%  bind_rows() -> wang
    
  # Wap model
    print("Calculating WAP model")
    sensitivity %>% 
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        calc_wap(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
      })  %>%  bind_rows() -> wap
    
    ## PMODEL Ecrit ##
    print("Calculating PMODEL Ecrit")
    pmodel_ecrit <- calc_pmodel_ecrit(df, PHYDRO_TRUE, par_plant_std, soil)%>% as_tibble()
    
    df_res <- pmodel %>% 
      bind_rows(pmodel_swc) %>% 
      bind_rows(phydro) %>% 
      bind_rows(sperry) %>% 
      bind_rows(wang) %>% 
      bind_rows(wap) %>%
      suppressMessages()
# 
#     df_res %>% 
#       ggplot()+
#       geom_ribbon(aes(x=TIMESTAMP, ymin = (E_sapflow - E_sapflow_sd), ymax = (E_sapflow + E_sapflow_sd)), fill = "grey" )+
#       geom_line(aes(TIMESTAMP,E_sapflow))+
#       geom_line(aes(TIMESTAMP,E_pmodel), color="brown")+
#       geom_line(aes(TIMESTAMP,E_pmodel_swc), color="purple")+
#       geom_line(aes(TIMESTAMP,E_phydro), color="green")+
#       geom_line(aes(TIMESTAMP,E_sperry), color="orange")+
#       geom_line(aes(TIMESTAMP,E_wang), color="yellow")+
#       geom_line(aes(TIMESTAMP,E_wap), color="red")+
#       geom_line(aes(TIMESTAMP,aet), color="blue")+
#       theme_bw()
    
  #### SAVE DATA ####
    if(nrow(df_res)>0){
      save(df_res, file = paste0("DATA/OUTCOME_WEEKLY/",(gsub(".*/","",x = x))))
  }
  
  
  
  #### PMODEL MONTHLY ####
    df <- data_prep(sfn_aggr_month, env_month, opt_swc)
  
    if(nrow(df)>1 && any(!is.na(df$E_stand))){

  #### MODEL CALCULATION ####
      ## PMODEL ##
      pmodel <- calc_pmodel(df)%>% as_tibble()
      
      # PMODEL SWC limitation if there is AET/PET then calculate it
      pmodel_swc <- calc_pmodel_swc(df, meanalpha, soil)%>% as_tibble()
      
      # PHydro if there is species information 
      sensitivity %>% 
        split(seq(nrow(.)))%>%
        purrr::map(function(x){
          calc_phydro(df, PHYDRO_TRUE, par_plant_std, soil, x) %>% as_tibble()
        })  %>%  bind_rows() -> phydro 
      
      # Sperry model
      sensitivity %>% 
        split(seq(nrow(.)))%>%
        purrr::map(function(x){
          calc_sperry(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
        }) %>%  bind_rows() -> sperry
      
      # Wang model
      sensitivity %>% 
        split(seq(nrow(.)))%>%
        purrr::map(function(x){
          calc_wang(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
        })  %>%  bind_rows() -> wang
      
      # Wap model
      sensitivity %>% 
        split(seq(nrow(.)))%>%
        purrr::map(function(x){
          calc_wap(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
        })  %>%  bind_rows() -> wap
  
      df_res <- pmodel %>% 
        bind_rows(pmodel_swc) %>% 
        bind_rows(phydro) %>% 
        bind_rows(sperry) %>% 
        bind_rows(wang) %>% 
        bind_rows(wap) %>%
        suppressMessages()

    if(nrow(df_res)>0){
      save(df_res, file = paste0("DATA/OUTCOME_MONTHLY/",(gsub(".*/","",x = x))))
    }
  #return(df_res)
  }
})#->DF
