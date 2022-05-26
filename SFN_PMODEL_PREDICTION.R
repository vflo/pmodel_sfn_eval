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
    sfn_aggr_day <- calc_sfn_aggr_E(sfn_agr,"1 day")
    env_day <- calc_sfn_aggr_env(sfn,"1 day")
    env_day <- env_day %>% mutate(CO2 = oce::fillGap(CO2))

  # CALCULATE MEANALPHA
  
    meanalpha <- mean_alpha_calc(sfn, sw_in_monthly_average, temp_monthly_average, precip_monthly_average, env_month, soil)
  
  #### PMODEL WEEKLY ####
  #DATA preparation
    df <- data_prep(sfn_aggr_week, env_week, opt_swc) %>% bind_cols(si_code = si_code)


  #### MODEL CALCULATION ####
  ## PMODEL ##
    pmodel <- calc_pmodel(df)

  # PMODEL SWC limitation if there is AET/PET then calculate it
    pmodel_swc <- calc_pmodel_swc(df, meanalpha, soil)

  # PHydro if there is species information  
    phydro <- calc_phydro(df, PHYDRO_TRUE, par_plant_std, soil)
  
  # Sperry model
    sperry <- calc_sperry(df, PHYDRO_TRUE, par_plant_std, soil)
  
  # Wang model
    wang <- calc_wang(df, PHYDRO_TRUE, par_plant_std, soil)

  # Wap model
    wap <- calc_wap(df, PHYDRO_TRUE, par_plant_std, soil)
  
    df_res <- pmodel %>% 
      left_join(pmodel_swc) %>% 
      left_join(phydro) %>% 
      left_join(sperry) %>% 
      left_join(wang) %>% 
      left_join(wap) %>%
      suppressMessages()


  #### SAVE DATA ####
    if(nrow(df_res)>0){
      save(df_res, file = paste0("DATA/OUTCOME_WEEKLY/",(gsub(".*/","",x = x))))
  }
  
  
  
  #### PMODEL MONTHLY ####
    df <- data_prep(sfn_aggr_month, env_month, opt_swc)
  
    if(nrow(df)>1 && any(!is.na(df$E_stand))){

  #### MODEL CALCULATION ####
  ## PMODEL ##
    pmodel <- calc_pmodel(df)
  
  # PMODEL SWC limitation if there is AET/PET then calculate it
    pmodel_swc <- calc_pmodel_swc(df, meanalpha, soil)
  
  # PHydro if there is species information  
    phydro <- calc_phydro(df, PHYDRO_TRUE, par_plant_std, soil)
  
  # Sperry model
    sperry <- calc_sperry(df, PHYDRO_TRUE, par_plant_std, soil)
  
  # Wang model
    wang <- calc_wang(df, PHYDRO_TRUE, par_plant_std, soil)
  
  # Wap model
    wap <- calc_wap(df, PHYDRO_TRUE, par_plant_std, soil)
  
    df_res <- pmodel %>% 
      left_join(pmodel_swc) %>% 
      left_join(phydro) %>% 
      left_join(sperry) %>% 
      left_join(wang) %>% 
      left_join(wap) %>%
      suppressMessages()
  

    if(nrow(df_res)>0){
      save(df_res, file = paste0("DATA/OUTCOME_MONTHLY/",(gsub(".*/","",x = x))))
    }
  #return(df_res)
  }
})#->DF
