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
    
    #subset 
    df <- df[c(10:100),]
    # df$swvl <- df$swc_shallow
    # df <- df %>% filter(vpd>=1)
    
    #### MODEL CALCULATION ####
    ## PMODEL ##
    print("Calculating PMODEL")
    pmodel <- calc_pmodel(df, soil, par_plant_std)%>% as_tibble()
    
    # PMODEL SWC limitation if there is AET/PET then calculate it
    print("Calculating PMODEL with swc limitation")
    pmodel_swc <- calc_pmodel_swc(df, meanalpha, soil, par_plant_std)%>% as_tibble()
    
    # PHydro if there is species information 
    print("Calculating PHYDRO")
    sensitivity[6,]  %>% 
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        calc_phydro(df, PHYDRO_TRUE, par_plant_std, soil, x) %>% as_tibble()
      })  %>%  bind_rows() -> phydro 
    
    # Sperry model
    print("Calculating SPERRY model")
    sensitivity[6,] %>% 
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        calc_sperry(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
      }) %>%  bind_rows() -> sperry
    
    # Wang model
    print("Calculating WANG model")
    sensitivity[6,]%>% 
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        calc_wang(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
      }) %>% bind_rows() -> wang
    
    # Wap model
    print("Calculating WAP model")
    sensitivity[6,] %>% 
      split(seq(nrow(.)))%>%
      purrr::map(function(x){
        calc_wap(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
      }) %>%  bind_rows() -> wap
    
    ## PMODEL Ecrit ##
    # print("Calculating PMODEL Ecrit")
    # sensitivity[6,] %>% 
    #   split(seq(nrow(.)))%>%
    #   purrr::map(function(x){
    #     calc_pmodel_ecrit(df, PHYDRO_TRUE, par_plant_std, soil,x)%>% as_tibble()
    #   }) %>% bind_rows() -> pmodel_ecrit
    
    
    df_res <- pmodel %>% 
      bind_rows(pmodel_swc) %>% 
      bind_rows(phydro) %>% 
      bind_rows(sperry) %>% 
      bind_rows(wang) %>% 
      bind_rows(wap) %>%
      # bind_rows(pmodel_ecrit) %>%
      suppressMessages()
    # 
    
    df_res$grp <- format(df_res$TIMESTAMP, "%Y")
    df_res %>%
      filter(low_swp == FALSE) %>%
      # filter(vpd>0.5) %>% # dplyr::select(p_leaf, psi_soil,model_type, gs) %>% View()
      mutate(model_type = as.factor(model_type)) %>% 
      ggplot()+
      # geom_point(aes(log(-psi_soil),gs,color=model_type), shape = 21, alpha = 0.5)+
      # geom_point(aes(log(-psi_soil),Gs_sapflow),fill = "grey40",color = "grey40", shape = 21, alpha = 0.5)+
      # geom_smooth(aes(log(-psi_soil),gs,color=model_type),method="gam",se = FALSE)+
      # geom_smooth(aes(log(-psi_soil),Gs_sapflow),color = "grey40",method="lm", formula = 'y~poly(x,4)')+
      # geom_point(aes(log(-psi_soil),E,color=model_type), shape = 21, alpha = 0.5)+
      # geom_point(aes(log(-psi_soil),E_sapflow),fill = "grey40",color = "grey40", shape = 21, alpha = 0.5)+
      # geom_smooth(aes(log(-psi_soil),E,color=model_type),method="gam",se = FALSE)+
      # geom_smooth(aes(log(-psi_soil),E_sapflow),color = "grey40",method="lm", formula = 'y~poly(x,3)')+
      # geom_point(aes(psi_soil, p_leaf,color=model_type), shape = 21, alpha = 0.5)+
      # geom_smooth(aes(psi_soil,p_leaf,color=model_type),method="gam",se = FALSE)+
      # geom_abline(intercept = 0, slope = 1, linetype = 2)+
      geom_ribbon(aes(x=TIMESTAMP, ymin = (Gs_sapflow - Gs_sapflow_sd), ymax = (Gs_sapflow + Gs_sapflow_sd),
                      group = grp), fill = "grey" )+
      geom_line(aes(TIMESTAMP,Gs_sapflow,group = grp),color = "black")+
      geom_line(aes(TIMESTAMP,gs, color=model_type,group = interaction(grp,model_type)))+
      # geom_line(aes(TIMESTAMP,LAI,group = grp), color = "green")+
      # geom_ribbon(aes(x=TIMESTAMP, ymin = (E_sapflow -E_sapflow_sd), ymax = (E_sapflow + E_sapflow_sd),
      #                 group = grp), fill = "grey" )+
      # geom_line(aes(TIMESTAMP,E_sapflow,group = grp),color = "black")+
      # geom_line(aes(TIMESTAMP,E, color=model_type,group = interaction(grp,model_type)))+
      theme_bw()+
      # ylim(-4,0)+
      NULL
    
    #### SAVE DATA ####
    if(nrow(df_res)>0){
      save(df_res, file = paste0("DATA/OUTCOME_WEEKLY/",(gsub(".*/","",x = x))))
    }
    
    
    # 
    # #### PMODEL MONTHLY ####
    #   df <- data_prep(sfn_aggr_month, env_month, opt_swc)
    # 
    #   if(nrow(df)>1 && any(!is.na(df$E_stand))){
    # 
    # #### MODEL CALCULATION ####
    #     ## PMODEL ##
    #     pmodel <- calc_pmodel(df, soil)%>% as_tibble()
    #     
    #     # PMODEL SWC limitation if there is AET/PET then calculate it
    #     pmodel_swc <- calc_pmodel_swc(df, meanalpha, soil)%>% as_tibble()
    #     
    #     # PHydro if there is species information 
    #     sensitivity %>% 
    #       split(seq(nrow(.)))%>%
    #       purrr::map(function(x){
    #         calc_phydro(df, PHYDRO_TRUE, par_plant_std, soil, x) %>% as_tibble()
    #       })  %>%  bind_rows() -> phydro 
    #     
    #     # Sperry model
    #     sensitivity %>% 
    #       split(seq(nrow(.)))%>%
    #       purrr::map(function(x){
    #         calc_sperry(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
    #       }) %>%  bind_rows() -> sperry
    #     
    #     # Wang model
    #     sensitivity %>% 
    #       split(seq(nrow(.)))%>%
    #       purrr::map(function(x){
    #         calc_wang(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
    #       })  %>%  bind_rows() -> wang
    #     
    #     # Wap model
    #     sensitivity %>% 
    #       split(seq(nrow(.)))%>%
    #       purrr::map(function(x){
    #         calc_wap(df, PHYDRO_TRUE, par_plant_std, soil,x) %>% as_tibble()
    #       })  %>%  bind_rows() -> wap
    # 
    #     df_res <- pmodel %>% 
    #       bind_rows(pmodel_swc) %>% 
    #       bind_rows(phydro) %>%  
    #       bind_rows(sperry) %>% 
    #       bind_rows(wang) %>% 
    #       bind_rows(wap) %>%
    #       suppressMessages()
    # 
    #   if(nrow(df_res)>0){
    #     save(df_res, file = paste0("DATA/OUTCOME_MONTHLY/",(gsub(".*/","",x = x))))
    #   }
    # #return(df_res)
    # }
  })#->DF
