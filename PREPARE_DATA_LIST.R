source("init_PREPARE_DATA_LIST.R")

#### MODEL CALCULATION ####
as.list(flx_files$sfn_sites) %>%
  furrr::future_map(function(x){
    #load SAPFLUXNET site and aggregation at daily level
    sfn <- read_sfn_data(x, folder = path) %>%
      sfn_metrics(period = "1 day", .funs = list(~ mean(., na.rm = TRUE)),
                  solar = TRUE, interval = "general") %>%
      sapfluxnetr::metrics_tidyfier(metadata = sfn_meta, interval = "general")
    
    if (!is.na(sfn$st_soil_depth %>% unique())){sensors_d = sfn$st_soil_depth %>% unique()} else {sensors_d = 0.5}
    #load FLUXNET and extract simplified variables if there is EC and SFN overlap
    
    flx_file <- flx_files[which(flx_files$sfn_sites == x),"flx_sites"]
    if(!is.na(flx_file)){
      flx <- splashTools::readFluxdata(paste0("DATA/EC/", flx_file$flx_sites), sfn$si_elev %>% unique(), sensors_d = sensors_d)
      
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
      left_join(df_swc%>%
                  group_by(si_code) %>% 
                  mutate(max_swvl = max(swvl2, na.rm = TRUE)), 
                by=c("TIMESTAMP","si_code")) 
    
    #if there is no column add it filled with NA
    sfn <- add_miss_var(sfn)
    
    #use in site swc if present
    if(all(is.na(sfn$swc_shallow))){
      sfn <- sfn %>% 
        mutate(swc_shallow = swvl2,
               is_st_swc_shallow = FALSE)
    }else{
      swc_model <- lm(swc_shallow ~ swvl2,data=sfn)
      sfn <- sfn %>% 
        mutate(is_st_swc_shallow = case_when(is.na(swc_shallow) ~ FALSE,
                                             TRUE ~ TRUE),
               swc_shallow = case_when(!is_st_swc_shallow ~ swc_model$coefficients[1]+swc_model$coefficients[2]*swvl2,
                                       TRUE ~ swc_shallow))
    }
    
    #join soil properties
    soil <- splashTools::getSoilSite(unique(sfn$si_lat), unique(sfn$si_long)) %>% as.list %>% as_tibble
      #if the cell don't have data, get the 9 adjacent cells and do the mean
    if(any(is.na(soil))){
      print("No soil data in SoilGrids cell, obtainning 9 adjacent cells")
      cells <- matrix(c(1,1,1,0,1,-1,0,1,0,0,0,-1,-1,1,-1,0,-1,-1),nrow=9,byrow = TRUE)
      cells %>% 
        split(rep(1:nrow(cells), each = ncol(cells))) %>% 
        purrr::map_df(function(x){
          splashTools::getSoilSite(unique(sfn$si_lat)+x[1]*0.0025, unique(sfn$si_long)+x[2]*0.0025) %>% as.list %>% as_tibble
        }) -> soil_grid
      soil <- soil_grid %>% summarise_all(mean,na.rm = TRUE)
    }
      #if the 9 cells don't have data, get the 9 adjacent cells of the 9 cells and do the mean
    if(any(is.na(soil))){
      print("No soil data in SoilGrids cell, obtainning 27 adjacent cells")
      cells <- matrix(c(1,1,1,0,1,-1,0,1,0,0,0,-1,-1,1,-1,0,-1,-1),nrow=9,byrow = TRUE)
      cells %>% 
        split(rep(1:nrow(cells), each = ncol(cells))) %>% 
        purrr::map_df(function(x){      
          lat <- unique(sfn$si_lat)+x[1]*0.005
          lon <- unique(sfn$si_long)+x[2]*0.005
          cells %>%
            split(rep(1:nrow(cells), each = ncol(cells))) %>% 
            purrr::map_df(function(y){
              splashTools::getSoilSite(lat+y[1]*0.0025, lon+y[2]*0.0025) %>% as.list %>% as_tibble
              }) -> soil_grid
            soil_grid %>% summarise_all(mean,na.rm = TRUE)
        }) -> soil_grid_cell
      soil <- soil_grid_cell %>% summarise_all(mean, na.rm = TRUE)
      }
    
    sfn <- sfn %>% bind_cols(soil) %>% 
      mutate(is_st_clay = case_when(is.na(st_clay_perc)~ FALSE,
                                    !is.na(st_clay_perc)~ TRUE),
             is_st_sand = case_when(is.na(st_sand_perc)~ FALSE,
                                    !is.na(st_sand_perc)~ TRUE),
             is_st_soil_depth = case_when(is.na(st_soil_depth)~ FALSE,
                                          !is.na(st_soil_depth)~ TRUE),
             st_clay_perc = case_when(is.na(st_clay_perc)~ clay,
                                      !is.na(st_clay_perc)~ st_clay_perc),
             st_sand_perc = case_when(is.na(st_sand_perc)~ sand,
                                      !is.na(st_sand_perc)~ st_sand_perc),
             st_soil_depth = case_when(is.na(st_soil_depth)~ depth*100,
                                      !is.na(st_soil_depth)~ st_soil_depth))
    
    #join SW ERA5
    sfn <- sfn %>% 
      left_join(sw_ERA5,by = c("si_code", "TIMESTAMP_daylight")) %>% 
      mutate(sw_ERA5 = sw_ERA5_18/2) # transform into mean daily value (divide by two to account for night value)
    
    #join elevation and set it if NA
    sfn <- sfn %>% 
      left_join(elevation %>% dplyr::select(si_code, elevation),by = c("si_code")) %>% 
      mutate(is_si_elev = case_when(is.na(si_elev) ~ FALSE,
                                    TRUE ~ TRUE),
             si_elev = case_when(is.na(si_elev) ~ elevation,
                                 TRUE ~ si_elev)) %>% 
      dplyr::select(-elevation)
    
    #join tree height
    sfn <- sfn %>% 
      left_join(simard,by = c("si_code")) %>% 
      left_join(GEDI,by = c("si_code")) %>% 
      mutate(height = case_when((is.na(pl_height) & !is.na(st_height)) ~ st_height,
                                (is.na(pl_height) & is.na(st_height) & !is.na(height_GEDI)) ~ height_GEDI,
                                (is.na(pl_height) & is.na(st_height) & is.na(height_GEDI) & !is.na(height_simard)) ~ height_simard,
                                !is.na(pl_height) ~ pl_height))

    
    ##CALIBRATE SOIL POTENTIAL by adding SWC
    ## THE OPTIMUM SWC value to add
    psi_data_site <-  PSI_DF[PSI_DF$id_sfn == unique(sfn$si_code) & PSI_DF$time_psi == "pre-dawn",]
    is_psi_data <- any(lubridate::date(psi_data_site$TIMESTAMP) %in% lubridate::date(sfn$TIMESTAMP)) 
    
    if(is_psi_data){
      temp <- optim_swc(sfn, psi_data_site, soil, type ="swc")
      opt_swc <- temp[[1]]
      sfn$swp_corrected <- TRUE
      sfn <- sfn %>% 
        mutate(date = lubridate::date(TIMESTAMP)) %>% 
        left_join(temp[[2]],by="date") %>% 
        rename(psi_predawn = psi) %>% 
        select(-date)
      
    }else{
      opt_swc <- c(1,0)
      sfn$swp_corrected <- FALSE
      sfn$psi_predawn <- NA
    }

    
    #### FAPAR ####
    sfn <- include_fapar_lai(sfn, fapar_noaa, fapar)
    
    #### PHYDRO PREPARATION ####
    PHYDRO_TRUE <- FALSE
    par_plant_std <- sp_params(sfn, SpParams)
    if(!is.null(par_plant_std)){
      if(!is.na(par_plant_std$psi50) & !is.na(par_plant_std$K)){
        PHYDRO_TRUE <- TRUE
      }
    }
    
    #### JOIN ALL and SAVE ####
    
    sfn_list <- LIST(sfn, PHYDRO_TRUE, par_plant_std, opt_swc, soil)
    
    save(sfn_list, file = paste0("DATA/PREPARED_DATA_LIST/",x,".RData"))
    
  })
