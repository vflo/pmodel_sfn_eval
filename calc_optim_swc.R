#' Optimization of swc addition
#'
#' Calculate the optim swc to be added to the swc measured at the site to match estimated soil water potential using "Saxton" with actual soil water potential
#'
#' @param sfn_object tibble object, a sfn_object as the output of sapfluxnetr::metrics_tidifier()
#' @param psi_data tibble object, using PSI SAPFLUXNET data base format.
#' @param min_swc numeric, the minimum value of swc that can be optimised.
#' @param max_swc numeric, the maximum value of swc that can be optimised. Should be less than 1.
#' @param type character, "swc" or "swvl" depending on the soil water content source.
#'
#' @return numeric, optimal swc value to add to the site's swc
#'
#' @references  
#'
#' @export
#'
optim_swc <- function(sfn_object, psi_data, soil=soil, min_swc = 0, max_swc = 0.6, type = "swvl"){
  #check if sfn_object is a tibble. If not, stop calculation.
  if(!is_tibble(sfn_object)){stop("sfn_object should be a tibble")}
  #check if psi_data is a tibble. If not, stop calculation.
  if(!is_tibble(psi_data)){stop("psi_data should be a tibble")}
  #check if type is a valid value. If not, stop calculation.
  if(!type %in% c("swvl","swc")){stop("type variable must be a valid value")}
  
  #extract the site code
  si_code <- sfn_object$si_code %>% unique()
  
  #extract predawn psi data of the site 
  psi <- psi_data %>% 
    filter(time_psi == 'pre-dawn')
  
  #check if there is data of the site
  if(!any(unique(psi$id_sfn) == si_code)){
    df <- NA
    print("There is no psi data for this site")
  }else{
    #psi data preparation
    psi <- psi %>% 
      mutate(date = lubridate::date(TIMESTAMP)) %>% 
      rename(n_psi = 'psi_N') %>% 
      dplyr::select(date,psi,n_psi) %>% 
      group_by(date) %>% 
      summarise(psi = weighted.mean(psi,n_psi)) # summarize actual soil water potential as the weighted mean of each day (the number of measures are the weights)
    
    #sf_object preparation
    if(type == "swvl"){
      sfn_swvl <- sfn_object %>%
        mutate(hour = lubridate::hour(TIMESTAMP)) %>%
        # filter(hour == 5) %>%
        dplyr::select(TIMESTAMP, st_soil_depth,st_sand_perc,st_clay_perc,swvl1:swvl4) %>%
        unique() %>%
        mutate(swvl = swvl_aggregation(st_soil_depth,swvl1,swvl2,swvl3,swvl4)) %>% 
        group_by(TIMESTAMP,st_soil_depth,st_sand_perc,st_clay_perc) %>%
        summarise(swvl = swvl, .groups ="drop" ) %>% 
        # ungroup() %>% 
        filter(!is.na(swvl)) 
    }
    
    if(type == "swc"){
      sfn_swvl <- sfn_object %>%
        mutate(hour = lubridate::hour(TIMESTAMP)) %>%
        # filter(hour == 5) %>%
        dplyr::select(TIMESTAMP, st_soil_depth,st_sand_perc,st_clay_perc,swvl1:swvl4,swc_shallow) %>%
        unique() %>%
        mutate(swvl = swvl_aggregation(st_soil_depth,swvl1,swvl2,swvl3,swvl4)) %>% 
        group_by(TIMESTAMP,st_soil_depth,st_sand_perc,st_clay_perc) %>%
        summarise(swvl = swc_shallow,.groups = "drop" ) %>% 
        # ungroup() %>% 
        filter(!is.na(swvl)) 
    }

    
    # objective function
    swc_increment <- function(par){
      b <- par[1]
      s_x <- par[2]
      psi_sfn <- sfn_swvl %>% 
        mutate(psi_sfn = medfate::soil_psi(medfate::soil(tibble(widths = st_soil_depth*10, clay = st_clay_perc, 
                                                                sand = st_sand_perc, om = soil$OM, bd = soil$bd, rfc = 70),
                                                         W = (b * swvl+s_x)/
                                                           soil_thetaSATSX(clay = st_clay_perc, sand = st_sand_perc, om = soil$OM), model = "SX"))) %>%
        dplyr::select(TIMESTAMP, psi_sfn) %>%
        mutate(date = lubridate::date(TIMESTAMP))
    
    df_psi <- psi %>% left_join(psi_sfn, by="date") %>% na.omit()
    
    rmse <- sqrt(sum((df_psi$psi_sfn-df_psi$psi)^2)/nrow(df_psi))
    if(any((b * sfn_swvl$swvl+s_x)>=max_swc)){rmse <- 1e3}
    if(any((b * sfn_swvl$swvl+s_x)<min_swc)){rmse <- 1e3}
    return(rmse)
    }
    
  df <- optim(c(1.1,0.1), swc_increment)$par
  }
  return(df)
}
