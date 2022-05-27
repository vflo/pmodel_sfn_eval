#' Summarise environmental data
#'
#' Summarise environmental data at the desired time aggregation level.
#'
#' @param sfn_object tibble object, a sfn_object as the output of sapfluxnetr::metrics_tidifier() with a period of 1 hour
#' @param time_aggregation character, time aggregation as "1 hour", "3 hour", "1 day", "1 week", "1 month"
#'
#' @return tibble, time aggregated environmental variables as mean values. Note that even precipitation is the mean hourly value of the aggregation period!!
#'
#' @importFrom (tidyverse, lubridate, tsibble)
#'
#' @export
#'
calc_sfn_aggr_env <- function(sfn_object, time_aggregation = "1 hour"){
  #check if sfn_object is a tibble. If not, stop calculation.
  if(!is_tibble(sfn_object)){stop("sfn_object should be a tibble")}
  #check if time_aggregation is a character string. If not, stop calculation.
  if(!is.character(time_aggregation)){stop("time_aggregation should be a character string")}
  #check if the aggregation time value is valid. If not, stop calculation.
  valid_time <- c("1 hour", "3 hour", "1 day", "1 week", "1 month")
  if(!time_aggregation %in% valid_time){
    stop(paste0("The time aggregation '", time_aggregation, 
                "' is not a valid entry."))}
  
  env_data <- sfn_object %>% 
    mutate(TIMESTAMP = lubridate::as_date(TIMESTAMP))%>%
    dplyr::select(TIMESTAMP,ta,rh,ppfd_in,sw_in,vpd,ext_rad,netrad,ws,swc_shallow,ext_rad,
                  swvl1:swvl4,st_soil_depth,st_sand_perc,st_clay_perc,OM,gravel,bd,depth,
                  si_elev,FAPAR,Fapar,LAI,CO2,aet, netr, VPD, PPFD, tc,sw_ERA5) %>% 
    unique() %>% 
    group_by(TIMESTAMP) %>%
    summarise_all(mean) %>% 
    as_tsibble(index = TIMESTAMP)
  
  
  if(time_aggregation == "3 hour" ){
    env_data <- env_data %>% 
      group_by_key() %>%
      index_by(timestamp_aggr = ~ lubridate::round_date(., "3 hour")) %>%
      summarise_all(~mean(.,na.rm=TRUE))
  }
  
  if(time_aggregation == "1 day" ){
    env_data <- env_data %>% 
      group_by_key() %>%
      index_by(timestamp_aggr = ~ lubridate::as_date(.)) %>%
      summarise_all(~mean(.,na.rm=TRUE))
  }
  
  if(time_aggregation == "1 week" ){
    env_data <- env_data %>% 
      group_by_key() %>%
      index_by(timestamp_aggr = ~ yearweek(.)) %>%
      summarise_all(~mean(.,na.rm=TRUE))
  }
  
  if(time_aggregation == "1 month" ){
    env_data <- env_data %>% 
      group_by_key() %>%
      index_by(timestamp_aggr = ~ yearmonth(.)) %>%
      summarise_all(~mean(.,na.rm=TRUE))
  }
  
  return(env_data)
  
}
