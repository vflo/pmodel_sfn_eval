#' Summarise transpiration data
#'
#' Summarise transpiration data at the desired time aggregation level.
#'
#' @param sfn_scale tibble object, an object as the output of calc_sfn_scale()
#' @param time_aggregation character, time aggregation as "1 hour", "3 hour"
#'
#' @return tibble, time aggregated Ttranspiration as mean values.
#'
#' @importFrom (tidyverse, lubridate, tsibble)
#'
#' @export
#'
calc_sfn_aggr_E_hourly <- function(sfn_scale, time_aggregation = "1 hour"){
  #check if sfn_object is a tibble. If not, stop calculation.
  if(!is_tibble(sfn_scale)){stop("sfn_scale should be a tibble object as the output of calc_sfn_scale function")}
  #check if time_aggregation is a character string. If not, stop calculation.
  if(!is.character(time_aggregation)){stop("time_aggregation should be a character string")}
  #check if the aggregation time value is valid. If not, stop calculation.
  valid_time <- c("1 hour", "3 hour")
  if(!time_aggregation %in% valid_time){
    stop(paste0("The time aggregation '", time_aggregation, 
                "' is not a valid entry."))}
  
  E_data <- sfn_scale %>% 
    mutate(TIMESTAMP = lubridate::ymd_hms(TIMESTAMP)) %>%
    # unique() %>% 
    # mutate(key  = row_number()) %>%
    as_tsibble(index = TIMESTAMP)# %>%
    # dplyr::select(-key)
  
  
  if(time_aggregation == "3 hour" ){
    E_data <- E_data %>% 
      group_by_key() %>%
      index_by(timestamp_aggr = ~ lubridate::round_date(., "3 hour")) %>%
      summarise_all(~mean(.,na.rm=TRUE))
  }

  return(E_data)
  
}
