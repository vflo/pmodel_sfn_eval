#' Stand transpiration upscaling
#'
#' Calculate SAPFLUXNET stand upscaling of transpiration
#'
#' @param sfn_object tibble object, a sfn_object as the output of sapfluxnetr::metrics_tidifier()
#' @param min_basal_area numeric, minimum percentage of basal area of the stand represented by measured trees in SAPFLUXNET to allow upscaling, %
#'
#' @return tibble, stand aggregated transpiration
#'
#' @references  Poyatos, R., Granda, V., Flo, V., Molowny-Horas, R., Steppe, K., 
#' Mencuccini, M., & Martínez-Vilalta, J. (2019). SAPFLUXNET: A
#' global database of sap flow measurements. Zenodo. https://doi.
#' org/10.5281/zenodo.2530798
#'
#' @export
#'
calc_sfn_scale <- function(sfn_object, min_basal_area = 70){
  #check if sfn_object is a tibble. If not, stop calculation.
  if(!is_tibble(sfn_object)){stop("sfn_object should be a tibble")}
  #check if stand basal area is greater than min_basal_area. If not, stop calculation.
  ba_perc <- sfn_object$sp_basal_area_perc %>% unique() %>% sum(na.rm = TRUE) #percentage of total basal area of plants measured in SAPFLUXNET
  if(ba_perc < min_basal_area){
    warning(paste0("The basal area of the forest stand is less than ", as.character(min_basal_area), 
                "%, consider reducing the minimum basal area."))}
  
  st_ba <- sfn_object$st_basal_area %>% unique() #stand basal area m2 ha-1
  
  sfn_object %>% 
    group_by(TIMESTAMP,pl_code,pl_species) %>% 
    summarise(pl_basal_area = pi * (pl_dbh/2)^2 * 1e-4, #individual basal area m2
              E_plant_ba = sapflow_ / pl_basal_area, #plant transpiration per unit basal area cm3 m-2basalarea h-1
              sp_basal_area_perc = unique(sp_basal_area_perc),
              .groups = "drop"
              ) %>% 
    group_by(TIMESTAMP, pl_species) %>% 
    summarise(E_sp_ba = mean(E_plant_ba, na.rm = TRUE), #mean species transpiration per unit basal area
              E_sp_ba_variance = var(E_plant_ba, na.rm = TRUE),
              sp_basal_area_perc = unique(sp_basal_area_perc),
              E_sp_perc = E_sp_ba * sp_basal_area_perc/100, #species transpiration normalized by species basal area
              .groups = "drop") %>% 
    group_by(TIMESTAMP) %>% 
    summarise(E_stand = sum(E_sp_perc)*100/ba_perc*st_ba/1e4, #sum of total stand transpiration cm3 m-2(ground) h-1
              E_stand_sd = sqrt(sum(E_sp_ba_variance)/dplyr::n_distinct(pl_species))*100/ba_perc*st_ba/1e4) #sum of total stand transpiration cm3 m-2(ground) h-1
  
}


