do_phydro_params <- function(sp, sfn = sfn, visco = visco, density = density){
  v25 = sp$v25
  j25=sp$j25
  if(is.na(v25)){v25 = 50}
  if(is.na(j25)){j25 = 100}
  cond_Kl <- sp$Kl
  c <- sp$c
  d <- sp$d
  psi <- seq(0,-9,-0.01)
  cond_Ks <- sp$Ks
  if(is.na(cond_Ks)){cond_Ks <- NA}
  if(is.na(cond_Kl)){cond_Kl <- NA}
  if(is.na(c)){c <- NA}
  if(is.na(d)){d <- NA}
  vulne_curve <- cond_Ks*exp(-(psi/d)^c)
  
  opt_curve_param <- function(par){
    p <- par[1]
    b <- par[2]
    res <- (1/2)^((psi/p)^b)*cond_Ks
    rmse <- sqrt(sum((vulne_curve-res)^2)/length(vulne_curve))
    return(rmse)
  }
  par <- tibble(par = c(NA,NA))
  if(any(!is.na(vulne_curve))){par <- optim(c(-3,2), opt_curve_param)}
  
  #res <- (1/2)^((psi/par$par[1])^par$par[2])
  #plot(vulne_curve)
  #lines(cond*res)
  sfn %>% 
    filter(pl_species == sp$pl_species) %>% 
    group_by(pl_code) %>%
    dplyr::select(pl_code, pl_sapw_area, pl_leaf_area, pl_height, st_lai, st_basal_area, sp_basal_area_perc, st_height) %>% 
    summarise_all(unique,
                  .groups = "drop") %>% 
    dplyr::select(-pl_code) %>% 
    summarise_all(mean, na.rm=TRUE) -> plant_sp_data
  
  par_plant_std = list(
    #Ks0=1e-12, # m2
    st_basal_area = plant_sp_data$st_basal_area,
    sp_basal_area = plant_sp_data$sp_basal_area_perc,
    v_huber = plant_sp_data$pl_sapw_area/plant_sp_data$pl_leaf_area/10^4, # m2sapwood m-2leaf
    height = plant_sp_data$pl_height, # m
    LAI = plant_sp_data$st_lai, # m2leaf m-2soil
    conductivity = cond_Ks, #Kg m-1 MPa-1 s-1
    conductivity_Kl = cond_Kl, #Kg m-1 MPa-1 s-1
    # conductivity = cond*visco/(1e9*density*55.5), # mmol s-1 m-2 MPa-1; it should be provided as m3/m2/s/Pa: to back transform -> kmax*viscosity_water/(1e3*1e6*density*55.5)
    #conductivity_scalar=1,  # Scales (Ks0*v_huber/H)
    psi50 = par$par[1], # MPa
    b=par$par[2],
    c = c,
    d = d,
    v25 = v25,
    j25 = j25
  )
  
  return(par_plant_std)
}
