

###########################################
## HYDRAULIC FUNCTIONS
###########################################
huber_value <- function(par_env, par_plant){
  par_plant$pl_sapw_area / par_plant$pl_ba * par_plant$st_basal_area * 1e-4 / par_env$LAI
}

## Returns conductivity in mol/m2/s/Mpa
scale_conductivity_ks = function(K, par_env, par_plant, do_backtransform = TRUE){
  if(do_backtransform){
    #from Kg m-1 s-1 MPa-1 to m3 m-2
    K = (K*par_env$viscosity_water*huber_value(par_env, par_plant))/(par_env$density_water*1e6*par_plant$height) #back transform to m3 m-2
  }
  
  # Flow rate in m3/m2/s/Pa
  K2 = K / par_env$viscosity_water
  
  # Flow rate in mol/m2/s/Pa
  mol_h20_per_kg_h20 = 55.5
  K3 = K2 * par_env$density_water * mol_h20_per_kg_h20
  
  # Flow rate in mol/m2/s/Mpa
  K4 = K3*1e6
  
  return(K4)  
}

calc_conductivity_m = function(sapwood_perm, hv, height){
  sapwood_perm*hv/height
}

## Returns conductivity in mol/m2/s/Mpa
scale_conductivity = function(K, par_env){
  # Flow rate in m3/m2/s/Pa
  K2 = K / par_env$viscosity_water
  
  # Flow rate in mol/m2/s/Pa
  mol_h20_per_kg_h20 = 55.5
  K3 = K2 * par_env$density_water * mol_h20_per_kg_h20
  
  # Flow rate in mol/m2/s/Mpa
  K4 = K3*1e6
  
  return(K4)  
}

calc_gs_phydro <- function (dpsi, psi_soil, par_plant, par_env, ...){
  K = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (leaf) s-1 MPa-1
  # K = K * par_plant$pl_sapw_area/par_plant$pl_ba*par_plant$st_basal_area*1e-4 #kg m-2 (ground) s-1 MPa-1
  # K = K * par_env$LAI #kg m-2 (leaf) s-1 MPa-1
  D = (par_env$vpd/par_env$patm)
  # D = par_env$vpd
  K/1.6/D * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (leaf) s-1
}

calc_gs <- function (dpsi, psi_soil, par_plant, par_env, ...){
  K = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE) #mol m-2 (leaf) s-1 MPa-1
  # K = K * par_plant$pl_sapw_area/par_plant$pl_ba*par_plant$st_basal_area*1e-4 #kg m-2 (ground) s-1 MPa-1
  # K = K * par_env$LAI #kg m-2 (leaf) s-1 MPa-1
  D = (par_env$vpd/par_env$patm)
  # D = par_env$vpd
  K/1.6/D * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (leaf) s-1
}

integral_P <- function(dpsi, psi_soil, psi50, b, ...){
  ps = psi_soil/psi50
  pl = (psi_soil-dpsi)/psi50
  l2 = log(2)
  pl_ = l2*pl^b
  ps_ = l2*ps^b
  -(psi50/b)*(l2^(-1/b))*(expint::gammainc(a = 1/b, x = pl_)-expint::gammainc(a = 1/b, x = ps_))
}

integral_P_e_ecrit <- function(dpsi, psi_soil, psi50, b, ...){
  ps = psi_soil/psi50
  pl = (psi_soil-dpsi)/psi50
  l2 = log(2)
  pl_ = l2*pl^b
  ps_ = l2*ps^b
  -1*((expint::gammainc(a = 1/b, x = pl_)-expint::gammainc(a = 1/b, x = ps_)))/
    (-1*expint::gammainc(a = 1/b, x = ps_))
}

integral_P_ecrit <- function(psi_soil, psi50, b, ...){
  ps = psi_soil/psi50
  l2 = log(2)
  ps_ = l2*ps^b
  (psi50/b)*(l2^(-1/b))*(expint::gammainc(a = 1/b, x = ps_))
}

integral_P_soil <- function(psi_soil, psi50, b, ...){
  ps = psi_soil/psi50
  l2 = log(2)
  ps_ = l2*ps^b
  
  expint::gammainc(a = 1/b, x = ps_)/expint::gammainc(a = 1/b, x = 0)
}

integral_P_leaf <- function(psi_soil, psi50, b, ...){
  pl = (psi_soil-dpsi)/psi50
  l2 = log(2)
  pl_ = l2*pl^b
  expint::gammainc(a = 1/b, x = pl_)/expint::gammainc(a = 1/b, x = 0)
}

get_p <- function(psi_soil, e, K, d, c, h, dens_water) {
  dp = 0.0
  psi_soil = -psi_soil
  tension = psi_soil
  # 3.2 iterate through each layer to get stem dp including gravity
  for (i in 1:50){
    p = tension
    f = exp( -1.0 * (p/-d)^c )
    layer_k = K *50* f
    layer_k = max(layer_k, 1E-12)
    dp <-  dp + e / layer_k + dens_water*9.8*h*0.02*1e-6
    tension = psi_soil + dp
  }
  
  return(-tension)
}


get_e_crit <- function(psi_soil, K, d, c, h, dens_water,precision=6){
  p_crit = d * log(1000.0) ** (1.0/c)
  e_min  = 0.0
  e_max  = 100.0
  e_crit = 50.0
  continue <- TRUE
  while(continue){
    p = get_p(psi_soil,e_max, K, d, c, h, dens_water)
    if (p>p_crit){
      e_max <-  e_max*2
    }else{
      continue <- FALSE
    }
  }
  continue <- TRUE
  while(continue){
    e = 0.5 * (e_max+e_min)
    p = get_p(psi_soil,e, K, d, c, h, dens_water)
    if (abs(p-p_crit)<10^(-1*precision) | (e_max-e_min)<10^(-1*precision)){
      e_crit = e
      continue <- FALSE
    }
    if (p<p_crit){ 
      e_max  = e
    }else{
      e_min  = e
    }
    
  }
  return (e_crit)
}
