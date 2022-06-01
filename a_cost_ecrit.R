get_p <- function(psi_soil, e, K, d, c, h) {
  dp = 0.0
  tension = 0
  # 3.2 iterate through each layer to get stem dp including gravity
  for (i in 1:20){
    p = min(psi_soil, tension)
    f = exp( -1.0 * (p/d)**c )
    layer_k = K *20* f
    layer_k = max(layer_k, 1E-12)
    dp <-  dp + e / layer_k + 998.0*9.8*h*0.05*1E-6
    tension = psi_soil - dp
  }
  
  return(tension)
}

get_e_crit <- function(psi_soil, K, d, c, h){
  p_crit = d * log(1000.0) ** (1.0/c)
  e_min  = 0.0
  e_max  = 100.0
  e_crit = 50.0
  continue <- TRUE
  while(continue){
    p = get_p(psi_soil,e_max, K, d, c, h)
    if (p>p_crit){
      e_max <-  e_max*2
    }else{
      continue <- FALSE}}
  continue <- TRUE
  while(continue){
    e = 0.5 * (e_max+e_min)
    p = get_p(psi_soil,e, K, d, c, h)
    if (abs(p-p_crit)<1E-3 | (e_max-e_min)<1E-3){
      e_crit = e
      continue <- FALSE
    }
    if (p<p_crit){ 
      e_max  = e
    }else{e_min  = e}
    
  }
  return (e_crit)
}


ecrit <- get_e_crit(psi_soil, K, d, c, h)*LAI*1e-3 #micromols??
resp <- rs * rpmodel::calc_ftemp_inst_rd(tc)
a_cost = resp/ecrit