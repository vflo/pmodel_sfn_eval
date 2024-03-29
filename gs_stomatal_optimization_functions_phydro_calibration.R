###########################################
## PHYDRO SCHEMES
###########################################
fn_profit_gs <- function(par, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                      par_plant, par_env, do_optim = FALSE, stomatal_model){
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  gs = par[2]#      # delta Psi in MPa
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  dpsi    = calc_dpsi_phydro(gs, psi_soil, par_plant, par_env) # dpsi in MPa
  psi_leaf = psi_soil-dpsi  # leaf water potential
  # if(jmax<=0|jmax>1e5|
  #    gs > 10|gs<=0|
  #    psi_leaf<p_crit){
  #   out <- -1e20
  #   }else{
  ## light-limited assimilation
    a_j   <- calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
    a     = a_j$a
    ci    = a_j$ci
    vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
  
  ## dummy cost
  dummy_costs = 0*exp(20*(-abs(dpsi/4)-abs(jmax/1))) # ONLY added near (0,0) for numerical stability.
  
  ## phydro
  if(stomatal_model == "phydro"){
    costs = par_cost$alpha * jmax + par_cost$gamma * dpsi^2 #((abs((-dpsi)/par_plant$psi50)))^2  
    out <- a - costs - dummy_costs
  }
  
  ## wang
  if(stomatal_model == "phydro_wang"){
    cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    if(is.na(e_e_crit)){e_e_crit = 1}
    out_e = (a- cost) * e_e_crit
    out = exp((a - cost) - out_e)
  }
  
  ## wang modificated
  if(stomatal_model == "phydro_wang_mod"){
    cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    if(is.na(e_e_crit)){e_e_crit = 1}
    out_e = a * e_e_crit
    out = exp(a - out_e - cost)
  }
  
  # Least-cost
  if(stomatal_model == "phydro_least_cost"){
    cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
    out = exp((a-cost)/(par_cost$gamma*e+jmax))
  }
  
  #CGain
  if(stomatal_model == "phydro_cgain"){
    K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl     = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    ks     = K*(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    out = exp(a - par_cost$alpha*jmax - par_cost$gamma*(K-kl)/K)
  }
  
  # WUE
  if(stomatal_model == "phydro_wue"){
    out = exp(a - par_cost$alpha*jmax - par_cost$gamma*e)
  }
  
  ## cmax
  if(stomatal_model == "phydro_cmax"){
    aa         = par_cost$gamma
    bb         = 1
    p = psi_leaf
    out <- exp((a - aa*p**2.0 - bb*p) - (par_cost$alpha * jmax))
  }
  
  ## sperry
  if(stomatal_model == "phydro_sperry"){
    K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    g   <- e_crit/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    amax <-  c_a$a
    ks     = K *(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    kl     = K *(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    k_crit = 0
    risk   = amax*(1-kl/ks)
    cost   = (par_cost$alpha * jmax)
    # amax   = max(unlist(list_a))
    out <- exp(a - risk - cost)
  }
  
  ## SOX
  if(stomatal_model == "phydro_sox"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    out <- exp((a- (par_cost$alpha * jmax))*reduction_factor )
  }
  
  ## SOX mod
  if(stomatal_model == "phydro_sox_mod"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    out <- exp((a)*reduction_factor - (par_cost$alpha * jmax))
    # if (a<(par_cost$alpha * jmax)){out <- 0}
  }
  # }
  if (do_optim){
    return(-out)
  } else {
    return(out)
  }
}

fn_profit_inst_schemes_gs <- function(par, jmax, vcmax, psi_soil, e_crit, p_crit, par_cost, 
                                   par_photosynth, par_plant, par_env, 
                                   stomatal_model, do_optim){
  gs = par[1]#      # gs in mol m-2 (ground) s-1
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)$a
  dpsi    = calc_dpsi_phydro(gs, psi_soil, par_plant, par_env) # dpsi in MPa
  psi_leaf = psi_soil-dpsi  # leaf water potential

  
  ## phydro
  if(stomatal_model == "phydro"){
    profit = A - par_cost$gamma * dpsi^2 
  }
  
  ## wang
  if(stomatal_model == "phydro_wang"){
    cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    if(is.na(e_e_crit)){e_e_crit = 1}
    out_e = (A- cost) * e_e_crit
    profit = exp((A - cost) - out_e)
  }
  
  ## wang modificated
  if(stomatal_model == "phydro_wang_mod"){
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    if(is.na(e_e_crit)){e_e_crit = 1}
    profit = A*(1-e_e_crit)
    profit = exp(profit)
  }
  
  ## Least-cost #Doesn't work because of jmax in A and as cost.
  if(stomatal_model == "phydro_least_cost"){
    profit = exp(A/(par_cost$gamma*e+jmax))
  }
  
  ## CGain
  if(stomatal_model == "phydro_cgain"){
    K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl     = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    profit = exp(A - par_cost$gamma*(K-kl)/K)
  }
  
  ## WUE
  if(stomatal_model == "phydro_wue"){
    profit = exp(A - (par_cost$gamma*e))
  }
  
  ## cmax
  if(stomatal_model == "phydro_cmax"){
    aa         = par_cost$gamma
    bb         = 1
    p = psi_leaf
    profit <- exp(A - aa*p**2.0 - bb*p)
  }
  
  ## sperry
  if(stomatal_model == "phydro_sperry"){
    K          = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (ground) s-1 MPa-1
    g   <- e_crit/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    amax <-  c_a$a
    ks     = K *(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    kl     = K *(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    k_crit = 0
    risk   = amax*(1-kl/ks)
    # amax   = max(unlist(list_a))
    profit <- A - risk
    profit = exp(profit)
  }
  
  ## SOX
  if(stomatal_model == "phydro_sox"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    profit <- exp((A - (par_cost$alpha * jmax))*reduction_factor )
  }
  
  ## SOX mod
  if(stomatal_model == "phydro_sox_mod"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    profit <- (A)*reduction_factor
    profit = exp(profit)
  }
  
  if (do_optim){
    return(-profit)
  } else {
    return(profit)
  }
}


optimise_stomata_phydro_schemes_gs <- function(fn_profit_gs, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                                            par_plant, par_env, jmax_lim, return_all = FALSE, do_optim=TRUE, stomatal_model){
  if(is.null(jmax_lim)){jmax_lim = 7}
  
  jmax_ini = 0
  gs_ini = 1e-1
  continue = TRUE
  while(continue){
    out_optim <- optimr::optimr(
      par       = c(logjmax=jmax_ini, gs=gs_ini),  
      lower     = c(-10, 1e-10),
      upper     = c(jmax_lim, 10),
      fn             = fn_profit_gs,
      psi_soil       = psi_soil,
      e_crit         = e_crit,
      p_crit         = p_crit,
      par_cost       = par_cost,
      par_photosynth = par_photosynth,
      par_plant      = par_plant,
      par_env        = par_env,
      do_optim       = do_optim, 
      stomatal_model = stomatal_model,
      method         = "L-BFGS-B",
      control        = list(maxit = 500, maximize = TRUE, fnscale = 1e2)
    )
    
    gs_prov <- out_optim$par[2]
    e  = 1.6*gs_prov*(par_env$vpd/par_env$patm)
    dpsi    = calc_dpsi_phydro(gs_prov, psi_soil, par_plant, par_env) # dpsi in MPa
    psi_leaf = psi_soil-dpsi  # leaf water potential
     
    if(all(e>=e_crit, (psi_leaf)>=p_crit)){continue = FALSE
    }else{
      gs_ini = gs_ini/10
    }
    if(gs_ini < 1e-10){continue = FALSE}
  }
  out_optim$value <- -out_optim$value
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}

optimise_shortterm_schemes_gs <- function(fn_profit_inst_gs, jmax, vcmax, psi_soil, e_crit, p_crit,
                               par_cost, par_photosynth, par_plant, par_env, 
                               stomatal_model, return_all = FALSE, do_optim){
  gs_ini = 1e-1
  continue = TRUE
  while(continue){
    out_optim <- optimr::optimr(
      par       = c(gs=gs_ini),  
      lower     = c(1e-10),
      upper     = c(10),
      fn        = fn_profit_inst_schemes_gs,
      psi_soil  = psi_soil,
      jmax      = jmax,
      vcmax     = vcmax,
      e_crit    = e_crit,
      p_crit    = p_crit,
      par_cost  = par_cost,
      par_photosynth = par_photosynth,
      par_plant = par_plant,
      par_env   = par_env,
      do_optim  = do_optim, 
      stomatal_model = stomatal_model,
      method    = "L-BFGS-B",
      control   = list(maxit = 500, maximize = TRUE) 
    )
    
    gs_prov <- out_optim$par[1]
    e  = 1.6*gs_prov*(par_env$vpd/par_env$patm)
    dpsi    = calc_dpsi_phydro(gs_prov, psi_soil, par_plant, par_env) # dpsi in MPa
    psi_leaf = psi_soil-dpsi  # leaf water potential
    
    if(all(e>=e_crit, (psi_leaf)>=p_crit)){continue = FALSE}else{gs_ini = gs_ini/10}
    if(gs_ini < 1e-10){continue = FALSE}
  }
  
  out_optim$value <- -out_optim$value
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}


###########################################
## MODEL OPTIMIZATION AND DATA PREPARATION
###########################################
model_numerical_gs <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, 
                            rdark, par_plant, par_cost, stomatal_model){

  patm = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm       = rpmodel::calc_kmm(tc, patm),
    gammastar = rpmodel::calc_gammastar(tc, patm),
    phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
    Iabs      = ppfd * fapar,
    ca        = co2 * patm * 1e-6,  # Convert to partial pressure
    patm      = patm,
    delta     = rdark
    # delta     = ftemp_inst_rd( tc )
  )
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    patm            = patm,
    tc              = tc,
    vpd             = vpd
  )
  

    K      = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (ground) s-1 MPa-1
    p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
    e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1

    # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit 
    # low_swp <- FALSE
    # if(psi_soil <= p_crit){
    #   psi_soil = p_crit*0.95
    #   low_swp  = TRUE
    #   print("WARNING: soil water potential is lower than critical plant water potential.")
    # }

    # 3. Optimizer
    lj_gs = optimise_stomata_phydro_schemes_gs(fn_profit_gs, 
                                             psi_soil = psi_soil,
                                             e_crit = e_crit,
                                             p_crit = p_crit,
                                             par_cost  = par_cost, 
                                             par_photosynth = par_photosynth, 
                                             par_plant = par_plant, 
                                             par_env = par_env,
                                             jmax_lim = 7, 
                                             return_all = FALSE, 
                                             do_optim=TRUE, 
                                             stomatal_model = stomatal_model)
  
    
    jmax  = exp(lj_gs[1]) %>% unname()
    gs  = lj_gs[2] %>% unname()

    dpsi    = calc_dpsi_phydro(gs, psi_soil, par_plant, par_env) # dpsi in MPa
    psi_l = psi_soil-dpsi
    a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
    a     = a_j$a
    ci    = a_j$ci
    vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
    E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
    gs    = gs #transform to umol m-2(ground) s-1 Pa-1

    return(list(
      jmax         = jmax,
      dpsi         = dpsi,
      p_leaf       = psi_l,
      gs           = gs,
      E            = E,
      a            = a,
      ci           = ci,
      chi          = ci/par_photosynth$ca,
      vcmax        = vcmax,
      chi_jmax_lim = 0
    ))
}


########################################################
## INSTANTANEOUS MODEL OPTIMIZATION AND DATA PREPARATION
########################################################
model_numerical_instantaneous_gs <- function(vcmax, jmax, tc, ppfd, vpd, co2, elv, fapar, 
                                          kphio, psi_soil,  rdark, par_plant, 
                                          par_cost, stomatal_model){
  
  patm = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm       = rpmodel::calc_kmm(tc, patm),
    gammastar = rpmodel::calc_gammastar(tc, patm),
    phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
    Iabs      = ppfd * fapar,
    ca        = co2 * patm * 1e-6,  # Convert to partial pressure
    patm      = patm,
    delta     = rdark
    # delta     = ftemp_inst_rd( tc )
  )
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    patm            = patm,
    tc              = tc,
    vpd             = vpd
  )

  K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
  p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
  e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1
    
  
  # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit
  # low_swp <- FALSE
  # if(psi_soil <= p_crit){
  #   psi_soil = p_crit*0.95
  #   low_swp  = TRUE
  #   print("WARNING: soil water potential is lower than critical plant water potential.")
  # }

  # 3. Optimizer
  lj_gs = optimise_shortterm_schemes_gs(fn_profit_inst_schemes_gs,
                                      jmax = jmax, 
                                      vcmax = vcmax,
                                      psi_soil = psi_soil,
                                      e_crit = e_crit,
                                      p_crit= p_crit,
                                      par_cost  = par_cost, 
                                      par_photosynth = par_photosynth, 
                                      par_plant = par_plant, 
                                      par_env = par_env,
                                      return_all = FALSE, 
                                      do_optim = TRUE, 
                                      stomatal_model = stomatal_model)
  
  gs  = lj_gs[1] %>% unname() # delta Psi in MPa

  dpsi    = calc_dpsi_phydro(gs, psi_soil, par_plant, par_env) # dpsi in MPa
  psi_l = psi_soil-dpsi  # leaf water potential
  a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
  a     = a_j$a
  ci    = a_j$ci
  E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
  gs    = gs #transform to umol m-2(ground) s-1 Pa-1
    
  return(list(
      jmax         = jmax,
      dpsi         = dpsi,
      p_leaf       = psi_l,
      gs           = gs,
      E            = E,
      a            = a,
      ci           = ci,
      chi          = ci/par_photosynth$ca,
      vcmax        = vcmax,
      chi_jmax_lim = 0
  ))

}
