###########################################
## ACCLIMATION OPTIMIZATION SCHEMES
###########################################
fn_profit <- function(par, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                      par_plant, par_env, do_optim = FALSE, stomatal_model){
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]#      # delta Psi in MPa
  psi_leaf = psi_soil-dpsi #MPa
  
  # if(stomatal_model %in% c("phydro","phydro_wang","phydro_wang_mod", "phydro_sperry")){
  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)+1e-270   # gs in mol/m2ground/s
  # }else{
  #   K = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (leaf) s-1 MPa-1
  #   D = (par_env$vpd/par_env$patm)
  #   gs = K/1.6/D * dpsi #mol m-2 (leaf) s-1
  # }
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  
  ## light-limited assimilation
  a_j   <- calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
  a     = a_j$a
  ci    = a_j$ci
  vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth) # vcmax based on coordination theory
  
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
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b) #Calculate e over e_crit
    if(is.na(e_e_crit)){e_e_crit = 1}
    out_e = a * e_e_crit
    out = exp(a - out_e - cost)
  }
  
  # Least-cost
  if(stomatal_model == "phydro_least_cost"){
    cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
    out = exp((a-cost)/(par_cost$gamma*e*1e3+vcmax))
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
    out = exp(a - par_cost$alpha*jmax - par_cost$gamma*e*1e3)
  }

  ## cmax
  if(stomatal_model == "phydro_cmax"){
    aa         = par_cost$gamma
    bb         = 1
    p = psi_leaf
    out <- ((a - aa*p^2 - bb*p) - (par_cost$alpha * jmax))
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

  if (do_optim){
    return(-out)
  } else {
    return(out)
  }
}


##############################################
# INSTANTANEOUS OPTIMIZATION SCHEMES
##############################################
fn_profit_inst_schemes <- function(par, jmax, vcmax, psi_soil, e_crit, p_crit, par_cost, 
                                   par_photosynth, par_plant, par_env, 
                                   stomatal_model, do_optim){
  dpsi = par[1]#      # delta Psi in MPa
  psi_leaf = psi_soil-dpsi #MPa

  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)+1e-270  # gs in mol_co2/m2/s/Mpa
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol_h2o/m2_ground/s
  A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)$a

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
    profit = exp(A/(par_cost$gamma*e*1e3+vcmax))
  }
  
  ## CGain
  if(stomatal_model == "phydro_cgain"){
    K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl     = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    profit = exp(A - par_cost$gamma*(K-kl)/K)
  }
  
  ## WUE
  if(stomatal_model == "phydro_wue"){
    profit = exp(A - (par_cost$gamma*e*1e3))
  }

  ## cmax
  if(stomatal_model == "phydro_cmax"){
    aa         = par_cost$gamma
    bb         = 1
    p = psi_leaf
    profit <- (A - aa*p^2 - bb*p)
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


##########################################
# ACCLIMATED OPTIMIZATION PROCESS
##########################################
optimise_stomata_phydro_schemes <- function(fn_profit, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                                            par_plant, par_env, jmax_lim, return_all = FALSE, do_optim=TRUE, stomatal_model){
  if(is.null(jmax_lim)){jmax_lim = 7}
  
  jmax_ini = 0
  dpsi_ini = 1
  # count = 0
  # continue = TRUE
  # while(continue){
  out_optim <- optimr::optimr(
    par       = c(logjmax=jmax_ini, dpsi=dpsi_ini),  
    lower     = c(-10, .000001),
    upper     = c(jmax_lim, 20),
    fn             = fn_profit,
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
#   
#   dpsi_prov <- out_optim$par[2]
#   gs = calc_gs_phydro(dpsi_prov, psi_soil, par_plant, par_env)  # gs in mol_co2/m2/s/Mpa
#   e  = 1.6*gs*(par_env$vpd/par_env$patm) 
#   count <- count + 1 
#   if(all(#(psi_soil-dpsi_prov)>=p_crit, 
#          dpsi_ini != dpsi_prov)){
#     continue = FALSE
#   }else{
#     dpsi_ini = dpsi_ini*runif(1,0.1,10) 
#     if(dpsi_ini >= 20){dpsi_ini <- 19.99}
#   }
#   if(count > 20){continue = FALSE}
# }
  out_optim$value <- -out_optim$value
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}


#############################################
# INSTANTANEOUS OPTIMIZATION PROCESS
#############################################
optimise_shortterm_schemes <- function(fn_profit_inst, jmax, vcmax, psi_soil, e_crit, p_crit,
                               par_cost, par_photosynth, par_plant, par_env, 
                               stomatal_model, return_all = FALSE, do_optim){
  
  dpsi_ini = 1
  # count = 0
  # continue = TRUE
  # while(continue){
  out_optim <- optimr::optimr(
    par       = c(dpsi=dpsi_ini),  
    lower     = c(.000001),
    upper     = c(20),
    fn        = fn_profit_inst,
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
  
  # dpsi_prov <- out_optim$par[1]
  # gs = calc_gs_phydro(dpsi_prov, psi_soil, par_plant, par_env)  # gs in mol_co2/m2/s/Mpa
  # e  = 1.6*gs*(par_env$vpd/par_env$patm) 
  # count <- count + 1 
  # if(all(#(psi_soil-dpsi_prov)>=p_crit,
  #        dpsi_ini != dpsi_prov)){
  #   continue = FALSE
  # }else{
  #   dpsi_ini = dpsi_ini*runif(1,0.1,10) 
  #   if(dpsi_ini >= 20){dpsi_ini <- 19.99}
  # }
  # if(count > 20){continue = FALSE}
  # }
  
  out_optim$value <- -out_optim$value
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}


###########################################################
## MODEL OPTIMIZATION AND DATA PREPARATION WITH ACCLIMATION
###########################################################
# This function has as input the environmental variables, 
# as well as the hydraulic and cost parameters and the stomata model to be calculated. 
# The optimisation process of jmax and dpsi is then carried out. 
# The function returns jmax, dpsi, psi_leaf, gs, E, A, ci, chi and vcmax.
model_numerical <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, 
                            rdark, par_plant, par_cost, stomatal_model){
  # 1. input preparation
  patm = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm       = rpmodel::calc_kmm(tc, patm),
    gammastar = rpmodel::calc_gammastar(tc, patm),
    phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
    Iabs      = ppfd * fapar,
    ca        = co2 * patm * 1e-6,  # Convert to partial pressure
    patm      = patm,
    delta     = rdark
  )
  
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    patm            = patm,
    tc              = tc,
    vpd             = vpd
  )
  
    # Calculate psi and critical E which are based on the species parameters and 
    # for critical E also on the water potential in the soil.
    K      = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (ground) s-1 MPa-1
    p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
    e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b)+1e-270 #mol m-2 (ground) s-1

    # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit 
    # low_swp <- FALSE
    # if(psi_soil <= p_crit){
    #   psi_soil = p_crit*0.95
    #   low_swp  = TRUE
    #   print("WARNING: soil water potential is lower than critical plant water potential.")
    # }

    # 3. Optimizer
    lj_dps = optimise_stomata_phydro_schemes(fn_profit, 
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
  
    # 4. Calculate output variables
    jmax  = exp(lj_dps[1]) %>% unname()
    dpsi  = lj_dps[2] %>% unname()
    psi_l = psi_soil-dpsi
    gs    = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env) # gs in mol m-2 (ground) s-1
    a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
    a     = a_j$a
    ci    = a_j$ci
    vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
    E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
    gs    = gs #transform to umol m-2(ground) s-1 Pa-1

    # 5. Prepare output list and return
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
# This function has as input the environmental variables, 
# as well as the hydraulic and cost parameters and the stomata model to be calculated. 
# The optimisation process of jmax and dpsi is then carried out. 
# The function returns jmax, dpsi, psi_leaf, gs, E, A, ci, chi and vcmax.
model_numerical_instantaneous <- function(vcmax, jmax, tc, ppfd, vpd, co2, elv, fapar, 
                                          kphio, psi_soil,  rdark, par_plant, 
                                          par_cost, stomatal_model){
  
  # 1. input preparation
  patm = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm       = rpmodel::calc_kmm(tc, patm),
    gammastar = rpmodel::calc_gammastar(tc, patm),
    phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
    Iabs      = ppfd * fapar,
    ca        = co2 * patm * 1e-6,  # Convert to partial pressure
    patm      = patm,
    delta     = rdark
  )
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    patm            = patm,
    tc              = tc,
    vpd             = vpd
  )

  # Calculate psi and critical E which are based on the species parameters and 
  # for critical E also on the water potential in the soil.
  K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
  p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
  e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b)+1e-270 #mol m-2 (ground) s-1
    
  
  # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit
  # low_swp <- FALSE
  # if(psi_soil <= p_crit){
  #   psi_soil = p_crit*0.95
  #   low_swp  = TRUE
  #   print("WARNING: soil water potential is lower than critical plant water potential.")
  # }

  # 3. Optimizer
  lj_dps = optimise_shortterm_schemes(fn_profit_inst_schemes,
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
  
  # 4. Calculate output variables
  dpsi  = lj_dps[1] %>% unname() # delta Psi in MPa
  psi_l = psi_soil-dpsi  # leaf water potential
  gs    = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env) # gs in mol m-2 (ground) s-1
  a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
  a     = a_j$a
  ci    = a_j$ci
  E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
  gs    = gs #transform to umol m-2(ground) s-1 Pa-1
  
  # 5. Prepare output list and return
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
