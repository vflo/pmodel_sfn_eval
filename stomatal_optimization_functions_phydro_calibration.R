###########################################
## PHYDRO SCHEMES
###########################################
fn_profit <- function(par, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                      par_plant, par_env, do_optim = FALSE, stomatal_model){
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]#      # delta Psi in MPa
  psi_leaf = psi_soil-dpsi #MPa
  
  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)   # gs in mol/m2ground/s
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  
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
    out_ecrit = (a-cost)
    out_e = (a-cost)*e_e_crit
    out = out_ecrit-out_e
    if(out_e>=out_ecrit){out = 0}
  }
  
  ## wang instantaneous simplificated
  if(stomatal_model == "phydro_wang_inst_simple"){
    cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    out_ecrit = (a-cost)
    out_e = (a-cost)*e_e_crit
    out = out_ecrit-out_e
    if(out_e>=out_ecrit){out = 0}
  }
  
  ## wang modificated
  if(stomatal_model == "phydro_wang_mod"){
    cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    out_e = a * e_e_crit
    out = a - out_e - cost + 1
  }
  
  # Least-cost
  if(stomatal_model == "phydro_least_cost"){
    out = (a-par_cost$alpha*jmax)/(par_cost$gamma*e+jmax)
  }
  
  #CGain
  if(stomatal_model == "phydro_cgain"){
    K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl     = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    ks     = K*(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    out = (a- par_cost$alpha*jmax) - ((a- par_cost$alpha*jmax)*par_cost$gamma*(K-kl)/K)  #modified version
    # out = a - (par_cost$gamma*((K-kl)/K)) - par_cost$alpha*jmax
    # if (A<(par_cost$alpha * jmax)){profit <- 0}
  }
  
  # WUE
  if(stomatal_model == "phydro_wue"){
    out = ((a - par_cost$alpha*jmax)-(par_cost$gamma*e))
    if (a<(par_cost$alpha * jmax)){out <- 0}
  }

  ## cmax
  if(stomatal_model == "phydro_cmax"){
    aa         = par_cost$aa
    bb         = par_cost$bb
    p = psi_leaf
    out <- (a - aa*p**2.0 - bb*p) - (par_cost$alpha * jmax)
  }
  
  ## sperry
  if(stomatal_model == "phydro_sperry"){
    K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    # list_a = vector("list", 200)
    # for(i in 1:200){
    #   e   <- i * 0.005 * e_crit
    #   g   <- e/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    #   c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    #   list_a[[i]] <- c_a$a
    # }
    g   <- e_crit/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    amax <-  c_a$a
    ks     = K*(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    kl     = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    k_crit = 0
    risk   = (ks-kl)/(ks-k_crit)
    cost   = (par_cost$alpha * jmax)
    # amax   = max(unlist(list_a))
    out <- a/amax - risk - cost 
  }
  
  ## SOX
  if(stomatal_model == "phydro_sox"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    out <- (a-(par_cost$alpha * jmax))*reduction_factor
    if (a<(par_cost$alpha * jmax)){out <- 0}
  }
  
  ## SOX instantaneous simple
  if(stomatal_model == "phydro_sox_inst_simple"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    out <- (a-(par_cost$alpha * jmax))*reduction_factor
    if (a<(par_cost$alpha * jmax)){out <- 0}
  }
  
  ## SOX mod
  if(stomatal_model == "phydro_sox_mod"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    out <- (a)*reduction_factor - (par_cost$alpha * jmax)
    if (a<(par_cost$alpha * jmax)){out <- 0}
  }
  
  ## SOX alternative
  if(stomatal_model == "phydro_sox_alt"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(ks - k_crit)
    out <- (a)*reduction_factor-(par_cost$alpha * jmax)
  }
  
  if (do_optim){
    return(-out)
  } else {
    return(out)
  }
}

fn_profit_inst_schemes <- function(par, jmax, vcmax, psi_soil, e_crit, p_crit, par_cost, 
                                   par_photosynth, par_plant, par_env, 
                                   stomatal_model, do_optim){
  dpsi = par[1]#      # delta Psi in MPa
  psi_leaf = psi_soil-dpsi #MPa

  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)  # gs in mol/m2/s/Mpa
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)$a

  ## phydro
  if(stomatal_model == "phydro"){
    profit = A - par_cost$gamma * dpsi^2
  }
  
  ## wang
  if(stomatal_model == "phydro_wang"){
    cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    out_ecrit = (A-cost)
    out_e = (A-cost)*e_e_crit
    out = out_ecrit-out_e
    if(out_e>=out_ecrit){out = 0}
    profit <- out
  }
  
  ## wang instantaneous simplificated
  if(stomatal_model == "phydro_wang_inst_simple"){
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    profit = A*(1-e_e_crit)
  }
  
  ## wang modificated
  if(stomatal_model == "phydro_wang_mod"){
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    profit = A*(1-e_e_crit)
  }
  
  ## Least-cost #Doesn't work because of jmax in A and as cost.
  if(stomatal_model == "phydro_least_cost"){
    profit = (A-par_cost$alpha*jmax)/(par_cost$gamma*e+jmax)
  }
  
  ## CGain
  if(stomatal_model == "phydro_cgain"){
    K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl     = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    ks     = K*(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    profit = (A- par_cost$alpha*jmax) - ((A - par_cost$alpha*jmax)*par_cost$gamma*(K-kl)/K) 
    # profit = A - (A*par_cost$gamma*(K-kl)/K)  #modified version
    # profit = A - (par_cost$gamma*((K-kl)/K))
  }

  ## WUE
  if(stomatal_model == "phydro_wue"){
    profit = ((A - par_cost$alpha*jmax)- (par_cost$gamma*e))
    if (A<(par_cost$alpha * jmax)){profit <- 0}
  }

  ## cmax
  if(stomatal_model == "phydro_cmax"){
    aa         = par_cost$aa
    bb         = par_cost$bb
    p = psi_leaf
    profit <- (A - aa*p**2.0 - bb*p)
  }
  
  ## sperry
  if(stomatal_model == "phydro_sperry"){
    K          = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (ground) s-1 MPa-1
    # list_a = vector("list", 100)
    # for(i in 1:100){
    #   e   <- i * 0.01 * e_crit
    #   g   <- e/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    #   c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    #   list_a[[i]] <- c_a$a
    # }
    g   <- e_crit/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    amax <-  c_a$a
    ks     = K*(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    kl     = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    k_crit = 0
    risk   = (ks-kl)/(ks-k_crit)
    # cost   = (par_cost$alpha * jmax)
    # amax   = max(unlist(list_a))
    profit <- A/amax - risk
  }
  
  ## SOX
  if(stomatal_model == "phydro_sox"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    profit <- (A-(par_cost$alpha * jmax))*reduction_factor
    if (A<(par_cost$alpha * jmax)){profit <- 0}
  }
  
  ## SOX instantaneous simple
  if(stomatal_model == "phydro_sox_inst_simple"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    profit <- (A)*reduction_factor+100
  }
  
  ## SOX mod
  if(stomatal_model == "phydro_sox_mod"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    profit <- (A)*reduction_factor
  }
  
  ## SOX alternative
  if(stomatal_model == "phydro_sox_alt"){
    K          = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(ks - k_crit)
    profit <-(A)*reduction_factor
  }
  
  
  if (do_optim){
    return(-profit)
  } else {
    return(profit)
  }
}


optimise_stomata_phydro_schemes <- function(fn_profit, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                                            par_plant, par_env, jmax_lim, return_all = FALSE, do_optim=TRUE, stomatal_model){
  if(is.null(jmax_lim)){jmax_lim = 7}
  
  out_optim <- optimr::optimr(
    par       = c(logjmax=0, dpsi=1),  
    lower     = c(-10, .0001),
    upper     = c(10, 10),
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
  
  out_optim$value <- -out_optim$value
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}

optimise_shortterm_schemes <- function(fn_profit_inst, jmax, vcmax, psi_soil, e_crit, p_crit,
                               par_cost, par_photosynth, par_plant, par_env, 
                               stomatal_model, return_all = FALSE, do_optim){
  
  # dpsi_ini = (psi_soil - p_crit)/2
  # dpsi_ini = 0.01
  # dpsi_ini = 0.5
  dpsi_ini = 1
  out_optim <- optimr::optimr(
    par       = c(dpsi=dpsi_ini),  
    lower     = c(.000001),
    upper     = c(10),
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
model_numerical <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, 
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
model_numerical_instantaneous <- function(vcmax, jmax, tc, ppfd, vpd, co2, elv, fapar, 
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
  
  dpsi  = lj_dps[1] %>% unname() # delta Psi in MPa
  psi_l = psi_soil-dpsi  # leaf water potential
  gs    = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env) # gs in mol m-2 (ground) s-1
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
