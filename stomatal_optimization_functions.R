###########################################
## SPERRY
###########################################
get_optima_sperry <- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env){
  # 1. calculate the p_crit @ layer_f = 1E-6
  LAI        = par_env$LAI
  K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
  p_crit     = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)

  # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit 
  low_swp <- FALSE
  if(psi_soil <= p_crit){
    psi_soil = p_crit * 0.95
    low_swp  = TRUE
    print("WARNING: soil water potential is lower than critical plant water potential.")
  } #(adjust values if psi soil is less than psi critical)
  
  # 3. Calculate E critical to set the loop max
  e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1
  de     = 0.00000001
  
  # 4. Calculate a for each e step
  list_e = vector("list", 200)
  list_k = vector("list", 200)
  list_a = vector("list", 200)
  list_p = vector("list", 200)
  list_g = vector("list", 200)
  list_ci = vector("list", 200)
  for(i in 1:200){
    e   <- i * 0.005 * e_crit
    p   <- get_p(psi_soil, e, K, d, c, h, dens_water)
    g   <- e/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    e_de = e + de
    p_de = get_p(psi_soil, e_de, K, d, c, h, dens_water)
    k = de / (abs(p_de)-abs(p))
    list_e[[i]] <- e
    list_k[[i]] <- k
    list_a[[i]] <- c_a$a
    list_ci[[i]] <- c_a$ci
    list_p[[i]] <- p
    list_g[[i]] <- g
  }
  
  # 5. Calculate Maxim Profit
  gain = unlist(list_a)/max(unlist(list_a))
  risk = 1.0 - unlist(list_k)/max(unlist(list_k))
  prof = gain - risk
  opt_site = which(prof == max(prof))
  opt_a = list_a[[opt_site]]
  opt_e = list_e[[opt_site]]
  opt_g = list_g[[opt_site]]
  opt_p = list_p[[opt_site]]
  
  # 6. return optimized values
  return(tibble(a = opt_a, E = opt_e, p_leaf =  opt_p, gs =  opt_g, e_crit = e_crit, K = K, low_swp = low_swp ))
}


###########################################
## WANG
###########################################
get_optima_wang<- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env){
  # 1. calculate the p_crit @ layer_f = 1E-6
  LAI        = par_env$LAI
  K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
  p_crit     = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)

  # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit 
  low_swp <- FALSE
  if(psi_soil <= p_crit){
    psi_soil = p_crit * 0.95
    low_swp  = TRUE
    print("WARNING: soil water potential is lower than critical plant water potential.")
  }
  
  # 3. Optimizer
  e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1
  e_min  = 0.0
  e_max  = e_crit
  e_opt  = 0.0
  a_opt  = 0.0
  de     = 0.00000001
  continue <- TRUE
  while(continue){
    e         = 0.5 * (e_min+e_max)
    g         = e/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a       <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    e_de      = e + de
    g_de      = e_de/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a_de    <- calc_assimilation_limiting(vcmax, jmax, g_de, par_photosynth)
    optimizer = c_a_de$a * (e_crit-e_de) - c_a$a * (e_crit-e)
    if ((e_max-e_min)<0.00000000001){
      e_opt = e
      a_opt = c_a$a
      p_opt = get_p(psi_soil, e, K, d, c, h, dens_water)
      g_opt =  g
      continue <- FALSE
      }
    if (optimizer>0){
      e_min = e
      }else{
        e_max = e
      }
  }
  # 4. return optimized values
  return(tibble(a = a_opt, E = e_opt, p_leaf =  p_opt, gs = g_opt, e_crit = e_crit, low_swp = low_swp))
}


###########################################
## WOLF-ANDEREGG-PACALA (CMAX)
###########################################
get_optima_cmax <- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env, aa=0.1, bb=0.1){
  # 1. calculate the p_crit @ layer_f = 1E-6
  LAI        = par_env$LAI
  K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
  p_crit     = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
  
  # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit 
  low_swp <- FALSE
  if(psi_soil <= p_crit){
    psi_soil = p_crit * 0.95
    low_swp  = TRUE
    print("WARNING: soil water potential is lower than critical plant water potential.")
  }
  
  # 3. Optimizer
  e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1
  e_min  = 0.0
  e_max  = e_crit
  e_opt  = 0.0
  a_opt  = 0.0
  de     = 0.00000001
  continue = TRUE
  while(continue){
    e          = 0.5 * (e_min+e_max)
    p          = get_p(psi_soil, e, K, d, c, h, dens_water)
    g          = e/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a        <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    e_de       = e + de
    p_de       = get_p(psi_soil, e_de, K, d, c, h, dens_water)
    g_de       = e_de/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a_de     <- calc_assimilation_limiting(vcmax, jmax, g_de, par_photosynth)
    optimizer  = c_a_de$a - aa*p_de**2.0 - bb*p_de - c_a$a + aa*p**2.0 + bb*p
    if ((e_max-e_min)<0.00000000001){
      e_opt    = e
      a_opt    = c_a$a
      p_opt    = get_p(psi_soil, e, K, d, c, h, dens_water)
      g_opt    =  g
      continue = FALSE
    }
    if (optimizer>0){
      e_min = e
    }else{
      e_max = e
    }
  }
  
  # 4. return optimized values
  return(tibble(a = a_opt, E = e_opt, p_leaf =  p_opt, gs = g_opt, e_crit = e_crit, low_swp = low_swp))
}



###########################################
## SOX
###########################################
get_optima_sox <- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env){
  # 1. calculate the p_crit @ layer_f = 1E-6
  LAI        = par_env$LAI
  K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
  p_crit     = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
  
  # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit 
  low_swp <- FALSE
  if(psi_soil <= p_crit){
    psi_soil = p_crit * 0.95
    low_swp  = TRUE
    print("WARNING: soil water potential is lower than critical plant water potential.")
  }
  
  # 3. Optimizer
  e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1
  e_min  = 0.0
  e_max  = e_crit
  e_opt  = 0.0
  a_opt  = 0.0
  de     = 0.0000001
  continue = TRUE
  while(continue){
    e          = 0.5 * (e_min+e_max)
    p          = get_p(psi_soil, e, K, d, c, h, dens_water)
    g          = e/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a        <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    e_de       = e + de
    p_de       = get_p(psi_soil, e_de, K, d, c, h, dens_water)
    g_de       = e_de/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a_de     <- calc_assimilation_limiting(vcmax, jmax, g_de, par_photosynth)
    e_df       = e_de + de
    p_df       = get_p(psi_soil, e_df, K, d, c, h, dens_water)
    m          = de / (p_de - p)
    m_de       = de / (p_df - p_de)
    optimizer  = c_a_de$a * m_de - c_a$a * m
    if ((e_max-e_min)<0.00000000001){
      e_opt    = e
      a_opt    = c_a$a
      p_opt    = get_p(psi_soil, e, K, d, c, h, dens_water)
      g_opt    =  g
      continue = FALSE
    }
    if (optimizer>0){
      e_min = e
    }else{
      e_max = e
    }
  }
  
  # 4. return optimized values
  return(tibble(a = a_opt, E = e_opt, p_leaf =  p_opt, gs = g_opt, e_crit = e_crit, low_swp = low_swp))
}


###########################################
## PHYDRO SCHEMES
###########################################
fn_profit <- function(par, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                      par_plant, par_env, do_optim = FALSE, stomatal_model){
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]#      # delta Psi in MPa
  psi_leaf = psi_soil-dpsi #MPa
  
  gs = calc_gs(dpsi, psi_soil, par_plant, par_env) * par_env$LAI  # gs in mol/m2ground/s
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
  }
  
  # ## cmax
  # if(stomatal_model == "phydro_cmax"){
  #   aa         = par_cost$aa
  #   bb         = par_cost$bb
  #   p = psi_leaf
  #   out <- (a - aa*p**2.0 - bb*p) - (par_cost$alpha * jmax)
  # }
  
  ## sperry
  if(stomatal_model == "phydro_sperry"){
    LAI        = par_env$LAI
    K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
    list_a = vector("list", 200)
    for(i in 1:200){
      e   <- i * 0.005 * e_crit
      g   <- e/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
      c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
      list_a[[i]] <- c_a$a
    }
    ks     = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl     = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit = 0
    risk   = (ks-kl)/(ks-k_crit)
    cost   = (par_cost$alpha * jmax)
    amax   = max(unlist(list_a))
    out <- a - amax*risk - cost
  }
  
  ## SOX
  if(stomatal_model == "phydro_sox"){
    LAI        = par_env$LAI
    K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(ks - k_crit)
    out <- a*reduction_factor - (par_cost$alpha * jmax)
  }
  
  if (do_optim){
    return(-out)
  } else {
    return(out)
  }
}

fn_profit_instantaneous <- function(par, jmax, vcmax, psi_soil, e_crit, p_crit, par_cost, 
                                   par_photosynth, par_plant, par_env, 
                                   stomatal_model, do_optim=F){
  dpsi = par[1]#      # delta Psi in MPa
  psi_leaf = psi_soil-dpsi #MPa

  gs = calc_gs(dpsi, psi_soil, par_plant, par_env)  # gs in mol/m2/s/Mpa
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)$a

  ## phydro
  if(stomatal_model == "phydro"){
    profit = A - par_cost$gamma * dpsi^2
  }
  
  ## wang
  if(stomatal_model == "phydro_wang"){
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    profit = A*e_e_crit
  }
  
  # ## cmax
  # if(stomatal_model == "phydro_cmax"){
  #   aa         = par_cost$aa
  #   bb         = par_cost$bb
  #   p = psi_leaf
  #   profit <- (A - aa*p**2.0 - bb*p)
  # }
  
  ## sperry
  if(stomatal_model == "phydro_sperry"){
    LAI        = par_env$LAI
    K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
    list_a = vector("list", 200)
    for(i in 1:200){
      e   <- i * 0.005 * e_crit
      g   <- e/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
      c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
      list_a[[i]] <- c_a$a
    }
    ks     = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl     = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit = 0
    risk   = (ks-kl)/(ks-k_crit)
    amax   = max(unlist(list_a))
    profit = A - amax*risk
  }
  
  ## SOX
  if(stomatal_model == "phydro_sox"){
    LAI        = par_env$LAI
    K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
    ks         = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl         = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit     = 0
    reduction_factor = (kl - k_crit)/(ks - k_crit)
    profit <- A*reduction_factor
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
    par            = c(logjmax=-4, dpsi=(psi_soil-p_crit)/2),  
    lower          = c(-10, 0.001),
    upper          = c(jmax_lim,(psi_soil-p_crit)),
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
                               stomatal_model, return_all = FALSE){
  
  out_optim <- optimr::optimr(
    par       = c(dpsi=(psi_soil-p_crit)/2),  
    lower     = c(0.001),
    upper     = c((psi_soil-p_crit)),
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
    do_optim  = TRUE, 
    stomatal_model = stomatal_model,
    method    = "L-BFGS-B",
    control   = list() 
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
model_numerical <- function(tc, ppfd, vpd, co2, elv, LAI, fapar, kphio, psi_soil, r_dark, par_plant, par_cost = NULL, stomatal_model = "sperry"){

  patm = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm       = rpmodel::calc_kmm(tc, patm),
    gammastar = rpmodel::calc_gammastar(tc, patm),
    phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
    Iabs      = ppfd * fapar,
    ca        = co2 * patm * 1e-6,  # Convert to partial pressure
    patm      = patm,
    delta     = r_dark
    # delta     = ftemp_inst_rd( tc )
  )
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    patm            = patm,
    tc              = tc,
    vpd             = vpd,
    LAI             = LAI
  )
  if(is.null(par_cost)){
    par_cost = list(
      alpha = 0.1,       # cost of Jmax
      gamma = 1,          # cost of hydraulic repair
      aa    = 0.1,       # CMax parabola coefficient
      bb    = 0.1        # CMax parabola coefficient
    )
  }
  
  #vcmax and jmax calculation (Arhenius)
  j25      = par_plant$j25
  v25      = par_plant$v25
  tem      = par_env$tc
  jmax_ar  = GetPhotosyntheticJmax(j25, tem)
  vcmax_ar = GetPhotosyntheticVcmax(v25, tem)

  #Sperry
  if(stomatal_model == "sperry"){
    dps <- get_optima_sperry(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, 
                             psi_soil, par_photosynth, par_plant, par_env)
    return(list(
      p_leaf  = dps$p_leaf,
      jmax    = jmax_ar,
      vcmax   = vcmax_ar,
      E       = dps$E,
      a       = dps$a,
      gs      = dps$gs,
      e_crit  = dps$e_crit,
      K       = dps$K,
      low_swp = dps$low_swp
    ))
  }
  
  #Wang
  if(stomatal_model == "wang"){
    dps <- get_optima_wang(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, 
                           psi_soil, par_photosynth, par_plant, par_env)
    return(list(
      p_leaf  = dps$p_leaf,
      jmax    = jmax_ar,
      vcmax   = vcmax_ar,
      E       = dps$E,
      a       = dps$a,
      gs      = dps$gs,
      e_crit  = dps$e_crit,
      low_swp = dps$low_swp
    ))
  }
  
  # #cmax
  # if(stomatal_model == "cmax"){
  #   dps <- get_optima_cmax(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, 
  #                         psi_soil, par_photosynth, par_plant, par_env, 
  #                         aa = par_cost$aa, bb = par_cost$bb)
  #   return(list(
  #     p_leaf  = dps$p_leaf,
  #     jmax    = jmax_ar,
  #     vcmax   = vcmax_ar,
  #     E       = dps$E,
  #     a       = dps$a,
  #     gs      = dps$gs,
  #     e_crit  = dps$e_crit,
  #     low_swp = dps$low_swp
  #   ))
  # }
  
  #SOX
  if(stomatal_model == "sox"){
    dps <- get_optima_sox(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, 
                          psi_soil, par_photosynth, par_plant, par_env)
    return(list(
      p_leaf  = dps$p_leaf,
      jmax    = jmax_ar,
      vcmax   = vcmax_ar,
      E       = dps$E,
      a       = dps$a,
      gs      = dps$gs,
      e_crit  = dps$e_crit,
      low_swp = dps$low_swp
    ))
  }
  
  #Phydro/phydro_cmax/phydro_sperry/phydro_sox
  if(stomatal_model %in% c("phydro", "phydro_cmax", 'phydro_sperry',"phydro_sox","phydro_wang")){
    K      = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
    p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
    e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1

    # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit 
    low_swp <- FALSE
    if(psi_soil <= p_crit){
      psi_soil = p_crit*0.95
      low_swp  = TRUE
      print("WARNING: soil water potential is lower than critical plant water potential.")
    }
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
    
    profit = fn_profit(par = lj_dps, psi_soil = psi_soil, par_cost  = par_cost,
                       e_crit = e_crit, p_crit = p_crit, par_photosynth = par_photosynth, 
                       par_plant = par_plant, par_env = par_env, do_optim = FALSE, 
                       stomatal_model = stomatal_model)
    
    jmax  = exp(lj_dps[1]) %>% unname()
    dpsi  = lj_dps[2] %>% unname()
    psi_l = psi_soil-dpsi
    gs    = calc_gs(dpsi, psi_soil, par_plant, par_env) * par_env$LAI # gs in mol m-2 (ground) s-1
    a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
    a     = a_j$a
    ci    = a_j$ci
    vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
    E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
    gs    = gs/patm*1e6 #transform to umol m-2(ground) s-1 Pa-1

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
      chi_jmax_lim = 0,
      profit       = profit,
      low_swp      = low_swp
    ))
  }
}


########################################################
## INSTANTANEOUS MODEL OPTIMIZATION AND DATA PREPARATION
########################################################
model_numerical_instantaneous <- function(tc, ppfd, vpd, co2, elv, LAI, fapar, 
                                          kphio, r_dark, psi_soil, jmax, vcmax, par_plant, 
                                          par_cost = NULL, stomatal_model = "phydro_sperry"){
  
  patm = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm       = rpmodel::calc_kmm(tc, patm),
    gammastar = rpmodel::calc_gammastar(tc, patm),
    phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
    Iabs      = ppfd * fapar,
    ca        = co2 * patm * 1e-6,  # Convert to partial pressure
    patm      = patm,
    delta     = r_dark,
    # delta     = ftemp_inst_rd( tc )
  )
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    patm            = patm,
    tc              = tc,
    vpd             = vpd,
    LAI             = LAI
  )
  if(is.null(par_cost)){
    par_cost = list(
      alpha = 0.1,       # cost of Jmax
      gamma = 1,          # cost of hydraulic repair
      aa    = 0.1,       # CMax parabola coefficient
      bb    = 0.1        # CMax parabola coefficient
    )
  }
  
  #vcmax and jmax calculation (Arhenius)
  j25      = par_plant$j25
  v25      = par_plant$v25
  tem      = par_env$tc
  jmax_ar  = GetPhotosyntheticJmax(j25, tem)
  vcmax_ar = GetPhotosyntheticVcmax(v25, tem)
  K = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*par_env$LAI #mol m-2 (ground) s-1 MPa-1
  p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
  e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1
    
  # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit 
  low_swp <- FALSE
  if(psi_soil <= p_crit){
     psi_soil = p_crit*0.95
     low_swp  = TRUE
     print("WARNING: soil water potential is lower than critical plant water potential.")
     }
    
  # 3. Optimizer
  lj_dps = optimise_shortterm_schemes(fn_profit_instantaneous,
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
    
  profit = fn_profit_instantaneous(par = lj_dps, jmax = jmax, vcmax = vcmax,
                                   psi_soil = psi_soil, e_crit = e_crit, 
                                   p_crit= p_crit, par_cost  = par_cost, 
                                   par_photosynth = par_photosynth, 
                                   par_plant = par_plant, par_env = par_env, 
                                   do_optim = TRUE, stomatal_model = stomatal_model)
  
  dpsi  = lj_dps[1] %>% unname() # delta Psi in MPa
  psi_l = psi_soil-dpsi  # leaf water potential
  gs    = calc_gs(dpsi, psi_soil, par_plant, par_env) * par_env$LAI # gs in mol m-2 (ground) s-1
  a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
  a     = a_j$a
  ci    = a_j$ci
  E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
  gs    = gs/patm*1e6 #transform to umol m-2(ground) s-1 Pa-1
    
  return(list(
      jmax         = jmax,
      profit       = profit,
      dpsi         = dpsi,
      p_leaf       = psi_l,
      gs           = gs,
      E            = E,
      a            = a,
      ci           = ci,
      chi          = ci/par_photosynth$ca,
      vcmax        = vcmax,
      chi_jmax_lim = 0,
      low_swp      = low_swp
  ))

}
