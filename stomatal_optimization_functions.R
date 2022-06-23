#Some functions are from yugie Wang code

###########################################
## PHOTOSYNTHETIC FUNCTIONS
###########################################

get_a_ci <- function(v25,j25,gamma,gc,ca,tem,par){
  tar_p  = 0.0
  tar_a  = 0.0
  max_p  = ca
  min_p  = gamma
  adjust = 0.98
  r_day  = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
  # r_day = 0
  continue <- TRUE
  while(continue){
    tar_p = 0.5 * (max_p+min_p)
    vmax = GetPhotosyntheticVcmax(v25,tem)
    jmax = GetPhotosyntheticJmax(j25,tem)
    j = GetPhotosyntheticJ(jmax,par)
    kc = 41.01637 * 2.1**(0.1*(tem-25.0))
    ko = 28201.92 * 1.2**(0.1*(tem-25.0))
    km = kc * (1.0+21000.0/ko)
    aj = j * (tar_p-gamma) / (4.0*(tar_p+2*gamma))
    ac = vmax * (tar_p-gamma) / (tar_p+km)
    af = (aj + ac - sqrt((aj+ac)**2.0 - 4*adjust*aj*ac) ) / adjust * 0.5
    af = af - r_day
    #af = min(aj,ac) - r_day
    tmp_g = af / (ca-tar_p)
    if(abs(tmp_g-gc)/gc < 1E-12){
      tar_a = af
      continue <- FALSE
    }else if(tmp_g<gc){
      min_p = tar_p
    }else{
      max_p = tar_p
    }
    if(abs(max_p-min_p) < 1E-12){
      tar_a = af
      continue <- FALSE
    }
  }
  return(list(tar_p, tar_a))
}

# get one-point measurement Vcmax
GetOnePointVcmax <- function(ci, an, tem, gamma=2.5){
  v_min = 1.0
  v_max = 200.0
  v25   = 100.0
  while(v_max-v_min>1){
    v25   = 0.5 * (v_max + v_min)
    r_day = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
    vmax = GetPhotosyntheticVcmax(v25,tem)
    kc = 41.01637 * 2.1**(0.1*(tem-25.0))
    ko = 28201.92 * 1.2**(0.1*(tem-25.0))
    km = kc * (1.0+21000.0/ko)
    ac = vmax * (ci-gamma) / (ci+km)
    af = ac - r_day
    if ((af>an & an>0) | (an>af & an<0)){
      v_max = v25
      }else{
        v_min = v25
      }
  }
  return (v25)
}  
  
# calculate j from light
GetPhotosyntheticJ <- function(jmax, light){
  a = 0.9
  b = -0.3*light - jmax
  c = 0.3*light*jmax
  j = ( -b - sqrt(b*b-4*a*c) ) / a * 0.5
  return(j)
  }
  
# calculate jmax from temperature
GetPhotosyntheticJmax <- function(jmax25, tem){
  ha=50300.0
  hd=152044.0
  sv=495.0
  t0=298.15
  r=8.315
  c = 1.0 + exp((sv*t0 -hd)/(r*t0))
  t1 = tem + 273.15
  factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
  jmax = jmax25 * factor
  return(jmax)
  }
  
# calculate vcmax from temperature
GetPhotosyntheticVcmax <- function(vcmax25, tem){
  ha=73637.0
  hd=149252.0
  sv=486.0
  t0=298.15
  r=8.315
  c = 1.0 + exp((sv*t0 -hd)/(r*t0))
  t1 = tem + 273.15
  factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
  vcmax = vcmax25 * factor
  return(vcmax)
  }


calc_assimilation_limiting = function(vcmax, jmax, gs, par_photosynth){
  # gs = calc_gs(dpsi, psi_soil, par_plant = par_plant, par_env = par_env)
  
  # We need not employ numerical root-finding. calculate chi independently assuming Ac and Aj, and bigger of the two will be the limiting one. Accordingly return Ac or Aj
  Ac = calc_assim_rubisco_limited(gs = gs, vcmax = vcmax, par_photosynth = par_photosynth)
  Aj = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
  
  A = list(ac = Ac$a, aj=Aj$a)
  
  if (Ac$ci > Aj$ci ){
    A$ci = Ac$ci
    A$a  = Ac$a
  } else {
    A$ci = Aj$ci
    A$a  = Aj$a
  }
  
  A
}


calc_assim_rubisco_limited <- function(gs, vcmax, par_photosynth){
  
  ca = par_photosynth$ca             # ca is in Pa
  gs = gs * 1e6/par_photosynth$patm  # convert to umol/m2/s/Pa
  
  d = par_photosynth$delta
  ## with
  A <- -1.0 * gs
  B <- gs * ca - gs * par_photosynth$kmm - vcmax*(1-d)
  C <- gs * ca * par_photosynth$kmm + vcmax * (par_photosynth$gammastar + par_photosynth$kmm*d)
  
  ci <- QUADM(A, B, C)
  # a_c <- vcmax * (ci - par$gammastar) / (ci + par$kmm)
  a_c <- gs*(ca-ci) 
  
  return(list(a=a_c, ci=ci))
}


calc_assim_light_limited <- function(gs, jmax, par_photosynth){
  
  ## Only light is limiting
  ## Solve Eq. system
  ## A = gs (ca- ci)
  ## A = phi0 * Iabs * jlim * (ci - gammastar)/(ci + 2*gamma_star)
  
  ## This leads to a quadratic equation:
  ## A * ci^2 + B * ci + C  = 0
  ## 0 = a + b*x + c*x^2
  
  ca = par_photosynth$ca             # ca is in Pa
  gs = gs * 1e6/par_photosynth$patm  # convert to umol/m2/s/Pa
  
  phi0iabs = par_photosynth$phi0 * par_photosynth$Iabs
  jlim = phi0iabs / sqrt(1+ (4*phi0iabs/jmax)^2)
  
  d = par_photosynth$delta 
  ## with
  A <- -1.0 * gs
  B <- gs * ca - gs * 2 * par_photosynth$gammastar - jlim*(1-d)
  C <- gs * ca * 2*par_photosynth$gammastar + jlim * (par_photosynth$gammastar + d*par_photosynth$kmm)
  
  ci <- QUADM(A, B, C)
  aj <- gs*(ca-ci)
  # vcmax_pot <- a*(ci + par$kmm)/(ci - par$gammastar)
  
  return(list(a=aj, ci=ci))
}


calc_vcmax_coordinated_numerical <-  function(aj, ci, par_photosynth){
  d = par_photosynth$delta
  vcmax_coord = aj*(ci + par_photosynth$kmm)/(ci*(1-d) - (par_photosynth$gammastar+par_photosynth$kmm*d))
  return(vcmax_coord)
}


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


get_e_crit <- function(psi_soil, K, d, c, h, dens_water){
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
    if (abs(p-p_crit)<1e-6 | (e_max-e_min)<1e-6){
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


###########################################
## SPERRY
###########################################
get_optima_sperry <- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env){
  # 1. calculate the p_crit @ layer_f = 1E-6
  LAI        = par_env$LAI
  K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
  d          = par_plant$d
  c          = par_plant$c
  h          = par_plant$height
  dens_water = par_env$density_water
  p_crit     = d * log(1000.0) ^ (1.0/c)

  # 2. if soil psi is lower than p_crit, set soil psi to 98% p_crit 
  low_swp <- FALSE
  if(psi_soil <= p_crit){
    psi_soil = p_crit * 0.98
    low_swp  = TRUE
    print("WARNING: soil water potential is lower than critical plant water potential.")
  } #(adjust values if psi soil is less than psi critical)
  
  # 3. Calculate E critical to set the loop max
  E_crit = get_e_crit(psi_soil, K, d, c, h, dens_water) #mol m-2 (ground) s-1
  de     = 0.00001
  
  # 4. Calculate a for each e step
  list_e = vector("list", 200)
  list_k = vector("list", 200)
  list_a = vector("list", 200)
  list_p = vector("list", 200)
  list_g = vector("list", 200)
  list_ci = vector("list", 200)
  for(i in 1:200){
    e   <- i * 0.005 * E_crit
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
  return(tibble(a = opt_a, E = opt_e, p_leaf =  opt_p, gs =  opt_g, e_crit = E_crit, K = K, low_swp = low_swp ))
}


###########################################
## WANG
###########################################
get_optima_wang<- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env){
  # 1. calculate the p_crit @ layer_f = 1E-6
  LAI        = par_env$LAI
  K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
  d          = par_plant$d
  c          = par_plant$c
  h          = par_plant$height
  dens_water = par_env$density_water
  p_crit     = d * log(1000.0) ^ (1.0/c)

  # 2. if soil psi is lower than p_crit, set soil psi to 98% p_crit 
  low_swp <- FALSE
  if(psi_soil <= p_crit){
    psi_soil = p_crit * 0.98
    low_swp  = TRUE
    print("WARNING: soil water potential is lower than critical plant water potential.")
  }
  
  # 3. Optimizer
  e_crit = get_e_crit(psi_soil, K, d, c, h, dens_water) #mol m-2 (ground) s-1
  e_min  = 0.0
  e_max  = e_crit
  e_opt  = 0.0
  a_opt  = 0.0
  de     = 0.0001
  continue <- TRUE
  while(continue){
    e         = 0.5 * (e_min+e_max)
    g         = e/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a       <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    e_de      = e + de
    g_de      = e_de/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a_de    <- calc_assimilation_limiting(vcmax, jmax, g_de, par_photosynth)
    optimizer = c_a_de$a * (e_crit-e_de) - c_a$a * (e_crit-e)
    if ((e_max-e_min)<0.0000001){
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
## WOLF-ANDEREGG-PACALA
###########################################
get_optima_wap <- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env, aa=0.1, bb=0.1){
  # 1. calculate the p_crit @ layer_f = 1E-6
  LAI        = par_env$LAI
  K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
  d          = par_plant$d
  c          = par_plant$c
  h          = par_plant$height
  dens_water = par_env$density_water
  p_crit     = d * log(1000.0) ^ (1.0/c)
  
  # 2. if soil psi is lower than p_crit, set soil psi to 98% p_crit 
  low_swp <- FALSE
  if(psi_soil <= p_crit){
    psi_soil = p_crit * 0.98
    low_swp  = TRUE
    print("WARNING: soil water potential is lower than critical plant water potential.")
  }
  
  # 3. Optimizer
  e_crit = get_e_crit(psi_soil, K, d, c, h, dens_water) #mol m-2 (ground) s-1
  e_min  = 0.0
  e_max  = e_crit
  e_opt  = 0.0
  a_opt  = 0.0
  de     = 0.0001
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
    if ((e_max-e_min)<0.000001){
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
## PHYDRO
###########################################
fn_profit <- function(par, psi_soil, par_cost, e_crit, par_photosynth, par_plant, par_env, do_optim = FALSE, opt_hypothesis){
  
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]       # delta Psi in MPa
  
  gs = calc_gs(dpsi, psi_soil, par_plant, par_env) * par_env$LAI  # gs in mol/m2ground/s
  # E  = 1.6*gs*(par_env$vpd/par_env$patm)*1e6         # E in umol/m2/s
  
  ## light-limited assimilation
  a_j   <- calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
  a     = a_j$a
  # ci    = a_j$ci
  # vcmax = calc_vcmax_coordinated_numerical(a,ci, par_photosynth)
  
  costs = par_cost$alpha * jmax + par_cost$gamma * dpsi^2 #((abs((-dpsi)/par_plant$psi50)))^2  
  
  benefit = 1 #(1+1/(par_photosynth$ca/40.53))/2
  
  dummy_costs = 0*exp(20*(-abs(dpsi/4)-abs(jmax/1))) # ONLY added near (0,0) for numerical stability. 
  
  if (opt_hypothesis == "PM"){
    ## Profit Maximisation
    out <- a*benefit - costs - dummy_costs
  } else if (opt_hypothesis == "LC"){
    ## Least Cost
    out <- -(costs+dummy_costs) / (a+1e-4)
  }
  
  if (do_optim){
    return(-out)
  } else {
    return(out)
  }
}

optimise_stomata_phydro <- function(fn_profit, psi_soil, par_cost, e_crit, par_photosynth, par_plant, par_env, return_all = FALSE, opt_hypothesis){
  
  out_optim <- optimr::optimr(
    par            = c(logjmax=0, dpsi=1),  
    lower          = c(-10, .0001),
    upper          = c(10, 1e6),
    fn             = fn_profit,
    psi_soil       = psi_soil,
    e_crit         = e_crit,
    par_cost       = par_cost,
    par_photosynth = par_photosynth,
    par_plant      = par_plant,
    par_env        = par_env,
    do_optim       = TRUE,
    opt_hypothesis = opt_hypothesis,
    method         = "L-BFGS-B",
    control        = list(maxit = 500, maximize = TRUE, fnscale = 1e4)
  )
  
  out_optim$value <- -out_optim$value
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}


fn_profit_instantaneous = function(par, jmax, vcmax, psi_soil, par_cost, par_photosynth, par_plant, par_env, gs_approximation = gs_approximation, do_optim=F){
  dpsi = par

  gs = calc_gs(dpsi, psi_soil, par_plant, par_env)  # gs in mol/m2/s/Mpa
  A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)$a
  profit = A - par_cost$gamma * dpsi^2  
  
  if (do_optim){
    return(-profit)
  } else {
    return(profit)
  }
}


optimise_shortterm <- function(fn_profit_inst, jmax, vcmax, psi_soil, par_cost, par_photosynth, par_plant, par_env, gs_approximation = gs_approximation, return_all = FALSE){
  
  out_optim <- optimr::optimr(
    par       = c(dpsi=1),  
    lower     = c(0.001),
    upper     = c(20),
    fn        = fn_profit_inst,
    psi_soil  = psi_soil,
    jmax      = jmax,
    vcmax     = vcmax,
    par_cost  = par_cost,
    par_photosynth = par_photosynth,
    par_plant = par_plant,
    par_env   = par_env,
    do_optim  = TRUE,
    # opt_hypothesis = opt_hypothesis,
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
## PHYDRO + WANG
###########################################
fn_profit_phydro_wang <- function(par, psi_soil, e_crit, par_cost, par_photosynth, par_plant, par_env, do_optim = FALSE, opt_hypothesis){
  
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]       # delta Psi in MPa
  
  gs = calc_gs(dpsi, psi_soil, par_plant, par_env) * par_env$LAI  # gs in mol/m2ground/s
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  
  ## light-limited assimilation
  a_j   <- calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
  a     = a_j$a
  ci    = a_j$ci
  vcmax = calc_vcmax_coordinated_numerical(a,ci, par_photosynth)
  
  costs = (0.1*vcmax + 0.1*jmax) #umolco2 umolh2o m-4 s-2
  # out <- jmax - (a*(e_crit*1e6-e*1e6) - costs)/(e_crit*1e6) #umolco2 m-2 s-1
  out <- (a*(e_crit*1e6-e*1e6) - costs*(e_crit*1e6-e*1e6))/(e_crit*1e6) #umolco2 m-2 s-1

  if (do_optim){
    return(-out)
  } else {
    return(out)
  }
}

###########################################
## PHYDRO + WAP
###########################################
fn_profit_phydro_wap <- function(par, psi_soil, e_crit, par_cost, par_photosynth, par_plant, par_env, do_optim = FALSE, opt_hypothesis, aa=0.1, bb=0.1){
  
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]       # delta Psi in MPa
  
  LAI        = par_env$LAI
  K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*LAI #mol m-2 (ground) s-1 MPa-1
  d          = par_plant$d
  c          = par_plant$c
  h          = par_plant$height
  dens_water = par_env$density_water
  gs = calc_gs(dpsi, psi_soil, par_plant, par_env) * par_env$LAI  # gs in mol/m2ground/s
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  
  ## light-limited assimilation
  a_j   <- calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
  a     = a_j$a
  ci    = a_j$ci
  
  p = get_p(psi_soil, e, K, d, c, h, dens_water)
  out <- (a - aa*p**2.0 - bb*p) - (par_cost$alpha * jmax)
  
  if (do_optim){
    return(-out)
  } else {
    return(out)
  }
}


###########################################
## MODEL OPTIMIZATION AND DATA PREPARATION
###########################################
model_numerical <- function(tc, ppfd, vpd, nR, co2, elv, LAI, fapar, kphio, psi_soil, par_plant, stomatal_model = "sperry"){
  
  patm = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm       = rpmodel::calc_kmm(tc, patm),
    gammastar = rpmodel::calc_gammastar(tc, patm),
    phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
    Iabs      = ppfd * fapar,
    ca        = co2 * patm * 1e-6,  # Convert to partial pressure
    patm      = patm,
    delta     = 0
    # delta     = ftemp_inst_rd( tc )
  )
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    patm            = patm,
    tc              = tc,
    vpd             = vpd,
    nR              = nR,
    LAI             = LAI
  )
  
  #vcmax and jmax calculation (Arhenius)
  j25      = par_plant$j25
  v25      = par_plant$v25
  tem      = par_env$tc
  jmax_ar  = GetPhotosyntheticJmax(j25, tem)
  vcmax_ar = GetPhotosyntheticVcmax(v25, tem)
  par_plant$conductivity = par_plant$conductivity * 0.2
  # par_plant$psi50 = par_plant$psi50 / 3
  # par_plant$d = par_plant$d / 3
  # par_plant$c = 1
  # par_plant$b = 1
  K = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE)*par_env$LAI #mol m-2 (ground) s-1 MPa-1
  e_crit = get_e_crit(psi_soil, K, par_plant$d, par_plant$c, par_plant$height, par_env$density_water) #mol m-2 (ground) s-1
  
  #Sperry
  if(stomatal_model == "sperry"){
    dps <- get_optima_sperry(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env)
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
    dps <- get_optima_wang(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env)
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
  
  #Wap
  if(stomatal_model == "wap"){
    dps <- get_optima_wap(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env, aa=0.1, bb=0.1)
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
  
  #Phydro
  if(stomatal_model == "phydro"){
    par_cost = list(
      alpha = 0.1,       # cost of Jmax
      gamma = 1          # cost of hydraulic repair
    )
    p_crit = par_plant$d * log(1000.0) ^ (1.0/par_plant$c)
    # 2. if soil psi is lower than p_crit, set soil psi to 98% p_crit 
    low_swp <- FALSE
    if(psi_soil <= p_crit){
      psi_soil = p_crit*0.98
      low_swp  = TRUE
      print("WARNING: soil water potential is lower than critical plant water potential.")
    }
    # 3. Optimizer
    lj_dps = optimise_stomata_phydro(fn_profit, 
                                    psi_soil = psi_soil,
                                    e_crit = e_crit,
                                    par_cost  = par_cost, 
                                    par_photosynth = par_photosynth, 
                                    par_plant = par_plant, 
                                    par_env = par_env, 
                                    opt_hypothesis = "PM")
    
    profit = fn_profit(par = lj_dps, psi_soil = psi_soil, par_cost  = par_cost, par_photosynth = par_photosynth, par_plant = par_plant, par_env = par_env, opt_hypothesis = "PM")
    
    jmax  = exp(lj_dps[1]) %>% unname()
    dpsi  = lj_dps[2] %>% unname()
    gs    = calc_gs(dpsi, psi_soil, par_plant, par_env) * par_env$LAI # gs in mol m-2 (ground) s-1
    a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
    a     = a_j$a
    ci    = a_j$ci
    vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
    E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
    gs    = gs/patm*1e6 #transform to umol m-2(ground) s-1 Pa-1
    psi_l = psi_soil - dpsi

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
  
  ## Phydro Wang ##
  if(stomatal_model == "phydro_wang"){
    par_cost = list(
      alpha = 0.1     # cost of Jmax
    )
    p_crit = par_plant$d * log(1000.0) ^ (1.0/par_plant$c)
    # 2. if soil psi is lower than p_crit, set soil psi to 98% p_crit 
    low_swp <- FALSE
    if(psi_soil <= p_crit){
      psi_soil = p_crit*0.98
      low_swp  = TRUE
      print("WARNING: soil water potential is lower than critical plant water potential.")
    }

    # 3. Optimizer
    lj_dps = optimise_stomata_phydro(fn_profit_phydro_wang, 
                                     psi_soil = psi_soil, 
                                     e_crit = e_crit,
                                     par_cost  = par_cost, 
                                     par_photosynth = par_photosynth, 
                                     par_plant = par_plant, 
                                     par_env = par_env, 
                                     opt_hypothesis = "PM")
    
    profit = fn_profit_phydro_wang(par = lj_dps, psi_soil = psi_soil, e_crit = e_crit,
                                   par_cost  = par_cost, par_photosynth = par_photosynth, 
                                   par_plant = par_plant, par_env = par_env, opt_hypothesis = "PM")
    
    jmax  = exp(lj_dps[1]) %>% unname()
    dpsi  = lj_dps[2] %>% unname()
    gs    = calc_gs(dpsi, psi_soil, par_plant, par_env) * par_env$LAI # gs in mol m-2 (ground) s-1
    a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
    a     = a_j$a
    ci    = a_j$ci
    vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
    E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
    gs    = gs/patm*1e6 #transform to umol m-2(ground) s-1 Pa-1
    psi_l = psi_soil - dpsi
    
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
  
  ## Phydro Wap ##
  if(stomatal_model == "phydro_wap"){
    par_cost = list(
      alpha = 0.1      # cost of Jmax
    )
    p_crit = par_plant$d * log(1000.0) ^ (1.0/par_plant$c)
    # 2. if soil psi is lower than p_crit, set soil psi to 98% p_crit 
    low_swp <- FALSE
    if(psi_soil <= p_crit){
      psi_soil = p_crit*0.98
      low_swp  = TRUE
      print("WARNING: soil water potential is lower than critical plant water potential.")
    }
    
    # 3. Optimizer
    lj_dps = optimise_stomata_phydro(fn_profit_phydro_wap, 
                                     psi_soil = psi_soil, 
                                     e_crit = e_crit,
                                     par_cost  = par_cost, 
                                     par_photosynth = par_photosynth, 
                                     par_plant = par_plant, 
                                     par_env = par_env, 
                                     opt_hypothesis = "PM")
    
    profit = fn_profit_phydro_wap(par = lj_dps, psi_soil = psi_soil, e_crit = e_crit,
                                  par_cost  = par_cost, par_photosynth = par_photosynth, 
                                  par_plant = par_plant, par_env = par_env, opt_hypothesis = "PM",
                                  aa=0.1, bb=0.1)
    
    jmax  = exp(lj_dps[1]) %>% unname()
    dpsi  = lj_dps[2] %>% unname()
    gs    = calc_gs(dpsi, psi_soil, par_plant, par_env) * par_env$LAI # gs in mol m-2 (ground) s-1
    a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
    a     = a_j$a
    ci    = a_j$ci
    vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
    E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
    gs    = gs/patm*1e6 #transform to umol m-2(ground) s-1 Pa-1
    psi_l = psi_soil - dpsi
    
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
  
}
