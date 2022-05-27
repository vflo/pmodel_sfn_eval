#Some functions are from Jugie Wang code

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


calc_assimilation_limiting_max = function(vcmax, jmax, v25, gs_min, gs_crit, par_photosynth, par_env){
  A <- tibble::tibble(ci=NA,a=NA)
  gs <- seq(gs_crit, gs_min, length.out = 100)
  for (i in 1:100){
    
    # We need not employ numerical root-finding. calculate chi independently assuming Ac and Aj, and bigger of the two will be the limiting one. Accordingly return Ac or Aj
    Ac = calc_assim_rubisco_limited(gs = gs[i], vcmax = vcmax, par_photosynth = par_photosynth, v25 = v25, par_env = par_env)
    Aj = calc_assim_light_limited(gs = gs[i], jmax = jmax, par_photosynth = par_photosynth, v25 = v25, par_env = par_env)
    A[i,'ac'] = Ac$a
    A[i,'aj'] = Aj$a
    
    if (Ac$ci > Aj$ci ){
      A[i,'ci'] = Ac$ci
      A[i,'a']  = Ac$a
    } else {
      A[i,'ci'] = Aj$ci
      A[i,'a']  = Aj$a
    }
  }
  A[which(A$a == max(A$a)),]
}


calc_assimilation_limiting <- function (vcmax, jmax, gs, par_photosynth, v25, par_env) {
  Ac = calc_assim_rubisco_limited(gs = gs, vcmax = vcmax, 
                                  par_photosynth = par_photosynth,
                                  v25 = v25, par_env = par_env)
  Aj = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth, v25 = v25, par_env = par_env)
  A = list(ac = Ac$a, aj = Aj$a)
  if (Ac$ci > Aj$ci) {
    A$ci = Ac$ci
    A$a = Ac$a
  }
  else {
    A$ci = Aj$ci
    A$a = Aj$a
  }
  A
}


calc_assim_rubisco_limited <- function (gs, vcmax, par_photosynth, v25, par_env){
  tem <- par_env$tc
  ca = par_photosynth$ca
  gs = gs * 1e+06/par_photosynth$patm
  # d = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0))) # day respiration as f(temperature)
  d = par_photosynth$delta
  A <- -1 * gs
  B <- gs * ca - gs * par_photosynth$kmm - vcmax * (1 - d)
  C <- gs * ca * par_photosynth$kmm + vcmax * (par_photosynth$gammastar + 
                                                 par_photosynth$kmm * d)
  ci <- QUADM(A, B, C)
  a_c <- gs * (ca - ci)
  return(list(a = a_c, ci = ci))
}


calc_assim_light_limited <- function (gs, jmax, par_photosynth, v25, par_env){
  tem <- par_env$tc
  ca = par_photosynth$ca
  gs = gs * 1e+06/par_photosynth$patm
  phi0iabs = par_photosynth$phi0 * par_photosynth$Iabs
  jlim = phi0iabs/sqrt(1 + (4 * phi0iabs/jmax)^2)
  # d = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0))) # day respiration as f(temperature)
  d = par_photosynth$delta
  A <- -1 * gs
  B <- gs * ca - gs * 2 * par_photosynth$gammastar - jlim * (1 - d)
  C <- gs * ca * 2 * par_photosynth$gammastar + jlim * (par_photosynth$gammastar + d * par_photosynth$kmm)
  ci <- QUADM(A, B, C)
  aj <- gs * (ca - ci)
  return(list(a = aj, ci = ci))
}

calc_vcmax_coordinated_numerical <-  function(aj, ci, par_photosynth, v25, par_env){
  tem <- par_env$tc
  # d = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0))) # day respiration as f(temperature)
  d = par_photosynth$delta
  vcmax_coord = aj*(ci + par_photosynth$kmm)/(ci*(1-d) - (par_photosynth$gammastar+par_photosynth$kmm*d))
  return(vcmax_coord)
}


###########################################
## HYDRAULIC FUNCTIONS
###########################################
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

## Returns conductivity in mol/m2/s/Mpa
scale_conductivity_Ks = function(K, par_env, par_plant){
  # Flow rate in Kg/m/s/MPa
  K2 = K / par_plant$height
  # Flow rate in mol/m2/s/MPa
  mol_h20_per_kg_h20 = 55.5
  K3 = K2  * mol_h20_per_kg_h20

  return(K3)  
}

calc_gs <- function (dpsi, psi_soil, par_plant, par_env, ...){
  # K = scale_conductivity(par_plant$conductivity, par_env)
  K = par_plant$conductivity
  D = (par_env$vpd/par_env$patm)
  K/1.6/D * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b, 
                        ...)
}

integral_P <- function(dpsi, psi_soil, psi50, b, ...){
  ps = psi_soil/psi50
  pl = (psi_soil-dpsi)/psi50
  l2 = log(2)
  -(psi50/b)*(l2^(-1/b))*(expint::gammainc(a = 1/b, x = l2*pl^b)-expint::gammainc(a = 1/b, x = l2*ps^b))
}


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


###########################################
## SPERRY
###########################################
get_optima_sperry <- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env){
  # 1. calculate the p_crit @ layer_f = 1E-6
  # K = scale_conductivity_Ks(par_plant$conductivity, par_env,par_plant)
  K = par_plant$conductivity
  d = par_plant$d
  c = par_plant$c
  h = par_plant$height
  p_crit = d * log(1000.0) ^ (1.0/c)
  low_swp <- FALSE
  if(p_crit>=psi_soil){psi_soil <- p_crit*0.99
  low_swp <- TRUE
  print("WARNING: soil water potential is lower than critical plant water potential.")} #(adjust values if psi soil is less than psi critical)
  E_crit = get_e_crit(psi_soil, K, d, c, h)
  de = 0.1
  # 2. increase the e stepwise
  list_e = vector("list", 100)
  list_k = vector("list", 100)
  list_a = vector("list", 100)
  list_p = vector("list", 100)
  list_ci = vector("list", 100)
  for(i in 1:100){
    e   <- i * 0.01 * E_crit
    p   <- get_p(psi_soil, e, K, d, c, h)
    g   <- e/1.6/(par_env$vpd/par_env$patm)/1e3
    c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth,v25, par_env)
    # c_a <- get_a_ci(v25,j25,par_photosynth$gammastar,g,par_photosynth$ca,par_env$tc,par_photosynth$Iabs)
    # k <- K*exp(-(p/d)^c)
    e_de = e + de
    p_de = get_p(psi_soil, e_de, K, d, c, h)
    k = de / (abs(p_de)-abs(p))
    list_e[[i]] <- e
    list_k[[i]] <- k
    list_a[[i]] <- c_a$a
    list_ci[[i]] <- c_a$ci
    # list_a[[i]] <- c_a[[2]]
    # list_ci[[i]] <- c_a[[1]]
    list_p[[i]] <- p
  }
# 3. extend the lists
  gain = unlist(list_a)/max(unlist(list_a))
  risk = 1.0 - unlist(list_k)/max(unlist(list_k))
  prof = gain - risk
  opt_site = which(prof == max(prof))
  opt_a = list_a[[opt_site]]
  opt_e = list_e[[opt_site]]
  opt_p = list_p[[opt_site]]
  return(tibble(a = opt_a, E = opt_e, p_leaf =  opt_p, low_swp = low_swp ))
}


###########################################
## WANG
###########################################
get_optima_wang<- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env){
  K = par_plant$conductivity
  d = par_plant$d
  c = par_plant$c
  h = par_plant$height
  p_crit = d * log(1000.0) ^ (1.0/c)
  low_swp <- FALSE
  if(p_crit>=psi_soil){psi_soil <- p_crit*0.99
  low_swp <- TRUE
  print("WARNING: soil water potential is lower than critical plant water potential.")} #(adjust values if psi soil is less than psi critical)
  e_crit = get_e_crit(psi_soil, K, d, c, h)
  e_min  = 0.0
  e_max  = e_crit
  e_opt  = 0.0
  a_opt  = 0.0
  de     = 1.0
  continue <- TRUE
  while(continue){
    e   = 0.5 * (e_min+e_max)
    g   <- e/1.6/(par_env$vpd/par_env$patm)/1e3
    c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth,v25, par_env)
    # c_a <- get_a_ci(v25,j25,par_photosynth$gammastar,g,par_photosynth$ca,par_env$tc,par_photosynth$Iabs)
    e_de = e + de
    g_de   <- e_de/1.6/(par_env$vpd/par_env$patm)/1e3
    c_a_de <- calc_assimilation_limiting(vcmax, jmax, g_de, par_photosynth,v25, par_env)
    # c_a_de <- get_a_ci(v25,j25,par_photosynth$gammastar,g_de,par_photosynth$ca,par_env$tc,par_photosynth$Iabs)
    optimizer = c_a_de$a * (e_crit-e_de) - c_a$a * (e_crit-e)
    # optimizer = c_a_de[[2]] * (e_crit-e_de) - c_a[[2]] * (e_crit-e)
    if ((e_max-e_min)<0.001){
      e_opt = e
      a_opt = c_a$a
      # a_opt = c_a[[2]]
      p_opt = get_p(psi_soil, e, K, d, c, h)
      continue <- FALSE
      }
    if (optimizer>0){
      e_min = e
      }else{
        e_max = e
      }
  }
  return(tibble(a = a_opt, E = e_opt, p_leaf =  p_opt, low_swp = low_swp))
}


###########################################
## WOLF-ANDEREGG-PACALA
###########################################
get_optima_wap <- function(jmax = jmax, vcmax = vcmax, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env, aa=0.1, bb=0.1){
  K = par_plant$conductivity
  d = par_plant$d
  c = par_plant$c
  h = par_plant$height
  p_crit = d * log(1000.0) ^ (1.0/c)
  low_swp <- FALSE
  if(p_crit>=psi_soil){psi_soil <- p_crit*0.99
  low_swp <- TRUE
  print("WARNING: soil water potential is lower than critical plant water potential.")} #(adjust values if psi soil is less than psi critical)
  e_crit = get_e_crit(psi_soil, K, d, c, h)
  e_min  = 0.0
  e_max  = e_crit
  e_opt  = 0.0
  a_opt  = 0.0
  de     = 1.0
  continue <- TRUE
  while(continue){
    e   = 0.5 * (e_min+e_max)
    p = get_p(psi_soil, e, K, d, c, h)
    g   <- e/1.6/(par_env$vpd/par_env$patm)/1e3
    c_a <- calc_assimilation_limiting(vcmax, jmax, g, par_photosynth,v25, par_env)
    # c_a <- get_a_ci(v25,j25,par_photosynth$gammastar,g,par_photosynth$ca,par_env$tc,par_photosynth$Iabs)
    e_de = e + de
    p_de = get_p(psi_soil, e_de, K, d, c, h)
    g_de   <- e_de/1.6/(par_env$vpd/par_env$patm)/1e3
    c_a_de <- calc_assimilation_limiting(vcmax, jmax, g_de, par_photosynth,v25, par_env)
    # c_a_de <- get_a_ci(v25,j25,par_photosynth$gammastar,g_de,par_photosynth$ca,par_env$tc,par_photosynth$Iabs)
    optimizer = c_a_de$a - aa*p_de**2.0 - bb*p_de - c_a$a + aa*p**2.0 + bb*p
    # optimizer = c_a_de[[2]] - aa*p_de**2.0 - bb*p_de - c_a[[2]] + aa*p**2.0 + bb*p
    if ((e_max-e_min)<0.001){
      e_opt = e
      a_opt = c_a$a
      # a_opt = c_a[[2]]
      p_opt = get_p(psi_soil, e, K, d, c, h)
      continue <- FALSE
    }
    if (optimizer>0){
      e_min = e
    }else{
      e_max = e
    }
  }
  return(tibble(a = a_opt, E = e_opt, p_leaf =  p_opt, low_swp = low_swp))
}

###########################################
## PHYDRO
###########################################
optimise_stomata_phydro <- function(fn_profit, v25, psi_soil, par_photosynth, par_plant, par_env, par_cost = par_cost, return_all = FALSE){
  
  out_optim <- optimr::optimr(
    par       = c(logjmax=0, dpsi=1),  
    lower     = c(-10, .0001),
    upper     = c(10, 20),
    fn        = fn_profit,
    psi_soil  = psi_soil,
    v25 = v25,
    par_photosynth = par_photosynth,
    par_cost = par_cost,
    par_plant = par_plant,
    par_env   = par_env,
    method    = "L-BFGS-B",
    control   = list( maxit = 5000, maximize = TRUE, fnscale=1e1 )
  )
  
  out_optim$value <- -out_optim$value

  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}


phydro_optim = function(par, jmax=jmax, vcmax=vcmax, v25 = v25, psi_soil, par_cost, par_photosynth, par_plant, par_env){
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]       # delta Psi in MPa
  gs = calc_gs(dpsi, psi_soil, par_plant, par_env)/1e3  # gs in mol/m2/s/Mpa
  
  ## light-limited assimilation
  a_j <- calc_assim_light_limited(gs, jmax, par_photosynth,v25,par_env) # Aj in umol/m2/s
  a = a_j$a
  ci = a_j$ci
  vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth, v25, par_env)
  costs = par_cost$alpha * jmax + par_cost$gamma * dpsi^2 #((abs((-dpsi)/par_plant$psi50)))^2  
  a = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth, v25, par_env)$a
  out <- costs - a
  
  return(out)
  
}



###########################################
## MODEL OPTIMIZATION AND DATA PREPARATION
###########################################
model_numerical <- function(tc, ppfd, vpd, nR, co2, elv, LAI, fapar, kphio, psi_soil, par_plant, stomatal_model = "sperry"){
  
  p = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm = rpmodel::calc_kmm(tc, p),
    gammastar = rpmodel::calc_gammastar(tc, p),
    phi0 = kphio*rpmodel::calc_ftemp_kphio(tc),
    Iabs = ppfd*fapar,
    ca = co2*p*1e-6,  # Convert to partial pressure
    patm = p,
    delta = 0
  )
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, p),  # Needs to be imported from rpmodel.R
    density_water = rpmodel::calc_density_h2o(tc, p),  # Needs to be imported from rpmodel.R
    patm = p,
    tc = tc,
    vpd = vpd,
    nR = nR,
    LAI = LAI
  )
  #vcmax and jmax calculation (Arhenius)
  j25 = par_plant$j25
  v25 = par_plant$v25
  tem = par_env$tc
  jmax_ar = GetPhotosyntheticJmax(j25,tem)
  vcmax_ar = GetPhotosyntheticVcmax(v25,tem)
  
  if(stomatal_model == "sperry"){
    dps = get_optima_sperry(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env)
    return(list(
      p_leaf=dps$p_leaf,
      jmax=jmax_ar,
      vcmax=vcmax_ar,
      E=dps$E,
      a=dps$a
    ))
  }
  
  if(stomatal_model == "wang"){
    dps = get_optima_wang(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env)
    return(list(
      p_leaf=dps$p_leaf,
      jmax=jmax_ar,
      vcmax=vcmax_ar,
      E=dps$E,
      a=dps$a
    ))
  }
  
  if(stomatal_model == "wap"){
    dps = get_optima_wap(jmax = jmax_ar, vcmax = vcmax_ar, j25 = j25, v25 = v25, psi_soil, par_photosynth, par_plant, par_env, aa=0.1, bb=0.1)
    return(list(
      p_leaf=dps$p_leaf,
      jmax=jmax_ar,
      vcmax=vcmax_ar,
      E=dps$E,
      a=dps$a
    ))
  }
  
  if(stomatal_model == "phydro"){
    par_cost = list(
      alpha = 0.1,       # cost of Jmax
      gamma = 1          # cost of hydraulic repair
    )
    p_crit = par_plant$d * log(1000.0) ^ (1.0/par_plant$c)
    low_swp <- FALSE
    if(p_crit>=psi_soil){psi_soil <- p_crit*0.99
    low_swp <- TRUE
    print("WARNING: soil water potential is lower than critical plant water potential.")} #(adjust values if psi soil is less than psi critical)
    lj_dps = optimise_stomata_phydro(phydro_optim,
                                     v25 = v25,
                                     psi_soil = psi_soil,
                                     par_photosynth = par_photosynth,
                                     par_plant = par_plant,
                                     par_env = par_env,
                                     par_cost = par_cost)
    jmax = exp(lj_dps[1])
    dpsi = lj_dps[2]
    gs = calc_gs(dpsi, psi_soil, par_plant, par_env)/1e3  # gs in mol/m2/s/Mpa
    E = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2/s
    a_j = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth, v25=v25, par_env = par_env)
    a = a_j$a
    ci = a_j$ci
    vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth, v25=v25, par_env = par_env)
    return(list(
      jmax=jmax,
      dpsi=dpsi,
      gs=gs,
      E=E,
      a=a,
      ci=ci,
      chi = ci/par_photosynth$ca,
      vcmax=vcmax,
      chi_jmax_lim = 0
    ))
  }
  
  
}
