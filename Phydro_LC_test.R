source("stomatal_optimization_functions.R")
library(tidyverse)
library(rpmodel)

par_plant = list(
  # Ks0=1e-12,              # m2
  # v_huber=1e-4,           #
  height=10,              # m
  # conductivity_scalar=3,  # Scales (Ks0*v_huber/H)
  conductivity = 3e-17,
  psi50 = -3.46,             # Mpa
  d = -3.66,
  c = 6.55,
  b=6.55,                    # 
  rs = 0.1,
  rr =5,
  dens_sapwood = 700
)

tc=25
p = rpmodel::calc_patm(0)
vpd = 10000

par_env = list(
  viscosity_water = rpmodel::calc_viscosity_h2o(tc, p),  # Needs to be imported from rpmodel.R
  density_water = rpmodel::calc_density_h2o(tc, p),  # Needs to be imported from rpmodel.R
  patm = p,
  tc = tc,
  vpd = vpd,
  LAI = 1.5
)

calc_gs_test <- function (dpsi, psi_soil, par_plant, par_env, ...){
  K = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = FALSE) #mol m-2 (leaf) s-1 MPa-1
  # K = K * par_plant$pl_sapw_area/par_plant$pl_ba*par_plant$st_basal_area*1e-4 #kg m-2 (ground) s-1 MPa-1
  # K = K * par_env$LAI #kg m-2 (leaf) s-1 MPa-1
  D = (par_env$vpd/par_env$patm)
  # D = par_env$vpd
  K/1.6/D * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (leaf) s-1
}

get_std_photosythnesis_params = function(){ 
  ## Set P-model parameters
  kphio <- 0.05        # quantum yield efficiency
  # c_molmass <- 12.0107 # molar mass, g / mol
  
  ## Define environmental conditions
  tc <- 25             # temperature, deg C
  ppfd <- 400          # umol/m2/s
  # vpd  <- 1000         # Pa
  co2  <- 400          # ppm
  elv  <- 0            # m.a.s.l.
  fapar <- 0.7         # fraction
  
  p = rpmodel::calc_patm(elv)
  return (list(
    kmm = rpmodel::calc_kmm(tc, p),  # Why does this use std. atm pressure, and not p(z)?
    gammastar = rpmodel::calc_gammastar(tc, p),
    phi0 = kphio*rpmodel::calc_ftemp_kphio(tc),
    Iabs = ppfd*fapar,
    ca = co2*p*1e-6,  # Convert to partial pressure
    patm = p,
    delta = 0.00
  ))
}

calc_gs_LC <- function (dpsi, psi_soil, par_plant, par_env, ...){
  K = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = FALSE) #mol m-2 (leaf) s-1 MPa-1
  # K = K * par_plant$pl_sapw_area/par_plant$pl_ba*par_plant$st_basal_area*1e-4 #kg m-2 (ground) s-1 MPa-1
  # K = K * par_env$LAI #kg m-2 (leaf) s-1 MPa-1
  D = (par_env$vpd/par_env$patm)
  # D = par_env$vpd
  K/1.6/D * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (leaf) s-1
}


###########################################
## PHYDRO least cost
###########################################
# 
# 
# 
# optimise_stomata_phydro_LC <- function(fn_profit, psi_soil, par_cost, par_photosynth, par_plant, par_env, return_all = FALSE, opt_hypothesis){
#   
#   out <- capture.output(out_optim <- optimr::optimr(
#     par            = c(logjmax=0, dpsi=1),  
#     lower          = c(-10, .0001),
#     upper          = c(10, 1e6),
#     fn             = fn_profit_LC,
#     psi_soil       = psi_soil,
#     par_cost       = par_cost,
#     par_photosynth = par_photosynth,
#     par_plant      = par_plant,
#     par_env        = par_env,
#     do_optim       = TRUE,
#     opt_hypothesis = opt_hypothesis,
#     method         = "L-BFGS-B",
#     control        = list(maxit = 500, maximize = TRUE, fnscale = 1e4, trace = TRUE)
#   ))
#   
#   out_optim$value <- -out_optim$value
#   
#   if (return_all){
#     out_optim
#   } else {
#     return(out_optim$par)
#   }
# }
# 
# par_photosynth = get_std_photosythnesis_params()
# 

# psi_soil = -1
# p_crit = par_plant$d * log(1000.0) ^ (1.0/par_plant$c)
# # 2. if soil psi is lower than p_crit, set soil psi to 98% p_crit 
# low_swp <- FALSE
# if(psi_soil <= p_crit){
#   psi_soil = p_crit*0.98
#   low_swp  = TRUE
#   print("WARNING: soil water potential is lower than critical plant water potential.")
# }
# 
# # 3. Optimizer
# lj_dps = optimise_stomata_phydro_LC(fn_profit_LC, 
#                                  psi_soil = psi_soil, 
#                                  par_cost  = par_cost, 
#                                  par_photosynth = par_photosynth, 
#                                  par_plant = par_plant, 
#                                  par_env = par_env, 
#                                  opt_hypothesis = "LC")
# 
# profit = fn_profit(par = lj_dps, psi_soil = psi_soil, par_cost  = par_cost, par_photosynth = par_photosynth, par_plant = par_plant, par_env = par_env, opt_hypothesis = "PM")
# 
# jmax  = exp(lj_dps[1]) %>% unname()
# dpsi  = lj_dps[2] %>% unname()
# gs    = calc_gs_LC(dpsi, psi_soil, par_plant, par_env) * par_env$LAI # gs in mol m-2 (ground) s-1
# a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
# a     = a_j$a
# ci    = a_j$ci
# vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
# E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
# gs    = gs/patm*1e6 #transform to umol m-2(ground) s-1 Pa-1
# psi_l = psi_soil - dpsi
# 
# 



# jmax <- seq(-2,10,0.1)
# pl <- seq(0.1,5,0.1)
# par <- expand.grid(logjmax = jmax, dpsi = pl)
# 
# par %>% 
#   split(seq(nrow(.))) %>% 
#   purrr::map(function(x){
#     fn_profit_LC(x, 
#                  psi_soil = psi_soil, 
#                  par_cost = par_cost, 
#                  par_photosynth = par_photosynth, 
#                  par_plant = par_plant, 
#                  par_env = par_env, 
#                  do_optim = FALSE, 
#                  opt_hypothesis = "LC")->res
#     return(res)
#   }) %>% bind_rows()->res
# 
# df_res <- res %>% cbind(par)
# 
# df_res %>% 
#   ggplot(aes(dpsi,out))+
#   geom_point()
# 
# df_res %>% 
#   ggplot(aes(dpsi,a_cost))+
#   geom_point()
# 
# df_res %>% 
#   ggplot(aes(dpsi,E_cost))+
#   geom_point()
# 
# df_res %>% 
#   ggplot(aes(dpsi,vcmax_cost, color = E_cost))+
#   geom_point()
# 
# df_res %>% 
#   ggplot(aes(dpsi,jmax_cost, color = ))+
#   geom_point()
# 
# df_res[which(df_res$out == df_res$out %>% max(na.rm = TRUE)),]
# 
# 

p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)

par_plant = list(
  # Ks0=1e-12,              # m2
  # v_huber=1e-4,           #
  height=10,              # m
  # conductivity_scalar=3,  # Scales (Ks0*v_huber/H)
  conductivity = 3e-17,
  psi50 = -3.46,             # Mpa
  d = -3.66,
  c = 6.55,
  b=6.55,                    # 
  rs = 0.1*1e-6,
  rr =2,
  dens_sapwood = 700
)

par_cost = list(
  b_cost = 0.031,       # cost of vcmax
  c_cost = 0.41*1e-1          # cost of jmax
)


fn_profit_LC <- function(par, psi_soil, par_cost, par_photosynth, par_plant, par_env, do_optim = FALSE, opt_hypothesis){
  
  jmax = exp(par[1]$logjmax)  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]$dpsi      # delta Psi in MPa
  
  gs = calc_gs_LC(dpsi, psi_soil, par_plant, par_env) * par_env$LAI  # gs in mol/m2ground/s
  E  = 1.6*gs*(par_env$vpd/par_env$patm)        # E in umol/m2/s
  h = par_plant$height
  dw = par_env$density_water
  dsw = par_plant$dens_sapwood
  rs = par_plant$rs
  rr = par_plant$rr
  nu = par_env$viscosity_water
  d = par_plant$d
  c = par_plant$c
  K = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = FALSE)
  ks = par_plant$conductivity
  p = par_plant$psi50
  b = par_plant$b
  ecrit = get_e_crit(psi_soil,K,d,c,h,par_env$density_water)
  # a_cost = (rs+(rr*(1-(1/2)^(((psi_soil-dpsi)/p)^b))))/(E) 
  a_cost = (rs/E)+(rr*(1-(1/2)^(((psi_soil-dpsi)/p)^b)))/(ecrit) 
  
  
  ## light-limited assimilation
  a_j   <- calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2/s
  a     = a_j$a
  ci    = a_j$ci
  vcmax = calc_vcmax_coordinated_numerical(a,ci, par_photosynth)
  
  costs = par_cost$alpha * jmax + par_cost$gamma * dpsi^2 #((abs((-dpsi)/par_plant$psi50)))^2  
  
  benefit = 1 #(1+1/(par_photosynth$ca/40.53))/2
  
  dummy_costs = 0*exp(20*(-abs(dpsi/4)-abs(jmax/1))) # ONLY added near (0,0) for numerical stability.
  
  if (opt_hypothesis == "PM"){
    ## Profit Maximisation
    out <- a*benefit - costs - dummy_costs
  } else if (opt_hypothesis == "LC"){
    ## Least Cost
    # b_cost = a_cost * 146
    b_cost = par_cost$b_cost
    out <- tibble(out_MP = a - (a_cost * E + b_cost * vcmax + par_cost$c_cost * jmax ),
                  # out_LC = (a_cost * E + b_cost * vcmax + par_cost$c_cost * jmax )/a,
                  a_cost = a_cost,
                  E_cost =  a_cost * E,
                  vcmax_cost = b_cost * vcmax,
                  jmax_cost = par_cost$c_cost*jmax,
                  vcmax = vcmax,
                  jmax = jmax,
                  a = a,
                  E = E,
                  gs = gs,
                  psi_soil = psi_soil)
  }
  
  
  return(out)
}

par_cost = list(
  alpha = 0.1     # cost of jmax
)

par_photosynth = get_std_photosythnesis_params()

fn_profit_phydro_wang <- function(par, psi_soil, e_crit, par_cost, par_photosynth, par_plant, par_env, do_optim = FALSE, opt_hypothesis){
  
  jmax = exp(par[1])$logjmax  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]$dpsi       # delta Psi in MPa
  K = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = FALSE)
  gs = calc_gs_test(dpsi, psi_soil, par_plant, par_env) * par_env$LAI  # gs in mol/m2ground/s
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  # e = K*dpsi
  # gs_crit <- e_crit/(1.6*(par_env$vpd/par_env$patm))
  ## light-limited assimilation
  a_j   <- calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
  a     = a_j$a
  ci    = a_j$ci
  vcmax = calc_vcmax_coordinated_numerical(a,ci, par_photosynth)
  
  # a_j_crit   <- calc_assim_light_limited(gs_crit, jmax, par_photosynth) # Aj in umol/m2ground/s
  # a_crit     = a_j$a
  # ci_crit    = a_j$ci
  # vcmax_crit = calc_vcmax_coordinated_numerical(a_crit,ci_crit, par_photosynth)
  
  cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
  # out <- jmax - (a*(e_crit*1e6-e*1e6) - costs)/(e_crit*1e6) #umolco2 m-2 s-1
  # out = (a*(e_crit-e) - costs*(e_crit-e))/(e_crit), #umolco2 m-2 s-1
  
  out <- tibble(out_ecrit = (a-cost),
                out_e = (a-cost)*e/e_crit,
                out = out_ecrit-out_e, #umolco2 m-2 s-1
                vcmax = vcmax,
                jmax = jmax,
                a = a,
                E = e,
                gs = gs,
                psi_soil = psi_soil)
  return(out)
}


tc=25
p = rpmodel::calc_patm(0)
vpd = 100
psi_soil = -0.01

par_env = list(
  viscosity_water = rpmodel::calc_viscosity_h2o(tc, p),  # Needs to be imported from rpmodel.R
  density_water = rpmodel::calc_density_h2o(tc, p),  # Needs to be imported from rpmodel.R
  patm = p,
  tc = tc,
  vpd = vpd,
  LAI = 1.5
)

LAI        = par_env$LAI
K          = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = FALSE)*LAI #mol m-2 (ground) s-1 MPa-1
d          = par_plant$d
c          = par_plant$c
h          = par_plant$height
dens_water = par_env$density_water
p_crit     = d * log(1000.0) ^ (1.0/c)
e_crit = get_e_crit(psi_soil, K, d, c, h, dens_water) #mol m-2 (ground) s-1

jmax <- 4
pl <- seq(0.02,5,0.02)
par <- expand.grid(logjmax = jmax, dpsi = pl)

par %>% 
  split(seq(nrow(.))) %>% 
  purrr::map_df(function(x){
    fn_profit_phydro_wang(x,
                 psi_soil = psi_soil, 
                 e_crit = e_crit,
                 par_cost = par_cost, 
                 par_photosynth = par_photosynth, 
                 par_plant = par_plant, 
                 par_env = par_env, 
                 do_optim = FALSE, 
                 opt_hypothesis = "LC")->res
    return(res)
  }) ->res

df_res <- res %>% cbind(par)

df_res[which(df_res$out == df_res$out %>% max(na.rm = TRUE)),]
# df_res[which(df_res$out_LC == df_res$out_LC %>% min(na.rm = TRUE)),]

df_res %>% 
  ggplot()+
  geom_point(aes(dpsi,out_ecrit),color = "green")+
  geom_point(aes(dpsi,out_e),color = "purple")+
  geom_point(aes(dpsi,out, color = logjmax))+
  # ylim(1,3 %>% )
  NULL

df_res %>% 
  ggplot()+
  geom_point(aes(gs,out_ecrit),color = "green")+
  geom_point(aes(gs,out_e),color = "purple")+
  geom_point(aes(gs,out, color = logjmax))+
  # ylim(1,3 %>% )
  NULL

df_res %>% 
  ggplot(aes(dpsi,out, color = vcmax))+
  geom_point()+
  # ylim(1,3 %>% )
  NULL

df_res %>% 
  ggplot(aes(dpsi,out, color = E))+
  geom_point()+
  # ylim(1,3 %>% )
  NULL


df_res %>% 
  ggplot(aes(dpsi,out, color = gs))+
  geom_point()+
  # ylim(1,3 %>% )
  NULL
# 
# df_res %>% 
#   ggplot(aes(dpsi,a_cost))+
#   geom_point()
# 
# df_res %>% 
#   ggplot(aes(dpsi,E_cost))+
#   geom_point()
# 
# df_res %>% 
#   ggplot(aes(dpsi,vcmax_cost, color = E_cost))+
#   geom_point()
# 
# df_res %>% 
#   ggplot(aes(dpsi,jmax_cost ))+
#   geom_point()
# 
# df_res %>% 
#   ggplot(aes(vcmax,exp(logjmax), color = a))+geom_point()
# 

