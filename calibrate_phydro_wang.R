library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)
library(rphydro)
source("stomatal_optimization_functions_phydro_calibration.R")
source('hydraulic_functions.R')
source('photosynthetic_functions.R')
source("QUADP.R")

dat = read.csv("DATA/drying_experiments_meta-analysis_Joshi_et_al_2022.csv")
dat$obs_id = 1:nrow(dat)
dat$Ciest = dat$ca..ppm. - dat$A..umol.m.2.s.1./dat$gC..mol.m.2.s.1.
dat$a_pred = NA
dat$g_pred = NA
dat$c_pred = NA

dpsi_df = read.csv("DATA/drying_experiments_dpsi_Joshi_et_al_2022.csv")
dpsi_df <- dpsi_df %>% dplyr::rename(Predawn.LWP..MPa. = SWP)

dat1 <- dat %>%  full_join(dpsi_df)

cpar = read.csv("DATA/fitted_params_Joshi_et_al_2022.csv") 


optimise_stomata_phydro_schemes_calib <- function(fn_profit, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                                            par_plant, par_env, jmax_lim, return_all = FALSE, do_optim=TRUE, stomatal_model){
  
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
    control        = list(maxit = 500,maximize = TRUE,  fnscale = 1e2)
  )
  
  
  if (out_optim$convergence == 9999){out_optim$par <- c(0,0)}
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}


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
  
  #Phydro/phydro_cmax/phydro_sperry/phydro_sox
    K      = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (ground) s-1 MPa-1
    p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
    e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1
    
    # 3. Optimizer
    lj_dps = optimise_stomata_phydro_schemes_calib(fn_profit, 
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


fn_calib <- function(par, data_temp, stomatal_model){
  K.scalar <- par[1]
  P50 <- -1*par[2]
  b <- par[3]
  alpha <- par[4]
  print(paste("K.scalar: ", K.scalar, "; P50: ", P50, "; b: ", b,
              "; alpha: ", alpha))
  
  
  if(K.scalar<=0 | P50 >=0 | b<1 | alpha <= 0|
     K.scalar>1e4 | P50 < -10 | b>20 | alpha > 0.4){ 
    print(paste("Over-limits"))
    return(1e6)

    }else{
  
  par_plant_now = list(
    conductivity = K.scalar*1e-16, # in Joshi et al 2021, K is in  10^-16 m 
    psi50 = P50,  
    b = b
  )
  
  par_cost_now = list(
    alpha  = alpha          # cost of Jmax
  )
  
  if (cpar[which(cpar$Species == species), "inst"] == 0){
    cat("Acclimated response: days = ", unique(data_temp$Drydown.days), "\n")
    calib_phydro = tibble(var = data_temp$Predawn.LWP..MPa.) %>% 
      mutate(pmod = map(var, ~ model_numerical(tc = mean(data_temp$T..deg.C., na.rm = TRUE), 
                                               ppfd = mean(data_temp$Iabs.in.measurement.chamber..liCor...umol.m.2.s.1., na.rm = TRUE), 
                                               vpd = mean(data_temp$D..unitless...Pa...Pa.atm..*101325, na.rm = TRUE),
                                               co2 = mean(data_temp$ca..ppm., na.rm = TRUE), 
                                               elv = 0, 
                                               fapar = .99, 
                                               kphio = 0.087, 
                                               psi_soil = ., 
                                               rdark = 0.02, 
                                               par_plant = par_plant_now, 
                                               par_cost = par_cost_now,
                                               stomatal_model = stomatal_model))) %>% 
      unnest_wider(pmod)
    
  } else {
    ndays = mean(data_temp$Drydown.days,na.rm=TRUE)
    psi_min = min(data_temp$Predawn.LWP..MPa.,na.rm=TRUE)
    psi_max = 0 #max(data$LWP)
    # cat(ndays,"\n")
    
    # Assuming the rate of (linear) drydown, calc the days required to reach each lwp value in [-6,0]
    # lwp = seq(-6,0, length.out=30)
    lwp = data_temp$Predawn.LWP..MPa.
    day = ndays * (lwp-psi_max)/(psi_min-psi_max)
    
    # Function (interpolate) to back-calculate lwp on any given day
    # Thus lwp_day(day) = lwp
    lwp_day = function(day_num){
      psi_max + day_num/ndays * (psi_min-psi_max)
    } 
    
    # Calculate week-average LWP corresponding to days from 0 to last day (corresponding to -6 MPa). Days before start (-ve days) are considered well-watered
    k = 7
    lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
    
    cat("Instantaneous response with slow acclimation: days = ", unique(data_temp$Drydown.days), ", k = ", k, "\n")
    
    # To back-calculate LWP-avg for the desired instantaneous LWP, we need to use an interpolator. 
    # desired lwp --> decimal day assuming linear drydown ---> week_avg lwp for decimal day
    spl = splinefun(x = max(day):0, y=lwp_week)
    
    dat_acc = tibble(var = spl(day)) %>% 
      mutate(pmod = map(var,
                        ~model_numerical(
                          tc = mean(data_temp$T..deg.C., na.rm = TRUE),
                          ppfd = mean(data_temp$Iabs.growth.in.growth.chamber..umol.m.2.s.1., na.rm = TRUE), 
                          vpd = mean(data_temp$D..unitless...Pa...Pa.atm..*101325, na.rm = TRUE), 
                          co2 = mean(data_temp$ca..ppm., na.rm = TRUE), 
                          elv = 0, 
                          fapar = .99, 
                          kphio = 0.087, 
                          psi_soil = ., 
                          rdark = 0.02, 
                          par_plant = par_plant_now, 
                          par_cost = par_cost_now,
                          stomatal_model = stomatal_model))) %>% 
      unnest_wider(pmod)
    
  calib_wang = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>% 
    mutate(pmod = purrr::pmap(list(var, jmax_a, vcmax_a),
                              ~ model_numerical_instantaneous(
                                  vcmax = ..3,
                                  jmax = ..2,
                                  tc = mean(data_temp$T..deg.C., na.rm = TRUE), 
                                  ppfd = mean(data_temp$Iabs.in.measurement.chamber..liCor...umol.m.2.s.1., na.rm = TRUE), 
                                  vpd = mean(data_temp$D..unitless...Pa...Pa.atm..*101325, na.rm = TRUE),
                                  co2 = mean(data_temp$ca..ppm., na.rm = TRUE),
                                  elv = 0, 
                                  fapar = .99,
                                  kphio = 0.087,
                                  psi_soil = ..1, 
                                  rdark = 0.02,
                                  par_plant = par_plant_now, 
                                  par_cost = par_cost_now,
                                  stomatal_model = stomatal_model))) %>% 
    unnest_wider(pmod)
  }
  
  data_temp1 <- data_temp %>%
    left_join(calib_wang %>% dplyr::rename('Predawn.LWP..MPa.'=var),
              by = 'Predawn.LWP..MPa.')%>% 
    mutate(Chi = Ciest/ca..ppm.)
  
  
  Error <- data_temp1 %>% 
      mutate(a_act = A..umol.m.2.s.1.,
             a_act_mean = mean(a_act, na.rm = TRUE),
             a_pred = a,
             gs_act = gC..mol.m.2.s.1.,
             gs_act_mean = mean(gs_act, na.rm = TRUE),
             gs_pred = gs,
             chi_act = Chi,
             chi_act_mean = mean(chi_act, na.rm = TRUE),
             chi_pred = chi, 
             dpsi_act = Dpsi,
             dpsi_act_mean = mean(dpsi_act, na.rm = TRUE),
             dpsi_pred = dpsi
             ) %>% 
      summarise(E_a = sum(((a_pred - a_act)/a_act_mean)^2, na.rm = TRUE),
                E_gs = sum(((gs_pred - gs_act)/gs_act_mean)^2, na.rm = TRUE),
                E_chi = sum(((chi_pred - chi_act)/chi_act_mean)^2, na.rm = TRUE),
                E_dpsi = sum(((dpsi_pred - dpsi_act)/dpsi_act_mean)^2, na.rm = TRUE),
                E_dpsi = case_when(is.na(E_dpsi)~0,
                                   TRUE~E_dpsi),
                Error = E_a+E_gs+E_chi+E_dpsi
                ) %>% 
    dplyr::select(Error)
      

  # print(paste("RMSE A:",rmse_a))
  # print(paste("RMSE gs:",rmse_gs))
  print(paste("Total Error:", Error$Error))

  return(Error$Error)
    }
  
}



calibrate_wang_parameters <- function(data_temp, cpar, stomatal_model){
  print(unique(data_temp$Species))
  if (any(cpar$Species == unique(data_temp$Species))){
    init <- cpar[which(cpar$Species == unique(data_temp$Species)),]
    if(is.na(init$K.scalar)){init$K.scalar <- 1}
    if(is.na(init$P50)){init$P50 <- -1}
    if(is.na(init$b)){init$b <- 1}
    if(is.na(init$alpha)){init$alpha <- 0.1}
  
  
    param <- optimr::optimr(
      par = c(init$K.scalar, -init$P50, init$b, init$alpha),
      fn = fn_calib,
      data_temp = data_temp,
      # K.scalar = init$K.scalar,
      # P50 = init$P50,
      stomatal_model = stomatal_model,
      # method = "CG",
      control = list(maxit = 10000,
                     all.methods = TRUE)
    )  

  return(tibble(Species = unique(data_temp$Species),
                K.scalar = param$par[1],
                P50 = -1* param$par[2],
                b = param$par[3],
                alpha = param$par[4]))
  }else{
    return(tibble(Species = unique(data_temp$Species),
                  K.scalar = NA,
                  P50 = NA,
                  b = NA,
                  alpha = NA))
}

}

cpar_wang <- dat1 %>% 
  # filter(!is.na(K.scalar),!is.na(P50),!is.na(b),!is.na(alpha),) %>% 
  group_by(Species) %>% 
  dplyr::group_split()%>% 
  purrr::map_df(~calibrate_wang_parameters(., cpar,"phydro_wang")) %>% 
  dplyr::bind_cols()


save(cpar_wang, file = "DATA/cpar_wang.RData")
