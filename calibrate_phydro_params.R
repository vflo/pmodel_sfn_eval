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

fn_calib <- function(par, data_temp, stomatal_model){
  K.scalar <- par[1]
  P50 <- -1*par[2]
  gamma <- par[3]
  alpha <- par[4]
  print(paste("K.scalar: ", K.scalar, "; P50: ", P50, "; gamma: ", gamma,
              "; alpha: ", alpha))
  
  
  # if(K.scalar<=0 | P50 >=0 | gamma<=0 | alpha <= 0|
  #    K.scalar>1e4 | P50 < -10 | gamma>10 | alpha > 0.4){ 
  #   print(paste("Over-limits"))
  #   return(1e6)
  # 
  #   }else{
  if(unique(data_temp$Species) == "Helianthus annuus"){b = 1.4}else{b=1}
  
  par_plant_now = list(
    conductivity = K.scalar*1e-16, # in Joshi et al 2021, K is in  10^-16 m 
    psi50 = P50,  
    b = b
  )

  par_cost_now = list(
    alpha  = alpha, # cost of Jmax
    gamma = gamma             
  )
  
  if (cpar[which(cpar$Species == species), "inst"] == 0){
    cat("Acclimated response: days = ", unique(data_temp$Drydown.days), "\n")
    calib_phydro = tibble(var = data_temp$Predawn.LWP..MPa.) %>% 
      mutate(pmod = map(var, ~ rphydro_analytical(tc = mean(data_temp$T..deg.C., na.rm = TRUE), 
                                                  ppfd = mean(data_temp$Iabs.in.measurement.chamber..liCor...umol.m.2.s.1., na.rm = TRUE), 
                                                  vpd = mean(data_temp$D..unitless...Pa...Pa.atm..*101325, na.rm = TRUE),
                                                  co2 = mean(data_temp$ca..ppm., na.rm = TRUE), 
                                                  elv = 0, 
                                                  fapar = .99, 
                                                  kphio = 0.087, 
                                                  psi_soil = ., 
                                                  rdark = 0.02, 
                                                  par_plant = par_plant_now, 
                                                  par_cost = par_cost_now))) %>% 
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
                      ~rphydro_analytical(tc = mean(data_temp$T..deg.C., na.rm = TRUE),
                                          ppfd = mean(data_temp$Iabs.growth.in.growth.chamber..umol.m.2.s.1., na.rm = TRUE), 
                                          vpd = mean(data_temp$D..unitless...Pa...Pa.atm..*101325, na.rm = TRUE), 
                                          co2 = mean(data_temp$ca..ppm., na.rm = TRUE), 
                                          elv = 0, 
                                          fapar = .99, 
                                          kphio = 0.087, 
                                          psi_soil = ., 
                                          rdark = 0.02, 
                                          par_plant = par_plant_now, 
                                          par_cost = par_cost_now))) %>% 
    unnest_wider(pmod)
  
  calib_phydro = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>% 
    mutate(pmod = purrr::pmap(list(var, jmax_a, vcmax_a),
                      ~ rphydro_instantaneous_analytical(jmax = ..2,
                                                         vcmax = ..3, 
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
                                                par_cost = par_cost_now))) %>% 
    unnest_wider(pmod)
  
  }
  
  data_temp1 <- data_temp %>%
    left_join(calib_phydro %>% dplyr::rename('Predawn.LWP..MPa.'=var),
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
    # }
  
}



calibrate_phydro_parameters <- function(data_temp, cpar, stomatal_model){
  print(unique(data_temp$Species))
  if (any(cpar$Species == unique(data_temp$Species))){
    init <- cpar[which(cpar$Species == unique(data_temp$Species)),]
    # if(is.na(init$K.scalar)){init$K.scalar <- 1}
    # if(is.na(init$P50)){init$P50 <- -1}
    # if(is.na(init$b)){init$b <- 1}
    # if(is.na(init$alpha)){init$alpha <- 0.1}
  
  
    param <- optimr::optimr(
      par = c(init$K.scalar, -init$P50, init$gamma, init$alpha),
      # par = c(1,1,1,0.1),
      fn = fn_calib,
      data_temp = data_temp,
      # K.scalar = init$K.scalar,
      # P50 = init$P50,
      stomatal_model = stomatal_model,
      # method = "BFGS",
      control = list(maxit = 500, fnscale = 1e2,
                     all.methods = TRUE)
    )  

  return(tibble(Species = unique(data_temp$Species),
                K.scalar = param$par[1],
                P50 = -1* param$par[2],
                gamma = param$par[3],
                alpha = param$par[4]))
  }else{
    return(tibble(Species = unique(data_temp$Species),
                  K.scalar = NA,
                  P50 = NA,
                  gamma = NA,
                  alpha = NA))
}

}

cpar_phydro <- dat1 %>% 
  # filter(!is.na(K.scalar),!is.na(P50),!is.na(b),!is.na(alpha),) %>% 
  group_by(Species) %>% 
  dplyr::group_split()%>% 
  purrr::map_df(~calibrate_phydro_parameters(., cpar,"phydro")) %>% 
  dplyr::bind_cols()


save(cpar_phydro, file = "DATA/cpar_phydro_victor.RData")
