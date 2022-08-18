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

# 
# model_numerical <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, 
#                             rdark, par_plant, par_cost, stomatal_model){
#   
#   patm = rpmodel::calc_patm(elv)
#   par_photosynth <- list(
#     kmm       = rpmodel::calc_kmm(tc, patm),
#     gammastar = rpmodel::calc_gammastar(tc, patm),
#     phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
#     Iabs      = ppfd * fapar,
#     ca        = co2 * patm * 1e-6,  # Convert to partial pressure
#     patm      = patm,
#     delta     = rdark
#     # delta     = ftemp_inst_rd( tc )
#   )
#   par_env = list(
#     viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
#     density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
#     patm            = patm,
#     tc              = tc,
#     vpd             = vpd
#   )
#   
#   #Phydro/phydro_cmax/phydro_sperry/phydro_sox
#     K      = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (ground) s-1 MPa-1
#     p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
#     e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1
#     
#     # 3. Optimizer
#     lj_dps = optimise_stomata_phydro_schemes_calib(fn_profit, 
#                                              psi_soil = psi_soil,
#                                              e_crit = e_crit,
#                                              p_crit = p_crit,
#                                              par_cost  = par_cost, 
#                                              par_photosynth = par_photosynth, 
#                                              par_plant = par_plant, 
#                                              par_env = par_env,
#                                              jmax_lim = 7, 
#                                              return_all = FALSE, 
#                                              do_optim=TRUE, 
#                                              stomatal_model = stomatal_model)
#     
#     
#     jmax  = exp(lj_dps[1]) %>% unname()
#     dpsi  = lj_dps[2] %>% unname()
#     psi_l = psi_soil-dpsi
#     gs    = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env) # gs in mol m-2 (ground) s-1
#     a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
#     a     = a_j$a
#     ci    = a_j$ci
#     vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
#     E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
#     gs    = gs #transform to umol m-2(ground) s-1 Pa-1
#     
#     return(list(
#       jmax         = jmax,
#       dpsi         = dpsi,
#       p_leaf       = psi_l,
#       gs           = gs,
#       E            = E,
#       a            = a,
#       ci           = ci,
#       chi          = ci/par_photosynth$ca,
#       vcmax        = vcmax,
#       chi_jmax_lim = 0
#     ))
#   }


fn_calib <- function(par, data_temp, stomatal_model){
  K.scalar <- par[1]
  P50 <- -1*par[2]
  b <- par[3]
  alpha <- par[4]
  print(paste("K.scalar: ", K.scalar, "; P50: ", P50, "; b: ", b,
              "; alpha: ", alpha))
  
  
  if(K.scalar<=0 | P50 >=0 | b<1 | alpha <= 0|
     K.scalar>1e4 | P50 < -10 | b>3 | alpha > 0.4){ 
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
                Error = E_a+E_gs#+E_chi#+E_dpsi
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
      control = list(maxit = 500, fnscale = 1e2,
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








################## Jaideep optimization code ###################################


rm(list=ls())

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)
library(stringr)
library(rphydro)
source("stomatal_optimization_functions_phydro_calibration.R")
source('hydraulic_functions.R')
source('photosynthetic_functions.R')
source("QUADP.R")

dat = read.csv(file="DATA/drying_experiments_meta-analysis_Joshi_et_al_2022.csv")
names = colnames(dat)
# [4] "Predawn.LWP..MPa."                                 
# [5] "A..umol.m.2.s.1."                                  
# [6] "gC..mol.m.2.s.1."                                  
# [7] "ca..ppm."                                          
# [8] "D..unitless...Pa...Pa.atm.."                       
# [9] "T..deg.C."                                                          
names[4:9] = c("LWP", "A", "gC", "ca", "D", "T")
names[13] = "Iabs_growth"
names[12] = "Iabs_used"
colnames(dat) = names

dpsi_df = read.csv(file = "DATA/drying_experiments_dpsi_Joshi_et_al_2022.csv")
par_data = read.csv("DATA/fitted_params_Joshi_et_al_2022.csv") 


plot_all = function(df_w_vol, varname, species, data, dpsi_data=NULL, analytical=F){
  # df_w_vol = df_w_vol[complete.cases(df_w_vol),]
  # View(df_w_vol)
  
  gx = log(df_w_vol$gs+1e-20)
  gy = df_w_vol$var
  # gx1 = seq(max(min(gx), -20), max(gx), length.out=100)
  f = splinefun(x = gx, y=gy)
  gs0 = df_w_vol$gs[which(df_w_vol$var==0)]
  psi88S = f(log(gs0*0.12))
  psi50S = f(log(gs0*0.50))
  psi12S = f(log(gs0*0.88))
  cat("psi88S = ", psi88S, "\n")
  cat("psi50S = ", psi50S, "\n")
  cat("psi12S = ", psi50S, "\n")
  
  dpx = df_w_vol$var
  dpy = df_w_vol$dpsi
  f1 = splinefun(dpy~dpx)
  dpx1 = seq(min(dpx), max(dpx), length.out=100)
  dp88S = f1(psi88S)
  dp50S = f1(psi50S)
  dp12S = f1(psi12S)
  cat("psiL88S = ", psi88S-dp88S, "\n")
  cat("psiL50S = ", psi50S-dp50S, "\n")
  cat("psiL12S = ", psi12S-dp12S, "\n")
  
  cat(psi88S, "\t", psi50S, "\t", psi12S, "\t", psi88S-dp88S, "\t", psi50S-dp50S, psi12S-dp12S, "\n")
  
  subdata = data %>% filter(LWP > -6)
  
  p1 <- df_w_vol %>% 
    ggplot() +
    geom_line(aes(x = var, y = vcmax), col="green3", size=1) +
    # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Vcmax))+
    # geom_vline(xintercept = par_plant_std$psi50, col="orange") +
    expand_limits(y=0)+
    xlab(varname)
  if (analytical) p1 = p1 + geom_line(aes(x = var, y = vcmax ), col="grey", size=1) 
  
  p2 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = dpsi), col="blue", size=1)+
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    expand_limits(y=0)+
    xlab(varname)
  if (!is.null(dpsi_data)){
    cat("Adding points...\n")
    p2 <- p2 + 
      geom_point(data=dpsi_data, aes(x=SWP, y=Dpsi))
  }
  
  p3 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = gs), col="cyan2", size=1)+
    geom_point(data=subdata, aes(x=LWP, y=gC))+
    # geom_vline(xintercept = psi88S, col="grey")+
    expand_limits(y=0)+
    xlab(varname)
  if (analytical) p3 = p3 + geom_line(aes(x = var, y = gs), col="grey", size=1)
  
  p4 <- df_w_vol %>%
    # mutate(chi = ci/out_hydraulics_ca) %>% 
    ggplot() +
    geom_line(aes(x = var, y = chi), col="magenta", size=1)+
    geom_point(data=subdata, aes(x=LWP, y=1-A/gC/ca))+
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    expand_limits(y=0)+
    xlab(varname)
  if (analytical) p4 = p4 + geom_line(aes(x = var, y = chi), col="grey", size=1)
  
  p5 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = jmax), col="goldenrod1", size=1) +
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Jmax))+
    expand_limits(y=0)+
    xlab(varname)
  
  
  p6 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = a), col="green4", size=1) +
    geom_point(data=subdata, aes(x=LWP, y=A))+
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    expand_limits(y=0)+
    xlab(varname)
  if (analytical) p6 = p6 + geom_line(aes(x = var, y = gpp), col="grey",size=1) 
  
  grid.arrange(p3,p6,p1,p2,p5,p4, ncol=2)
}



error_fun = function(x, data, dpsi_file="", plot=F, scale = 1, 
                     dpsi_calib=T, inst=F, k=7, stomatal_model = stomatal_model){
  x = x*scale
  # cat("start: ",x, "\n")
  
  data$Ciest = data$ca-data$A/data$gC
  
  par_plant_now = list(
    conductivity= x[1]*1e-16, #.7, 
    psi50 = x[2], #-.5, 
    b= x[3] #1.2
  )
  
  par_cost_now = list(
    alpha  = x[4], #0.13,       # cost of Jmax
    gamma = x[5]
  )
  if(x[1]<=0 | x[2] >=0 | x[3]<1 | x[4] <= 0| #x[5] < 1|
     x[1]>1e4 | x[2] < -15 | x[3]>10 | x[4] > 0.4 #|x[5] > 1
     ){ 
    print(paste("Over-limits"))
    return(1e6)
    
  }else{
  if (dpsi_file != ""){
    dpsi_data = dpsi_df %>% filter(Species == dpsi_file) #read.csv(paste0("drying_experiments_meta/",dpsi_file,".csv"))
  }
  else{
    dpsi_data=NULL
  }
  
  if (inst == F){
    cat("Acc resp\n")
    lwp = seq(min(data$LWP), 0, length.out=20)
    dat1 = tibble(var = lwp) %>%
      mutate(p = purrr::map(var,
                            ~model_numerical(tc = mean(data$T), 
                                             ppfd = mean(data$Iabs_growth), 
                                             vpd = mean(data$D*101325), 
                                             co2 = mean(data$ca), elv = 0,
                                             fapar = .99, kphio = 0.087, psi_soil = ., 
                                             rdark = 0.02, par_plant=par_plant_now, 
                                             par_cost = par_cost_now,stomatal_model = stomatal_model)) ) %>% 
      unnest_wider(p)
  } else{
    cat("inst resp\n")
    ndays = mean(data$Drydown.days)
    psi_min = -5 #min(data$LWP)
    # psi_min = psi_crit
    psi_max = 0 #max(data$LWP)
    # cat(ndays,"\n")
    
    lwp = seq(psi_min,0, length.out=10)
    day = ndays * (lwp-psi_max)/(psi_min-psi_max)
    
    lwp_day = function(day_num){
      psi_max + day_num/ndays * (psi_min-psi_max)
    }
    
    # k = 7
    lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
    
    spl = splinefun(x = max(day):0, y=lwp_week)
    
    # plot(x=day, y=lwp_day(day))
    # points(x=day, y=spl(day), col="red", type="o")
    # points(x=max(day):0, y=lwp_week, col="blue", type="l")
    
    dat_acc = tibble(var = spl(day)) %>% 
      mutate(pmod = map(var, ~model_numerical(tc = mean(data$T), ppfd = mean(data$Iabs_growth), 
                                              vpd = mean(data$D*101325), co2 = mean(data$ca), 
                                              elv = 0, fapar = .99, kphio = 0.087, 
                                              psi_soil = ., rdark = 0.02, par_plant=par_plant_now, 
                                              par_cost = par_cost_now, stomatal_model = stomatal_model))) %>% 
      unnest_wider(pmod)
    dat1 = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>%
      mutate(p = purrr::pmap(list(var, jmax_a, vcmax_a), 
                             ~model_numerical_instantaneous(tc = mean(data$T), 
                                                            ppfd = mean(data$Iabs_used), 
                                                            vpd = mean(data$D*101325), 
                                                            co2 = mean(data$ca), elv = 0, 
                                                            fapar = .99, kphio = 0.087, 
                                                            psi_soil = ..1, rdark = 0.02, 
                                                            par_plant=par_plant_now, 
                                                            par_cost = par_cost_now, 
                                                            jmax = ..2, vcmax = ..3, 
                                                            stomatal_model = stomatal_model)) ) %>% 
      unnest_wider(p)
  }
  
  # if(plot==T) dat_acc %>% plot_all(varname = "psi_soil", species=species, data = data, dpsi_data=dpsi_data, analytical = F)
  if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=species, data = data, dpsi_data=dpsi_data)
  
  a_spl = splinefun(x = lwp, y=dat1$a)
  g_spl = splinefun(x = lwp, y=dat1$gs)
  c_spl = splinefun(x = lwp, y=dat1$chi)
  
  w = rep(1, length(data$LWP))
  # w[data$LWP < -1.75] = 0
  y1 = mean((a_spl(data$LWP) - data$A)^2*w,na.rm  = TRUE)/mean(data$A,na.rm  = TRUE)^2
  y2 = mean((g_spl(data$LWP) - data$gC)^2*w,na.rm  = TRUE)/mean(data$gC,na.rm  = TRUE)^2 #*4000
  
  # w = exp(-data$LWP/par_plant_now$psi50)
  y4 = mean((c_spl(data$LWP) - data$Ciest/data$ca)^2*w,na.rm  = TRUE)/mean(data$Ciest/data$ca,na.rm  = TRUE)^2 #*800
  cat("c_spl:", c_spl(data$LWP), "\n")
  
  if (dpsi_calib == T){
    d_spl = splinefun(lwp, y=dat1$dpsi)
    w = 1 #0.1 + exp(dpsi_data$SWP)
    # dpsi_swp <- dpsi_data %>% filter(SWP>=-1)
    y3 = mean((d_spl(dpsi_data$SWP) - dpsi_data$Dpsi)^2*w,na.rm  = TRUE)/mean(dpsi_data$Dpsi,na.rm  = TRUE)^2 #*40
    # y3 = mean((d_spl(dpsi_swp$SWP) - dpsi_swp$Dpsi)^2*w)/mean(dpsi_swp$Dpsi)^2 #*40
    cat("d_spl:", d_spl(dpsi_data$SWP), "\n")
  }else{
    # w=1
    # p50 = -1.5 #par_plant_now$psi50
    # y3 = mean((dat1$dpsi[lwp>p50] - -p50)^2*w)/mean(dat1$dpsi[lwp>p50])^2 #*20
    y3=0
  }
  y = (y1+y2+y3+y4)*100
  cat(x, "|", y1, "/ ", y2, " / ", y3, " / ", y4, " / ", y, "\n")
  cat(x, "|", y, "\n")
  
  
  x0 <<- x
  y
  }
}




stomatal_model <- "phydro_cgain"

species = "Eucalyptus pilularis"
data1=filter(dat, Species==species)

# data1 = data1 %>% arrange(desc(LWP)) %>% mutate(day = abs(LWP)/(max(abs(LWP)))*Drydown.days )
# data1 = data1 %>% filter(LWP > -1.5)

x0 = as.numeric((par_data %>% filter(Species == species) )[2:6])
x0 = c(x0[1], x0[2],x0[3], x0[4], 1e6)

# error_fun(x0, data1,plot=T, dpsi_file=species, dpsi_calib = T, inst=T, stomatal_model = stomatal_model)


# x0= c(9, -0.655592, 1.2, 0.1070148)
# x0 = c(1.545316, -1.826197, 1.70648, 0.1028241 )
error_fun(x0, data1 %>% mutate(Drydown.days=38), plot=T, dpsi_file=species, dpsi_calib = T, inst=T, stomatal_model = stomatal_model)

optimr::optimr(par = x0/abs(x0),
               fn = error_fun,
               data=data1 %>% mutate(Drydown.days=38),
               scale=abs(x0),
               dpsi_file = species,  # Set to species name or "" if dpsi data is not available
               dpsi_calib = T,       # Set to F if dpsi data is not available
               inst=T,
               stomatal_model = stomatal_model,
               control = list(par_scale = 1, maxit=700, all.methods = TRUE)
)
# error_fun(x0, data1, plot=T, dpsi_file="")

error_fun(x0, data1,#%>% mutate(Drydown.days=38), 
          plot=T, dpsi_file=species, inst=F,
          stomatal_model = stomatal_model)


### Cross Validation
species_master = unique(dat$Species) %>% str_subset(pattern = "Mediterranean", negate = T) %>% str_subset("Glycine max2", negate=T)
res_master = data.frame()

K=5

species = species_master[1]
species_short = paste(substr(strsplit(species, " ")[[1]][-(5:8)], 0,2), collapse = ".")
data1=filter(dat, Species==species)

data1 = data1[sample(1:nrow(data1), replace = F), ]
N = ceiling(nrow(data1)/K)
if (nrow(data1) < 2*K-1) N=1
start_ids = seq(1, nrow(data1), N)
start_ids = c(start_ids, nrow(data1)+1)
cat("splits = ", start_ids, "\n")


ids_list = list()
for (i in 1:K){
  ids_list[[i]] = start_ids[i]:(start_ids[i+1]-1)
}

for (iter in 1:5){
  ids_test = unlist(ids_list[iter])
  ids_train = unlist(ids_list[-iter])
  
  dat_train = data1[ids_train,]
  dat_test  = data1[ids_test,]
  
  x0 = as.numeric((par_data %>% filter(Species == species) )[2:6])
  names(x0) = c("K.scalar", "P50", "b", "alpha" ,"yxp50")
  error_fun(x0, dat_train, plot=F, dpsi_file="")
  
  out = optimr::optimr(par = x0,
                       fn = error_fun,
                       data=dat_train,
                       # dpsi_file = species
                       control = list(par_scale = scale, maxit=100)
  )
  
  res = data.frame(as.list(out$par))
  res$height = mean(dat_train$Ht_used)
  res$species = species
  res$iter = iter
  res$err_train = error_fun(out$par, dat_train, plot=F, dpsi_file="")
  res$err_test  = error_fun(out$par, dat_test, plot=F, dpsi_file="")
  
  res_master = rbind(res_master, res)
}

write.csv(res_master, file = paste0("cross_validation/cross_validation_", species_short,".csv"))


# error_fun(out$par, dat_train, plot=T, dpsi_file="")
# error_fun(out$par, dat_test, plot=T, dpsi_file="")

# #### Analysis ####
# 
# library(MASS)
# library(tidyverse)
# d = read.csv("C:/Users/Jaideep/Documents/codes/rpmodel/vignettes_add/drying_experiments_meta/fitted_params_brd_all_3_p50pmin.csv")
# # d = d %>% filter(!(Species %in% c("Helianthus annuus", "Glycine max")))
# d = d %>% mutate(gamma = yxp50/P50^2)
# d = d %>% mutate(P88 = (log(0.12)/log(0.5))^(1/b)*P50 )
# d = d %>% mutate(K.scalar = K.scalar/H)
# d = d %>% dplyr::select(-Species, -Ref.1, -Ref2, -Ref, -b..Cavit., -b..SC., -Pclose, -Slope, -P12, -P50..stem., -PgS88, -PgS50, -pg12, -H)
# d = d %>% dplyr::select(-SLA_nopetiole, -SLA_wpetiole, -SLA1, -A.G)
# d = d %>% dplyr::select(-Pmin, -dpsi_int, -P88)
# d = d %>% dplyr::select(-pgL50, -Pgs90, -b)
# d = d %>% mutate(SM_skelton = -P50 + pgL12)
# # d = d %>% mutate(SM_choat = P50 - P88)
# # d1 = sapply(d, FUN = function(x){ifelse(x<0, -x, x)}, simplify = T)
# # d = d %>% mutate(logb = log(b)) %>% dplyr::select(-b)
# # d = d %>% mutate(logY = log(gamma)) %>% dplyr::select(-gamma)
# # d1 = d1 %>% mutate(logb = log(b)) %>% dplyr::select(-b)
# # d = d %>% mutate(logP50 = log(-P50)) %>% dplyr::select(-P50)
# # d = d %>% mutate(logK = log(K.scalar)) %>% dplyr::select(-K.scalar)
# 
# panel_func = function(x,y,...){
#   points(y~x,...);
#   mod = rlm(y~x, weights = rep(1, length(x)) );
#   if (abs(summary(mod)$coefficients[2,3]) > 1.96) {
#     abline(lm(y~x), col="black")
#   }
#   else {
#     if (abs(summary(mod)$coefficients[2,3]) > .9){
#       abline(lm(y~x), col="grey65")
#     }
#   }
# }
# 
# d %>% pairs(panel= panel_func, pch=20, cex=1.2, col="cyan4") #
# 
# d %>% dplyr::select(P88, pgL88, Ptlp) %>% pairs(panel= panel_func)
# 
# plot(d$pgL12~d$P50)
# abline(0,1)
# 
# par(mfrow=c(2,3), mar=c(4,4,1,1), cex.lab=1.2)
# 
# with(d1, plot(y=alpha, x=gamma, xlab="log(gamma)", ylab="log(alpha)"))
# abline(rlm((d1$alpha)~(d1$gamma), method="MM"), col="cyan3")
# 
# with(d1, plot(log(K.scalar)~log(P50)))
# abline(rlm(log(d1$K.scalar)~log(d1$P50), method="MM"), col="cyan3")
# 
# with(d1, plot((P50)~(P88)))
# abline(rlm((P50)~log(P88), method="MM"), col="cyan3")
# 
# with(d1, plot((P50)~(pgL12)))
# abline(rlm((d1$P50)~(d1$pgL12), method="MM"), col="cyan3")
# 
# with(d1, plot(log(K.scalar)~(alpha), xlab="log(alpha)"))
# abline(rlm(log(d1$K.scalar)~(d1$alpha), method="MM"), col="cyan3")
# 
# with(d1, plot((pgL50)~(X.tlp)))
# abline(rlm((d1$pgL50)~(d1$X.tlp), method="MM"), col="cyan3")
# 
