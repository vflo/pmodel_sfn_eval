################## optimization code ###################################

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
source("gs_stomatal_optimization_functions_phydro_calibration.R")
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
# par_data = read.csv("DATA/fitted_params_Joshi_et_al_2022.csv") 
# par_data_new = read.csv("DATA/fitted_params_flo_acclimated.csv") 
# template = read.csv("DATA/fitted_params_template.csv")
# template <- template[-c(1:160),]
path_par <- "DATA/parameters/"
par_data <- list.files(path_par) %>% 
  purrr::map_df(function(x){
    readr::read_csv(paste0(path_par,x))
  })
  
par_data$alpha <- 0.1


################################################################################

fun_no_accl = function(data, dpsi_calib=T,  k=7, 
                       stomatal_model = stomatal_model_now,
                       par_plant_acclimation = par_plant_acclimation,
                       par_cost_acclimation = par_cost_acclimation,
                       par_plant_now = par_plant,
                       par_cost_now = par_cost,
                       Species_now = species){

  data$Ciest = data$ca-data$A/data$gC
  dpsi_data = dpsi_df %>% filter(Species == Species_now)

    cat("inst resp\n")
    ndays = mean(data$Drydown.days)
    psi_min = min(data$LWP)# -6
    # psi_min = psi_crit
    psi_max = 0 #max(data$LWP)
    # cat(ndays,"\n")
    
    lwp = seq(psi_min,0, length.out=20)
    day = ndays * (lwp-psi_max)/(psi_min-psi_max)
    
    lwp_day = function(day_num){
      psi_max + day_num/ndays * (psi_min-psi_max)
    }
    
    # k = 7
    lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
    
    spl = splinefun(x = max(day):0, y=lwp_week)

    if(stomatal_model == "phydro"){
      dat_pmodel = tibble(var = 0) %>%
        mutate(p = purrr::map(var,
                              ~rphydro_analytical(tc = mean(data$T), 
                                                  ppfd = mean(data$Iabs_growth), 
                                                  vpd = mean(data$D*101325), 
                                                  co2 = mean(data$ca), elv = 0,
                                                  fapar = .99, kphio = 0.087, psi_soil = ., 
                                                  rdark = 0.02, par_plant=par_plant_acclimation, 
                                                  par_cost = par_cost_acclimation)) ) %>% 
        unnest_wider(p)
      dat1 = tibble(var = lwp, jmax_a=dat_pmodel$jmax, vcmax_a=dat_pmodel$vcmax) %>%
        mutate(var = case_when(var>0~0,
                               TRUE~var),,
               p = purrr::pmap(list(var, jmax_a, vcmax_a), 
                               ~ rphydro_instantaneous_analytical(vcmax = ..3,jmax = ..2,
                                                                  tc = mean(data$T), 
                                                                  ppfd = mean(data$Iabs_used), 
                                                                  vpd = mean(data$D*101325), 
                                                                  co2 = mean(data$ca), elv = 0, 
                                                                  fapar = .99, kphio = 0.087, 
                                                                  psi_soil = ..1, rdark = 0.02, 
                                                                  par_plant=par_plant_now, 
                                                                  par_cost = par_cost_now)) ) %>% 
        unnest_wider(p)
    }else{
      dat_pmodel = tibble(var = 0) %>%
        mutate(p = purrr::map(var,
                              ~model_numerical(tc = mean(data$T), 
                                               ppfd = mean(data$Iabs_growth), 
                                               vpd = mean(data$D*101325), 
                                               co2 = mean(data$ca), elv = 0,
                                               fapar = .99, kphio = 0.087, psi_soil = ., 
                                               rdark = 0.02, par_plant=par_plant_acclimation, 
                                               par_cost = par_cost_acclimation,
                                               stomatal_model = stomatal_model)) ) %>% 
        unnest_wider(p)
      dat1 = tibble(var = lwp, jmax_a=dat_pmodel$jmax, vcmax_a=dat_pmodel$vcmax) %>%
        mutate(var = case_when(var>0~0,
                               TRUE~var),
               p = purrr::pmap(list(var, jmax_a, vcmax_a), 
                               ~model_numerical_instantaneous(tc = mean(data$T), 
                                                              ppfd = mean(data$Iabs_used), 
                                                              vpd = mean(data$D*101325), 
                                                              co2 = mean(data$ca), elv = 0, 
                                                              fapar = .99, kphio = 0.087, 
                                                              psi_soil = ..1, rdark = 0.02, 
                                                              par_plant = par_plant_now, 
                                                              par_cost = par_cost_now, 
                                                              jmax = ..2, vcmax = ..3, 
                                                              stomatal_model = stomatal_model)) ) %>% 
        unnest_wider(p)

  }

  
    dat2 <- dat1 %>% filter(gs>=1e-40)
    gx = log(dat2$gs)
    gy = dat2$var
    fpsi = splinefun(x = gx, y=gy, method = "natural")
    gs0 = dat2$gs[which(dat2$var==0)]
    psi88S = fpsi(log(gs0*0.12))
    dpx = dat2$var
    dpy = dat2$dpsi
    f1 = splinefun(dpy~dpx)
    dp88S = f1(psi88S)
    psiL88S = psi88S-dp88S
  
  
  # A, G, CHI
  a_spl = splinefun(x = lwp, y=dat1$a)
  g_spl = splinefun(x = lwp, y=dat1$gs)
  c_spl = splinefun(x = lwp, y=dat1$chi)
  
  data_f <- data

  a_pred = a_spl(data_f$LWP)
  g_pred = g_spl(data_f$LWP)
  c_pred = c_spl(data_f$LWP)
  
  data_f <- data_f %>% 
    cbind(a_pred = a_pred,
          g_pred = g_pred,
          c_pred = c_pred,
          psi88S = psi88S,
          psiL88S = psiL88S,
          dp88S = dp88S,
          gs0 = gs0) 
  
  #DPSI
  # if (dpsi_calib){
    d_spl = splinefun(lwp, y=dat1$dpsi)
    d_pred = d_spl(dpsi_data$SWP)
  # }else{
  #   d_pred = NA
  # }
  dpsi_data_f <- dpsi_data %>% 
    dplyr::select(Species,DAY,MDWP,SWP,Dpsi) %>%
    dplyr::rename(LWP = SWP) %>% 
    cbind(d_pred = d_pred) %>% 
    dplyr::bind_cols(data_f %>% 
                       dplyr::select(Source,Drydown.days) %>%
                       dplyr::summarise(Source = unique(Source),
                                        Drydown.days = unique(Drydown.days)))%>% 
    
    dplyr::bind_cols(data_f %>% 
                       dplyr::select(8:13) %>%
                       dplyr::summarise_all(mean,na.rm = TRUE))
  
  # JOIN EVERYTHING
  res <- full_join(data_f, dpsi_data_f)
    
  return(res)
  
}


################################################################################

fun_accl = function(data, dpsi_calib=T, inst=F, k=7, 
                    stomatal_model = stomatal_model_now,
                    par_plant_now = par_plant,
                    par_cost_now = par_cost,
                    Species_now = species){

  data$Ciest = data$ca-data$A/data$gC
  dpsi_data = dpsi_df %>% filter(Species == Species_now)

    if (inst == F){
      cat("Acc resp\n")
      lwp = seq(min(data$LWP), 0, length.out=20)
      if(stomatal_model == "phydro"){
        dat1 = tibble(var = lwp) %>%
          mutate(var = case_when(var>0~0,
                                 TRUE~var),
                 p = purrr::map(var,
                                ~rphydro_analytical(tc = mean(data$T),
                                                    ppfd = mean(data$Iabs_growth),
                                                    vpd = mean(data$D*101325),
                                                    co2 = mean(data$ca), elv = 0,
                                                    fapar = .99, kphio = 0.087, psi_soil = .,
                                                    rdark = 0.02, par_plant=par_plant_now,
                                                    par_cost = par_cost_now)) ) %>%
          unnest_wider(p)
      }else{
        dat1 = tibble(var = lwp) %>%
          mutate(var = case_when(var>0~0,
                                 TRUE~var),
                 p = purrr::map(var,
                                ~model_numerical(tc = mean(data$T), 
                                                 ppfd = mean(data$Iabs_growth), 
                                                 vpd = mean(data$D*101325), 
                                                 co2 = mean(data$ca), elv = 0,
                                                 fapar = .99, kphio = 0.087, psi_soil = ., 
                                                 rdark = 0.02, par_plant=par_plant_now, 
                                                 par_cost = par_cost_now,
                                                 stomatal_model = stomatal_model))) %>% 
          unnest_wider(p)
      }
    }else{
      cat("inst resp\n")
      ndays = mean(data$Drydown.days)
      psi_min = min(data$LWP) #-6
      # psi_min = psi_crit
      psi_max = 0 #max(data$LWP)
      # cat(ndays,"\n")
      
      lwp = seq(psi_min,0, length.out=20)
      day = ndays * (lwp-psi_max)/(psi_min-psi_max)
      
      lwp_day = function(day_num){
        psi_max + day_num/ndays * (psi_min-psi_max)
      }
      
      # k = 7
      lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
      
      spl = splinefun(x = max(day):0, y=lwp_week)

      if(stomatal_model == "phydro"){
        dat_acc = tibble(var = spl(day)) %>%
          mutate(var = case_when(var>0~0,
                                 TRUE~var),
                 pmod = map(var, ~rphydro_analytical(tc = mean(data$T), ppfd = mean(data$Iabs_growth),
                                                     vpd = mean(data$D*101325), co2 = mean(data$ca),
                                                     elv = 0, fapar = .99, kphio = 0.087,
                                                     psi_soil = ., rdark = 0.02, par_plant = par_plant_now,
                                                     par_cost = par_cost_now))) %>%
          unnest_wider(pmod)
        dat1 = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>%
          mutate(p = purrr::pmap(list(var, jmax_a, vcmax_a),
                                 ~ rphydro_instantaneous_analytical(vcmax = ..3,jmax = ..2,
                                                                    tc = mean(data$T),
                                                                    ppfd = mean(data$Iabs_used),
                                                                    vpd = mean(data$D*101325),
                                                                    co2 = mean(data$ca), elv = 0,
                                                                    fapar = .99, kphio = 0.087,
                                                                    psi_soil = ..1, rdark = 0.02,
                                                                    par_plant=par_plant_now,
                                                                    par_cost = par_cost_now)) ) %>%
          unnest_wider(p)
      }else{
        dat_acc = tibble(var = spl(day)) %>% 
          mutate(var = case_when(var>0~0,
                                 TRUE~var),
                 pmod = map(var, ~model_numerical(tc = mean(data$T), ppfd = mean(data$Iabs_growth), 
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
    }

    dat2 <- dat1 %>% filter(gs>=1e-40)
    gx = log(dat2$gs)
    gy = dat2$var
    fpsi = splinefun(x = gx, y=gy, method = "natural")
    gs0 = dat2$gs[which(dat2$var==0)]
    psi88S = fpsi(log(gs0*0.12))
    dpx = dat2$var
    dpy = dat2$dpsi
    f1 = splinefun(dpy~dpx)
    dp88S = f1(psi88S)
    psiL88S = psi88S-dp88S
    
    
    # A, G, CHI
    a_spl = splinefun(x = lwp, y=dat1$a)
    g_spl = splinefun(x = lwp, y=dat1$gs)
    c_spl = splinefun(x = lwp, y=dat1$chi)
    
    data_f <- data
    
    a_pred = a_spl(data_f$LWP)
    g_pred = g_spl(data_f$LWP)
    c_pred = c_spl(data_f$LWP)
    
    data_f <- data_f %>% 
      cbind(a_pred = a_pred,
            g_pred = g_pred,
            c_pred = c_pred,
            psi88S = psi88S,
            psiL88S = psiL88S,
            dp88S = dp88S,
            gs0 = gs0) 
    
    #DPSI
    # if(dpsi_calib){
      d_spl = splinefun(lwp, y=dat1$dpsi)
      d_pred = d_spl(dpsi_data$SWP)
    # }else{
    #   d_pred = NA
    # }
    
    dpsi_data_f <- dpsi_data %>% 
      dplyr::select(Species,DAY,MDWP,SWP,Dpsi) %>%
      dplyr::rename(LWP = SWP) %>% 
      cbind(d_pred = d_pred) %>% 
      dplyr::bind_cols(data_f %>% 
                         dplyr::select(Source,Drydown.days) %>%
                         dplyr::summarise(Source = unique(Source),
                                          Drydown.days = unique(Drydown.days)))%>% 
      
      dplyr::bind_cols(data_f %>% 
                         dplyr::select(8:13) %>%
                         dplyr::summarise_all(mean,na.rm = TRUE))
    
    # JOIN EVERYTHING
    res <- full_join(data_f, dpsi_data_f)
    
    return(res)
    }

par_scheme <- list("phydro","phydro_cgain","phydro_wue", "phydro_cmax")
par_scheme_no_alpha <- list("phydro_wang_mod","phydro_sox_mod","phydro_sperry")



##### SIMULATION #####
get_simulations <- function(x){
  species = x$Species %>% unique()
  stomatal_model_now = x$scheme %>% unique()
  print(species)
  print(stomatal_model_now)
  dpsi_calib = x$dpsi %>% unique()
  inst = x$inst %>% unique()
  data1=filter(dat, Species==species)
  
  ##### SIMULATION WITH ACCLIMATION #####
  par_plant_acclimation = as.list(tibble(
    conductivity = x[1,6][[1]]*1e-16 ,
    psi50 = x[1,7][[1]], 
    b= x[1,8][[1]]
    ))
  par_cost_acclimation = as.list(data.frame(
    alpha  = x[1,9][[1]], 
    gamma = x[1,10][[1]]
    ))

  
  accl <- fun_accl(data1,
           Species_now=species,
           dpsi_calib = dpsi_calib, 
           inst = inst,
           par_plant_now = par_plant_acclimation,
           par_cost_now = par_cost_acclimation,
           stomatal_model = stomatal_model_now)
  
  accl <- accl %>% left_join(x[1,])

  ##### SIMULATION WITHOUT ACCLIMATION ####
  par_plant = list(
    conductivity= x[2,6][[1]]*1e-16,
    psi50 = x[2,7][[1]], 
    b= x[2,8][[1]]
  )
  par_cost = list(
    alpha  = x[2,9][[1]], 
    gamma = x[2,10][[1]]
  )

  no_accl <- fun_no_accl(data1,
                    dpsi_calib = dpsi_calib,  stomatal_model = stomatal_model_now,
                    par_plant_acclimation = par_plant_acclimation,
                    par_cost_acclimation = par_cost_acclimation,
                    par_plant_now = par_plant,
                    par_cost_now = par_cost,
                    Species_now = species)
  
  no_accl <- no_accl %>% left_join(x[2,])
  
  ##### RETURN #####
  
  return(bind_rows(accl,no_accl))
  
}

# par_data_extra <- par_data %>% 
#   filter(scheme %in% c("phydro_wue","phydro_cgain"
#                        )) %>% 
#   filter(Species %in% c( "Allocasuarina luehmannii","Olea europaea var. Meski",
#                  "Pseudotzuga menziesii","Eucalyptus pilularis","Eucalyptus populnea",
#                  "Glycine max","Quercus coccifera",
#                  "Quercus ilex","Quercus suber"
#            ))

df <- par_data %>%
  # filter(dpsi == TRUE) %>%
  # rbind(par_data_extra) %>%
  group_split(Species,scheme,dpsi) %>%
  purrr::map(get_simulations) %>%
  bind_rows()
# 
save(df, file = "DATA/simulations_16_09_2022.RData")



