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
# library(furrr)
# plan('multisession', workers = 6)
# options('future.global.maxsize'=2*1024*1024^2)
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
par_data = read.csv("DATA/fitted_params_Joshi_et_al_2022.csv") 
# par_data_new = read.csv("DATA/fitted_params_flo_acclimated.csv") 
template = read.csv("DATA/fitted_params_template.csv")
template <- template%>%
filter(Species %in% c("Allocasuarina luehmannii","Olea europaea var. Meski",
                    "Pseudotzuga menziesii","Eucalyptus pilularis","Eucalyptus populnea",
                    "Glycine max","Quercus coccifera","Quercus ilex","Quercus suber"
                     )#,
       # !scheme %in% c("phydro_wang", "phydro_sox")
       )
# template <- template[96,]

plot_all = function(df_w_vol, varname, species, data, dpsi_data=NULL, analytical=F){
  # df_w_vol = df_w_vol[complete.cases(df_w_vol),]
  # View(df_w_vol)
  
  gx = log(df_w_vol$gs+1e-20)
  gy = df_w_vol$var
  # gx1 = seq(max(min(gx), -20), max(gx), length.out=100)
  f = splinefun(x = gx, y=gy, method = "monoH.FC")
  gs0 = df_w_vol$gs[which(df_w_vol$var==0)]
  psi88S = f(log(gs0*0.12))
  psi50S = f(log(gs0*0.50))
  psi12S = f(log(gs0*0.88))
  cat("psi88S = ", psi88S, "\n")
  cat("psi50S = ", psi50S, "\n")
  cat("psi12S = ", psi50S, "\n")
  
  dpx = df_w_vol$var
  dpy = df_w_vol$dpsi
  f1 = splinefun(dpy~dpx, method = "monoH.FC")
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
    xlab(varname)+
    theme_bw()
  if (analytical) p1 = p1 + geom_line(aes(x = var, y = vcmax ), col="grey", size=1) 
  
  p2 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = dpsi), col="blue", size=1)+
    geom_vline(xintercept = psi88S, col="grey", size=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (!is.null(dpsi_data)){
    cat("Adding points...\n")
    p2 <- p2 + 
      geom_point(data=dpsi_data, aes(x=SWP, y=Dpsi))
  }
  
  p3 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = gs), col="cyan2", size=1)+
    geom_point(data=subdata, aes(x=LWP, y=gC))+
    geom_vline(xintercept = psi88S, col="grey", size=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (analytical) p3 = p3 + geom_line(aes(x = var, y = gs), col="grey", size=1)
  
  p4 <- df_w_vol %>%
    # mutate(chi = ci/out_hydraulics_ca) %>% 
    ggplot() +
    geom_line(aes(x = var, y = chi), col="magenta", size=1)+
    geom_point(data=subdata, aes(x=LWP, y=1-A/gC/ca))+
    geom_vline(xintercept = psi88S, col="grey", size=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (analytical) p4 = p4 + geom_line(aes(x = var, y = chi), col="grey", size=1)
  
  p5 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = jmax), col="goldenrod1", size=1) +
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Jmax))+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  
  
  p6 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = a), col="green4", size=1) +
    geom_point(data=subdata, aes(x=LWP, y=A))+
    geom_vline(xintercept = psi88S, col="grey", size=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (analytical) p6 = p6 + geom_line(aes(x = var, y = gpp), col="grey",size=1) 
  
  grid.arrange(p3,p6,p1,p2,p5,p4, ncol=2)
}


pmodel_wang17 = function(tc, ppfd, vpd, co2, elv, fapar, kphio, ...){
  
  out_analytical <- rpmodel::rpmodel(
    tc             = tc,
    vpd            = vpd,
    co2            = co2,
    elv            = elv,
    kphio          = kphio,
    beta           = 146,
    fapar          = fapar,
    ppfd           = ppfd*1e-6*86400, # p-model requires in mol/m2/day
    ...
  )
  
  # Convert some outputs to facilitate comparison
  return(list_modify(out_analytical,
                     gs = out_analytical$gs/86400*rpmodel::calc_patm(elv), # mol/m2/day/Pa --> mol/m2/s
                     gpp = out_analytical$gpp/1.03772448,  # gC/m2/day --> umol/m2/s
                     vcmax = out_analytical$vcmax/0.0864,   # mol/m2/day --> umol/m2/s
                     jmax = out_analytical$jmax/0.0864   # mol/m2/day --> umol/m2/s
  )
  )
}


################################################################################

error_fun_no_accl = function(x, data, plot=F, scale = 1, 
                             dpsi_calib=T,  k=7, stomatal_model = stomatal_model_now,
                             vcmax = vcmax,
                             jmax = jmax, 
                             Species_now = species){
  x = x*scale
  cat("start: ",x, "\n")
  
  data$Ciest = data$ca-data$A/data$gC
  
  if(stomatal_model %in% par_scheme){
    if(Species_now == "Helianthus annuus"){
      par_plant_now = list(
        conductivity = x[1]*1e-16,
        psi50 = x[2],
        b = 1.4
      )
    }else{
      par_plant_now = list(
        conductivity = x[1]*1e-16,
        psi50 = x[2],
        b = 1
      )
    }
    par_cost_now = list(
      alpha = 0.1,
      gamma = x[3]
    )
    if(x[1]<=0| x[2]>=0 |  x[3] < 0|
       x[1]>1e4 | x[2]<(-15)
    ){over <- TRUE}else{over <- FALSE} #set boundaries
  }else{
    par_plant_now = list(
      conductivity= x[1]*1e-16,
      psi50 = x[2],
      b= x[3]
    )
    par_cost_now = list(
      alpha  = 0.1
    )
    if(x[1]<=0| x[2]>=0 | x[3]<1 | 
       x[1]>1e4 | x[2]<(-15)
    ){over <- TRUE}else{over <- FALSE} #set boundaries
  }
  
  if(over){
    print(paste("Over-limits"))
    return(1e6)
  }else{
    if (Species_now != ""){
      dpsi_data = dpsi_df %>% filter(Species == Species_now) #read.csv(paste0("drying_experiments_meta/",dpsi_file,".csv"))
    }
    else{
      dpsi_data=NULL
    }
    
    
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
    
    # dat_pmodel = pmodel_wang17(tc = mean(data$T),
    #                            ppfd = mean(data$Iabs_used),
    #                            vpd = mean(data$D*101325),
    #                            co2 = mean(data$ca), elv = 0,
    #                            fapar = .99, kphio = 0.087)
    if(stomatal_model == "phydro"){
      dat1 = tibble(var = lwp, jmax_a=jmax, vcmax_a=vcmax) %>%
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
      dat1 = tibble(var = lwp, jmax_a=jmax, vcmax_a=vcmax) %>%
        mutate(var = case_when(var>0~0,
                               TRUE~var),,
               p = purrr::pmap(list(var, jmax_a, vcmax_a), 
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
  if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=Species_now, data = data, dpsi_data=dpsi_data)
  
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
  
  a_spl = splinefun(x = lwp, y=dat1$a)
  g_spl = splinefun(x = lwp, y=dat1$gs)
  c_spl = splinefun(x = lwp, y=dat1$chi)
  
  data_f <- data #%>% filter(LWP >= psi88S) #use only values over Psi88S
  w = rep(1, length(data_f$LWP))
  y1 = mean((a_spl(data_f$LWP) - data_f$A)^2*w,na.rm  = TRUE)/mean(data_f$A,na.rm  = TRUE)^2
  y2 = mean((g_spl(data_f$LWP) - data_f$gC)^2*w,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2#*4000
  
  # w = exp(-data$LWP/par_plant_now$psi50)
  y4 = mean((c_spl(data_f$LWP) - data_f$Ciest/data_f$ca)^2*w,na.rm  = TRUE)/mean(data_f$Ciest/data_f$ca,na.rm  = TRUE)^2#*800
  cat("c_spl:", c_spl(data_f$LWP), "\n")
  
  if (dpsi_calib ){
    d_spl = splinefun(lwp, y=dat1$dpsi)
    # w = 1 #0.1 + exp(dpsi_data$SWP)
    dpsi_data_f <- dpsi_data #%>% filter(SWP >= psi88S) #use only values over Psi88S
    y3 = mean((d_spl(dpsi_data_f$SWP) - dpsi_data_f$Dpsi)^2,na.rm  = TRUE)/mean(dpsi_data_f$Dpsi,na.rm  = TRUE)^2 #*40
    # y3 = mean((d_spl(dpsi_swp$SWP) - dpsi_swp$Dpsi)^2*w)/mean(dpsi_swp$Dpsi)^2 #*40
    cat("d_spl:", d_spl(dpsi_data_f$SWP), "\n")
  }else{
    # w=1
    # p50 = -1.5 #par_plant_now$psi50
    # y3 = mean((dat1$dpsi[lwp>p50] - -p50)^2*w)/mean(dat1$dpsi[lwp>p50])^2 #*20
    y3=0
  }
  # y = (y1+y2+y3+y4)
  
  y=y2+y3
  # if(stomatal_model %in% c("phydro", "phydro_cgain")){
  #   y = y2 + y3 + y4
  # }else{y = y2+y4}

  cat(x, "|", y1, "/ ", y2, " / ", y3, " / ", y4, " / ", y, "\n")
  cat(x, "|", y, "\n")
  
  
  x0 <<- x
  y
  }
}


################################################################################

error_fun = function(x, data,  plot=F, scale = 1, 
                     dpsi_calib=T, inst=F, k=7, stomatal_model = stomatal_model_now, 
                     Species_now = species, par_plant = par_plant_no_acclimated,
                     par_cost = par_cost_no_acclimation){
  x = x*scale
  #cat("start: ",x, "\n")
  
  data$Ciest = data$ca-data$A/data$gC
  
  # if(stomatal_model %in% par_scheme){
  #   if(Species_now == "Helianthus annuus"){
  #     par_plant_now = list(
  #       conductivity= par_plant$conductivity,
  #       psi50 = par_plant$psi50,
  #       b= 1.4
  #     )
  #   }else{
  #     par_plant_now = list(
  #       conductivity= par_plant$conductivity,
  #       psi50 = par_plant$psi50,
  #       b= par_plant$b
  #     )
  #   }
  # 
  #   par_cost_now = list(
  #     alpha  = x[1], 
  #     gamma = par_cost$gamma# x[2]
  #   )
  #   
  #   if(x[1] <= 0 | x[1] > 0.4 #| x[2] < 0 
  #      ){over <- TRUE}else{over <- FALSE} #set boundaries
  # }else{
  #   par_plant_now = list(
  #     conductivity= par_plant$conductivity,
  #     psi50 = par_plant$psi50,
  #     b= par_plant$b
  #   )
  #   par_cost_now = list(
  #     alpha  = x[1]
  #   )
  #   if(x[1] <= 0 | x[1] > 0.4){over <- TRUE}else{over <- FALSE} #set boundaries
  # }
  if(stomatal_model %in% par_scheme){
    if(Species_now == "Helianthus annuus"){
      par_plant_now = list(
        conductivity = x[1]*1e-16,
        psi50 = x[2],
        b = 1.4
      )
    }else{
      par_plant_now = list(
        conductivity = x[1]*1e-16,
        psi50 = x[2],
        b = 1
      )
    }
    
    par_cost_now = list(
      alpha  = 0.1, 
      gamma = x[3]
    )
    
    if(x[1]<=0| x[2]>=0 | x[3] < 0|
       x[1]>1e4 | x[2]<(-15)#| x[2] < 0 
    ){over <- TRUE}else{over <- FALSE} #set boundaries
  }else{
    par_plant_now = list(
      conductivity = x[1]*1e-16,
      psi50 = x[2],
      b = x[3]
    )
    par_cost_now = list(
      alpha  =0.1
    )
    if(x[1]<=0| x[2]>=0 | x[3]<1 | 
       x[1]>1e4 | x[2]<(-15) ){over <- TRUE}else{over <- FALSE} #set boundaries
  }
  if(over){
    print(paste("Over-limits"))
    return(1e6)
  }else{
    if (Species_now != ""){
      dpsi_data = dpsi_df %>% filter(Species == Species_now) #read.csv(paste0("drying_experiments_meta/",dpsi_file,".csv"))
    }
    else{
      dpsi_data=NULL
    }
    
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
      
      # plot(x=day, y=lwp_day(day))
      # points(x=day, y=spl(day), col="red", type="o")
      # points(x=max(day):0, y=lwp_week, col="blue", type="l")
      if(stomatal_model == "phydro"){
        dat_acc = tibble(var = spl(day)) %>%
          mutate(var = case_when(var>0~0,
                                 TRUE~var),
                 pmod = map(var, ~rphydro_analytical(tc = mean(data$T), ppfd = mean(data$Iabs_growth),
                                                  vpd = mean(data$D*101325), co2 = mean(data$ca),
                                                  elv = 0, fapar = .99, kphio = 0.087,
                                                  psi_soil = ., rdark = 0.02, par_plant=par_plant_now,
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
    # if(plot==T) dat_acc %>% plot_all(varname = "psi_soil", species=species, data = data, dpsi_data=dpsi_data, analytical = F)
    if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=species, data = data, dpsi_data=dpsi_data)
    
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
    
    a_spl = splinefun(x = lwp, y=dat1$a)
    g_spl = splinefun(x = lwp, y=dat1$gs)
    c_spl = splinefun(x = lwp, y=dat1$chi)
    
    
    data_f <- data #%>% filter(LWP >= psi88S) #use only values over Psi88S
    w = rep(1, length(data_f$LWP))
    y1 = mean((a_spl(data_f$LWP) - data_f$A)^2*w,na.rm  = TRUE)/mean(data_f$A,na.rm  = TRUE)^2
    y2 = mean((g_spl(data_f$LWP) - data_f$gC)^2*w,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2#*4000
    
    # w = exp(-data$LWP/par_plant_now$psi50)
    y4 = mean((c_spl(data_f$LWP) - data_f$Ciest/data_f$ca)^2*w,na.rm  = TRUE)/mean(data_f$Ciest/data_f$ca,na.rm  = TRUE)^2#*800
    cat("c_spl:", c_spl(data_f$LWP), "\n")
    
    if (dpsi_calib == T){
      d_spl = splinefun(lwp, y=dat1$dpsi)
      # w = 1 #0.1 + exp(dpsi_data$SWP)
      dpsi_data_f <- dpsi_data #%>% filter(SWP >= psi88S) #use only values over Psi88S
      y3 = mean((d_spl(dpsi_data_f$SWP) - dpsi_data_f$Dpsi)^2,na.rm  = TRUE)/mean(dpsi_data_f$Dpsi,na.rm  = TRUE)^2 #*40
      # y3 = mean((d_spl(dpsi_swp$SWP) - dpsi_swp$Dpsi)^2*w)/mean(dpsi_swp$Dpsi)^2 #*40
      cat("d_spl:", d_spl(dpsi_data_f$SWP), "\n")
    }else{
      # w=1
      # p50 = -1.5 #par_plant_now$psi50
      # y3 = mean((dat1$dpsi[lwp>p50] - -p50)^2*w)/mean(dat1$dpsi[lwp>p50])^2 #*20
      y3=0
    }
    
    y=y2+y3
    # y = (y1+y2+y3+y4)
  #     if(stomatal_model  %in% c("phydro", "phydro_cgain")){
  #   y = y2 + y3 + y4
  # }else{y = y2+y4}
    cat(x, "|", y1, "/ ", y2, " / ", y3, " / ", y4, " / ", y, "\n")
    cat(x, "|", y, "\n")
    
    
    x0 <<- x
    y
  }
}




# x0 <- c(3.060048, -2.166255, 6.514309)
# species <- "Allocasuarina luehmannii"
# stomatal_model_now = "phydro_cgain"
# x0 <- c(0.541, -0.789,   NA,      0.120)
# list( "Allocasuarina luehmannii","Olea europaea var. Meski",
#       "Pseudotzuga menziesii","Eucalyptus pilularis","Eucalyptus populnea",
#       "Glycine max","Quercus coccifera",
#       "Quercus ilex","Quercus suber"
# ) %>%
# list("Broussonetia papyrifera (Linnaeus) L_He ritier ex Ventenat","Ficus tikoua",
#      "Cinnamomum bodinieri H. Leveille", "Platycarya longipes Wu",
#      "Pteroceltis tatarinowii Maximowicz","Rosa cymosa Trattinnick"
# ) %>%
# list("Helianthus annuus") %>%
# list("Olea europaea var. Chemlali","Cedrus atlantica") %>%
# list("Broussonetia papyrifera (Linnaeus) L_He ritier ex Ventenat","Ficus tikoua",
#       "Cinnamomum bodinieri H. Leveille", "Platycarya longipes Wu",
#       "Pteroceltis tatarinowii Maximowicz","Rosa cymosa Trattinnick"
#      )

par_scheme <- list("phydro","phydro_cgain","phydro_wue", "phydro_cmax")
par_scheme_no_alpha <- list("phydro_wang_mod","phydro_sox_mod","phydro_sperry")

##### PARAMETERIZATION #####
get_parameters <- function(x){
    species = x$Species
    stomatal_model_now = x$scheme
    dpsi_calib = x$dpsi
    inst = x$inst
    data1=filter(dat, Species==species)
    
    ##### PARAMETERIZATION WITHOUT ACCLIMATION #####

    data_ww <- data1 %>% 
      mutate(LWP_q90 = quantile(LWP, 0.9),
             ci = ca-A/gC) %>% 
      filter(LWP>=LWP_q90) %>% 
      dplyr::select(LWP,A,gC,T,ci,Iabs_growth) %>% 
      dplyr::summarise_all(mean, na.rm = TRUE)
      
      
    vcmax <- calc_vcmax_no_acclimated_ww(A = data_ww$A,
                                         ci = data_ww$ci,
                                         tc = data_ww$T,
                                         patm = calc_patm(0,data_ww$T),
                                         rdark = 0.02
    )
    jmax <- calc_jmax_no_acclimated_ww(A = data_ww$A,
                                       vcmax = vcmax,
                                       ci = data_ww$ci,
                                       I = data_ww$Iabs_growth,
                                       tc = data_ww$T,
                                       patm = calc_patm(0,data_ww$T),
                                       kphio = 0.087
                                       )
    

                                        
    if(stomatal_model_now %in% par_scheme){x0 <- c(1,-1,1)}else{x0 <- c(1,-1,1)}

    
    conv <- TRUE
    count <- 0
    while(all(conv , count < 2)){
      print(stomatal_model_now)
      print(species)
      optimr::optimr(par = x0/abs(x0),
                     fn = error_fun_no_accl,
                     data=data1,
                     scale=abs(x0),
                     dpsi_calib = dpsi_calib,       # Set to F if dpsi data is not available
                     stomatal_model = stomatal_model_now,
                     vcmax = vcmax,
                     jmax = jmax, 
                     Species = species,
                     control = list(par_scale = 1, maxit=2000, all.methods = TRUE)
      ) -> opt
      convergence_no_accl <- opt$convergence
      count <- count + 1
      if(convergence_no_accl==0){conv = FALSE}else{conv = TRUE}
    }
    
    # error_fun_no_accl(x0, data1,plot=T,dpsi_calib = dpsi_calib,
    #                   stomatal_model = stomatal_model_now, ,
    #                   vcmax = vcmax,
    #                   jmax = jmax,
    #                   Species_now = species)
    # x_no_accl <- x0
    x_no_accl <- abs(opt$par) * x0
    if(stomatal_model_now %in% par_scheme){
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=x_no_accl[1],
                            psi50=x_no_accl[2],
                            b = 1,
                            alpha=NA,
                            gamma=x_no_accl[3],
                            convergence = convergence_no_accl)
    }else if(stomatal_model_now %in% par_scheme_no_alpha){
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=x_no_accl[1],
                            psi50=x_no_accl[2],
                            b = x_no_accl[3],
                            alpha=NA,
                            gamma=NA,
                            convergence = convergence_no_accl)
    }else{
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=x_no_accl[1],
                            psi50=x_no_accl[2],
                            b = x_no_accl[3],
                            alpha=0.1,
                            gamma=NA,
                            convergence = convergence_no_accl)
    }
    

##### PARAMETERIZATION WITH ACCLIMATION #####

    
    
    if(stomatal_model_now %in% par_scheme){
      par_plant_no_acclimation = list(
        conductivity= x_no_accl[1]*1e-16,
        psi50 = x_no_accl[2], 
        b= 1
      )
      par_cost_no_acclimation = list(
        gamma =  x_no_accl[3]
      )
    }else{
      par_plant_no_acclimation = list(
        conductivity=  x_no_accl[1]*1e-16,
        psi50 =  x_no_accl[2], 
        b=  x_no_accl[3]
      )
      par_cost_no_acclimation = list(
        gamma  = NA
      )
    }
    
    
    if(stomatal_model_now %in% par_scheme){x0 <- c(1,-1,1)}else{x0 <- c(1,-1,1)}

# error_fun(x0, data1,# %>% mutate(Drydown.days=38),
#           plot=T, dpsi_file=species,
#           dpsi_calib = FALSE, inst=T, stomatal_model = stomatal_model_now)

    conv <- TRUE
    count <- 0
    while(all(conv , count < 2)){
      print(stomatal_model_now)
      print(species)
      optimr::optimr(par = x0/abs(x0),
                     fn = error_fun,
                     data=data1,# %>% mutate(Drydown.days=38),
                     scale=abs(x0),
                     dpsi_calib = dpsi_calib,       # Set to F if dpsi data is not available
                     inst=inst,
                     stomatal_model = stomatal_model_now, 
                     Species_now = species,
                     par_plant = par_plant_no_acclimation,
                     par_cost = par_cost_no_acclimation,
            # method = "Brent",
            # lower = 0,
            # upper = 0.4
                     control = list(par_scale = 1, maxit=2000, all.methods = TRUE)
                     ) -> opt_accl
      convergence_accl <- opt_accl$convergence
      count <- count + 1
      if(convergence_accl==0){conv = FALSE}else{conv = TRUE}
      }

# error_fun(x0, data1,# %>% mutate(Drydown.days=38),
#           plot=T, dpsi_calib = dpsi_calib,inst=T,
#           stomatal_model = stomatal_model_now, Species_now = species,
#           par_plant = par_plant_no_acclimation,
#           par_cost = par_cost_no_acclimation)

    x_accl <- abs(opt_accl$par) * x0
    if(stomatal_model_now %in% par_scheme){
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=x_accl[1],
                         psi50=x_accl[2],
                         b = 1,#x_accl[3],
                         alpha=0.1,#x_accl[4],
                         gamma=x_accl[3],
                         convergence = convergence_accl)
    }else{
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=x_accl[1],
                         psi50=x_accl[2],
                         b = x_accl[3],
                         alpha=0.1,#x_accl[4],
                         gamma=NA,
                         convergence = convergence_accl)
      }


  df <- bind_rows(res_accl,res_no_accl)
  
  readr::write_csv(df,file=paste0("DATA/parameters/",stomatal_model_now,"_",species,"_",dpsi_calib,"_",inst,".csv"))
  
  return(bind_rows(res_accl,res_no_accl))

}

##### COMPUTE PARAMETERS #####
template[96,] %>% 
  # filter(Species != "Helianthus annuus",
  #        !scheme %in% c("phydro_wang", "phydro_sox"),
  #        dpsi == TRUE) %>% 
  split(seq(nrow(.))) %>%  
  purrr::map_df(get_parameters)->res

View(res)
