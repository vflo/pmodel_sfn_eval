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
library(DEoptim)
# library(furrr)
# plan('multisession', workers = 6)
# options('future.global.maxsize'=2*1024*1024^2)
source("stomatal_optimization_functions_phydro_calibration.R")
source("gs_stomatal_optimization_functions_phydro_calibration.R")
source('hydraulic_functions.R')
source('photosynthetic_functions.R')
source("QUADP.R")

#Convenience functions
opt_curve_param_vc <- function(par,p88,p50,p12){
  p <- par[1]
  b <- par[2]
  if(p12>=-1){
    psi <- c(p88,p50)
    res <- (1/2)^((psi/p)^b)
    vulne_curve <- c(0.12,0.5)
    rmse <- sqrt(sum((vulne_curve-res)^2)/length(vulne_curve))
    
  }else{
    psi <- c(p88,p50,p12)
    res <- (1/2)^((psi/p)^b)
    vulne_curve <- c(0.12,0.5,0.88)
    rmse <- sqrt(sum((vulne_curve-res)^2)/length(vulne_curve))
  }
  
  if(b<1){rmse <- 1e6}
  return(rmse)
}

#LOAD DATA
dat <-  read.csv(file="DATA/drying_experiments_meta-analysis_extended.csv")
names = colnames(dat)
# [4] "Predawn.LWP..MPa."                                 
# [5] "A..umol.m.2.s.1."                                  
# [6] "gC..mol.m.2.s.1."
# [7] "vcmax_obs"
# [8] "jmax_obs"
# [9] "ca..ppm."                                          
# [10] "D..unitless...Pa...Pa.atm.."                       
# [11] "T..deg.C."                                                          
names[5:12] = c("LWP", "A", "gC", "vcmax_obs","jmax_obs","ca", "D", "T")
names[16] = "Iabs_growth"
names[15] = "Iabs_used"
colnames(dat) = names

dat <- dat %>% filter(!is.na(D),!is.na(LWP),!is.na(T),!is.na(ca),!is.na(Iabs_used))

dpsi_df <-  read.csv(file = "DATA/drying_experiments_dpsi_extended.csv")

traits <- read.csv(file="DATA/imputation_result_plus_subsps.csv")
traits <- traits %>% separate(taxon,c("genus","species", "subsp"),sep="_")

template <-  read.csv("DATA/fitted_params_template.csv")
template <- template %>% 
  left_join(traits) %>% 
  filter(!is.na(P50..MPa.)) %>%
  rowwise() %>%
  mutate(P50 = optimr::optimr(c(p=-2,b=2),
                              fn=opt_curve_param_vc,
                              p88=P88..MPa.,
                              p50=P50..MPa.,
                              p12=P12..MPa.)$par[1],
         b = optimr::optimr(c(p=-2,b=2),
                            fn=opt_curve_param_vc,
                            p88=P88..MPa.,
                            p50=P50..MPa.,
                            p12=P12..MPa.)$par[2]) %>% 
  ungroup()
  

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
    # geom_vline(xintercept = psi88S, col="grey", size=0.8)+
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
    # geom_vline(xintercept = psi88S, col="grey", size=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (analytical) p3 = p3 + geom_line(aes(x = var, y = gs), col="grey", size=1)
  
  p4 <- df_w_vol %>%
    # mutate(chi = ci/out_hydraulics_ca) %>% 
    ggplot() +
    geom_line(aes(x = var, y = chi), col="magenta", size=1)+
    geom_point(data=subdata, aes(x=LWP, y=1-A/gC/ca))+
    # geom_vline(xintercept = psi88S, col="grey", size=0.8)+
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
    # geom_vline(xintercept = psi88S, col="grey", size=0.8)+
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

error_fun_no_accl = function(x, data, data_template, plot=F, scale = 1, 
                             dpsi_calib=T,  k=7, 
                             stomatal_model = stomatal_model_now,
                             vcmax = vcmax,
                             jmax = jmax, 
                             Species_now = species,
                             K_sperry = K_sperry_no_acclimate){
  data$Ciest = data$ca-data$A/data$gC
  if(stomatal_model %in% par_scheme){
    par_plant_now = list(
      conductivity = K_sperry$K_sperry*1e-16,
      psi50 = data_template$P50 %>% unique(),
      b = data_template$b%>% unique()
    )
    par_cost_now = list(
      alpha = 0.1,
      gamma = x[1]
    )
  }else{
    par_plant_now = list(
      conductivity = x[1]*1e-16,
      psi50 = data_template$P50%>% unique(),
      b = data_template$b%>% unique()
    )
    par_cost_now = list(
      alpha  = 0.1
    )
  }
    if (Species_now != ""){
      dpsi_data = dpsi_df %>% filter(Species == Species_now) #read.csv(paste0("drying_experiments_meta/",dpsi_file,".csv"))
    }
    else{
      dpsi_data=NULL
    }
      lwp = data$LWP
      dat1 = tibble(var = lwp, jmax_a=jmax, vcmax_a=vcmax) %>% 
        cbind(data %>% select(t=T,Iabs_used, D,ca)) %>% 
        mutate(var = case_when(var>0~0,
                               TRUE~var),
               p = purrr::pmap(list(var, jmax_a, vcmax_a,t,Iabs_used,D,ca), 
                               ~model_numerical_instantaneous(tc = ..4, 
                                                              ppfd = ..5, 
                                                              vpd = ..6*101325, 
                                                              co2 = ..7, elv = 0, 
                                                              fapar = .99, kphio = 0.087, 
                                                              psi_soil = ..1, rdark = 0.02, 
                                                              par_plant=par_plant_now, 
                                                              par_cost = par_cost_now, 
                                                              jmax = ..2, vcmax = ..3, 
                                                              stomatal_model = stomatal_model)) ) %>% 
        unnest_wider(p)
    # }
  
  
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
  
  # a_spl = splinefun(x = lwp, y=dat1$a)
  # g_spl = splinefun(x = lwp, y=dat1$gs)
  # c_spl = splinefun(x = lwp, y=dat1$chi)
  
  data_f <- data #%>% filter(LWP >= psi88S) #use only values over Psi88S
  # w = rep(1, length(data_f$LWP))
  y2 = mean((dat1$gs - data_f$gC)^2,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2
  # y1 = mean((a_spl(data_f$LWP) - data_f$A)^2*w,na.rm  = TRUE)/mean(data_f$A,na.rm  = TRUE)^2
  # y2 = mean((g_spl(data_f$LWP) - data_f$gC)^2*w,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2#*4000
  
  # w = exp(-data$LWP/par_plant_now$psi50)
  # y4 = mean((c_spl(data_f$LWP) - data_f$Ciest/data_f$ca)^2*w,na.rm  = TRUE)/mean(data_f$Ciest/data_f$ca,na.rm  = TRUE)^2#*800
  # cat("c_spl:", c_spl(data_f$LWP), "\n")
  
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

  cat(x, "|", y2, " / ", y3, " / ", y, "\n")
  cat(x, "|", y, "\n")
  
  y
}


################################################################################

error_fun = function(x, data, data_template,  plot=F, inst=TRUE,
                     dpsi_calib=T, k=7, stomatal_model = stomatal_model_now, 
                     Species_now = species,K_sperry = K_sperry_acclimate){
  data$Ciest = data$ca-data$A/data$gC
  if(stomatal_model %in% par_scheme){
    par_plant_now = list(
      conductivity = K_sperry$K_sperry*1e-16,
      psi50 = data_template$P50%>% unique(),
      b = data_template$b%>% unique()
    )
    par_cost_now = list(
      alpha = 0.1,
      gamma = x[1]
    )
  }else{
    par_plant_now = list(
      conductivity = x[1]*1e-16,
      psi50 = data_template$P50 %>% unique(),
      b = data_template$b%>% unique()
    )
    par_cost_now = list(
      alpha  = 0.1
    )
  }
    if (Species_now != ""){
      dpsi_data = dpsi_df %>% filter(Species == Species_now) #read.csv(paste0("drying_experiments_meta/",dpsi_file,".csv"))
    }
    else{
      dpsi_data=NULL
    }
  ndays = mean(data$Drydown.days)
  psi_crit = data_template$P50 * (log(1000)/log(2)) ^ ( 1/data_template$b)
  if(min(data$LWP,na.rm = TRUE)<psi_crit){
    psi_min = psi_crit
  }else{
    psi_min = min(data$LWP,na.rm = TRUE) #-6
  }
  psi_max = 0 #max(data$LWP)
  # cat(ndays,"\n")
  
  lwp = seq(psi_min,0, length.out=20)
  day = ndays * (lwp-psi_max)/(psi_min-psi_max)
  actual_day = ndays * (data$LWP-psi_max)/(psi_min-psi_max)
  lwp_day = function(day_num){
    psi_max + day_num/ndays * (psi_min-psi_max)
  }
  
  # k = 7

  
    if (inst == F){
      cat("Acc resp\n")
      ndays = 20
      actual_day = ndays * (data$LWP-psi_max)/(psi_min-psi_max)
      lwp_week = rollmean(x = lwp_day(c(20:0, rep(0,k-1))), k = k, align = "right")
      spl = splinefun(x = 20:0, y=lwp_week)
      # lwp = seq(min(data$LWP), 0, length.out=20)
      dat_acc = tibble(var = spl(actual_day)) %>%
        mutate(var = case_when(var>0~0,
                               TRUE~var),
               p = purrr::map(var,
                              ~model_numerical(tc = mean(data$T,na.rm = TRUE),
                                               ppfd = mean(data$Iabs_growth,na.rm = TRUE),
                                               vpd = mean(data$D*101325,na.rm = TRUE),
                                               co2 = mean(data$ca,na.rm = TRUE), elv = 0,
                                               fapar = .99, kphio = 0.087, psi_soil = .,
                                               rdark = 0.02, par_plant=par_plant_now,
                                               par_cost = par_cost_now,
                                               stomatal_model = stomatal_model))) %>%
        unnest_wider(p)
      lwp = data$LWP
      dat1 = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>% 
        cbind(data %>% select(t=T,Iabs_used, D,ca)) %>% 
        mutate(var = case_when(var>0~0,
                               TRUE~var),
               p = purrr::pmap(list(var, jmax_a, vcmax_a,t,Iabs_used,D,ca), 
                               ~model_numerical_instantaneous(tc = ..4, 
                                                              ppfd = ..5, 
                                                              vpd = ..6*101325, 
                                                              co2 = ..7, elv = 0, 
                                                              fapar = .99, kphio = 0.087, 
                                                              psi_soil = ..1, rdark = 0.02, 
                                                              par_plant=par_plant_now, 
                                                              par_cost = par_cost_now, 
                                                              jmax = ..2, vcmax = ..3, 
                                                              stomatal_model = stomatal_model)) ) %>% 
        unnest_wider(p)
    }else{
      cat("inst resp\n")
      lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
      spl = splinefun(x = max(day):0, y=lwp_week)
      dat_acc = tibble(var = spl(actual_day)) %>% 
        mutate(var = case_when(var>0~0,
                               TRUE~var),
               pmod = map(var, ~model_numerical(tc = mean(data$T,na.rm = TRUE), ppfd = mean(data$Iabs_growth,na.rm = TRUE), 
                                                vpd = mean(data$D*101325,na.rm = TRUE), co2 = mean(data$ca,na.rm = TRUE), 
                                                elv = 0, fapar = .99, kphio = 0.087, 
                                                psi_soil = ., rdark = 0.02, par_plant=par_plant_now, 
                                                par_cost = par_cost_now, stomatal_model = stomatal_model))) %>% 
        unnest_wider(pmod)
      lwp = data$LWP
      dat1 = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>% 
        cbind(data %>% select(t=T,Iabs_used, D,ca)) %>% 
        mutate(var = case_when(var>0~0,
                               TRUE~var),
               p = purrr::pmap(list(var, jmax_a, vcmax_a,t,Iabs_used,D,ca), 
                               ~model_numerical_instantaneous(tc = ..4, 
                                                              ppfd = ..5, 
                                                              vpd = ..6*101325, 
                                                              co2 = ..7, elv = 0, 
                                                              fapar = .99, kphio = 0.087, 
                                                              psi_soil = ..1, rdark = 0.02, 
                                                              par_plant=par_plant_now, 
                                                              par_cost = par_cost_now, 
                                                              jmax = ..2, vcmax = ..3, 
                                                              stomatal_model = stomatal_model)) ) %>% 
        unnest_wider(p)
      }
    # }
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
      
    # a_spl = splinefun(x = lwp, y=dat1$a)
    # g_spl = splinefun(x = lwp, y=dat1$gs)
    # c_spl = splinefun(x = lwp, y=dat1$chi)
      
    data_f <- data #%>% filter(LWP >= psi88S) #use only values over Psi88S
    # w = rep(1, length(data_f$LWP))
    y2 = mean((dat1$gs - data_f$gC)^2,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2
    # y1 = mean((a_spl(data_f$LWP) - data_f$A)^2*w,na.rm  = TRUE)/mean(data_f$A,na.rm  = TRUE)^2
    # y2 = mean((g_spl(data_f$LWP) - data_f$gC)^2*w,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2#*4000
      
    # w = exp(-data$LWP/par_plant_now$psi50)
    # y4 = mean((c_spl(data_f$LWP) - data_f$Ciest/data_f$ca)^2*w,na.rm  = TRUE)/mean(data_f$Ciest/data_f$ca,na.rm  = TRUE)^2#*800
    # cat("c_spl:", c_spl(data_f$LWP), "\n")
      
    if (dpsi_calib ){
      d_spl = splinefun(lwp, y=dat1$dpsi)
      # w = 1 #0.1 + exp(dpsi_data$SWP)
      dpsi_data_f <- dpsi_data #%>% filter(SWP >= psi88S) #use only values over Psi88S
      y3 = mean((d_spl(dpsi_data_f$SWP) - dpsi_data_f$Dpsi)^2,na.rm  = TRUE)/mean(dpsi_data_f$Dpsi,na.rm  = TRUE)^2 #*40
      # y3 = mean((d_spl(dpsi_swp$SWP) - dpsi_swp$Dpsi)^2*w)/mean(dpsi_swp$Dpsi)^2 #*40
      cat("d_spl:", d_spl(dpsi_data_f$SWP), "\n")
    }else{
      y3=0
    }
    # y = (y1+y2+y3+y4)
      
    y=y2+y3
      
      cat(x, "|", y2, " / ", y3, " / ", y, "\n")
      cat(x, "|", y, "\n")
      
      y
}


par_scheme <- list("phydro","phydro_cgain","phydro_wue", "phydro_cmax")
par_scheme_no_alpha <- list("phydro_wang_mod","phydro_sox_mod","phydro_sperry")

##### PARAMETERIZATION #####
get_parameters <- function(x){
    species = x$Species %>% unique()
    stomatal_model_now = x$scheme %>% unique()
    dpsi_calib_now = x$dpsi %>% unique()
    inst = x$inst
    data_template_now = x
    data1 = filter(dat, Species==species, Source == unique(x$source))
    if(!is.null(K_sperry)){
    K_sperry_no_acclimate = K_sperry %>% 
      filter(Species == species,
             acclimation == FALSE,
             dpsi == dpsi_calib_now,
             source == x$source) %>% 
      select(K_sperry) %>% 
      unique()
    K_sperry_acclimate = K_sperry %>% 
      filter(Species == species,
             acclimation == TRUE,
             dpsi == dpsi_calib_now,
             source == x$source) %>% 
      select(K_sperry) %>% 
      unique()
    }
    ##### PARAMETERIZATION WITHOUT ACCLIMATION #####

    data_ww <- data1 %>% 
      filter(!is.na(gC)) %>% 
      mutate(LWP_q90 = quantile(LWP, 0.9,na.rm = TRUE),
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
    # 
    # jmax <- vcmax * 1.67 #Medlyn et al. 2002

                                        
    print(stomatal_model_now)
    print(species)
    parameter_max <- 10
    if(stomatal_model_now %in% c("phydro_cgain")){
      parameter_max <- 50}
    if(stomatal_model_now %in% c("phydro_cmax")){
      parameter_max <- 6}

      optimise(error_fun_no_accl,
               interval = c(0,parameter_max),
               data=data1,
               data_template = data_template_now,
               dpsi_calib = dpsi_calib_now,
               stomatal_model = stomatal_model_now,
               vcmax = vcmax,
               jmax = jmax, 
               Species = species,
               K_sperry = K_sperry_no_acclimate
      ) -> opt_no_accl
      x_no_accl <- opt_no_accl$minimum
    
    error_fun_no_accl(x_no_accl, data1, 
                      data_template = data_template_now, plot=T, 
                      dpsi_calib = dpsi_calib_now,
                      stomatal_model = stomatal_model_now, ,
                      vcmax = vcmax,
                      jmax = jmax,
                      Species_now = species,
                      K_sperry = K_sperry_no_acclimate)
    
    if(stomatal_model_now %in% par_scheme){
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=K_sperry_no_acclimate$K_sperry,
                            alpha=NA,
                            gamma=x_no_accl[1])
    }else if(stomatal_model_now %in% par_scheme_no_alpha){
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=x_no_accl[1],
                            alpha=NA,
                            gamma=NA)
    }else{
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=x_no_accl[1],
                            alpha=0.1,
                            gamma=NA)
    }
    
    
    ##### PARAMETERIZATION WITH ACCLIMATION #####
    print(stomatal_model_now)
    print(species)

    parameter_max <- 10
    if(stomatal_model_now %in% c("phydro")){
      if(species == "Broussonetia papyrifera"){
        parameter_max <- 1000
      }else{
        parameter_max <- 10
      }
    }
    if(stomatal_model_now %in% c("phydro_cmax")){
      if(species %in% c("Broussonetia papyrifera","Cinnamomum bodinieri")){
        parameter_max <- 1 #B papyrifera = 10
      }else{
        parameter_max <- 6
      }
    }
    if(stomatal_model_now %in% c("phydro_sox_mod")){
      if(species %in% c("Broussonetia papyrifera","Platycarya longipes")){
        parameter_max <- 30000
      }else{
        parameter_max <- 6
      }
    }
    if(stomatal_model_now %in% c("phydro_cgain")){
      if(species != "Broussonetia papyrifera"){
        parameter_max <- 50
      }else{
        parameter_max <- 30
        }
      }
    # if(stomatal_model_now %in% c("phydro_cmax")){
    #   parameter_max <- 6}
    
      optimise(error_fun,
               interval = c(0,parameter_max),
               data=data1,
               data_template = data_template_now,
               dpsi_calib = dpsi_calib_now,
               inst = inst,
               stomatal_model = stomatal_model_now, 
               Species_now = species,
               K_sperry = K_sperry_acclimate
      ) -> opt_accl
      x_accl <- opt_accl$minimum

    error_fun(x_accl, data1, data_template = data_template_now,
              plot=T, dpsi_calib = dpsi_calib_now, inst = inst,
              stomatal_model = stomatal_model_now, Species_now = species,
              K_sperry = K_sperry_acclimate)
    
    if(stomatal_model_now %in% par_scheme){
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=K_sperry_acclimate$K_sperry,
                         alpha=0.1,
                         gamma=x_accl[1])
    }else{
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=x_accl[1],
                         alpha=0.1,
                         gamma=NA)
    }
    

  df <- bind_rows(res_accl,res_no_accl)
  
  readr::write_csv(df,file=paste0("DATA/parameter_kmax_vpd/",stomatal_model_now,"_",species,"_",dpsi_calib_now,"_",x$source,".csv"))
  
  return(bind_rows(res_accl,res_no_accl))

}

##### COMPUTE PARAMETERS #####
#First compute sperry model to obtain Kmax for CMAX. CGAIN, WUE and PHYDRO models
K_sperry <- NULL
template %>% filter(scheme == "phydro_sperry",dpsi == FALSE
                    # Species %in% c("Rosa cymosa",
                    #                "Broussonetia papyrifera",
                    #                "Cinnamomum bodinieri",
                    #                "Platycarya longipes",
                    #                "Pteroceltis tatarinowii")
                    ) %>%
  group_split(scheme, dpsi, Species,source) %>%
  purrr::map_df(get_parameters)->res

save(res,file = "DATA/K_sperry_meta-analysis_kmax_vpd.RData")

load(file = "DATA/K_sperry_meta-analysis_kmax_vpd.RData")

K_sperry <- res %>% 
  select(Species,K_sperry = K.scale,dpsi,acclimation,source) %>% 
  group_by(Species,dpsi, acclimation,source) %>% 
  summarise_all(unique)

#Compute the rest of the models
template %>% 
  filter(!scheme %in% c("phydro_sperry"),
         dpsi == FALSE
         # Species %in% c(#"Rosa cymosa",
         #                "Broussonetia papyrifera",
         #                 #"Cinnamomum bodinieri"
         #                 "Platycarya longipes"
         #                 # "Pteroceltis tatarinowii"
         #                )
         ) %>%
  # filter(scheme %in% c("phydro")) %>%
  # filter(scheme %in% c("phydro_cmax"), Species == "Quercus ilex") %>%
  group_split(scheme, dpsi, Species,source) %>% 
  purrr::map_df(get_parameters)

View(res)
