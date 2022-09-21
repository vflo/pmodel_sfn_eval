################## optimization code ###################################

rm(list=ls())

library(tidyverse)
library(purrr)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)
library(stringr)
library(rphydro)
library(rstatix)
# library(furrr)
# plan('multisession', workers = 6)
# options('future.global.maxsize'=2*1024*1024^2)
source("stomatal_optimization_functions_phydro_calibration.R")
source("gs_stomatal_optimization_functions_phydro_calibration.R")
source('hydraulic_functions.R')
source('photosynthetic_functions.R')
source("QUADP.R")
#Convenience functions
opt_curve_param <- function(par,psi50,d,c){
  # p <- par[1]
  b <- par[1]
  psi <- seq(0,-9,-0.01)
  res <- (1/2)^((psi/psi50)^b)
  vulne_curve <- exp(-(psi/d)^c)
  rmse <- sqrt(sum((vulne_curve-res)^2)/length(vulne_curve))
  return(rmse)
}

opt_curve_param_2 <- function(par,psi50,d,c){
  p <- par[1]
  b <- par[2]
  psi <- seq(0,-9,-0.01)
  res <- (1/2)^((psi/p)^b)
  vulne_curve <- exp(-(psi/d)^c)
  rmse <- sqrt(sum((vulne_curve-res)^2)/length(vulne_curve))
  return(rmse)
}


#Import meteo data from Caldes de Montbui Automatic Meteorological Station https://www.meteo.cat/observacions/xema/dades?codi=X9&dia=2021-09-15T00:00Z
dat <-  read_csv("DATA/FR_Pue_HH_EC+Sap.csv")
#Summarise central meteo data at the daily scale
meteo <- meteo_caldes_2019 %>% 
  dplyr::select(-hour) %>% 
  mutate(DATE = lubridate::dmy(date),
         VPD = bigleaf::rH.to.VPD(HRM/100,TM)*1000,
         D = VPD/(PM*100),
         PPFD = bigleaf::Rg.to.PPFD(RS),
         ca = 400) %>% 
  group_by(date) %>% 
  summarise_all(mean,na.rm = TRUE) %>% 
  dplyr::select(DATE,
                "T" = TM,
                D,
                ca,
                "Iabs_used" = PPFD
                )

#Calculate the average meteo values for the long and short acclimation period
long_period_meteo <- meteo %>% 
  dplyr::select(-DATE) %>% 
  summarise_all(mean,na.rm = TRUE) %>% 
  cbind(Drydown.days = 34)

medium_period_meteo <- meteo %>%
  filter(DATE<=lubridate::dmy("18/07/2019")) %>% 
  dplyr::select(-DATE) %>% 
  summarise_all(mean,na.rm = TRUE) %>% 
  cbind(Drydown.days = 30)

short_period_meteo <- meteo %>%
  filter(DATE<=lubridate::dmy("15/07/2019")) %>% 
  dplyr::select(-DATE) %>% 
  summarise_all(mean,na.rm = TRUE) %>% 
  cbind(Drydown.days = 27)

all_period_meteo <- rbind(long_period_meteo,medium_period_meteo,short_period_meteo)

#Import DRESS dataset
Dataset_TM <- read_csv("DATA/Dataset_Modelos_Victor.csv")

#Transform variables and select the species
foo <- Dataset_TM %>% mutate(DATE = lubridate::dmy(DATE_CALEND)) %>% 
  filter(TREAT_2019 == "D", !is.na(GS_LEAF_TIME), CAMPAIGN %in% c("i", "d"#, "r", "dd"
                                                                  )
  ) %>% 
  arrange(SPECIES, DATE) %>% 
  select(DATE,DATE_CONT, gs_leaf_time_hour, CAMPAIGN, RH, TEMP, RAD, VPD, SPECIES,
         PTLP_Caract:RGR) %>% 
  mutate(gC = GS*1e-3/1.6,
         LWP = -WP_PD/10,
         K = Kl_Caract *18 * 1e3,
         Dpsi = MD_PD.WP/10,
         Dpsi = case_when(Dpsi<0 ~ NA_real_,
                          TRUE ~ Dpsi),
         PPFD = bigleaf::Rg.to.PPFD(RAD), 
         Species = plyr::mapvalues(SPECIES, 
                                   from = c("VT","SN","AG","PA","SC","BP","AC","QP","PS","QI","FA",
                                            "AM","PH","PL","AU","MC","FL","CM","QC","BS"),
                                   to = c("Viburnum spp.","Sambucus nigra","Alnus glutinosa",
                                          "Populus spp.","Salix cinerea","Betula pubescens",
                                          "Acer campestre","Quercus petraea","Pinus sylvestris",
                                          "Quercus ilex","Fraxinus spp.","Acer monspessulanum",
                                          "Pinus halepensis","Pistacia lentiscus","Arbutus unedo",
                                          "Myrtus communis","Phillyrea latifolia","Crataegus monogyna",
                                          "Quercus coccifera","Buxus sempervirens"))) 

foo %>%
  group_by(Species) %>%
  pairwise_t_test(LWP~CAMPAIGN,alternative = "greater")

foo <- foo %>% 
  filter(Species %in% c('Viburnum spp.', #genus level
                        # 'Fraxinus spp.', #genus level
                        'Myrtus communis',
                        'Quercus ilex',
                        'Pinus sylvestris',
                        'Pinus halepensis',
                        'Phillyrea latifolia',
                        'Acer campestre',
                        'Arbutus unedo'
                        # 'Populus spp.' #genus level #no long enough
                        )) %>%  # select species with enough separation in soil water potential between irrigation and drought determined by T-test
  group_by(Species) %>% 
  mutate(Drydown.days = max(DATE_CONT,na.rm = TRUE))

foo %>% ggplot(aes(LWP,gC, color = CAMPAIGN))+
  geom_point()+
  facet_wrap(~SPECIES)

foo %>% ggplot(aes(LWP,Dpsi, color = CAMPAIGN))+
  geom_point()+
  facet_wrap(~SPECIES)

#Import species parameters data
medpar <- medfate::SpParamsMED %>% 
  dplyr::select(Name,Kmax_stemxylem,VCstem_c,VCstem_d,Vmax298, Jmax298) %>% 
  dplyr::rename(Species= Name)

dat <- left_join(foo,medpar)



b_param <- dat %>% 
  group_by(Species) %>% 
  summarise(P50 = unique(P50),
            VCstem_d = unique(VCstem_d),
            VCstem_c = unique(VCstem_c)) %>% 
  split(seq(nrow(.))) %>% 
  purrr::map_df(function(x){
    
    b = optimize(opt_curve_param, c(1,30), psi50 = -x$P50, d=x$VCstem_d,c=x$VCstem_c)
    alt = optim(c(-2,3), opt_curve_param_2, psi50 = -x$P50, d=x$VCstem_d,c=x$VCstem_c)
    return(tibble(Species=x$Species,b = b$minimum,P50_alt = alt$par[1], b_alt = alt$par[2]))
})

dat <- dat %>% left_join(b_param) %>% left_join(all_period_meteo)

dat <- full_join(tibble(dpsi_calib = c(TRUE, FALSE)),dat,by = character()) 
dat <- full_join(tibble(scheme = c("phydro","phydro_wang","phydro_wang_mod",
                                   "phydro_sox","phydro_sox_mod","phydro_sperry",
                                   "phydro_wue", "phydro_cgain", "phydro_cmax")),
                 dat,
                 by = character())

#####################################################################################
# Plot function
plot_all = function(df_w_vol, varname, species, data, analytical=F){
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
  
  subdata = data %>% filter(LWP > -6.5)
  
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
  if (any(!is.na(data$Dpsi))){
    cat("Adding points...\n")
    p2 <- p2 + 
      geom_point(data=data, aes(x=LWP, y=Dpsi))
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
    ggplot() +
    geom_line(aes(x = var, y = jmax), col="goldenrod1", size=1) +
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Jmax))+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()

  
  grid.arrange(p3,p2,p1,p4, ncol=2)
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
                             dpsi_calib=T,  k=7, 
                             stomatal_model = stomatal_model_now,
                             vcmax = vcmax,
                             jmax = jmax, 
                             Species_now = species){
  x = x*scale
  
  if(stomatal_model %in% par_scheme){
      par_plant_now = list(
        conductivity = x[1]*1e-16,
        psi50 = -data$P50 %>% unique(),
        b = data$b_alt%>% unique()
      )
    par_cost_now = list(
      alpha = 0.1,
      gamma = x[2]
    )
    if(x[1]<=0 | x[1]>1e4 | x[2] < 0
    ){over <- TRUE}else{over <- FALSE} #set boundaries
  }else{
    par_plant_now = list(
      conductivity = x[1]*1e-16,
      psi50 = -data$P50%>% unique(),
      b = data$b_alt%>% unique()
    )
    par_cost_now = list(
      alpha  = 0.1
    )
    if(x[1]<=0| x[1]>1e4 
    ){over <- TRUE}else{over <- FALSE} #set boundaries
  }
  
  if(over){
    print(paste("Over-limits"))
    return(1e6)
  }else{
    p_crit = par_plant_now$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant_now$b)
    ndays = mean(data$Drydown.days)
    # psi_min = min(data$LWP)# -6
    psi_min = p_crit
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
                               ~model_numerical_instantaneous_gs(tc = mean(data$T), 
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
  
  if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=Species_now, data = data)
  
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
  
    g_spl = splinefun(x = lwp, y=dat1$gs)
    
    data_f <- data %>%  filter(CAMPAIGN == "i")#%>% filter(LWP >= psi88S) #use only values over Psi88S
    w = rep(1, length(data_f$LWP))
    y1 = mean((g_spl(data_f$LWP) - data_f$gC)^2*w,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2
    
    if (dpsi_calib == T){
      d_spl = splinefun(lwp, y=dat1$dpsi)
      y2 = mean((d_spl(data_f$LWP) - data_f$Dpsi)^2,na.rm  = TRUE)/mean(data_f$Dpsi,na.rm  = TRUE)^2 #*40
      cat("d_spl:", d_spl(data_f$LWP), "\n")
    }else{
      y2=0
    }
    
    y=y1+y2
    
    cat(x, "|", y1, "/ ", y2, " / ", y, "\n")
    cat(x, "|", y, "\n")
  
  
  x0 <<- x
  y
  }
}


################################################################################

error_fun = function(x, data,  plot=F, scale = 1, 
                     dpsi_calib=T, k=7, stomatal_model = stomatal_model_now, 
                     Species_now = species){
  x = x*scale
  
  if(stomatal_model %in% par_scheme){
    par_plant_now = list(
      conductivity = x[1]*1e-16,
      psi50 = -data$P50%>% unique(),
      b = data$b_alt%>% unique()
    )
    par_cost_now = list(
      alpha = 0.1,
      gamma = x[2]
    )
    if(x[1]<=0 | x[1]>1e4 | x[2] < 0
    ){over <- TRUE}else{over <- FALSE} #set boundaries
  }else{
    par_plant_now = list(
      conductivity = x[1]*1e-16,
      psi50 = -data$P50 %>% unique(),
      b = data$b_alt%>% unique()
    )
    par_cost_now = list(
      alpha  = 0.1
    )
    if(x[1]<=0| x[1]>1e4 
    ){over <- TRUE}else{over <- FALSE} #set boundaries
  }
  if(over){
    print(paste("Over-limits"))
    return(1e6)
  }else{
      p_crit = par_plant_now$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant_now$b)
      ndays = mean(data$Drydown.days)
      # psi_min = min(data$LWP)# -6
      psi_min = p_crit
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
                 pmod = map(var, ~rphydro_analytical(tc = mean(data$T), ppfd = mean(data$Iabs_used),
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
               pmod = map(var, ~model_numerical_gs(tc = mean(data$T), ppfd = mean(data$Iabs_used), 
                                                vpd = mean(data$D*101325), co2 = mean(data$ca), 
                                                elv = 0, fapar = .99, kphio = 0.087, 
                                                psi_soil = ., rdark = 0.02, par_plant=par_plant_now, 
                                                par_cost = par_cost_now, stomatal_model = stomatal_model))) %>% 
        unnest_wider(pmod)
      dat1 = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>%
        mutate(p = purrr::pmap(list(var, jmax_a, vcmax_a), 
                               ~model_numerical_instantaneous_gs(tc = mean(data$T), 
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
      
    if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=species, data = data)
    
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

    g_spl = splinefun(x = lwp, y=dat1$gs)
    
    data_f <- data %>%  filter(CAMPAIGN == "i") #%>% filter(LWP >= psi88S) #use only values over Psi88S
    w = rep(1, length(data_f$LWP))
    y1 = mean((g_spl(data_f$LWP) - data_f$gC)^2*w,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2

    if (dpsi_calib == T){
      d_spl = splinefun(lwp, y=dat1$dpsi)
      y2 = mean((d_spl(data_f$LWP) - data_f$Dpsi)^2,na.rm  = TRUE)/mean(data_f$Dpsi,na.rm  = TRUE)^2 #*40
      cat("d_spl:", d_spl(data_f$LWP), "\n")
    }else{
      y2=0
    }
    
    y=y1+y2

    cat(x, "|", y1, "/ ", y2, " / ", y, "\n")
    cat(x, "|", y, "\n")
    
    
    x0 <<- x
    y
  }
}


par_scheme <- list("phydro","phydro_cgain","phydro_wue", "phydro_cmax")
par_scheme_no_alpha <- list("phydro_wang_mod","phydro_sox_mod","phydro_sperry")

##### PARAMETERIZATION #####
get_parameters <- function(x){
    species = x$Species %>% unique()
    stomatal_model_now = x$scheme %>% unique()
    dpsi_calib = x$dpsi_calib %>% unique()
    tc = x$T %>% mean(na.rm = TRUE)
    PPFD = x$PPFD %>% mean(na.rm = TRUE)
    VPD = x$VPD %>% mean(na.rm = TRUE)
    vc25 = x$Vmax298 %>% mean(na.rm = TRUE)
    j25 = x$Jmax298 %>% mean(na.rm = TRUE)
    data1 = x
    
    ##### PARAMETERIZATION WITHOUT ACCLIMATION #####

    #VCmax and Jmax parameters are not determinant on this analysis since we are
    # only interested in gs.
    vcmax <- rpmodel::calc_ftemp_inst_vcmax(tc)*vc25
    jmax <- vcmax*1.67 #Medlyn et al. 2002
    if(is.na(jmax)){
      foo <- rpmodel::rpmodel(tc = tc,vpd = VPD,co2 = 400,fapar = 0.99,ppfd = PPFD,patm = 101325,elv = 203)
      jmax = foo$jmax
      vcmax = foo$vcmax} 

                                        
    if(stomatal_model_now %in% par_scheme){x0 <- c(1,1)}else{x0 <- c(1)}
    # error_fun_no_accl(x0, data1, #%>% mutate(Drydown.days=38), 
    #                   plot=T, dpsi_file=species,
    #                   dpsi_calib = dpsi_calib,  stomatal_model = stomatal_model_now,par_plant_no_acclimation = par_plant_no_acclimation,
    #                   par_cost_no_acclimation = par_cost_no_acclimation)
    
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

    # x_no_accl <- x0
    x_no_accl <- abs(opt$par) * x0
    
    error_fun_no_accl(x_no_accl, data1, plot=T, dpsi_calib = dpsi_calib,
                      stomatal_model = stomatal_model_now, ,
                      vcmax = vcmax,
                      jmax = jmax,
                      Species_now = species)
    
    if(stomatal_model_now %in% par_scheme){
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=x_no_accl[1],
                            psi50= -data1$P50 %>% unique(),
                            b_used = data1$b_alt %>% unique(),
                            alpha=NA,
                            gamma=x_no_accl[2],
                            convergence = convergence_no_accl)
    }else if(stomatal_model_now %in% par_scheme_no_alpha){
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=x_no_accl[1],
                            psi50= -data1$P50 %>% unique(),
                            b_used = data1$b_alt %>% unique(),
                            alpha=NA,
                            gamma=NA,
                            convergence = convergence_no_accl)
    }else{
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            psi50= -data1$P50 %>% unique(),
                            b_used = data1$b_alt %>% unique(),
                            alpha=0.1,
                            gamma=NA,
                            convergence = convergence_no_accl)
    }


##### PARAMETERIZATION WITH ACCLIMATION #####
    
    if(stomatal_model_now %in% par_scheme){x0 <- c(1,1)}else{x0 <- c(1)}

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
                     data=data1 ,# %>% mutate(Drydown.days=38),
                     scale=abs(x0),
                     dpsi_calib = dpsi_calib,       # Set to F if dpsi data is not available
                     stomatal_model = stomatal_model_now, 
                     Species_now = species,
                     control = list(par_scale = 1, maxit=2000, all.methods = TRUE)
                     ) -> opt_accl
      convergence_accl <- opt_accl$convergence
      count <- count + 1
      if(convergence_accl==0){conv = FALSE}else{conv = TRUE}
      }


    x_accl <- abs(opt_accl$par) * x0
    
    error_fun(x_accl, data1,# %>% mutate(Drydown.days=38),
              plot=T, dpsi_calib = dpsi_calib,
              stomatal_model = stomatal_model_now, Species_now = species)
    
    if(stomatal_model_now %in% par_scheme){
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=x_accl[1],
                         psi50= -data1$P50 %>% unique(),
                         b_used = data1$b_alt %>% unique(),
                         alpha=0.1,#x_accl[4],
                         gamma=x_accl[2],
                         convergence = convergence_accl)
    }else{
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=x_accl[1],
                         psi50= -data1$P50 %>% unique(),
                         b_used = data1$b_alt %>% unique(),
                         alpha=0.1,#x_accl[4],
                         gamma=NA,
                         convergence = convergence_accl)
      }


  df <- bind_rows(res_accl,res_no_accl)
  
  readr::write_csv(df,file=paste0("DATA/parameters_torre_marimon/",stomatal_model_now,"_",species,"_",dpsi_calib,".csv"))
  
  return(bind_rows(res_accl,res_no_accl))

}

##### COMPUTE PARAMETERS #####
dat %>% filter(scheme != "phydro") %>% group_split(scheme, dpsi_calib, Species) %>% 
  purrr::map_df(get_parameters)->res

View(res)
