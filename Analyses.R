library(sapfluxnetr)
library(rpmodel)
library(tidyverse)
library(lmerTest)
library(lubridate)
library(Metrics)
library(ghibli)
library(xts)
library(wCorr)
library(pROC)
library(jcolors)
library(emmeans)

ghibli_pal <- ghibli_palette("MononokeMedium", 7, direction = -1, type = c("discrete"))
#### ANALISIS ####


# LOAD METADATA -----------------------------------------------------------

load(file = "DATA/sfn_meta.RData")


# LOAD WEEKLY ------------------------------------------------------------
path_weekly_folder <- "DATA/OUTCOME_WEEKLY"
paths_weekly <- list.files(path = path_weekly_folder,full.names = T)

paths_weekly %>% 
  purrr::map_df(function(x){
    load(file=paste0(x))
    df_res
  })->df

df <- df %>% 
  left_join(sfn_meta$site_md %>% dplyr::select(-si_elev), by = "si_code") %>% 
  mutate(sensitivity_K = case_when(is.na(sensitivity_K)~1,
                                   TRUE~sensitivity_K),
         sensitivity_psi = case_when(is.na(sensitivity_psi)~1,
                                   TRUE~sensitivity_psi),
         gpp = case_when(model_type %in% c("pmodel")~gpp*1e6,
                         TRUE~gpp),
         a = case_when(model_type %in% c("phydro","sperry", "wang", "wap")~a*LAI,
                         TRUE~a),
         gpp_all = case_when(is.na(gpp)~a*3600*1e-6, #to mol m-2 h-1
                             TRUE ~ gpp*3600*1e-6),  #to mol m-2 h-1
         log_psi = log(-psi_soil)) %>% 
  filter(!is.na(E_sapflow),!is.na(E), E>0, E_sapflow>0
         ,low_swp==FALSE, model_type != "pmodel_swc"
         ) %>% 
  group_by(si_code, model_type,sensitivity_K,sensitivity_psi) %>%
  mutate(max_E = max(E),
         max_E_sapflow = max(E_sapflow),
         E_rel = E/max_E,
         E_sapflow_rel = E_sapflow/max_E_sapflow,
         E_dif = E - E_sapflow,
         E_dif_rel = E_rel - E_sapflow_rel) %>% 
  as_tibble()



# WEEKLY DATA ------------------------------------------------------
df_lm <- df %>% 
  group_by(si_code, model_type,sensitivity_K,sensitivity_psi) %>% 
  do(lm(E_sapflow ~ E, data = .) %>% broom::glance())


df_summary <- df %>% 
  group_by(si_code, model_type,sensitivity_K,sensitivity_psi) %>% 
  summarise(n = n(),
            cor_value = weightedCorr(E_sapflow,E, method ="Pearson"),
            rmse = rmse(E_sapflow, E),
            mae = mae(E_sapflow, E),
            bias = bias(E_sapflow, E),
            si_map = unique(si_map),
            si_mat = unique(si_mat),
            si_elev = unique(si_elev),
            si_biome = unique(si_biome),
            si_igbp = unique(si_igbp),
            sand = unique(st_sand_perc),
            clay = unique(st_clay_perc),
            FAPAR = mean(FAPAR,na.rm = TRUE),
            VPD = mean(vpd,na.rm = TRUE),
            Ta = mean(ta,na.rm = TRUE),
            swvl = mean(swvl,na.rm = TRUE),
            log_psi = mean(log_psi,na.rm = TRUE),
            ppfd = mean(ppfd_in,na.rm = TRUE),
            LAI = mean(LAI,na.rm = TRUE),
            gpp_all = mean(gpp_all, na.rm = TRUE)) 

df_summary <- df_summary %>% left_join(df_lm)

sensitivity_references <- df_summary %>% 
  ungroup() %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  mutate(ref_bias = bias,
         ref_cor = cor_value,
         ref_r2 = r.squared,
         ref_mae = mae,
         ref_rmse = rmse
         ) %>% 
  dplyr::select(si_code,model_type,ref_bias,ref_cor,ref_r2,ref_mae,ref_rmse)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

df %>%
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(log_E = log(E),
         log_E_sapflow = log(E_sapflow),
         dens = get_density(E,E_sapflow,n = 100)) %>% 
  ggplot(aes(E, E_sapflow, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type, scales = "free_x")+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab("E actual [mol m-2 h-1]")+
  xlab("E modeled [mol m-2 h-1]")+
  theme_bw()

df %>%
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(log_E = log(E),
         log_E_sapflow = log(E_sapflow),
         dens = get_density(log_E,log_E_sapflow,n = 100)) %>% 
  ggplot(aes(log(E), log(E_sapflow), group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type, scales = "free_x")+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab(expression(paste("ln(E actual) [ln(mol ", m^{-2}, h^{-1},")]")))+
  xlab(expression(paste("ln(E modeled) [ln(mol ", m^{-2}, h^{-1},")]")))+
  theme_bw()

df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(E_rel,E_sapflow_rel,n = 200)) %>% 
  ggplot(aes(E_rel*100, E_sapflow_rel*100, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab("Relative E actual [%]")+
  xlab("Relative E modeled [%]")+
  theme_bw()


mod_E_E <- lmer(log(E_sapflow)~log(E)*model_type+(log(E)|si_code), data = df%>% 
                  filter(sensitivity_K == 1, sensitivity_psi == 1))
summary(mod_E_E)

emt1 <- emtrends(mod_E_E, "model_type", var = "E")
emt1         
pairs(emt1)


mod_E_E_rel <- lmer(E_sapflow_rel~E_rel*model_type+(1|si_code), data = df%>% 
                      filter(sensitivity_K == 1, sensitivity_psi == 1))
summary(mod_E_E_rel)

emt1_rel <- emtrends(mod_E_E_rel, "model_type", var = "E_rel")
emt1_rel   
pairs(emt1_rel)

#### SWC all data ####
df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(swvl,E_dif,n = 200)) %>% 
  ggplot(aes(swvl, E_dif, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 2)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("SWC [v/v]")+
  theme_bw()

df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(swvl,E_dif_rel,n = 200)) %>% 
  ggplot(aes(swvl, E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab("SWC [v/v]")+
  theme_bw()


#### SOIL PSI all data ####
df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  filter(!is.na(log_psi)) %>% 
  mutate(dens = get_density(log_psi,E_dif,n = 200)) %>% 
  ggplot(aes(log_psi, E_dif, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 2)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab(expression(paste("log(|",psi[soil],"|)")))+
  theme_bw()

df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>%
  filter(!is.na(log_psi)) %>%
  mutate(dens = get_density(log_psi,E_dif_rel,n = 200)) %>% 
  ggplot(aes(log_psi, E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab(expression(paste("log(|",psi[soil],"|)")))+
  theme_bw()


#### VPD all data ####
df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(vpd,E_dif,n = 200)) %>% 
  ggplot(aes(log(vpd), E_dif, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("log(VPD) [log(kPa)]")+
  theme_bw()

df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(vpd,E_dif_rel,n = 200)) %>% 
  ggplot(aes(log(vpd), E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab("log(vpd) [log(kPa)]")+
  theme_bw()


#### PPFD all data ####
df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1, ppfd_in<1100) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(ppfd_in,E_dif,n = 200)) %>% 
  ggplot(aes(ppfd_in, E_dif, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab(expression(paste("PPFD [", mu,"mol ", m^{-2},s^{-1},"]")))+
  theme_bw()

df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1, ppfd_in<1100) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(ppfd_in,E_dif_rel,n = 200)) %>% 
  ggplot(aes(ppfd_in, E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab(expression(paste("PPFD [", mu,"mol ", m^{-2},s^{-1},"]")))+
  theme_bw()


#### Temperature all data ####
df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(ta,E_dif,n = 200)) %>% 
  ggplot(aes(ta, E_dif, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("T [ºC]")+
  theme_bw()

df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(ta,E_dif_rel,n = 200)) %>% 
  ggplot(aes(ta, E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab("T [ºC]")+
  theme_bw()

#### LAI all data ####
df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(LAI,E_dif,n = 200)) %>% 
  ggplot(aes(LAI, E_dif, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("LAI")+
  theme_bw()

df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(LAI,E_dif_rel,n = 200)) %>% 
  ggplot(aes(LAI, E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab("LAI")+
  theme_bw()


#### MODELS weekly ####

model_pmodel <- lmer(log(E_sapflow) ~ log(E)*(swvl+log(vpd)+ppfd_in+ta+LAI)+(1|si_code), 
                     data = df %>% filter(model_type == "pmodel") %>% mutate(ppfd_in = ppfd_in/1000 ))
model_pmodel <- lmer(E_dif ~ swvl+log(vpd)+ppfd_in+ta+LAI +(1|si_code), 
                     data = df %>% filter(model_type == "pmodel") %>% mutate(ppfd_in = ppfd_in/1000 ))
summary(model_pmodel)
MuMIn::r.squaredGLMM(model_pmodel)


model_pmodel_swc <- lmer(log(E_sapflow) ~ log(E)*(swvl+log(vpd)+ppfd_in+ta+LAI)+(1|si_code), 
                     data = df %>% filter(model_type == "pmodel_swc") %>% mutate(ppfd_in = ppfd_in/1000 ))
model_pmodel_swc <- lmer(E_dif ~ swvl+log(vpd)+ppfd_in+ta + (1|si_code), 
                         data = df %>% filter(model_type == "pmodel_swc") %>% mutate(ppfd_in = ppfd_in/1000 ))
summary(model_pmodel_swc)
MuMIn::r.squaredGLMM(model_pmodel_swc)


model_phydro <- lmer(log(E_sapflow) ~ log(E)*(swvl+log(vpd)+ppfd_in+ta+LAI)+(1|si_code), 
                     data = df %>% 
                       filter(model_type == "phydro",sensitivity_K == 1, sensitivity_psi == 1) %>%
                       mutate(ppfd_in = ppfd_in/1000 ))
model_phydro <- lmer(E_dif ~ swvl+log(vpd)+ppfd_in+ta+LAI + (1|si_code),
                     data = df %>% 
                       filter(model_type == "phydro",sensitivity_K == 1, sensitivity_psi == 1)%>% 
                       mutate(ppfd_in = ppfd_in/1000 ))
summary(model_phydro)
MuMIn::r.squaredGLMM(model_phydro)


model_sperry <- lmer(log(E_sapflow) ~ log(E)*(swvl+log(vpd)+ppfd_in+ta+LAI)+(1|si_code), 
                     data = df %>% 
                       filter(model_type == "sperry",sensitivity_K == 1, sensitivity_psi == 1) %>% 
                       mutate(ppfd_in = ppfd_in/1000 ))
model_sperry <- lmer(E_dif ~ swvl+log(vpd)+ppfd_in+ta+LAI + (1|si_code),
                     data = df %>% 
                       filter(model_type == "sperry",sensitivity_K == 1, sensitivity_psi == 1)%>%
                       mutate(ppfd_in = ppfd_in/1000 ))
summary(model_sperry)
MuMIn::r.squaredGLMM(model_sperry)


model_wang <- lmer(log(E_sapflow) ~ log(E)*(swvl+log(vpd)+ppfd_in+ta+LAI)+(1|si_code), 
                     data = df %>% 
                     filter(model_type == "wang",sensitivity_K == 1, sensitivity_psi == 1) %>% 
                     mutate(ppfd_in = ppfd_in/1000 ))
model_wang <- lmer(E_dif ~ swvl+log(vpd)+ppfd_in+ta+LAI + (1|si_code),
                     data = df %>% 
                     filter(model_type == "wang",sensitivity_K == 1, sensitivity_psi == 1)%>% 
                     mutate(ppfd_in = ppfd_in/1000 ))
summary(model_wang)
MuMIn::r.squaredGLMM(model_wang)


model_wap <- lmer(log(E_sapflow) ~ log(E)*(swvl+log(vpd)+ppfd_in+ta)+(1|si_code), 
                     data = df %>% 
                    filter(model_type == "wap",sensitivity_K == 1, sensitivity_psi == 1) %>% 
                    mutate(ppfd_in = ppfd_in/1000 ))
model_wap <- lmer(E_dif ~ swvl+log(vpd)+ppfd_in+ta+LAI + (1|si_code),
                   data = df %>% 
                    filter(model_type == "wap",sensitivity_K == 1, sensitivity_psi == 1)%>% 
                    mutate(ppfd_in = ppfd_in/1000 ))
summary(model_wap)
MuMIn::r.squaredGLMM(model_wap)


#### MAP ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()



#### MAT ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()


#### ELEVATION ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Site elevation [m]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Site elevation [m]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Site elevation [m]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Site elevation [m]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Site elevation [m]")+
  theme_bw()

#### VPD AVERAGE ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=VPD, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Average VPD [kPa]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=VPD, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Average VPD [kPa]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=VPD, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Average VPD [kPa]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=VPD, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Average VPD [kPa]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=VPD, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Average VPD [kPa]")+
  theme_bw()


#### T AVERAGE ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=Ta, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Average T [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=Ta, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Average T [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=Ta, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Average T [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=Ta, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Average T [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=Ta, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Average T [ºC]")+
  theme_bw()


#### PPFD AVERAGE ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=ppfd, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Average PPFD")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=ppfd, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Average PPFD")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=ppfd, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Average PPFD")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=ppfd, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Average PPFD")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=ppfd, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Average PPFD")+
  theme_bw()


#### SWC AVERAGE ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=swvl, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Average SWC")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=swvl, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Average SWC")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=swvl, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Average SWC")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=swvl, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Average SWC")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=swvl, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Average SWC")+
  theme_bw()



#### SOIL PSI AVERAGE ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=log_psi, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab(expression(paste("log(|",psi[soil],"|)")))+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=log_psi, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab(expression(paste("log(|",psi[soil],"|)")))+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=log_psi, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab(expression(paste("log(|",psi[soil],"|)")))+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=log_psi, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab(expression(paste("log(|",psi[soil],"|)")))+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=log_psi, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab(expression(paste("log(|",psi[soil],"|)")))+
  theme_bw()


#### SAND PERCENTAGE ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=sand, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Soil sand percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=sand, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Soil sand percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=sand, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Soil sand percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=sand, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Soil sand percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=sand, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Soil sand percentage")+
  theme_bw()


#### CLAY PERCENTAGE ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=clay, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Soil clay percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=clay, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Soil clay percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=clay, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Soil clay percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=clay, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Soil clay percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=clay, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Soil clay percentage")+
  theme_bw()


#### LAI PERCENTAGE ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=LAI, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("LAI")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=LAI, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("LAI")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=LAI, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("LAI")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=LAI, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("LAI")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=LAI, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("LAI")+
  theme_bw()



#### GPP MODELED ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>%
  filter(gpp_all>=0) %>% 
  ggplot(aes(x=gpp_all, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Average GPP modeled")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  filter(gpp_all>=0) %>% 
  ggplot(aes(x=gpp_all, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Average GPP modeled")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>%
  filter(gpp_all>=0) %>% 
  ggplot(aes(x=gpp_all, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Average GPP modeled")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  filter(gpp_all>=0) %>% 
  ggplot(aes(x=gpp_all, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Average GPP modeled")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  filter(gpp_all>=0) %>% 
  ggplot(aes(x=gpp_all, y = mae, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Average GPP modeled")+
  theme_bw()

#### BIOMES ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_biome, y = cor_value, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Biomes")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_biome, y = r.squared, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Biomes")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_biome, y = bias, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("BIAS")+
  xlab("Biomes")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_biome, y = rmse, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Biomes")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_biome, y = mae, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Biomes")+
  theme_bw()


#### IGBP ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_igbp, y = cor_value, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("IGBP land cover classification")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_igbp, y = r.squared, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("IGBP land cover classification")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_igbp, y = bias, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("BIAS")+
  xlab("IGBP land cover classification")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_igbp, y = rmse, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("IGBP land cover classification")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_igbp, y = mae, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("IGBP land cover classification")+
  theme_bw()



#### METHODS ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=model_type, y = cor_value, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Model")+
  theme_bw()


df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=model_type, y = r.squared, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Model")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=model_type, y = bias, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Model")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=model_type, y = rmse, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Model")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=model_type, y = mae, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("MAE")+
  xlab("Model")+
  theme_bw()



#### MODELS SUMMARY ####
bias <- lmer(bias~model_type +(1|si_code), data = df_summary%>% 
               filter(sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(bias)
emmeans::emmeans(bias,'model_type') %>% test()
emmeans::emmeans(bias,'model_type') %>% pairs()

r2 <- lmer(r.squared~model_type +(1|si_code), data = df_summary%>% 
               filter(sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(r2)
emmeans::emmeans(r2,'model_type')
emmeans::emmeans(r2,'model_type') %>% pairs()

corr <- lmer(cor_value~model_type +(1|si_code), data = df_summary%>% 
             filter(sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(corr)
emmeans::emmeans(corr,'model_type')
emmeans::emmeans(corr,'model_type') %>% pairs()

rmse <- lmer(rmse~model_type +(1|si_code), data = df_summary%>% 
               filter(sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(rmse)
emmeans::emmeans(rmse,'model_type')
emmeans::emmeans(rmse,'model_type') %>% pairs()

mae <- lmer(mae~model_type +(1|si_code), data = df_summary%>% 
               filter(sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(mae)
emmeans::emmeans(mae,'model_type')
emmeans::emmeans(mae,'model_type') %>% pairs()


#PMODEL
model_pmodel <- lm(bias ~ si_map + si_mat + si_elev + clay + sand + VPD + Ta + swvl + ppfd + si_igbp + si_biome , 
                     data = df_summary %>% filter(model_type == "pmodel"), weights = n)
summary(model_pmodel)
step(model_pmodel) %>% summary()

model_pmodel <- lm(cor_value ~ si_map + si_mat + si_elev + clay + sand + VPD + Ta + swvl + ppfd + si_igbp + si_biome, 
                   data = df_summary %>% filter(model_type == "pmodel"), weights = n)
summary(model_pmodel)
step(model_pmodel) %>% summary()

model_pmodel <- lm(r.squared ~ si_map + si_mat + si_elev + clay + sand + VPD + Ta + swvl + ppfd + si_igbp + si_biome, 
                   data = df_summary %>% filter(model_type == "pmodel"), weights = n)
summary(model_pmodel)
step(model_pmodel) %>% summary()

model_pmodel <- lm(rmse ~ si_map + si_mat + si_elev + clay + sand + VPD + Ta + swvl + ppfd + si_igbp + si_biome, 
                   data = df_summary %>% filter(model_type == "pmodel"), weights = n)
summary(model_pmodel)
step(model_pmodel) %>% summary()


#### SENSITIVITY ANALYSIS ####
df_summary %>% left_join(sensitivity_references) %>% 
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_psi ==1) %>% 
  group_by(model_type) %>% 
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_K, bias-ref_bias, group = model_type))+
  geom_line(aes(sensitivity_K, bias-ref_bias, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab("Bias - Reference Bias")+
  xlab("K Sensitivity multiplier")+
  theme_bw()

df_summary %>% left_join(sensitivity_references) %>%
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_K ==1) %>% 
  group_by(model_type) %>%
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_psi, bias-ref_bias, group = model_type))+
  geom_line(aes(sensitivity_psi, bias-ref_bias, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab("Bias - Reference Bias")+
  xlab("P50 Sensitivity multiplier")+
  theme_bw()



df_summary %>% left_join(sensitivity_references) %>% 
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_psi ==1) %>% 
  group_by(model_type) %>% 
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_K, rmse-ref_rmse, group = model_type))+
  geom_line(aes(sensitivity_K, rmse-ref_rmse, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab("RMSE - Reference RMSE")+
  xlab("K Sensitivity multiplier")+
  theme_bw()

df_summary %>% left_join(sensitivity_references) %>%
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_K ==1) %>% 
  group_by(model_type) %>%
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_psi, rmse-ref_rmse, group = model_type))+
  geom_line(aes(sensitivity_psi, rmse-ref_rmse, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab("RMSE - Reference RMSE")+
  xlab("P50 Sensitivity multiplier")+
  theme_bw()


df_summary %>% left_join(sensitivity_references) %>% 
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_psi ==1) %>% 
  group_by(model_type) %>% 
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_K, r.squared - ref_r2, group = model_type))+
  geom_line(aes(sensitivity_K, r.squared - ref_r2, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab(expression(paste(R^{2}, " - Reference ", R^{2})))+
  xlab("K Sensitivity multiplier")+
  theme_bw()

df_summary %>% left_join(sensitivity_references) %>%
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_K ==1) %>% 
  group_by(model_type) %>%
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_psi, r.squared - ref_r2, group = model_type))+
  geom_line(aes(sensitivity_psi, r.squared - ref_r2, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab(expression(paste(R^{2}, " - Reference ", R^{2})))+
  xlab("P50 Sensitivity multiplier")+
  theme_bw()


df_summary %>% left_join(sensitivity_references) %>% 
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_psi ==1) %>% 
  group_by(model_type) %>% 
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_K, cor_value-ref_cor, group = model_type))+
  geom_line(aes(sensitivity_K,cor_value-ref_cor, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab("r Pearson's Correlation - Reference r Pearson's Correlation")+
  xlab("K Sensitivity multiplier")+
  theme_bw()

df_summary %>% left_join(sensitivity_references) %>% 
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_K ==1) %>% 
  group_by(model_type) %>% 
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_psi, cor_value-ref_cor, group = model_type))+
  geom_line(aes(sensitivity_psi, cor_value-ref_cor, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab("r Pearson's Correlation - Reference r Pearson's Correlation")+
  xlab("P50 Sensitivity multiplier")+
  theme_bw()


  
  
  
  
#### GPP ####

df %>%
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  filter(model_type %in% c("pmodel","pmodel_ecrit")) %>% 
  mutate(dens = get_density(gpp,E,n = 100)) %>% 
  ggplot(aes(E, gpp, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab(expression(paste("GPP [mol ", m^{-2},h^{-1},"]")))+
  xlab("E modeled [mol m-2 h-1]")+
  theme_bw()
  
  
df %>%
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  filter(model_type %in% c("pmodel","pmodel_ecrit")) %>% 
  mutate(dens = get_density(gpp,swvl,n = 100)) %>% 
  ggplot(aes(swvl, gpp, group = model_type))+
  geom_point(aes(color = dens))+
  # geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab(expression(paste("GPP [mol ", m^{-2},h^{-1},"]")))+
  xlab("SWC [v/v]")+
  theme_bw()

  
  
df %>%
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  filter(gpp_all >= 0, gpp_all<200) %>% 
  # filter(model_type %in% c("pmodel","pmodel_ecrit")) %>% 
  mutate(dens = get_density(gpp_all,E,n = 100)) %>% 
  ggplot(aes(E, gpp_all, group = model_type))+
  geom_point(aes(color = dens))+
  # geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab(expression(paste("GPP [mol ", m^{-2},h^{-1},"]")))+
  xlab("E modeled [mol m-2 h-1]")+
  theme_bw() 

df %>%
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  filter(gpp_all >= 0, gpp_all<200) %>% 
  # filter(model_type %in% c("pmodel","pmodel_ecrit")) %>% 
  mutate(dens = get_density(gpp_all,swvl,n = 100)) %>% 
  ggplot(aes(swvl, gpp_all, group = model_type))+
  geom_point(aes(color = dens))+
  # geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab(expression(paste("GPP [mol ", m^{-2},h^{-1},"]")))+
  xlab("SWC [v/v]")+
  theme_bw() 
  
  
  
  
  
  
  
  
  
  
