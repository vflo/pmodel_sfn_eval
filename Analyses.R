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
                                   TRUE~sensitivity_psi)) %>% 
  filter(!is.na(E_sapflow),!is.na(E), E>0, E_sapflow>0,low_psi==FALSE) %>% 
  group_by(si_code, model_type,sensitivity_K,sensitivity_psi) %>%
  mutate(max_E = max(E),
         max_E_sapflow = max(E_sapflow),
         E_rel = E/max_E,
         E_sapflow_rel = E_sapflow/max_E_sapflow) %>% 
  as_tibble()



# CORRELATION WEEKLY ------------------------------------------------------
df_lm <- df %>% 
  group_by(si_code, model_type,sensitivity_K,sensitivity_psi) %>% 
  do(lm(E_sapflow~ E, data = .) %>% broom::glance())


df_summary <- df %>% 
  group_by(si_code, model_type,sensitivity_K,sensitivity_psi) %>% 
  summarise(n = n(),
            cor_value = weightedCorr(E_sapflow,E, method ="Pearson"),
            bias = mean(E - E_sapflow, na.rm = TRUE),
            bias_rel = mean(E_rel - E_sapflow_rel, na.rm = TRUE),
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
            ppfd = mean(ppfd_in,na.rm = TRUE)) 

df_summary <- df_summary %>% left_join(df_lm)

sensitivity_references <- df_summary %>% 
  ungroup() %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  mutate(ref_bias = bias,
         ref_rel_bias = bias_rel,
         ref_cor = cor_value,
         ref_r2 = r.squared) %>% 
  dplyr::select(si_code,model_type,ref_bias,ref_rel_bias,ref_cor,ref_r2)

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
    dens = get_density(log_E,log_E_sapflow,n = 100)) %>% 
  ggplot(aes(log(E), log(E_sapflow), group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab("E actual [mol m-2 h-1]")+
  xlab("E modeled [mol m-2 h-1]")+
  theme_bw()

df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(E_rel,E_sapflow_rel,n = 200)) %>% 
  ggplot(aes(E_rel, E_sapflow_rel, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab("Relative E actual")+
  xlab("Relative E modeled")+
  theme_bw()


mod_E_E <- lmer(E_sapflow~E*model_type+(1|si_code), data = df)
summary(mod_E_E)

emt1 <- emtrends(mod_E_E, "model_type", var = "E")
emt1         
pairs(emt1)


mod_E_E_rel <- lmer(E_sapflow_rel~E_rel*model_type+(1|si_code), data = df)
summary(mod_E_E_rel)

emt1_rel <- emtrends(mod_E_E_rel, "model_type", var = "E_rel")
emt1_rel   
pairs(emt1_rel)

#### SWC all data ####
df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(bias = E - E_sapflow,
         bias_rel = E_rel - E_sapflow_rel,
         dens = get_density(swvl,bias,n = 200)) %>% 
  ggplot(aes(swvl, bias, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("SWC [v/v]")+
  theme_bw()

df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(bias = E - E_sapflow,
         bias_rel = E_rel - E_sapflow_rel,
         dens = get_density(swvl,bias_rel,n = 200)) %>% 
  ggplot(aes(swvl, bias_rel, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab("SWC [v/v]")+
  theme_bw()


#### VPD all data ####
df %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(bias = E - E_sapflow,
         bias_rel = E_rel - E_sapflow_rel,
         dens = get_density(vpd,bias,n = 200)) %>% 
  ggplot(aes(log(vpd), bias, group = model_type))+
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
  mutate(bias = E - E_sapflow,
         bias_rel = E_rel - E_sapflow_rel,
         dens = get_density(vpd,bias_rel,n = 200)) %>% 
  ggplot(aes(log(vpd), bias_rel, group = model_type))+
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
  mutate(bias = E - E_sapflow,
         bias_rel = E_rel - E_sapflow_rel,
         dens = get_density(ppfd_in,bias,n = 200)) %>% 
  ggplot(aes(ppfd_in, bias, group = model_type))+
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
  mutate(bias = E - E_sapflow,
         bias_rel = E_rel - E_sapflow_rel,
         dens = get_density(ppfd_in,bias_rel,n = 200)) %>% 
  ggplot(aes(ppfd_in, bias_rel, group = model_type))+
  geom_point(aes(color = dens))+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab(expression(paste("PPFD [", mu,"mol ", m^{-2},s^{-1},"]")))+
  theme_bw()



#### MAP ###
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = bias_rel*100, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = bias_rel*100, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Site elevation [m]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = bias_rel*100, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average VPD [kPa]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=VPD, y = bias_rel*100, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average T [ºC]")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=Ta, y = bias_rel*100, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average PPFD")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=ppfd, y = bias_rel*100, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average SWC")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=swvl, y = bias_rel*100, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab("Average SWC")+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Soil sand percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=sand, y = bias_rel*100, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Soil clay percentage")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=clay, y = bias_rel*100, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab("Soil clay percentage")+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Biomes")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_biome, y = bias_rel*100, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab("Biomes")+
  theme_bw()


#### IGBP ####
df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_igbp, y = cor_value, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE, position = position_dodge(width = .75))+
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
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("IGBP land cover classification")+
  theme_bw()

df_summary %>% 
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_igbp, y = bias_rel*100, color = model_type, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(paste("Relative ", E[Estimate], " - Relative ", E[Sapflow], " [%]")))+
  xlab("IGBP land cover classification")+
  theme_bw()



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
  ggplot(aes(sensitivity_K, bias_rel-ref_rel_bias, group = model_type))+
  geom_line(aes(sensitivity_K, bias_rel-ref_rel_bias, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab("Relative Bias - Reference Relative Bias")+
  xlab("K Sensitivity multiplier")+
  theme_bw()

df_summary %>% left_join(sensitivity_references) %>%
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_K ==1) %>% 
  group_by(model_type) %>%
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_psi, bias_rel-ref_rel_bias, group = model_type))+
  geom_line(aes(sensitivity_psi, bias_rel-ref_rel_bias, group = interaction(si_code, model_type)),alpha = 0.5)+
  stat_summary(fun = mean, fun.min = sd, fun.max = sd, colour = ghibli_pal[[3]])+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  ylab("Relative Bias - Reference Relative Bias")+
  xlab("P50 Sensitivity multiplier")+
  theme_bw()



df_summary %>% left_join(sensitivity_references) %>% 
  filter(!model_type %in% c("pmodel", "pmodel_swc"),!is.na(cor_value),sensitivity_psi ==1) %>% 
  group_by(model_type) %>% 
  # mutate(dens = get_density(E,E_sapflow,n = 100)) %>%
  ggplot(aes(sensitivity_K, r.squared-ref_r2, group = model_type))+
  geom_line(aes(sensitivity_K, r.squared-ref_r2, group = interaction(si_code, model_type)),alpha = 0.5)+
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
  ggplot(aes(sensitivity_psi, r.squared-ref_r2, group = model_type))+
  geom_line(aes(sensitivity_psi, r.squared-ref_r2, group = interaction(si_code, model_type)),alpha = 0.5)+
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


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

# LOAD MONTHLY ------------------------------------------------------------
path_monthly_folder <- "DATA/OUTCOME_MONTHLY"
paths_monthly <- list.files(path = path_monthly_folder,full.names = T)

paths_monthly %>% 
  purrr::map_df(function(x){
    load(file=paste0(x))
    if(!any(colnames(df_res) == "PMODEL swc limitation")){df_res <- df_res %>% mutate(`PMODEL swc limitation` = c(NA))}
    if(!any(colnames(df_res) == 'PHydro')){df_res <- df_res %>% mutate(PHydro = c(NA))}
    if(!any(colnames(df_res) == 'PHydro Penman-Monteith')){df_res <- df_res %>% mutate(`PHydro Penman-Monteith` = c(NA))}
    df_res %>% as_tibble()
  }) ->df

df <- df %>% left_join(sfn_meta$site_md)

# CORRELATION MONTHLY ------------------------------------------------------

df_si <- df %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation")) %>%
  as_tibble() %>%
  dplyr::select(si_code,name, value, `E Sap flow`, si_map, si_mat, si_elev, si_biome, si_igbp, st_sand_perc, st_clay_perc, FAPAR, ta, vpd, swvl, ppfd_in) %>% 
  na.omit() %>% 
  group_by(si_code, name) %>% 
  mutate() %>% 
  summarise(n = n(),
            cor_value = weightedCorr(`E Sap flow`,value, method ="Pearson"),
            bias = mean(value - `E Sap flow`, na.rm = TRUE),
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
            ppfd = mean(ppfd_in,na.rm = TRUE)) 

df_pmodel <- df %>% 
  pivot_longer(cols = c("PMODEL")) %>%
  as_tibble() %>%
  dplyr::select(si_code,name, value, `E Sap flow`, vcmax_pmodel, jmax_pmodel, chi_pmodel, gammastar_pmodel) %>% 
  na.omit() %>% 
  group_by(si_code, name) %>% 
  summarise(n = n(),
            cor_value = weightedCorr(`E Sap flow`,value, method ="Pearson"),
            bias = mean(value - `E Sap flow`, na.rm = TRUE),
            vcmax = mean(vcmax_pmodel,na.rm = TRUE),
            jmax = mean(jmax_pmodel,na.rm = TRUE),
            chi = mean(chi_pmodel,na.rm = TRUE),
            gammastar = mean(gammastar_pmodel,na.rm = TRUE)) 

df_pmodel_swc <- df %>% 
  pivot_longer(cols = c("PMODEL swc limitation")) %>%
  as_tibble() %>%
  dplyr::select(si_code,name, value, `E Sap flow`, vcmax_pmodel_swc, jmax_pmodel_swc, chi_pmodel_swc, gammastar_pmodel_swc) %>% 
  na.omit() %>% 
  group_by(si_code, name) %>% 
  summarise(n = n(),
            cor_value = weightedCorr(`E Sap flow`,value, method ="Pearson"),
            bias = mean(value - `E Sap flow`, na.rm = TRUE),
            vcmax = mean(vcmax_pmodel_swc,na.rm = TRUE),
            jmax = mean(jmax_pmodel_swc,na.rm = TRUE),
            chi = mean(chi_pmodel_swc,na.rm = TRUE),
            gammastar = mean(gammastar_pmodel_swc,na.rm = TRUE)) 

df_pmodel <- df_pmodel %>% bind_rows(df_pmodel_swc)
# 
# df_si %>% 
#   ggplot(aes(x=name, y = cor_value, fill =name))+
#   stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE)+
#   geom_abline(intercept = 0,slope=0)+
#   coord_flip()+
#   scale_fill_manual(values = ghibli_pal[c(2,5)])+
#   ylab("r Pearson's Correlation")+
#   xlab("MODEL")+
#   facet_wrap(~si_code)+
#   theme_bw()

df_si %>% 
  ggplot(aes(x=si_map, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=si_mat, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_pmodel %>% 
  ggplot(aes(x=vcmax, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Average Vcmax")+
  theme_bw()

df_pmodel %>% 
  ggplot(aes(x=jmax, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Average Jmax")+
  theme_bw()

df_pmodel %>% 
  ggplot(aes(x=chi, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Average Chi")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=si_elev, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Site elevation")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=FAPAR, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Average FAPAR")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=VPD, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Average VPD [kPa]")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=Ta, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Average T [ºC]")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=ppfd, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Average PPFD")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=swvl, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Average SWC")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=sand, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Soil sand percentage")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=clay, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Soil clay percentage")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=si_biome, y = cor_value, color = name, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Biomes")+
  theme_bw()

df_si %>% 
  ggplot(aes(x=si_igbp, y = cor_value, color = name, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("IGBP land cover classification")+
  theme_bw()


# BIAS MONTHLY -------------------------------------------------------------
# df_si %>%
#   ggplot(aes(x=name, y = bias,fill =name))+
#   ggbeeswarm::geom_quasirandom(shape =21,size=2, dodge.width = .75,alpha=.5,show.legend = F)+
#   stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE)+
#   geom_abline(intercept = 0,slope=0)+
#   coord_flip()+
#   scale_fill_manual(values =ghibli_pal[c(1:3)])+
#   ylab(expression(paste(E[Estimate], " - ", E[Sapflow])))+
#   xlab("MODEL")+
#   facet_wrap(~si_code)+
#   theme_classic()

df %>% 
  filter(si_code !="AUS_WOM") %>% 
  pivot_longer(names_to = 'model type', values_to = "Estimate", cols = c("PMODEL","PMODEL swc limitation"#,'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
  )) %>%
  pivot_longer(names_to = 'Environmental variable', values_to = "value", cols = c("vpd","ppfd_in",'FAPAR',"swvl","REW","ta")) %>%
  as_tibble() %>%
  dplyr::select(si_code,`model type`, Estimate,`E Sap flow`,`Environmental variable`, value) %>% 
  na.omit() %>% 
  ggplot(aes(x=value, y = (Estimate -`E Sap flow`), color = `model type`))+
  geom_point( shape = 21, alpha = 0.2)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow])))+
  xlab("")+
  facet_wrap(~`Environmental variable`, scales = "free_x")+
  theme_bw()

df_si %>%
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=si_map, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_si%>%
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=si_mat, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_pmodel %>% 
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=vcmax, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average Vcmax")+
  theme_bw()

df_pmodel %>% 
  filter(si_code !="AUS_WOM") %>% 
  ggplot(aes(x=jmax, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average Jmax")+
  theme_bw()

df_si %>% 
  filter(si_code !="AUS_WOM") %>% 
  ggplot(aes(x=si_elev, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Site elevation")+
  theme_bw()

df_si %>%
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=FAPAR, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average FAPAR")+
  theme_bw()

df_si %>% 
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=VPD, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average VPD [kPa]")+
  theme_bw()

df_si %>% 
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=Ta, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average T [ºC]")+
  theme_bw()

df_si %>%
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=ppfd, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average PPFD")+
  theme_bw()

df_si %>% 
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=swvl, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Average SWC")+
  theme_bw()

df_si %>% 
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=sand, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Soil sand percentage")+
  theme_bw()

df_si %>% 
  filter(si_code !="AUS_WOM") %>%
  ggplot(aes(x=clay, y = bias, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Soil clay percentage")+
  theme_bw()

df_si %>% 
  filter(si_code !="AUS_WOM")%>% 
  ggplot(aes(x=si_biome, y = bias, color = name, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("Biomes")+
  theme_bw()

df_si %>% 
  filter(si_code !="AUS_WOM") %>% 
  ggplot(aes(x=si_igbp, y = bias, color = name, weight = n^0.5))+
  ggbeeswarm::geom_quasirandom(size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE, position = position_dodge(width = .75))+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("IGBP land cover classification")+
  theme_bw()



