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
library(padr)
source('stomatal_optimization_functions.R')



ghibli_pal <- ghibli_palette("MononokeMedium", 7, direction = -1, type = c("discrete"))
#### ANALISIS ####


# LOAD METADATA -----------------------------------------------------------

load(file = "DATA/sfn_meta.RData")
# 
# path <- "DATA/PREPARED_DATA_LIST"
# list_files <- list.files(path, full.names = TRUE)
# 
# par_plant <- as.list(list_files) %>%
#   purrr::map_df(function(x){
#     load(x)
#     sfn_list$par_plant_std %>% cbind(si_code = sfn_list$sfn$si_code %>% unique())
#   })
# 
# save(par_plant, file = "DATA/par_plant.RData")
load(file = "DATA/par_plant.RData")

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
  left_join(par_plant %>% dplyr::select(-LAI), by = "si_code") %>% 
  mutate(
         # sensitivity_K = case_when(is.na(sensitivity_K)~1,
         #                           TRUE~sensitivity_K),
         # sensitivity_psi = case_when(is.na(sensitivity_psi)~1,
         #                           TRUE~sensitivity_psi),
         gpp = case_when(model_type %in% c("pmodel")~gpp*1e6,
                         TRUE~gpp),
         a = case_when(model_type %in% c("phydro","phydro_wang","sperry", "wang", "wap")~a*LAI,
                         TRUE~a),
         gpp_all = case_when(is.na(gpp)~a*3600*1e-6, #to mol m-2 h-1
                             TRUE ~ gpp*3600*1e-6),  #to mol m-2 h-1
         log_psi = log(-psi_soil)
         ) %>% 
  filter(!is.na(E_sapflow),!is.na(E), E>0, E_sapflow>0
         ,low_swp==FALSE, model_type != "pmodel_swc",
         E_stand_QF==TRUE, is_st_swc_shallow == 1) %>% 
  group_by(si_code, model_type,sensitivity_K,sensitivity_psi) %>%
  mutate(max_E = max(E),
         max_E_sapflow = max(E_sapflow),
         E_rel = E/max_E,
         E_sapflow_rel = E_sapflow/max_E_sapflow,
         E_dif = E - E_sapflow,
         E_dif_rel = E_rel - E_sapflow_rel,
         max_gs = max(gs),
         max_Gs_sapflow = max(Gs_sapflow),
         gs_rel = gs/max_gs,
         Gs_sapflow_rel = Gs_sapflow/max_Gs_sapflow,
         gs_dif = gs - Gs_sapflow,
         gs_dif_rel = gs_rel - Gs_sapflow_rel) %>% 
  as_tibble()# %>% 
  # dplyr::select(-c(1,3:4,8,14:18,74:85,90:92))

rm("par_plant","sfn_meta")
gc(reset = TRUE)



# WEEKLY DATA ------------------------------------------------------
df_lm <- df %>% 
  group_by(si_code, model_type#,sensitivity_K,sensitivity_psi
           ) %>% 
  do(lm(E_sapflow ~ E, data = .) %>% broom::glance())

df_gs_lm <- df %>% 
  group_by(si_code, model_type#,sensitivity_K,sensitivity_psi
  ) %>% 
  do(lm(Gs_sapflow ~ gs, data = .) %>% broom::glance())

df_summary <- df %>% 
  group_by(si_code, model_type#,sensitivity_K,sensitivity_psi
           ) %>% 
  summarise(n = n(),
            cor_value = weightedCorr(E_sapflow,E, method ="Pearson"),
            rmse = rmse(E_sapflow, E),
            mae = mae(E_sapflow, E),
            bias = bias(E_sapflow, E),
            gs_cor_value = weightedCorr(Gs_sapflow,gs, method ="Pearson"),
            gs_rmse = rmse(Gs_sapflow, gs),
            gs_mae = mae(Gs_sapflow, gs),
            gs_bias = bias(Gs_sapflow, gs),
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
            height = mean(height, na.rm = TRUE),
            gpp_all = mean(gpp_all, na.rm = TRUE)) 

df_summary <- df_summary %>% 
  left_join(df_lm) %>% 
  left_join(df_gs_lm%>% 
              rename_with(~paste0(., "_gs"), 3:14))

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


#### TIMESERIES example ####
df$grp <- format(df$TIMESTAMP, "%Y")
df %>% 
  filter(#sensitivity_K == 1, 
         #sensitivity_psi == 1, 
         model_type %in% c("pmodel", 'phydro','phydro_wang',"sperry","wap","wang"), 
         si_code == "PRT_LEZ_ARN") %>%
  dplyr::select(TIMESTAMP, grp, model_type, E, E_sapflow,E_sapflow_sd) %>% 
  pivot_wider(names_from = model_type, values_from = E) %>% 
  ggplot()+
  geom_ribbon(mapping = aes(TIMESTAMP, 
                            ymax = (E_sapflow + E_sapflow_sd), 
                            ymin = (E_sapflow - E_sapflow_sd),
                            group = grp),
              fill = "grey70")+
  geom_line(aes(TIMESTAMP,wang,group = grp), color = ghibli_pal[[6]])+
  geom_line(aes(TIMESTAMP,wap,group = grp), color = ghibli_pal[[5]])+
  geom_line(aes(TIMESTAMP,sperry,group = grp), color = ghibli_pal[[4]])+
  geom_line(aes(TIMESTAMP,phydro,group = grp), color = ghibli_pal[[3]])+
  geom_line(aes(TIMESTAMP,phydro_wang,group = grp), color = ghibli_pal[[2]])+
  geom_line(aes(TIMESTAMP,pmodel,group = grp), color = ghibli_pal[[1]])+
  geom_path(aes(TIMESTAMP,E_sapflow,group = grp), color = "grey20")+
  # annotate("text",x = lubridate::ymd("2014/01/01"),y = 7, label = "FRA FON")+
  ylab(expression(paste(E[modelled], " [mol ", m[ground]^{-2}, h^{-1},"]")))+
  theme_bw()+
  theme(strip.background = element_blank())


df %>% 
  filter(#sensitivity_K == 1, 
    #sensitivity_psi == 1, 
    model_type %in% c("pmodel", 'phydro','phydro_wang',"sperry","wap","wang")) %>%  
  dplyr::select(TIMESTAMP, grp, model_type, E, E_sapflow,E_sapflow_sd) %>% 
  pivot_wider(names_from = model_type, values_from = E) %>% 
  # lm(pmodel~sperry,data=.) %>% summary()
  filter(!is.na(pmodel), !is.na(phydro)) %>% 
  mutate(dens = get_density(pmodel,phydro,n = 100)) %>% 
  ggplot()+
  geom_point(aes(pmodel,phydro,color = dens))+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  theme_bw()+
  theme(strip.background = element_blank())


df %>% 
  filter(#sensitivity_K == 1, 
    #sensitivity_psi == 1, 
    model_type %in% c("pmodel", 'phydro','phydro_wang',"sperry","wap","wang")) %>%  
  dplyr::select(TIMESTAMP, grp, model_type, E, E_sapflow,E_sapflow_sd) %>% 
  pivot_wider(names_from = model_type, values_from = E) %>% 
  # lm(pmodel~sperry,data=.) %>% summary()
  filter(!is.na(phydro_wang), !is.na(phydro)) %>% 
  mutate(dens = get_density(phydro_wang,phydro,n = 100)) %>% 
  ggplot()+
  geom_point(aes(phydro_wang,phydro,color = dens))+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  theme_bw()+
  theme(strip.background = element_blank())

df %>% 
  filter(#sensitivity_K == 1, 
    #sensitivity_psi == 1, 
    model_type %in% c("pmodel", 'phydro','phydro_wang',"sperry","wap","wang")) %>%  
  dplyr::select(TIMESTAMP, grp, model_type, E, E_sapflow,E_sapflow_sd) %>% 
  pivot_wider(names_from = model_type, values_from = E) %>% 
  # lm(pmodel~sperry,data=.) %>% summary()
  filter(!is.na(pmodel), !is.na(phydro_wang)) %>% 
  mutate(dens = get_density(pmodel,phydro_wang,n = 100)) %>% 
  ggplot()+
  geom_point(aes(pmodel,phydro_wang,color = dens))+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  theme_bw()+
  theme(strip.background = element_blank())


#Gs
df %>% 
  filter(#sensitivity_K == 1, 
    #sensitivity_psi == 1, 
    model_type %in% c("pmodel", 'phydro','phydro_wang',"sperry","wap","wang")) %>%  
  dplyr::select(TIMESTAMP, grp, model_type, gs, Gs_sapflow,Gs_sapflow_sd) %>% 
  pivot_wider(names_from = model_type, values_from = gs) %>% 
  # lm(pmodel~sperry,data=.) %>% summary()
  filter(!is.na(pmodel), !is.na(phydro)) %>% 
  mutate(dens = get_density(pmodel,phydro,n = 100)) %>% 
  ggplot()+
  geom_point(aes(pmodel,phydro,color = dens))+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  theme_bw()+
  theme(strip.background = element_blank())


df %>% 
  filter(#sensitivity_K == 1, 
    #sensitivity_psi == 1, 
    model_type %in% c("pmodel", 'phydro','phydro_wang',"sperry","wap","wang")) %>%  
  dplyr::select(TIMESTAMP, grp, model_type, gs, Gs_sapflow,Gs_sapflow_sd) %>% 
  pivot_wider(names_from = model_type, values_from = gs) %>% 
  # lm(pmodel~sperry,data=.) %>% summary()
  filter(!is.na(phydro_wang), !is.na(phydro)) %>% 
  mutate(dens = get_density(phydro_wang,phydro,n = 100)) %>% 
  ggplot()+
  geom_point(aes(phydro_wang,phydro,color = dens))+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  theme_bw()+
  theme(strip.background = element_blank())

df %>% 
  filter(#sensitivity_K == 1, 
    #sensitivity_psi == 1, 
    model_type %in% c("pmodel", 'phydro','phydro_wang',"sperry","wap","wang")) %>%  
  dplyr::select(TIMESTAMP, grp, model_type, gs, Gs_sapflow,Gs_sapflow_sd) %>% 
  pivot_wider(names_from = model_type, values_from = gs) %>% 
  # lm(pmodel~sperry,data=.) %>% summary()
  filter(!is.na(pmodel), !is.na(phydro_wang)) %>% 
  mutate(dens = get_density(pmodel,phydro_wang,n = 100)) %>% 
  ggplot()+
  geom_point(aes(pmodel,phydro_wang,color = dens))+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  theme_bw()+
  theme(strip.background = element_blank())


#### E vs E modelled ####

mod_E_E <- lmer(E_sapflow~E*model_type+(1|si_code), data = df#%>% 
                  # filter(sensitivity_K == 1, sensitivity_psi == 1)
                  )
summary(mod_E_E)

emt1 <- emtrends(mod_E_E, "model_type", var = "E") %>%
  as_tibble() %>% 
  rename(trend = `E.trend`) %>% 
  mutate(model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG")),
         trend = format(round(trend, 3), nsmall = 3),
         trend = paste("beta"," == ", trend)) %>% 
  cbind(E = 45, E_sapflow =29)
emt1         


mod_E_E_log <- lmer(log_E_sapflow~log_E*model_type+(1|si_code), data = df%>% 
                # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
                  mutate(log_E = log(E), log_E_sapflow = log(E_sapflow)))
summary(mod_E_E)

emt2 <- emtrends(mod_E_E_log, "model_type", var = "log_E")%>%
  as_tibble() %>% 
  rename(trend = `log_E.trend`) %>% 
  mutate(model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG")),
         trend = format(round(trend, 3), nsmall = 3),
         trend = paste("beta"," == ", trend)) %>% 
  cbind(E = 0.01, E_sapflow =0.00001)
emt2         


mod_E_E_rel <- lmer(E_sapflow_rel~E_rel*model_type+(1|si_code), data = df#%>% 
                      # filter(sensitivity_K == 1, sensitivity_psi == 1)
                    )
summary(mod_E_E_rel)

emt3 <- emtrends(mod_E_E_rel, "model_type", var = "E_rel")%>%
  as_tibble() %>% 
  rename(trend = `E_rel.trend`) %>% 
  mutate(model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG")),
         trend = format(round(trend, 3), nsmall = 3),
         trend = paste("beta"," == ", trend)) %>% 
  cbind(E_rel = 0.5, E_sapflow_rel =0.85)
emt3   


p1 <- df %>%
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(log_E = log(E),
         log_E_sapflow = log(E_sapflow),
         dens = get_density(E,E_sapflow,n = 100),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(E, E_sapflow, group = model_type))+
  geom_point(aes(color = dens),show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  geom_text(data = emt1, mapping = aes(label = trend),parse = TRUE)+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab(expression(paste(E[actual], " [mol ", m[ground]^{-2}, h^{-1},"]")))+
  xlab(expression(paste(E[modelled]," [mol ", m[ground]^{-2}, h^{-1},"]")))+
  theme_bw()+
  theme(strip.background = element_blank())

p2 <- df %>%
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(log_E = log(E),
         log_E_sapflow = log(E_sapflow),
         dens = get_density(log_E,log_E_sapflow,n = 100),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(log(E), log(E_sapflow), group = model_type))+
  geom_point(aes(color = dens),show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  geom_text(data = emt2, mapping = aes(label = trend),parse = TRUE)+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab(expression(paste("ln(",E[actual],") [ln(mol ", m[ground]^{-2}, h^{-1},")]")))+
  xlab(expression(paste("ln(",E[modelled],") [ln(mol ", m[ground]^{-2}, h^{-1},")]")))+
  theme_bw()+
  theme(strip.background = element_blank())

p3 <- df %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(E_rel,E_sapflow_rel,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>%  
  ggplot(aes(E_rel*100, E_sapflow_rel*100, group = model_type))+
  geom_point(aes(color = dens),show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  geom_text(data = emt3, mapping = aes(label = trend),parse = TRUE, color = "white")+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ",E[actual]," [%]")))+
  xlab(expression(paste("Relative ",E[modelled]," [%]")))+
  theme_bw()+
  theme(strip.background = element_blank())

ggpubr::ggarrange(p1,p2,p3, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)


#### SWC all data ####
# df %>% 
#   filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
#   group_by(model_type) %>% 
#   mutate(dens = get_density(swvl,E_dif,n = 200),
#          model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
#   ggplot(aes(swvl, E_dif, group = model_type))+
#   geom_point(aes(color = dens),show.legend = FALSE)+
#   geom_abline(intercept = 0, slope = 0, linetype = 2)+
#   geom_smooth(method= "lm", color = ghibli_pal[[3]])+
#   facet_wrap(. ~ model_type)+
#   scale_color_jcolors_contin(palette = "pal2",bias = 2)+
#   ylab(expression(paste(E[modelled], " - ", E[actual], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
#   xlab("SWC [v/v]")+
#   theme_bw()+
#   theme(strip.background = element_blank())
# 
# df %>% 
#   filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
#   group_by(model_type) %>% 
#   mutate(dens = get_density(swvl,E_dif_rel,n = 200),
#          model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
#   ggplot(aes(swvl, E_dif_rel*100, group = model_type))+
#   geom_point(aes(color = dens),show.legend = FALSE)+
#   geom_abline(intercept = 0, slope = 0, linetype = 2)+
#   geom_smooth(method= "lm", color = ghibli_pal[[3]])+
#   facet_wrap(. ~ model_type)+
#   scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
#   ylab(expression(paste("Relative ", E[modelled], " - Relative ", E[actual], " [%]")))+
#   xlab("SWC [v/v]")+
#   theme_bw()+
#   theme(strip.background = element_blank())


#### SOIL PSI all data ####

mod_gs_psi <- lmer((gs-Gs_sapflow)~log_psi*model_type+(1|si_code), data = df#%>% 
                  # filter(sensitivity_K == 1, sensitivity_psi == 1)
                  )
summary(mod_gs_psi)

emt1 <- emtrends(mod_gs_psi, "model_type", var = "log_psi") %>%
  as_tibble() %>% 
  rename(trend = `log_psi.trend`) %>% 
  mutate(model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG")),
         trend = format(round(trend, 3), nsmall = 3),
         trend = paste("beta"," == ", trend)) %>% 
  cbind(E = 45, E_sapflow =29)
emt1         


mod_gs_psi_rel <- lmer((gs_rel-Gs_sapflow_rel)~log_psi*model_type+(1|si_code), data = df#%>% 
                    # filter(sensitivity_K == 1, sensitivity_psi == 1)
                    )
summary(mod_gs_psi_rel)

emt2 <- emtrends(mod_gs_psi_rel, "model_type", var = "log_psi") %>%
  as_tibble() %>% 
  rename(trend = `log_psi.trend`) %>% 
  mutate(model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG")),
         trend = format(round(trend, 3), nsmall = 3),
         trend = paste("beta"," == ", trend)) %>% 
  cbind(E = 45, E_sapflow =29)
emt2 

df %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  filter(!is.na(log_psi)) %>% 
  mutate(dens = get_density(log_psi,gs_dif,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(log_psi, gs_dif, group = model_type))+
  geom_point(aes(color = dens),show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 2)+
  ylab(expression(paste(G[modelled], " - ", G[actual], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab(expression(paste("log(|",psi[soil],"|)")))+
  theme_bw()+
  theme(strip.background = element_blank())

df %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>%
  filter(!is.na(log_psi)) %>%
  mutate(dens = get_density(log_psi,E_dif_rel,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(log_psi, E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens),show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[modelled], " - Relative ", E[actual], " [%]")))+
  xlab(expression(paste("log(|",psi[soil],"|)")))+
  theme_bw()+
  theme(strip.background = element_blank())

# df %>% 
#   filter(sensitivity_K == 1, sensitivity_psi == 1, psi_soil>-20) %>% 
#   group_by(model_type) %>% 
#   filter(!is.na(log_psi)) %>% 
#   mutate(dens = get_density(psi_soil,E_dif,n = 200),
#          model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
#   ggplot(aes(psi_soil, E_dif, group = model_type))+
#   geom_point(aes(color = dens),show.legend = FALSE)+
#   geom_abline(intercept = 0, slope = 0, linetype = 2)+
#   geom_smooth(method= "lm", color = ghibli_pal[[3]])+
#   facet_wrap(. ~ model_type)+
#   scale_color_jcolors_contin(palette = "pal2",bias = 2)+
#   ylab(expression(paste(E[modelled], " - ", E[actual], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
#   xlab(expression(paste("log(|",psi[soil],"|)")))+
#   theme_bw()+
#   theme(strip.background = element_blank())
# 
# df %>% 
#   filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
#   group_by(model_type) %>%
#   filter(!is.na(log_psi)) %>%
#   mutate(dens = get_density(log_psi,E_dif_rel,n = 200),
#          model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
#   ggplot(aes(log_psi, E_dif_rel*100, group = model_type))+
#   geom_point(aes(color = dens),show.legend = FALSE)+
#   geom_abline(intercept = 0, slope = 0, linetype = 2)+
#   geom_smooth(method= "lm", color = ghibli_pal[[3]])+
#   facet_wrap(. ~ model_type)+
#   scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
#   ylab(expression(paste("Relative ", E[modelled], " - Relative ", E[actual], " [%]")))+
#   xlab(expression(paste("log(|",psi[soil],"|)")))+
#   theme_bw()+
#   theme(strip.background = element_blank())

#### VPD all data ####
df %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(vpd,E_dif,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(log(vpd), E_dif, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste(E[modelled], " - ", E[actual], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("log(VPD) [log(kPa)]")+
  theme_bw()+
  theme(strip.background = element_blank())

df %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(vpd,E_dif_rel,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(log(vpd), E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[modelled], " - Relative ", E[actual], " [%]")))+
  xlab("log(vpd) [log(kPa)]")+
  theme_bw()+
  theme(strip.background = element_blank())


#### PPFD all data ####
df %>% 
  filter(#sensitivity_K == 1, sensitivity_psi == 1, 
    ppfd_in<1100) %>%
  group_by(model_type) %>% 
  mutate(dens = get_density(ppfd_in,E_dif,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(ppfd_in, E_dif, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste(E[modelled], " - ", E[actual], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab(expression(paste("PPFD [", mu,"mol ", m^{-2},s^{-1},"]")))+
  theme_bw()+
  theme(strip.background = element_blank())

df %>% 
  filter(#sensitivity_K == 1, sensitivity_psi == 1, 
         ppfd_in<1100) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(ppfd_in,E_dif_rel,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(ppfd_in, E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[modelled], " - Relative ", E[actual], " [%]")))+
  xlab(expression(paste("PPFD [", mu,"mol ", m^{-2},s^{-1},"]")))+
  theme_bw()+
  theme(strip.background = element_blank())


#### Temperature all data ####
df %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(ta,E_dif,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(ta, E_dif, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste(E[modelled], " - ", E[actual], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("T [ºC]")+
  theme_bw()+
  theme(strip.background = element_blank())

df %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(ta,E_dif_rel,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>%
  ggplot(aes(ta, E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[modelled], " - Relative ", E[actual], " [%]")))+
  xlab("T [ºC]")+
  theme_bw()+
  theme(strip.background = element_blank())

#### LAI all data ####
df %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(LAI,E_dif,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(LAI, E_dif, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste(E[modelled], " - ", E[actual], " [", mol," ", m[soil]^{-2}," ", h^{-1},"]")))+
  xlab("LAI")+
  theme_bw()+
  theme(strip.background = element_blank())

df %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  mutate(dens = get_density(LAI,E_dif_rel,n = 200),
         model_type = factor(model_type, c("pmodel","phydro_wang","phydro","sperry","wap", "wang"), c("PMODEL", "PHYDRO CRIT", "PHYDRO", "SPERRY", "WAP", "WANG"))) %>% 
  ggplot(aes(LAI, E_dif_rel*100, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  geom_abline(intercept = 0, slope = 0, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 1.5)+
  ylab(expression(paste("Relative ", E[modelled], " - Relative ", E[actual], " [%]")))+
  xlab("LAI")+
  theme_bw()+
  theme(strip.background = element_blank())


#### MODELS weekly ####

model_pmodel <- lmer(log(E_sapflow) ~ log(E)*(swvl+log(vpd)+ppfd_in+ta+LAI)+(1|si_code), 
                     data = df %>% filter(model_type == "pmodel") %>% mutate(ppfd_in = ppfd_in/1000 ))
model_pmodel <- lmer(E_dif ~ swvl+log(vpd)+ppfd_in+ta+LAI +(1|si_code), 
                     data = df %>% filter(model_type == "pmodel") %>% mutate(ppfd_in = ppfd_in/1000 ))
summary(model_pmodel)
MuMIn::r.squaredGLMM(model_pmodel)


model_phydro_wang <- lmer(log(E_sapflow) ~ log(E)*(swvl+log(vpd)+ppfd_in+ta+LAI)+(1|si_code),
                     data = df %>% filter(model_type == "phydro_wang") %>% mutate(ppfd_in = ppfd_in/1000 ))
model_phydro_wang <- lmer(E_dif ~ swvl+log(vpd)+ppfd_in+ta + (1|si_code),
                         data = df %>% filter(model_type == "phydro_wang") %>% mutate(ppfd_in = ppfd_in/1000 ))
summary(model_phydro_wang)
MuMIn::r.squaredGLMM(model_phydro_wang)


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
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_map, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Mean annual precipitation [mm]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
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
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_mat, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Mean annual temperature [ºC]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
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
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Site elevation [m]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = r.squared, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab(expression(R^{2}))+
  xlab("Site elevation [m]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = bias, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("Bias")+
  xlab("Site elevation [m]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=si_elev, y = rmse, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("RMSE")+
  xlab("Site elevation [m]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
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
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  ggplot(aes(x=VPD, y = cor_value, color = model_type, weight = n))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(1:6)])+
  ylab("r Pearson's Correlation")+
  xlab("Average VPD [kPa]")+
  theme_bw()

df_summary %>% 
  # filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
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
  scale_color_manual(values = ghibli_pal[c(1:7)])+
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
model_pmodel <- lm(bias ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                     data = df_summary %>% filter(model_type == "pmodel"), weights = n)
summary(model_pmodel)
step(model_pmodel) %>% summary()

model_pmodel <- lm(cor_value ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                   data = df_summary %>% filter(model_type == "pmodel"), weights = n)
summary(model_pmodel)
step(model_pmodel) %>% summary()

model_pmodel <- lm(r.squared ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                   data = df_summary %>% filter(model_type == "pmodel"), weights = n)
summary(model_pmodel)
step(model_pmodel) %>% summary()

model_pmodel <- lm(rmse ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                   data = df_summary %>% filter(model_type == "pmodel"), weights = n)
summary(model_pmodel)
step(model_pmodel) %>% summary()
caret::varImp(model_pmodel)



#PMODEL Ecrit
model_phydro_wang <- lm(bias ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                   data = df_summary %>% filter(model_type == "phydro_wang",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_phydro_wang)
step(model_phydro_wang) %>% summary()

model_phydro_wang <- lm(cor_value ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                   data = df_summary %>% filter(model_type == "phydro_wang",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_phydro_wang)
step(model_phydro_wang) %>% summary()

model_phydro_wang <- lm(r.squared ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                   data = df_summary %>% filter(model_type == "phydro_wang",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_phydro_wang)
step(model_phydro_wang) %>% summary()

model_phydro_wang <- lm(rmse ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                   data = df_summary %>% filter(model_type == "phydro_wang",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_phydro_wang)
step(model_phydro_wang) %>% summary()
caret::varImp(model_phydro_wang)


#Sperry
model_sperry <- lm(bias ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                         data = df_summary %>% filter(model_type == "sperry",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_sperry)
step(model_sperry) %>% summary()

model_sperry <- lm(cor_value ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                         data = df_summary %>% filter(model_type == "sperry",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_sperry)
step(model_sperry) %>% summary()

model_sperry <- lm(r.squared ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                         data = df_summary %>% filter(model_type == "sperry",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_sperry)
step(model_sperry) %>% summary()

model_sperry <- lm(rmse ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                         data = df_summary %>% filter(model_type == "sperry",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_sperry)
step(model_sperry) %>% summary()
caret::varImp(model_sperry)


#Sperry
model_phydro <- lm(bias ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                   data = df_summary %>% filter(model_type == "phydro",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_phydro)
step(model_phydro) %>% summary()

model_phydro <- lm(cor_value ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                   data = df_summary %>% filter(model_type == "phydro",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_phydro)
step(model_phydro) %>% summary()

model_phydro <- lm(r.squared ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                   data = df_summary %>% filter(model_type == "phydro",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_phydro)
step(model_phydro) %>% summary()

model_phydro <- lm(rmse ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                   data = df_summary %>% filter(model_type == "phydro",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_phydro)
step(model_phydro) %>% summary()
caret::varImp(model_phydro)


#Wang
model_wang <- lm(bias ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                   data = df_summary %>% filter(model_type == "wang",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_wang)
step(model_wang) %>% summary()

model_wang <- lm(cor_value ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                   data = df_summary %>% filter(model_type == "wang",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_wang)
step(model_wang) %>% summary()

model_wang <- lm(r.squared ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                   data = df_summary %>% filter(model_type == "wang",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_wang)
step(model_wang) %>% summary()

model_wang <- lm(rmse ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                   data = df_summary %>% filter(model_type == "wang",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_wang)
step(model_wang) %>% summary()
caret::varImp(model_wang)


#WAP
model_wap <- lm(bias ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                 data = df_summary %>% filter(model_type == "wap",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_wap)
step(model_wap) %>% summary()

model_wap <- lm(cor_value ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                 data = df_summary %>% filter(model_type == "wap",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_wap)
step(model_wap) %>% summary()

model_wap <- lm(r.squared ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height, 
                 data = df_summary %>% filter(model_type == "wap",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_wap)
step(model_wap) %>% summary()

model_wap <- lm(rmse ~ si_map * si_mat + si_elev + clay + sand + ppfd + LAI + height,
                 data = df_summary %>% filter(model_type == "wap",sensitivity_K == 1, sensitivity_psi == 1), weights = n)
summary(model_wap)
step(model_wap) %>% summary()
caret::varImp(model_wap)



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
  filter(model_type %in% c("pmodel","phydro_wang")) %>% 
  mutate(dens = get_density(gpp,E,n = 100)) %>% 
  ggplot(aes(E, gpp, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
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
  filter(model_type %in% c("pmodel","phydro_wang")) %>% 
  mutate(dens = get_density(gpp,swvl,n = 100)) %>% 
  ggplot(aes(swvl, gpp, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
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
  # filter(model_type %in% c("pmodel","phydro_wang")) %>% 
  mutate(dens = get_density(gpp_all,E,n = 100)) %>% 
  ggplot(aes(E, gpp_all, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
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
  # filter(model_type %in% c("pmodel","phydro_wang")) %>% 
  mutate(dens = get_density(gpp_all,swvl,n = 100)) %>% 
  ggplot(aes(swvl, gpp_all, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  # geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab(expression(paste("GPP [mol ", m^{-2},h^{-1},"]")))+
  xlab("SWC [v/v]")+
  theme_bw() 
  
  
  


#### PSI soil vs PSI leaf
  
df %>%
  filter(sensitivity_K == 1, sensitivity_psi == 1) %>% 
  group_by(model_type) %>% 
  summarise(psi_l = mean(psi_l, na.rm = TRUE))
  filter(psi_soil <= 0, psi_l <= 0, !is.na(psi_soil), !is.na(psi_l)) %>% 
  # filter(model_type %in% c("pmodel","phydro_wang")) %>% 
  # mutate(dens = get_density(psi_l,psi_soil,n = 100)) %>% 
  ggplot(aes(psi_soil, psi_l, group = model_type))+
  geom_point(aes(color = dens), show.legend = FALSE)+
  # geom_abline(intercept = 0, slope = 1, linetype = 2)+
  geom_smooth(method= "lm", color = ghibli_pal[[3]])+
  facet_wrap(. ~ model_type)+
  scale_color_jcolors_contin(palette = "pal2",bias = 3)+
  ylab(expression(paste("GPP [mol ", m^{-2},h^{-1},"]")))+
  xlab("SWC [v/v]")+
  theme_bw() 
  
  
  
  
