library(sapfluxnetr)
library(rpmodel)
library(tidyverse)
library(lmerTest)
library(lubridate)
library(Metrics)
library(ghibli)
library(xts)
library(wCorr)

ghibli_pal <- ghibli_palette("MononokeMedium", 7, direction = -1, type = c("discrete"))
#### ANALISIS ####


# LOAD METADATA -----------------------------------------------------------

load(file = "DATA/sfn_meta.RData")

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




# LOAD WEEKLY ------------------------------------------------------------
path_weekly_folder <- "DATA/OUTCOME_WEEKLY"
paths_weekly <- list.files(path = path_weekly_folder,full.names = T)

paths_weekly %>% 
  purrr::map_df(function(x){
    load(file=paste0(x))
    if(!any(colnames(df_res) == "PMODEL swc limitation")){df_res <- df_res %>% mutate(`PMODEL swc limitation` = c(NA))}
    if(!any(colnames(df_res) == 'PHydro')){df_res <- df_res %>% mutate(PHydro = c(NA))}
    if(!any(colnames(df_res) == 'PHydro Penman-Monteith')){df_res <- df_res %>% mutate(`PHydro Penman-Monteith` = c(NA))}
    df_res %>% as_tibble
  })->df

df <- df %>% left_join(sfn_meta$site_md)



# CORRELATION WEEKLY ------------------------------------------------------

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
  ggplot(aes(x=gammastar, y = cor_value, color = name, weight = n^0.5))+
  geom_point(aes(size = n), show.legend = FALSE)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0,slope=0)+
  scale_color_manual(values = ghibli_pal[c(2,5)])+
  ylab("r Pearson's Correlation")+
  xlab("Average Gammastar")+
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


# BIAS WEEKLY -------------------------------------------------------------
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
