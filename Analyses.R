library(sapfluxnetr)
library(rpmodel)
library(tidyverse)
library(lmerTest)
library(lubridate)
library(Metrics)
library(ghibli)
library(xts)


#### ANALISIS ####

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

df %>% 
  filter(si_code !="AUS_WOM") %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
                        )) %>%
  ggplot(aes(x=name, y = (value-`E Sap flow`),fill =name))+
  ggbeeswarm::geom_quasirandom(shape =21,size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  coord_flip()+
  ghibli::scale_fill_ghibli_d("MononokeMedium",direction = -1)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow])))+
  xlab("MODEL")+
  facet_wrap(~si_code)+
  theme_classic()

df %>% 
  filter(si_code !="AUS_WOM") %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
  )) %>%
  as_tibble() %>%
  dplyr::select(si_code,name,value,`E Sap flow`) %>% 
  na.omit() %>% 
  group_by(si_code, name) %>% 
  summarise(cor_value = cor(`E Sap flow`,value)) %>% 
  ggplot(aes(x=name, y = cor_value, fill =name))+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  coord_flip()+
  ghibli::scale_fill_ghibli_d("MononokeMedium",direction = -1)+
  ylab("r Pearson's Correlation")+
  xlab("MODEL")+
  facet_wrap(~si_code)+
  theme_bw()

df %>% 
  filter(si_code !="AUS_WOM") %>% 
  pivot_longer(names_to = 'model type', values_to = "Estimate", cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
  )) %>%
  pivot_longer(names_to = 'Environmental variable', values_to = "value", cols = c("vpd","ppfd_in",'FAPAR',"swvl","ta"
  )) %>%
  as_tibble() %>%
  dplyr::select(si_code,`model type`, Estimate,`E Sap flow`,`Environmental variable`, value) %>% 
  na.omit() %>% 
  ggplot(aes(x=value, y = (Estimate -`E Sap flow`), color = `model type`))+
  geom_point()+
  geom_smooth()+
  geom_abline(intercept = 0,slope=0)+
  ghibli::scale_color_ghibli_d("MononokeMedium",direction = -1)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow])))+
  xlab("")+
  facet_wrap(~`Environmental variable`, scales = "free_x")+
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

df %>% 
  filter(si_code !="AUS_WOM") %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
  )) %>%
  ggplot(aes(x=name, y = (value-`E Sap flow`),fill =name))+
  ggbeeswarm::geom_quasirandom(shape =21,size=2, dodge.width = .75,alpha=.5,show.legend = F)+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  coord_flip()+
  ghibli::scale_fill_ghibli_d("MononokeMedium",direction = -1)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow])))+
  xlab("MODEL")+
  facet_wrap(~si_code)+
  theme_classic()

df %>% 
  filter(si_code !="AUS_WOM") %>% 
  pivot_longer(cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
  )) %>%
  as_tibble() %>%
  dplyr::select(si_code,name,value,`E Sap flow`) %>% 
  na.omit() %>% 
  group_by(si_code, name) %>% 
  summarise(cor_value = cor(`E Sap flow`,value)) %>% 
  ggplot(aes(x=name, y = cor_value, fill =name))+
  stat_summary(fun = mean, geom = "point", shape=23, size = 5, show.legend = FALSE)+
  geom_abline(intercept = 0,slope=0)+
  coord_flip()+
  ghibli::scale_fill_ghibli_d("MononokeMedium",direction = -1)+
  ylab("r Pearson's Correlation")+
  xlab("MODEL")+
  facet_wrap(~si_code)+
  theme_bw()

df %>% 
  filter(si_code !="AUS_WOM") %>% 
  pivot_longer(names_to = 'model type', values_to = "Estimate", cols = c("PMODEL","PMODEL swc limitation",'PMODEL Smith'#,'PHydro','PHydro Penman-Monteith'
  )) %>%
  pivot_longer(names_to = 'Environmental variable', values_to = "value", cols = c("vpd","ppfd_in",'FAPAR',"swvl","ta"
  )) %>%
  as_tibble() %>%
  dplyr::select(si_code,`model type`, Estimate,`E Sap flow`,`Environmental variable`, value) %>% 
  na.omit() %>% 
  ggplot(aes(x=value, y = (Estimate -`E Sap flow`), color = `model type`))+
  geom_point( shape = 21)+
  geom_smooth()+
  geom_abline(intercept = 0,slope=0)+
  ghibli::scale_color_ghibli_d("MononokeMedium",direction = -1)+
  ylab(expression(paste(E[Estimate], " - ", E[Sapflow])))+
  xlab("")+
  facet_wrap(~`Environmental variable`, scales = "free_x")+
  theme_bw()
