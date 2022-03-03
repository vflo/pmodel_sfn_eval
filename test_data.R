library(sapfluxnetr)
library(rpmodel)
library(tidyverse)
library(lubridate)
library(tsibble)

sapply(list('calc_sfn_aggr_E.R','calc_sfn_aggr_env.R','calc_sfn_scale.R'), source, .GlobalEnv)

#### LOAD SAPFLUXNET DATA ####
# path <- "../0.1.5/RData/plant"
# sfn_meta <- sapfluxnetr::read_sfn_metadata(folder = path, .write_cache = TRUE)
# save(sfn_meta, file = "sfn_meta.RData")
# load(file = "sfn_meta.RData")

# 
# FRA_FON <- read_sfn_data("FRA_FON", folder = path) %>%
#   sfn_metrics(period = "1 hour", .funs = list(~ mean(., na.rm = TRUE)),
#               solar = TRUE, interval = "general") %>%
#   sapfluxnetr::metrics_tidyfier(metadata = sfn_meta,interval = "general")
  
# 
# save(FRA_FON, file = "FRA_FON.RData")
load(file = "FRA_FON.RData")

#### SPATIAL AND TEMPORAL DATA AGGREGATION ####
sfn <- calc_sfn_scale(FRA_FON)
env_hour <- calc_sfn_aggr_env(FRA_FON,"1 hour")
sfn_aggr <- calc_sfn_aggr_E(sfn,"1 month")
env <- calc_sfn_aggr_env(FRA_FON,"1 month")

df_hour <- sfn %>% left_join(env_hour)
df <- sfn_aggr %>% 
  left_join(env) %>% 
  mutate(ppfd = ppfd_in/1e6,
         patm = rpmodel::calc_patm(FRA_FON$si_elev %>% unique()))


pm_month <- df %>% 
  split(seq(nrow(.)))%>%
  purrr::map_df(function(x){
  rpmodel::rpmodel(tc=x[["ta"]], vpd=x[["vpd"]]*1000, co2=400, fapar=1, ppfd = x[["ppfd"]], elv=FRA_FON$si_elev %>% unique() )
  })

R =  8.31446261815324 #universal gas constant m3 Pa K−1 mol−1
pm_E <- 1.6*pm_month$gs*(df$vpd*1e3)*3600 #mol m-2 h-1
df <- df %>% bind_cols(pm_E = pm_E) %>% 
  mutate(E_stand_mol= E_stand/18.2) #mol m-2 h-1

df %>% 
  ggplot(aes(E_stand_mol,pm_E))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_smooth(method= "lm", color = "grey30")+
  theme_classic()

#He de comprobar unidades!!!! mirar bien si es gas o 
