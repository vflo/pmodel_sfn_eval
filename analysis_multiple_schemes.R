################################################################################
# PRELIMINARY ANALYSIS
################################################################################
library(lmerTest)
library(emmeans)
library(ggpubr)
library(tidyverse)
library(scales)
library(ggalt)
library(grid)
library(ggConvexHull)
library(rphydro)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)
library(stringr)
source("stomatal_optimization_functions_phydro_calibration.R")
source("gs_stomatal_optimization_functions_phydro_calibration.R")
source('hydraulic_functions.R')
source('photosynthetic_functions.R')
source("QUADP.R")
makeTransparent = function(col, alpha=0.7){
  rgb(t(col2rgb(col)/255), alpha = alpha)
}

show_col(brewer_pal(palette = "BrBG")(5))

mytheme = function(){
  theme_classic()+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          legend.text=element_text(size=10),
          plot.tag.position = "topleft") 
  
}

mytheme2 = function(){
  theme_classic()+
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=10),
          legend.text=element_text(size=10),
          plot.tag.position = "topleft") 
  
}

col_df <- tibble(scheme = factor(c(#"phydro_wue",
  "phydro_cmax",
  "phydro_cgain", "phydro",
  "phydro_wang_mod","phydro_sox_mod",
  "phydro_sperry"),
  levels = c(#"phydro_wue",
    "phydro_cmax",
    "phydro_cgain", "phydro",
    "phydro_wang_mod","phydro_sox_mod",
    "phydro_sperry"),
  labels = c(#"WUE",
    'CMAX','CGAIN','PHYDRO','PROFITMAX2','SOX','PROFITMAX')),
  col = brewer_pal(palette = "Dark2")(6)
)

traits <- read.csv(file="DATA/imputation_result_plus_subsps.csv")
traits <- traits %>% separate(taxon,c("genus","species", "subsp"),sep="_")
olea_hv <- traits %>% 
  filter(genus == "Olea",species=="europaea",is.na(subsp)) %>% 
  select(Huber.value)
olea_kl <- traits %>% 
  filter(genus == "Olea",species=="europaea",is.na(subsp)) %>% 
  select(KL..kg.m.1.MPa.1.s.1.)

load(file = "DATA/simulations_kmax_vpd.RData")

#### SIMULATIONS RESULTS ####

df <- df %>% 
  filter(dpsi == FALSE) %>%
  mutate(chi = Ciest/ca,
         scheme = factor(scheme, 
                         levels = c(
                           "phydro_cmax",
                           "phydro_cgain", "phydro",
                           "phydro_sox_mod",
                           "phydro_wang_mod",
                           "phydro_sperry"),
                         labels = c(
                           'CMAX','CGAIN','PHYDRO',
                           'SOX','PROFITMAX2','PROFITMAX')),
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "No acclimated"))
  ) %>% 
  filter(
    !Species %in% c('Betula pendula','Pinus sylvestris')
  )


#### PARAMETERS DATA ####
path_par <- "DATA/parameter_kmax_vpd/"
par_data <- list.files(path_par) %>% 
  purrr::map_df(function(x){
    readr::read_csv(paste0(path_par,x))
  })

df_param <- par_data %>% 
  filter(dpsi == FALSE) %>%
  mutate(scheme = factor(scheme, 
                         levels = c(#"phydro_wue",
                                    "phydro_cmax",
                                    "phydro_cgain", "phydro",
                                    # "phydro_wang","phydro_sox",
                                    "phydro_wang_mod","phydro_sox_mod",
                                    "phydro_sperry"),
                         labels = c(#"WUE",
                                    'CMAX','CGAIN','PHYDRO',
                                    # 'PROFITMAX2net','SOXnet',
                                    'PROFITMAX2','SOX','PROFITMAX')) ,
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "No acclimated"))
  )%>% 
  filter(
    # !scheme %in% c('PROFITMAX2net',
    #                     'SOXnet',
    #                     'WUE'),
         !Species %in% c('Betula pendula','Pinus sylvestris'))
df_param[which(df_param$Species %in% 
                 c("Olea europaea var. Chemlali","Olea europaea var. Meski")),
         "Huber.value"] <- olea_hv$Huber.value
df_param[which(df_param$Species %in% 
                 c("Olea europaea var. Chemlali","Olea europaea var. Meski")),
         "KL..kg.m.1.MPa.1.s.1."] <- olea_kl$KL..kg.m.1.MPa.1.s.1.

path_par_kmax_alpha <- "DATA/parameter_kmax_alpha_vpd/"
par_data_kmax_alpha <- list.files(path_par_kmax_alpha) %>% 
  purrr::map_df(function(x){
    readr::read_csv(paste0(path_par_kmax_alpha,x))
  })


df_param_kmax_alpha <- par_data_kmax_alpha %>% 
  filter(dpsi == FALSE) %>%
  mutate(stomatal_model = scheme,
         scheme = factor(scheme, 
                         levels = c("phydro_wue","phydro_cmax",
                                    "phydro_cgain", "phydro",
                                    # "phydro_wang","phydro_sox",
                                    "phydro_wang_mod","phydro_sox_mod",
                                    "phydro_sperry"),
                         labels = c("WUE",'CMAX','CGAIN','PHYDRO',
                                    # 'PROFITMAX2net','SOXnet',
                                    'PROFITMAX2','SOX','PROFITMAX')) ,
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "No acclimated"))
  )%>% 
  filter(
    # !scheme %in% c('PROFITMAX2net',
    #                     'SOXnet',
    #                     'WUE'),
    !Species %in% c('Betula pendula','Pinus sylvestris'))
df_param_kmax_alpha[which(df_param_kmax_alpha$Species %in% 
                 c("Olea europaea var. Chemlali","Olea europaea var. Meski")),
         "Huber.value"] <- olea_hv$Huber.value
df_param_kmax_alpha[which(df_param_kmax_alpha$Species %in% 
                 c("Olea europaea var. Chemlali","Olea europaea var. Meski")),
         "KL..kg.m.1.MPa.1.s.1."] <- olea_kl$KL..kg.m.1.MPa.1.s.1.
df_param_kmax_alpha <- df_param_kmax_alpha %>% 
  left_join(df %>% 
              select(scheme,source,Species,Iabs_growth,D,T,ca) %>%
              group_by(scheme,source,Species) %>%  
              summarise_all(mean,na.rm=TRUE)
            ) %>% 
  mutate(var = 0,
         pmod = purrr::pmap(list(var, T,Iabs_growth,D,ca,K.scale,P50,b,alpha,gamma,stomatal_model), 
                            ~model_numerical(tc =..2, ppfd = ..3, 
                                          vpd = ..4*101325, co2 = ..5, 
                                          elv = 0, fapar = .99, kphio = 0.087, 
                                          psi_soil = ..1, rdark = 0.02, 
                                          par_plant=list(conductivity = ..6*1e-16,
                                                         psi50 = ..7%>% unique(),
                                                         b = ..8%>% unique()), 
                                          par_cost = list(
                                            alpha = ..9,
                                            gamma = ..10
                                          ), stomatal_model = ..11 %>% unique())))%>% 
  unnest_wider(pmod,names_sep = ".")





#### ALPHA ANALYSIS ####
sla_alpha_mod <- lmerTest::lmer(alpha~SLA..cm2.g.1.+ (1|Species)+ (1|scheme) + (1|source), data = df_param_kmax_alpha)
anova(sla_alpha_mod)
sla_alpha <- summary(sla_alpha_mod)$coefficients
sla_sim <- seq(min(df_param_kmax_alpha$SLA..cm2.g.1.,na.rm = TRUE),
    max(df_param_kmax_alpha$SLA..cm2.g.1.,na.rm = TRUE),
    1)
alpha_sim_sla <- sla_alpha[1,1]+sla_alpha[2,1]*sla_sim
alpha_sd_max_sla <- sla_alpha[1,1] + sla_alpha[1,2] + (sla_alpha[2,1]+sla_alpha[2,2])*sla_sim
alpha_sd_min_sla <- sla_alpha[1,1] - sla_alpha[1,2] + (sla_alpha[2,1]-sla_alpha[2,2])*sla_sim
p1 <- ggplot()+
  stat_summary(data = df_param_kmax_alpha,
               mapping=aes(SLA..cm2.g.1.,alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="grey20")+
  geom_line(aes(x = sla_sim, y = alpha_sim_sla),color = "grey20")+
  geom_line(aes(x = sla_sim, y = alpha_sd_max_sla),color = "grey20",linetype=2)+
  geom_line(aes(x = sla_sim, y = alpha_sd_min_sla),color = "grey20",linetype=2)+
  annotate(geom = "text",x=5,y=0.13,hjust=0,vjust=0,label="***",size=6)+
  mytheme2()+
  xlab(expression("SLA  (cm"^2~"g"^-1*")"))+
  ylab(expression(alpha))

p50_alpha_mod <- lmerTest::lmer(alpha~P50 + (1|Species)+ (1|scheme) + (1|source), data = df_param_kmax_alpha)
anova(p50_alpha_mod)
p50_alpha <- summary(p50_alpha_mod)$coefficients
p50_sim <- seq(min(df_param_kmax_alpha$P50,na.rm = TRUE),
               max(df_param_kmax_alpha$P50,na.rm = TRUE),
               1)
alpha_sim_p50 <- p50_alpha[1,1]+p50_alpha[2,1]*p50_sim
alpha_sd_max_p50 <- p50_alpha[1,1] + p50_alpha[1,2] + (p50_alpha[2,1]+p50_alpha[2,2])*p50_sim
alpha_sd_min_p50 <- p50_alpha[1,1] - p50_alpha[1,2] + (p50_alpha[2,1]-p50_alpha[2,2])*p50_sim
p2 <- ggplot()+
  stat_summary(data = df_param_kmax_alpha,
               mapping=aes(P50,alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="grey20")+
  geom_line(aes(x = p50_sim, y = alpha_sim_p50),color = "grey20")+
  geom_line(aes(x = p50_sim, y = alpha_sd_max_p50),color = "grey20",linetype=2)+
  geom_line(aes(x = p50_sim, y = alpha_sd_min_p50),color = "grey20",linetype=2)+
  annotate(geom = "text",x=-7.5,y=0.13,hjust=0,vjust=0,label="**",size=6)+
  mytheme2()+
  xlab(expression(psi[50]~"(MPa)"))+
  ylab(expression(alpha))


hv_alpha_mod <- lmerTest::lmer(alpha~log(Huber.value)+ (1|Species)+ (1|scheme) + (1|source), data = df_param_kmax_alpha)
anova(hv_alpha_mod)
hv_alpha <- summary(hv_alpha_mod)$coefficients
hv_sim <- seq(min(log(df_param_kmax_alpha$Huber.value),na.rm = TRUE),
               max(log(df_param_kmax_alpha$Huber.value),na.rm = TRUE),
               0.1)
alpha_sim_hv<- hv_alpha[1,1]+hv_alpha[2,1]*hv_sim
alpha_sd_max_hv <- hv_alpha[1,1] + hv_alpha[1,2] + (hv_alpha[2,1]+hv_alpha[2,2])*hv_sim
alpha_sd_min_hv <- hv_alpha[1,1] - hv_alpha[1,2] + (hv_alpha[2,1]-hv_alpha[2,2])*hv_sim
p3 <- ggplot()+
  stat_summary(data = df_param_kmax_alpha,
               mapping=aes(log(Huber.value),alpha,group=interaction(Species,source)),fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="grey20")+
  geom_line(aes(x = hv_sim, y = alpha_sim_hv),color = "grey20")+
  geom_line(aes(x = hv_sim, y = alpha_sd_max_hv),color = "grey20",linetype=2)+
  geom_line(aes(x = hv_sim, y = alpha_sd_min_hv),color = "grey20",linetype=2)+
  annotate(geom = "text",x=-7,y=0.13,hjust=0,vjust=0,label="*",size=6)+
  mytheme2()+
  xlab(expression("Huber value (ln("*cm[sw]^2~m[leaf]^-2*"))"))+
  ylab(expression(alpha))

k_alpha_mod <- lmerTest::lmer(alpha~log(K.scale)+ (1|Species)+ (1|scheme) + (1|source), data = df_param_kmax_alpha)
anova(k_alpha_mod) 
k_alpha <- summary(k_alpha_mod)$coefficients
k_sim <- seq(min(log(df_param_kmax_alpha$K.scale),na.rm = TRUE),
              max(log(df_param_kmax_alpha$K.scale),na.rm = TRUE),
              0.1)
alpha_sim_k <- k_alpha[1,1]+k_alpha[2,1]*k_sim
alpha_sd_max_k <- k_alpha[1,1] + k_alpha[1,2] + (k_alpha[2,1]+k_alpha[2,2])*k_sim
alpha_sd_min_k <- k_alpha[1,1] - k_alpha[1,2] + (k_alpha[2,1]-k_alpha[2,2])*k_sim
p4 <- ggplot()+
  stat_summary(data = df_param_kmax_alpha,
               mapping=aes(log(K.scale),alpha,group=interaction(species,source)),
               fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="grey20")+
  # stat_summary(data = df_param_kmax_alpha,
  #              aes(log(55.55*KL..kg.m.1.MPa.1.s.1.),alpha,
  #                  group=interaction(species,source)),
  #              fun.data=mean_sdl, fun.args = list(mult=1),
  #              geom="pointrange", shape=1)+
  # geom_line(aes(x = k_sim, y = alpha_sim_k),color = "grey20")+
  # geom_line(aes(x = k_sim, y = alpha_sd_max_k),color = "grey20",linetype=2)+
  # geom_line(aes(x = k_sim, y = alpha_sd_min_k),color = "grey20",linetype=2)+
  mytheme2()+
  xlab(expression(K[leaf]*" (ln("*mol[H2O]~m^-1~MPa^-1~s^-1*"))"))+
  ylab(expression(alpha))

b_alpha_mod <- lmerTest::lmer(alpha~log(b) + (1|Species)+ (1|scheme) + (1|source), data = df_param_kmax_alpha)
anova(b_alpha_mod)
b_alpha <- summary(b_alpha_mod)$coefficients
b_sim <- seq(min(log(df_param_kmax_alpha$b),na.rm = TRUE),
               max(log(df_param_kmax_alpha$b),na.rm = TRUE),
               1)
alpha_sim_b <- b_alpha[1,1]+b_alpha[2,1]*b_sim
alpha_sd_max_b <- b_alpha[1,1] + b_alpha[1,2] + (b_alpha[2,1]+b_alpha[2,2])*b_sim
alpha_sd_min_b <- b_alpha[1,1] - b_alpha[1,2] + (b_alpha[2,1]-b_alpha[2,2])*b_sim
p5 <- ggplot()+
  stat_summary(data = df_param_kmax_alpha,
               mapping=aes(log(b),alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1),
               geom="pointrange", color="grey20")+
  geom_line(aes(x = b_sim, y = alpha_sim_b),color = "grey20")+
  geom_line(aes(x = b_sim, y = alpha_sd_max_b),color = "grey20",linetype=2)+
  geom_line(aes(x = b_sim, y = alpha_sd_min_b),color = "grey20",linetype=2)+
  # annotate(geom = "text",x=0.1,y=0.13,hjust=0,vjust=0,label="***",size=6)+
  mytheme2()+
  xlab(expression("b"))+
  ylab(expression(alpha))


Iabs_alpha_mod <- lmerTest::lmer(alpha~ Iabs_growth + (1|Species)+ (1|scheme) + (1|source), data = df_param_kmax_alpha)
anova(Iabs_alpha_mod)
Iabs_alpha <- summary(Iabs_alpha_mod)$coefficients
Iabs_sim <- seq(min(df_param_kmax_alpha$Iabs_growth,na.rm = TRUE),
             max(df_param_kmax_alpha$Iabs_growth,na.rm = TRUE),
             1)
alpha_sim_Iabs <- Iabs_alpha[1,1]+Iabs_alpha[2,1]*Iabs_sim
alpha_sd_max_Iabs <- Iabs_alpha[1,1] + Iabs_alpha[1,2] + (Iabs_alpha[2,1]+Iabs_alpha[2,2])*Iabs_sim
alpha_sd_min_Iabs <- Iabs_alpha[1,1] - Iabs_alpha[1,2] + (Iabs_alpha[2,1]-Iabs_alpha[2,2])*Iabs_sim
p6 <- ggplot()+
  stat_summary(data = df_param_kmax_alpha,
               mapping=aes(Iabs_growth,alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1),
               geom="pointrange", color="grey20")+
  geom_line(aes(x = Iabs_sim, y = alpha_sim_Iabs),color = "grey20")+
  geom_line(aes(x = Iabs_sim, y = alpha_sd_max_Iabs),color = "grey20",linetype=2)+
  geom_line(aes(x = Iabs_sim, y = alpha_sd_min_Iabs),color = "grey20",linetype=2)+
  # annotate(geom = "text",x=0.1,y=0.13,hjust=0,vjust=0,label="***",size=6)+
  mytheme2()+
  xlab(expression(I[growth]~"("*mu*"mol"~m^-2~s^-1*")"))+
  ylab(expression(alpha))


sla_jmax0_mod <- lmerTest::lmer(pmod.jmax~SLA..cm2.g.1.+ (1|Species)+ (1|scheme) + (1|source), data = df_param_kmax_alpha)
anova(sla_jmax0_mod)
sla_alpha <- summary(sla_jmax0_mod)

I_jmax0_mod <- lmerTest::lmer(pmod.jmax~Iabs_growth+ (1|Species)+ (1|scheme) + (1|source), data = df_param_kmax_alpha)
anova(I_jmax0_mod)
summary(I_jmax0_mod)

alpha_jmax0_mod <- lmerTest::lmer(pmod.jmax~alpha+ (1|Species)+ (1|scheme) + (1|source), data = df_param_kmax_alpha)
anova(alpha_jmax0_mod)
summary(alpha_jmax0_mod)

ggarrange(p1,p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'),
          ncol=2, nrow = 2)

alpha_mod <- lmerTest::lmer(alpha~ SLA..cm2.g.1.+Iabs_growth+D+ l_b + lk.scale+P50 + (1|Species)+ (1|scheme)+ (1|source), 
                            data = df_param_kmax_alpha %>%
                              mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
                            hv = log(Huber.value),
                            lk.scale =log(K.scale),l_b = log(b))
                            )
summary(alpha_mod)
step(alpha_mod)
alpha_mod <- lmerTest::lmer(#alpha~ SLA..cm2.g.1.+D+ l_b + lk.scale+P50 + (1|scheme)+ (1|source), 
                            alpha ~ SLA..cm2.g.1. + Iabs_growth + l_b + (1 | Species) + (1 | scheme),
                            data = df_param_kmax_alpha %>%
                              mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
                                     hv = log(Huber.value),
                                     lk.scale =log(K.scale),l_b = log(b))
)
summary(alpha_mod)
MuMIn::r.squaredGLMM(alpha_mod)

# library("factoextra")
# df_k_a_pca <- df_param_kmax_alpha %>% 
#   select(alpha,Huber.value,SLA..cm2.g.1.,K.scale,P50,b) %>% 
#   drop_na 
# res.pca <- prcomp(df_k_a_pca, scale = TRUE)
# fviz_eig(res.pca)
# fviz_pca_ind(res.pca,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )
# fviz_pca_var(res.pca,
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )
# 
# groups <- df_param_kmax_alpha %>% 
#   select(Species,alpha,Huber.value,SLA..cm2.g.1.,K.scale,P50,b) %>% 
#   drop_na %>% 
#   select(Species)
# fviz_pca_ind(res.pca,
#              col.ind = groups$Species, # color by groups
#              # palette = c("#00AFBB",  "#FC4E07"),
#              addEllipses = TRUE, # Concentration ellipses
#              ellipse.type = "confidence",
#              legend.title = "Groups",
#              repel = TRUE
# )

# 
# df_param %>% 
#   ggplot(aes(scheme, log(K.scale), color = acclimation))+
#   geom_boxplot()+
#   mytheme()
# 
# K_par_mod <- lmerTest::lmer(log(K.scale)~scheme*acclimation + (1|Species), data = df_param)
# 
# summary(K_par_mod)
# emmeans(K_par_mod,'acclimation', by='scheme')
# emmeans(K_par_mod,'scheme', by='acclimation')
# pwpp(emmeans(K_par_mod, ~ scheme*acclimation), type = "response")
# plot(emmeans(K_par_mod, ~ scheme*acclimation), comparisons = TRUE)
# pwpp(emmeans(K_par_mod, ~ scheme*acclimation), by = 'acclimation', type = "response")
# pwpp(emmeans(K_par_mod, ~ scheme*acclimation), by = 'scheme', type = "response")
# 






df_param %>% 
  filter(scheme %in% c('PROFITMAX2',
                        'SOX',
                        'PROFITMAX'))%>%
  ggplot(aes(scheme, K.scale, color = acclimation))+
  geom_boxplot()+
  mytheme()+
  scale_color_manual(values = c("#A6611A","#018571"))

df_param %>%
  filter(!scheme %in% c('PROFITMAX2',
                       'SOX',
                       'PROFITMAX'))%>%
  ggplot(aes(scheme, gamma, color = acclimation))+
  geom_boxplot()+
  mytheme()+
  scale_color_manual(values = c("#A6611A","#018571"))

# df_param %>%
#   ggplot(aes(scheme, -psi50, color = acclimation))+
#   geom_boxplot()+
#   mytheme()+
#   scale_color_manual(values = c("#A6611A","#018571"))


# p50_par_mod <- lmerTest::lmer(psi50~scheme*acclimation + (1|Species), data = df_param)
# 
# summary(p50_par_mod)
# emmeans::emmeans(p50_par_mod,'acclimation', by='scheme')
# emmeans::emmeans(p50_par_mod,'scheme', by='acclimation')
# emmeans::pwpp(emmeans::emmeans(p50_par_mod, ~ scheme*acclimation), type = "response")
# 
# df_param %>% 
#   filter(scheme %in% c("phydro_wang_mod", "phydro_sox_mod","phydro_sperry")) %>% 
#   ggplot(aes(scheme, log(b), color = acclimation))+
#   geom_boxplot()+
#   mytheme()
# 
# df_param %>% 
#   filter(acclimation == "Acclimated") %>% 
#   ggplot(aes(scheme, alpha))+
#   geom_boxplot()+
#   mytheme()
# 
# df_param %>% 
#   filter(!scheme %in% c("phydro_wang_mod", "phydro_sox_mod","phydro_sperry")) %>% 
#   ggplot(aes(scheme, gamma, color = acclimation))+
#   geom_boxplot()+
#   mytheme()





################################################################################
#### A
################################################################################
df_a <- df %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(!is.na(A)) %>% 
  mutate(diff_a = A - a_pred) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2])# %>% 
  # filter(beta>= -40)



r_a <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_a, weights = n_dist
                      )
r_a_p <- emmeans::contrast(emmeans(r_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
r_a <- emmeans(r_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  left_join(r_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_a %>%
  ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Assimilation rate r Pearson's correlation"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



beta_a <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_a, weights = n_dist
                         )
beta_a_p <- emmeans::contrast(emmeans(beta_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
beta_a <- emmeans(beta_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(beta_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_a%>% 
  ggplot(aes(scheme,beta, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 1, slope = 0, color = "grey20")+
  geom_pointrange(data = beta_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylim(-10,30)+
  ylab(expression(' Assimilation rate'~beta))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()+
  # scale_y_log10()+
  NULL



rmse_a <- lmerTest::lmer(rmse~scheme*acclimation + (1|Species), data = df_a, weights = n_dist
                         )
rmse_a_p <- emmeans::contrast(emmeans(rmse_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
rmse_a <- emmeans(rmse_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(rmse_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p3 <- df_a %>%
  ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_pointrange(data = rmse_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Assimilation rate RMSE  ("*mu*"mol m"^-2~"s"^-1*")"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()

bias_a <- lmerTest::lmer(bias~scheme*acclimation + (1|Species), data = df_a, weights = n_dist
                         )
bias_a_p <- emmeans::contrast(emmeans(bias_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
bias_a <- emmeans(bias_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(bias_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_a %>%
  ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = bias_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Assimilation rate BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()


ggarrange(p1,p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'),
          common.legend = T,ncol=2, nrow = 2)



df %>% 
  filter(scheme == "PROFITMAX") %>% 
  ggplot(aes(a_pred,A,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression(textstyle("Predicted A")))+
  ylab(expression(textstyle("Observed A")))+
  ggtitle(expression(atop("PROFITMAX Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))

df %>% 
  filter(scheme == "PROFITMAX2") %>% 
  ggplot(aes(a_pred,A,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression(textstyle("Predicted A")))+
  ylab(expression(textstyle("Observed A")))+
  ggtitle(expression(atop("PROFITMAX2 Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))

df %>% 
  filter(scheme == "SOX") %>% 
  ggplot(aes(a_pred,A,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression(textstyle("Predicted A")))+
  ylab(expression(textstyle("Observed A")))+
  ggtitle(expression(atop("SOX Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))

df %>% 
  filter(scheme == "PHYDRO") %>% 
  ggplot(aes(a_pred,A,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression(textstyle("Predicted A")))+
  ylab(expression(textstyle("Observed A")))+
  ggtitle(expression(atop("PHYDRO Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))


df %>% 
  filter(scheme == "CMAX") %>% 
  ggplot(aes(a_pred,A,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression(textstyle("Predicted A")))+
  ylab(expression(textstyle("Observed A")))+
  ggtitle(expression(atop("CMAX Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))


df %>% 
  filter(scheme == "CGAIN") %>% 
  ggplot(aes(a_pred,A,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression(textstyle("Predicted A")))+
  ylab(expression(textstyle("Observed A")))+
  ggtitle(expression(atop("CMAX Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))




################################################################################
#### G
################################################################################
df_g <- df %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(!is.na(gC)) %>% 
  mutate(diff_g = gC - g_pred) %>% 
  summarise(n_dist = n(),
            r = cor(gC, g_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_g,na.rm = TRUE)/mean(gC,na.rm = TRUE),
            rmse = Metrics::rmse(gC,g_pred),
            beta = lm(gC~g_pred)$coefficients[2]) %>% 
  filter(beta>-10)



r_g <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_g, weights = n_dist)
r_g_p <- emmeans::contrast(emmeans(r_g, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
r_g <- emmeans(r_g,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  left_join(r_g_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_g %>%
  ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_g, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Stomatal conductance r Pearson's correlation"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



beta_g <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_g, weights = n_dist)
beta_g_p <- emmeans::contrast(emmeans(beta_g, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
beta_g <- emmeans(beta_g,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(beta_g_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_g%>% 
  ggplot(aes(scheme, beta, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 1, slope = 0, color = "grey20")+
  geom_pointrange(data = beta_g, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(0,2)+
  ylab(expression('Stomatal conductance'~beta))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



rmse_g <- lmerTest::lmer(rmse~scheme*acclimation + (1|Species), data = df_g, weights = n_dist)
rmse_g_p <- emmeans::contrast(emmeans(rmse_g, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
rmse_g <- emmeans(rmse_g,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(rmse_g_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p3 <- df_g %>%
  ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_pointrange(data = rmse_g, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Stomatal conductance RMSE  (mol m"^-2~"s"^-1*")"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()

bias_g <- lmerTest::lmer(bias~scheme*acclimation + (1|Species), data = df_g, weights = n_dist)
bias_g_p <- emmeans::contrast(emmeans(bias_g, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
bias_g <- emmeans(bias_g,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(bias_g_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_g %>%
  ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = bias_g, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Stomatal conductance BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()


ggarrange(p1,p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'),
          common.legend = T,ncol=2, nrow = 2)



df %>% 
  filter(scheme == "PROFITMAX") %>% 
  ggplot(aes(g_pred,gC,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted g"[s]))+
  ylab(expression("Observed g"[s]))+
  ggtitle(expression(atop("Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))

df %>% 
  filter(scheme == "PROFITMAX2") %>% 
  ggplot(aes(g_pred,gC,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted g"[s]))+
  ylab(expression("Observed g"[s]))+
  ggtitle(expression(atop(" PROFITMAX2 Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))

df %>% 
  filter(scheme == "SOX") %>% 
  ggplot(aes(g_pred,gC,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted g"[s]))+
  ylab(expression("Observed g"[s]))+
  ggtitle(expression(atop("SOX Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))

df %>% 
  filter(scheme == "PHYDRO") %>% 
  ggplot(aes(g_pred,gC,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted g"[s]))+
  ylab(expression("Observed g"[s]))+
  ggtitle(expression(atop("PHYDRO Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))


df %>% 
  filter(scheme == "CGAIN") %>% 
  ggplot(aes(g_pred,gC,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted g"[s]))+
  ylab(expression("Observed g"[s]))+
  ggtitle(expression(atop("CGAIN Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))

df %>% 
  filter(scheme == "CMAX") %>% 
  ggplot(aes(g_pred,gC,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted g"[s]))+
  ylab(expression("Observed g"[s]))+
  ggtitle(expression(atop("CMAX Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))





################################################################################
#### CHI
################################################################################
df_c <- df %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(!is.na(chi)) %>% 
  mutate(diff_c = chi - c_pred) %>% 
  summarise(n_dist = n(),
            r = cor(chi, c_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_c,na.rm = TRUE)/mean(chi,na.rm = TRUE),
            rmse = Metrics::rmse(chi,c_pred),
            beta = lm(chi~c_pred)$coefficients[2]) %>% 
  filter(beta<10)



r_c <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_c, weights = n_dist)
r_c_p <- emmeans::contrast(emmeans(r_c, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
r_c <- emmeans(r_c,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  left_join(r_c_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_c %>%
  ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_c, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio r Pearson's correlation"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



beta_c <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_c, weights = n_dist)
beta_c_p <- emmeans::contrast(emmeans(beta_c, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
beta_c <- emmeans(beta_c,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(beta_c_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_c%>% 
  ggplot(aes(scheme, beta, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 1, slope = 0, color = "grey20")+
  geom_pointrange(data = beta_c, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(0,2)+
  ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio"~beta))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



rmse_c <- lmerTest::lmer(rmse~scheme*acclimation + (1|Species), data = df_c, weights = n_dist)
rmse_c_p <- emmeans::contrast(emmeans(rmse_c, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
rmse_c <- emmeans(rmse_c,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(rmse_c_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p3 <- df_c %>%
  ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_pointrange(data = rmse_c, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio RMSE"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()

bias_c <- lmerTest::lmer(bias~scheme*acclimation + (1|Species), data = df_c, weights = n_dist)
bias_c_p <- emmeans::contrast(emmeans(bias_c, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
bias_c <- emmeans(bias_c,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(bias_c_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_c %>%
  ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = bias_c, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()


ggarrange(p1,p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'),
          common.legend = T,ncol=2, nrow = 2)



df %>% 
  filter(scheme == "PROFITMAX") %>% 
  ggplot(aes(c_pred,chi,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted "~chi))+
  ylab(expression("Observed "~chi))+
  ggtitle(expression(atop("PROFITMAX Leaf internal-to-ambient CO"[2]~"ratio")))

df %>% 
  filter(scheme == "PROFITMAX2") %>% 
  ggplot(aes(c_pred,chi,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted "~chi))+
  ylab(expression("Observed "~chi))+
  ggtitle(expression(atop("PROFITMAX2 Leaf internal-to-ambient CO"[2]~"ratio")))

df %>% 
  filter(scheme == "SOX") %>% 
  ggplot(aes(c_pred,chi,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted "~chi))+
  ylab(expression("Observed "~chi))+
  ggtitle(expression(atop("SOX Leaf internal-to-ambient CO"[2]~"ratio")))

df %>% 
  filter(scheme == "PHYDRO") %>% 
  ggplot(aes(c_pred,chi,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted "~chi))+
  ylab(expression("Observed "~chi))+
  ggtitle(expression(atop("PHYDRO Leaf internal-to-ambient CO"[2]~"ratio")))




################################################################################
#### DPSI
################################################################################
df_d <- df %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(!is.na(Dpsi)) %>% select(Dpsi,d_pred)%>% 
  mutate(diff_d = Dpsi - d_pred) %>% 
  summarise(n_dist = n(),
            r = cor(Dpsi, d_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_d,na.rm = TRUE)/mean(Dpsi,na.rm = TRUE),
            rmse = Metrics::rmse(Dpsi,d_pred),
            beta = lm(Dpsi~d_pred)$coefficients[2]) %>% 
  filter(beta> (-50),beta <50, rmse<100)



r_d <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_d, weights = n_dist)
r_d_p <- emmeans::contrast(emmeans(r_d, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
r_d <- emmeans(r_d,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  left_join(r_d_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_d %>%
  ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_d, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Soil-leaf water-potential difference"~italic(Delta*psi)~"ratio r Pearson's correlation"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



beta_d <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_d, weights = n_dist)
beta_d_p <- emmeans::contrast(emmeans(beta_d, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
beta_d <- emmeans(beta_d,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(beta_d_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_d%>% 
  ggplot(aes(scheme, beta, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 1, slope = 0, color = "grey20")+
  geom_pointrange(data = beta_d, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(0,2)+
  ylab(expression("Soil-leaf water-potential difference"~italic(Delta*psi)~beta))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



rmse_d <- lmerTest::lmer(rmse~scheme*acclimation + (1|Species), data = df_d, weights = n_dist)
rmse_d_p <- emmeans::contrast(emmeans(rmse_d, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
rmse_d <- emmeans(rmse_d,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(rmse_d_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p3 <- df_d %>%
  ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_pointrange(data = rmse_d, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Soil-leaf water-potential difference"~italic(Delta*psi)~"(MPa) RMSE"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()

bias_d <- lmerTest::lmer(bias~scheme*acclimation + (1|Species), data = df_d, weights = n_dist)
bias_d_p <- emmeans::contrast(emmeans(bias_d, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
bias_d <- emmeans(bias_d,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(bias_d_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_d %>%
  ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = bias_d, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Soil-leaf water-potential difference"~italic(Delta*psi)~" BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()


ggarrange(p1,p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'),
          common.legend = T,ncol=2, nrow = 2)



df %>% 
  filter(scheme == "PROFITMAX") %>% 
  ggplot(aes(d_pred,Dpsi,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted "~Delta*psi))+
  ylab(expression("Observed "~Delta*psi))+
  ggtitle(expression(atop("PROFITMAX Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))

df %>% 
  filter(scheme == "PROFITMAX2") %>% 
  ggplot(aes(d_pred,Dpsi,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted "~Delta*psi))+
  ylab(expression("Observed "~Delta*psi))+
  ggtitle(expression(atop("PROFITMAX2 Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))

df %>% 
  filter(scheme == "SOX") %>% 
  ggplot(aes(d_pred,Dpsi,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted "~Delta*psi))+
  ylab(expression("Observed "~Delta*psi))+
  ggtitle(expression(atop("SOX Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))

df %>% 
  filter(scheme == "PHYDRO") %>% 
  ggplot(aes(d_pred,Dpsi,color = acclimation))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, color = "grey20")+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~Species)+
  mytheme2()+
  theme(legend.title = element_blank(),
        legend.position="top",
        plot.title = element_text(vjust = -10))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  xlab(expression("Predicted "~Delta*psi))+
  ylab(expression("Observed "~Delta*psi))+
  ggtitle(expression(atop("PHYDRO Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))




################################################################################
#### SIMULATION + PARAMETERS ANALYSIS
################################################################################
df_a_diff <- df %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(!is.na(A)) %>% 
  select(-c(dpsi,genus,species,subsp,P88..MPa.,P50..MPa.,P12..MPa.,
             SLA..cm2.g.1.,Height.max..m.,Ks..kg.m.1.MPa.1.s.1.,Huber.value,
             KL..kg.m.1.MPa.1.s.1.)) %>% 
  mutate(diff_a = A - a_pred) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2]) %>%
  left_join(df_param%>%
              select(- c(dpsi, inst,genus,species,subsp,P88..MPa.,P50..MPa.,P12..MPa.,
                         SLA..cm2.g.1.,Height.max..m.,Ks..kg.m.1.MPa.1.s.1.,Huber.value,
                         KL..kg.m.1.MPa.1.s.1.))) %>%
  ungroup() %>%
  pivot_wider(names_from = acclimation, values_from = c(5:8,10:14)) %>% 
  mutate(r_diff = r_Acclimated - `r_No acclimated`,
         bias_diff = bias_Acclimated - `bias_No acclimated`,
         rmse_diff = rmse_Acclimated - `rmse_No acclimated`,
         beta_diff = beta_Acclimated - `beta_No acclimated`,
         K.scale_diff = K.scale_Acclimated - `K.scale_No acclimated`,
         psi50_diff = P50_Acclimated - `P50_No acclimated`,
         b_diff = b_Acclimated - `b_No acclimated`,
         gamma_diff = gamma_Acclimated - `gamma_No acclimated`
  ) 

# A - r Pearson's
p1 <- df_a_diff %>% 
  # filter(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")) %>%
  mutate(r_diff = case_when(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")~r_diff,
                            TRUE~NA_real_)) %>% 
  ggplot(aes(K.scale_diff, r_diff, color = scheme, group = scheme))+
  geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
  geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
  # geom_encircle()+
  geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
  geom_point()+
  # geom_density_2d()+
  mytheme2()+
  theme(legend.title = element_blank())+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
  # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
  ylab("")+
  xlab(expression(Delta*"K"))

p2 <- df_a_diff %>% 
  # filter(scheme %in% c("WUE","CGAIN","CMAX",'PHYDRO')) %>% 
  ggplot(aes(gamma_diff, r_diff, color = scheme, group = scheme))+
  geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
  geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
  # geom_encircle()+
  geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
  geom_point()+
  # geom_density_2d()+
  mytheme2()+
  theme(legend.title = element_blank())+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
  # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
  ylab("")+
  xlab(expression(Delta*gamma))


annotate_figure(
  ggarrange(p1,p2, 
          align='hv', labels=c('a', 'b'),
          common.legend = T,ncol=2, nrow = 1),
  left = textGrob("(Acclimated - No acclimated) A rate r Pearson's correlation", 
                  rot = 90, vjust = 1, gp = gpar(cex = 1)))

# A - BETA
p1 <- df_a_diff %>% 
  # filter(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")) %>%
  mutate(beta_diff = case_when(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")~beta_diff,
                               TRUE~NA_real_)) %>% 
  ggplot(aes(K.scale_diff, beta_diff, color = scheme, group = scheme))+
  geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
  geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
  # geom_encircle()+
  geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
  geom_point()+
  # geom_density_2d()+
  mytheme2()+
  theme(legend.title = element_blank())+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
  # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
  ylab("")+
  xlab(expression(Delta*"K"))

p2 <- df_a_diff %>% 
  # filter(scheme %in% c("WUE","CGAIN","CMAX",'PHYDRO')) %>% 
  ggplot(aes(gamma_diff, beta_diff, color = scheme, group = scheme))+
  geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
  geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
  # geom_encircle()+
  geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
  geom_point()+
  # geom_density_2d()+
  mytheme2()+
  theme(legend.title = element_blank())+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
  # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
  ylab("")+
  xlab(expression(Delta*gamma))

annotate_figure(
  ggarrange(p1,p2,
            align='hv', labels=c('a', 'b'),
            common.legend = T,ncol=2, nrow = 1),
  left = textGrob("(Acclimated - No acclimated) A rate \u03B2", 
                  rot = 90, vjust = 1, gp = gpar(cex = 1)))


# A - RMSE
p1 <- df_a_diff %>% 
  # filter(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")) %>%
  mutate(rmse_diff = case_when(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")~rmse_diff,
                            TRUE~NA_real_)) %>% 
  ggplot(aes(K.scale_diff, rmse_diff, color = scheme, group = scheme))+
  geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
  geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
  # geom_encircle()+
  geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
  geom_point()+
  # geom_density_2d()+
  mytheme2()+
  theme(legend.title = element_blank())+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
  # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
  ylab("")+
  xlab(expression(Delta*"K"))

p2 <- df_a_diff %>% 
  # filter(scheme %in% c("WUE","CGAIN","CMAX",'PHYDRO')) %>% 
  ggplot(aes(gamma_diff, rmse_diff, color = scheme, group = scheme))+
  geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
  geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
  # geom_encircle()+
  geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
  geom_point()+
  # geom_density_2d()+
  mytheme2()+
  theme(legend.title = element_blank())+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
  # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
  ylab("")+
  xlab(expression(Delta*gamma))

annotate_figure(
  ggarrange(p1,p2,
            align='hv', labels=c('a', 'b'),
            common.legend = T,ncol=2, nrow = 1),
  left = textGrob("(Acclimated - No acclimated) A rate RMSE", 
                  rot = 90, vjust = 1, gp = gpar(cex = 1)))


# A - BIAS
p1 <- df_a_diff %>% 
  # filter(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")) %>%
  mutate(bias_diff = case_when(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")~bias_diff,
                               TRUE~NA_real_)) %>% 
  ggplot(aes(K.scale_diff, bias_diff, color = scheme, group = scheme))+
  geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
  geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
  # geom_encircle()+
  geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
  geom_point()+
  # geom_density_2d()+
  mytheme2()+
  theme(legend.title = element_blank())+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
  # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
  ylab("")+
  xlab(expression(Delta*"K"))

p2 <- df_a_diff %>% 
  # filter(scheme %in% c("WUE","CGAIN","CMAX",'PHYDRO')) %>% 
  ggplot(aes(gamma_diff, bias_diff, color = scheme, group = scheme))+
  geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
  geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
  # geom_encircle()+
  geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
  geom_point()+
  # geom_density_2d()+
  mytheme2()+
  theme(legend.title = element_blank())+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
  # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
  ylab("")+
  xlab(expression(Delta*gamma))

annotate_figure(
  ggarrange(p1,p2,
            align='hv', labels=c('a', 'b'),
            common.legend = T,ncol=2, nrow = 1),
  left = textGrob("(Acclimated - No acclimated) A rate BIAS", 
                  rot = 90, vjust = 1, gp = gpar(cex = 1)))






################################################################################
get_no_acclimated_response <- function(x,tc_now,ppfd_now,vpd_now,co2_now){
  stomatal_model = x$stomatal_model
  par_plant_now = list(
    conductivity = x$K.scale*1e-16,
    psi50 = x$P50%>% unique(),
    b = x$b%>% unique()
  )
  par_cost_now = list(
    alpha = 0.1,
    gamma = x$gamma
  )
  ndays = 30
  psi_max = 0
  psi_min = -6
  lwp = seq(psi_min, 0, length.out=50)

  dat_a <- rpmodel::rpmodel(tc =tc_now, ppfd = ppfd_now, 
                            vpd = vpd_now, co2 = co2_now, 
                            elv = 0, fapar = .99)
  
  dat1 = tibble(var = lwp, jmax_a=dat_a$jmax, vcmax_a=dat_a$vcmax) %>%
    mutate(p = purrr::pmap(list(var, jmax_a, vcmax_a), 
                           ~model_numerical_instantaneous(tc =tc_now, ppfd = ppfd_now, 
                                                          vpd = vpd_now, co2 = co2_now, 
                                                          elv = 0, fapar = .99, 
                                                          kphio = 0.087, 
                                                          psi_soil = ..1, rdark = 0.02, 
                                                          par_plant=par_plant_now, 
                                                          par_cost = par_cost_now, 
                                                          jmax = ..2, vcmax = ..3, 
                                                          stomatal_model = stomatal_model)) ) %>% 
    unnest_wider(p)
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
  return(dat1 %>% cbind(scheme = x$scheme, Species = x$Species, psi88S=psi88S,psiL88S=psiL88S,gs0 = gs0))
  
}

Qi_params_no <- df_param %>% 
  filter(Species=="Quercus ilex",
         acclimation == "No acclimated",
         source == "Epron and Dreyer (1990)",
         dpsi == FALSE) %>% 
  mutate(stomatal_model = fct_recode(scheme,
                                     phydro_cgain = 'CGAIN',
                                     phydro_cmax = 'CMAX',
                                     phydro = 'PHYDRO',
                                     phydro_sox_mod = 'SOX',
                                     phydro_sperry = 'PROFITMAX',
                                     phydro_wang_mod = 'PROFITMAX2'))

Ep_params_no <- df_param %>% 
  filter(Species=="Eucalyptus pilularis",
         acclimation == "No acclimated",
         dpsi == FALSE) %>% 
  mutate(stomatal_model = fct_recode(scheme,
                                     phydro_cgain = 'CGAIN',
                                     phydro_cmax = 'CMAX',
                                     phydro = 'PHYDRO',
                                     phydro_sox_mod = 'SOX',
                                     phydro_sperry = 'PROFITMAX',
                                     phydro_wang_mod = 'PROFITMAX2'))

Epo_params_no <- df_param %>% 
  filter(Species=="Eucalyptus populnea",
         acclimation == "No acclimated",
         dpsi == FALSE) %>% 
  mutate(stomatal_model = fct_recode(scheme,
                                     phydro_cgain = 'CGAIN',
                                     phydro_cmax = 'CMAX',
                                     phydro = 'PHYDRO',
                                     phydro_sox_mod = 'SOX',
                                     phydro_sperry = 'PROFITMAX',
                                     phydro_wang_mod = 'PROFITMAX2'))


df_sim_Qi_no <- Qi_params_no %>%
  split(seq(nrow(.))) %>% 
  purrr::map_df(get_no_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 400)

df_sim_Ep_no <- Ep_params_no %>% 
  split(seq(nrow(.))) %>% 
  purrr::map_df(get_no_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 400)

df_sim_Epo_no <- Epo_params_no %>% 
  split(seq(nrow(.))) %>% 
  purrr::map_df(get_no_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 400)

df_summary_no <- df_sim_Qi_no  %>% 
  rbind(df_sim_Ep_no) %>% 
  rbind(df_sim_Epo_no) %>% 
  # filter(a>0.2) %>% 
  group_by(Species,scheme) %>% 
  select(psi88S) %>% 
  summarise_all(mean,na.rm=TRUE)



get_acclimated_response <- function(x,tc_now,ppfd_now,vpd_now,co2_now){
  stomatal_model = x$stomatal_model
  par_plant_now = list(
    conductivity = x$K.scale*1e-16,
    psi50 = x$P50%>% unique(),
    b = x$b%>% unique()
  )
  par_cost_now = list(
    alpha = x$alpha,
    gamma = x$gamma
  )
  ndays = 30
  psi_max = 0
  psi_min = -7
  lwp = seq(psi_min, 0, length.out=100)
  day = ndays * (lwp-psi_max)/(psi_min-psi_max)
  
  lwp_day = function(day_num){
    psi_max + day_num/ndays * (psi_min-psi_max)
  }
  
  k = 7
  lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
  
  spl = splinefun(x = max(day):0, y=lwp_week)
  
  dat_acc = tibble(var = spl(day)) %>% 
    mutate(var = case_when(var>0~0,
                           TRUE~var),
           pmod = map(var, ~model_numerical(tc =tc_now, ppfd = ppfd_now, 
                                            vpd = vpd_now, co2 = co2_now, 
                                            elv = 0, fapar = .99, kphio = 0.087, 
                                            psi_soil = ., rdark = 0.02, par_plant=par_plant_now, 
                                            par_cost = par_cost_now, stomatal_model = stomatal_model))) %>% 
    unnest_wider(pmod)
  dat1 = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>%
    mutate(p = purrr::pmap(list(var, jmax_a, vcmax_a), 
                           ~model_numerical_instantaneous(tc =tc_now, ppfd = ppfd_now, 
                                                          vpd = vpd_now, co2 = co2_now, 
                                                          elv = 0, fapar = .99, 
                                                          kphio = 0.087, 
                                                          psi_soil = ..1, rdark = 0.02, 
                                                          par_plant=par_plant_now, 
                                                          par_cost = par_cost_now, 
                                                          jmax = ..2, vcmax = ..3, 
                                                          stomatal_model = stomatal_model)) ) %>% 
    unnest_wider(p)
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
  return(dat1 %>% cbind(scheme = x$scheme, Species = x$Species, psi88S=psi88S,psiL88S=psiL88S,gs0 = gs0))
  
}


Q_params <- df_param %>% 
  filter(Species%in%c("Quercus ilex","Quercus petraea","Quercus pubescens"),
         acclimation == "Acclimated",
         source == "Epron and Dreyer (1990)",
         dpsi == FALSE) %>% 
  mutate(stomatal_model = fct_recode(scheme,
                             phydro_cgain = 'CGAIN',
                             phydro_cmax = 'CMAX',
                             phydro = 'PHYDRO',
                             phydro_sox_mod = 'SOX',
                             phydro_sperry = 'PROFITMAX',
                             phydro_wang_mod = 'PROFITMAX2'))

Q_params_actual <- df_param %>% 
  left_join(df %>% 
              select(Species,ca,T,D,Iabs_used,gC) %>% 
              group_by(Species) %>% 
              summarise_all(mean, na.rm=TRUE)) %>% 
  filter(Species%in%c("Quercus ilex","Quercus petraea","Quercus pubescens"),
         acclimation == "Acclimated",
         source == "Epron and Dreyer (1990)",
         dpsi == FALSE) %>% 
  mutate(stomatal_model = fct_recode(scheme,
                                     phydro_cgain = 'CGAIN',
                                     phydro_cmax = 'CMAX',
                                     phydro = 'PHYDRO',
                                     phydro_sox_mod = 'SOX',
                                     phydro_sperry = 'PROFITMAX',
                                     phydro_wang_mod = 'PROFITMAX2'))

Q_params_actual_kmax_alpha <- df_param_kmax_alpha %>% 
  left_join(df %>% 
              select(Species,ca,T,D,Iabs_used,gC) %>% 
              group_by(Species) %>% 
              summarise_all(mean, na.rm=TRUE)) %>% 
  filter(Species%in%c("Quercus ilex","Quercus petraea","Quercus pubescens"),
         acclimation == "Acclimated",
         source == "Epron and Dreyer (1990)",
         dpsi == FALSE) %>% 
  mutate(stomatal_model = fct_recode(scheme,
                                     phydro_cgain = 'CGAIN',
                                     phydro_cmax = 'CMAX',
                                     phydro = 'PHYDRO',
                                     phydro_sox_mod = 'SOX',
                                     phydro_sperry = 'PROFITMAX',
                                     phydro_wang_mod = 'PROFITMAX2'))
# 
# Ep_params <- df_param %>% 
#   filter(Species=="Eucalyptus pilularis",
#          acclimation == "Acclimated",
#          dpsi == FALSE) %>% 
#   mutate(stomatal_model = fct_recode(scheme,
#                                      phydro_cgain = 'CGAIN',
#                                      phydro_cmax = 'CMAX',
#                                      phydro = 'PHYDRO',
#                                      phydro_sox_mod = 'SOX',
#                                      phydro_sperry = 'PROFITMAX',
#                                      phydro_wang_mod = 'PROFITMAX2'))
# 
# Epo_params <- df_param %>% 
#   filter(Species=="Eucalyptus populnea",
#          acclimation == "Acclimated",
#          dpsi == FALSE) %>% 
#   mutate(stomatal_model = fct_recode(scheme,
#                                      phydro_cgain = 'CGAIN',
#                                      phydro_cmax = 'CMAX',
#                                      phydro = 'PHYDRO',
#                                      phydro_sox_mod = 'SOX',
#                                      phydro_sperry = 'PROFITMAX',
#                                      phydro_wang_mod = 'PROFITMAX2'))


# df_sim_Qi <- Qi_params %>%
#   split(seq(nrow(.))) %>% 
#   purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 4000)
# 
# df_sim_Ep <- Ep_params %>% 
#   split(seq(nrow(.))) %>% 
#   purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 4000)
# 
# df_sim_Epo <- Epo_params%>% 
#   split(seq(nrow(.))) %>% 
#   purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 4000)

# df_test <- df %>% 
#   filter(Species%in%c("Eucalyptus populnea",
#                       "Eucalyptus pilularis",
#                       "Quercus ilex"))
# 
# df_summary <- df_sim_Qi  %>% 
#   rbind(df_sim_Ep) %>% 
#   rbind(df_sim_Epo) %>% 
#   # filter(a>0.2) %>% 
#   group_by(Species,scheme) %>% 
#   select(psi88S) %>% 
#   summarise_all(mean,na.rm=TRUE)


df_sim_Q_low <- Q_params %>%
  split(seq(nrow(.))) %>%
  purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 400)
df_summary_low <- df_sim_Q_low  %>% 
  group_by(Species,scheme) %>% 
  select(psi88S) %>% 
  summarise_all(mean,na.rm=TRUE)

df_sim_Q_high <- Q_params %>%
  split(seq(nrow(.))) %>%
  purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=3000,co2_now = 400)
df_summary_high <- df_sim_Q_high  %>% 
  group_by(Species,scheme) %>% 
  select(psi88S) %>% 
  summarise_all(mean,na.rm=TRUE)

df_sim_Q_actual <- Q_params_actual %>%
  rowwise() %>%
  do(mapped = get_acclimated_response(.,tc_now=.$T,ppfd_now=.$Iabs_used,vpd_now=.$D*101325,co2_now = .$ca)) %>% 
  unnest(cols = c(mapped))
df_summary_actual<- df_sim_Q_actual  %>% 
  group_by(Species,scheme) %>% 
  select(psi88S) %>% 
  summarise_all(mean,na.rm=TRUE)

df_sim_Q_actual_kmax_alpha <- Q_params_actual_kmax_alpha %>%
  rowwise() %>%
  do(mapped = get_acclimated_response(.,tc_now=.$T,ppfd_now=.$Iabs_used,vpd_now=.$D*101325,co2_now = .$ca)) %>% 
  unnest(cols = c(mapped))
df_summary_actual_kmax_alpha <- df_sim_Q_actual_kmax_alpha  %>% 
  group_by(Species,scheme) %>% 
  select(psi88S) %>% 
  summarise_all(mean,na.rm=TRUE)

df_test <- df %>% 
  filter(Species%in%c("Quercus ilex","Quercus petraea","Quercus pubescens"),
         acclimation == "Acclimated",
         source == "Epron and Dreyer (1990)")
df_summary_sp <- df_test %>% 
  group_by(Species) %>% 
  select(P50) %>% 
  summarise_all(mean,na.rm=TRUE)



p1 <- df_sim_Q_low  %>%
  # rbind(df_sim_Ep) %>% 
  # rbind(df_sim_Epo) %>%
  # filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,jmax,color = scheme),size = 1)+
  geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
            linetype = 2,size = 1)+
  # geom_segment(aes(x=psi88S,xend=psi88S, y= 190, yend= 170,color=scheme),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #            data=df_summary_low,
  #            linetype = 1,size=1,
  #            show.legend = FALSE)+
  # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=190,yend=170),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #            linetype = 1,color="grey20",size=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
  scale_colour_manual(breaks = col_df$scheme, 
                    values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~Species)+
  xlim(-7,0)+
  ylim(0,170)+
  NULL

p2 <- df_sim_Q_high  %>%
  # rbind(df_sim_Ep) %>% 
  # rbind(df_sim_Epo) %>%
  # filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,jmax,color = scheme),size = 1)+
  geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
            linetype = 2,size = 1)+
  # geom_segment(aes(x=psi88S,xend=psi88S, y= 190, yend= 170,color=scheme),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #              data=df_summary_high,
  #              linetype = 1,size=1,
  #              show.legend = FALSE)+
  # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=190,yend=170),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #              linetype = 1,color="grey20",size=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~Species)+
  xlim(-7,0)+
  ylim(0,170)+
  NULL

p3 <- df_sim_Q_actual  %>%
  # rbind(df_sim_Ep) %>% 
  # rbind(df_sim_Epo) %>%
  # filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,jmax,color = scheme),size = 1)+
  geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
            linetype = 2,size = 1)+
  # geom_segment(aes(x=psi88S,xend=psi88S, y= 190, yend= 170,color=scheme),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #              data=df_summary_high,
  #              linetype = 1,size=1,
  #              show.legend = FALSE)+
  # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=190,yend=170),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #              linetype = 1,color="grey20",size=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~Species)+
  xlim(-7,0)+
  ylim(0,60)+
  NULL

p4 <- df_sim_Q_actual_kmax_alpha %>%
  # rbind(df_sim_Ep) %>% 
  # rbind(df_sim_Epo) %>%
  # filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,jmax,color = scheme),size = 1)+
  geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
            linetype = 2,size = 1)+
  # geom_segment(aes(x=psi88S,xend=psi88S, y= 190, yend= 170,color=scheme),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #              data=df_summary_high,
  #              linetype = 1,size=1,
  #              show.legend = FALSE)+
  # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=190,yend=170),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #              linetype = 1,color="grey20",size=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~Species)+
  xlim(-7,0)+
  ylim(0,70)+
  NULL

ggarrange(p1,p2, 
          align='hv', labels=c('a', 'b'),
          common.legend = T,ncol=2, nrow = 1)



p1 <- df_sim_Q_low  %>%
  # rbind(df_sim_Ep) %>% 
  # rbind(df_sim_Epo) %>%
  # filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,gs,color = scheme),size = 1)+
  geom_line(aes(p_leaf,gs,color = scheme,group = interaction(Species,scheme)), 
            linetype = 2,size = 1)+
  geom_segment(aes(x=psi88S,xend=psi88S, y= 0.5, yend= 0.45,color=scheme),
               arrow = arrow(length = unit(0.3, "cm")),
               lineend = c('round'),linejoin = c('mitre'),
               data=df_summary_low,
               linetype = 1,size=1,
               show.legend = FALSE)+
  geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.5,yend=0.45),
               arrow = arrow(length = unit(0.3, "cm")),
               lineend = c('round'),linejoin = c('mitre'),
               linetype = 1,color="grey20",size=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~Species)+
  xlim(-7,0)+
  ylim(0,0.6)+
  NULL

p2 <- df_sim_Q_high  %>%
  # rbind(df_sim_Ep) %>% 
  # rbind(df_sim_Epo) %>%
  # filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,gs,color = scheme),size = 1)+
  geom_line(aes(p_leaf,gs,color = scheme,group = interaction(Species,scheme)), 
            linetype = 2,size = 1)+
  geom_segment(aes(x=psi88S,xend=psi88S, y= 0.5, yend= 0.45,color=scheme),
               arrow = arrow(length = unit(0.3, "cm")),
               lineend = c('round'),linejoin = c('mitre'),
               data=df_summary_high,
               linetype = 1,size=1,
               show.legend = FALSE)+
  geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.5,yend=0.45),
               arrow = arrow(length = unit(0.3, "cm")),
               lineend = c('round'),linejoin = c('mitre'),
               linetype = 1,color="grey20",size=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~Species)+
  xlim(-7,0)+
  ylim(0,0.6)+
  NULL

p3 <- df_sim_Q_actual  %>%
  # rbind(df_sim_Ep) %>% 
  # rbind(df_sim_Epo) %>%
  # filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,gs,color = scheme),size = 1)+
  geom_line(aes(p_leaf,gs,color = scheme,group = interaction(Species,scheme)), 
            linetype = 2,size = 1)+
  geom_point(aes(LWP,gC),data=df_test)+
  geom_segment(aes(x=psi88S,xend=psi88S, y= 0.24, yend= 0.22,color=scheme),
               arrow = arrow(length = unit(0.3, "cm")),
               lineend = c('round'),linejoin = c('mitre'),
               data=df_summary_high,
               linetype = 1,size=1,
               show.legend = FALSE)+
  geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.24,yend=0.22),
               arrow = arrow(length = unit(0.3, "cm")),
               lineend = c('round'),linejoin = c('mitre'),
               linetype = 1,color="grey20",size=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~Species)+
  xlim(-7,0)+
  ylim(0,0.31)+
  NULL

ggarrange(p1,p2, 
          align='hv', labels=c('a', 'b'),
          common.legend = T,ncol=2, nrow = 1)

p_dpsi_accl <- df_sim_Qi %>%
  rbind(df_sim_Ep) %>% 
  rbind(df_sim_Epo) %>%
  # filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,dpsi,color = scheme),size = 1)+
  # geom_line(aes(p_leaf,dpsi,color = scheme,group = interaction(Species,scheme)), 
  #           linetype = 2,size = 1)+
  # geom_segment(aes(x=psi88S,xend=psi88S, y= 0.5, yend= 0.45,color=scheme),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #              data=df_summary,
  #              linetype = 1,size=1,
  #              show.legend = FALSE)+
  # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.5,yend=0.45),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #              linetype = 1,color="grey20",size=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(Delta*psi*" (MPa)"))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  facet_wrap(~Species)+
  xlim(-6,0)+
  ylim(0,1.8)+
  NULL


p_dpsi_no_accl <- df_sim_Qi_no %>%
  rbind(df_sim_Ep_no) %>% 
  rbind(df_sim_Epo_no) %>%
  # filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,dpsi,color = scheme),size = 1)+
  # geom_line(aes(p_leaf,dpsi,color = scheme,group = interaction(Species,scheme)), 
  #           linetype = 2,size = 1)+
  # geom_segment(aes(x=psi88S,xend=psi88S, y= 0.5, yend= 0.45,color=scheme),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
  #              data=df_summary,
  #              linetype = 1,size=1,
  #              show.legend = FALSE)+
  # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.5,yend=0.45),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              lineend = c('round'),linejoin = c('mitre'),
#              linetype = 1,color="grey20",size=1)+
mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(Delta*psi*" (MPa)"))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  facet_wrap(~Species)+
  xlim(-6,0)+
  ylim(0,2)+
  NULL

ggarrange(p_dpsi_accl,p_dpsi_no_accl, 
          align='hv', labels=c('a', 'b'),
          common.legend = T,ncol=1, nrow = 2)



df_sim_Qi %>% 
  rbind(df_sim_Ep) %>% 
  rbind(df_sim_Epo) %>% 
  filter(a>0.2) %>% 
  ggplot()+
  geom_point(data=df_test,
             aes(gC,A),
             color = "grey20",shape=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
  ylab(expression(atop("A ("*mu*"mol m"^-2~"s"^-1*")")))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  facet_wrap(~Species)+
  NULL






df_sim_Qi  %>% 
  rbind(df_sim_Ep) %>% 
  rbind(df_sim_Epo) %>% 
  filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(var,a/ci,color = scheme),size = 1)+
  # geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
  #           linetype = 2,size = 1)+
  geom_point(data=df_test,
             aes(LWP,A/(Ciest*101325/1e6)),
             color = "grey20",shape=1)+
  geom_line(aes(gs,a,color = scheme,group = interaction(Species,scheme)),size = 1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(atop("A/"*c[i]~" ("*mu*"mol m"^-2~"s"^-1*"Pa"^-1~")")))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  facet_wrap(~Species)+
  NULL


df_sim_Qi  %>% 
  rbind(df_sim_Ep) %>% 
  rbind(df_sim_Epo) %>% 
  filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(gs,a/ci,color = scheme),size = 1)+
  # geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
  #           linetype = 2,size = 1)+
  geom_point(data=df %>% 
               filter(Species%in%c("Eucalyptus populnea",
                                   "Eucalyptus pilularis",
                                   "Quercus ilex")),
             aes(gC,A/(Ciest*101325/1e6)),
             color = "grey20",shape=1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
  ylab(expression(atop("A/"*c[i]~" ("*mu*"mol m"^-2~"s"^-1*"Pa"^-1~")")))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  facet_wrap(~Species)+
  NULL

df_sim_Qi  %>% 
  rbind(df_sim_Ep) %>% 
  rbind(df_sim_Epo) %>% 
  filter(a>0.2) %>%
  ggplot()+
  geom_line(aes(a,ci,color = scheme),size = 1)+
  # geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
  #           linetype = 2,size = 1)+
  geom_point(data=df %>% 
               filter(Species%in%c("Eucalyptus populnea",
                                   "Eucalyptus pilularis",
                                   "Quercus ilex")),
             aes(A,(Ciest*101325/1e6)),
             color = "grey20",shape=1)+
  mytheme()+
  labs(color = "")+
  ylab(expression(c[i]*" (Pa)"))+
  xlab(expression(atop("A ("*mu*"mol m"^-2~"s"^-1~")")))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  facet_wrap(~Species)+
  NULL

################################################################################

df_sim_Q_low %>%
  # rbind(df_sim_Ep) %>% 
  # rbind(df_sim_Epo) %>% 
  group_by(scheme,Species) %>% 
  mutate(
    rprima1 = a/(lead(ci)-lag(ci)),
    r1 = 1/gs,
    r2 = 1/lead(gs),
    rprima2 = lead(rprima1),
    A1 = a*1e-6,
    A2 = lead(a)*1e-6,
    L = (r2-r1)/(A2-A1),
    sens_g = -(0.5*((A1/(rprima1+r1))+(A2/(rprima2+r2)))*L)) %>% 
  filter(sens_g<=1,gs>=gs0*0.12,sens_g>=0) %>% 
  ggplot()+
  geom_line(aes(var,sens_g,color = scheme,group= scheme),size = 1)+
  geom_vline(data=df_summary_sp, aes(xintercept=P50),linetype=3)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression("Relative stomatal limitation"))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~Species)+
  NULL





################################################################################


vcmax_jmax_params <- df_param %>% 
  left_join(df %>% 
              select(Species,ca,T,D,Iabs_used) %>% 
              group_by(Species) %>% 
              summarise_all(mean, na.rm=TRUE)) %>% 
  filter(Species%in%c("Picea abies","Populus tremula"),
         acclimation == "Acclimated",
         dpsi == FALSE) %>% 
  mutate(stomatal_model = fct_recode(scheme,
                                     phydro_cgain = 'CGAIN',
                                     phydro_cmax = 'CMAX',
                                     phydro = 'PHYDRO',
                                     phydro_sox_mod = 'SOX',
                                     phydro_sperry = 'PROFITMAX',
                                     phydro_wang_mod = 'PROFITMAX2'))


df_sim_jmax <- vcmax_jmax_params %>%
  rowwise() %>%
  do(mapped = get_acclimated_response(.,tc_now=.$T,ppfd_now=.$Iabs_used,vpd_now=.$D*101325,co2_now = .$ca)) %>% 
  unnest(cols = c(mapped))
df_summary_jmax <- df_sim_jmax  %>% 
  group_by(Species,scheme) %>% 
  select(psi88S) %>% 
  summarise_all(mean,na.rm=TRUE)

df_test_jmax <- df %>% 
  filter(Species%in%c("Picea abies","Populus tremula"))

df_summary_jmax <- df_test_jmax %>% 
  group_by(Species) %>% 
  select(P50) %>% 
  summarise_all(mean,na.rm=TRUE)


df_summary_jmax_error <- df %>% filter(Species%in%c("Picea abies","Populus tremula"),
              acclimation == "Acclimated") %>% 
  select(LWP,Species,week,jmax_obs) %>% 
  group_by(Species,week) %>% 
  drop_na() %>% 
  summarise(
    jmax_sd = sd(jmax_obs,na.rm=TRUE),
    LWP_sd = sd(LWP,na.rm=TRUE),
    jmax_obs=mean(jmax_obs,na.rm=TRUE),
    LWP=mean(LWP,na.rm=TRUE)) %>% 
  drop_na()

df_sim_jmax %>%
  ggplot()+
  geom_line(aes(var,jmax,color = scheme),size = 1)+
  geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
            linetype = 2,size = 1)+
  geom_point(data = df %>% filter(Species%in%c("Picea abies","Populus tremula"),
                                  acclimation == "Acclimated") %>%
               group_by(Species,LWP) %>% summarise(jmax_obs=mean(jmax_obs,na.rm=TRUE)),
             mapping=aes(LWP,jmax_obs),
             color="grey50")+
  geom_point(data=df_summary_jmax_error,mapping = aes(LWP,jmax_obs))+
  geom_errorbarh(aes(y=jmax_obs,
                     xmin=LWP-LWP_sd,
                     xmax=LWP+LWP_sd),
                 data=df_summary_jmax_error,
                 height=2)+
  geom_errorbar(aes(x=LWP,
                    ymin=jmax_obs-jmax_sd,
                    ymax=jmax_obs+jmax_sd),
                data=df_summary_jmax_error,
                width=0.1)+
  mytheme()+
  labs(color = "")+
  xlab(expression(psi*" (MPa)"))+
  ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  facet_wrap(~Species)+
  xlim(-7,0)+
  guides(colour = guide_legend(nrow = 1))+
  NULL


################################################################################
#### A ONLY DROUGHT
################################################################################
df_a <- df %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  mutate(LWP_q50 = quantile(LWP, 0.5,na.rm = TRUE)) %>% 
  filter(!is.na(A),LWP<=LWP_q50) %>% 
  mutate(diff_a = A - a_pred) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2]) #%>% 
# filter(beta<40)



r_a <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_a#, weights = n_dist
)
r_a_p <- emmeans::contrast(emmeans(r_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
r_a <- emmeans(r_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  left_join(r_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_a %>%
  ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Assimilation rate r Pearson's correlation"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



beta_a <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_a#, weights = n_dist
)
beta_a_p <- emmeans::contrast(emmeans(beta_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
beta_a <- emmeans(beta_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(beta_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_a%>% 
  ggplot(aes(scheme,beta, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 1, slope = 0, color = "grey20")+
  geom_pointrange(data = beta_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(-10,40)+
  ylab(expression(' Assimilation rate'~beta))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()+
  # scale_y_log10()+
  NULL



rmse_a <- lmerTest::lmer(rmse~scheme*acclimation + (1|Species), data = df_a#, weights = n_dist
)
rmse_a_p <- emmeans::contrast(emmeans(rmse_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
rmse_a <- emmeans(rmse_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(rmse_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p3 <- df_a %>%
  ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_pointrange(data = rmse_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Assimilation rate RMSE  ("*mu*"mol m"^-2~"s"^-1*")"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()

bias_a <- lmerTest::lmer(bias~scheme*acclimation + (1|Species), data = df_a#, weights = n_dist
)
bias_a_p <- emmeans::contrast(emmeans(bias_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
bias_a <- emmeans(bias_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(bias_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_a %>%
  ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = bias_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Assimilation rate BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()


ggarrange(p1,p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'),
          common.legend = T,ncol=2, nrow = 2)




################################################################################
#### CHI
################################################################################
df_c <- df %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  mutate(LWP_q75 = quantile(LWP, 0.75,na.rm = TRUE)) %>% 
  filter(!is.na(chi),LWP<=LWP_q75) %>% 
  mutate(diff_c = chi - c_pred) %>% 
  summarise(n_dist = n(),
            r = cor(chi, c_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_c,na.rm = TRUE)/mean(chi,na.rm = TRUE),
            rmse = Metrics::rmse(chi,c_pred),
            beta = lm(chi~c_pred)$coefficients[2])



r_c <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_c, weights = n_dist)
r_c_p <- emmeans::contrast(emmeans(r_c, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
r_c <- emmeans(r_c,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  left_join(r_c_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_c %>%
  ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_c, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio r Pearson's correlation"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



beta_c <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_c, weights = n_dist)
beta_c_p <- emmeans::contrast(emmeans(beta_c, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
beta_c <- emmeans(beta_c,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(beta_c_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_c%>% 
  ggplot(aes(scheme, beta, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 1, slope = 0, color = "grey20")+
  geom_pointrange(data = beta_c, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(0,2)+
  ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio"~beta))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



rmse_c <- lmerTest::lmer(rmse~scheme*acclimation + (1|Species), data = df_c, weights = n_dist)
rmse_c_p <- emmeans::contrast(emmeans(rmse_c, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
rmse_c <- emmeans(rmse_c,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(rmse_c_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p3 <- df_c %>%
  ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_pointrange(data = rmse_c, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio RMSE"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()

bias_c <- lmerTest::lmer(bias~scheme*acclimation + (1|Species), data = df_c, weights = n_dist)
bias_c_p <- emmeans::contrast(emmeans(bias_c, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
bias_c <- emmeans(bias_c,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(bias_c_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_c %>%
  ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = bias_c, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.3),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()


ggarrange(p1,p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'),
          common.legend = T,ncol=2, nrow = 2)



################################################################################
library(gridExtra)
library(grid)
df_a %>% 
  ungroup() %>% 
  filter(acclimation == "Acclimated") %>% 
  select(Species,scheme,r,source) %>% 
  # pivot_wider(values_from = 'r', names_from = 'scheme') %>% 
  # grid.table(theme = ttheme_minimal(
  #   core=list(bg_params = list(fill = blues9[1:4], col=NA),
  #             fg_params=list(fontface=3))))
ggpubr::ggballoonplot( 
  # aes(x = scheme, y = Species, fill = r)) +
  # geom_tile(color = "grey20") +
  # scale_fill_viridis_b(name = "r Pearson's") +
  # scale_x_discrete(name = "") +
  # scale_y_discrete(name = "")+
  # mytheme2()
x = "scheme", y = "Species",
size = "r", fill = "r") +
  scale_fill_viridis_b(name = "r Pearson's") +
  guides(size = FALSE)
    




################################################################################

df %>% 
  ggplot()+
  geom_point(aes(gC,A),shape = 1)+
  geom_smooth(aes(g_pred,a_pred,color=scheme, group=interaction(Species,scheme)),
              se = FALSE,method="lm",linetype = 2, size=0.5)+
  geom_smooth(aes(g_pred,a_pred,color=scheme),se = TRUE,method="lm")+
  facet_wrap(~acclimation)+
  mytheme2()+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(
    legend.title = element_blank(),
    legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  ylab(expression(atop("A ("*mu*"mol m"^-2~"s"^-1~")")))+
  xlab(expression(atop("g"[s]~"(mol m"^-2~"s"^-1*")")))+
  guides(colour = guide_legend(nrow = 1))+
  NULL



df %>% 
  ggplot(aes(LWP,gC/A,color=Species))+
  geom_point(shape = 1)+
  # geom_smooth(aes(g_pred,a_pred,color=scheme, group=interaction(Species,scheme)),
  #             se = FALSE,method="lm",linetype = 2, size=0.5)+
  geom_smooth(se = FALSE, method="gam",formula = y ~ s(x, bs = "cs",k=3))+
  # facet_wrap(~acclimation)+
  mytheme2()+
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.background = element_rect(fill="transparent"))+
  ylab(expression(atop("WUE"[i]~"("*mu*"mol mol"^{-1}*")")))+
  xlab(expression(psi*" (MPa)"))+
  guides(colour = guide_legend(nrow = 2))+
  NULL
