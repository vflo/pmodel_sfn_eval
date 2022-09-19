
################################################################################
# PRELIMINARY ANALYSIS
################################################################################
library(lmerTest)
library(emmeans)
library(ggpubr)
library(tidyverse)
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

load(file = "DATA/simulations_16_09_2022.RData")

path_par <- "DATA/parameters/"
par_data <- list.files(path_par) %>% 
  purrr::map_df(function(x){
    readr::read_csv(paste0(path_par,x))
  })

# par_data_extra <- par_data %>% 
#   filter(scheme %in% c("phydro_wue","phydro_cgain"
#   )) %>% 
#   filter(Species %in% c( "Allocasuarina luehmannii","Olea europaea var. Meski",
#                          "Pseudotzuga menziesii","Eucalyptus pilularis","Eucalyptus populnea",
#                          "Glycine max","Quercus coccifera",
#                          "Quercus ilex","Quercus suber"
#   ))
df_param <- par_data %>% 
  filter(dpsi == TRUE & scheme %in% c("phydro")|
           dpsi == FALSE & !scheme %in% c("phydro")) %>%
  # filter(dpsi == FALSE) %>%
  # filter(dpsi == TRUE & !scheme %in% c("phydro_wue","phydro_cgain",'phydro_wang_mod','phydro_sox')|
  #          dpsi == FALSE & scheme %in% c("phydro_wue","phydro_cgain",'phydro_wang_mod','phydro_sox')) %>%
  # rbind(par_data_extra) %>%
  mutate(scheme = factor(scheme, 
                         levels = c("phydro_wue","phydro_cmax",
                                    "phydro_cgain", "phydro",
                                    "phydro_wang","phydro_sox",
                                    "phydro_wang_mod","phydro_sox_mod",
                                    "phydro_sperry"),
                         labels = c("WUE",'CMAX','CGAIN','PHYDRO','PROFITMAX2net','SOXnet','PROFITMAX2','SOX','PROFITMAX')) 
  )%>% 
  filter(!scheme %in% c('PROFITMAX2_alt',
                        'SOX_alt'))
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
  ggplot(aes(scheme, log(K.scale), color = acclimation))+
  geom_boxplot()+
  mytheme()+
  scale_color_manual(values = c("#A6611A","#018571"))

df_param %>%
  ggplot(aes(scheme, -psi50, color = acclimation))+
  geom_boxplot()+
  mytheme()+
  scale_color_manual(values = c("#A6611A","#018571"))


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


#### SIMULATIONS RESULTS ####

df <- df %>% 
  filter(dpsi == TRUE & scheme %in% c("phydro")|
           dpsi == FALSE & !scheme %in% c("phydro")) %>%
  # filter(dpsi == FALSE) %>%
  mutate(chi = Ciest/ca,
         scheme = factor(scheme, 
                         levels = c("phydro_wue","phydro_cmax",
                                    "phydro_cgain", "phydro",
                                    "phydro_wang","phydro_sox","phydro_sox_mod",
                                    "phydro_wang_mod",
                                    "phydro_sperry"),
                         labels = c("WUE",'CMAX','CGAIN','PHYDRO','PROFITMAX2net',
                                    'SOXnet','SOX','PROFITMAX2','PROFITMAX')),
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "No acclimated"))
  ) %>% 
  filter(!scheme %in% c('PROFITMAX2net',
                        'SOXnet'))


################################################################################
#### A
################################################################################
df_a <- df %>% 
  group_by(scheme,acclimation,Species) %>% 
  filter(!is.na(A)) %>% 
  mutate(diff_a = A - a_pred) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2]) %>% 
  filter(beta>0)



r_a <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_a, weights = n_dist)
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
  ylab(expression(" Assimilation rate r Pearson's correlation"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



beta_a <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_a, weights = n_dist)
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
  ylim(0,5)+
  ylab(expression(' Assimilation rate'~beta))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



rmse_a <- lmerTest::lmer(rmse~scheme*acclimation + (1|Species), data = df_a, weights = n_dist)
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

bias_a <- lmerTest::lmer(bias~scheme*acclimation + (1|Species), data = df_a, weights = n_dist)
bias_a_p <- emmeans::contrast(emmeans(bias_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
bias_a <- emmeans(bias_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(rmse_a_p) %>% 
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
  scale_shape_manual(values = c(19))+ #since all the pairs have significant difference, select only sig shape
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





################################################################################
#### G
################################################################################
df_g <- df %>% 
  group_by(scheme,acclimation,Species) %>% 
  filter(!is.na(gC)) %>% 
  mutate(diff_g = gC - g_pred) %>% 
  summarise(n_dist = n(),
            r = cor(gC, g_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_g,na.rm = TRUE)/mean(gC,na.rm = TRUE),
            rmse = Metrics::rmse(gC,g_pred),
            beta = lm(gC~g_pred)$coefficients[2]) %>% 
  filter(beta>0, beta<10)



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
  left_join(rmse_g_p) %>% 
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





################################################################################
#### CHI
################################################################################
df_c <- df %>% 
  group_by(scheme,acclimation,Species) %>% 
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
  left_join(rmse_c_p) %>% 
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
  group_by(scheme,acclimation,Species) %>% 
  filter(!is.na(Dpsi)) %>% select(Dpsi,d_pred)%>% 
  mutate(diff_d = Dpsi - d_pred) %>% 
  summarise(n_dist = n(),
            r = cor(Dpsi, d_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_d,na.rm = TRUE)/mean(Dpsi,na.rm = TRUE),
            rmse = Metrics::rmse(Dpsi,d_pred),
            beta = lm(Dpsi~d_pred)$coefficients[2]) %>% 
  filter(beta> (-50),beta <50)



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
  left_join(rmse_d_p) %>% 
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









#########################

df_A_annot_t <- df %>% 
  filter(acclimation == "Acclimated") %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(A, a_pred, use = "pairwise.complete.obs")^2,
            label = paste('Acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

df_A_annot_n_t <- df %>% 
  filter(acclimation == "No acclimated") %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(A, a_pred, use = "pairwise.complete.obs")^2,
            label = paste('No~acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

p1 =
  df %>% 
  # filter(acclimation == "Acclimated") %>%  
  ggplot(mapping = aes(x=a_pred, y=A, 
                       fill = acclimation, 
                       color = acclimation,
                       group = interaction(scheme,acclimation))) + 
  mytheme()+
  geom_point( shape=21, color=makeTransparent("black",0)) + 
  geom_smooth(method = "lm", se = F)+
  geom_abline(slope = 1, intercept=0, color="grey") +
  # scale_fill_gradient(oob=squish, limits = c(-4,0), low = makeTransparent("goldenrod1"), high = makeTransparent("blue"))+
  # scale_fill_viridis_c(oob=squish, limits = c(-1,0), direction = -1, begin=.25)+
  # labs(fill = expression(psi[s]*"/|"*psi[g88]*"|"))+
  geom_text(data = df_A_annot_t,
            aes(x = 1, y = 20, label = label), parse=TRUE,hjust = 0)+
  geom_text(data = df_A_annot_n_t,
            aes(x = 1, y = 17, label = label), parse=TRUE,hjust = 0)+
  facet_wrap(~scheme)+
  xlab(expression(textstyle("Predicted A")))+
  ylab(expression(textstyle("Observed A")))+
  ggtitle(expression(atop("Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))

p1





df_g_annot_t <- df %>% 
  filter(acclimation == "Acclimated") %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(chi, c_pred, use = "pairwise.complete.obs")^2,
            label = paste('Acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

df_g_annot_n_t <- df %>% 
  filter(acclimation == "No acclimated") %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(chi, c_pred, use = "pairwise.complete.obs")^2,
            label = paste('No~acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

p2 =
  df %>% 
  # filter(acclimation == "Acclimated") %>%  
  ggplot(mapping = aes(x=c_pred, y=chi, 
                       fill = acclimation, 
                       color = acclimation,
                       group = interaction(scheme,acclimation))) + 
  mytheme()+
  geom_point( shape=21, color=makeTransparent("black",0)) + 
  geom_smooth(method = "lm", se = F)+
  geom_abline(slope = 1, intercept=0, color="grey") +
  # scale_fill_gradient(oob=squish, limits = c(-4,0), low = makeTransparent("goldenrod1"), high = makeTransparent("blue"))+
  # scale_fill_viridis_c(oob=squish, limits = c(-1,0), direction = -1, begin=.25)+
  # labs(fill = expression(psi[s]*"/|"*psi[g88]*"|"))+
  geom_text(data = df_g_annot_t,
            aes(x = 0, y = 0.3, label = label), parse=TRUE,hjust = 0)+
  geom_text(data = df_g_annot_n_t,
            aes(x = 0, y = 0.27, label = label), parse=TRUE,hjust = 0)+
  facet_wrap(~scheme)+
  xlab(expression("Predicted g"[s]))+
  ylab(expression("Observed g"[s]))+
  ggtitle(expression(atop("Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))

p2





df_c_annot_t <- df %>% 
  filter(acclimation == "Acclimated") %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(chi, c_pred, use = "pairwise.complete.obs")^2,
            label = paste('Acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

df_c_annot_n_t <- df %>% 
  filter(acclimation == "No acclimated") %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(chi, c_pred, use = "pairwise.complete.obs")^2,
            label = paste('No~acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

p3 =
  df %>% 
  # filter(acclimation == "Acclimated") %>%  
  ggplot(mapping = aes(x=c_pred, y=chi, 
                       fill = acclimation, 
                       color = acclimation,
                       group = interaction(scheme,acclimation))) + 
  mytheme()+
  geom_point( shape=21, color=makeTransparent("black",0)) + 
  geom_smooth(method = "lm", se = F)+
  geom_abline(slope = 1, intercept=0, color="grey") +
  # scale_fill_gradient(oob=squish, limits = c(-4,0), low = makeTransparent("goldenrod1"), high = makeTransparent("blue"))+
  # scale_fill_viridis_c(oob=squish, limits = c(-1,0), direction = -1, begin=.25)+
  # labs(fill = expression(psi[s]*"/|"*psi[g88]*"|"))+
  geom_text(data = df_c_annot_t,
            aes(x = 0.7, y = 0.43, label = label), parse=TRUE,hjust = 0)+
  geom_text(data = df_c_annot_n_t,
            aes(x = 0.7, y = 0.37, label = label), parse=TRUE,hjust = 0)+
  facet_wrap(~scheme)+
  xlab(expression("Predicted"~chi))+
  ylab(expression("Observed"~chi))+
  ggtitle(expression(atop("Leaf internal-to-","ambient CO"[2]~"ratio,"~italic(chi))))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))

p3







df_d_annot_t <- df %>% 
  filter(acclimation == "Acclimated", !scheme%in%c('WUE',"CGAIN"),Dpsi>=0) %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(Dpsi, d_pred, use = "pairwise.complete.obs")^2,
            label = paste('Acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

df_d_annot_n_t <- df %>% 
  filter(acclimation == "No acclimated", !scheme%in%c('WUE',"CGAIN"),Dpsi>=0) %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(Dpsi, d_pred, use = "pairwise.complete.obs")^2,
            label = paste('No~acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

p4 =
  df %>% 
  filter(Dpsi>=0, !scheme%in%c('WUE',"CGAIN")) %>%  
  ggplot(mapping = aes(x=d_pred, y=Dpsi, 
                       fill = acclimation, 
                       color = acclimation,
                       group = interaction(scheme,acclimation))) + 
  mytheme()+
  geom_point( shape=21, color=makeTransparent("black",0)) + 
  geom_smooth(method = "lm", se = F)+
  geom_abline(slope = 1, intercept=0, color="grey") +
  # scale_fill_gradient(oob=squish, limits = c(-4,0), low = makeTransparent("goldenrod1"), high = makeTransparent("blue"))+
  # scale_fill_viridis_c(oob=squish, limits = c(-1,0), direction = -1, begin=.25)+
  # labs(fill = expression(psi[s]*"/|"*psi[g88]*"|"))+
  geom_text(data = df_d_annot_t,
            aes(x = 1, y = 0.55, label = label), parse=TRUE,hjust = 0)+
  geom_text(data = df_d_annot_n_t,
            aes(x = 1, y = 0.4, label = label), parse=TRUE,hjust = 0)+
  facet_wrap(~scheme)+
  xlab(expression("Predicted"~Delta*psi))+
  ylab(expression("Observed"~Delta*psi)) + 
  ggtitle(expression(atop("Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))

p4





# png(filename = "calibration_validation_viridis_wtrend_g88_2.png", width = 820*3, height = 650*3, res=300)
cowplot::plot_grid(p1,p5, labels=LETTERS, label_size = 16, label_colour = "grey50", label_x = 0.05, align = "hv", rel_widths = 1, ncol=2)
cowplot::plot_grid(p2,p6,labels=LETTERS, label_size = 16, label_colour = "grey50", label_x = 0.05, align = "hv", rel_widths = 1, ncol=2)
cowplot::plot_grid( p3,p7, labels=LETTERS, label_size = 16, label_colour = "grey50", label_x = 0.05, align = "hv", rel_widths = 1, ncol=2)
cowplot::plot_grid( p4,p8, labels=LETTERS, label_size = 16, label_colour = "grey50", label_x = 0.05, align = "hv", rel_widths = 1, ncol=2)
