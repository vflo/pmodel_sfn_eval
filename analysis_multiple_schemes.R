
################################################################################
# PRELIMINARY ANALYSIS
################################################################################
library(lmerTest)
library(emmeans)
library(ggpubr)
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
  filter(dpsi == TRUE) %>%
  # rbind(par_data_extra) %>%
  mutate(scheme = factor(scheme, 
                         levels = c("phydro_wue","phydro_cmax",
                                    "phydro_cgain", "phydro",
                                    "phydro_wang","phydro_sox",
                                    "phydro_wang_mod","phydro_sox_mod",
                                    "phydro_sperry"),
                         labels = c("WUE",'CMAX','CGAIN','PHYDRO','PROFITMAX2net','SOXnet','PROFITMAX2','SOX','PROFITMAX')) 
  )# %>% 
  # filter(!scheme %in% c('PROFITMAX2_alt',
  #                       'SOX_alt'))
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
#   filter(acclimation == TRUE) %>% 
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
  filter(dpsi == TRUE) %>% 
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
  ) #%>% 
  # filter(!scheme %in% c('PROFITMAX2_alt',
  #                       'SOX_alt'))



#### R2 A
df_a <- df %>% 
  group_by(scheme,acclimation,Species) %>% 
  filter(!is.na(A)) %>% 
  mutate(diff_a = A - a_pred) %>% 
  summarise(r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias_scale = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            bias = Metrics::bias(A,a_pred))

r_a <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_a)
emmeans(r_a,'scheme', by='acclimation')
p1 <- df_a %>%
  ggplot(aes(scheme,r,color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  ylab(expression("A r pearson"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  coord_flip()

p2 <- df_a %>%
  ggplot(aes(scheme,bias,color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  ylab(expression("A BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  coord_flip()



p3 <- df_a %>%
  ggplot(aes(scheme,rmse,color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  ylab(expression("A RMSE"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  coord_flip()

p4 <- df_a %>%
  ggplot(aes(scheme,bias_scale,color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  theme(legend.title = element_blank())+
  ylab(expression("A BIAS SCALED"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  coord_flip()


ggarrange(p1,p2,p3,p4,
          align='hv', labels=c('a', 'b','c','d'),
          common.legend = T,ncol=2, nrow = 2)




summary(r_sqrt_a)
car::Anova(r_sqrt_a)

bias_a <- lmerTest::lmer(bias~acclimation + (1|Species) + (1|scheme), data = df_a)

summary(bias_a)
car::Anova(bias_a)
emmeans::contrast(emmeans::emmeans(bias_a,specs = "acclimation",by='scheme'))
emmeans::contrast(emmeans::emmeans(bias_a,specs = "scheme",by='acclimation'))
difflsmeans(bias_a, test.effs = "acclimation", ddf="Kenward-Roger")
pwpm(emmeans(bias_a,  ~ scheme),means = FALSE, flip = TRUE,reverse = TRUE)



#### R2 G
df_g <- df %>% 
  group_by(scheme,acclimation,Species) %>% 
  filter(!is.na(gC)) %>% 
  summarise(r_sqrt = cor(gC, g_pred, use = "pairwise.complete.obs")^2 ,
            bias = Metrics::mape(gC, g_pred))

p3 <- df_g %>%
  ggplot(aes(scheme,r_sqrt,color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange",position=position_dodge(width = 0.3))+
  mytheme2()+
  # xlab("Optimization scheme")+
  ylab(expression("G"[CO2]~"R"^2))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))
p4 <- df_g %>%
  ggplot(aes(scheme,bias,color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  # xlab("Optimization scheme")+
  ylab(expression("G"[CO2]~"BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))


#### R2 CHI
df_c <- df %>% 
  group_by(scheme,acclimation,Species) %>% 
  filter(!is.na(chi)) %>% 
  summarise(r_sqrt = cor(chi, c_pred, use = "pairwise.complete.obs")^2,
            bias = Metrics::mape(chi, c_pred))

p5 <- df_c %>%
  ggplot(aes(scheme,r_sqrt,color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  # xlab("Optimization scheme")+
  ylab(expression(chi~"R"^2))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))

p6 <- df_c %>%
  ggplot(aes(scheme,bias,color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  # xlab("Optimization scheme")+
  ylab(expression(chi~"BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))

#### R2 DPSI
df_d <- df %>% 
  filter(!is.na(Dpsi)) %>% 
  group_by(scheme,acclimation,Species) %>% 
  summarise(r_sqrt = cor(Dpsi, d_pred, use = "pairwise.complete.obs")^2,
            bias = Metrics::mape(Dpsi, d_pred))

p7 <- df_d %>%
  ggplot(aes(scheme,r_sqrt,color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  # xlab("Optimization scheme")+
  ylab(expression(Delta*psi~"R"^2))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))

p8 <- df_d %>%
  ggplot(aes(scheme,log(bias),color = acclimation, group = acclimation))+
  geom_point(shape= 1,position=position_dodge(width = 0.3))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", position=position_dodge(width = 0.3))+
  mytheme2()+
  # xlab("Optimization scheme")+
  ylab(expression(Delta*psi~"BIAS"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))



ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,
          align='hv', labels=c('a', 'b','c','d','e','f','g','h'),
          common.legend = T,ncol=2, nrow = 4)



#########################

df_A_annot_t <- df %>% 
  filter(acclimation == TRUE) %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(A, a_pred, use = "pairwise.complete.obs")^2,
            label = paste('Acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

df_A_annot_n_t <- df %>% 
  filter(acclimation == FALSE) %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(A, a_pred, use = "pairwise.complete.obs")^2,
            label = paste('No~acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

p1 =
  df %>% 
  # filter(acclimation == TRUE) %>%  
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
  filter(acclimation == TRUE) %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(gC, g_pred, use = "pairwise.complete.obs")^2,
            label = paste('Acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

df_g_annot_n_t <- df %>% 
  filter(acclimation == FALSE) %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(gC, g_pred, use = "pairwise.complete.obs")^2,
            label = paste('No~acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

p2 =
  df %>% 
  # filter(acclimation == TRUE) %>%  
  ggplot(mapping = aes(x=g_pred, y=gC, 
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
  filter(acclimation == TRUE) %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(chi, c_pred, use = "pairwise.complete.obs")^2,
            label = paste('Acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

df_c_annot_n_t <- df %>% 
  filter(acclimation == FALSE) %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(chi, c_pred, use = "pairwise.complete.obs")^2,
            label = paste('No~acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

p3 =
  df %>% 
  # filter(acclimation == TRUE) %>%  
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
  filter(acclimation == TRUE, !scheme%in%c('WUE',"CGAIN"),Dpsi>=0) %>%
  group_by(scheme, acclimation) %>% 
  summarise(r_sqrt = cor(Dpsi, d_pred, use = "pairwise.complete.obs")^2,
            label = paste('Acclimated~',"italic(R)^2 ==", sprintf("%.2f",r_sqrt))) %>%
  summarise_all(unique)

df_d_annot_n_t <- df %>% 
  filter(acclimation == FALSE, !scheme%in%c('WUE',"CGAIN"),Dpsi>=0) %>%
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
