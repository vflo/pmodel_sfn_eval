stomatal_model <- "phydro_cgain"
par_plant <- list(conductivity=4.693275e-17,psi50=-2.269907,b=1.895929)
par_cost <- list(alpha=0.1,gamma=14.45469)
par_env <- list(viscosity_water=0.0008900227,density_water=997.047,patm=101325,tc=25,vpd=1000)
par_photosynth <- list(kmm=70.84225,gammastar=4.332,phi0=0.0599865,Iabs=1188,ca=40.53,patm=101325,delta=0.02)
jmax <- seq(-14,4,0.1)
pl <- seq(-30,-0.5,0.1)
par <- expand.grid(jmax, pl)
psi_soil <- -1.71
# psi_soil <- -1.85
# e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (ground) s-1
par %>% 
  split(seq(nrow(.))) %>% 
  purrr::map(function(x){
    x <- as.numeric(x)
    fn_profit(x,
              psi_soil       = psi_soil,
              par_cost       = par_cost,
              par_photosynth = par_photosynth,
              par_plant      = par_plant,
              par_env        = par_env,
              do_optim       = FALSE, 
              stomatal_model = stomatal_model)->res
    return(tibble(jmax = x[1],dpsi = x[2],out = res[1]))
    
  }) %>% bind_rows() -> res



res %>% 
  ggplot()+
  geom_contour_filled(aes(x=exp(dpsi), y = exp(jmax), z = out))+#, show.legend = FALSE)+
  # geom_contour_filled(aes(x=dpsi, y = jmax, z = out))+#, show.legend = FALSE)+
  # geom_point(data = max_out, aes(x=dpsi,y=jmax),shape=3)+
  # facet_wrap(~ psi_soil)+
  theme_bw()+
  theme(strip.background = element_rect(fill="transparent"))+
  ylab(expression(atop(J["max"],"("*mu*"mol/m"^2*"/s)")))+
  xlab(expression(atop(Delta*psi,"(MPa)")))+
  viridis::scale_fill_viridis(discrete=TRUE)
