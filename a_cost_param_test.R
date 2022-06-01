library(tidyverse)


get_p <- function(psi_soil, e, K, d, c, h) {
  dp = 0.0
  tension = 0
  # 3.2 iterate through each layer to get stem dp including gravity
  for (i in 1:20){
    p = min(psi_soil, tension)
    f = exp( -1.0 * (p/d)**c )
    layer_k = K *20* f
    layer_k = max(layer_k, 1E-12)
    dp <-  dp + e / layer_k + 998.0*9.8*h*0.05*1E-6
    tension = psi_soil - dp
  }
  
  return(tension)
}

get_e_crit <- function(psi_soil, K, d, c, h){
  p_crit = d * log(1000.0) ** (1.0/c)
  e_min  = 0.0
  e_max  = 100.0
  e_crit = 50.0
  continue <- TRUE
  while(continue){
    p = get_p(psi_soil,e_max, K, d, c, h)
    if (p>p_crit){
      e_max <-  e_max*2
    }else{
      continue <- FALSE}}
  continue <- TRUE
  while(continue){
    e = 0.5 * (e_max+e_min)
    p = get_p(psi_soil,e, K, d, c, h)
    if (abs(p-p_crit)<1E-3 | (e_max-e_min)<1E-3){
      e_crit = e
      continue <- FALSE
    }
    if (p<p_crit){ 
      e_max  = e
    }else{e_min  = e}
    
  }
  return (e_crit)
}


#### SIMULATION values ####

psi_soil <- seq(0,-5,-0.1)
psi_leaf <- seq(0,-5,-0.1)
tc <- seq(0,40,0.1)
K <- 1
d <- -1
c <- 1
h <- 10
p_crit = d * log(1000.0) ^ (1.0/c)
rs <- 0.1
resp <- rs * rpmodel::calc_ftemp_inst_rd(tc)

ecrit <- psi_soil %>% as.list() %>% 
  map(function(x){get_e_crit(psi_soil = x,K = K, d = d, c = c, h = h)}) %>% unlist()

combinations <- expand.grid(resp = resp, ecrit = ecrit) %>% as_tibble()
comb_psi_tc <- expand.grid(tc = tc, psi_soil = psi_soil) %>% as_tibble()
comb_psi_leaf_tc <- expand.grid(tc_leaf = tc, psi_leaf = psi_leaf) %>% as_tibble()

df <- combinations %>% cbind(comb_psi_tc) %>% cbind(comb_psi_leaf_tc) %>% 
  mutate(alpha = resp/ecrit,
         e_0_1 = K*(-0.1-psi_leaf),
         e_1 = K*(-1-psi_leaf),
         e_2 = K*(-2-psi_leaf),
         alpha_prentice_0_1 = resp/e_0_1,
         alpha_prentice_1 = resp/e_1,
         alpha_prentice_2 = resp/e_2)

df %>% 
  ggplot(aes(psi_soil,tc, fill = log(alpha)))+
  geom_raster()+
  viridis::scale_fill_viridis()

df %>% filter(alpha_prentice_0_1 >=0) %>% 
  ggplot(aes(psi_leaf,tc, fill = alpha_prentice_0_1))+
  geom_raster()+
  xlim(-5,-0.1)+
  viridis::scale_fill_viridis()

df %>% filter(alpha_prentice_1 >=0) %>% 
  ggplot(aes(psi_leaf,tc, fill = alpha_prentice_1))+
  geom_raster()+
  xlim(-5,-1)+
  viridis::scale_fill_viridis()

df %>% filter(alpha_prentice_2 >=0) %>% 
  ggplot(aes(psi_leaf,tc, fill = alpha_prentice_2))+
  geom_raster()+
  xlim(-5,-2)+
  viridis::scale_fill_viridis()
