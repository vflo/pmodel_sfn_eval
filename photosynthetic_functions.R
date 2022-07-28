#Some functions are from yujie Wang code

###########################################
## PHOTOSYNTHETIC FUNCTIONS
###########################################

# get_a_ci <- function(v25,j25,gamma,gc,ca,tem,par){
#   tar_p  = 0.0
#   tar_a  = 0.0
#   max_p  = ca
#   min_p  = gamma
#   adjust = 0.95
#   r_day  = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
#   # r_day = 0
#   continue <- TRUE
#   while(continue){
#     tar_p = 0.5 * (max_p+min_p)
#     vmax = GetPhotosyntheticVcmax(v25,tem)
#     jmax = GetPhotosyntheticJmax(j25,tem)
#     j = GetPhotosyntheticJ(jmax,par)
#     kc = 41.01637 * 2.1**(0.1*(tem-25.0))
#     ko = 28201.92 * 1.2**(0.1*(tem-25.0))
#     km = kc * (1.0+21000.0/ko)
#     aj = j * (tar_p-gamma) / (4.0*(tar_p+2*gamma))
#     ac = vmax * (tar_p-gamma) / (tar_p+km)
#     af = (aj + ac - sqrt((aj+ac)**2.0 - 4*adjust*aj*ac) ) / adjust * 0.5
#     af = af - r_day
#     #af = min(aj,ac) - r_day
#     tmp_g = af / (ca-tar_p)
#     if(abs(tmp_g-gc)/gc < 1E-12){
#       tar_a = af
#       continue <- FALSE
#     }else if(tmp_g<gc){
#       min_p = tar_p
#     }else{
#       max_p = tar_p
#     }
#     if(abs(max_p-min_p) < 1E-12){
#       tar_a = af
#       continue <- FALSE
#     }
#   }
#   return(list(tar_p, tar_a))
# }
# 
# # get one-point measurement Vcmax
# GetOnePointVcmax <- function(ci, an, tem, gamma=2.5){
#   v_min = 1.0
#   v_max = 200.0
#   v25   = 100.0
#   while(v_max-v_min>1){
#     v25   = 0.5 * (v_max + v_min)
#     r_day = v25 * 0.01 *  2.0**(0.1*(tem-25.0)) / (1.0+exp(1.3*(tem-55.0)))
#     vmax = GetPhotosyntheticVcmax(v25,tem)
#     kc = 41.01637 * 2.1**(0.1*(tem-25.0))
#     ko = 28201.92 * 1.2**(0.1*(tem-25.0))
#     km = kc * (1.0+21000.0/ko)
#     ac = vmax * (ci-gamma) / (ci+km)
#     af = ac - r_day
#     if ((af>an & an>0) | (an>af & an<0)){
#       v_max = v25
#       }else{
#         v_min = v25
#       }
#   }
#   return (v25)
# }  
#   
# # calculate j from light
# GetPhotosyntheticJ <- function(jmax, light){
#   a = 0.9
#   b = -0.3*light - jmax
#   c = 0.3*light*jmax
#   j = ( -b - sqrt(b*b-4*a*c) ) / a * 0.5
#   return(j)
#   }
#   
# calculate jmax from temperature
GetPhotosyntheticJmax <- function(jmax25, tem){
  ha=50300.0
  hd=152044.0
  sv=495.0
  t0=298.15
  r=8.315
  c = 1.0 + exp((sv*t0 -hd)/(r*t0))
  t1 = tem + 273.15
  factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
  jmax = jmax25 * factor
  return(jmax)
  }

# calculate vcmax from temperature
GetPhotosyntheticVcmax <- function(vcmax25, tem){
  ha=73637.0
  hd=149252.0
  sv=486.0
  t0=298.15
  r=8.315
  c = 1.0 + exp((sv*t0 -hd)/(r*t0))
  t1 = tem + 273.15
  factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
  vcmax = vcmax25 * factor
  return(vcmax)
  }


calc_assimilation_limiting = function(vcmax, jmax, gs, par_photosynth){
  # gs = calc_gs(dpsi, psi_soil, par_plant = par_plant, par_env = par_env)
  
  # We need not employ numerical root-finding. calculate chi independently assuming Ac and Aj, and bigger of the two will be the limiting one. Accordingly return Ac or Aj
  Ac = calc_assim_rubisco_limited(gs = gs, vcmax = vcmax, par_photosynth = par_photosynth)
  Aj = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
  
  A = list(ac = Ac$a, aj=Aj$a)
  
  if (Ac$ci > Aj$ci ){
    A$ci = Ac$ci
    A$a  = Ac$a
  } else {
    A$ci = Aj$ci
    A$a  = Aj$a
  }
  
  A
}


calc_assim_rubisco_limited <- function(gs, vcmax, par_photosynth){
  
  ca = par_photosynth$ca             # ca is in Pa
  gs = gs * 1e6/par_photosynth$patm  # convert to umol/m2/s/Pa
  
  d = par_photosynth$delta
  ## with
  A <- -1.0 * gs
  B <- gs * ca - gs * par_photosynth$kmm - vcmax*(1-d)
  C <- gs * ca * par_photosynth$kmm + vcmax * (par_photosynth$gammastar + par_photosynth$kmm*d)
  
  ci <- QUADM(A, B, C)
  # a_c <- vcmax * (ci - par$gammastar) / (ci + par$kmm)
  a_c <- gs*(ca-ci) 
  
  return(list(a=a_c, ci=ci))
}


calc_assim_light_limited <- function(gs, jmax, par_photosynth){
  
  ## Only light is limiting
  ## Solve Eq. system
  ## A = gs (ca- ci)
  ## A = phi0 * Iabs * jlim * (ci - gammastar)/(ci + 2*gamma_star)
  
  ## This leads to a quadratic equation:
  ## A * ci^2 + B * ci + C  = 0
  ## 0 = a + b*x + c*x^2
  
  ca = par_photosynth$ca             # ca is in Pa
  gs = gs * 1e6/par_photosynth$patm  # convert to umol/m2/s/Pa
  
  phi0iabs = par_photosynth$phi0 * par_photosynth$Iabs
  jlim = phi0iabs / sqrt(1+ (4*phi0iabs/jmax)^2)
  
  d = par_photosynth$delta 
  ## with
  A <- -1.0 * gs
  B <- gs * ca - gs * 2 * par_photosynth$gammastar - jlim*(1-d)
  C <- gs * ca * 2*par_photosynth$gammastar + jlim * (par_photosynth$gammastar + d*par_photosynth$kmm)
  
  ci <- QUADM(A, B, C)
  aj <- gs*(ca-ci)
  # vcmax_pot <- a*(ci + par$kmm)/(ci - par$gammastar)
  
  return(list(a=aj, ci=ci))
}


calc_vcmax_coordinated_numerical <-  function(aj, ci, par_photosynth){
  d = par_photosynth$delta
  vcmax_coord = aj*(ci + par_photosynth$kmm)/(ci*(1-d) - (par_photosynth$gammastar+par_photosynth$kmm*d))
  return(vcmax_coord)
}
