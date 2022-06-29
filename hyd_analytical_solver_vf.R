#### Analytical Solver ####

calc_J = function(gs, x, par_photosynth){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  d = par_photosynth$delta
  4*gs*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k)) 
}

calc_jmax_from_J = function(J, par_photosynth){
  p = par_photosynth$phi0 * par_photosynth$Iabs
  4*p/((4*p/J)^2-1)^(1/2)
}


calc_djmax_dJ = function(J, par_photosynth){
  p = par_photosynth$phi0 * par_photosynth$Iabs
  (4*p)^3/((4*p)^2-J^2)^(3/2)
}

calc_dJ_dchi = function(gs, x, par_photosynth){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  d = par_photosynth$delta
  # gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x^2) + 2*g^2 + g*(2*x - 3) - x^2)/(d*(k + x) + g - x)^2)
  4*gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x^2) - ((x-g)^2+3*g*(1-g)))/(d*(k + x) + g - x)^2)
  #gs*ca*(3*(g-1)*g/(g-x)^2 - 1)
}

calc_dJ_ddpsi = function(gsprime, x, par_photosynth){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  d = par_photosynth$delta
  4*gsprime*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k))
}

calc_x_from_dpsi = function(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost){
  gstar = par_photosynth$gammastar/par_photosynth$patm*1e6
  Km = par_photosynth$kmm/par_photosynth$patm*1e6
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  br = par_photosynth$delta
  y = par_cost$gamma
  
  gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)
  
  x = (-2*ca*dpsi*(gstar + br*Km)*y + 
         ca^2*((3 - 2*br)*gstar + br*Km)*gsprime + 
         -sqrt(2)*sqrt(
           ca^2*dpsi*((-3 + 2*br)*gstar - br*Km)*((-1 + br)*ca + gstar + 
                                                    br*Km)*y*
             (-2*dpsi*y + (ca + 2*gstar)*
                gsprime)))/
    (ca^2*(2*(-1 + br)*dpsi*y + ((3 - 2*br)*gstar + br*Km)*
             gsprime))  
  
  x[x<(gstar + br*Km)/(ca - br*ca)]=(gstar + br*Km)/(ca - br*ca)+1e-12
  x
}

calc_delta_from_dpsi = function(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost){
  gstar = par_photosynth$gammastar/par_photosynth$patm*1e6
  Km = par_photosynth$kmm/par_photosynth$patm*1e6
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  br = par_photosynth$delta
  y = par_cost$gamma
  
  gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)
  
  delt = (-2*dpsi*y + (ca + 2*gstar)*gsprime)
  delt
}


calc_dpsi_bound = function(psi_soil, par_plant, par_env, par_photosynth, par_cost){
  gstar = par_photosynth$gammastar/par_photosynth$patm*1e6
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  LAI = par_env$LAI
  y = par_cost$gamma
  K = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant,do_backtransform = TRUE)*LAI
  K = K/(1.6*par_env$vpd/par_env$patm)
  Pox = P(psi_soil, par_plant$psi50, par_plant$b)
  Ppox = Pprime(psi_soil, par_plant$psi50, par_plant$b)
  Pppox = Pprimeprime(psi_soil, par_plant$psi50, par_plant$b)
  
  f1 = function(dpsi){
    gs = calc_gs(dpsi, psi_soil, par_plant, par_env)
    x = calc_x_from_dpsi(dpsi,psi_soil,par_plant, par_env, par_photosynth, par_cost)
    J=calc_J(gs, x, par_photosynth)-4*par_photosynth$phi0*par_photosynth$Iabs
    J
  }
  
  a = (ca + 2*gstar)*K*Pppox*4/8
  b = -(2*y + (ca + 2*gstar)*K*Ppox)
  c = (ca + 2*gstar)*K*Pox
  del = b^2-4*a*c
  
  approx_O2 = (-b-sqrt(b^2-4*a*c))/2/a
  exact = uniroot(f = function(dpsi){(-2*dpsi*y + (ca + 2*gstar)*
                                        calc_gsprime(dpsi, psi_soil, par_plant, par_env))}, interval = c(0,10))$root
  
  use_bound = exact
  
  # if (del>=0){
  #   approx_O2 = (-b-sqrt(b^2-4*a*c))/2/a
  # 
  #   del_x = (-2*approx_O2*y + (ca + 2*gstar)*calc_gsprime(approx_O2, psi_soil, par_plant, par_env))
  #   ajmax_m_phi0iabs = f1(approx_O2)*f1(approx_O2*0.0001)
  #   cat("signs of ajmax-phi0Iabs = ", ajmax_m_phi0iabs,"\n")
  #   if (del_x < 0 | ajmax_m_phi0iabs < 0){
  #     cat("Warning: Need exact calc: delx = ", del_x, ", ajm_by_phi0iabs = ", ajmax_m_phi0iabs,"\n")
  #     exact = uniroot(f = function(dpsi){(-2*dpsi*y + (ca + 2*gstar)*
  #                                           calc_gsprime(dpsi, psi_soil, par_plant, par_env))}, interval = c(0,10))$root
  #     use_bound = exact
  #   }
  #   else{
  #     use_bound = approx_O2
  #     ## v--- Only for debug sake
  #     exact = uniroot(f = function(dpsi){(-2*dpsi*y + (ca + 2*gstar)*
  #                                            calc_gsprime(dpsi, psi_soil, par_plant, par_env))}, interval = c(0,10))$root
  #   }
  # 
  # }
  # else{
  #   approx_O2 = NA
  #   exact = uniroot(f = function(dpsi){(-2*dpsi*y + (ca + 2*gstar)*
  #                                calc_gsprime(dpsi, psi_soil, par_plant, par_env))}, interval = c(0,10))$root
  #   cat("Warning: Need exact calc: approx_O2 = NA \n")
  #   # approx = (ca*K*Pox + 2*gstar*K*Pox)/(ca*K*Ppox +
  #   #                               2*gstar*K*Ppox + 2*y)
  #   use_bound= exact
  # }
  
  
  # cat(psi_soil, ":", exact, " ", approx_O2, " ", use_bound, "\n")
  Iabs_bound = uniroot(f = f1, interval = c(use_bound*0.001,use_bound*0.99))$root
  
  # dpsi=seq(exact*0.001,exact*0.99, length.out=200)
  # plot(y=sapply(X = dpsi, FUN = f1), x=dpsi, type="l")
  
  list(exact=exact, Iabs_bound=Iabs_bound, approx_O2 = approx_O2)
}

derivatives = function(x, psi_soil, par_photosynth, par_plant, par_env, par_cost){
  X = x[1]
  dpsi = x[2]
  
  gs = calc_gs(dpsi, psi_soil, par_plant, par_env)#* 1e6/par_photosynth$patm
  gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)#* 1e6/par_photosynth$patm
  J = calc_J(gs, X, par_photosynth)
  
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  g = par_photosynth$gammastar/par_photosynth$ca
  
  dP_dx = -gs*ca - par_cost$alpha * calc_djmax_dJ(ajmax, par_photosynth) * calc_dJ_dchi(gs, X, par_photosynth)
  
  dP_ddpsi = gsprime*ca*(1-X) - par_cost$alpha * calc_djmax_dJ(ajmax, par_photosynth) * calc_dJ_ddpsi(gsprime, X, par_photosynth) - 2*par_cost$gamma*dpsi #/par_plant$psi50^2
  # cat(c(dP_dx, dP_ddpsi), "\n")
  c(dP_dx, dP_ddpsi)
}


dFdx = function(dpsi, psi_soil, par_photosynth, par_plant, par_env, par_cost){
  gs = calc_gs(dpsi, psi_soil, par_plant, par_env)#* 1e6/par_photosynth$patm
  gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)#* 1e6/par_photosynth$patm
  
  X =  calc_x_from_dpsi(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost)
  
  J = calc_J(gs, X, par_photosynth)
  
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  g = par_photosynth$gammastar/par_photosynth$ca
  
  djmax_dJ = calc_djmax_dJ(J, par_photosynth)
  dJ_dchi = calc_dJ_dchi(gs, X, par_photosynth)
  
  dP_dx = -gs*ca - par_cost$alpha * djmax_dJ * dJ_dchi
  
  # dP_ddpsi = gsprime*ca*(1-X) - par_cost$alpha * calc_djmax_dAjmax(ajmax, par_photosynth) * calc_dAjmax_ddpsi(gsprime, X, par_photosynth) - 2*par_cost$gamma*dpsi #/par_plant$psi50^2
  # # cat(c(dP_dx, dP_ddpsi), "\n")
  # c(dP_dx, dP_ddpsi)
  list(dP_dx=dP_dx, J=J, djmax_dJ=djmax_dJ, dJ_dchi=dJ_dchi)
}

#' Electron-transport limited \eqn{\chi} 
#'
#' Analytically calculate \eqn{\chi} in the case of strong Jmax limitation
#'
#' @param par_photosynth numeric, A list of photosynthesis parameters, which must include ...... . All concentrations must be converted to partial pressures (Pa)
#' @param par_cost numeric, A list of cost parameters
#'
#' @return \eqn{\chi}
#'
#' @export
#'
# Analytical chi in the case of strong Jmax limitation
chi_jmax_limited = function(par_photosynth, par_cost){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  b = par_photosynth$delta
  a = par_cost$alpha
  
  #(g*(1-4*a) + 2*sqrt(3)*sqrt(a*(1-4*a)*(1-g)*g))/(1-4*a)
  # (2*sqrt(-a*(4*a + b - 1)*(2*b^2*g*k + 2*b^2*g + b^2*(-k^2) - b^2*k + 2*b*g^2 - 4*b*g*k - 5*b*g + b*k - 3*g^2 + 3*g)) - 4*a*b*k - 4*a*g + b^2*(-k) - b*g + b*k + g)/((b - 1)*(4*a + b - 1))
  (2*sqrt(-a*(4*a + b - 1)*(-3*g + 2*b*g - b*k)*(-1 + b + g + b*k)) - (4*a + b - 1)*(b*k + g))/((b - 1)*(4*a + b - 1))
}


## Derivative of gs wrt dpsi (mol/m2ground/s/Mpa)
calc_gsprime = function(dpsi, psi_soil, par_plant, par_env){
  K = scale_conductivity_ks(par_plant$conductivity, par_env,par_plant, do_backtransform = TRUE)
  D = (par_env$vpd/par_env$patm)
  K/1.6/D*P(psi_soil-dpsi, par_plant$psi50, par_plant$b)*par_env$LAI
}


#' Hydraulic vulnerability curve
#'
#' Calculates the conductivity of the hydraulic pathway (stem or leaf) as a function 
#' of water potential, using a 2-parameter Weibull function to describe the hydraulic 
#' vulnerability curve:
#' \deqn{P(\psi) = (1/2)^((\psi/\psi_50)^b) }
#'
#' @param psi numeric, Water potential (Mpa)
#' @param psi50 numeric, Water potential at which 50\% conductivity is lost (Mpa)
#' @param b numeric, Slope of the vulnerability curve
#'
#' @return Fractional conductivity
#'
#' @export
#' 
P <- function(psi, psi50, b){
  0.5^((psi/psi50)^b)
}

#' @describeIn P First derivative of the vulnerability curve
Pprime = function(psi, psi50, b){
  log(0.5)*P(psi, psi50, b)*b*(psi/psi50)^(b-1)/psi50
}

#'@describeIn P Second derivative of the vulnerability curve
Pprimeprime = function(psi, psi50, b){
  log(0.5)*b*(psi/psi50)^(b-1)/psi50*Pprime(psi, psi50, b) + log(0.5)*P(psi, psi50, b)/psi50^2*b*(b-1)*(psi/psi50)^(b-2)
}

