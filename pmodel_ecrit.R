rpmodel_ecrit <- function (tc, vpd, co2, fapar, ppfd, patm = NA, elv = NA, 
                     kphio = ifelse(do_ftemp_kphio, ifelse(do_soilmstress, 0.087182, 0.081785), 0.049977), beta = 146, b_cost = 0.03,
                     psi_soil, K, d, c, h, rs, LAI = LAI, c_cost = 0.41, soilm = 1, meanalpha = 1, apar_soilm = 0, 
                     bpar_soilm = 0.733, c4 = FALSE, method_optci = "prentice14", 
                     method_jmaxlim = "wang17", do_ftemp_kphio = TRUE, do_soilmstress = FALSE, 
                     returnvar = NULL, verbose = FALSE){
  if (identical(NA, elv) && identical(NA, patm)) {
    rlang::abort("Aborted. Provide either elevation (arugment elv) or atmospheric pressure (argument patm).")
  }
  else if (!identical(NA, elv) && identical(NA, patm)) {
    if (verbose) 
      rlang::warn("Atmospheric pressure (patm) not provided. Calculating it as a function of elevation (elv), assuming standard atmosphere (101325 Pa at sea level).")
    patm <- calc_patm(elv)
  }
  c_molmass <- 12.0107
  kPo <- 101325
  kTo <- 25
  rd_to_vcmax <- 0.015
  if (do_ftemp_kphio) {
    ftemp_kphio <- calc_ftemp_kphio(tc)
  }
  else {
    ftemp_kphio <- 1
  }
  if (do_soilmstress) {
    soilmstress <- calc_soilmstress(soilm, meanalpha, apar_soilm, 
                                    bpar_soilm)
  }
  else {
    soilmstress <- 1
  }
  ca <- co2_to_ca(co2, patm)
  gammastar <- rpmodel::calc_gammastar(tc, patm)
  kmm <- rpmodel::calc_kmm(tc, patm)
  ecrit <- get_e_crit(psi_soil, K, d, c, h)*LAI*1e-3 #micromols??
  resp <- rs * rpmodel::calc_ftemp_inst_rd(tc)
  a_cost = resp/ecrit

  if (c4) {
    out_optchi <- calc_chi_c4()
  } else if (method_optci=="prentice14"){
    
    ## Full formualation (Gamma-star not zero), analytical solution
    ##-----------------------------------------------------------------------
    out_optchi <- calc_optimal_chi_ecrit( kmm, gammastar, a_cost, ca, vpd, beta )
    
  } else if (method_optci=="prentice14_num"){
    out_optim_num <- calc_optim_num_ecrit(kmm = kmm, gammastar = gammastar, ca = ca, vpd = vpd, ppfd = ppfd, fapar = fapar,
                                          kphio = kphio, beta = beta, b_cost = b_cost, c_cost = c_cost, ecrit = ecrit, resp = resp,
                                          vcmax_start = 30, gs_start = 0.8, jmax_start = 40)
    gamma <- gammastar/ca
    kappa <- kmm/ca
    out_optchi <- list(chi = out_optim_num$chi, 
                       mj = (out_optim_num$chi - gamma)/(out_optim_num$chi + 2 * gamma), 
                       mc = (out_optim_num$chi - gamma)/(out_optim_num$chi + kappa), 
                       mjoc = (out_optim_num$chi + kappa)/(out_optim_num$chi + 2 * gamma))
    method_jmaxlim <- "prentice14_num"
  }
  else {
    rlang::abort("rpmodel(): argument method_optci not idetified.")
  }
  ci <- out_optchi$chi * ca
  iwue = (ca - ci)/1.6
  if (c4) {
    out_lue_vcmax <- calc_lue_vcmax_c4(kphio, ftemp_kphio, c_molmass, soilmstress)
  } else {
    out_lue_vcmax <- calc_lue_vcmax_wang17(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress, c_cost)
  }
  
  ftemp25_inst_vcmax <- calc_ftemp_inst_vcmax(tc, tc, tcref = 25)
  vcmax25_unitiabs <- out_lue_vcmax$vcmax_unitiabs/ftemp25_inst_vcmax
  ftemp_inst_rd <- calc_ftemp_inst_rd(tc)
  rd_unitiabs <- rd_to_vcmax * (ftemp_inst_rd/ftemp25_inst_vcmax) * out_lue_vcmax$vcmax_unitiabs
  len <- length(out_lue_vcmax[[1]])
  iabs <- rep(fapar * ppfd, len)
  gpp <- ifelse(!is.na(iabs), iabs * out_lue_vcmax$lue, rep(NA,len))
  vcmax <- ifelse(!is.na(iabs), iabs * out_lue_vcmax$vcmax_unitiabs, rep(NA, len))
  vcmax25 <- ifelse(!is.na(iabs), iabs * vcmax25_unitiabs, rep(NA, len))
  rd <- ifelse(!is.na(iabs), iabs * rd_unitiabs, rep(NA, len))
  fact_jmaxlim <- ifelse(!is.na(iabs), vcmax * (ci + 2 * gammastar)/(kphio * iabs * (ci + kmm)), rep(NA, len))
  jmax <- ifelse(!is.na(iabs), 4 * kphio * iabs/sqrt((1/fact_jmaxlim)^2 - 1), rep(NA, len))
  gs <- (gpp/c_molmass)/(ca - ci)
  eactual <- gs * vpd * 1.6
  a_cost <- resp/ecrit
  out <- list(ca = rep(ca, len), gammastar = rep(gammastar, len), kmm = rep(kmm, len), a_cost = rep(a_cost, len), 
              chi = out_optchi$chi, mj = out_optchi$mj, mc = out_optchi$mc, ci = ci, lue = out_lue_vcmax$lue, gpp = gpp, 
              iwue = iwue, E = eactual, ecrit = ecrit, resp = resp, gs = gs, vcmax = vcmax, vcmax25 = vcmax25, 
              jmax = jmax, rd = rd)
  if (!is.null(returnvar)) 
    out <- out[returnvar]
  return(out)
}

calc_optimal_chi_ecrit <- function (kmm, gammastar, a_cost, ca, vpd, beta){
  vpd <- ifelse(vpd < 0, 0, vpd)
  xi <- sqrt((beta * (kmm + gammastar))/(1.6 * a_cost))
  chi <- gammastar/ca + (1 - gammastar/ca) * xi/(xi + sqrt(vpd))
  vdcg <- ca - gammastar
  vacg <- ca + 2 * gammastar
  vbkg <- beta * (kmm + gammastar)
  calc_mj <- function(a_cost, vpd, vbkg) {
    vsr <- sqrt(1.6 * a_cost * vpd/vbkg)
    mj <- vdcg/(vacg + 3 * gammastar * vsr)
    return(mj)
  }
  mj <- ifelse(a_cost > 0 & vpd > 0 & vbkg > 0, calc_mj(a_cost,vpd, vbkg), rep(NA, length(vpd)))
  gamma <- gammastar/ca
  kappa <- kmm/ca
  mc <- (chi - gamma)/(chi + kappa)
  mjoc <- (chi + kappa)/(chi + 2 * gamma)
  out <- list(chi = chi, mc = mc, mj = mj, mjoc = mjoc)
  return(out)
}

calc_optim_num_ecrit <- function (kmm, gammastar, ecrit, ca, vpd, resp, ppfd, fapar, kphio, beta, b_cost, c_cost, vcmax_start, gs_start, jmax_start){
  
  optimise_this_gs_vcmax_jmax <- function(par, args, iabs, kphio, beta, b_cost, c_cost, maximize = FALSE, return_all = FALSE) {
    kmm <- args[1]
    gammastar <- args[2]
    ecrit <- args[3]
    ca <- args[4]
    vpd <- args[5]
    resp <- args[6]
    vcmax <- par[1]
    gs <- par[2]
    jmax <- par[3]
    L <- 1/sqrt(1 + ((4 * kphio * iabs)/jmax)^2)
    A <- -gs
    B <- gs * ca - 2 * gammastar * gs - L * kphio * iabs
    C <- 2 * gammastar * gs * ca + L * kphio * iabs * gammastar
    ci_j <- QUADM(A, B, C)
    a_j <- kphio * iabs * (ci_j - gammastar)/(ci_j + 2 * gammastar) * L
    A <- -1 * gs
    B <- gs * ca - gs * kmm - vcmax
    C <- gs * ca * kmm + vcmax * gammastar
    ci_c <- QUADM(A, B, C)
    a_c <- vcmax * (ci_c - gammastar)/(ci_c + kmm)
    assim <- min(a_j, a_c)
    ci <- max(ci_c, ci_j)
    eactual <- gs * vpd * 1.6
    a_cost <- resp/ecrit
    cost_transp <- 1.6 * a_cost * gs * vpd
    cost_vcmax <- b_cost * vcmax
    cost_jmax <- c_cost * jmax
    if (assim <= 0) {
      net_assim <- (999999999.9)
    }
    else if(a_cost <= 0) {
      net_assim <- (999999999.9)
    }
    else {
      net_assim <- (cost_transp + cost_vcmax + cost_jmax)/assim
    }
    if (maximize) 
      net_assim <- -net_assim
    if (return_all) {
      return(list(vcmax = vcmax, jmax = jmax, gs = gs, 
                  ci = ci, chi = ci/ca, a_c = a_c, a_j = a_j, 
                  assim = assim, ci_c = ci_c, ci_j = ci_j, cost_transp = cost_transp, 
                  cost_vcmax = cost_vcmax, cost_jmax = cost_jmax, 
                  net_assim = net_assim))
    }
    else {
      return(net_assim)
    }
  }
  
  out_optim <- optimr::optimr(par = c(vcmax_start, gs_start, jmax_start), 
                              lower = c(vcmax_start * 0.001, gs_start * 0.001, jmax_start * 0.001), 
                              upper = c(vcmax_start * 1000, gs_start * 1000, jmax_start * 1000), 
                              fn = optimise_this_gs_vcmax_jmax, 
                              args = c(kmm, gammastar, ecrit, ca, vpd, resp), 
                              iabs = (ppfd * fapar), 
                              kphio = kphio, 
                              beta = beta, 
                              b_cost = b_cost,
                              c_cost = c_cost/4, 
                              method = "L-BFGS-B", 
                              maximize = FALSE,
                              control = list(maxit = 1e+05))
  
  varlist <- optimise_this_gs_vcmax_jmax(par = out_optim$par, 
                                         args = c(kmm, gammastar, ecrit, ca, vpd, resp), 
                                         iabs = (fapar * ppfd), 
                                         kphio, beta,b_cost, c_cost/4, 
                                         maximize = FALSE, 
                                         return_all = TRUE)
  return(varlist)
}



co2_to_ca <- function( co2, patm ){
  #-----------------------------------------------------------------------
  # Input:    - float, annual atm. CO2, ppm (co2)
  #           - float, monthly atm. pressure, Pa (patm)
  # Output:   - ca in units of Pa
  # Features: Converts ca (ambient CO2) from ppm to Pa.
  #-----------------------------------------------------------------------
  ca   <- ( 1.0e-6 ) * co2 * patm         # Pa, atms. CO2
  return( ca )
}


calc_lue_vcmax_wang17 <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress, c_cost){
  
  ## Include effect of Jmax limitation
  len <- length(out_optchi[[1]])
  mprime <- calc_mprime( out_optchi$mj, c_cost )
  
  out <- list(
    
    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * ftemp_kphio * mprime * c_molmass * soilmstress,
    
    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = kphio * ftemp_kphio * out_optchi$mjoc * mprime / out_optchi$mj * soilmstress,
    
    ## complement for non-smith19 
    omega      = rep(NA, len),
    omega_star = rep(NA, len)
    
  )
  
  return(out)
}


calc_mprime <- function( mc, kc ){
  #-----------------------------------------------------------------------
  # Input:  mc   (unitless): factor determining LUE
  # Output: mpi (unitless): modiefied m accounting for the co-limitation
  #                         hypothesis after Prentice et al. (2014)
  #-----------------------------------------------------------------------
  mpi <- mc^2 - kc^(2.0/3.0) * (mc^(4.0/3.0))
  
  # Check for negatives:
  mpi <- ifelse(mpi>0, sqrt(mpi), NA)
  
  return(mpi)
}