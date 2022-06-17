#' Invokes a P-model function call for sub-daily estimations accounting for acclimation
#'
#' R implementation of the P-model and its 
#' corollary predictions (Prentice et al., 2014; Han et al., 2017; Mengoli et al. 2022).
#'
#' @param TIMESTAMP timestamp in a ymd_hms format (as obtained from the package lubridate)
#' @param tc Temperature, relevant for photosynthesis (deg C). Numeric vector.
#' @param vpd Vapour pressure deficit (Pa)
#' @param co2 Atmospheric CO2 concentration (ppm)
#' @param fapar (Optional) Fraction of absorbed photosynthetically active
#'  radiation (unitless, defaults to \code{NA})
#' @param ppfd Incident photosynthetic photon flux density 
#'  (mol m-2 d-1, defaults to \code{NA}). Note that the units of 
#'  \code{ppfd} (per area and per time) determine the units of outputs 
#'  \code{lue}, \code{gpp}, \code{vcmax}, and \code{rd}. For example, 
#'  if \code{ppfd} is provided in units of mol m-2 month-1, then
#'  respective output variables are returned as per unit months.
#' @param sw_in (optional) needed when upscaling_method is "max_rad". 
#' Numeric vector that should be provided in W m-2.
#' @param patm Atmospheric pressure (Pa). When provided, overrides
#'  \code{elv}, otherwise \code{patm} is calculated using standard
#'  atmosphere (101325 Pa), corrected for elevation (argument \code{elv}),
#'  using the function \link{calc_patm}.
#' @param elv Elevation above sea-level (m.a.s.l.). Is used only for 
#'  calculating atmospheric pressure (using standard atmosphere (101325 Pa),
#'  corrected for elevation (argument \code{elv}), using the function
#' \link{calc_patm}), if argument \code{patm} is not provided. If argument
#' \code{patm} is provided, \code{elv} is overridden.
#' @param kphio Apparent quantum yield efficiency (unitless). Defaults to
#'  0.081785 for \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, 
#'  do_soilmstress=FALSE}, 0.087182 for \code{method_jmaxlim="wang17",
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}, and 0.049977 for 
#'  \code{method_jmaxlim="wang17", do_ftemp_kphio=FALSE, do_soilmstress=FALSE},
#'  corresponding to the empirically fitted value as presented in Stocker et al.
#'  (2019) Geosci. Model Dev. for model setup 'BRC', 'FULL', and 'ORG' 
#'  respectively, corresponding to \eqn{(a_L b_L)/4} in 
#'  Eq.20 in Stocker et al. (2020) for C3 photosynthesis. For C4 photosynthesis
#'  (\code{c4 = TRUE}), \code{kphio} defaults to 1.0, corresponding to the 
#'   parametrisation by  Cai & Prentice (2020).
#' @param beta Unit cost ratio. Defaults to 146.0 (see Stocker et al., 2019) for
#'   C3 plants and 146/9 for C4 plants.
#' @param soilm (Optional, used only if \code{do_soilmstress==TRUE}) Relative 
#'  soil moisture as a fraction of field capacity (unitless). Defaults to 1.0 
#'  (no soil moisture stress). This information is used to calculate
#'  an empirical soil moisture stress factor (\link{calc_soilmstress}) whereby
#'  the sensitivity is determined by average aridity, defined by the local 
#'  annual mean ratio of actual over potential evapotranspiration, supplied by
#'  argument \code{meanalpha}.
#' @param meanalpha (Optional, used only if \code{do_soilmstress==TRUE}) Local 
#'  annual mean ratio of actual over potential evapotranspiration, measure for 
#'  average aridity. Defaults to 1.0. Only scalar numbers are accepted. If 
#'  a vector is provided, only the first element will be used.
#' @param apar_soilm (Optional, used only if \code{do_soilmstress==TRUE}) 
#'  Parameter determining the sensitivity of the empirical soil moisture stress 
#'  function. Defaults to 0.0, the empirically fitted value as presented in 
#'  Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' 
#'  (corresponding to a setup with \code{method_jmaxlim="wang17", 
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#' @param bpar_soilm (Optional, used only if \code{do_soilmstress==TRUE}) 
#'  Parameter determining the sensitivity of the empirical soil moisture stress
#'  function. Defaults to 0.7330, the empirically fitted value as presented in 
#'  Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' 
#'  (corresponding to a setup with \code{method_jmaxlim="wang17", 
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#' @param c4 (Optional) A logical value specifying whether the C3 or C4
#'  photosynthetic pathway is followed.Defaults to \code{FALSE}. If \code{TRUE},
#'  the leaf-internal CO2 concentration is still estimated using beta but
#'  \eqn{m} (returned variable \code{mj}) tends to 1, and \eqn{m'} tends to
#'  0.669 (with \code{c = 0.41}) to represent CO2 concentrations within the leaf.
#'  With \code{do_ftemp_kphio = TRUE}, a C4-specific temperature dependence of
#'  the quantum yield efficiency is used (see \link{ftemp_kphio}).
#' @param method_jmaxlim (Optional) A character string specifying which method 
#'  is to be used for factoring in Jmax limitation. Defaults to \code{"wang17"},
#'  based on Wang Han et al. 2017 Nature Plants and (Smith 1937). Available is 
#'  also \code{"smith19"}, following the method by Smith et al., 2019 Ecology 
#'  Letters, and \code{"none"} for ignoring effects of Jmax limitation.
#' @param do_ftemp_kphio (Optional) A logical specifying whether 
#'  temperature-dependence of quantum yield efficiency is used. See \link{ftemp_kphio}
#'  for details. Defaults to \code{TRUE}. Only scalar numbers are accepted. If 
#'  a vector is provided, only the first element will be used.
#' @param do_soilmstress (Optional) A logical specifying whether an empirical 
#' soil moisture stress factor is to be applied to down-scale light use 
#' efficiency (and only light use efficiency). Defaults to \code{FALSE}.
#' @param returnvar (Optional) A character string of vector of character strings
#'  specifying which variables are to be returned (see return below).
#' @param verbose Logical, defines whether verbose messages are printed. 
#'  Defaults to \code{FALSE}.
#' @param upscaling_method String defining the method applied in the sub-daily down scaling. 
#'  Accepted values are "noon", "daily" and "max_rad". "noon" approach computes the 
#'  down scaling using the average of conditions around midday (hour_reference_t +/- nr_window).
#'  "daily" approach uses average daytime conditions. "max_rad" uses average of daytime 
#'  conditions around the point of maximum radiation.
#' @param hour_reference_T numeric from 0 to 23. Reference time for the upscaling process.
#'  
#'
#' @return A named list of numeric values (including temperature and pressure 
#' dependent parameters of the photosynthesis model, P-model predictions, 
#' including all its corollary). This includes :
#' 
#' \itemize{
#'  \item \code{ca}: Ambient CO2 expressed as partial pressure (Pa)
#'  
#'  \item \code{gammastar}: Photorespiratory compensation point \eqn{\Gamma*}, 
#'   (Pa), see \link{calc_gammastar}.
#'  
#'  \item \code{kmm}: Michaelis-Menten coefficient \eqn{K} for photosynthesis 
#'  (Pa), see \link{calc_kmm}.
#'  
#'  \item \code{ns_star}: Change in the viscosity of water, relative to its 
#'   value at 25 deg C (unitless).
#'   \deqn{\eta* = \eta(T) / \eta(25 deg C)}
#'   This is used to scale the unit cost of transpiration. 
#'   Calculated following Huber et al. (2009).
#'  
#'  \item \code{chi}: Optimal ratio of leaf internal to ambient CO2 (unitless). 
#'   Derived following Prentice et al.(2014) as:
#'  \deqn{
#'   \chi = \Gamma* / ca + (1- \Gamma* / ca) \xi / (\xi + \sqrt D )
#'   }
#'   with
#'   \deqn{
#'    \xi = \sqrt (\beta (K+ \Gamma*) / (1.6 \eta*))
#'   }
#'   \eqn{\beta} is given by argument \code{beta}, \eqn{K} is 
#'   \code{kmm} (see \link{calc_kmm}), \eqn{\Gamma*} is 
#'   \code{gammastar} (see \link{calc_gammastar}). \eqn{\eta*} is \code{ns_star}.
#'   \eqn{D} is the vapour pressure deficit (argument \code{vpd}), \eqn{ca} is 
#'   the ambient CO2 partial pressure in Pa (\code{ca}).
#'   
#'   \item \code{ci}: Leaf-internal CO2 partial pressure (Pa), calculated as \eqn{(\chi ca)}.
#'   
#'   \item \code{lue}: Light use efficiency (g C / mol photons), calculated as
#'                         \deqn{
#'                              LUE = \phi(T) \phi0 m' Mc
#'                         }
#'                         where \eqn{\phi(T)} is the temperature-dependent quantum yield efficiency modifier
#'                         (\link{ftemp_kphio}) if \code{do_ftemp_kphio==TRUE}, and 1 otherwise. \eqn{\phi 0}
#'                         is given by argument \code{kphio}.
#'                         \eqn{m'=m} if \code{method_jmaxlim=="none"}, otherwise
#'                         \deqn{
#'                                m' = m \sqrt( 1 - (c/m)^(2/3) )
#'                         }
#'                         with \eqn{c=0.41} (Wang et al., 2017) if \code{method_jmaxlim=="wang17"}. \eqn{Mc} is
#'                         the molecular mass of C (12.0107 g mol-1). \eqn{m} is given returned variable \code{mj}.
#'                         If \code{do_soilmstress==TRUE}, \eqn{LUE} is multiplied with a soil moisture stress factor,
#'                         calculated with \link{calc_soilmstress}.
#'         \item \code{mj}: Factor in the light-limited assimilation rate function, given by
#'                         \deqn{
#'                             m = (ci - \Gamma*) / (ci + 2 \Gamma*)
#'                        }
#'                        where \eqn{\Gamma*} is given by \code{calc_gammastar}.
#'         \item \code{mc}: Factor in the Rubisco-limited assimilation rate function, given by
#'                         \deqn{
#'                             mc = (ci - \Gamma*) / (ci + K)
#'                        }
#'                        where \eqn{K} is given by \code{calc_kmm}.
#'         \item \code{gpp}: Gross primary production (g C m-2), calculated as
#'                        \deqn{
#'                            GPP = Iabs LUE
#'                        }
#'                        where \eqn{Iabs} is given by \code{fapar*ppfd} (arguments), and is
#'                        \code{NA} if \code{fapar==NA} or \code{ppfd==NA}. Note that \code{gpp} scales with
#'                        absorbed light. Thus, its units depend on the units in which \code{ppfd} is given.
#'         \item \code{iwue}: Intrinsic water use efficiency (iWUE, Pa), calculated as
#'                        \deqn{
#'                              iWUE = ca (1-\chi)/(1.6)
#'                        }
#'         \item \code{gs}: Stomatal conductance (gs, in mol C m-2 Pa-1), calculated as
#'                        \deqn{
#'                             gs = A / (ca (1-\chi))
#'                        }
#'                        where \eqn{A} is \code{gpp}\eqn{/Mc}.
#'         \item \code{vcmax}: Maximum carboxylation capacity \eqn{Vcmax} (mol C m-2) at growth temperature (argument
#'                       \code{tc}), calculated as
#'                       \deqn{
#'                            Vcmax = \phi(T) \phi0 Iabs n
#'                       }
#'                       where \eqn{n} is given by \eqn{n=m'/mc}.
#'         \item \code{vcmax25}: Maximum carboxylation capacity \eqn{Vcmax} (mol C m-2) normalised to 25 deg C
#'                      following a modified Arrhenius equation, calculated as \eqn{Vcmax25 = Vcmax / fv},
#'                      where \eqn{fv} is the instantaneous temperature response by Vcmax and is implemented
#'                      by function \link{ftemp_inst_vcmax}.
#'         \item \code{jmax}: The maximum rate of RuBP regeneration () at growth temperature (argument
#'                       \code{tc}), calculated using
#'                       \deqn{
#'                            A_J = A_C
#'                       }
#'         \item \code{rd}: Dark respiration \eqn{Rd} (mol C m-2), calculated as
#'                      \deqn{
#'                          Rd = b0 Vcmax (fr / fv)
#'                      }
#'                      where \eqn{b0} is a constant and set to 0.015 (Atkin et al., 2015), \eqn{fv} is the
#'                      instantaneous temperature response by Vcmax and is implemented by function
#'                      \link{ftemp_inst_vcmax}, and \eqn{fr} is the instantaneous temperature response
#'                      of dark respiration following Heskel et al. (2016) and is implemented by function
#'                      \link{ftemp_inst_rd}.
#' }
#'
#' Additional variables are contained in the returned list if argument \code{method_jmaxlim=="smith19"}
#' \itemize{
#'  \item \code{omega}: Term corresponding to \eqn{\omega}, defined by Eq. 16 in
#'   Smith et al. (2019), and Eq. E19 in Stocker et al. (2019).
#'   
#'  \item \code{omega_star}: Term corresponding to \eqn{\omega^\ast}, defined by
#'   Eq. 18 in Smith et al. (2019), and Eq. E21 in Stocker et al. (2019).
#'  }patm
#'
#' @references  
#'  Bernacchi, C. J., Pimentel, C., and Long, S. P.:  In vivo temperature response func-tions  of  parameters
#'  required  to  model  RuBP-limited  photosynthesis,  Plant  Cell Environ., 26, 1419–1430, 2003
#'
#   Cai, W., and Prentice, I. C.: Recent trends in gross primary production 
#'  and their drivers: analysis and modelling at flux-site and global scales,
#'  Environ. Res. Lett. 15 124050 https://doi.org/10.1088/1748-9326/abc64e, 2020
#
#'  Heskel,  M.,  O’Sullivan,  O.,  Reich,  P.,  Tjoelker,  M.,  Weerasinghe,  L.,  Penillard,  A.,Egerton, J.,
#'  Creek, D., Bloomfield, K., Xiang, J., Sinca, F., Stangl, Z., Martinez-De La Torre, A., Griffin, K.,
#'  Huntingford, C., Hurry, V., Meir, P., Turnbull, M.,and Atkin, O.:  Convergence in the temperature response
#'  of leaf respiration across biomes and plant functional types, Proceedings of the National Academy of Sciences,
#'  113,  3832–3837,  doi:10.1073/pnas.1520282113,2016.
#'
#'  Huber,  M.  L.,  Perkins,  R.  A.,  Laesecke,  A.,  Friend,  D.  G.,  Sengers,  J.  V.,  Assael,M. J.,
#'  Metaxa, I. N., Vogel, E., Mares, R., and Miyagawa, K.:  New international formulation for the viscosity
#'  of H2O, Journal of Physical and Chemical ReferenceData, 38, 101–125, 2009
#'
#'  Mengoli, G., Agustí-Panareda, A., Boussetta, S., Harrison, S. P., Trotta, C., 
#'  & Prentice, I. C. (2022). Ecosystem photosynthesis in land-surface models: A first-principles approach 
#'  incorporating acclimation. Journal of Advances in Modeling Earth Systems, 14,
#'   e2021MS002767. https://doi.org/10.1029/2021MS002767 
#'  
#'  Prentice,  I. C.,  Dong,  N.,  Gleason,  S. M.,  Maire,  V.,  and Wright,  I. J.:  Balancing the costs
#'  of carbon gain and water transport:  testing a new theoretical frameworkfor  plant  functional  ecology,
#'  Ecology  Letters,  17,  82–91,  10.1111/ele.12211,http://dx.doi.org/10.1111/ele.12211, 2014.
#'
#'  Wang, H., Prentice, I. C., Keenan, T. F., Davis, T. W., Wright, I. J., Cornwell, W. K.,Evans, B. J.,
#'  and Peng, C.:  Towards a universal model for carbon dioxide uptake by plants, Nat Plants, 3, 734–741, 2017.
#'  Atkin, O. K., et al.:  Global variability in leaf respiration in relation to climate, plant func-tional
#'  types and leaf traits, New Phytologist, 206, 614–636, doi:10.1111/nph.13253,
#'  https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.13253.
#'
#'  Smith, N. G., Keenan, T. F., Colin Prentice, I. , Wang, H. , Wright, I. J., Niinemets, U. , Crous, K. Y.,
#'  Domingues, T. F., Guerrieri, R. , Yoko Ishida, F. , Kattge, J. , Kruger, E. L., Maire, V. , Rogers, A. ,
#'  Serbin, S. P., Tarvainen, L. , Togashi, H. F., Townsend, P. A., Wang, M. , Weerasinghe, L. K. and Zhou, S.
#'  (2019), Global photosynthetic capacity is optimized to the environment. Ecol Lett, 22: 506-517.
#'  doi:10.1111/ele.13210
#'
#'  Stocker, B. et al. Geoscientific Model Development Discussions (in prep.)
#'
#' @export
#'
#' @examples \dontrun{
#'  rpmodel(
#'   tc = 20,
#'   vpd = 1000,
#'   co2 = 400,
#'   ppfd = 30,
#'   elv = 0)
#' }
#'
rpmodel_subdaily <- function(
    TIMESTAMP, tc, vpd, co2, fapar, ppfd, sw_in, patm = NA, elv = NA,
    kphio = ifelse(do_ftemp_kphio, ifelse(do_soilmstress, 0.087182, 0.081785), 0.049977),
    beta = 146.0, c_cost = 0.41, soilm = 1.0, meanalpha = 1.0, apar_soilm = 0.0, bpar_soilm = 0.73300,
    c4 = FALSE, method_optci = "prentice14", method_jmaxlim = "wang17",
    do_ftemp_kphio = TRUE, do_soilmstress = FALSE, returnvar = NULL, verbose = FALSE,
    upscaling_method = c("noon","daily","max_rad"), hour_reference_T = 12
){
  
  # 1.0 Calculate P model without acclimation
  df_Or <- rpmodel::rpmodel(tc, vpd, co2, fapar, ppfd, patm, elv, kphio, beta, c_cost, 
                            soilm, meanalpha, apar_soilm, bpar_soilm, c4, method_optci, 
                            method_jmaxlim, do_ftemp_kphio, do_soilmstress , returnvar, 
                            verbose)
  
  # 2.0 Create df for upscaling
  dfIn <- tibble(YEAR = lubridate::year(TIMESTAMP),
                 MONTH = lubridate::month(TIMESTAMP),
                 DAY = lubridate::day(TIMESTAMP),
                 HOUR = lubridate::hour(TIMESTAMP),
                 MINUTE = lubridate::minute(TIMESTAMP),
                 fapar = fapar,
                 tc = tc,
                 ppfd = ppfd,
                 vpd = vpd,
                 co2 = co2,
                 sw_in = sw_in)
  
  # 2.1 apply the dailyUpscaling function 
  dataDaily <- dailyUpscaling(df = dfIn, 
                                nrWindow = 1, 
                                hour_reference_T = 12, 
                                upscaling_method = "noon")
  
  # 3.0 apply running mean
  dataDaily = runningMean(data_x_running_mean_t = dataDaily, daily_window = 15)
  
  
  # 4.0 Calculate P-model on daily upscaling values
  df_Mm <- rpmodel::rpmodel(dataDaily$tc, dataDaily$vpd, dataDaily$co2, dataDaily$fapar, 
                            dataDaily$ppfd, patm, unique(elv), kphio, beta, c_cost, 
                            soilm, meanalpha, apar_soilm, bpar_soilm, c4, method_optci, 
                            method_jmaxlim, do_ftemp_kphio, do_soilmstress , returnvar, 
                            verbose)
  
  
  # 5.0 Downscale from daily to subdaily
  dfReplicate = dfIn[,c("YEAR","MONTH","DAY","HOUR","MINUTE")]
  
  dfReplicate[,'KOpt'] = NA
  dfReplicate[,'GammaStarOpt'] = NA
  dfReplicate[,'viscosityOpt'] = NA
  
  dfReplicate[,'xiOpt'] = NA
  dfReplicate[,'ciOpt'] =NA
  
  dfReplicate[,'VPDOpt'] = NA
  dfReplicate[,'CO2Opt'] = NA
  dfReplicate[,'caOpt'] = NA
  dfReplicate[,'IabsOpt'] = NA
  dfReplicate[,'TempOpt'] = NA
  
  
  # xiAcclimated option: 'off'/'on'
  dfReplicate[,'vcmaxOpt'] = NA     
  dfReplicate[,'JmaxOpt'] = NA      
  
  dataDaily = dataDaily %>% cbind(df_Mm) 
  
  while (nrow(dataDaily) > 0 ) {
    posValues = which(
      dfReplicate[,'YEAR'] == dataDaily[,'YEAR'][1] &
        dfReplicate[,'MONTH'] == dataDaily[,'MONTH'][1] &
        dfReplicate[,'DAY'] == dataDaily[,'DAY'][1] &
        dfReplicate[,'HOUR'] == dataDaily[,'HOUR'][1] &
        dfReplicate[,'MINUTE'] == dataDaily[,'MINUTE'][1]
    )
    dfReplicate[,'KOpt'][posValues] = dataDaily[,'kmPa'][1]
    dfReplicate[,'GammaStarOpt'][posValues] = dataDaily[,'GammaStarM'][1]
    dfReplicate[,'viscosityOpt'][posValues] = dataDaily[,'viscosityH2oStar'][1]
    dfReplicate[,'xiOpt'][posValues] = dataDaily[,'xiPaM'][1]                    # xiAcclimated: 'off'
    # dfReplicate[,'xiOpt'][posValues] = dataDaily[,'xiPa'][1]
    dfReplicate[,'ciOpt'][posValues] = dataDaily[,'ciM'][1]                      # xiAcclimated: 'off'
    # dfReplicate[,'ciOpt'][posValues] = dataDaily[,'ci'][1]
    dfReplicate[,'VPDOpt'][posValues] = dataDaily[,'VPDPa'][1]
    dfReplicate[,'CO2Opt'][posValues] = dataDaily[,'CO2'][1]
    dfReplicate[,'caOpt'][posValues] = dataDaily[,'ca'][1]
    dfReplicate[,'IabsOpt'][posValues] = dataDaily[,'Iabs'][1]
    # if xiAcclimated is ON vcmaxOpt and JmaxOpt will be overwritten
    dfReplicate[,'vcmaxOpt'][posValues] = dataDaily[,'vcmaxPmodelM1'][1]  # xiAcclimated: 'off'
    dfReplicate[,'JmaxOpt'][posValues] = dataDaily[,'JmaxM1'][1]          # xiAcclimated: 'off'
    
    
    # AVERAGE TEMPERATURE
    dfReplicate[,'TempOpt'][posValues] = dataDaily[,'Ta'][1]
    
    dataDaily = dataDaily[-1,] 
  }
  
  # gapFilling function (set what type of approach 'i.e. constant' in the file setting)----
  for (col in colnames(dfReplicate)) {
    if ( col == 'YEAR') next
    if ( col == 'MONTH') next
    if ( col == 'DAY') next
    if ( col == 'HOUR') next
    if ( col == 'MINUTE') next
    cat(sprintf('%s gap = %d\n',col,sum(is.na(dfReplicate[,col]))))
    dfReplicate[,col] = gapFilling(v = dfReplicate[,col],
                                   vYear = dfReplicate[,'YEAR'],
                                   vMonth = dfReplicate[,'MONTH'],
                                   vDay = dfReplicate[,'DAY'],
                                   vHour = dfReplicate[,'HOUR'],
                                   vMinute = dfReplicate[,'MINUTE'],
                                   approccio = df_settaggio$method_gf_replicate)
    cat(sprintf('%s post gap = %d\n',col,sum(is.na(dfReplicate[,col]))))
  }
  rm(col)
  
  # OPTIMAL Vcmax and Jmax----
  # Vcmax with acclimated 'xiPa', 'ci' and 'phi0' (ON)
  if (df_settaggio$xi_acclimated == 'off') {
    cat(sprintf('vcmaxOpt and JmaxOpt replicated\n'))
    dfReplicate[,'ci'] = NA
    dfReplicate[,'xiPa'] = NA
  } else {
    
    cat(sprintf('vcmaxOpt and JmaxOpt calculated\n'))
    
    # Intrinsic quantum efficiency of photosynthesis (phi0)
    if (is.na(df_settaggio$phi0_method) ) {
      
      phi0 =(1/8) *(0.352+0.022*modelloOr[,'Ta'] -
                      0.00034*modelloOr[,'Ta']^(2) )      # Temperature dependence function of phi0 (Bernacchi et al.,2003)
    } else {
      phi0 = phi0Method
    }
    
    c=0.41                        # cost factor for electron transport capacity
    
    # acclimated xiPa (parameter that determines the sensitivity of ci/ca to VPD)
    betha = 146
    
    dfReplicate[,'xiPa'] = sqrt ( (betha * (dfReplicate[,'KOpt'] +
                                              dfReplicate[,'GammaStarOpt']) ) /
                                    (1.6 * dfReplicate[,'viscosityOpt']) )                      # [Pa^1/2]
    
    # acclimated ci (with acclimated xiPa, and adjusted with the actual VPD)
    dfReplicate[,'ci'] = ( dfReplicate[,'xiPa'] * dfReplicate[,'caOpt'] +
                             dfReplicate[,'GammaStarOpt']*sqrt(modelloOr[,'VPDPa']))/(
                               dfReplicate[,'xiPa']+sqrt(modelloOr[,'VPDPa']))
    
    
    # OPTIMAL Vcmax
    dfReplicate[,'vcmaxOpt']  = phi0 * dfReplicate[,'IabsOpt'] *
      ((dfReplicate[,'ci'] + dfReplicate[,'KOpt']) / (dfReplicate[,'ci'] +
                                                        2*dfReplicate[,'GammaStarOpt'])) *
      sqrt (1-(c * (dfReplicate[,'ci']+2*dfReplicate[,'GammaStarOpt'])/(
        dfReplicate[,'ci']-dfReplicate[,'GammaStarOpt']))^(2/3))                                        # [micromol/m2s]
    
    # OPTIMAL Jmax
    dfReplicate[,'JmaxOpt']  = ( 4 * phi0 * dfReplicate[,'IabsOpt'] ) / sqrt (
      1 / ( 1 - ( c * ( dfReplicate[,'ci'] + 2*dfReplicate[,'GammaStarOpt'] ) /
                    ( dfReplicate[,'ci'] - dfReplicate[,'GammaStarOpt'] ) ) ^ (2.0/3.0) ) - 1 )           #[micromol/m2s]
  }
  
  # Optimal and actual temperature conversion
  dfReplicate[,'topt'] = dfReplicate[,'TempOpt'] + 273.15                               # [K]   
  dfReplicate[,'tk'] = modelloOr[,'Ta'] + 273.15                                        # [K]
  
  # INSTANTANEOUS Vcmax and Jmax----
  # The Arrhenius equation constants:
  Ha = 65330                       # J mol-1 
  Haj = 43900
  Rgas = 8.314                     # J/mol*K
  
  
  dfReplicate[,'vcmaxAdjusted'] = dfReplicate[,'vcmaxOpt'] * exp( (Ha / Rgas)*
                                                                    (1/dfReplicate[,'topt'] - 1/dfReplicate[,'tk']))
  
  dfReplicate[,'JmaxAdjusted'] = dfReplicate[,'JmaxOpt'] * exp( (Haj / Rgas)*
                                                                  (1/dfReplicate[,'topt'] - 1/dfReplicate[,'tk']))
  
  rm(Rgas,Ha,Haj)
  
  
  # CALCULATE the assimilation rate: Ac-----
  if (df_settaggio$xi_acclimated == 'off') {
    Ac <- dfReplicate[,'vcmaxAdjusted']*(modelloOr[,'ciM'] - modelloOr[,'GammaStarM']) /
      (modelloOr[,'ciM'] + modelloOr[,'kmPa'])                                              #[micromol/m2s]
  } else {
    # acclimated Ac with the acclimated xiPa term
    Ac = dfReplicate[,'vcmaxAdjusted']*(
      dfReplicate[,'ci'] - modelloOr[,'GammaStarM']) /
      (dfReplicate[,'ci'] + modelloOr[,'kmPa'])                                             #[micromol/m2s]
  }
  
  
  # CALCULATE the assimilation rate: AJ----
  # Intrinsic quantum efficiency of photosynthesis [mol Co2/mol photons]
  if (is.na(df_settaggio$phi0_method) ) {
    phi0 =(1/8) *(0.352+0.022*modelloOr[,'Ta'] - 0.00034*modelloOr[,'Ta']^(2) )      
  } else {
    phi0 = phi0_method
  }
  
  # electron transport rate (Smith equation)
  dfReplicate[,'J'] =  (4 *phi0*modelloOr[,'Iabs'])/sqrt(1 + (
    (4*phi0*modelloOr[,'Iabs'])/(dfReplicate[,'JmaxAdjusted']))^(2))
  
  if ( df_settaggio$xi_acclimated == 'off' ) {
    AJ = (dfReplicate[,'J']/4)*(modelloOr[,'ciM'] - modelloOr[,'GammaStarM'])/(
      modelloOr[,'ciM'] + 2*modelloOr[,'GammaStarM'])    
  } else {# acclimated AJ with the acclimated xiPa term
    AJ = (dfReplicate[,'J']/4)*(dfReplicate[,'ci'] -
                                  modelloOr[,'GammaStarM'])/(dfReplicate[,'ci'] + 2*modelloOr[,'GammaStarM'])    #[micromol/m2s]
  }
  
  # Rename Ac, AJ and save variables in a df.
  dfReplicate[,'Ac1Opt'] = Ac
  dfReplicate[,'AJp1Opt'] = AJ
  rm(Ac, AJ, phi0)
  
  out = dfReplicate[,c(
    'Ac1Opt','AJp1Opt','caOpt','ci','ciOpt','CO2Opt',
    'GammaStarOpt','IabsOpt','J','JmaxOpt',
    'JmaxAdjusted','KOpt','TempOpt','vcmaxAdjusted',
    'vcmaxOpt','VPDOpt',
    'xiPa')]
  
  # COMPUTE GPPp as minimum value between Ac and AJ ----
  out$GPPpOpt = NA
  
  cat(sprintf('Calculate GPPpOpt, please wait ....\n'))
  for ( cy1 in seq(1,nrow(out)) ) {
    tmp = c(out[,'Ac1Opt'][cy1],out[,'AJp1Opt' ][cy1])
    posNa = which(is.na(tmp) == 1)
    if ( length(posNa) > 0 ) tmp = tmp[-1*posNa]
    rm(posNa)
    if (length(tmp) == 2 )
      out[,'GPPpOpt'][cy1] = min(tmp)
  }
  rm(cy1,tmp)
  
  cat(sprintf('Calculate GPPpOpt .... OK\n'))
  
  # CREATE and SAVE a data.frame with model outputs ----
  
  dataOr = cbind(modelloOr,out)
}




dailyUpscaling <- function(df = dfIn, nrWindow = 1, hour_reference_T = 12, 
                             upscaling_method = "noon") {
  
  #1.0 Header control
  if(is.character(upscaling_method)){
    stop(
      cat('The upscaling_method should be a character string either "noon", "daily" or "max_rad".')
    )
  }else  if(!upscaling_method %in% c("noon","daily","max_rad")){
    stop(
      cat('The upscaling_method provided is not a valid value. It should be either "noon", "daily" or "max_rad"')
      )
  }
  
  colMandatory = c('YEAR','MONTH','DAY','HOUR','MINUTE','NIGHT','TIMESTAMP_START')
  if (upscaling_method == "max_rad")
    colMandatory = c(colMandatory,"sw_in")
  
  if ( headerControl_dd(df = df, colMandatory = colMandatory))
    stop(headerControl_dd(df = df, colMandatory = colMandatory, showMsg = TRUE))
  
  #2.0 T calculation
  dfDayT = df[1,]    
  
  
  for (cicloAnno in sort(unique(df$YEAR))) {
    for (cicloMesi in seq(1,12)) {
      for (cicloGiorni in seq(1,31)) {
        
        posDay = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi &
                         df$DAY == cicloGiorni)
        if (length(posDay) == 0) next
        
        for (hourReference in hour_reference_T) {
          # upscaling_method noon
          if (upscaling_method == "noon") {
            posReference = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni &
                                   df$NIGHT == 0 & df$HOUR == hourReference & df$MINUTE == 0)
            if (length(posReference) == 0) next
            windowDay = seq(posReference - nrWindow,posReference + nrWindow)
          }
          # upscaling_method daily
          if (upscaling_method == "daily") {
            posReference = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni &
                                   df$NIGHT == 0 & df$HOUR == 12 & df$MINUTE == 0)
            windowDay = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni &
                                df$NIGHT == 0)
          }
          # upscaling_method max_rad
          if (upscaling_method == "max_rad") {
            posReferenceMax = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni &
                                      df$NIGHT == 0 & is.na(df$SWINPOT) == 0)
            if ( length(posReferenceMax) > 0 ) {
              
              posMaxSWINPOT = which(df$SWINPOT[posReferenceMax] == max(df$SWINPOT[posReferenceMax],na.rm = T)[1])
              
              posReference = posReferenceMax[posMaxSWINPOT[1]]
              windowDay = seq(posReference - nrWindow,posReference + nrWindow)
            } else {
              posReference = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni &
                                     df$NIGHT == 0 & df$HOUR == 12 & df$MINUTE == 0)
              windowDay = NA
            }
          }
          
          
          dfDay = df[1,]
          for (col in colnames(df)) {
            if ( is.na(sum(windowDay)) ) {
              dfDay[1,col] = NA  
              next
            }
            tmp = df[windowDay,col]
            if (!is.numeric(tmp) ) {
              dfDay[1,col] = NA
              next
            }
            
            posNa = which(is.na(tmp) == 1)
            if (length(posNa) > 0) tmp = tmp[-1*posNa]
            if (length(tmp) == 0) {
              dfDay[1,col] = NA
            } else {
              dfDay[1,col] = mean(tmp)
            }
            rm(tmp,posNa)
          }
          rm(col)
          dfDay$TIMESTAMP_START = df$TIMESTAMP_START[posReference]
          dfDayT = rbind(dfDayT,dfDay)     
          rm(dfDay,windowDay,posReference)
          if (upscaling_method == "daily") break;
          if (upscaling_method == "max_rad") break;
        }
      }
    }
  }
  dfDayT = dfDayT[-1,]
  
  dfDayT$YEAR = year(ymd_hm(as.character(dfDayT$TIMESTAMP_START)))
  dfDayT$MONTH = month(ymd_hm(as.character(dfDayT$TIMESTAMP_START)))
  dfDayT$DAY = day(ymd_hm(as.character(dfDayT$TIMESTAMP_START)))
  dfDayT$HOUR = hour(ymd_hm(as.character(dfDayT$TIMESTAMP_START)))
  dfDayT$MINUTE = minute(ymd_hm(as.character(dfDayT$TIMESTAMP_START)))
  
  dfDayT = dfDayT[order(dfDayT$TIMESTAMP_START),]
  
  return(dfDayT)
}




#' Auxiliar function for dailydounscaling function that check if the mandatory variables 
#' exist in the dataset (missing or redundant columns)
#' @param df dataframe to check the colnames
#' @param colMandatory list of mandatory headers
#' @param showMsg Logical. If TRUE, it shows the function messages. Defaults to \code{FALSE}.
#' @return TRUE or FALSE. 

headerControl_dd <- function(df = dfToCheck, colMandatory = listMandatoryToCheck,showMsg = F) {
  
  # string of mandatory variables, missing
  errorColMissing = c()
  
  # string of mandatory variables, duplicate
  errorColMultiple = c()
  
  for (col in colMandatory){
    ckMandatory = which(colnames(df) == col)
    if (length(ckMandatory) == 0) {
      errorColMissing = c(errorColMissing,col)
      next
    }
    if (length(ckMandatory) > 1) {
      errorColMultiple = c(errorColMultiple,col)
      next
    }
  }
  rm(col)
  if (showMsg) {
    cat(sprintf('missing mandatory variables: %d \n',length(errorColMissing)))   
    if (length(errorColMissing) > 0)
      cat(sprintf(' (%s)\n',paste(errorColMissing,collapse = ',')))
    
    cat(sprintf('multiple mandatory variables: %d \n',length(errorColMultiple)))
    if (length(errorColMultiple) > 0)
      cat(sprintf(' (%s)\n',paste(errorColMultiple,collapse = ',')))
  }
  if ( (length(errorColMissing) + length(errorColMultiple)) > 0) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}



# Function to apply the running mean method on the model inputs 
# param: data_x_running_mean_t, dataframe to use for the runningMean 
# param: daily_window, number of days to compute the running mean 
# param: showMsg (if T shows the function messages)
# return: data.frame with runningMean output
# rdname: runningMean


runningMean <- function(data_x_running_mean_t = df, daily_window = 10) {
  
  unique_hour = sort(unique(data_x_running_mean_t$HOUR))
  
  data_running_mean_t = data_x_running_mean_t[1,]
  
  for (cy_hour_ref in unique_hour) {
    pos_hour = which(data_x_running_mean_t$HOUR == cy_hour_ref)
    
    data_x_running_mean = data_x_running_mean_t[pos_hour,]
    rm(pos_hour)
    # window of time 
    lunghezza_finestra = daily_window - 1# remotion of one day since the computation starts from the 1st
    data_running_mean = data_x_running_mean[1,]
    
    for ( ciclo_inizio_mm in seq(2,nrow(data_x_running_mean)) ) {
      # definition of positions to compute the running mean by making a counter in reverse
      posizioni_rm = (ciclo_inizio_mm - lunghezza_finestra):ciclo_inizio_mm
      
      pos_meno1 = which(posizioni_rm < 1)
      if ( length(pos_meno1) > 0 ) posizioni_rm = posizioni_rm[-1*pos_meno1]
      rm(pos_meno1)
      
      # cycle for computing the mean of the dataset's variables
      media_1 = data_x_running_mean[1,]
      for ( ciclo_variabili in colnames(media_1) ) {
        # dataset creation with values to use
        data_1 = data_x_running_mean[posizioni_rm,ciclo_variabili]
        if ( !is.numeric(data_1) ) {
          media_1[ciclo_variabili] = NA
          next
        }
        # NA remotion
        pos_na = which(is.na(data_1) == 1 )
        if ( length(pos_na) > 0 ) data_1 = data_1[-1*pos_na]
        rm(pos_na)
        # if there are missing values it puts NA, otherwise it computes the mean
        if ( length(data_1) == 0 ) {
          media_1[ciclo_variabili] = NA
        } else {
          media_1[ciclo_variabili] = mean(data_1)
        }
        rm(data_1)
      }
      rm(ciclo_variabili)
      media_1$TIMESTAMP_START = data_x_running_mean$TIMESTAMP_START[ciclo_inizio_mm]
      data_running_mean = rbind(data_running_mean,media_1)
      rm(media_1)
      rm(posizioni_rm)
    }
    rm(ciclo_inizio_mm)
    rm(lunghezza_finestra)
    rm(data_x_running_mean)
    data_running_mean_t = rbind(data_running_mean_t,data_running_mean)
  }
  
  data_running_mean_t = data_running_mean_t[-1,]
  
  library(lubridate)
  data_running_mean_t$YEAR = year(ymd_hm(as.character(data_running_mean_t$TIMESTAMP_START)))
  data_running_mean_t$MONTH = month(ymd_hm(as.character(data_running_mean_t$TIMESTAMP_START)))
  data_running_mean_t$DAY = day(ymd_hm(as.character(data_running_mean_t$TIMESTAMP_START)))
  data_running_mean_t$HOUR = hour(ymd_hm(as.character(data_running_mean_t$TIMESTAMP_START)))
  data_running_mean_t$MINUTE = minute(ymd_hm(as.character(data_running_mean_t$TIMESTAMP_START)))
  
  data_running_mean_t = data_running_mean_t[order(data_running_mean_t$TIMESTAMP_START),]
  
  return(data_running_mean_t)
}



# Function to fill in the missing data of a vector according to the settings
# param: v, data vector to be filled
# param: vYear, vector of the years
# param: vMonth, vector of the months
# param: vDay, vector of the days
# param: vHour, vector of the hours
# param: vMinute, vector of the minutes
# param: approccio (approach) if 'constant' replicates measured data on missing values;
#                             if 'linear' it applies linear regression to fill in missing values among measured data
# param: showMsg (if T shows the function messages)
# return: a vector with required inputs for pmodelPlus
# rdname: gapFilling

gapFilling = function(v = vettore,v_name = NA,
                      vYear = vYear, vMonth = vMonth, vDay = vDay,
                      vHour = vHour, vMinute = vMinute, showMsg = F,
                      approccio = 'constant') {
  
  # check if there are values in v 
  if ( sum(is.na(v)) == length(v) ) {
    warning('variable not gap filled\n') 
    return(v)
  }
  str_msg = ''
  if (!is.na(v_name))
    str_msg = paste0(str_msg,sprintf('GAP FILLING variable: %s\n',v_name))
  
  str_msg = paste0(str_msg,
                   sprintf('approach: %s\nmissing value%s: %s\n',
                           approccio,
                           ifelse(sum(is.na(v)) == 1,'','s'),sum(is.na(v))))
  
  if (showMsg)
    cat(str_msg)
  
  if (approccio == 'constant') {
    # 
    # for (cyr in seq(2,length(v))) {
    # 
    #   if (is.na(v[cyr]))
    #     v[cyr] = v[cyr-1]
    # }
    # rm(cyr)
    posValori = which(is.na(v) == 0)
    for (contaV in posValori) {
      
      yearToUse  = vYear[contaV]
      monthToUse = vMonth[contaV]
      dayToUse   = vDay[contaV]
      vToUse = v[contaV]
      
      posToUse = which(
        vYear == yearToUse &
          vMonth == monthToUse &
          vDay == dayToUse)
      
      if ( length(posToUse) > 0 )
        v[posToUse] = vToUse
      
      rm(yearToUse,monthToUse,dayToUse)
      rm(vToUse,posToUse)
    }
    rm(contaV)
  }
  #approccio: 'linear' 
  if (approccio == 'linear') {
    posValori = which(is.na(v) == 0)
    for (contaV in seq(2,length(posValori))) {
      # slope computetation 
      x1 = posValori[contaV - 1]
      x2 = posValori[contaV]
      y1 = v[x1]
      y2 = v[x2]
      slope = (y2 - y1) / (x2-x1)
      intercept = y1 - (slope*x1)
      itp = (slope * seq(x1,x2)) + intercept
      v[seq(x1,x2)] = itp
    }
  }
  
  if (showMsg)
    cat(sprintf('GAP FILLING COMPLETED\nmissing value%s:%d\n',
                ifelse(sum(is.na(v)) == 1,'','s'),sum(is.na(v))))
  
  return(v)
  
}