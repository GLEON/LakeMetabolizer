#'@title Metabolism calculated from the maximum likelihood estimates of the parameters in a standard linear regression model
#'@description Process-error-only model with parameters fitted via maximum likelihood estimation (MLE). This function runs the maximum likelihood metabolism model on the supplied gas concentration and other supporting data.
#'@param do.obs Vector of dissolved oxygen concentration observations, \eqn{mg O_2 L^{-1}}{mg O2 / L}
#'@param do.sat Vector of dissolved oxygen saturation values based on water temperature. Calculate using \link{o2.at.sat}
#'@param k.gas Vector of kGAS values calculated from any of the gas flux models
#'(e.g., \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}
#'@param z.mix Vector of mixed-layer depths in meters. To calculate, see \link[rLakeAnalyzer]{ts.meta.depths}
#'@param irr Vector of photosynthetically active radiation in \eqn{\mu mol\ m^{-2} s^{-1}}{micro mols / m^2 / s}
#'@param wtr Vector of water temperatures in \eqn{^{\circ}C}{degrees C}. Used in scaling respiration with temperature
#'@param error.type Option specifying if model should assume pure Process Error 'PE' or Observation Error 'OE'. Defaults to observation error 'OE'.
#'@param n.boot Numeric of the number of bootstrap iterations. Set to 0 to bypass uncertainty quantification (default). Set to 2 or greater to quantify uncertainty via bootstrapping. To ensure reproducible confidence intervals, use \link{set.seed}.
#'@param ar1.resids Logical. If TRUE, the AR(1) structure of the residuals will be retained during randomization when bootstrapping.
#'@param constrain.sign Logical. If TRUE, estimates of GPP and R coefficients with be constrained to positive and negative values, respectively.
#'@param ... additional arguments; currently "datetime" is the only recognized argument passed through \code{...}
#'@return
#'A data.frame with columns corresponding to components of metabolism
#'\describe{
	#'\item{GPP}{numeric estimate of Gross Primary Production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
	#'\item{R}{numeric estimate of Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
	#'\item{NEP}{numeric estimate of Net Ecosystem production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'}
#'If n.boot is 2 or greater, the data frame will also contain the following columns:
#'\describe{
#'\item{GPP.lci}{lower 95\% confidence interval of the estimate of Gross Primary Production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'\item{GPP.uci}{upper 95\% confidence interval of the estimate of Gross Primary Production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'\item{GPP.cv}{coefficient of variation of the bootstrapped estimates of Gross Primary Production}
#'\item{R.lci}{lower 95\% confidence interval of the estimate of Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'\item{R.uci}{upper 95\% confidence interval of the estimate of Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'\item{R.cv}{coefficient of variation of the bootstrapped estimates of Respiration. As R is typically negative, the CV is also negative.}
#'\item{NEP.lci}{lower 95\% confidence interval of the estimate of Net Ecosystem Production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'\item{NEP.uci}{upper 95\% confidence interval of the estimate of Net Ecosystem Production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'\item{NEP.cv}{coefficient of variation of the bootstrapped estimates of Net Ecosystem Production. As the mean NEP can approach 0, the CV can become inflated.}
#'}
#' The maximum likelihood estimates of model parameters can be accessed via \code{attributes(metab.mle(...))[["params"]]}
#'
#'@details
#'The model has the three parameters, \eqn{c_1, c_2, epsilon}{c1, c2, epsilon}, and has the form
#'
#'\deqn{v=k.gas/z.mix}{v=k.gas/z.mix}
#'
#'\deqn{a_t = c_1*irr_{t-1} + c_2*log_e(wtr_{t-1}) + v_{t-1}*do.sat_{t-1}}{a[t] = c1*irr[t-1] + c2*log(wtr[t-1]) + v[t-1]*do.sat[t-1]}
#'
#'\deqn{beta = e^{-v}}{beta = exp(-v)}
#'
#'\deqn{do.obs_t = a_t/v_{t-1} + -e^{-v_{t-1}}*a_t/v_{t-1} + beta_{t-1}*do.obs_{t-1} + epsilon_t}{do.obs[t] = a[t]/v[t-1] + -exp(-v[t-1])*a[t]/v[t-1] + beta[t-1]*do.obs[t-1] + epsilon[t]}
#'
#'
#' The above model is used during model fitting, but if gas flux is not integrated between time steps, those equations simplify to the following:
#'
#' \deqn{F_{t-1} = k.gas_{t-1}*(do.sat_{t-1} - do.obs_{t-1})/z.mix_{t-1}}{F[t-1] = k.gas[t-1]*(do.sat[t-1] - do.obs[t-1])/z.mix[t-1]}
#'
#'\deqn{do.obs_t=do.obs_{t-1}+c_1*irr_{t-1}+c_2*log_e(wtr_{t-1}) + F_{t-1} + epsilon_t}{do.obs[t] = do.obs[t-1] + c1*irr[t-1] + c2*log(wtr[t-1]) + F[t-1] + epsilon[t]}
#'
#'
#'The parameters are fit using maximum likelihood, and the optimization (minimization of the negative log likelihood function) is performed by \code{optim} using default settings.
#'
#'GPP is then calculated as \code{mean(c1*irr, na.rm=TRUE)*freq}, where \code{freq} is the number of observations per day, as estimated from the typical size between time steps. Thus, generally \code{freq==length(do.obs)}.
#'
#'Similarly, R is calculated as \code{mean(c2*log(wtr), na.rm=TRUE)*freq}.
#'
#'NEP is the sum of GPP and R.
#'
#'@note Currently, missing values in any arguments will result in an error, so freq must always equal nobs.
#'@author Luke A Winslow, Ryan Batt, Bennett McAfee, GLEON Fellows
#'@references
#'Hanson, Paul C., Stephen R. Carpenter, Nobuaki Kimura, Chin Wu, Sean P. Cornelius, 
#'and Timothy K. Kratz. 2008. \emph{Evaluation of Metabolism Models for Free-Water Dissolved 
#'Oxygen Methods in Lakes}. Limnology and Oceanography: Methods 6 (9): 454–65. 
#'\url{https://doi.org/10.4319/lom.2008.6.454}.
#'
#'Solomon, Christopher T., Denise A. Bruesewitz, David C. Richardson, et al. 2013. 
#'\emph{Ecosystem Respiration: Drivers of Daily Variability and Background Respiration in Lakes around the Globe}. 
#'Limnology and Oceanography 58 (3): 849–66. \url{https://doi.org/10.4319/lo.2013.58.3.0849}.
#'
#'Dugan, Hilary A., R. Iestyn Woolway, Arianto B. Santoso, et al. 2016. 
#'\emph{Consequences of Gas Flux Model Choice on the Interpretation of Metabolic 
#'Balance across 15 Lakes}. Inland Waters 6 (4): 581–92. \url{https://doi.org/10.1080/IW-6.4.836}.
#'
#'Corman, Jessica R., Jacob A. Zwart, Jennifer Klug, et al. 2023. \emph{Response of 
#'Lake Metabolism to Catchment Inputs Inferred Using High-Frequency Lake and 
#'Stream Data from across the Northern Hemisphere}. Limnology and Oceanography 
#'68 (12): 2617–31. \url{https://doi.org/10.1002/lno.12449}.
#'
#'@importFrom stats dnorm optim rnorm quantile
#'
#'@seealso
#'\link{metab}, \link{metab.bookkeep}, \link{metab.ols}, \link{metab.kalman}, \link{metab.bayesian}
#'@examples
#'library(rLakeAnalyzer)
#'doobs = load.ts(system.file('extdata',
#'                            'sparkling.doobs', package="LakeMetabolizer"))
#'wtr = load.ts(system.file('extdata',
#'                          'sparkling.wtr', package="LakeMetabolizer"))
#'wnd = load.ts(system.file('extdata',
#'                          'sparkling.wnd', package="LakeMetabolizer"))
#'irr = load.ts(system.file('extdata',
#'                          'sparkling.par', package="LakeMetabolizer"))
#'
#'#Subset a day
#'mod.date = as.POSIXct('2009-07-08', 'GMT')
#'doobs = doobs[trunc(doobs$datetime, 'day') == mod.date, ]
#'wtr = wtr[trunc(wtr$datetime, 'day') == mod.date, ]
#'wnd = wnd[trunc(wnd$datetime, 'day') == mod.date, ]
#'irr = irr[trunc(irr$datetime, 'day') == mod.date, ]
#'z.mix = ts.thermo.depth(wtr)
#'
#'k600 = k.cole.base(wnd[,2])
#'k.gas = k600.2.kGAS.base(k600, wtr[,3], 'O2')
#'do.sat = o2.at.sat.base(wtr[,3], altitude=300)
#'
#'metab.mle(doobs[,2], do.sat, k.gas, z.mix[,2], irr[,2], wtr[,3])
#'@export
metab.mle <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, error.type="OE", n.boot = 0, ar1.resids = TRUE, constrain.sign = FALSE, ...){
  
  complete.inputs(do.obs=do.obs, do.sat=do.sat, k.gas=k.gas,
                  z.mix=z.mix, irr=irr, wtr=wtr, error=TRUE)

	match.arg(error.type, choices=c('OE', 'PE'))

	nobs <- length(do.obs)

	mm.args <- list(...)

  if(any(z.mix <= 0)){
    stop("z.mix must be greater than zero.")
  }
  if(any(wtr <= 0)){
    stop("all wtr must be positive.")
  }
  
  if("datetime"%in%names(mm.args)){ # check to see if datetime is in the ... args
    datetime <- mm.args$datetime # extract datetime
    freq <- calc.freq(datetime) # calculate sampling frequency from datetime
    if(nobs!=freq){ # nobs and freq should agree, if they don't issue a warning
      bad.date <- format.Date(datetime[1], format="%Y-%m-%d")
      warning("number of observations on ", bad.date, " (", nobs, ") ", "does not equal estimated sampling frequency", " (", freq, ")", sep="")
    }
  }else{ # if datetime is *not* in the ... args
    warning("datetime not found, inferring sampling frequency from # of observations") # issue a warning (note checks in addNAs)
    # NOTE: because of the checks in addNA's, it is unlikely a user would receive this warning via metab()
    # warning will only be seen through direct use of metab.bookkeep when datettime is not supplied
    freq <- nobs
  }
  
  chk.list <- list(do.obs, irr, do.sat, z.mix, k.gas, wtr)
  if(!all(sapply(chk.list, is.numeric)) || !all(sapply(chk.list, is.vector))){
    stop('All metab.mle inputs must be numeric vectors.')
  }
  
  if(!all(nobs==sapply(chk.list, length))){
    stop('All input data to metab.mle must be the same length')
  }
  
  Q0 <- ((diff(range(do.obs,na.rm=TRUE)) - mean(do.obs,na.rm=TRUE))^2 / length(do.obs))
  
  if (constrain.sign == TRUE){
    guesses <- c(log(1E-4), log(1E-4), log(Q0))
  } else {
    guesses <- c(1E-4, 1E-4, log(Q0))
  }
  
  
  #We have a different number of fitted parameters depending on error type of the model
  if(error.type=='OE'){
    guesses <- c(guesses, do.obs[1])
    
    fit <- optim(guesses, fn=mleNllOE, do.obs=do.obs, do.sat=do.sat, k.gas=(k.gas/freq), z.mix=z.mix, irr=irr, wtr=wtr, constrain.sign = constrain.sign)
    
    pars0 <- fit$par
    if (constrain.sign == TRUE){
      pars <- c("gppCoeff"=exp(pars0[1]), "rCoeff"=-exp(pars0[2]), "Q"=exp(pars0[3]), "nll"=fit$value, "doInit"=pars0[4])
    } else {
      pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]), "nll"=fit$value, "doInit"=pars0[4])
    }
    
    
  }else if(error.type=='PE'){
    guesses <- c(guesses)
    
    fit <- optim(guesses, fn=mleNllPE, do.obs=do.obs, do.sat=do.sat, k.gas=(k.gas/freq), z.mix=z.mix, irr=irr, wtr=wtr, constrain.sign = constrain.sign)
    
    pars0 <- fit$par
    if (constrain.sign == TRUE){
      pars <- c("gppCoeff"=exp(pars0[1]), "rCoeff"=-exp(pars0[2]), "Q"=exp(pars0[3]), "nll"=fit$value, "doInit"=pars0[4])
    } else {
      pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]), "nll"=fit$value, "doInit"=pars0[4])
    }
    
  }else{
    stop("error.type must be either 'OE' or 'PE', Observation Error or Process Error respectively.")
  }
  
  # ====================================
  # = Use fits to calculate metabolism =
  # ====================================
  GPP <- mean(pars[1]*irr, na.rm=TRUE) * freq
  R <- mean(pars[2]*log(wtr), na.rm=TRUE) * freq
  
  # =================================
  # = Bootstrapping for uncertainty =
  # =================================
  if (n.boot == 0){
    results <- list("params"=pars, "metab"=c("GPP"=GPP,"R"=R,"NEP"=GPP+R))
  } else if (n.boot >= 2){
    boot.results <- bootstrap.mle(n.boot, do.obs, pars, guesses, do.sat, k.gas, freq, z.mix, irr, wtr, error.type, ar1.resids, constrain.sign)
    ci.GPP <- stats::quantile(boot.results$GPP, c(0.025, 0.975), na.rm=TRUE) # 95% confidence interval
    ci.R   <- stats::quantile(boot.results$R,   c(0.025, 0.975), na.rm=TRUE) # 95% confidence interval
    ci.NEP <- stats::quantile(boot.results$NEP, c(0.025, 0.975), na.rm=TRUE) # 95% confidence interval
    cv.GPP <- stats::sd(boot.results$GPP, na.rm=TRUE)/mean(boot.results$GPP, na.rm=TRUE) # coefficient of variation
    cv.R   <- stats::sd(boot.results$R, na.rm=TRUE)/mean(boot.results$R, na.rm=TRUE) # coefficient of variation
    cv.NEP <- stats::sd(boot.results$NEP, na.rm=TRUE)/mean(boot.results$NEP, na.rm=TRUE) # coefficient of variation
    results <- list("params"=pars, 
                "metab"=c("GPP"=GPP,"R"=R,"NEP"=GPP+R,
                          "GPP.lci" = unname(ci.GPP[1]), "GPP.uci" = unname(ci.GPP[2]), "GPP.cv" = cv.GPP,
                          "R.lci" = unname(ci.R[1]), "R.uci" = unname(ci.R[2]), "R.cv" = cv.R,
                          "NEP.lci" = unname(ci.NEP[1]), "NEP.uci" = unname(ci.NEP[2]), "NEP.cv" = cv.NEP))
  } else {
    stop("n.boot must be 0 (no uncertainty quantification) or greater than 2 (uncertainty quantification via bootstrapping).")
  }
  
  return(results)
}

# ============================================
# = The R loop for Observation Error mle NLL =
# ============================================
mleLoopOE <- function(alpha, doobs, c1, c2, beta, irr, wtr, kz, dosat){
	nobs <- length(doobs)
	a.loop <- .C("mleLoopCoe", alpha=as.double(alpha), as.double(doobs), as.double(c1), as.double(c2), as.double(beta), as.double(irr), as.double(wtr), as.double(kz), as.double(dosat), as.integer(nobs), PACKAGE="LakeMetabolizer")
	return(a.loop[["alpha"]])
}

# ============================================
# = The R loop for Process Error mle NLL =
# ============================================
mleLoopPE <- function(alpha, doobs, c1, c2, beta, irr, wtr, kz, dosat){
	nobs <- length(doobs)
	a.loop <- .C("mleLoopCpe", alpha=as.double(alpha), as.double(doobs), as.double(c1), as.double(c2), as.double(beta), as.double(irr), as.double(wtr), as.double(kz), as.double(dosat), as.integer(nobs), PACKAGE="LakeMetabolizer")
	return(a.loop[["alpha"]])
}

# ====================
# = mle NLL function =
# ====================
mleNllPE <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr, constrain.sign){
	if (constrain.sign == TRUE){
	  c1 <- exp(Params[1]) #PAR coeff
	  c2 <- -exp(Params[2]) #log(Temp) coeff
	} else {
	  c1 <- Params[1] #PAR coeff
	  c2 <- Params[2] #log(Temp) coeff
	}
	Q <- exp(Params[3]) # Variance of the process error

	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	beta <- exp(-kz) # This beta is for using the differential equation form

	# Set first true value equal to first observation
	alpha <- rep(0, length(do.obs))
	alpha[1] <- do.obs[1]#Let's give this model some starting values

	#R version of C loop
	#for(i in 2:length(do.obs)){
	#	a1 <- c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
	#	alpha[i] <- a1/kz[i-1] + -exp(-kz[i-1])*a1/kz[i-1] + beta[i-1]*alpha[i-1] # NOTE: beta==exp(-kz); kz=K/Zmix
	#}
	alpha <- mleLoopPE(alpha=alpha, doobs=do.obs, c1=c1, c2=c2, beta=beta, irr=irr, wtr=wtr, kz=kz, dosat=do.sat)

	return(-sum(dnorm(do.obs, alpha, sd=sqrt(Q), log=TRUE), na.rm=TRUE))
}#End function

# ====================
# = mle NLL function =
# ====================
mleNllOE <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr, error.type, constrain.sign){
  if (constrain.sign == TRUE){
    c1 <- exp(Params[1]) #PAR coeff
    c2 <- -exp(Params[2]) #log(Temp) coeff
  } else {
    c1 <- Params[1] #PAR coeff
    c2 <- Params[2] #log(Temp) coeff
  }
	Q <- exp(Params[3]) # Variance of the process error

	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	beta <- exp(-kz) # This beta is for using the differential equation form

	# Set first true value equal to first observation
	alpha <- rep(0, length(do.obs))
	alpha[1] <- Params[4] #Free varying initial DO value

	alpha <- mleLoopOE(alpha=alpha, doobs=do.obs, c1=c1, c2=c2, beta=beta, irr=irr, wtr=wtr, kz=kz, dosat=do.sat)

	return(-sum(dnorm(do.obs, alpha, sd=sqrt(Q), log=TRUE), na.rm=TRUE))
}#End function

# =========================================
# = Function to bootstrap for uncertainty =
# =========================================
bootstrap.mle <- function(n.boot, do.obs, pars, guesses, do.sat, k.gas, freq, z.mix, irr, wtr, error.type = "OE", ar1.resids = FALSE, constrain.sign = FALSE){
  
  ## Calculate fitted values and residuals
  if (error.type == "OE"){
    alpha <- mleLoopOE(alpha = c(do.obs[1], rep(0, length(do.obs) - 1)),
                       doobs = do.obs,
                       c1 = pars[1],
                       c2 = pars[2],
                       beta = exp(-(k.gas / freq / z.mix)),
                       irr = irr,
                       wtr = wtr,
                       kz = k.gas / freq / z.mix,
                       dosat = do.sat)
  } else if (error.type == "PE"){
    alpha <- mleLoopPE(alpha = c(do.obs[1], rep(0, length(do.obs) - 1)),
                       doobs = do.obs,
                       c1 = pars[1],
                       c2 = pars[2],
                       beta = exp(-(k.gas / freq / z.mix)),
                       irr = irr,
                       wtr = wtr,
                       kz = k.gas / freq / z.mix,
                       dosat = do.sat)
  } else {
    stop("error.type must be either 'OE' or 'PE', Observation Error or Process Error respectively.")
  }
  
  
  resids <- do.obs - alpha
  n.obs <- length(do.obs)
  
  if (ar1.resids){
    ar1.lm <- lm(resids[2:n.obs] ~ resids[1:(n.obs - 1)] - 1)
    ar1.coeff <- as.numeric(ar1.lm$coefficients)
    ar1.sd <- sd(ar1.lm$residuals)
  }
  
  ## Initialize results object
  boot.results <- data.frame(boot.iter = seq_len(n.boot),
                             gppCoeff = rep(NA, n.boot),
                             rCoeff = rep(NA, n.boot),
                             Q = rep(NA, n.boot),
                             doInit = rep(NA, n.boot),
                             convergence = rep(NA, n.boot),
                             nll = rep(NA, n.boot),
                             GPP = rep(NA, n.boot),
                             R = rep(NA, n.boot),
                             NEP = rep(NA, n.boot))
  
  ## Bootstrap loop
  for (i in boot.results$boot.iter){
    
    if (ar1.resids) {
      # Residual randomization preserving AR(1)
      simRes <- rep(NA, n.obs)
      simRes[1] <- sample(resids[!is.na(resids)], 1)
      for (j in 2:n.obs){simRes[j] <- ar1.coeff * simRes[j - 1] + stats::rnorm(1, sd = ar1.sd)}
    } else {
      # Residual randomization without replacement
      simRes <- sample(resids[!is.na(resids)], length(resids), replace = FALSE)
    }
    
    doSim <- alpha + simRes
    
    if (error.type == "OE"){
      simFit <- optim(guesses, fn = mleNllOE, do.obs = doSim, do.sat = do.sat, k.gas = k.gas / freq, z.mix = z.mix, irr = irr, wtr = wtr, constrain.sign = constrain.sign)
    } else if (error.type == "PE"){
      simFit <- optim(guesses, fn = mleNllPE, do.obs = doSim, do.sat = do.sat, k.gas = k.gas / freq, z.mix = z.mix, irr = irr, wtr = wtr, constrain.sign = constrain.sign)
    }
    
    simFitPar <- simFit$par
    if (constrain.sign == TRUE){
      simFit.gppCoeff <- exp(simFitPar[1])
      simFit.rCoeff <- -exp(simFitPar[2])
    } else {
      simFit.gppCoeff <- simFitPar[1]
      simFit.rCoeff <- simFitPar[2]
    }
    simFitGPP <- mean(simFit.gppCoeff * irr, na.rm = TRUE) * freq
    simFitR <- mean(simFit.rCoeff * log(wtr), na.rm = TRUE) * freq
    
    boot.results[i, "gppCoeff"] <- simFit.gppCoeff
    boot.results[i, "rCoeff"] <- simFit.rCoeff
    boot.results[i, "Q"] <- exp(simFitPar[3])
    boot.results[i, "doInit"] <- simFitPar[4]
    boot.results[i, "convergence"] <- simFit$convergence
    boot.results[i, "nll"] <- simFit$value
    boot.results[i, "GPP"] <- simFitGPP
    boot.results[i, "R"] <- simFitR
    boot.results[i, "NEP"] <- simFitGPP + simFitR
  }
  
  return(boot.results)
}
