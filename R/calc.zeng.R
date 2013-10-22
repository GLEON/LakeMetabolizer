calc.zeng <- function(dateTime,Ts,airT,Uz,rh,atm.press,hu,ht,hq){
  
  # INPUTS
  #   dateTime = datetime, YYYY-mm-dd HH:MM
  #   Ts = surface water temperature, degC
  #   Uz = wind speed, m/s
  #   rh = relative humidity, %
  #   airT = air temperature, degC
  #   hu: height of wind measurement
  #   ht: height of temperature measurement
  #   hq: height of humidity measurement
  #   atm.press: atmospheric pressure (mb)
  #
  # OUTPUTS:
  #   mm: matrix of sensible and latent heat fluxes
  #                mm also contains other variables used in calculating 
  #                these fluxes. 
  
  # Author R.Iestyn Woolway

  # define functions used in script
  psi <- function(k,zeta){
    chik <- (1 - 16*zeta)^0.25
    if (k == 1){
      psi <- 2*log((1 + chik)*0.5) + log((1+chik*chik)*0.5) - 2*atan(chik) + (pi/2)
    } else{
      psi <- 2*log((1 + chik*chik)*0.5)
    }
  }
  
  # generate data.frame of all variables (to ensure consistent times)
  dat <- data.frame(dateTime = dateTime,
                    Ts = Ts,
                    airT = airT,
                    Uz = Uz,
                    rh = rh)
  
  # remove duplicated time stamps
  dat$dateTime <- as.POSIXct(strptime(dat$dateTime,"%Y-%m-%d %H:%M")) # ensure times are POSIXct
  dat <- subset(dat,!duplicated(dat$dateTime)) #remove duplicate time stamps
  
  # store original dates - used for final data frame
  original_dates <- data.frame(dateTime = dat$dateTime)  
  
  # remove missing data
  dat <-  dat[complete.cases(dat),]
  
  # re-define variables
  Ts <- dat$Ts
  airT <- dat$airT
  Uz <- dat$Uz
  rh <- dat$rh
  
  # if temperature and humidity height are missing, assume same as wind
  if (missing(ht)){ht <- wndZ}
  if (missing(hq)){hq <- wndZ}

  # define constants
  const_vonKarman <- 0.41 # von Karman constant
  const_gas <- 287.1 # gas constant for dry air J kg-1 K-1
  const_SpecificHeatAir <- 1005 # Specific heat capacity of air, J kg-1 K-1
  const_Charnock <- 0.013 # charnock constant
  const_Gravity <- 9.81 # gravitational acceleration, m/s2

  # ensure wind below 0.2 are changed to 0.1
  thres <- 0.2
  Uz[Uz < thres] <- thres/2 # m/s

  # calculate humidity values
  e_s <- 6.11*exp(17.27*airT/(237.3 + airT)) # saturated vapour pressure at airT, mb
  e_a <- rh*e_s/100 # vapour pressure, mb
  q_z <- 0.622*e_a/atm.press # specific humidity, kg kg-1 
  e_sat <- 6.11*exp(17.27*Ts/(237.3 + Ts)) # saturated vapour pressure at Ts, mb
  q_s <- 0.622*e_sat/atm.press # humidity at saturation, kg kg-1

  # calculate gas constant for moist air, J kg-1 K-1
  R_a <- 287*(1 + 0.608*q_z)

  # -- calculate latent heat of vaporization, J kg-1 
  xlv <- 2.501e6-2370*Ts    

  # calculate air density, kg/m3
  rho_a <- 100*atm.press/(R_a*(airT + 273.16))

  # -- kinematic viscosity, m2 s-1
  KinV <- (1/rho_a)*(4.94e-8*airT + 1.7184e-5)   

  # calculate virtual air temperature, K
  t_virt <- (airT + 273.16)*(1 + 0.61*q_z)

  # estimate initial values of u* and zo
  ustar <- Uz*sqrt(0.00104+0.0015/(1+exp((-Uz+12.5)/1.56)))
  zo <- (const_Charnock*ustar^2/const_Gravity) + (0.11*KinV/ustar)
  for (i in 1:20){
    ustar <- const_vonKarman*Uz/(log(hu/zo))
    zo <- (const_Charnock*ustar^2/const_Gravity) + (0.11*KinV/ustar)    
  }
  zo <- Re(zo)
  
  # calculate neutral transfer coefficients
  C_DN <- (ustar^2)/(Uz^2)
  re <- ustar*zo/KinV
  zot <- zo*exp(-2.67*(re)^(0.25) + 2.57)
  zot <- Re(zot)
  zoq <- Re(zot)
  C_HN <- const_vonKarman*sqrt(C_DN)/(log(hu/zot)) 
  C_EN <- Re(C_HN)

  # calculate neutral transfer coefficients at 10 m
  C_D10N <- (const_vonKarman/log(10/zo))*(const_vonKarman/log(10/zo)) 
  C_E10N <- (const_vonKarman*const_vonKarman)/(log(10/zo)*log(10/zoq))
  C_H10N <- C_E10N    

  # calculate neutral latent and sensible heat fluxes, W/m2
  alhN <- rho_a*xlv*C_EN*Uz*(q_s-q_z)
  ashN <- rho_a*const_SpecificHeatAir*C_HN*Uz*(Ts-airT) 

  # calculate initial monin obukhov length scale, m
  obu <- (-rho_a*t_virt*(ustar*ustar*ustar))/(const_vonKarman*const_Gravity*(ashN/const_SpecificHeatAir + 0.61*(airT + 273.16)*alhN/xlv)) 

  # iteration to compute corrections for atmospheric stability    
  zeta_thres <- 15
  
  # pre-define arrays
  tstar <- rep(0,length(Uz))
  qstar <- rep(0,length(Uz))
  wc <- rep(0,length(Uz))
  
  for (i in 1:20){ 
    
    # calulate roughness lengths
    zo <- (0.013*((ustar^2)/const_Gravity)) + (0.11*(KinV/ustar))                       
    re <- (ustar*zo)/KinV     
    xq <- 2.67*re^0.25 - 2.57
    xq[xq < 0] <- 0
    zoq <- zo/exp(xq)        
    zot <- zoq     

    # define zeta
    zetam <- -1.574
    zetat <- -0.465
  
    # calculate ustar
    zeta <- hu/obu
    
    zeta[zeta < -zeta_thres] <- -zeta_thres
    zeta[zeta > zeta_thres] <- zeta_thres
    
    idx <- zeta < zetam & !is.na(zeta) # very unstable conditions      
    ustar[idx] <- (Uz[idx]*const_vonKarman)/((log((zetam*obu[idx])/zo[idx]) - psi(1,zetam)) + 1.14*(((-zeta[idx])^0.333) - ((-zetam)^0.333))) 
    idx <- zeta < 0 & zeta >= zetam & !is.na(zeta) # unstable conditions
    ustar[idx] = (Uz[idx]*const_vonKarman)/(log(hu/zo[idx]) - psi(1,zeta[idx]))
    idx <- zeta > 0 & zeta <= 1 & !is.na(zeta) # stable conditions
    ustar[idx] <- (Uz[idx]*const_vonKarman)/(log(hu/zo[idx]) + 5*zeta[idx])
    idx <- zeta > 1 & !is.na(zeta) # very stable conditions
    ustar[idx] <- (Uz[idx]*const_vonKarman)/((log(obu[idx]/zo[idx])+ 5) + (5*log(zeta[idx]) + zeta[idx] - 1))
    
    # calculate tstar
    zeta <- ht/obu
    zeta[zeta < -zeta_thres] <- -zeta_thres
    zeta[zeta > zeta_thres] <- zeta_thres

    idx <- zeta < zetat & !is.na(zeta) # very unstable conditions
    tstar[idx] <- (const_vonKarman*(airT[idx] - Ts[idx]))/((log((zetat*obu[idx])/zot[idx]) - psi(2,zetat)) +  0.8*((-zetat)^-0.333 - ((-zeta[idx]))^-0.333))
    idx <- zeta >= zetat & zeta < 0 & !is.na(zeta) # unstable conditions
    tstar[idx] <- (const_vonKarman*(airT[idx] - Ts[idx]))/(log(ht/zot[idx]) - psi(2,zeta[idx]))
    idx <- zeta > 0 & zeta <= 1 & !is.na(zeta) # stable conditions
    tstar[idx] <- (const_vonKarman*(airT[idx] - Ts[idx]))/(log(ht/zot[idx]) + 5*zeta[idx])
    idx <- zeta > 1 & !is.na(zeta) # very stable conditions
    tstar[idx] <- (const_vonKarman*(airT[idx] - Ts[idx]))/((log(obu[idx]/zot[idx]) + 5) + (5*log(zeta[idx]) + zeta[idx] - 1))

    # calculate qstar
    zeta <- hq/obu
    zeta[zeta < -zeta_thres] <- -zeta_thres
    zeta[zeta > zeta_thres] <- zeta_thres

    idx <- zeta < zetat & !is.na(zeta) # very unstable conditions
    qstar[idx] <- (const_vonKarman*(q_z[idx] - q_s[idx]))/((log((zetat*obu[idx])/zoq[idx]) - psi(2,zetat)) + 0.8*((-zetat)^-0.333 - ((-zeta[idx]))^-0.333))
    idx <- zeta >= zetat & zeta < 0 & !is.na(zeta) # unstable conditions
    qstar[idx] <- (const_vonKarman*(q_z[idx] - q_s[idx]))/(log(hq/zoq[idx]) - psi(2,zeta[idx]))
    idx <- zeta > 0 & zeta <= 1 & !is.na(zeta) # stable conditions
    qstar[idx] <- (const_vonKarman*(q_z[idx] - q_s[idx]))/(log(hq/zoq[idx]) + 5*zeta[idx])
    idx <- zeta > 1 & !is.na(zeta) # very stable conditions
    qstar[idx] <- (const_vonKarman*(q_z[idx] - q_s[idx]))/((log(obu[idx]/zoq[idx]) + 5) + (5*log(zeta[idx]) + zeta[idx] - 1))
    
    # calculate zeta at 10 m
    zeta <- 10/obu
    zeta[zeta < -zeta_thres] <- -zeta_thres
    zeta[zeta > zeta_thres] <- zeta_thres

    # calculate transfer coefficients corrected for atmospheric stability
    C_H <- (-rho_a*const_SpecificHeatAir*ustar*tstar)/(rho_a*const_SpecificHeatAir*Uz*(Ts - airT))
    C_E <- C_H
    C_D <- (ustar*ustar)/(Uz*Uz)

    # calculate tau and sensible and latent heat fluxes
    tau <- C_D*rho_a*Uz*Uz
    ash <- rho_a*const_SpecificHeatAir*C_H*Uz*(Ts - airT)
    alh <- rho_a*xlv*C_E*Uz*(q_s - q_z)

    # calculate new monin obukhov length
    obu <- (-rho_a*t_virt*(ustar*ustar*ustar))/(const_Gravity*const_vonKarman*((ash/const_SpecificHeatAir) + (0.61*(airT + 273.16)*alh/xlv)))

    # alter zeta in stable cases
    zeta <- hu/obu
    idx <- zeta >= 1
    Uz[idx] <- max(Uz[idx],0.1)

    # avoid singularity at um = 0 for unstable conditions     
    idx <- zeta < 0 & !is.na(zeta)
    th <- (airT + 273.16)*(1000/atm.press)^(287.1/1004.67) # potential temperature  
    thvstar <- tstar*(1 + 0.61*q_z/1000) + 0.61*th*qstar # temperature scaling parameter
    thv <- th*(1 + 0.61*q_z/1000) #virtual potential temperature    
    wc[idx] <- 1*(-const_Gravity*ustar[idx]*thvstar[idx]/thv[idx])^0.333
    Uz[idx] <- sqrt(Uz[idx]*Uz[idx] + wc[idx]*wc[idx])
  }

  # take real values to remove any complex values that arise from missing data or NA.
  C_D <- Re(C_D)
  C_E <- Re(C_E)
  C_H <- Re(C_H)
  zo <- Re(zo)
  zoq <- Re(zoq)
  zot <- Re(zot)
  
  # store results in data.frame and merge with original dateTime
  mm <- data.frame(dateTime = dat$dateTime,
                   Ts = Ts,airT = airT,rh = rh,Uz = Uz,
                   C_D = C_D,C_E = C_E,C_H = C_H,zo = zo,zot = zot,
                   zoq = zoq,ustar = ustar,alh = alh,ash = ash)
  mm <- merge(mm,original_dates,all = TRUE)
  return(mm)
  
}
