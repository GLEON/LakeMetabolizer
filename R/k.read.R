#calculateK600 R script - conversion from MatLab 'calculateK600.m' by Jordan Read
#author: Hilary Dugan 
#Edits 2013-09-10: Luke Winslow
#updated: R.I.Woolway Oct 2013

k.read <- function(wndZ, Kd, lat, lake.area, atm.press, dateTime, wtr, depth, airT, Uz, RH, sw, lwnet, par, lw){ 
  
  require(rLakeAnalyzer)
  # define constants used in function
  dT <- 0.5   # change in temp for mixed layer depth
  C1 <- 114.278 # from Soloviev et al. 2007
  nu <- 0.29 # proportionality constant from Zappa et al. 2007, lower bounds
  KeCrit <- 0.18     # constant for wave age = 20 (Soloviev et al. 2007)
  albedo_SW <- 0.07
  swRat <- 0.46 # percentage of SW radiation that penetrates the water column
  g <- 9.81 # gravity
  C_w <- 4186 # J kg-1 ?C-1 (Lenters et al. 2005)
  mnWnd <- 0.2 # minimum wind speed
  Kelvin <- 273.15 # temp mod for deg K   
  emiss <- 0.972 # emissivity;
  S_B <- 5.67E-8 # Stefan-Boltzman constant (?K is used)
  
  
  # Get short wave radiation data 
  if(!missing(sw)){ 
    sw <- sw
  } else if (!missing(par)){
    #sw <- par
    #parMult <- 0.4957
    sw = par.to.sw(par)
    #sw <- sw*parMult
  } else {  
    stop("no SW equivalent file available\n")
  }
  
  # Get water temperature data
  if(!missing(wtr)){ 
    wtr <- wtr
    Ts <- wtr[1]
  } else {  
    stop("no wtr file available\n")
  }
  
  # Get air temperature
  if(!missing(airT)){ 
    airT <- airT
  } else {  
    stop("no air temp data available")
  }
  
  # Get relative humidity data
  if(!missing(RH)){ 
    RH <- RH
  } else {  
    stop("no relative humidity data available")
  }
  
  # Get long wave radiation data
  if(!missing(lwnet)){ 
    lwnet <- lwnet
  } else if(!missing(lw)){
    lw_in <- lw # long wave in
    Tk <- Ts+Kelvin # water temperature in Kelvin
    LWo <- S_B*emiss*Tk^4 # long wave out
    lwnet <- lw_in-LWo
  } else {  
    stop("no longwave radiation available")
  }
  
  # Get wind speed data
  if(!missing(Uz)){ 
    wnd <- Uz
  } else{  
    stop("no wind speed data available")
  }
  
  # impose limit on wind speed
  rpcI <- wnd < mnWnd
  wnd[rpcI] <- mnWnd
  
  # calculate sensible and latent heat fluxes
  mm <- calc.zeng(dateTime,Ts,airT,wnd,RH,atm.press,wndZ)
  C_D <- mm$C_D # drag coefficient for momentum
  E <- mm$alh # latent heat flux
  H <- mm$ash # sensible heat flux
  
  # calculate total heat flux
  dUdt <- sw*0.93 - E - H + lwnet
  Qo <- sw*(1-albedo_SW)*swRat
  
  # calculate water density
  rho_w <- water.density(Ts)
  
  # calculate u*
  vonK <- 0.41 # von Karman  constant
  if (wndZ != 10) {
    e1 <- sqrt(C_D)
    wnd <- wnd/(1-e1/vonK*log(10/wndZ))
  }
  rhoAir <- 1.2 #  air density
  tau <- C_D*wnd^2*rhoAir
  uSt <- sqrt(tau/rho_w)
  
  # find Z_aml
  if(!is.na(wtr[1] - wtr[length(wtr)]) && wtr[1] - wtr[length(wtr)] <= dT){
    z_aml <- depth[length(depth)]
  } else {
    zI <- depth[wtr[1] - dT > wtr]
    z_aml <- zI[1]
  }
  
  # calculate the effective heat flux
  q1 <- 2-2*exp(z_aml*-Kd)
  q2 <- z_aml*Kd
  q3 <- exp(z_aml*-Kd)
  H_star <- dUdt-Qo*(q1/q2-q3) # Kim 1976
  
  # calculate the thermal expansion coefficient 
  thermalExpFromTemp <- function(Ts) {    
    V <- 1       
    dT <- 0.001 
    T1 <- water.density(Ts)
    T2 <- water.density(Ts+dT)
    V2 <- T1/T2
    dv_dT <- (V2-V)/dT
    alpha <- dv_dT
    return (alpha)
  }
  tExp <- thermalExpFromTemp(Ts)
  
  # calculate buoyancy flux and w*
  B1 <- H_star*tExp*g
  B2 <- rho_w*C_w
  Bflx <- B1/B2
  ltI <- Bflx>0
  B1 <- Bflx
  B1[ltI] <- 0
  divi <- 1/3
  w1 <- -B1*z_aml
  wSt <- w1^divi
  
  # calculate kinematic viscosiy
  getKinematicVis <- function(Ts) {
    # from Mays 2005, Water Resources Engineering
    tempTable <- seq(0,100,by=5)
    # table in m2/s E-6
    visTable <- c(1.792,1.519,1.308,1.141,1.007,0.897,
                  0.804,0.727,0.661,0.605,0.556,0.513,0.477,0.444,
                  0.415,0.39,0.367,0.347,0.328,0.311,0.296)
    v <- data.frame(approx(tempTable,visTable,xout = Ts))[2]
    v <- v*1e-6
    return(v)
  }
  kinV <- getKinematicVis(Ts)
  
  KeDe <- (kinV*g)
  KeNm <- uSt^3
  Ke <- KeNm/KeDe
  tau <- tau    # estimate of total tau (includes wave stress)
  euPw <- (1+Ke/KeCrit)  # tau_shear = tau/(1+Ke/Kecr)
  tau_t <- tau/euPw      # tau_shear, Soloviev
  uTanS <- tau_t/rho_w   
  uTanS <- uTanS^0.5
  
  # calculate viscous sublayer
  Sv <- C1*kinV/uTanS
  eu_N <- uTanS^3      # e_u(0) = (tau_t/rho)^1.5/(vonK*Sv)
  eu_D <- vonK*Sv       # denominator
  eu_0 <- eu_N/eu_D    # in m2/s3
  ew_0 <- -1.0*B1       # buoyancy flux, but only when outward
  e_0 <- ew_0+eu_0     # e(0) from Soloviev (w/o wave effects)
  K1 <- e_0*kinV       # in units of m4/s4, want cm4/hr4
  K2 <- ew_0*kinV      # convective component (m4/s4)
  K1 <- K1*100^4*3600^4 # now in cm4/hr4  (Total)
  K2 <- K2*100^4*3600^4 # now in cm4/hr4  (Convective)
  K600 <- nu*600^(-0.5)*K1^(1/4)   # in cm/hr (Total)
  
  k600 <- as.numeric(K600)
  
  return(k600)
}
