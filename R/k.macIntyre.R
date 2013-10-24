# ---Author: Hilary Dugan, 2013-10-20 --- 
# modified from K600 code by Richard Woolway and Jordan Read

# INPUTS;
# wndz: Height of anemometer
# Kd: Diffuse attenuation coefficient
# atm.press: Atmospheric pressure in hPa or Mb. If unknown use '1013' for sea level 
# dateTime: date and time vector
# wtr: dataframe of water temperatures 
# airT: vector of air temperatures in *C
# Uz: vector of wind speed in m/s
# Rh: vector of relative humidity in %
# sw: vector of shortwave radiation in W/m2
# lw: vector of longwave radiation. If missing, run calc.lw.net function first
# par: vector of par data in umol m-2 s-1 

# OUTPUT: returns the gas exchange velocity for O2 in units of m/(timeStep*min) (i.e. 30 minute sampling 
#          interval will return kO2 in units of m/(1/48) - converts to fraction of day)

k.macIntyre <- function(wndZ, Kd, atm.press, dateTime, wtr, depth, airT, Uz, RH, sw, lwnet, par){
  
  #Constants
  dT <- 0.5   # change in temp for mixed layer depth. Step change in temperature from the surface
  #temperature is set equivalent to the accuracy of the loggers.
  albedo_SW <- 0.07
  vonK <- 0.41 #von Karman constant
  swRat <- 0.46 # percentage of SW radiation that penetrates the water column
  mnWnd <- 0.2 # minimum wind speed
  g <- 9.81 # gravity
  C_w <- 4186 # J kg-1 ?C-1 (Lenters et al. 2005)
 
  # Get short wave radiation data 
  if(!missing(sw)){ 
    sw <- sw
  } else if (!missing(par)){
    sw <- par
    parMult <- 0.4957
    sw <- sw*parMult
  } else {  
    stop("no SW equivalent file available\n")
  }
  
  # Get water temperature data
  if(!missing(wtr)){ 
    wtr <- wtr
    Ts <- wtr[[1]]
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
  mm <- calc.zeng(dateTime,Ts,airT,wnd,RH,atm.press,wndZ,2,2)
  C_D <- mm$C_D # drag coefficient for momentum
  E <- mm$alh # latent heat flux
  H <- mm$ash # sensible heat flux
  
  # calculate total heat flux
  dUdt <- sw*0.93 - E - H + lwnet
  Qo <- sw*(1-albedo_SW)*swRat
  
  # calculate water density
  rho_w <- water.density(Ts)
  
  # calculate u*
  if (wndZ != 10) {
    e1 <- sqrt(C_D)
    u10 <- wnd/(1-e1/vonK*log(10/wndZ))
  }
  rhoAir <- 1.2 #  air density
  vonK <- 0.41 # von Karman  constant
  tau <- C_D*u10^2*rhoAir
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
  
  B1 = H_star*tExp*g
  B2 = rho_w*C_w
  Bflx = B1/B2
  
  
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
  
  
  KeNm = uSt^3
  SmE   = 0.84*(-0.58*Bflx+1.76*KeNm/(vonK*z_aml))
  Sk   = SmE*kinV
  Sk   = Sk*100^4*3600^4 # Sally's K now in cm4/h4
  Sk600= 1.2*600^(-0.5)*Sk^(1/4) # in cm/hr (Total)
    
  k600 <- as.numeric(Sk600)
  k600 <- k600*24/100 #units in m d-1
  ScO2 <- 1800.6 - (120.1 * Ts) + (3.7818 *Ts^2) - (0.047608 * Ts^3) # Schmidt number for O2 (Waninkhof, 1992)
  Sc600 <- ScO2 / 600
  n <- ifelse(wnd < 3.7,0.5,0.67) # n = 0.5 if U10 is less than 3.7 ms^-1, and 0.67 for U10 > 3.7 ms^-1
  kO2 <- k600 * (Sc600)^(-n) # gas exchange for O2 m d-1 
  kO2 <- kO2*(timeStep/1440) #change kO2 to units of m/(timeStep*min)
  return(kO2)
}

# -- References 
#MACINTRYE, sALLY, ANDERS JONSSON, MATS JANSSON, JAN ABERG, DAMON E. TURNEY AND SCOTT D. MILLER.
#2010. Buoyancy flux, turbulence, and the gas transfer coefficient in a stratified lake.
#Geophysical Research Letters. 37: L24604. doi:10.1029/2010GL044164 
