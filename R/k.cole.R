# ---Author: Jake Zwart, 2013-09-30 --- 
# cleaned from previous script - original author: Arianto Santoso 

# wnd: wind value in m/s
# wtr: water temperature at dissolved oxygen sensor depth in degrees celsius 
# wndHeight: height in meters of wind measurement 

# OUTPUT: returns the gas exchange velocity for O2 in units of m/(timeStep*min) (i.e. 30 minute sampling 
#          interval will return kO2 in units of m/(1/48) - converts to fraction of day)

k.cole <- function(wnd,wtr,wndHeight){
  U10 <- wnd * (10/wndHeight)^1/7
  k600 <- 2.07 + (0.215 * (U10^1.7)) # units in cm h-1
  k600 <- k600*24/100 #units in m d-1
  ScO2 <- 1800.6 - (120.1 * wtr) + (3.7818 *wtr^2) - (0.047608 * wtr^3) # Schmidt number for O2 (Waninkhof, 1992)
  Sc600 <- ScO2 / 600
  n <- ifelse(U10 < 3.7,0.5,0.67) # n = 0.5 if U10 is less than 3.7 ms^-1, and 0.67 for U10 > 3.7 ms^-1
  kO2 <- k600 * (Sc600)^(-n) # gas exchange for O2 m d-1 
  kO2 <- kO2*(timeStep/1440) #change kO2 to units of m/(timeStep*min)
  return(kO2)
}

# -- References 
#** COLE, J.J. and CARACO, N.F. 1998. Atmospheric exchange of **
#** carbon dioxide in a low-wind oligotrophic lake measured by **
#** the addition of SF6. Limnology and Oceanography. 43(4): 647-656. 

#** WANNINKHOF, R. 1992. Relationship between gas exchange and **
#** wind speed over the ocean. J. Geophys. Res. 97: 7373-7382.
