# ---Author: Hilary Dugan, 2013-10-20 --- 
# adapted from k.cole.R function

# Uses power function from:
# -- References 
#CRUSIUS, JOHN, AND RIK WANNINKHOF. 2003
#Gas transfer velocities measured at low wind speed over a lak.
#Limnology and Oceanography. 48(3): 1010:1017.

# INPUTS
# wnd: wind value in m/s
# wtr: water temperature at dissolved oxygen sensor depth in degrees celsius 
# wndHeight: height in meters of wind measurement 

# OUTPUT: returns the gas exchange velocity for O2 in units of m/(timeStep*min) (i.e. 30 minute sampling 
# interval will return kO2 in units of m/(1/48) - converts to fraction of day)


k.crusius <- function(wnd,wtr){
  U10 = wnd  #This function uses just the wind speed it is supplied. 
  k600 <- 0.228*U10^2.2+0.168 # units in cm h-1
  k600 <- k600*24/100 #units in m d-1
  ScO2 <- 1800.6 - (120.1 * wtr) + (3.7818 *wtr^2) - (0.047608 * wtr^3) # Schmidt number for O2 (Waninkhof, 1992)
  Sc600 <- ScO2 / 600
  n <- ifelse(U10 < 3.7,0.5,0.67) # n = 0.5 if U10 is less than 3.7 ms^-1, and 0.67 for U10 > 3.7 ms^-1
  kO2 <- k600 * (Sc600)^(-n) # gas exchange for O2 m d-1 
  kO2 <- kO2*(timeStep/1440) #change kO2 to units of m/(timeStep*min)
  return(kO2)
}
