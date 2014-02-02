# ---Author: Hilary Dugan, 2013-10-20 --- 
# Last update: 2014-02-01

# Uses power function from:
# -- References 
#CRUSIUS, JOHN, AND RIK WANNINKHOF. 2003
#Gas transfer velocities measured at low wind speed over a lak.
#Limnology and Oceanography. 48(3): 1010:1017.

# INPUTS
# wnd: wind value in metres/second
# wtr: water temperature at dissolved oxygen sensor depth 
# in degrees celsius 
# timeStep: time interval between measurements, in minutes
# method: Either "linear", "bilinear", or "power"
# as defined in Crusius et al., (2003). Default is "power"

# OUTPUT: returns the gas exchange velocity for O2 in units 
# of m/(timeStep).
# If timeStep = 60, then k02 measured in cm/60 min)

k.crusius <- function(wnd,wtr,timeStep,method){
  if (nargs()==3) {
    method = "power"
  }
    
  U10 = wnd  #This function uses just the wind speed it is supplied. 
  
  if (method=="constant"){
    k600 <- ifelse(U10<3.7,1,5.14*U10-17.9)
  } else if (method=="bilinear") {
    k600 <- ifelse(U10<3.7,0.72*U10,4.33*U10-13.3)
  } else if (method=="power") {
    k600 <- 0.228*U10^2.2+0.168 # units in cm h-1
  }
  
  k600 <- k600*24/100 #units in m d-1
  # Schmidt number for O2 (Waninkhof, 1992)
  ScO2 <- 1800.6 - (120.1 * wtr) + (3.7818 *wtr^2) - (0.047608 * wtr^3) 
  Sc600 <- ScO2 / 600
  # n=0.5 if U10 < 3.7 ms^-1, and 0.67 for U10 > 3.7 ms^-1
  n <- ifelse(U10 < 3.7,0.5,0.67) 
  kO2 <- k600 * (Sc600)^(-n) # gas exchange for O2 m d-1 
  kO2 <- kO2*(timeStep/1440) #change kO2 to units of m/(timeStep*min)
  return(kO2)
}
