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

k.crusius <- function(wnd,method='power'){

  U10 = wnd  #This function uses just the wind speed it is supplied. 
  
  if (method=="constant"){
    k600 <- ifelse(U10<3.7,1,5.14*U10-17.9)
  } else if (method=="bilinear") {
    k600 <- ifelse(U10<3.7,0.72*U10,4.33*U10-13.3)
  } else if (method=="power") {
    k600 <- 0.228*U10^2.2+0.168 # units in cm h-1
  }
  
  k600 <- k600*24/100 #units in m d-1
  return(k600)
}
