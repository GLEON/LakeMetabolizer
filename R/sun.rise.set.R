
is.day = function(lat, datetime){
  sr.ss = sun.rise.set(lat, datetime)
  
  is.daytime = xor(sr.ss[,1] > datetime, sr.ss[,2] > datetime)
  return(is.daytime)
}


is.night = function(lat, datetime){
  return(!is.daytime(lat, datetime))
}



sun.rise.set = function(lat, datetimes){
#SUNRISESET Calculates the time of sunrise and sunset
#
# Input:
#   lat - The latitude for the desired location on earth.
#   day - Either Matlab datenum or day-of-year (Jan 1st is 1)
#
# takes lattitude (deg south should be negative, north positive) and a
# matlab datenum (days since year 0) and returns the time of sunrise as a
# day fraction (e.g., 0.5 would be noon, 0 would be midnight)

doy = as.POSIXlt(datetimes)$yday


#TODO: Add leap-year fix
dayAngle = 2*pi*(doy-1)/365;


degToRad = 2*pi/360;
radToDeg = 180/pi;

#Declination of the sun "delta" (radians). Iqbal 1983 Eq. 1.3.1
dec = 0.006918 - 0.399912*cos(dayAngle) + 0.070257*sin(dayAngle) - 0.006758*cos(2*dayAngle) +  0.000907*sin(2*dayAngle) - 0.002697*cos(3*dayAngle) + 0.00148*sin(3*dayAngle);

#Sunrise hour angle "omega" (degrees). Iqbal 1983 Eq. 1.5.4
latRad = lat*degToRad;
sunriseHourAngle = acos(-tan(latRad)*tan(dec))*radToDeg;

#If we don't have a sunrise, then sunriseHourAngle is imaginary, replace with NaN
sunriseHourAngle[is.complex(sunriseHourAngle)] = NA;

#Sunrise and sunset times (decimal hours, relative to solar time) Iqbal 1983 Ex. 1.5.1
sr = 12 - sunriseHourAngle/15;
ss = 12 + sunriseHourAngle/15;

#convert to seconds into day
rise = trunc(datetimes, 'day') + sr*60*60
set = trunc(datetimes, 'day') + ss*60*60

return(as.POSIXct(matrix(c(rise, set), ncol=2), origin='1970-01-01'))

}