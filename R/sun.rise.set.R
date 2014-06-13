
#'@title determines if measurement was taken during the daytime
#'@description 
#'determines if measurement was taken during the daytime
#'
#'@usage
#'is.day(datetimes, lat)
#'@param datetimes Vector of dates as \code{POSIXct} or \code{POSIXlt} (see \code{\link{DateTimeClasses}}) format
#'@param lat Single latitude value of site. South should be negative, north positive
#'
#'@return a boolean vector of same length as \code{datetimes} 
#'
#'@keywords methods

#'@author
#'Luke A. Winslow
#'@seealso 
#'\link{is.night}
#'\link{sun.rise.set}
#'@export
is.day <- function(datetimes, lat){
  sr.ss <- sun.rise.set(datetimes, lat)
  
  is.daytime <- xor(sr.ss[,1] > datetimes, sr.ss[,2] > datetimes)
  return(is.daytime)
}

#'@title determines if measurement was taken during the night
#'@description 
#'determines if measurement was taken during the nighttime 
#'
#'@usage
#'is.night(datetimes, lat)
#'@param datetimes Vector of dates as \code{POSIXct} or \code{POSIXlt} (see \code{\link{DateTimeClasses}}) format
#'@param lat Single latitude value of site. South should be negative, north positive
#'
#'@return a boolean vector of same length as \code{datetimes} 
#'
#'@keywords methods

#'@author
#'Luke A. Winslow
#'@seealso 
#'\link{is.day}
#'\link{sun.rise.set}
#'@export
is.night <- function(datetimes, lat){
  return(!is.day(datetimes, lat))
}



sun.rise.set <- function(datetimes, lat){
	#SUNRISESET Calculates the time of sunrise and sunset
	#
	# Input:
	#   lat - The latitude for the desired location on earth.
	#   day - date
	#
	# takes lattitude (deg south should be negative, north positive) and a
	# matlab datenum (days since year 0) and returns the time of sunrise as a
	# day fraction (e.g., 0.5 would be noon, 0 would be midnight)

	doy <- as.POSIXlt(datetimes)$yday+1 # POSIX functions treat January 1 as day of year 0, so add 1 to compensate


	#TODO: Add leap-year fix
	dayAngle <- 2*pi*(doy-1)/365;


	degToRad <- 2*pi/360;
	radToDeg <- 180/pi;

	#Declination of the sun "delta" (radians). Iqbal 1983 Eq. 1.3.1
	dec <- 0.006918 - 0.399912*cos(dayAngle) + 0.070257*sin(dayAngle) - 0.006758*cos(2*dayAngle) +  0.000907*sin(2*dayAngle) - 0.002697*cos(3*dayAngle) + 0.00148*sin(3*dayAngle);

	#Sunrise hour angle "omega" (degrees). Iqbal 1983 Eq. 1.5.4
	latRad <- lat*degToRad;
	sunriseHourAngle <- acos(-tan(latRad)*tan(dec))*radToDeg;

	#If we don't have a sunrise, then sunriseHourAngle is imaginary, replace with NaN
	sunriseHourAngle[is.complex(sunriseHourAngle)] = NA;

	#Sunrise and sunset times (decimal hours, relative to solar time) Iqbal 1983 Ex. 1.5.1
	sr <- 12 - sunriseHourAngle/15;
	ss <- 12 + sunriseHourAngle/15;

	#convert to seconds into day
	rise <- trunc(datetimes, 'day') + sr*60*60
	set <- trunc(datetimes, 'day') + ss*60*60

	##Note, this does weird things. It *is* a matrix, but it doesn't print like one because it is viewed 
	# as POSIXct. I will leave it this way for now, though if someone knows how to get it to show up as
	# a matrix *and* a POSIXct value, that would be super cool.

	return(as.POSIXct(matrix(c(rise, set), ncol=2), origin='1970-01-01'))

}
