#---jread-usgs 2013-05-13---
getSchmidt	<-	function(temperature, gas){
	# temperature can be a number or a vector of numbers
	# gas must be a valid name. Code for more, but currently only in the supported data frame.
	
	# __Wanninkhof, Rik. "Relationship between wind speed and gas exchange over the ocean."__ 
	# __Journal of Geophysical Research: Oceans (1978–2012) 97.C5 (1992): 7373-7382.__
	
	range.t	<-	c(4,35)
	
	Schmidt	<-	data.frame(
		"He"=c(368,-16.75,0.374,-0.0036),
		"O2"=c(1568,-86.04,2.142,-0.0216),
		"CO2"=c(1742,-91.24,2.208,-0.0219),
		"CH4"=c(1824,-98.12,2.413,-0.0241),
		"SF6"=c(3255,-217.13,6.837,-0.0861),
		"N2O"=c(2105,-130.08,3.486,-0.0365),
		"Ar"=c(1799,-106.96,2.797,-0.0289),
		"N2"=c(1615,-92.15,2.349,-0.0240))
		
	if (!is.character(gas)){stop(paste('gas must be a character. was given as',gas))}
	if (length(gas)>1){stop("only one gas can be specified for this version")}
	if (!any(names(Schmidt)==gas)){stop(paste(gas,'not found in list of coded gasses'))}
	if (any(temperature < range.t[1] | temperature > range.t[2])){
		warning("temperature out of range")
	}
	A	<-	unlist(Schmidt[gas])[1]
	B	<-	unlist(Schmidt[gas])[2]
	C	<-	unlist(Schmidt[gas])[3]
	D	<-	unlist(Schmidt[gas])[4]

	Sc = as.numeric(A+B*temperature+C*temperature^2+D*temperature^3)

	return(Sc)
}