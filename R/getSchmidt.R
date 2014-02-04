#---jread-usgs 2013-05-13---
getSchmidt	<-	function(temperature, gas){
	# temperature can be a number or a vector of numbers
	# gas must be a valid name. Code for more, but currently only in the supported data frame.
	
	# __Wanninkhof, Rik. "Relationship between wind speed and gas exchange over the ocean."__ 
	# __Journal of Geophysical Research: Oceans (1978–2012) 97.C5 (1992): 7373-7382.__
	
	Schmidt	<-	data.frame("He"=c(377.09,19.154,0.50137,0.005669),
		"O2"=c(1800.6,120.1,3.7818,0.047608),
		"CO2"=c(1911.1,118.11,3.4527,0.04132),
		"CH4"=c(1897.8,114.28,3.2902,0.039061),
		"SF6"=c(3255.3,217.13,6.8370,0.086070),
		"N2O"=c(2055.6,137.11,4.3173,0.05435),
		"Ar"=c(1759.7,117.37,3.6959,0.046527))
		
	if (!is.character(gas)){stop(paste('gas must be a character. was given as',gas))}
	if (length(gas)>1){stop("only one gas can be specified for this version")}
	if (!any(names(Schmidt)==gas)){stop(paste(gas,'not found in list of coded gasses'))}

	A	<-	unlist(Schmidt[gas])[1]
	B	<-	unlist(Schmidt[gas])[2]
	C	<-	unlist(Schmidt[gas])[3]
	D	<-	unlist(Schmidt[gas])[4]

	Sc = as.numeric(A-B*temperature+C*temperature^2-D*temperature^3)

	return(Sc)
}