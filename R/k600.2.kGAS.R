#---jread-usgs 2014-01-04---

# INPUTS:
# k600: gas transfer velocity in units Time/Distance
# temperature: water temperature in °C
# gas: name of a gas (e.g., "CO2" or "O2")

# OUTPUT:
# kGAS: gas transfer velocity for the specified gas (in same units as k600 input)
k600.2.kGAS.base	<-	function(k600,temperature,gas){
	
	n	<-	0.5
	schmidt	<-	getSchmidt(temperature,gas)
	Sc600	<-	schmidt/600
	
	kGAS	<-	k600*(Sc600^-n)
	return(kGAS)
}


k600.2.kGAS = function(k600, ts.data, gas="O2"){
  all.data = k600
  temperature = get.vars(ts.data,'wtr')[,1:2]
  
  all.data = merge(all.data, temperature, by='datetime')
  
  kGAS = data.frame(datetime=all.data$datetime)
  
  kGAS$kGAS = k600.2.kGAS(all.data[,2], all.data[,3], gas)
  
  return(kGAS)
}
