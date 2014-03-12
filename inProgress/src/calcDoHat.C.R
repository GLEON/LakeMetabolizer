##This needs to be converted to just a straight C function

library(inline)
calcDoHat_C = cxxfunction(signature(iota_R = "numeric", rho_R = "numeric", 
                  doInit_R = "numeric", irr_R = "numeric", doSat_R = "numeric", 
                  zMix_R = "numeric", kO2_R = "numeric"), body = '

  

  double iota = as<double>(iota_R);
  double doInit = as<double>(doInit_R);
  double rho = as<double>(rho_R);
  
  //printf("Iota:%f  Rho:%f  doInit:%f\\n", iota, rho, doInit);

  Rcpp::NumericVector irr(irr_R);
  Rcpp::NumericVector doSat(doSat_R);
  Rcpp::NumericVector zMix(zMix_R);
  Rcpp::NumericVector kO2(kO2_R);
  
  int nObs = irr.size();
  //#Set up output
  Rcpp::NumericVector DOHat(nObs);
  Rcpp::NumericVector atmFlux(nObs);
  //#Initialize DOHat
  DOHat[0] = doInit;
  
  //#Calculate atmFlux and predicted DO for each time point
  //#Fluxes out of lake have negative sign
  for( int i=0; i<(nObs-1); i++){
    atmFlux[i] =   -kO2[i] * (DOHat[i] - doSat[i]) / zMix[i];
    DOHat[i+1] = 	DOHat[i] + iota * irr[i] - rho + atmFlux[i];
  }
  
  return DOHat;


', plugin="Rcpp")