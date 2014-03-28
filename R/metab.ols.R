metab.ols = function(do.obs, do.sat, irr, k.gas, z.mix){
  
  n.obs = length(do.obs)
  
  do.diff = diff(do.obs)
  
  #basically the average of the flux at time T_0 and T_1
  #flux = -0.5 * k * (D1 + D2 - 2*S)
  #
  # can be re-arranged
  # flux = 0.5 * [k * (D1-S) + k * (d2-s)]
  # Average of each inst flux
  
  #Hmm, this will need to be fixed
  inst_flux = k.gas * (do.sat - do.obs)  # positive is into the lake
  
  flux = apply(matrix(c(inst_flux[1:(n.obs-1)], inst_flux[2:(n.obs)]), ncol=2, byrow=TRUE), 1, mean)
  
  noflux.do.diff = do.diff - flux/z.mix
  
  mod = lm(noflux.do.diff ~ irr)
  
  rho = mod[[1]][1]
  iota = mod[[1]][2]
  gpp = iota + sum(irr)
  
  return(list(rho=rho, gpp=gpp, iota=iota))
}