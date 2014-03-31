metab.ols <- function(do.obs, do.sat, irr, k.gas, z.mix){

	n.obs <- length(do.obs)

	do.diff <- diff(do.obs)

	#basically the average of the flux at time T_0 and T_1
	#flux <- -0.5 * k * (D1 + D2 - 2*S)
	#
	# can be re-arranged
	# flux <- 0.5 * [k * (D1-S) + k * (d2-s)]
	# Average of each inst flux

	#Hmm, this will need to be fixed
	inst_flux <- k.gas * (do.sat - do.obs)  # positive is into the lake

	# flux = apply(matrix(c(inst_flux[1:(n.obs-1)], inst_flux[2:(n.obs)]), ncol=2, byrow=TRUE), 1, mean)
	# above could be simpler
	# note that inst_flux[2:(n.obs)] == inst_flux[-1] 
	# note that rather than applying mean row-wise, you can just add the two vectors, then divide by two
	# this none-apply() approach takes advantage of R's 'vectorization' – should be faster unless I'm messing something up
	flux <- (inst_flux[1:(n.obs-1)] + inst_flux[-1])/2

	noflux.do.diff <- do.diff - flux/z.mix

	mod <- lm(noflux.do.diff ~ irr) # note that if we decied to use resp=rho*log(Temp), you would do lm(do~irr+lntemp-1) (mind the -1)

	# also note that this model has different structure than Bayes, mle, Kalman 
	# these have resp as rho*log(Temp), rather than just an intercept – no prob, just need to pick one, easy to change once we get there
	rho <- mod[[1]][1] * length(irr) # Resp <- Beta0*nSteps, right? B/c the intercept is the average resp b/w steps? ~RDB
	iota <- mod[[1]][2]
	gpp <- iota + sum(irr)
	nep <- sum(fitted(mod))
	#other ways to get nep:
	# nep = 
	# also note that NEP is gpp+rho (rho is negative by this convention, which is consistent w/ Kalman, Bayes, mle – unsure of BK)

	# return(list(rho=rho, gpp=gpp, iota=iota))
	results <- data.frame("GPP"=gpp, "R"=rho, "NEP"=nep) 
	attr(results, "lm.mod") <- mod
	return(results) # i've seen dates about whether return() should be used, or if 'results' should just be the last line. 
	# Personally, I like return() b/c it makes debugging easier. I don't know if it matters beyond that.
	# You can use return() to stop a function and have it return NULL (can't do that w/o call to return())
		# e.g.
		# test <- function () i = 1+1
		# is.null(test()) # FALSE
		# test2 <- function () {i = 1+1; return()}
		# is.null(test()) # TRUE
		# see SO post http://stackoverflow.com/questions/11738823/explicitly-calling-return-in-a-function-or-not
		# I think using return() also has aesthetic benefits, 
		# so arguments for use of return() may be similar to the discussion of '=' and '<-' (please use '<-'!! <3)
	
	# notes on using attr() to simplify the output of our functions (extended from convo b/w LAW and RDB):
	# this makes it so that the results appear as the simple data frame
	# atrributes(results) shows the names of the attributes
	# attr(results, "lm.mod") will access the lm object
	# class(results) is data.frame
	# rbind(results, data.frame("NEP"=1, "GPP"=5, "R"=-1)) will provide the reasonable data frame 
	# dim(results) shows the data.frame dimensions
	# ... point i'm making is that it look and acts like a data frame, until you specifically attempt to get at lm.mod.
	# appears simple, has the functionality that is consistent w/ the simple appearence, but hides a lot of info if that info is desired
}