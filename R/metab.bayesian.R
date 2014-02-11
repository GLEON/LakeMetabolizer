# library("rjags")
library("R2jags")
library("R2WinBUGS")

# ==================
# = Load Functions =
# ==================	
source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/R/scale.exp.wind.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/R/k.cole.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/R/k600.2.kGAS.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/R/getSchmidt.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/R/date2doy.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/R/o2.saturation.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/R/longestRun.R")

# ====================
# = Read in Raw Data =
# ====================
# wind data >> scale.exp.wind() >> k.cole() >> k600.2.kGAS() >>
wo <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inst/extdata/troutbog.wnd", sep="\t", header=TRUE, colClasses=c("POSIXct","numeric"))
po <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inst/extdata/troutbog.par", sep="\t", header=TRUE, colClasses=c("POSIXct","numeric"))
to <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inst/extdata/troutbog.wtr", sep="\t", header=TRUE, colClasses=c("POSIXct",rep("numeric",10)))[,1:2]
doo <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inst/extdata/troutbog.doobs", sep="\t", header=TRUE, colClasses=c("POSIXct","numeric"))

# ===========================
# = Combine & Organize Data =
# ===========================
#merge
d1 <- merge(to, doo, all=TRUE)
d2 <- merge(wo, po, all=TRUE)
d3 <- merge(d1,d2, all=TRUE)

#convert to DoY format (not actually needed in this case, but conventient to have)
d3[,1] <- date2doy(d3[,1])

#subset to the portion of the data set with the most consecutive observations (not necessary)
data0 <- d3[longestRun(d3),]
names(data0) <- c("DoY", "Temp", "DO", "Wind", "PAR") # rename columns while still data frame
data0 <- as.matrix(data0) # convert to matrix
row.names(data0) <- NULL # remove row names left over from d3
data0 <- data0[data0[,"DoY"]>=319&data0[,"DoY"]<320,] # subset for fast trial run

Freq <- median(diff(data0[,"DoY"])) # determine the sampling frequency; i have a function for mode if we are worried about it

wind <- scale.exp.wind(data0[,"Wind"], 2) # convert wind
Kvec <- k600.2.kGAS(k.cole(wind)*Freq, data0[,"Temp"], "O2") # calculate K for relevant sampling frequency


# =========================
# = Define the Jags Model =
# =========================
j2Mod <- function(){
	# model process
	for(i in 2:N){	
		Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
		K[i] ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distributin on K
		Uk[i] <- (K[i]*(satO[i-1] - a[i-1]))/Zmix[i-1] # exchange
		aHat[i] <- a[i-1] + U[i,]%*%C + Uk[i] #U[i,]%*%C #+ Uk # the process
		a[i] ~ dnorm(aHat[i], tauW) # true values have a mean equivalent to estimated values, but accompanied by process error (process precision is tauW)
	}
	
	# K <- kP[1:N,1] #distributin on K
	
	# Starting values
	# a[1] <- Y[1] # set the first true value equal to first observed value
	a[1] <- a0
	C[1] ~ dnorm(cP[1,1], 1/cP[1,2])
	C[2] ~ dnorm(cP[2,1], 1/cP[2,2])
	#Priors on regression coefficients
	# C[1,1] ~ dnorm(cP[1,1], 1/cP[1,2]) # GPP coefficient
	# C[2,1] ~ dnorm(cP[2,1], 1/cP[2,2]) # R coefficient
	# C[3,1] <- 1 # dummy 1 for exchange
	
	#Prior on errors
	tauV ~ dgamma(1.0E-3, 1.0E-3)
	tauW ~ dgamma(1.0E-3, 1.0E-3)
	sigmaV <- 1/sqrt(tauV)
	sigmaW <- 1/sqrt(tauW)	
}
write.model(j2Mod, "j2Mod")

# ===========================================
# = Define objects to be used in jags model =
# ===========================================
#Supply elements of U (PAR, log(temp), 0's)
U <- matrix(NA, nrow=nrow(data0), ncol=2)
U[,1] <- data0[,"PAR"] # filler PAR Values
U[,2] <- log(data0[,"Temp"]) # filler log(temp) values

#Supply kP (for K), and cP (for C)
# ** I am making up these values for now
kP <- matrix(NA, nrow=length(Kvec), ncol=2)
kP[,1] <- Kvec # means for K
kP[,2] <- sum(Kvec)/length(Kvec)*0.1 # variances for K
cP <- matrix(NA, nrow=2, ncol=2)
cP[1,1] <- 0 # prior mean of GPP coefficient (C[1,1]*PAR=GPP)
cP[1,2] <- 1E5 # prior variance of GPP coefficient
cP[2,1] <- 0 # prior mean of R coefficient (C[2,1]*log(Temp)=R)
cP[2,2] <- 1E5 # prior variance of R coefficient

#Supply observed DO
Y <- data0[,"DO"]

#Supply saturated DO
satO <- o2.at.sat(data0[,"Temp"], 716)

#Supply Zmix
Zmix <- rep(1, length(Kvec)) # filler values

# Put in final format supplied to jags
data <- list(Y=Y, N=length(Y), U=U, kP=kP, cP=cP, satO=satO, a0=Y[1], Zmix=Zmix)
params <- c("C", "sigmaV", "sigmaW")


# ================================
# = Supply Data and run bayesFit =
# ================================
bayesFit <- function(data, params, tend="median", ...){ #function that writes jags model, traces params, supplies data, etc

	jags.m <- jags(data, NULL, parameters.to.save=params, "j2Mod", n.chains=3, n.iter=1E4, n.burnin=5E3)
	
	tF <- function(x, tend){ # tendency function
		switch(tend,
			median=median(x), #median
			mean=mean(x), #mean
			mode1 = unique(x)[which.max(tabulate(match(x, unique(x))))], # mode --- most frequently observed value
			mode2 = { # mode --- based on the highest peak of the posterior probability (from density plot)
				xd <- density(x)
				xd$x[which.max(xd$y)]
			}
			)
	}
	# medSim <- matrix(apply(jags2.m$BUGSoutput$sims.matrix, 2, median)[-(1)], nrow=115, ncol=3)
	simOut <- jags.m$BUGSoutput$sims.matrix
	ctSim <- apply(simOut, 2, tF, tend) # the central tendency metric
	sdSim <- apply(simOut, 2, sd)
	
	# apply(jags.m$BUGSoutput$sims.matrix, 2, tF, tend="median")
	# apply(jags.m$BUGSoutput$sims.matrix, 2, tF, tend="mean")
	# apply(jags.m$BUGSoutput$sims.matrix, 2, tF, tend="mode1")
	# apply(jags.m$BUGSoutput$sims.matrix, 2, tF, tend="mode2")
	
	#Figure out the order of the sims.matrix columns ...
	GPP <- sum(ctSim[1]*data$U[,1]) # gpp coef * par, then sum
	R <- sum(ctSim[2]*data$U[,2]) # r coef * log(temp), then sum
	
	GPPsd <- sqrt(sum(sdSim[1]^2*data$U[,1]^2))
	Rsd <- sqrt(sum(sdSim[2]^2*data$U[,2]^2))
	
	return(list(jags.m, matrix(c(GPP,GPPsd,R,Rsd), nrow=2, dimnames=list(c("mu", "sd"), c("GPP", "R"))))) # need to clean up format, and maybe include a return of the sd's of the estimates
}


test <- bayesFit(data, params)
# dev.new()
# par(mfrow=c(3,2), mar=c(2,2,1.5,0.5), ps=10)
nP <- ncol(test[[1]]$BUGSoutput$sims.matrix)
parD <- c(ceiling(sqrt(nP)), floor(sqrt(nP)))
# R2jags::traceplot(test[[1]], mfrow=c(parD))

dev.new()
par(mfrow=c(3,2), mar=c(2,2,1.5,0.5), ps=10)
coda::traceplot(as.mcmc(test[[1]]))

dev.new()
par(mfrow=c(3,2), mar=c(2,2,1.5,0.5), ps=10)
coda::densplot(as.mcmc(test[[1]]))


