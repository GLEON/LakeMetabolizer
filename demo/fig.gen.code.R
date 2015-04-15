## Code to run models and generate figures for Manuscript
#  This will be displayed when "manuscript.code()" function is called
#  A great outreach method for the package

data.path = system.file('figdata', package="LakeMetabolizer")

sp.data = load.all.data('sparkling', data.path)

ts.data = sp.data$data #pull out just the timeseries data

plot(ts.data$datetime, ts.data$doobs_0.5, type='l', ylab='DO',xlab='date')

plot(ts.data$datetime, ts.data$par, type='l', ylab='PAR',xlab='date')

plot(ts.data$datetime, ts.data$wnd_2.0, type='l', ylab='PAR',xlab='date')

#calculate U10 and add it back onto the original
u10 = wind.scale(ts.data)
ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset

#show the scaled wind
lines(ts.data$datetime, ts.data$wnd_10, col='blue')

#Now calculate k600 using the Cole method
k600.cole = k.cole(ts.data)

ts.data = merge(ts.data, k600.cole)

kgas = k600.2.kGAS(ts.data)
ts.data = rmv.vars(merge(kgas, ts.data), 'k600')

o2.sat = o2.at.sat(ts.data[,c('datetime','wtr_0')])

ts.data = merge(o2.sat, ts.data)
z.mix = ts.meta.depths(get.vars(ts.data, 'wtr'), seasonal=TRUE)
names(z.mix) = c('datetime','z.mix', 'bottom')

#set z.mix to bottom of lake when undefined
z.mix[z.mix$z.mix <=0 | is.na(z.mix$z.mix), 'z.mix'] = 20 
ts.data = merge(ts.data, z.mix[,c('datetime','z.mix')])


#Plot the data
plot(ts.data$datetime, ts.data$z.mix, type='l', ylab='z.mix',xlab='date')



#OLS
ols.res = metab(ts.data, method='ols', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')
#write.csv(res, 'sp.metab.ols.csv', row.names=FALSE)

#MLE
mle.res = metab(ts.data, method='mle', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

#Kalman
kalman.res = metab(ts.data, method='kalman', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

#Bayesian
bayes.res = metab(ts.data, method='bayesian', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

#Bookkeep
ts.data[, names(ts.data)%in%'par'] = as.numeric(ts.data[, names(ts.data)%in%'par'] >=0)

book.res = metab(ts.data, method='bookkeep', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

#Plot the metabolism results
par(mfrow=c(2,1), mar=c(1,4,1,1), las=1)
plot(ols.res[,2:3], type='l', col='black', ylim=c(-1,1))
lines(mle.res[,2:3], type='l', col='green')
lines(kalman.res[,2:3], type='l', col='yellow')
lines(bayes.res[,2:3], type='l', col='red')
lines(book.res[,2:3], type='l', col='blue')

plot(ols.res[,c(2,4)], type='l', col='black', ylim=c(-1,1))
lines(mle.res[,c(2,4)], type='l', col='green')
lines(kalman.res[,c(2,4)], type='l', col='yellow')
lines(bayes.res[,c(2,4)], type='l', col='red')
lines(book.res[,c(2,4)], type='l', col='blue')

################################################################################
## Run trout bog
################################################################################
library(zoo)

data.path = system.file('extdata', package="LakeMetabolizer")

sp.data = load.all.data('troutbog', data.path)

ts.data = sp.data$data #pull out just the timeseries data

#Do a simple linear interpolation to fill the few missing temp values
wtr = get.vars(ts.data, 'wtr')
wtr_interp = data.frame(lapply(wtr[,-1], function(d)na.approx(d, na.rm=FALSE)))
wtr_interp$datetime = wtr$datetime

plot(ts.data$datetime, ts.data$doobs_0.5, type='l', ylab='DO',xlab='date')

plot(ts.data$datetime, ts.data$par, type='l', ylab='PAR',xlab='date')

plot(ts.data$datetime, ts.data$wnd_2.0, type='l', ylab='PAR',xlab='date')

#calculate U10 and add it back onto the original
u10 = wind.scale(ts.data)
ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset

#show the scaled wind
lines(ts.data$datetime, ts.data$wnd_10, col='blue')

#Now calculate k600 using the Cole method
k600.cole = k.cole(ts.data)

ts.data = merge(ts.data, k600.cole)

kgas = k600.2.kGAS(ts.data)
ts.data = rmv.vars(merge(kgas, ts.data), 'k600')

o2.sat = o2.at.sat(ts.data[,c('datetime','wtr_0')])

ts.data = merge(o2.sat, ts.data)
z.mix = ts.meta.depths(get.vars(ts.data, 'wtr'), seasonal=TRUE)
names(z.mix) = c('datetime','z.mix', 'bottom')

#set z.mix to bottom of lake when undefined
z.mix[z.mix$z.mix <=0 | is.na(z.mix$z.mix), 'z.mix'] = 20 
ts.data = merge(ts.data, z.mix[,c('datetime','z.mix')])


#Plot the data
plot(ts.data$datetime, ts.data$z.mix, type='l', ylab='z.mix',xlab='date')



#OLS
ols.res = metab(ts.data, method='ols', 
								wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')
#write.csv(res, 'sp.metab.ols.csv', row.names=FALSE)

#MLE
mle.res = metab(ts.data, method='mle', 
								wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

#Kalman
kalman.res = metab(ts.data, method='kalman', 
									 wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

#Bayesian
bayes.res = metab(ts.data, method='bayesian', 
									wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

#Bookkeep
ts.data[, names(ts.data)%in%'par'] = as.numeric(ts.data[, names(ts.data)%in%'par'] >=0)

book.res = metab(ts.data, method='bookkeep', 
								 wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

#Plot the metabolism results
par(mfrow=c(2,1), mar=c(1,4,1,1), las=1)
plot(ols.res[,2:3], type='l', col='black', ylim=c(-1,1))
lines(mle.res[,2:3], type='l', col='green')
lines(kalman.res[,2:3], type='l', col='yellow')
lines(bayes.res[,2:3], type='l', col='red')
lines(book.res[,2:3], type='l', col='blue')

plot(ols.res[,c(2,4)], type='l', col='black', ylim=c(-1,1))
lines(mle.res[,c(2,4)], type='l', col='green')
lines(kalman.res[,c(2,4)], type='l', col='yellow')
lines(bayes.res[,c(2,4)], type='l', col='red')
lines(book.res[,c(2,4)], type='l', col='blue')

################################################################################
## Run Mendota
################################################################################
library(zoo)
Sys.setenv(TZ='GMT')
data.path = system.file('extdata', package="LakeMetabolizer")

me.data = load.all.data('mendota', data.path)

ts.data = me.data$data #pull out just the timeseries data

ts_interp = data.frame(lapply(ts.data[,-1], function(d)na.approx(d, na.rm=FALSE)))
ts.data = cbind(ts.data$datetime, ts_interp)
names(ts.data)[1] = 'datetime'


#Do a simple linear interpolation to fill the few missing temp values
#wtr = get.vars(ts.data, 'wtr')

#wtr_interp = data.frame(lapply(wtr[,-1], function(d)na.approx(d, na.rm=FALSE)))
#wtr_interp$datetime = wtr$datetime

#ts.data = merge(wtr_interp, rmv.vars(ts.data,'wtr',ignore.offset=TRUE),by='datetime')

#calculate U10 and add it back onto the original
u10 = wind.scale(ts.data, wnd.z=2)
ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
ts.data = merge(ts.data, u10,by='datetime')                          #merge new u10 into big dataset

ts.data$k600 = k.cole.base(ts.data$wnd_10)
ts.data$k.gas = k600.2.kGAS.base(ts.data$k600, ts.data$wtr_0)

ts.data = rmv.vars(ts.data, 'k600')

o2.sat = o2.at.sat(ts.data[,c('datetime','wtr_0')])

ts.data = merge(o2.sat, ts.data,by='datetime')

ts.data$z.mix = 0.25

#OLS
res = metab(ts.data, method='ols', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

write.csv(res, 'me.metab.ols.csv', row.names=FALSE)

#MLE
res = metab(ts.data, method='mle', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

write.csv(res, 'me.metab.mle.csv', row.names=FALSE)

#Kalman
res = metab(ts.data, method='kalman', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

write.csv(res, 'me.metab.kalman.csv', row.names=FALSE)

#Bayesian
res = metab(ts.data, method='bayesian', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

write.csv(res, 'me.metab.bayesian.csv', row.names=FALSE)

#Bookkeep
ts.data[, names(ts.data)%in%'par'] = as.numeric(ts.data[, names(ts.data)%in%'par'] > 5)

res = metab(ts.data, method='bookkeep', 
            wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='par')

write.csv(res, 'me.metab.bookkeep.csv', row.names=FALSE)


################################################################################
# metab.bayes example
################################################################################
# ==================
# = Load Functions =
# ==================	
#source("scale.exp.wind.R")
#source("k.cole.R")
#source("k600.2.kGAS.R")
#source("getSchmidt.R")
#source("date2doy.R")
#source("o2.saturation.R")
# source("longestRun.R")

# ====================
# = Read in Raw Data =
# ====================
# wind data >> scale.exp.wind() >> k.cole() >> k600.2.kGAS() >>
wo <- read.table("../inst/extdata/troutbog.wnd", sep="\t", header=TRUE, colClasses=c("POSIXct","numeric"))
po <- read.table("../inst/extdata/troutbog.par", sep="\t", header=TRUE, colClasses=c("POSIXct","numeric"))
to <- read.table("../inst/extdata/troutbog.wtr", sep="\t", header=TRUE, colClasses=c("POSIXct",rep("numeric",10)))[,1:2]
doo <- read.table("../inst/extdata/troutbog.doobs", sep="\t", header=TRUE, colClasses=c("POSIXct","numeric"))

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



metab.bayes(irr=data0[,"PAR"], z.mix=rep(1, length(Kvec)), 
            do.sat=o2.at.sat(data0[,"Temp"], 716), wtr=data0[,'Temp'],
            k.gas=Kvec, do.obs=data0[,"DO"])


# dev.new()
# par(mfrow=c(3,2), mar=c(2,2,1.5,0.5), ps=10)
#nP <- ncol(test[[1]]$BUGSoutput$sims.matrix)
#parD <- c(ceiling(sqrt(nP)), floor(sqrt(nP)))
# R2jags::traceplot(test[[1]], mfrow=c(parD))

##dev.new()
#par(mfrow=c(3,2), mar=c(2,2,1.5,0.5), ps=10)
#coda::traceplot(as.mcmc(test[[1]]))

#dev.new()
#par(mfrow=c(3,2), mar=c(2,2,1.5,0.5), ps=10)
#coda::densplot(as.mcmc(test[[1]]))


################################################################################
# Example use of KF
################################################################################
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", repos=NULL, type="source")
# library("LakeMetabolizer")
# source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inProgress/ryanData.R")
# library("LakeMetabolizer")

KFans <- metab.kalman(do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], k.gas=data[,"k.gas"], z.mix=data[,"z.mix"], irr=data[,"irr"], wtr=data[,"wtr"])

# =========================================
# = Plot smoothed and orginal time series =
# =========================================
dev.new(height=6, width=3.5)
par(mfrow=c(2,1), mar=c(2,2.5,0.5,0.5), oma=c(0.5, 0, 0, 0), ps=9, mgp=c(1.5,0.2,0), tcl=-0.2, xpd=TRUE)
plot(data[,"do.obs"], type="l", lwd=5, col="gray", xlab="", ylab=bquote(DO))
legend("topleft", legend=c("obs", "smoothed"), col=c("gray","red"), lty=1, lwd=c(5, 1))
lines(KFans$smoothDO, col="red")
plot(diff(data[,"do.obs"]), type="l", lwd=5, col="gray", xlab="", ylab=bquote(italic(d)*DO))
mtext("time",side=1, line=-0.5, outer=TRUE)
lines(diff(KFans$smoothDO), col="red")


################################################################################
# misc code
################################################################################

## pDens.R code
# Density Plot Color Scheme
# Densities will be estimated and plottted for each column of val (rows are multiple observations of same variable)
pDen <- function(vals=NULL, mu=NULL, sig=NULL){
	if(is.null(vals) & (is.null(mu)|is.null(sig))) stop("Must provide vals or (mu and sig), but not both.")
	if(any(sig<=0)) stop("Need positive sig")
	
	N <- length(mu)
	
	if(is.null(vals)){
		vals0 <- rnorm(n=100*N, mean=mu, sd=sig)
		vals <- matrix(vals0, ncol=N, byrow=TRUE)	
	}
	limX <- range(vals)
	
	dens <- apply(vals, 2, function(x)density(x, from=limX[1], to=limX[2])[c("x","y")])
	
	dX <- lapply(dens, function(z)z$x)
	dY <- lapply(dens, function(z)z$y)
	limY <- range(dY)
	
	cLine <- rainbow(n=N)
	cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=35, maxColorValue=255)
	
	# Need to be sure that the smallest and largest densities are 0 so that the bottom border of polygons are at 0 line
	xF <- 0.01 * diff(limX) # a "factor" by which to extend the range of X
	xA <- limX + c(-1,1)*xF # "add" this "adjustment" to the start and end of the dX

	dev.new(width=5, height=3.5)
	par(mar=c(1.5, 1.5, 0.5, 0.5), ps=10, cex=1, mgp=c(1.5, 0.1, 0), tcl=0.25, las=1)
	plot(c(xA[1],dX[[1]],xA[2]), c(0,dY[[1]],0), type="l", col=cLine[1], xlab="", ylab="", xlim=limX, ylim=limY, lwd=2)
	polygon(c(xA[1],dX[[1]],xA[2]), c(0,dY[[1]],0), col=cFill[1], border=NA)
	for(i in 2:N){
		polygon(c(xA[1],dX[[1]],xA[2]), c(0,dY[[i]],0), col=cFill[i], border=cLine[i], lwd=2)
	}
}

# Example (that simulates fake data):
# set.seed(3)
# pDen(mu=runif(5, -100, 100), sig=runif(5, 1, 15))
