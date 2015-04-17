#figure generation code for a metabolism example
library(LakeMetabolizer)

data.path = system.file('extdata/', package="LakeMetabolizer")
sp.data = load.all.data('sparkling', data.path)
ts.data = sp.data$data #pull out just the timeseries data


#calculate U10 and add it back onto the original
u10 = wind.scale(ts.data)
ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset


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

#Bring year and DOY together to get R datetime again
ols.res$datetime = ISOdate(ols.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24
mle.res$datetime = ISOdate(mle.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24
kalman.res$datetime = ISOdate(kalman.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24
bayes.res$datetime = ISOdate(bayes.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24
book.res$datetime = ISOdate(book.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24


png('~/fig_metab.png', res=300, width=1200, height=1800)
#Plot the metabolism results
par(mfrow=c(3,1), mar=c(0.1,4.4,1,1), las=1, oma=c(3,0,0,0))
plot(ols.res[,c(6,3)], type='o', col='black', ylim=c(-1,1), xaxt = 'n', ylab=expression(GPP~(mg~L^-1~day^-1)))
lines(mle.res[,c(6,3)], type='o', col='green')
lines(kalman.res[,c(6,3)], type='o', col='yellow')
lines(bayes.res[,c(6,3)], type='o', col='red')
lines(book.res[,c(6,3)], type='o', col='blue')
abline(0,0)
legend('bottomright', legend=c('metab.ols', 'metab.mle', 'metab.kalman', 'metab.bayesian', 'metab.bookkeep'), 
			 col=c('black','green','yellow','red','blue'), lty=1, horiz=FALSE)

plot(ols.res[,c(6,4)], type='o', col='black', ylim=c(-1,1), xaxt='n', ylab=expression(R~(mg~L^-1~day^-1)))
lines(mle.res[,c(6,4)], type='o', col='green')
lines(kalman.res[,c(6,4)], type='o', col='yellow')
lines(bayes.res[,c(6,4)], type='o', col='red')
lines(book.res[,c(6,4)], type='o', col='blue')
abline(0,0)

plot(ols.res[,c(6,5)], type='o', col='black', ylim=c(-0.5,0.5), ylab=expression(NEP~(mg~L^-1~day^-1)))
lines(mle.res[,c(6,5)], type='o', col='green')
lines(kalman.res[,c(6,5)], type='o', col='yellow')
lines(bayes.res[,c(6,5)], type='o', col='red')
lines(book.res[,c(6,5)], type='o', col='blue')
abline(0,0)


dev.off()

#summary stats
res = ols.res
cat('metab.ols GPP:', mean(res[,3]), ' R:', mean(res[,4]), 'NEP:', mean(res[,5]), '\n')
res = mle.res
cat('metab.mle GPP:', mean(res[,3]), ' R:', mean(res[,4]), 'NEP:', mean(res[,5]), '\n')
res = kalman.res
cat('metab.kalman GPP:', mean(res[,3]), ' R:', mean(res[,4]), 'NEP:', mean(res[,5]), '\n')
res = bayes.res
cat('metab.bayesian GPP:', mean(res[,3]), ' R:', mean(res[,4]), 'NEP:', mean(res[,5]), '\n')
res = book.res
cat('metab.bookkeep GPP:', mean(res[,3]), ' R:', mean(res[,4]), 'NEP:', mean(res[,5]), '\n')


