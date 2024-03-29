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
z.mix[z.mix$z.mix <=0 | is.na(z.mix$z.mix), 'z.mix'] = sp.data$metadata$maxdepth
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
# add column for day (1) and night (0)
ts.data$day_night = as.integer(is.day(datetimes=ts.data$datetime, lat=sp.data$metadata$latitude))

book.res = metab(ts.data, method='bookkeep',
								 wtr.name='wtr_0.5', do.obs.name='doobs_0.5', irr.name='day_night')

#Bring year and DOY together to get R datetime again
ols.res$datetime = ISOdate(ols.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24
mle.res$datetime = ISOdate(mle.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24
kalman.res$datetime = ISOdate(kalman.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24
bayes.res$datetime = ISOdate(bayes.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24
book.res$datetime = ISOdate(book.res$year, 1, 1) + book.res$doy*3600*24 - 3600*24



add_axes <- function(xlim, ylim, ylabel = pretty(ylim,10), panel.txt, no.x=TRUE){
	prc_x = 0.14 # position for text relative to axes
	prc_y = 0.07
	tick_len <- 0.15
	ext_x <- c(xlim[1]-86400, pretty(xlim,3), xlim[2]+86400)
	ext_y <- c(ylim[1]-10, pretty(ylim,10), ylim[2]+10)
	ylab <- c("",ylabel,"")
	if (is.na(ylabel[1])) ylab = NA
	#if(no.x) ext_x=NA
	if(no.x){
		axis(side = 1, at = ext_x , labels = FALSE, tcl = tick_len)
	}else{
		axis(side = 1, at = ext_x , labels = strftime(ext_x,'%d %b'), tcl = tick_len)
	}
	axis(side = 2, at = ext_y, labels = ylab, tcl = tick_len)
	axis(side = 3, at = ext_x, labels = NA, tcl = tick_len)
	axis(side = 4, at = ext_y, labels = NA, tcl = tick_len)
	x_txt <- (xlim[2] - xlim[1])*prc_x+xlim[1]
	y_txt <- ylim[2]-(ylim[2] - ylim[1])*prc_y
	text(x = x_txt, y_txt,labels = panel.txt)
}

add_legend <- function(models, xlim, ylim, prc_x = 0.26, prc_y = 0.06){

	y_strt <- ylim[2]-(ylim[2] - ylim[1])*prc_y
	y_spc <- (ylim[2] - ylim[1])*0.05
	x_len <- (xlim[2] - xlim[1])*0.16
	x <- c((xlim[2] - xlim[1])*prc_x+xlim[1], (xlim[2] - xlim[1])*prc_x+xlim[1] + x_len)

	for (i in 1:length(models)){
		y = y_strt-(i-1)*y_spc
		lines(x, c(y,y),
					col =models[[i]]$col,
					lty = models[[i]]$lty,
					lwd = models[[i]]$lwd)
		text(x[2],y, models[[i]]$name, pos = 4, cex = 0.65)
	}
}

add_models <- function(models, cols){
	.empty = sapply(X = models, FUN = function(x) {
		lines(x$data[,cols], col=x$col, lty = x$lty, lwd = x$lwd, type='o')
	})
}

cols <- c("#1b9e77", "#d95f02", "black", "#e7298a", "DodgerBlue", "#e6ab02", "grey50")

models <- list(
	list('name'="ols", data = ols.res, col = cols[1], lty = 6, lwd = 1.7),
	list('name'="mle", data = mle.res, col = cols[2], lty = 1, lwd = 1.2),
	list('name'="kalman", data = kalman.res, col = cols[3], lty = 1, lwd = 1.1),
	list('name'="bayesian", data = bayes.res, col = cols[4], lty = 1, lwd = 1.2),
	list('name'="bookkeep", data = book.res, col = cols[5], lty = 6, lwd = 1.7))


width = 3.37 # single column width for journal
night_col = 'grey90'
height = 5
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels
ylim=range(c(models[[1]]$data[3:5],models[[2]]$data[3:5],models[[3]]$data[3:5],models[[4]]$data[3:5],models[[5]]$data[3:5]))
xlim = as.POSIXct(c('2009-07-01 16:00', '2009-07-11'))

default_par = par(no.readonly = TRUE)
#Create plot and save in temporary directory
png(file.path(tempdir(), 'fig_metab.png'), res=300, width=width, height=height, units = 'in')

#layout(matrix(c(rep(1,10),rep(2,9)),ncol=1)) # 55% on the left panel
par(mai=c(b_mar,l_mar,t_mar,0), omi = c(0.1,0,0,r_mar),xpd=FALSE,
		mgp = c(1.15,.05,0), mfrow=c(3,1))

#Plot the metabolism results
plot(kalman.res[,c(6,3)], type='o', col=cols[1], ylim=c(-1,1.2), xaxt = 'n', ylab=expression(GPP~(O[2]~'in'~mg~L^-1~day^-1)), xlab='', axes=FALSE)
abline(0, 0, col=rgb(0,0,0,0.5))
add_models(models, c(6,3))
add_axes(xlim, ylim, panel.txt='a)')

plot(kalman.res[,c(6,4)], type='o', col=cols[1], ylim=c(-1,1), xaxt='n', ylab=expression(R~(O[2]~'in'~mg~L^-1~day^-1)), xlab='', axes=FALSE)
abline(0,0, col=rgb(0,0,0,0.5))
add_models(models, c(6,4))
add_axes(xlim, ylim, panel.txt='b)')

plot(ols.res[,c(6,5)], type='o', col=cols[1], ylim=c(-0.5,0.6), ylab=expression(NEP~(O[2]~'in'~mg~L^-1~day^-1)), xlab='', axes=FALSE)
abline(0,0, col=rgb(0,0,0,0.5))
add_models(models, c(6,5))
add_axes(xlim, ylim, panel.txt='c)', no.x=FALSE)

add_legend(models, xlim, ylim)

dev.off()
par(default_par)

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


#impossible values removed
cat('Averages after removal of impossible GPP and R values\n')
res = ols.res
cat('metab.ols GPP:', mean(res[res[,3]>0,3]), ' R:', mean(res[res[,4]<0,4]), 'NEP:', mean(res[res[,3]>0,3]) + mean(res[res[,4]<0,4]), '\n')
res = mle.res
cat('metab.mle GPP:', mean(res[res[,3]>0,3]), ' R:', mean(res[res[,4]<0,4]), 'NEP:',mean(res[res[,3]>0,3]) + mean(res[res[,4]<0,4]), '\n')
res = kalman.res
cat('metab.kalman GPP:', mean(res[res[,3]>0,3]), ' R:', mean(res[res[,4]<0,4]), 'NEP:', mean(res[res[,3]>0,3]) + mean(res[res[,4]<0,4]), '\n')
res = bayes.res
cat('metab.bayesian GPP:', mean(res[res[,3]>0,3]), ' R:', mean(res[res[,4]<0,4]), 'NEP:', mean(res[res[,3]>0,3]) + mean(res[res[,4]<0,4]), '\n')
res = book.res
cat('metab.bookkeep GPP:', mean(res[res[,3]>0,3]), ' R:', mean(res[res[,4]<0,4]), 'NEP:', mean(res[res[,3]>0,3]) + mean(res[res[,4]<0,4]), '\n')





