
# Generates a figure showing the results of the different k600 models
library(LakeMetabolizer)
data.path = system.file('extdata', package="LakeMetabolizer")

sp.data = load.all.data('sparkling', data.path)

ts.data = sp.data$data #pull out just the timeseries data

#calculate U10 and add it back onto the original

u10 = wind.scale(ts.data)
ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset

#Calculate k600 using the k.cole wind-based model
k600_cole = k.cole(ts.data)

k600_crusius = k.crusius(ts.data)

kd        = sp.data$metadata$averagekd
wnd.z      = 10   #because we converted to u10
atm.press  = 1018
lat       = sp.data$metadata$latitude
lake.area = sp.data$metadata$lakearea

#for k.read and k.macIntyre, we need LW_net.
#Calculate from the observations we have available.
lwnet = calc.lw.net(ts.data, lat, atm.press)
ts.data = merge(ts.data, lwnet)

k600_read = k.read(ts.data, wnd.z=wnd.z, Kd=kd, atm.press=atm.press,
									 lat=lat, lake.area=lake.area)

k600_soloviev = k.read.soloviev(ts.data, wnd.z=wnd.z, Kd=kd,
																atm.press=atm.press, lat=lat, lake.area=lake.area)

k600_macIntyre = k.macIntyre(ts.data, wnd.z=wnd.z, Kd=kd, atm.press=atm.press)

k600_heiskanen = k.heiskanen(ts.data, wnd.z, kd, atm.press)


k600_vachon = k.vachon(ts.data, lake.area)


#Create plot and save in user home directory (on Windows, Documents folder, on Mac, Home folder)
png('~/k600_figure.png', res=300, width=2400, height=1200)

#rise.set = sun.rise.set(xlim[1], lat)

cols = rainbow(7)
par(mar=c(5.1,4.3,4.1,2.1), mfrow=c(1,2), xpd=FALSE)
xlim = as.POSIXct(c('2009-07-03', '2009-07-04', '2009-07-09', '2009-07-10'))
plot(k600_macIntyre, type='l', xlab='', col=cols[1], ylim=c(0,10), ylab=expression(k600~(m~day^-1)), xlim=xlim[1:2])
lines(k600_cole, col=cols[2])
lines(k600_crusius, col=cols[3])
lines(k600_read, col=cols[4])
lines(k600_soloviev, col=cols[5])
lines(k600_heiskanen, col=cols[6])
lines(k600_vachon, col=cols[7])


plot(k600_macIntyre, type='l', xlab='', col=cols[1], ylim=c(0,10), ylab=expression(k600~(m~day^-1)), xlim=xlim[3:4])
lines(k600_cole, col=cols[2])
lines(k600_crusius, col=cols[3])
lines(k600_read, col=cols[4])
lines(k600_soloviev, col=cols[5])
lines(k600_heiskanen, col=cols[6])
lines(k600_vachon, col=cols[7])

#par(xpd=TRUE)

legend('topleft', legend=c('k.macIntyre', 'k.cole', 'k.crusius', 'k.read', 'k.read.soloviev', 'k.heiskanen', 'k.vachon'), 
			 lty=1, col=cols, horiz=FALSE, inset=c(-2,-0.1), ncol=4, xpd=TRUE)
dev.off()

#quick summary of k600 values for results section

cat('Cole:', mean(k600_cole[,2]), '\n')
cat('Crusius:', mean(k600_crusius[,2]), '\n')
cat('MacIntyre:', mean(k600_macIntyre[,2]), '\n')
cat('Read:', mean(k600_read[,2]), '\n')
cat('Read_soloviev:', mean(k600_soloviev[,2]), '\n')
cat('heiskanen:', mean(k600_heiskanen[,2]), '\n')
cat('vachon:', mean(k600_vachon[,2]), '\n')


