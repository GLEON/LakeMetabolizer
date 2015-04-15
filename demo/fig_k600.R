
# Generates a figure showing the results of the different k600 models
library(LakeMetabolizer)
data.path = system.file('extdata', package="LakeMetabolizer")

tb.data = load.all.data('sparkling', data.path)

ts.data = tb.data$data #pull out just the timeseries data

#calculate U10 and add it back onto the original

u10 = wind.scale(ts.data)
ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset

#Calculate k600 using the k.cole wind-based model
k600_cole = k.cole(ts.data)

k600_crusius = k.crusius(ts.data)

kd        = tb.data$metadata$averagekd
wnd.z      = 10   #because we converted to u10
atm.press  = 1018
lat       = tb.data$metadata$latitude
lake.area = tb.data$metadata$lakearea

#for k.read and k.macIntyre, we need LW_net.
#Calculate from the observations we have available.
lwnet = calc.lw.net(ts.data, lat, atm.press)
ts.data = merge(ts.data, lwnet)
## Not run: 
k600_read = k.read(ts.data, wnd.z=wnd.z, Kd=kd, atm.press=atm.press,
									 lat=lat, lake.area=lake.area)

k600_soloviev = k.read.soloviev(ts.data, wnd.z=wnd.z, Kd=kd,
																atm.press=atm.press, lat=lat, lake.area=lake.area)

k600_macIntyre = k.macIntyre(ts.data, wnd.z=wnd.z, Kd=kd, atm.press=atm.press)


#Create plot and save in user home directory (on Windows, Documents folder, on Mac, Home folder)
png('~/k600_figure.png', res=300, width=3000, height=1500)
cols = rainbow(5)
plot(k600_macIntyre, type='l', xlab='', col=cols[1], ylim=c(0,10), ylab=expression(k600~(m^-1)))
lines(k600_cole, col=cols[2])
lines(k600_crusius, col=cols[3])
lines(k600_read, col=cols[4])
lines(k600_soloviev, col=cols[5])

legend('topleft', legend=c('k.macIntyre', 'k.cole', 'k.crusius', 'k.read', 'k.read.soloviev'), lty=1, col=cols, horiz=FALSE)
dev.off()
