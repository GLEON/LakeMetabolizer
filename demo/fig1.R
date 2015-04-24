
##Figure showing example diel DO signal

library(LakeMetabolizer)
library(lubridate)
library(plyr)
library(dplyr)

#load data
data.path = system.file('figdata', package="LakeMetabolizer")

me = load.ts(file.path(data.path, 'troutlake.doobs'))
sp = load.ts(file.path(data.path, 'sparkling.doobs'))

names(me)[2] = 'doobs'
names(sp)[2] = 'doobs'

me$dayfrac = hour(me$datetime)/24
sp$dayfrac = hour(sp$datetime)/24

me$day     = as.POSIXct(trunc(me$datetime, units='days'))
sp$day     = as.POSIXct(trunc(sp$datetime, units='days'))

me = ddply(me, 'day', function(df){
				df$do_anom = df$doobs - mean(df$doobs)
				return(df)
					})
sp = ddply(sp, 'day', function(df){
	df$do_anom = df$doobs - mean(df$doobs)
	return(df)
})

me_day = group_by(me, dayfrac) %>% summarise(avg_do=mean(do_anom, na.rm=TRUE))
sp_day = group_by(sp, dayfrac) %>% summarise(avg_do=mean(do_anom, na.rm=TRUE))

sp_day$dayfrac = as.POSIXct(sp_day$dayfrac*3600*24, origin='2011-01-01', tz="UTC")
me_day$dayfrac = as.POSIXct(me_day$dayfrac*3600*24, origin='2011-01-01', tz="UTC")

plot(me_day, type='l', col='green', ylab="DO Anomaly", xlab="Time")
lines(sp_day, col='blue', ylab="")

