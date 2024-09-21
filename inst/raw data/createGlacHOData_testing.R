library(solaR)
library(sf)
#set path shortcut
basedir = '/mnt/ExtraDrive1/Work/desktop_data/2022_papers/HSsurv/'

#read in the data
load(paste0(basedir,'inst/raw data/tbl_percent_glacial.rda'))
# copy data into working data.frame
dt = tbl_percent_glacial
# remove records with missing data_time
dt = dt[!is.na(dt$date_time),]
# remove records with missing percent_dry
dt = dt[!is.na(dt$percent_dry),]
load(paste0(basedir,'data/polylonlat.rda'))
# attach polyid centroids
dt = merge(dt, polylonlat, by = 'polyid')
# create solar hour
dt$solhr = as.POSIXlt(local2Solar(st_drop_geometry(dt)[,"date_time"], lon = dt$lon))$hour
# as sanity check, do time zone correction from GMT
dt$hr = (as.POSIXlt(st_drop_geometry(dt)[,"date_time"])$hour + 15)%%24
# clip data to dates between 15 July and 30 September (year doesn't matter)
dt$hour = as.POSIXlt(st_drop_geometry(dt)[,"date_time"])$hour
# clip data to dates between 15 July and 30 September (year doesn't matter)
dt_molt = dt[as.POSIXlt(dt$date_time)$yday > as.POSIXlt("2004-07-15")$yday & 
  as.POSIXlt(dt$date_time)$yday < as.POSIXlt("2004-09-30")$yday,]
dt_pup = dt[as.POSIXlt(dt$date_time)$yday > as.POSIXlt("2004-05-01")$yday & 
  as.POSIXlt(dt$date_time)$yday < as.POSIXlt("2004-07-16")$yday,]
# create a data.frame for plotting
plot(dt_molt$hr, dt_molt$solhr)
plot(dt_molt$hr, dt_molt$hour)


load('/mnt/ExtraDrive1/Work/desktop_data/2022_papers/HSsurv/data/akpvpolys_sf.rda')
spenos = sort(unique(dt_molt$speno))
pdf(file = paste0(basedir,'inst/raw data/testing_solhr_molt.pdf'))
	plot(st_geometry(dt_molt), pch = 19, col = 'red')
	plot(st_geometry(akpvpolys_sf), add = TRUE)
	for( i in 1:25) {
		layout(matrix(1:2, nrow = 2))
		tmp = dt_molt[dt_molt$speno == spenos[i],]
		plot1data = aggregate(tmp$percent_dry, by = list(tmp$solhr), FUN = mean)
		plot(plot1data, type = 'l', lwd = 2, xlab = 'Solar Hour', ylim = c(0, 110),
			main = paste(spenos[i], '; sample size = ', dim(tmp)[1]))
		points(plot1data, pch = 19, cex = 2)

		plot(st_geometry(akpvpolys_sf), xlim = c(923352, 1339737),
			ylim = c(950079, 1209891))
		plot(st_geometry(tmp), pch = 19, col = 'red', add = TRUE)
	}
dev.off()

pdf(file = paste0(basedir,'inst/raw data/testing_AKhr.pdf'))
	plot(st_geometry(dt_molt), pch = 19, col = 'red')
	plot(st_geometry(akpvpolys_sf), add = TRUE)
	for( i in 1:25) {
		layout(matrix(1:2, nrow = 2))
		tmp = dt_molt[dt_molt$speno == spenos[i],]
		plot1data = aggregate(tmp$percent_dry, by = list(tmp$hr), FUN = mean)
		plot(plot1data, type = 'l', lwd = 2, xlab = 'Alaska Hour', ylim = c(0, 110),
			main = paste(spenos[i], '; sample size = ', dim(tmp)[1]))
		points(plot1data, pch = 19, cex = 2)

		plot(st_geometry(akpvpolys_sf), xlim = c(923352, 1339737),
			ylim = c(950079, 1209891))
		plot(st_geometry(tmp), pch = 19, col = 'red', add = TRUE)
	}
dev.off()

pdf(file = paste0(basedir,'inst/raw data/testing_rawhr.pdf'))
	plot(st_geometry(dt_molt), pch = 19, col = 'red')
	plot(st_geometry(akpvpolys_sf), add = TRUE)
	for( i in 1:25) {
		layout(matrix(1:2, nrow = 2))
		tmp = dt_molt[dt_molt$speno == spenos[i],]
		plot1data = aggregate(tmp$percent_dry, by = list(tmp$hour), FUN = mean)
		plot(plot1data, type = 'l', lwd = 2, xlab = 'Alaska Hour', ylim = c(0, 110),
			main = paste(spenos[i], '; sample size = ', dim(tmp)[1]))
		points(plot1data, pch = 19, cex = 2)

		plot(st_geometry(akpvpolys_sf), xlim = c(923352, 1339737),
			ylim = c(950079, 1209891))
		plot(st_geometry(tmp), pch = 19, col = 'red', add = TRUE)
	}
dev.off()

pdf(file = paste0(basedir,'inst/raw data/testing_solhr_pup.pdf'))
	plot(st_geometry(dt_pup), pch = 19, col = 'red')
	plot(st_geometry(akpvpolys_sf), add = TRUE)
	for( i in 1:25) {
		layout(matrix(1:2, nrow = 2))
		if(any(dt_pup$speno == spenos[i])) {
			tmp = dt_pup[dt_pup$speno == spenos[i],]
			plot1data = aggregate(tmp$percent_dry, by = list(tmp$solhr), FUN = mean)
			plot(plot1data, type = 'l', lwd = 2, xlab = 'Solar Hour', ylim = c(0, 110),
				main = paste(spenos[i], '; sample size = ', dim(tmp)[1]))
			points(plot1data, pch = 19, cex = 2)

			plot(st_geometry(akpvpolys_sf), xlim = c(923352, 1339737),
				ylim = c(950079, 1209891))
			plot(st_geometry(tmp), pch = 19, col = 'red', add = TRUE)
		}
	}
dev.off()

# drop the sf geometry
dt = st_drop_geometry(dt)
dt = as.data.frame(dt[,1:12])
# load sf object from harbor seal polygon data
load(paste0(basedir,'data/akpvpolys_sf.rda'))
# drop the sf geometry
DF  = st_drop_geometry(akpvpolys_sf)
DFstock = DF[,2:3]
# attach stockid (factor with 12 levels labeled 1,...,12) by using polyid
dt = merge(dt, DFstock, by = 'polyid')
#data(polylonlat)
# load centroids of polyids
load(paste0(basedir,'data/polylonlat.rda'))
# attach polyid centroids
dt = merge(dt, polylonlat, by = 'polyid')
# make sure polyid is a factor
dt$polyid = as.factor(as.character(dt$polyid))
# remove records with missing data_time
dt = dt[!is.na(dt$date_time),]
# create columns for year, day-of-year, and solar hour-of-day
dt$yr = as.POSIXlt(dt[,"date_time"])$year + 1900
dt$dy = as.POSIXlt(dt[,"date_time"])$yday
dt$solhr = as.POSIXlt(local2Solar(dt[,"date_time"], lon = dt$lon))$hour
# as a sanity check, take UCT/GMT - 9 hours to see if solar time makes sense
dt$hr = (as.POSIXlt(dt[,"date_time"])$hour + 15)%%24
# clip data to dates between 15 July and 30 September
dt = dt[dt$dy > as.POSIXlt("2004-07-15")$yday & 
  dt$dy < as.POSIXlt("2004-09-30")$yday,]

# get list of unique seal IDs (speno)
spenoList = unique(dt$speno)
#order by datadatetime within speno, and remove duplicate times
  dupTF = NULL
  d2 = NULL
  for(i in 1:length(spenoList)) {
    tmp = dt[dt$speno == spenoList[i],]
    ntmp = length(tmp[,1])
    oindx = order(tmp$date_time)
    tmp = tmp[oindx,]
    dupTF = c(dupTF,
      any(duplicated(tmp$date_time)) )
    tmp = tmp[!duplicated(tmp$date_time),] 
    d2 = rbind(d2,tmp)
  }
dt = d2
dt$speno = as.factor(as.character(dt$speno))

# create hour-of-year variable (for temporal autocorrelation modeling)
# and start from 0 for each seal
dt$yrhr = (dt$yr - min(dt$yr))*365*24 + dt$dy*24 + 
	as.POSIXlt(dt[,"date_time"])$hour
dt$yrhr0 = dt$yrhr
for(i in 1:length(levels(dt$speno)))
	dt[dt$speno == levels(dt$speno)[i],'yrhr0'] = 
		dt[dt$speno == levels(dt$speno)[i],'yrhr'] -
    min(dt[dt$speno == levels(dt$speno)[i],'yrhr'])

# order by yrhr0 within seal
dt = dt[order(dt$speno,dt$yrhr0),]
# change percent to proportion
dt$y = dt$percent_dry/100
# standardize explanatory variables
dt$dystd = (dt$dy-as.POSIXlt("2004-08-15")$yday)/30
dt$hrstd = (dt$solhr - 12)/6

dHOglac = dt
save(dHOglac, file = paste0(basedir,'/data/dHOglac.rda'))

