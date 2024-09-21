library(solaR)
library(sf)
#set path shortcut
basedir = '/mnt/ExtraDrive1/Work/desktop_data/2022_papers/HSsurv/'

#read in the data
load(paste0(basedir,'inst/raw data/tbl_percent_glacial.rda'))
# copy data into working data.frame
dt = tbl_percent_glacial
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

