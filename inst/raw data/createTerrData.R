library(sf)
library(solaR)

#set path shortcut
basedir = '/mnt/ExtraDrive1/Work/desktop_data/2022_papers/HSsurv/'

#read in the data
d1 = read.csv(paste0(basedir,'/inst/raw data/', 
	'CoastalHarborSealCounts_Thru2023_20240502_smk.csv'))
# create total counts by adding pups and non-pups
d1$count = as.numeric(d1$non_pup + d1$pup)

# check for numbers of missing data
sum(is.na(d1$survey_dt))
sum(is.na(d1$count))
# remove missing dates
d2 = d1[!is.na(d1$survey_dt),]
# check for missing counts
d2 = d2[!is.na(d2$count),]
# check for missing low tide date/time
d2 = d2[!is.na(d2$nearest_low_dt),]

# convert character dates to POSIXct
d2$survey_dt = as.POSIXct(d2$survey_dt,
	format="%Y-%m-%d %H:%M:%S",tz='UTC')
d2$nearest_low_dt = as.POSIXct(d2$nearest_low_dt,
	format="%Y-%m-%d %H:%M:%S",tz='UTC')

# load sf object from harbor seal polygon data
load(paste0(basedir,'/data/akpvpolys_sf.rda'))
# drop the sf geometry
DF  = st_drop_geometry(akpvpolys_sf)
DFstock = DF[,2:4]
# attach stockid (factor with 12 levels labeled 1,...,12) by using polyid
d2 = merge(d2, DFstock, by = 'polyid')

# load centroids of polyids
load(paste0(basedir,'/data/polylonlat.rda'))
# attach polyid centroids
d2 = merge(d2, polylonlat, by = 'polyid')

# get earliest and latest day-of-year for analysis but using character dates
# (for any year) and using POSIXct and POSIXlt
# i.e., we will only use data between 15 July and 30 September
startday = as.POSIXlt(as.POSIXct("2005-07-15"))$yday
stopday = as.POSIXlt(as.POSIXct("2005-09-30"))$yday
# now trim the data to records that fit within the day-of-year constraints
d2 = d2[as.POSIXlt(d2$survey_dt)$yday >= startday &
	as.POSIXlt(d2$survey_dt)$yday <= stopday,]

# create year, day-of-survey as separate columns in dataset
d2$yr = as.POSIXlt(d2$survey_dt)$year + 1900
d2$yday = as.POSIXlt(d2$survey_dt)$yday
d2$solhr = as.POSIXlt(local2Solar(d2[,"survey_dt"], lon = d2$lon))$hour
# as a sanity check, take UCT/GMT - 9 hours to see if solar time makes sense
d2$hr = (as.POSIXlt(d2[,"survey_dt"])$hour + 15)%%24
# see if they are close.  Of course, solar hour will not be exactly noon.
plot(d2$hr, d2$solhr)
# how many sites with d2$hr < 2
sum(d2$hr < 2)
# delete this site as an outlier, probably wrong time entered
d2 = d2[d2$hr > 2]
# create minutes-from-low tide variable
d2$minutes_from_low = as.numeric(difftime(d2$nearest_low_dt, 
	d2$survey_dt, units = 'mins'))
# clip data to plus and minus 5 hours from low tide
d2 = d2[d2$minutes_from_low > -300 & 
  d2$minutes_from_low < 300,]
# order by stock, polyid, year, and day
d2 = d2[order(d2$stockid, d2$polyid, d2$yr, d2$yday), ]

# standardize explanatory variables
d2$dystd = (d2$yday - as.POSIXlt("2004-08-15")$yday)/30
# because we clip at plus/minus 5 hours (300 minutes), dividing by 150 will
# keep tdstd as a fraction between -2 and 2
d2$tdstd = d2$minutes_from_low/150
d2$hrstd = (d2$solhr - 12)/6

# make sure polyid is a factor
d2$polyid = as.factor(as.character(d2$polyid))

# save dataset
dTerr = d2
save(dTerr, file = paste0(basedir,'/data/dTerr.rda'))
