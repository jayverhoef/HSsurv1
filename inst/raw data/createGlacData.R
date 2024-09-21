library(sf)
library(solaR)

#set path shortcut
basedir = '/mnt/ExtraDrive1/Work/desktop_data/2022_papers/HSsurv/'

#read in the data
d1 = read.csv(paste0(basedir,'/inst/raw data/', 
	'GlacialHarborSealCounts_Thru2018_20240802_smk.csv'))
# change some column names to match terrestrial data name
names(d1)[names(d1) == 'num_seals' ] = 'count'
 
# check for numbers of missing data
sum(is.na(d1$survey_dt))
sum(is.na(d1$count))
# which have missing counts?
d1[is.na(d1$count),
	c('polyid', 'survey_dt_gmt','glacier_name','counting_method','count')]

# remove missing dates
d2 = d1[!is.na(d1$survey_dt),]
# check for missing counts
d2 = d2[!is.na(d2$count),]

# convert character dates to POSIXct
d2$survey_dt = as.POSIXct(d2$survey_dt_gmt,
	format="%Y-%m-%d %H:%M:%S",tz='UTC')

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
# see if they are close (of course, solar noon will not be exactly noon)
plot(d2$hr, d2$solhr)
# which have d2$hr > 20
d2[d2$hr > 20,c('polyid', 'survey_dt','glacier_name','counting_method','count')]

# order by stock, polyid, year, and day
d2 = d2[order(d2$stockid, d2$polyid, d2$yr, d2$yday), ]

# standardize explanatory variables
d2$dystd = (d2$yday - as.POSIXlt("2004-08-15")$yday)/30
d2$hrstd = (d2$solhr - 12)/6

# make sure polyid is a factor
d2$polyid = as.factor(as.character(d2$polyid))

# save dataset
dGlac = d2
save(dGlac, file = paste0(basedir, '/data/dGlac.rda'))

