basedir = '/mnt/ExtraDrive1/Work/desktop_data/2022_papers/HSsurv/'
load(paste0(basedir,'data/MHO_3to7.rda'))
load(paste0(basedir,'data/dTerr.rda'))
source(paste0(basedir,'R/MCMCabu.R'))
source(paste0(basedir,'R/plot_ana.R'))

# ------------------------------------------------------------------------------
#                      Set these for each Stock
# ------------------------------------------------------------------------------

#read in the MCMC results from fitting the haul-out model
MHO = MHO_3to7
# give a base name for the results of this analysis
file_name = 'Anafigs_6'
# select which stock from the data.frame of all terrestrial counts
dstk = dTerr[dTerr$stockid == 6,]
# what is the latest year in the database?
maxyr = max(dTerr$yr)
# if this is the maximum year for modeling, set the number of years
# beginning in 1996
nyr = maxyr - 1995

# ------------------------------------------------------------------------------
#                      Preliminaries
# ------------------------------------------------------------------------------

# from Haulout results, extract the MCMC results for fixed effects
HObeta = MHO[['beta']]
# get the MCMC results for just the intercept (average haulout for standardized
# covariates)
beta0 = unlist(lapply(MHO[['beta']], function(x){x[[1]]}))
# from Haulout results, extract the MCMC results for random effects for animal
gam = MHO[['gam']]

# select only years from 1996 onwards
dstk = dstk[dstk$yr > 1995,]
# make sure that the polyid is a factor
dstk$polyid = as.factor(as.character(dstk$polyid))
# make sure that stockid is a factor
dstk$stockid = as.factor(as.character(dstk$stockid))
# get a list of all of the polyids
IDs = levels(dstk$polyid)
#order by datadatetime within speno, and remove duplicate times
dstk = dstk[order(dstk$polyid, dstk$survey_dt),]

# get polyids that have greater than 0 counts, and those that don't
pidnames = levels(dstk$polyid)
polyid_gt0 = NULL
polyid_0 = NULL
for(i in 1:max(as.integer(dstk$polyid))){
 if(max(dstk[dstk$polyid == pidnames[i],'count']) > 0) {
 polyid_gt0 = c(polyid_gt0, pidnames[i])} else {
	polyid_0 = c(polyid_0, pidnames[i])}
}

#create a working data.frame with sites with greater than 0 counts
dstk_gt0 = dstk[dstk$polyid %in% polyid_gt0,]

# in case we removed the only record of a stock, reset the factors so that
# we don't have factors without any records
dstk_gt0$polyid = as.factor(as.character(dstk_gt0$polyid))

# create a matrix with integers for polyid, and year
Ayr = cbind(as.integer(dstk_gt0$polyid),dstk_gt0$yr)
# subtract 1995, so this matrix serves as an index for polyid and year
A = cbind(as.integer(dstk_gt0$polyid), 
	dstk_gt0$yr - 1995)
# label them
colnames(A) = c('polyint', 'yr')

# design matrix
Xt <- model.matrix(count ~ dystd + I(dystd^2) + tdstd + I(tdstd^2)
  + hrstd + I(hrstd^2), data = dstk_gt0)

# create a data.frame for Poisson modeling of site intercepts
Df = data.frame(polyid = as.factor(as.integer(dstk_gt0$polyid)),
	yr = dstk_gt0$yr - 1995, 
	count = dstk_gt0$count)
# get site intercepts from Poisson regression
# and adjust them upwards by dividing by average haulout probability 
# (for standardized covariates)
# these will be fixed during MCMC because the ICAR model does not have a mean
# so it adjusts automatically to whatever means we choose
tau = coef(glm(count ~ polyid - 1, data = Df, family = 'poisson')) - 
 (mean(beta0) + mean(unlist(gam))) + 
 log(1 + exp(mean(beta0) + mean(unlist(gam))))
		
# initialize MCMC parameters
# random effects for day
REday = rep(0, times = length(dstk_gt0$count))
# temporal random walk random effect model within polyid
# that is, a random walk within each row, along the columns
del = matrix(0, nrow = length(tau), ncol = nyr)
# start the loglam matrix with polyid intercepts (tau)
loglam = matrix(tau, nrow = length(tau), ncol = nyr)
# start the abundance matrix with exponential of loglam
Nmat = ceiling(exp(loglam))
# initial variance parameter of the random walk
sigmaAR1 = 0.2

# initialize tuning parameters
# initial variance of the proposal distributions for the random effects for day
REday_tune = rep(.2, times = dim(dstk_gt0)[1])
# initial variance parameters for the proposal distributions for random walks 
# for each site
sigAR1_tune = 0.005

#-------------------------------------------------------------------------------
#                    Initital MCMC
#-------------------------------------------------------------------------------

set.seed(4001)

start = Sys.time()
M = MCMCabu(
  niter = 100, HObeta = HObeta, gam = gam, co = dstk_gt0$count, 
  Xt = Xt, A = A, Nmat = Nmat, del = del, tau = tau, 
  REday = REday, 
  sigmaAR1 = sigmaAR1, REday_tune = REday_tune, 
  sigAR1_tune = sigAR1_tune)
end = Sys.time()
end - start

#-------------------------------------------------------------------------------
#                    Tune Acceptance Rates and Burn-in
#-------------------------------------------------------------------------------

# continue sampling by picking up with the last MCMC sample previosly,
# but now we are going to tune the acceptance rates
burnin_start = Sys.time()
undebug(MCMCabu)
for(burnin in 1:50) {
  nM = length(M[['sigmaAR1']])
  M = MCMCabu(
		niter = 100, thin = 1, 
    HObeta = HObeta, gam = gam, co = dstk_gt0$count, 
    Xt = Xt, A = A,
    Nmat = M[['Nmat']][[nM]],
    del = M[['del']][[nM]],
    tau = tau, 
    REday = M[['REday']][[nM]], 
    sigmaAR1 = M[['sigmaAR1']][[nM]], 
		REday_tune = REday_tune,
    sigAR1_tune = sigAR1_tune
  )

  REday_tune[M[['REday_acc_rate']] < .21] = 
		REday_tune[M[['REday_acc_rate']] < .21] - 0.1
  REday_tune[M[['REday_acc_rate']] > .55] = 
		REday_tune[M[['REday_acc_rate']] > .55] + 0.5
  REday_tune[REday_tune <= 0] = .1
  REday_tune[REday_tune >= 2] = 2.1

  if(M[['sAR1_acc_rate']] < .21) sigAR1_tune = sigAR1_tune - 0.001
  if(M[['sAR1_acc_rate']] > .55) sigAR1_tune = sigAR1_tune + 0.005
  if(sigAR1_tune <= 0) sigAR1_tune = .001
  
}
burnin_end = Sys.time()
burnin_time = burnin_end - burnin_start

# now let autocorrelated errors mix after tuning, still with burn-in priors
nMCMC = 20000
thin = 20
MCMC_start = Sys.time()
  nM = length(M[['sigmaAR1']])
  #undebug(MCMCabu)
  M = MCMCabu(
		niter = nMCMC, thin = thin, 
    HObeta = HObeta, gam = gam, co = dstk_gt0$count, 
    Xt = Xt, A = A,
    Nmat = M[['Nmat']][[nM]],
    del = M[['del']][[nM]],
    tau = tau, 
    REday = M[['REday']][[nM]], 
    sigmaAR1 = M[['sigmaAR1']][[nM]], 
		REday_tune = REday_tune,
    sigmaYr_tune = sigmaYr_tune, sigAR1_tune = sigAR1_tune
  )
MCMC_end = Sys.time()
MCMC_time = MCMC_end - MCMC_start


# check the results so far after burn-in period
path = paste0(basedir, 'inst/scripts/ana_counts/figures/', 
	file_name,'_burnin.pdf')
plot_ana(M, dstk_gt0, path)

#-------------------------------------------------------------------------------
#          Tune Acceptance Rates One More Time and the Final Sample
#-------------------------------------------------------------------------------

# continue sampling by picking up with the last MCMC sample previosly,
# but now we are going to tune the acceptance rates without burn-in priors
burnin_start = Sys.time()
for(burnin in 1:50) {
source(paste0(basedir,'R/MCMCabu1.R'))
  nM = length(M[['sigmaAR1']])
  #undebug(MCMCabu)
  M = MCMCabu(
		niter = 100, thin = 1, 
    HObeta = HObeta, gam = gam, co = dstk_gt0$count, 
    Xt = Xt, A = A, 
    Nmat = M[['Nmat']][[nM]],
    del = M[['del']][[nM]],
    tau = tau, 
    REday = M[['REday']][[nM]],
		sigmaAR1 = M[['sigmaAR1']][[nM]], 
		REday_tune = REday_tune,	
		sigAR1_tune = sigAR1_tune
  )

  REday_tune[M[['REday_acc_rate']] < .21] = 
		REday_tune[M[['REday_acc_rate']] < .21] - 0.1
  REday_tune[M[['REday_acc_rate']] > .55] = 
		REday_tune[M[['REday_acc_rate']] > .55] + 0.5
  REday_tune[REday_tune <= 0] = .1
  REday_tune[REday_tune >= 2] = 2

  if(M[['sAR1_acc_rate']] < .21) sigAR1_tune = sigAR1_tune - 0.001
  if(M[['sAR1_acc_rate']] > .55) sigAR1_tune = sigAR1_tune + 0.005
  if(sigAR1_tune <= 0) sigAR1_tune = .001

}
burnin_end = Sys.time()
burnin_time = burnin_end - burnin_start


# final MCMC sampling
nMCMC = 100000
thin = 100
MCMC_start = Sys.time()
  nM = length(M[['sigmaAR1']])
  #undebug(MCMCabu)
  M = MCMCabu(
		niter = nMCMC, thin = thin, 
    HObeta = HObeta, gam = gam, co = dstk_gt0$count, 
    Xt = Xt, A = A, 
    Nmat = M[['Nmat']][[nM]],
    del = M[['del']][[nM]],
    tau = tau, 
    REday = M[['REday']][[nM]], 
    sigmaAR1 = M[['sigmaAR1']][[nM]], 
		REday_tune = REday_tune,	
		sigAR1_tune = sigAR1_tune
  )
MCMC_end = Sys.time()
MCMC_time = MCMC_end - MCMC_start

# add row names (polyid) and column names (year) to Nmat, the 
# estimated abundance matrices to keep things labelled for later
# analysis.  Also give stockid as an attribute
for(i in 1:1000) {
	Nmat = M[['Nmat']][[i]]
	rownames(Nmat) = levels(dstk_gt0$polyid)
	colnames(Nmat) = 1996:(1995 + dim(Nmat)[2])
	attr(Nmat,'stockid') = as.integer(as.character(
		dstk_gt0[!duplicated(dstk_gt0$polyid),c('stockid')]))
	M[['Nmat']][[i]] = Nmat
}

#-------------------------------------------------------------------------------
#          Set and Store for each Stock
#-------------------------------------------------------------------------------

Mana_6 = M
save(Mana_6, file = paste0(basedir,'data/Mana_6.rda'))

#-------------------------------------------------------------------------------
#          Diagnostic Graphics
#-------------------------------------------------------------------------------

path = paste0(basedir, 'inst/scripts/ana_counts/figures/', 
	file_name,'.pdf')
plot_ana(M, dstk_gt0, path)


