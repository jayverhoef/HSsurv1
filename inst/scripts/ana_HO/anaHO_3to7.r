# set directories and load haulout and terrestria count data
basedir = '/mnt/ExtraDrive1/Work/desktop_data/2022_papers/HSsurv/'
load(paste0(basedir,'data/dHOterr.rda'))
table(dHOterr$stockid)
load(paste0(basedir,'data/dTerr.rda'))
table(dTerr$stockid)
source(paste0(basedir,'R/MCMCHO.R'))
source(paste0(basedir,'R/plot_HO.R'))

# ------------------------------------------------------------------------------
#                      Set these for each run
# ------------------------------------------------------------------------------

file_name = 'HOfigs_3to7'

# get the aerial survey data for only stocks 3 through 7
dstk = dTerr[dTerr$stockid == 3 | dTerr$stockid == 4 | dTerr$stockid == 5 | 
	dTerr$stockid == 6 | dTerr$stockid == 7,]
# get the haul-out data for only stocks 3 through 7
dHO = dHOterr[dHOterr$stockid == 3 | dHOterr$stockid == 4 | 
	dHOterr$stockid == 5 | dHOterr$stockid == 6 | dHOterr$stockid == 7,]

# ------------------------------------------------------------------------------
#                      Preliminaries
# ------------------------------------------------------------------------------

# only use data after 1995
dstk = dstk[dstk$yr > 1995,]
# make sure polyid and stockid are factors in count data
dstk$polyid = as.factor(as.character(dstk$polyid))
dstk$stockid = as.factor(as.character(dstk$stockid))
IDs = levels(dstk$polyid)
#order by datadatetime within speno, and remove duplicate times
dstk = dstk[order(dstk$polyid, dstk$survey_dt),]

# make sure speno is factor in haulout data
dHO$speno = as.factor(as.character(dHO$speno))
HOIDs = levels(dHO$speno)
#order by datadatetime within speno
dHO = dHO[order(dHO$speno, dHO$date_time),]
# check for any duplicated times (must be increasing within animal)
any(dHO$yrhr0[2:dim(dHO)[1]] - dHO$yrhr0[1:(dim(dHO)[1] - 1)] == 0)

# inspect average haulout proportions by seal ID
y = dHO$y
# create binary data from binned data
y[dHO$y < 0.5] = 0
y[y != 0] = 1
# get the mean binned data to compare to mean binary data, by seal ID
DrybyID = NULL
HObyID = NULL
for(i in 1:length(levels(dHO$speno))){
  DrybyID = cbind(DrybyID, mean(dHO$y[dHO$speno == levels(dHO$speno)[i]]))
  HObyID = cbind(HObyID, mean(y[dHO$speno == levels(dHO$speno)[i]]))
}
# create a data.frame for comparing the two
HOdry = data.frame(ID = levels(dHO$speno), 
	DrybyID = as.vector(DrybyID), HObyID = as.vector(HObyID))
HOdry
# plot the average dry value for binned data versus the binary data
plot(c(0,1),c(0,1), type = 'l', lwd = 2, xlab = 'Recorded Dry Means by ID', 
	ylab = 'Binary Data Means by ID')
points(HOdry[,2:3], pch = 19)

# create design matrix
X <- model.matrix(y ~ dystd + I(dystd^2) + tdstd + I(tdstd^2)
  + hrstd + I(hrstd^2), data = dHO)
cbind(X[,6],dHO$hrstd)
# use Bernoulli distribution to get some initial parameters assuming independence
# create likelihood function as function of parameters and data
LLho = function(theta,y,X)
{
  # theta[1] is phi
	# rest of theta is part of linear model
	Xb = X %*% theta
	p = exp(Xb)/(1 + exp(Xb))
  -sum(dbern(y, p, log = TRUE))
}
# optimize the likelihood
optout1 = optim(c(.5,.5,-.5,.5,-.5,.5,-.5), LLho, y = y, X = X)

# starting values of fixed effects coefficients 
beta = optout1$par
# starting value for autocorrelation parameter
alpha = .8
# starting value for variance of random effects for animal
sigGam = 1
# starting value for variance of temporally-autocorrelated random effect
sigEps = 1
# starting values for vector of integers that identify an animal
Zi = as.integer(dHO$speno)
# initialize vector of random effects for animal
gam = rep(0, times = max(Zi))
# empty list of time differences
tdif = list()
# empty list of temporally-autocorrelated random effects
eps = list()
# fill in time differences by animal
for(i in 1:length(levels(dHO$speno))) {
	hrs = dHO[dHO$speno == levels(dHO$speno)[i],'yrhr0']
	tdif[[i]] = hrs[2:length(hrs)] -  hrs[1:(length(hrs) - 1)]
}
# initialize temporally-autocorrelated random effects by animal
for(i in 1:length(tdif)){
  eps[[i]] = rep(0, times = sum(dHO$speno == levels(dHO$speno)[i]))
}

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

# MCMC tuning parameter for Metropolis-Hasting
# these control the variance of the proposal distribution for each parameter
beta_tune = rep(0.02, times = length(beta))
alpha_tune = .02
sigGam_tune = 0.02
sigEps_tune = 0.02
gam_tune = rep(.05, times = length(gam))
eps_tune = rep(.02, times = length(gam))

#-------------------------------------------------------------------------------
#                    Initital MCMC
#-------------------------------------------------------------------------------

# set the random number seed so results are reproducible
set.seed(4001)
# an initial run to populate the outputs, all contained in M
start = Sys.time()
#undebug(MCMCHO)
M = MCMCHO(niter = 100, beta = beta, alpha = alpha, 
	sigGam = sigGam, sigEps = sigEps, gam = gam, eps = eps, 
	y = y, X = X, Zi = Zi, tdif = tdif,
	beta_tune = beta_tune, alpha_tune = alpha_tune, 
	sigGam_tune = sigGam_tune, sigEps_tune = sigEps_tune,
	gam_tune = gam_tune, eps_tune = eps_tune)
end = Sys.time()
end - start

#-------------------------------------------------------------------------------
#                    Tune Acceptance Rates and Burn-in
#-------------------------------------------------------------------------------

# do a series of burn-ins (50), as a loop, where we fine tune the tuning 
# parameters by monitoring the acceptance rates, and then adjusting accordingly
# We constantly update M and pass on current parameters to next iteration
start = Sys.time()
for(burnin in 1:50) {
	nM = length(M[['alpha']])
	#undebug(MCMCHO)
	M = MCMCHO(niter = 100, 
		beta = unlist(M[['beta']][nM]), 
		alpha = as.numeric(M[['alpha']][nM]),
		sigGam = as.numeric(M[['sigGam']][nM]), 
		sigEps = as.numeric(M[['sigEps']][nM]), 
		gam = unlist(M[['gam']][nM]), 
		eps = M[['eps']],
		y = y, X = X, Zi = Zi, tdif = tdif,
		beta_tune = beta_tune, alpha_tune = alpha_tune, 
		sigGam_tune = sigGam_tune, sigEps_tune = sigEps_tune,
		gam_tune = gam_tune, eps_tune = eps_tune)

    # in all the follows, adjust the tuning parameter (proposal variance) up
    # or down depending on whether the acceptance rate is less than 0.21 or
    # greater than 0.55
		beta_tune[M[['beta_acc_rate']] < .21] = 
			beta_tune[M[['beta_acc_rate']] < .21] - 0.001
		beta_tune[M[['beta_acc_rate']] > .55] = 
			beta_tune[M[['beta_acc_rate']] > .55] + 0.005
	  beta_tune[beta_tune <= 0] = .001
	  
		if(M[['alpha_acc_rate']] < .21) alpha_tune = alpha_tune - .01
		if(M[['alpha_acc_rate']] > .55) alpha_tune = alpha_tune + .01
		if(alpha_tune <= 0) alpha_tune = .01

		if(M[['sigGam_acc_rate']] < .21) sigGam_tune = sigGam_tune - .01
		if(M[['sigGam_acc_rate']] > .55) sigGam_tune = sigGam_tune + .05
		if(sigGam_tune <= 0) sigGam_tune = .01

		if(M[['sigEps_acc_rate']] < .21) sigEps_tune = sigEps_tune - .01
		if(M[['sigEps_acc_rate']] > .55) sigEps_tune = sigEps_tune + .05
		if(sigEps_tune <= 0) sigEps_tune = .01
		
		gam_tune[M[['gam_acc_rate']] < .21] = 
			gam_tune[M[['gam_acc_rate']] < .21] - 0.01
		gam_tune[M[['gam_acc_rate']] > .55] = 
			gam_tune[M[['gam_acc_rate']] > .55] + 0.02
		gam_tune[gam_tune <= 0] = .01
		
		eps_tune[M[['eps_acc_rate']] < .21] = 
			eps_tune[M[['eps_acc_rate']] < .21] - 0.01
		eps_tune[M[['eps_acc_rate']] > .55] = 
			eps_tune[M[['eps_acc_rate']] > .55] + 0.02
		eps_tune[eps_tune <= 0] = .01
								
}
end = Sys.time()
end - start

# Assuming that we have good proposal distributions now,
# do a burn-in run of 20,000 samples where we thin by 20,
# so we end up storing only 1,000 samples
start = Sys.time()
	nM = length(M[['alpha']])
	#undebug(MCMCHO)
	M = MCMCHO(niter = 20000, thin = 20,
		beta = unlist(M[['beta']][nM]), 
		alpha = as.numeric(M[['alpha']][nM]),
		sigGam = as.numeric(M[['sigGam']][nM]), 
		sigEps = as.numeric(M[['sigEps']][nM]), 
		gam = unlist(M[['gam']][nM]), 
		eps = M[['eps']],
		y = y, X = X, Zi = Zi, tdif = tdif,
		beta_tune = beta_tune, alpha_tune = alpha_tune, 
		sigGam_tune = sigGam_tune, sigEps_tune = sigEps_tune,
		gam_tune = gam_tune, eps_tune = eps_tune)
end = Sys.time()
end - start

#-------------------------------------------------------------------------------
#          Tune Acceptance Rates One More Time and the Final Sample
#-------------------------------------------------------------------------------

# do a series of burn-ins (50), as a loop, where we fine tune the tuning 
# parameters by monitoring the acceptance rates, and then adjusting accordingly
# We constantly update M and pass on current parameters to next iteration
start = Sys.time()
for(burnin in 1:50) {
	nM = length(M[['alpha']])
	#undebug(MCMCHO)
	M = MCMCHO(niter = 100, 
		beta = unlist(M[['beta']][nM]), 
		alpha = as.numeric(M[['alpha']][nM]),
		sigGam = as.numeric(M[['sigGam']][nM]), 
		sigEps = as.numeric(M[['sigEps']][nM]), 
		gam = unlist(M[['gam']][nM]), 
		eps = M[['eps']],
		y = y, X = X, Zi = Zi, tdif = tdif,
		beta_tune = beta_tune, alpha_tune = alpha_tune, 
		sigGam_tune = sigGam_tune, sigEps_tune = sigEps_tune,
		gam_tune = gam_tune, eps_tune = eps_tune)

    # in all the follows, adjust the tuning parameter (proposal variance) up
    # or down depending on whether the acceptance rate is less than 0.21 or
    # greater than 0.55
		beta_tune[M[['beta_acc_rate']] < .21] = 
			beta_tune[M[['beta_acc_rate']] < .21] - 0.001
		beta_tune[M[['beta_acc_rate']] > .55] = 
			beta_tune[M[['beta_acc_rate']] > .55] + 0.005
	  beta_tune[beta_tune <= 0] = .001
	  
		if(M[['alpha_acc_rate']] < .21) alpha_tune = alpha_tune - .01
		if(M[['alpha_acc_rate']] > .55) alpha_tune = alpha_tune + .01
		if(alpha_tune <= 0) alpha_tune = .01

		if(M[['sigGam_acc_rate']] < .21) sigGam_tune = sigGam_tune - .01
		if(M[['sigGam_acc_rate']] > .55) sigGam_tune = sigGam_tune + .05
		if(sigGam_tune <= 0) sigGam_tune = .01

		if(M[['sigEps_acc_rate']] < .21) sigEps_tune = sigEps_tune - .01
		if(M[['sigEps_acc_rate']] > .55) sigEps_tune = sigEps_tune + .05
		if(sigEps_tune <= 0) sigEps_tune = .01
		
		gam_tune[M[['gam_acc_rate']] < .21] = 
			gam_tune[M[['gam_acc_rate']] < .21] - 0.01
		gam_tune[M[['gam_acc_rate']] > .55] = 
			gam_tune[M[['gam_acc_rate']] > .55] + 0.02
		gam_tune[gam_tune <= 0] = .01
		
		eps_tune[M[['eps_acc_rate']] < .21] = 
			eps_tune[M[['eps_acc_rate']] < .21] - 0.01
		eps_tune[M[['eps_acc_rate']] > .55] = 
			eps_tune[M[['eps_acc_rate']] > .55] + 0.02
		eps_tune[eps_tune <= 0] = .01
								
}
end = Sys.time()
end - start


# Assuming that we have good proposal distributions now,
# do a final run of 100,000 samples where we thin by 100,
# so we end up storing only 1,000 samples
start = Sys.time()
	nM = length(M[['alpha']])
	#undebug(MCMCHO)
	M = MCMCHO(niter = 100000, thin = 100,
		beta = unlist(M[['beta']][nM]), 
		alpha = as.numeric(M[['alpha']][nM]),
		sigGam = as.numeric(M[['sigGam']][nM]), 
		sigEps = as.numeric(M[['sigEps']][nM]), 
		gam = unlist(M[['gam']][nM]), 
		eps = M[['eps']],
		y = y, X = X, Zi = Zi, tdif = tdif,
		beta_tune = beta_tune, alpha_tune = alpha_tune, 
		sigGam_tune = sigGam_tune, sigEps_tune = sigEps_tune,
		gam_tune = gam_tune, eps_tune = eps_tune)
end = Sys.time()
end - start

#-------------------------------------------------------------------------------
#          Set and Store 
#-------------------------------------------------------------------------------

# save the results to disk
 MHO_3to7 = M
 save(MHO_3to7, file = paste0(basedir,'data/MHO_3to7.rda'))

#-------------------------------------------------------------------------------
#          Diagnostic Graphics
#-------------------------------------------------------------------------------

path = paste0(basedir, 'inst/scripts/ana_HO/figures/', 
	file_name,'.pdf')

plot_HO(M, dHO, dstk, X, y, path) 

