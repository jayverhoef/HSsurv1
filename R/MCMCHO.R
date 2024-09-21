#-------------------------------------------------------------------------------
#
#          MCMCHO
#
#-------------------------------------------------------------------------------

#' MCMC sampling for haulout data using beta regression with AR1 errors
#'
#' MCMC sampling for haulout data using beta regression with AR1 errors
#'
#' @param niter the number of MCMC iterations
#' @param thin thinning of the MCMC chain
#' @param kappa upper truncation point
#' @param beta vector of regression parameters on logit scale
#' @param alpha AR1 autocorrelation parameter on logit scale
#' @param sigGam variance of random effects per animal
#' @param sigEps variance of autocorrelated random effects
#' @param gam random effects (intercepts) for animals
#' @param eps list with autocorrelated random effects per animal
#' @param y response variable as a proportion (0 < y < 1) (strict inequalities)
#' @param X  linear model design matrix
#' @param Zi vector of integers that indicate which random effect per observation
#' @param tdif list with sequential time differences within animal
#' @param phi_tune tuning parameter (variance of proposal) for phi
#' @param beta_tune tuning parameters (variance of proposal) for beta
#' @param alpha_tune tuning parameter (variance of proposal) for alpha
#' @param sigGam_tune tuning parameter (variance of proposal) for sigGam
#' @param sigEps_tune tuning parameter (variance of proposal) for sigEps
#' @param gam_tune tuning parameters (variance of proposal) for gam
#' @param eps_tune tuning parameter (variance of proposal) for eps
#'
#' @return a list of the MCMC chains values. Only the last iteration of the autocorrelated errors eps is kept.
#'
#' @author Jay Ver Hoef
#' @export

MCMCHO <- function(niter = 1000, thin = 1, nepsiter = 10,
	beta, alpha, sigGam, sigEps, gam, eps, y, X, Zi, tdif, 
  beta_tune, alpha_tune, sigGam_tune, sigEps_tune,
  gam_tune, eps_tune, sample_beta = TRUE, sample_alpha = TRUE, 
  sample_gam = TRUE, sample_sigGam = TRUE, sample_eps = TRUE,
  sample_sigEps = TRUE)
{

  cat("\n")
  
  # ----------------------------------------------
  #          ITERATIONS WITH THINNING
  # ----------------------------------------------

	M <- vector("list", 16)
  M[[1]] = list()
	M[[6]] = list()
	M[[7]] = list()
	beta.accept = rep(0, times = length(beta))
	alpha.accept = 0
	sigGam.accept = 0
	sigEps.accept = 0
	gam.accept = rep(0, times = length(gam))
	eps.accept = rep(0, times = length(gam))
  ikeep = 0
	for(iter in 1:niter) {
		
	# each chunk of the MCMC code works in the same way as we scroll through
	# sampling for all parameters.  All parameters are sampled by Metropolis/
	# Hasting, whether or not we might be able to do Gibbs sampling.  This keeps
	# the sampling structure exactly the same for all parameters.
	# The Metropolis/Hastings sampling structure is described here, and no 
	# further details are given for other parameters as they follow the same
	# structure
	
		if(sample_beta == TRUE) {
	    # ---------- sampling for beta -----------

			# loop through each of the beta parameters
      for(i in 1:length(beta)) {
				U <- log(runif(1))
				betai.try = beta
				# here, the proposal is plus/minus a random normal centered on current
				# value with standard deviation given by beta_tune.
				betai.try[i] <- rnorm(1, beta[i], beta_tune[i])
				Zge = gam[Zi] + unlist(eps)
				LLdif1 = LLHObeta(betai.try, y, X, Zge) - 
					LLHObeta(beta, y, X, Zge)
				if(LLdif1 > U) {
					beta <- betai.try
					beta.accept[i] = beta.accept[i] + 1
				}
			}
		}
#	  browser()
		if(sample_gam == TRUE) {
	    # ---------- sampling for gam -----------

      Xbe = X %*% beta + unlist(eps)
      # loop through each of the gam values
		  for(j in 1:length(gam)) {
			  U <- log(runif(1))
		    gami.try = gam
		    # the proposal is a truncated normal at +- 6, so still symmetric
		    # (so no need for Hastings)
        gami.try[j] <- min(max(-6, rnorm(1, gam[j], gam_tune[j])), 6)
				LLdif2 = LLHOgam(gami.try[j], y[Zi == j], sigGam, 
						Xbe[Zi == j]) -
					LLHOgam(gam[j], y[Zi == j], sigGam, Xbe[Zi == j])
        if(LLdif2 > U) {
					gam <- gami.try
					gam.accept[j] = gam.accept[j] + 1
				}				
		  }
		 }

		if(sample_sigGam == TRUE) {
	    # ---------- sampling for sigGam -----------

		  U <- log(runif(1))
      sigGam.try <-  min(max(sigGam + runif(1)*sigGam_tune - 
				sigGam_tune/2, 0.00001), 20)
      LLdif3 = LLHOsgam(sigGam.try, gam) - LLHOsgam(sigGam, gam)
      if(LLdif3 > U) {
				sigGam <- sigGam.try
				sigGam.accept = sigGam.accept + 1
			}
		}

		if(sample_eps == TRUE) {
	    # ------------------- sampling for eps ---------------
	    # ---because we are using batch sampling (changes more slowly) ---
	    # ---- do it nepsiter times per sampling for other parameters ----

			XZg = X %*% beta + gam[Zi]
      for(epsiter in 1:nepsiter) {
				for(j in 1:length(eps)) {								
					U <- log(runif(1))
					nj = length(eps[[j]])
					epsj = eps[[j]]
					epsj.try= pmin(pmax(-6, rnorm(nj, epsj, eps_tune[j])), 6)
					LLdif4 = LLHOeps(alpha, sigEps, epsj.try, y[Zi == j], 
							tdif[[j]], XZg[Zi == j]) -
						LLHOeps(alpha, sigEps, epsj, y[Zi == j], 
							tdif[[j]], XZg[Zi == j])
					if(LLdif4 > U) { 
						eps[[j]] <- epsj.try
						eps.accept[j] = eps.accept[j] + 1
					}
				}
			}
		}	
		
		if(sample_alpha == TRUE) {
	    # ---------- sampling for alpha -----------

		  U <- log(runif(1))
      alpha.try <-  min(max(alpha + runif(1)*alpha_tune - alpha_tune/2, 
				-0.99999),0.99999)
      LLdif5 = LLHOalpha(alpha.try, eps, sigEps, tdif) -
				LLHOalpha(alpha, eps, sigEps, tdif)
      if(LLdif5 > U) {
				alpha <- alpha.try
				alpha.accept = alpha.accept + 1
			}
		}

		if(sample_sigEps == TRUE) {
	    # ---------- sampling for sigEps -----------
		  U <- log(runif(1))
      sigEps.try <-  min(max(sigEps + runif(1)*sigEps_tune - 
				sigEps_tune/2, 0.00001),20)
      LLdif6 = LLHOsigEps(alpha, eps, sigEps.try, tdif, length(y)) -
				LLHOsigEps(alpha, eps, sigEps, tdif, length(y))
       if(LLdif6 > U) {
				sigEps <- sigEps.try
				sigEps.accept = sigEps.accept + 1
			}
		}
		
		# evaluate the whole loglikelihood to keep track of its values during
		# MCMC
		LLHOout =  LLHO(phi, kappa, beta, alpha, sigGam, sigEps,
			gam, eps, y = y, X = X, Zi = Zi, tdif = tdif)
			
    if(iter%%thin == 0) {
      ikeep = ikeep + 1
      cat("\r", "Iter Number: ", iter)
      M[[1]][[ikeep]] <- beta
	    M[[2]] <- c(M[[2]],alpha)
		  M[[3]] <- c(M[[3]],sigGam)
      M[[4]] <- c(M[[4]],sigEps)
		  M[[5]] <- c(M[[5]],LLHOout$logLike)
	    M[[6]][[ikeep]] = gam
	    M[[14]] <- c(M[[14]],LLHOout$LLbern)
	    M[[15]] <- c(M[[15]],LLHOout$LLgam)
	    M[[16]] <- c(M[[16]],LLHOout$LLeps)
		}

  }

	# keep the last value of eps
  M[[7]] = eps
  # compute acceptance rates
  M[[8]] = beta.accept/niter	
  M[[9]] = alpha.accept/niter	
  M[[10]] = sigGam.accept/niter	
  M[[11]] = sigEps.accept/niter	
  M[[12]] = gam.accept/niter
  M[[13]] = eps.accept/(nepsiter*niter)
  cat("\n")
  names(M) <- c("beta", "alpha", "sigGam", 
    "sigEps", "logLike", "gam", "eps", 
    "beta_acc_rate","alpha_acc_rate","sigGam_acc_rate",
    "sigEps_acc_rate", "gam_acc_rate", "eps_acc_rate",
    "LLHObern", "LLHOgam", "LLHOeps")

  M
}




#-------------------------------------------------------------------------------
#
#          LLHObeta
#
#-------------------------------------------------------------------------------

#' portion of loglikelihood containing beta for haul-out model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param phi variance parameter
#' @param beta regression parameters
#' @param y response variable as a proportion (0 < y < 1) (strict inequalities)
#' @param X  linear model design matrix
#' @param Zge random part of mean Z gam + unlist(eps)
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export

LLHObeta = function(beta, y, X, Zge)
{
	Xb = X %*% beta + Zge
	p = exp(Xb)/(1 + exp(Xb))
  sum(dbern(y, p, log = TRUE))
}

#-------------------------------------------------------------------------------
#
#          LLHOgam
#
#-------------------------------------------------------------------------------

#' portion of loglikelihood containing gam for haul-out model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param gami random effect (intercepts) for ith animal
#' @param y response variable, must be a 0 or 1
#' @param phi variance parameter of beta distribution
#' @param sigGam variance of random effects per animal
#' @param Xbe part of the mean unaffected by gam, X beta + unlist(eps)
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export

LLHOgam = function(gami, y, sigGam, Xbe)
{ 
  Xb = Xbe + gami 
	p = exp(Xb)/(1 + exp(Xb))
  sum(dbern(y, p, log = TRUE)) + 
	  dnorm(gami, 0, sigGam, log = TRUE)
}

#-------------------------------------------------------------------------------
#
#          LLHOsgam
#
#-------------------------------------------------------------------------------

#' portion of loglikelihood containing sigGam for haul-out model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param sigGam variance of random effects per animal
#' @param gam random effects for animals
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export

LLHOsgam = function(sigGam, gam)
{
		sum(dnorm(gam, 0, sigGam, log = TRUE)) - 
			length(gam)*log(1 - 2*pnorm(-6, 0, sigGam))
}

#-------------------------------------------------------------------------------
#
#          LLHOeps
#
#-------------------------------------------------------------------------------

#' portion of loglikelihood containing eps for haul-out model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param phi variance parameter of beta distribution
#' @param alpha the autocorrelation parameter
#' @param sigEps the standard error parameter
#' @param epsj vector with autocorrelated errors per animal
#' @param y response variable per animal
#' @param tdif the time differences per animal
#' @param XZg part of the mean unaffected by eps, X beta + gam[Zi]
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export

LLHOeps = function(alpha, sigEps, epsj, y, tdif, XZg)
{ 
	Xb = XZg + epsj
	nj = length(epsj)
	p = exp(Xb)/(1 + exp(Xb))
  sum(dbern(y, p, log = TRUE)) +  
		sum(dnorm(epsj[2:nj], alpha^tdif*epsj[1:(nj - 1)], sigEps, log = TRUE)) +
		dnorm(epsj[1], 0, sigEps/(1 - alpha^2), log = TRUE)
}

#-------------------------------------------------------------------------------
#
#          LLHOalpha
#
#-------------------------------------------------------------------------------

#' portion of loglikelihood containing alpha for haul-out model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param alpha AR1 autocorrelation parameter on logit scale
#' @param sigEps variance of autocorrelated random effects
#' @param eps list with autocorrelated random effects per animal
#' @param tdif list with sequential time differences within animal
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export

LLHOalpha = function(alpha, eps, sigEps, tdif)
{ 
	sumall = 0
	for(j in 1:length(eps)) {
		epsj = eps[[j]]
		tdifj = tdif[[j]]
		nj = length(epsj)
		sumall = sumall + sum(dnorm(epsj[2:nj], alpha^tdifj*epsj[1:(nj - 1)], 
			sigEps, log = TRUE)) + 
			dnorm(epsj[1], 0, sigEps/(1 - alpha^2), log = TRUE)
	}
	sumall
}

#-------------------------------------------------------------------------------
#
#          LLHOsigEps
#
#-------------------------------------------------------------------------------

#' portion of loglikelihood containing sigEps for haul-out model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param alpha AR1 autocorrelation parameter on logit scale
#' @param sigEps variance of autocorrelated random effects
#' @param eps list with autocorrelated random effects per animal
#' @param tdif list with sequential time differences within animal
#' @ny the number of observations per animal
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export

LLHOsigEps = function(alpha, eps, sigEps, tdif, ny)
{ 
	sumall = 0
	for(j in 1:length(eps)) {
		epsj = eps[[j]]
		tdifj = tdif[[j]]
		nj = length(epsj)
		sumall = sumall + sum(dnorm(epsj[2:nj], alpha^tdifj*epsj[1:(nj - 1)], 
			sigEps, log = TRUE)) + 
			dnorm(epsj[1], 0, sigEps/(1 - alpha^2), log = TRUE)
	}
	sumall - ny*log(1 - 2*pnorm(-6, 0, sigEps))

}
#-------------------------------------------------------------------------------
#
#          LLHO
#
#-------------------------------------------------------------------------------

#' evaluation of whole loglikelihood for haul-out model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param phi variance parameter of beta distribution
#' @param beta vector of regression parameters on logit scale
#' @param alpha AR1 autocorrelation parameter on logit scale
#' @param sigGam variance of random effects per animal
#' @param sigEps variance of autocorrelated random effects
#' @param gam random effects (intercepts) for animals
#' @param eps list with autocorrelated random effects per animal
#' @param y response variable as a proportion (0 < y < 1) (strict inequalities)
#' @param X  linear model design matrix
#' @param Zi vector of integers that indicate which random effect per observation
#' @param tdif list with sequential time differences within animal
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export

LLHO = function(phi, kappa, beta, alpha, sigGam, sigEps, gam, eps, 
	y, X, Zi, tdif)
{
	Xb = X %*% beta + gam[Zi] + unlist(eps)
	p = exp(Xb)/(1 + exp(Xb))
	sumall = 0
	for(j in 1:length(eps)) {
		epsj = eps[[j]]
		tdifj = tdif[[j]]
		nj = length(epsj)
		sumall = sumall + sum(dnorm(epsj[2:nj], alpha^tdifj*epsj[1:(nj - 1)], 
			sigEps, log = TRUE)) + 
			dnorm(epsj[1], 0, sigEps/(1 - alpha^2), log = TRUE)
	}
	LLbern = sum(dbern(y, p, log = TRUE))
	LLgamma = sum(dnorm(gam, 0, sigGam, log = TRUE)) -
			length(gam)*log(1 - 2*pnorm(-6, 0, sigGam))
	LLeps = sumall - length(y)*log(1 - 2*pnorm(-6, 0, sigEps)) 
  list(LLbern = LLbern, LLgamma = LLgamma, LLeps = LLeps,
		logLike = LLbern + LLgamma + LLeps 
	) 
}

#-------------------------------------------------------------------------------
#
#          Bernoulli distribution function
#
#-------------------------------------------------------------------------------

#' Bernoulli distribution function
#'
#' Bernoulli distribution function
#'
#' @param y response variable (must be 0 or 1)
#' @param p probability parameter
#' @param log true of false for the log of the density
#'
#' @return the (log)likelihood
#'
#' @author Jay Ver Hoef
#' @export

dbern = function(x, p, log = FALSE)
{
	if(log == FALSE) out = dbinom(x, 1, p)
	if(log == TRUE) out = dbinom(x, 1, p, log = TRUE)
	out
}

