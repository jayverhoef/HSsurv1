#-------------------------------------------------------------------------------
#
#          MCMCabu
#
#-------------------------------------------------------------------------------

#' MCMC sampling for abundance
#'
#' MCMC sampling for abundance

#' @param niter the number of MCMC iterations
#' @param thin thinning of the MCMC chain
#' @param HObeta beta (explanatory variable) MCMC fits from haul-out model. These are for the design matrix created with model.matrix(count ~ daystd + I(daystd^2) + tdstd + I(tdstd^2) + hrstd + I(hrstd^2)
#' @param gam random effects for animals from the MCMC fits from haul-out model.
#' @param beta0 overall intercept in linear model on logit scale
#' @param beta1 linear regression coefficient on hour-of-summer
#' @param beta2 quadratic regression coefficient on hour-of-summer
#' @param alpha autocorrelation range parameter
#' @param sigGam variance of random effects per animal
#' @param eps list with autocorrelated errors per animal
#' @param y response variable as a proportion (0 < y < 1) (strict inequalities)
#' @param X  linear model design matrix
#' @param Z random effects design matrix
#' @param tdif list with sequential time differences within animal
#'
#' @return a list of the MCMC chains values. Only the last iteration of the autocorrelated errors eps is kept.
#'
#' @author Jay Ver Hoef
#' @export

MCMCabu <- function(niter = 1000, thin = 1, HObeta, gam, co, vhat = NULL,
  Xt, A, Nmat, tau, del, REday, sigmaAR1, beta_noflir = NULL,
  prday_sd = 1, maxramult = 4, REday_tune, sigAR1_tune)
{
  
  cat("\n")
    
  # ----------------------------------------------
  #          SETUP
  # ----------------------------------------------

	M <- vector("list", 7)
  M[[1]] = list()
	M[[2]] = list()
	M[[3]] = list()
  del.accept = matrix(0, nrow = dim(del)[1], ncol = dim(del)[2])
	REday.accept = rep(0, times = length(co))
	sAR1.accept = 0
  ikeep = 0

  # ----------------------------------------------
  #          ITERATIONS WITH THINNING
  # ----------------------------------------------

	# round all counts.  Glacial surveys may have decimals, and binomial and betabinomial
	# distributions need integers
	co = round(co, 0)
	# if there are no glacial samples with variance estimates, set all vhat to NA
	if(is.null(vhat)) vhat = rep(NA, times = length(co))

	for(iter in 1:niter) {
			# randomly choose a beta vector and random intercepts from haulout model
			sampi = sample(1:length(HObeta), 1)
      betai = HObeta[[sampi]]
      gami = gam[[sampi]]
      # this will be nonzero only for stock 1 (Aleutians), where we adjust
      # all counts by the noFLIR effect
      if(!is.null(beta_noflir)) {
				# choose one of the noFLIR effects at random if nonNULL
				beta_nofliri = beta_noflir[[sample(1:length(beta_noflir), 1)]]} else
				# otherwise set to 0
				{beta_nofliri = 0}
      # cycle through sites and years for those 2-D arrays, Nmat and del

      for(i in 1:dim(Nmat)[1]) {
        for(j in 1:dim(Nmat)[2]) {  
					# create a temporally-autocorrelated linear mixed model for the
					# latent abundance values
			# Nmat 
			
					loglamij = tau[i] + del[i,j] 

					# if a sampled year    
          if(any(A[,1] == i & A[,2] == j)){
						# any counts (may be more than 1) from that year
            cnt = co[A[,1] == i & A[,2] == j]
            cvar = vhat[A[,1] == i & A[,2] == j]            
            # create a range of possible abundance values from the maximum
            # count up to maxramult times the maximum of the count
						range = seq.int(max(cnt),maxramult*max(co[A[,1]==i]))
						# if there are more than 100 values in the range, spread the 
						# possibe abundance values out evenly, on the integers, so
						# that there are only 100 values to evaluate.  Otherwise,
						# it takes too long		
						# uncomment to examine example for manuscript if(i == 28) browser()				
						if(length(range) > 100)
							range = ceiling(min(range) + (0:100)*(max(range) - 
								min(range))/100)
            # create linear mixed model for the probability of
            # the chosen count
            lmmij = Xt[A[,1] == i & A[,2] == j, , drop = FALSE] %*% 
							betai + mean(gami) + REday[A[,1] == i & A[,2] == j]
            # create the expected value (as a probability) for the
            # randomly chosen observed count
            p_i = exp(beta_nofliri)*exp(lmmij)/(1 + exp(lmmij))
            # choose which value in the range of abundance values to keep
            # randomly by using normalizing the discrete probability among
            # all the values, and then using an inverse CDF technique
						cntpdf = dmultcount(range, cvec = cnt, pvec = p_i, vhat = cvar,
							lambda = exp(loglamij))
						# turn discrete pdf into cdf
						cntcdf = cumsum(cntpdf)
						# now use inverse cdf technique to sample randomly from 
						# the dmultcount distribution
						Nmat[i,j] = range[min(which(runif(1) < cntcdf))]
          } else {
						# for years without data, create a range of possible abundance
						# values for ith site in jth year, starting with the minimum count in any
						# year, up to maxramult times the maximum count in any year
						range = seq.int(min(co[A[,1]==i]),maxramult*max(co[A[,1]==i]))	
            if(length(range) > 100) 							
							range = ceiling(min(range) + (0:100)*(max(range) - 
								min(range))/100)
            # sample from the predictive distribution using inverse cdf 
            # technique
            Nmat[i,j] = range[min(
              which(
                runif(1) < cumsum(
									# create the cumulative sum after subtracting the maximum
									# value (on the log scale, which is equivalent to dividing
									# by the maximum on the exponential scale) so that values do
									# not get too small.  The maximum values cancel when taking
									# dividing by the normalizing constant
                  exp(dpois(range, lambda = exp(loglamij),
										log = TRUE)  - 
                  max(dpois(range, lambda = exp(loglamij), 
										log = TRUE) ))
                )/
                sum(
                  exp(dpois(range, lambda = exp(loglamij), 
										log = TRUE)  - 
                  max(dpois(range, lambda = exp(loglamij), 
										log = TRUE) ))
                )
              )
            )]
          }
          
      #delta
      
          U <- log(runif(1))
					# create the proposal distribution boundaries based on changes from
					# neighboring values, then propose uniformly within those
					# boundaries
					minset = -4
					if(j > 1) minset = c(minset, del[i,j-1] - log(1.5))
					if(j < dim(Nmat)[2])	minset = c(minset, del[i,j+1]- log(1.5))	
					maxset = 4
					if(j > 1) maxset = c(maxset, del[i,j-1] + log(1.5))
					if(j < dim(Nmat)[2])	maxset = c(maxset, del[i,j+1] + log(1.5))	
					lb = max(minset)
					ub = min(maxset)
					delij.try = lb + runif(1)*(ub - lb)
					# do Metropolis sampling
          if(j > 1 & j < dim(Nmat)[2]) {
            LLdif =  -((delij.try - del[i,j-1])^2 + 
									(del[i,j+1] - delij.try)^2)/(2*sigmaAR1^2) +
								dpois(Nmat[i,j], lambda = exp(tau[i] + delij.try), 
									log = TRUE) +
              ((del[i,j] - del[i,j-1])^2 +
									(del[i,j+1] - del[i,j])^2)/(2*sigmaAR1^2) -
								dpois(Nmat[i,j], lambda = exp(tau[i] + del[i,j]), 
									log = TRUE)		
            if(LLdif > U) {
              del[i,j] <- delij.try
              del.accept[i,j] = del.accept[i,j] + 1
            }
          }
          # only have 1 neighbor if j = 1
          if(j == 1) {								
            LLdif = 
							-((del[i,j+1] - delij.try)^2)/(2*sigmaAR1^2) +
								dpois(Nmat[i,j], lambda = exp(tau[i] +  delij.try),
									log = TRUE) + 
							((del[i,j+1] - del[i,j])^2)/(2*sigmaAR1^2) -
								dpois(Nmat[i,j], lambda = exp(tau[i] + del[i,j]), 
									log = TRUE)
            if(LLdif > U) {
              del[i,j] <- delij.try
              del.accept[i,j] = del.accept[i,j] + 1
            }
          }
          # only have 1 neighbor if j is last year
          if(j == dim(Nmat)[2]) {
            LLdif =
              -((delij.try - del[i,j-1])^2)/(2*sigmaAR1^2) +
								dpois(Nmat[i,j], lambda = exp(tau[i] + delij.try), 
									log = TRUE) +
              ((del[i,j] - del[i,j-1])^2)/(2*sigmaAR1^2) -
								dpois(Nmat[i,j], lambda = exp(tau[i] + del[i,j]), 
									log = TRUE)
            if(LLdif > U) {
              del[i,j] <- delij.try
              del.accept[i,j] = del.accept[i,j] + 1
            }
          }
        }
      }
      

			#sigmaAR1
			
				U <- log(runif(1))			
				sigmaAR1.try <-  min(max(sigmaAR1 + runif(1)*sigAR1_tune - 
					sigAR1_tune/2, 0.00001),0.50)
				f = function(x, sar1, nt){sum(dnorm(x[2:nt], x[1:(nt - 1)], 
					sar1, log = TRUE))}
				LLdif = 
					sum(apply(del, 1, f, sar1 = sigmaAR1.try, nt = dim(Nmat)[2])) - 
					sum(apply(del, 1, f, sar1 = sigmaAR1, nt = dim(Nmat)[2]))
				if(LLdif > U) {
					sigmaAR1 <- sigmaAR1.try
					sAR1.accept = sAR1.accept + 1
				}

			
			# REday
      for(i in 1:length(REday)) {
        if(Nmat[A][i] > 0) {
          U <- log(runif(1))			
				  REdayi.try <-  min(max(REday[i] + runif(1)*REday_tune[i] - 
					REday_tune[i]/2, -5),5)
          pi.try = exp(beta_nofliri)*exp(Xt[i,] %*% betai + mean(gami) + REdayi.try)/
            (1 + exp(Xt[i,] %*% betai + mean(gami) + REdayi.try))
          p_i = exp(beta_nofliri)*exp(Xt[i,] %*% betai + mean(gami) + REday[i])/
            (1 + exp(Xt[i,] %*% betai + mean(gami) + REday[i]))
          if(!is.na(vhat[i])) {
						LLdif = dbebivhat(co[i], Nmat[A][i], pi.try, 
							vhat[i], log = TRUE) +
							dnorm(REdayi.try, 0, prday_sd) - 
						dbebivhat(co[i], Nmat[A][i], p_i, 
							vhat[i], log = TRUE) -
							dnorm(REday[i], 0, prday_sd)
					}
          if(is.na(vhat[i])) {
          LLdif = sum(dbinom(co[i], Nmat[A][i], pi.try, log = TRUE)) +
							dnorm(REdayi.try, 0, prday_sd, log = TRUE) - 
            sum(dbinom(co[i], Nmat[A][i], p_i, log = TRUE)) -
							dnorm(REday[i], 0, prday_sd, log = TRUE)
					}
          if(LLdif > U) {
              REday[i] <- REdayi.try
              REday.accept[i] = REday.accept[i] + 1
          }
        }
 
      }
    		  
      if(iter%%thin == 0) {
        ikeep = ikeep + 1
        cat("\r", "Iter Number: ", iter)
        M[[1]][[ikeep]] <- Nmat
        M[[2]][[ikeep]] = del
        M[[3]][[ikeep]] = REday
        M[[4]][[ikeep]] = sigmaAR1
     }

  }

  M[[5]] = del.accept/niter	
  M[[6]] = REday.accept/niter	
  M[[7]] = sAR1.accept/niter	

  cat("\n")
  names(M) <- c("Nmat", "del", "REday", "sigmaAR1",
		"del_acc_rate", "REday_acc_rate", "sAR1_acc_rate")
    
  M
  
}

# a function that reparameterizes the betabinomial distribution as described in the 
# manuscript
dbebivhat = function(x, size, prob, vhat, log = FALSE) {
	# compute kappa as given in manuscript
	kappa = pmax(0, size*prob*(1-prob)*(size -1)/vhat - 1)
	# in VGAM, dbetabinom is parameterized with rho rather than kappa
	# VGAM documentation says alpha = prob*(1 - rho)/rho and 
	# beta = (1 - prob)(1 - rho)/rho, so that means that
	# kappa = (1 - rho)/rho which implies rho = 1/(kappa + 1)
	rho = 1/(kappa + 1)
	# make sure that rho is between 0 and 1 
	rho = pmax(rho, rep(0.00000001, times = length(rho)))
	rho = pmin(rho,rep(0.9999999999, times = length(rho)))
	# abundance estimates have decimals so change them to nearest integer
  VGAM::dbetabinom(round(x,0), size, prob, rho, log = TRUE)
}

# a function for the product of a beta binomial (if estimated variances 
# are provided) plus binomial (if no estimated variances provided) 
# and Poisson probability distribution function
# for constant Nhat but with multivariate observed counts, 
# probabilities, and Poisson means, so returns a vector of 
# distribution values (on the log scale) for a given Nhat
dmc = function(Nhat, cvec, pvec, vhat, lambda) 
{
	LLsum = 0
	if(any(!is.na(vhat))) LLsum = LLsum + 
		sum(dbebivhat(cvec[!is.na(vhat)], Nhat, pvec[!is.na(vhat)], 
		vhat[!is.na(vhat)], log = TRUE))
	if(any(is.na(vhat))) LLsum = LLsum + 
		sum(dbinom(round(cvec[is.na(vhat)],0), 
			Nhat, pvec[is.na(vhat)], log = TRUE))	 
	LLsum + dpois(Nhat, lambda = lambda, log = TRUE)
}
# same as dmc, but now Nhat can be a vector as well
dmultcount = function(Nhat, cvec, pvec, vhat, lambda)
{
	# use dmc to compute product probability for a range of Nhat values 
	pdfvec = unlist(lapply(Nhat, dmc, cvec = cvec, pvec = pvec, vhat = vhat,
		lambda = lambda))
	# the log of the binomial and Poisson products for a vector of Nhat values
	# and subtract the maximum value so, when we exponentiate, we don't get values
	# that are too small or large.  Here, pdfstd at the maximum of pdfvec will 
	# be 0, so exponential of that will be 1.  We will normalize, so the addition
	# of this constant cancels
	pdfstd = pdfvec - max(pdfvec)
	# now normalize and return - a pdf that sums to 1
	exp(pdfstd)/sum(exp(pdfstd))
}


