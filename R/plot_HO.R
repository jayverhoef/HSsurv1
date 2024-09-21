#------------------------------------------------------
#        plot_HO
#------------------------------------------------------

plot_HO = function(M, dHO, dstk, X, y, path) {
pdf(path, width = 11, height = 8.5)

lenbeta = length(M[['beta']][[1]])
# Acceptance Rates
layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))

	nseals = length(M[['eps_acc_rate']])
	par(mar = c(5,5,1,1))
	plot(1:nseals, M[['eps_acc_rate']], ylim = c(0,1), pch = 19, cex = 2,
		xlab = 'year', ylab = 'eps[i] Acceptance Rates', 
		cex.lab = 2, cex.axis = 1.5)
	lines(c(1,nseals), c(0.21,0.21), lty = 2, lwd = 4)
	lines(c(1,nseals), c(0.55,0.55), lty = 2, lwd = 4)

	plot(1:nseals, M[['gam_acc_rate']], ylim = c(0,1), pch = 19, cex = 2,
		xlab = 'year', ylab = 'gam[i] Acceptance Rates', 
		cex.lab = 2, cex.axis = 1.5)
	lines(c(1,nseals), c(0.21,0.21), lty = 2, lwd = 4)
	lines(c(1,nseals), c(0.55,0.55), lty = 2, lwd = 4)

	plot(1:lenbeta, M[['beta_acc_rate']], ylim = c(0,1), pch = 19, cex = 2,
		xlab = 'index', ylab = 'beta[i] Acceptance Rates', 
		cex.lab = 2, cex.axis = 1.5)
	lines(c(1,nseals), c(0.21,0.21), lty = 2, lwd = 4)
	lines(c(1,nseals), c(0.55,0.55), lty = 2, lwd = 4)

	plot(c(1,2,3), c(M[['alpha_acc_rate']], M[['sigGam_acc_rate']],
		M[['sigEps_acc_rate']]),
		pch = 19, cex = 3, xlim = c(0.5, 3.5), ylim = c(0,1), xaxt = 'n',  
		xlab = 'alpha           sigGam       sigEps', 
		ylab = 'Acceptance Rates', cex.lab = 2, cex.axis = 1.5)
	lines(c(0.5, 3.5), c(0.21,0.21), lty = 2, lwd = 4)
	lines(c(0.5, 3.5), c(0.55,0.55), lty = 2, lwd = 4)

layout(1)

# fixed effects
beta0 = unlist(lapply(M[['beta']], function(x){x[[1]]}))
beta1 = unlist(lapply(M[['beta']], function(x){x[[2]]}))
beta2 = unlist(lapply(M[['beta']], function(x){x[[3]]}))
beta3 = unlist(lapply(M[['beta']], function(x){x[[4]]}))
beta4 = unlist(lapply(M[['beta']], function(x){x[[5]]}))
if(lenbeta > 5) {
beta5 = unlist(lapply(M[['beta']], function(x){x[[6]]}))
beta6 = unlist(lapply(M[['beta']], function(x){x[[7]]}))
}


# trace plots of beta elements
	par(mar = c(0,5,0,0))
	layout(matrix(1:7, ncol = 1, byrow = TRUE), 
		heights = c(rep(1, times = 6),1.7))
	if(lenbeta <= 5)
		layout(matrix(1:5, ncol = 1, byrow = TRUE), 
			heights = c(rep(1, times = 4),1.5))
	plot(beta0, type = 'l', xaxt = 'n')
	plot(beta1, type = 'l', xaxt = 'n')
	plot(beta2, type = 'l', xaxt = 'n')
	plot(beta3, type = 'l', xaxt = 'n')
	if(lenbeta <= 5) {
		par(mar = c(5,5,0,0))
		plot(beta4, type = 'l')
	}
	if(lenbeta > 5) {
		plot(beta4, type = 'l', xaxt = 'n')
		plot(beta5, type = 'l', xaxt = 'n')
		par(mar = c(5,5,0,0))
		plot(beta6, type = 'l')
	}
	layout(1)

# beta - date effects

  par(mar = c(5,5,0,0))
  plot(c(1,80),c(0,1), type = 'n',
    ylab = 'Haul-Out Probability', xlab = "Days since 15 July",
    cex.lab = 2, cex.axis = 1.5)
  x = ((1:75) - 30)/30
  HOvec = matrix(0, nrow = length(beta1), ncol = 75)
  for(k in 1:length(beta1)) {
    Xb = beta0[k] + mean(unlist(M[['gam']][k])) + beta1[k]*x + beta2[k]*x^2
    lines(1:75, exp(Xb)/(1+exp(Xb)), col = rgb(0,0,0,.03))
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  lines(1:75, apply(HOvec,2,mean), lwd = 4, col = 'blue')
  lines(1:75, apply(HOvec,2,quantile, probs = .05), 
    lty = 2, lwd = 4, col = 'blue')
  lines(1:75, apply(HOvec,2,quantile, probs = .95), 
    lty = 2, lwd = 4, col = 'blue')
  hist((dHO$dy-as.POSIXlt("2004-07-15")$yday), add = TRUE, freq = FALSE,
    col = rgb(.8,.1,.1,.4))
  hist(dstk$yday-as.POSIXlt(as.POSIXct("2005-07-15"))$yday, freq = FALSE, 
    add = TRUE, col = rgb(.1,.8,.1,.4))
  legend(50, 1, legend = c('Haulout histogram','Survey histogram'), 
		lty = c(1,1), lwd = c(8,8), col = c(rgb(.8,.1,.1,.4), rgb(.1,.8,.1,.4)))

# beta - tide effects (or, hours from solar noon for glacial haulout)

	xlabel = "Hours from Low Tide"
	if(lenbeta < 7) xlabel = "Hours from Solar Noon"
  par(mar = c(5,5,0,0))
	xlimit = c(-5,5)
	if(lenbeta < 7) xlimit = c(-12,12)
  plot(xlimit, c(0,1), type = 'n',
    ylab = 'Haul-Out Probability', xlab = xlabel,
    cex.lab = 2, cex.axis = 1.5)
  x = (-50:50)/25
  xtran = 2.5*x
	if(lenbeta < 7) xtran = 6*x
  HOvec = matrix(0, nrow = length(beta3), ncol = length(x))
  for(k in 1:length(beta3)) {
    Xb = beta0[k] + mean(unlist(M[['gam']][k])) + beta3[k]*x + beta4[k]*x^2
    lines(xtran, exp(Xb)/(1+exp(Xb)), col = rgb(0,0,0,.03))
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  lines(xtran, apply(HOvec,2,mean), lwd = 4, col = 'blue')
  lines(xtran, apply(HOvec,2,quantile, probs = .05), lty = 2, 
    lwd = 4, col = 'blue')
  lines(xtran, apply(HOvec,2,quantile, probs = .95), lty = 2, 
    lwd = 4, col = 'blue')
	if(lenbeta == 7) {  
		hist(dHO$minutes_from_low/60, add = TRUE, freq = FALSE,
			col = rgb(.8,.1,.1,.4))
		hist(as.numeric(dstk$minutes_from_low)/60, add = TRUE, freq = FALSE, 
			na.rm = TRUE, col = rgb(.1,.8,.1,.4))
		legend(1, 1, legend = c('Haulout histogram','Survey histogram'), 
			lty = c(1,1), lwd = c(8,8), col = c(rgb(.8,.1,.1,.4), rgb(.1,.8,.1,.4)))
	}
	if(lenbeta < 7) {  
		hist(dHO$hrstd*6, add = TRUE, freq = FALSE,
			col = rgb(.8,.1,.1,.4))
		hist(dstk$hrstd*6, add = TRUE, freq = FALSE, 
			na.rm = TRUE, col = rgb(.1,.8,.1,.4))
		legend(1, 1, legend = c('Haulout histogram','Survey histogram'), 
			lty = c(1,1), lwd = c(8,8), col = c(rgb(.8,.1,.1,.4), rgb(.1,.8,.1,.4)))
	}
	
# beta - hour-of-day effect

	if(lenbeta ==7 ) {
		par(mar = c(5,5,0,0))
		plot(c(-12,12),c(0,1), type = 'n',
			ylab = 'Haul-Out Probability', 
			xlab = "Hours from Solar Noon",
			cex.lab = 2, cex.axis = 1.5)
		x = (-60:60)/30
		HOvec = matrix(0, nrow = length(beta1), ncol = length(x))
		for(k in 1:length(beta1)) {
			Xb = beta0[k] + mean(unlist(M[['gam']][k])) + 
				beta5[k]*x + beta6[k]*x^2
			lines((-60:60)/5, exp(Xb)/(1+exp(Xb)), 
				col = rgb(0,0,0,.03))
			HOvec[k,] = exp(Xb)/(1+exp(Xb))
		}
		lines((-60:60)/5, apply(HOvec,2,mean), 
			lwd = 4, col = 'blue')
		lines((-60:60)/5, apply(HOvec,2,quantile, probs = .05), 
			lty = 2, lwd = 4, col = 'blue')
		lines((-60:60)/5, apply(HOvec,2,quantile, probs = .95), 
			lty = 2, lwd = 4, col = 'blue')
		hist(dHO$solhr - 12, add = TRUE, freq = FALSE,
			col = rgb(.8,.1,.1,.4))
		hist(dstk$hr - 12, add = TRUE, freq = FALSE,
			col = rgb(.1,.8,.1,.4))
		legend(1, 1, 
			legend = c('Haulout histogram','Survey histogram'), 
			lty = c(1,1), lwd = c(8,8), 
			col = c(rgb(.8,.1,.1,.4), rgb(.1,.8,.1,.4)))
}

# posterior densities for animals random intercepts
  par(mar = c(5,5,.1,.1))
  bw = 1
  maxden = 0
  maxx = 0
  minx = 0
  for(i in 1:length(M[['gam']][[1]])) {
    gam0 = unlist(lapply(M[['gam']], function(x){x[[i]]}))
    maxden = max(maxden, max(density(gam0, bw = bw)$y))
    maxx = max(maxx, max(density(gam0, bw = bw)$x))
    minx = min(minx, min(density(gam0, bw = bw)$x))
  }
  plot(c(minx,maxx),c(0,maxden*1.1), type = 'n', main = '',
    xlab = 'Random Effect Value', ylab = 'Posterior Density',
    cex.lab = 2, cex.axis = 1.5)
  for(i in 1:length(M[['gam']][[1]])) {
    gam0 = unlist(lapply(M[['gam']], function(x){x[[i]]}))
    lines(density(gam0, bw = bw), col = rgb(0,0,1,.4), lwd = 3)
  }

# posteriors of various variance parameters
  layout(matrix(1:4, nrow = 2, byrow = TRUE))
  #posterior density of hourly Haul-out autocorrelation
  par(mar = c(5,5,5,1))
    plot(density(M[['alpha']],bw = .01),
    ylab = 'Posterior Density', xlab = 'Lag-1 Autocorrelation',
    main = 'Haul-out Model', cex.main = 2,
    cex.lab = 2, cex.axis = 1.5, lwd = 3)
  #posterior density of variance parameters
  par(mar = c(5,5,5,1))
  plot(density(M[['sigGam']], bw = mean(M[['sigGam']])/10),
    ylab = 'Posterior Density', xlab = 'sigGam Parameter',
    main = 'Haul-out Model', cex.main = 2,
    cex.lab = 2, cex.axis = 1.5, lwd = 3)
  par(mar = c(5,5,5,1))
  plot(density(M[['sigEps']], bw = mean(M[['sigEps']])/10),
    ylab = 'Posterior Density', xlab = 'sigEps Parameter',
    main = 'Haul-out Model', cex.main = 2,
    cex.lab = 2, cex.axis = 1.5, lwd = 3)
  layout(1)

# trace plots of variance parameters, plus overall loglikelihood
	layout(matrix(1:4, ncol = 2, byrow = TRUE))
	plot(M[['alpha']], type = 'l')
	plot(M[['sigGam']], type = 'l')
	plot(M[['sigEps']], type = 'l')
	plot(M[['logLike']], type = 'l')

# trace plots of the first 6 epsilon values
	layout(matrix(1:6, ncol = 2, byrow = TRUE))
	plot(M[['eps']][[1]], type = 'l')
	plot(M[['eps']][[2]], type = 'l')
	plot(M[['eps']][[3]], type = 'l')
	plot(M[['eps']][[4]], type = 'l')
	plot(M[['eps']][[5]], type = 'l')
	plot(M[['eps']][[6]], type = 'l')

# autocorrelation and partial autocorrelation plots for epsilon
	layout(matrix(1:6, ncol = 2, byrow = TRUE))
	acf(M[['eps']][[1]])
	pacf(M[['eps']][[1]])
	acf(M[['eps']][[2]])
	pacf(M[['eps']][[2]])
	acf(M[['eps']][[3]])
	pacf(M[['eps']][[3]])
	layout(1)

	layout(matrix(1:6, ncol = 2, byrow = TRUE))
	acf(M[['eps']][[4]])
	pacf(M[['eps']][[4]])
	acf(M[['eps']][[5]])
	pacf(M[['eps']][[5]])
	acf(M[['eps']][[6]])
	pacf(M[['eps']][[6]])
	layout(1)
	
	# last MCMC fitted mean (conditional on random effects) versus observed data
	# vector of integers that identify an animal
	Zi = as.integer(dHO$speno)
	mu = X %*% M[['beta']][[1000]] + M[['gam']][[1000]][Zi] + 
		unlist(M[['eps']])
	fit  = exp(mu)/(1 + exp(mu))
	# add some jitter to y to make it easier to see
	par(mar = c(5,5,0,0))
	plot(y+rnorm(length(y),0,0.02), fit, 
		pch = 19, col = rgb(0,0,0,.01), cex = 2,
		xlab = 'Observed Value', ylab = 'Fitted Value',
		cex.lab = 2, cex.axis = 1.5)

	dev.off()

}

