#------------------------------------------------------
#        plot_ana
#------------------------------------------------------

plot_ana = function(M, dstk_gt0, path) {
pdf(path, width = 11, height = 8.5)

nyr = dim(M$Nmat[[1]])[2]
# Acceptance Rates
layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))

	nsites = dim(M[['del_acc_rate']])[1]
	par(mar = c(5,5,1,1))
	plot(c(1,nyr), c(0,1), ylim = c(0,1), type = 'n', 
		xlab = 'year', ylab = 'del[i,j] Acceptance Rates', cex.lab = 2, cex.axis = 1.5)
	lines(c(1,nyr), c(0.21,0.21), lty = 2, lwd = 4)
	lines(c(1,nyr), c(0.55,0.55), lty = 2, lwd = 4)
	for(i in 1:nsites)
		lines(c(1:nyr), M[['del_acc_rate']][i,], lwd = 2)

	plot(c(1:length(M[['REday_acc_rate']])), M[['REday_acc_rate']], 
		type = 'l', lwd = 2, ylim = c(0,1),  
		xlab = 'Index', ylab = 'REday Acceptance Rates', cex.lab = 2, cex.axis = 1.5)
	lines(c(1,length(M[['REday_acc_rate']])), c(0.21,0.21), lty = 2, lwd = 4)
	lines(c(1,length(M[['REday_acc_rate']])), c(0.55,0.55), lty = 2, lwd = 4)

	plot(1, M[['sAR1_acc_rate']],
		pch = 19, cex = 3, xlim = c(0.5, 2.5), ylim = c(0,1), xaxt = 'n',  
		xlab = 'sigmaAR1', 
		ylab = 'Acceptance Rates', cex.lab = 2, cex.axis = 1.5)
	lines(c(0.5, 2.5), c(0.21,0.21), lty = 2, lwd = 4)
	lines(c(0.5, 2.5), c(0.55,0.55), lty = 2, lwd = 4)

	plot(density(unlist(M[['sigmaAR1']])),
		ylab = 'Posterior Density', xlab = 'Delta Autocorrelated Error Variance',
		cex.lab = 2, cex.axis = 1.5, lwd = 3, main = '')


layout(1)


# Densities for each of the random effects for haulout by day
	maxy = 0
  for(i in 1:length(M[['REday']][[1]])) {
  REday = unlist(lapply(M[['REday']], 
      function(x){x[[i]]}))
	maxy = max(maxy, density(REday, bw = .2)$y)
	}
  par(mar = c(5,5,5,1))
  logitHOprob = rep(NA, times = length(beta0))
  for(i in 1:length(beta0))
  logitHOprob[i] = beta0[i] + mean(gam[[i]])
  plot(c(-5,3),c(0,1.5*maxy), type = 'n',
    xlab = 'Posterior Random Day Effect for Haulout', ylab = 'Posterior Density',
    main = '', 
    cex.lab = 1.5, cex.axis = 1.2)
  for(i in 1:length(M[['REday']][[1]])) {
  REday = unlist(lapply(M[['REday']], 
      function(x){x[[i]]}))
  lines(density(REday, bw = .2), col = rgb(0,.5,.5,.1), lwd = 3)
  }
  lines(density(logitHOprob, bw = .1), lwd = 3)
  legend(-5,1.5*maxy, legend = (c('Haul-out Model', 'Count Model')),
    col = c(rgb(0,0,0), rgb(0,.5,.5)), lwd = c(5,5),
    cex = 2)
  
# posterior binomial probalities for 15 August, low tide, solar noon
# for each day
	maxy = 0
  for(i in 1:length(M[['REday']][[1]])) {
		REday = unlist(lapply(M[['REday']], 
				function(x){x[[i]]})) + mean(beta0) + mean(unlist(gam))
		REday = exp(REday)/(1 + exp(REday))
		maxy = max(maxy, density(REday, bw = .01)$y)
  }
  par(mar = c(5,5,5,1))
  HOprob = rep(NA, times = length(beta0))
  for(i in 1:length(beta0))
  HOprob[i] = exp(beta0[i] + mean(gam[[i]]))/
		(1+exp(beta0[i] + mean(gam[[i]])))
  plot(c(0,1),c(0,1.5*maxy), type = 'n',
    xlab = 'Haul-out Probability', ylab = 'Posterior Density',
    main = '', 
    cex.lab = 2, cex.axis = 1.5)
  for(i in 1:length(M[['REday']][[1]])) {
  REday = unlist(lapply(M[['REday']], 
      function(x){x[[i]]})) + mean(beta0) + mean(unlist(gam))
  REday = exp(REday)/(1 + exp(REday))
  lines(density(REday, bw = .02), col = rgb(0,.5,.5,.05), lwd = 3)
  }
  lines(density(HOprob, bw = .02), lwd = 3)
  legend(2,1.5*maxy, legend = (c('Haul-out Model', 'Count Model')),
    col = c(rgb(0,0,0), rgb(0,.5,.5)), lwd = c(5,5),
    cex = 2)
  legend(0,1.5*maxy, legend = (c('Haul-out Model', 'Count Model')),
    col = c(rgb(0,0,0), rgb(0,.5,.5)), lwd = c(5,5),
    cex = 2)

# posterior of corrections factors for 15 August, low tide, solar noon
	maxy = 0
  for(i in 1:length(M[['REday']][[1]])) {
		REday = unlist(lapply(M[['REday']], 
				function(x){x[[i]]})) + mean(beta0) + mean(unlist(gam))
		REday = (1 + exp(REday))/exp(REday)
		maxy = max(maxy, density(REday, bw = .15)$y)
  }
  par(mar = c(5,5,5,1))
  plot(c(1,6),c(0,1.5*maxy), type = 'n',
    xlab = 'Standardized Correction Factor', ylab = 'Posterior Density',
    main = '', 
    cex.lab = 2, cex.axis = 1.5)
  for(i in 1:length(M[['REday']][[1]])) {
  REday = unlist(lapply(M[['REday']], 
      function(x){x[[i]]})) + mean(beta0) + mean(unlist(gam))
  REday = (1 + exp(REday))/exp(REday)
  lines(density(REday, bw = .15), col = rgb(0,.5,.5,.05), lwd = 3)
  }
  cf = (1+exp(beta0 + unlist(lapply(gam, mean))))/
		exp(beta0 + unlist(lapply(gam, mean)))
  lines(density(cf, bw = .15), lwd = 3)
  legend(2,1.5*maxy, legend = (c('Haul-out Model', 'Count Model')),
    col = c(rgb(0,0,0), rgb(0,.5,.5)), lwd = c(5,5),
    cex = 2)
  
# MCMC traces of first elements of del and Nmat
	all.cex.lab = 2
	all.cex.axis = 1.5
  layout(matrix(1:2, nrow = 2, byrow = TRUE))

		plot(unlist(lapply(M[['del']], function(x){x[1,1]})), type = 'l',
			ylab = 'First site and year for del', 
			cex.lab = all.cex.lab, cex.axis = all.cex.axis)
		plot(unlist(lapply(M[['Nmat']], function(x){x[1,1]})), type = 'l',
			ylab = 'First site and year for Nmat', 
			cex.lab = all.cex.lab, cex.axis = all.cex.axis)

	layout(1)

	# MCMC delta values by site and year
	npages = ceiling(dim(M[['Nmat']][[1]])[1]/12)
	maxsite = dim(M[['Nmat']][[1]])[1]
	polyids = levels(dstk_gt0$polyid)
	pg = 1
	for(pg in 1:npages) {
		sitestart = (pg-1)*12 + 1
		siteend = pg*12
		nMCMC = length(M$sigmaAR1)
		layout(matrix(1:12, nrow = 3, byrow = TRUE))
		par(mar = c(5,5,5,1))
		for(site in sitestart:min(maxsite,siteend)) {
			allval = unlist(lapply(1:nMCMC, function(x){M$del[[x]][site,]}))
			plot(M$del[[1]][site,], type = 'l', col = rgb(0,0,0,.04), 
			ylim = c(min(allval),max(allval)), 
			main = paste('polyid:', polyids[site]),
			ylab = 'Delta Values', cex.lab = 2, cex.axis = 1.5, cex.main = 2)
		for(j in 1:nMCMC)
			lines(1:nyr, M$del[[j]][site,], type = 'l', col = rgb(0,0,0,.04))
		}
		layout(1)
	}


# MCMC sampling of random day effect for haul-out
	npages = min(10, ceiling(dim(A)[1]/12))
	pg = 1
	i = 532
	for(pg in 1:npages) {
		obsstart = (pg-1)*12 + 1
		obsend = pg*12
		layout(matrix(1:12, nrow = 4, ncol = 3, byrow = TRUE))
		for(i in obsstart:obsend) {
			par(mar = c(5,5,1,1))
			if(i <= dim(A)[1] & i > 12*floor(dim(A)[1]/12)) {
				plot(unlist(lapply(M[['REday']], function(x){x[i]})), type = 'l',
					main = paste(polyid_gt0[A[i,1]], A[i,2] + 1995, 
						"count:", dstk_gt0$count[i]),
					ylab = 'REday Value',  xlab = 'MCMC iteration', 
					cex.lab = 1.5, cex.axis = 1.2)
			}
		}
		layout(1)
	}
	
	# Site by site MCMC sampling of abundance estimates
	npages = ceiling(dim(M[['Nmat']][[1]])[1]/12)
	maxsite = dim(M[['Nmat']][[1]])[1]
	polyids = levels(dstk_gt0$polyid)
	pg = 1
	site = 29
	for(pg in 1:npages) {
		sitestart = (pg-1)*12 + 1
		siteend = pg*12
			layout(matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE))
			for(site in sitestart:min(maxsite,siteend)) {
				# condition on site, matrix of year and MCMC value
				popsite = NULL
				cnts = dstk_gt0[dstk_gt0$polyid == 
					levels(dstk_gt0$polyid)[site],c('yr','count')]
				for(i in 1:length(M[['Nmat']])){
					tmp = M[['Nmat']][[i]][site,]
					popsite = rbind(popsite, tmp)
				}
				par(mar = c(5,5,5,1))
				plot(c(1,dim(M[['Nmat']][[1]])[2]), 
					c(min(popsite,cnts[,2]),max(popsite,cnts[,2])), type = 'n',
					xaxt = 'n', ylab = 'Abundance Estimate', cex.main = 2,
					xlab = '', main = paste('polyid:', polyids[site]), cex.lab = 2, cex.axis = 1.5,)
				axis(1, at = c(5,10,15,20,25), labels = c(2000, 2005, 2010, 2015, 2020))
				for(i in 1:length(M[['Nmat']]))
					lines(1:dim(M[['Nmat']][[1]])[2],popsite[i,], lty = 1, lwd = 2, col = rgb(0,0,0,.01))
				points(cnts[,1]-1995, cnts[,2], pch = 19, cex = 1.5, col = 'blue')
				}
	layout(1)
	}


# Abundance estimates with 95% intervals

  pop = NULL
  stockname = as.character(dstk_gt0$stockname[1])
  stockid = as.character(dstk_gt0$stockid[1])
  for(i in 1:length(M[['Nmat']])){
  junk = M[['Nmat']][[i]]
    pop = rbind(pop,apply(junk,2,sum))
  }
  bot = apply(pop,2,quantile, prob = .025)
  top = apply(pop,2,quantile, prob = .975)
  par(mar = c(5,5,5,1))
  plot(c(1,dim(M[['del']][[1]])[2]), c(min(bot),max(top)), type = 'n',
    xaxt = 'n', ylab = 'Estimated Abundance', cex.main = 1.5,
    xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
    main = paste0('Stock ', stockid, ': ',stockname))
  axis(1, at = c(5,10,15,20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
   cex.axis = 1.5)
  points(apply(pop,2,mean), pch = 19, cex = 2)
  for(i in 1:dim(M[['del']][[1]])[2])
    lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

  # moving trailing 8-yr trend estimates linear
  maxi = nyr
  trendlen = 8
  linTrendMat = NULL
  for(i in 1:(maxi - trendlen + 1))
    linTrendMat = cbind(linTrendMat,
      apply(pop,1,function(v){coef(lm(y~x, 
        data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))
  bot = apply(linTrendMat,2,quantile, prob = .025)
  top = apply(linTrendMat,2,quantile, prob = .975)
  par(mar = c(5,5,5,1))
    plot(c(1,length(top)), c(min(bot),max(top)), type = 'n',
      xaxt = 'n', ylab = 'Trend (Seals/Year)', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
      main = paste0('Trailing ', trendlen, '-Year Trend by Year'))
    axis(1, at = c(1,5,10,15,20), labels = c(2003, 2007, 2012, 2017, 2022),
     cex.axis = 1.5)
    lines(c(1,22),c(0,0), lty = 2, lwd = 3, col ='red')
    points(apply(linTrendMat,2,mean), pch = 19, cex = 2)
    for(i in 1:length(top))
      lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

# moving trailing 8-yr trend estimates multiplicative        
  maxi = nyr
  trendlen = 8
  propTrendMat = NULL
  for(i in 1:(maxi - trendlen + 1))
    propTrendMat = cbind(propTrendMat,
      100*(exp(apply(pop,1,function(v){coef(lm(I(log(y))~x, 
        data.frame(x=1:8,y = v[i:(i + trendlen - 1)])))[2]}))-1))
  bot = apply(propTrendMat,2,quantile, prob = .025)
  top = apply(propTrendMat,2,quantile, prob = .975)
  par(mar = c(5,5,5,1))
    plot(c(1,length(top)), c(min(bot),max(top)), type = 'n',
      xaxt = 'n', ylab = 'Trend (%/Year)', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
      main = paste0('Trailing ', trendlen, '-Year Trend by Year'))
    axis(1, at = c(1,5,10,15,20), labels = c(2003, 2007, 2012, 2017, 2022),
     cex.axis = 1.5)
    lines(c(1,22),c(0,0), lty = 2, lwd = 3, col ='red')
    points(apply(propTrendMat,2,mean), pch = 19, cex = 2)
    for(i in 1:length(top))
      lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

# CV by year

  pop = NULL
  stockname = as.character(dstk_gt0$stockname[1])
  stockid = as.character(dstk_gt0$stockid[1])
  for(i in 1:length(M[['Nmat']])){
  tmpp = M[['Nmat']][[i]]
    pop = rbind(pop,apply(tmpp,2,sum))
  }
  CV = sqrt(apply(pop,2,var))/apply(pop,2,mean)
  plot(1:length(CV), CV, 
      xaxt = 'n', ylab = 'CV', cex.main = 1.5,
      xlab = 'Year', cex.lab = 2, cex.axis = 1.5,
      main = 'CV by Year', type = 'l', lwd = 3)
  axis(1, at = c(5,10,15,20, 25), labels = c(2000, 2005, 2010, 2015, 2020),
   cex.axis = 1.5)

dev.off()
}

