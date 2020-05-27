##### 
##### 
##### 
##### 
##### 
# testing functions for the covariance of mixture data set

# tests the fitLKINLAStandard function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
# thetas: originally was c(.1, 3), but 3 is too large
testLKINLAModelMixture = function(seed=1, nLayer=3, nx=20, ny=nx, assumeMeanZero=TRUE, nu=1, 
                                  nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=14, testCovs=FALSE, 
                                  printVerboseTimings=FALSE, latInfo=NULL, n=900, thetas=NULL, 
                                  testfrac=1/9, plotNameRoot="", sigma2=.1^2, useKenya=FALSE, 
                                  effRangeRange=NULL, urbanOverSamplefrac=0, nHyperSamples=1000, 
                                  intStrategy="ccd", strategy="gaussian", separateRanges=FALSE, 
                                  leaveOutRegion=TRUE, gscratch=FALSE, 
                                  savePrecomputationResults=FALSE, loadPrecomputationResults=FALSE, 
                                  precomputationFileNameRoot="precomputationResults") {
  set.seed(seed)
  clusterEffect=TRUE
  
  if(useKenya)
    distanceBreaks = seq(0, 300, l=20)
  else
    distanceBreaks = seq(0, 0.5, l=20)
  
  # set plotNameRoot
  if(length(NC) == 1) {
    if(separateRanges)
      NC = c(14, 126) # by default, use two layers with the finest layer having resolution equal to 10km
  }
  if(separateRanges)
    nLayer = length(NC)
  
  # set plotNameRoot
  # > 2/.08 * 5
  # [1] 125
  # > 2/.8 * 5
  # [1] 12.5
  ncText = ""
  if(length(NC) == 1) {
    if(separateRanges)
      ncText = "_NC14_126"
    else {
      ncText = paste0("_NC", NC)
    }
  } else {
    tempText = do.call("paste0", as.list(c(NC[1], paste0("_", NC[-1]))))
    ncText = paste0("_NC", tempText)
  }
  plotNameRoot = paste0(plotNameRoot, "_L", nLayer, ncText, "_sepRange", separateRanges, "_n", n, "_nu", nu, "_nugV", 
                        round(sigma2, 2), "_Kenya", useKenya, "_noInt", assumeMeanZero, "_urbOversamp", round(urbanOverSamplefrac, 4))
  
  # set true parameter values
  if(useKenya) {
    if(is.null(thetas))
      thetas=c(.08, .8) * (1000/2) / sqrt(8)
  } else {
    if(is.null(thetas))
      thetas=c(.08, .8) / sqrt(8)
  }
  rho = 1
  effectiveRange = thetas * sqrt(8)
  
  # load data set if necessary
  if(is.null(n)) {
    out = load("mixtureDataSet.RData")
  } else {
    spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
        0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
    mixtureCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
        0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
    nTest = round(testfrac * n)
    if(leaveOutRegion) {
      simulationData = getSimulationDataSetsGivenCovarianceTest(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                                nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                                saveDataSetPlot=FALSE, doPredGrid=TRUE)
    } else {
      simulationData = getSimulationDataSetsGivenCovariance(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                            nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                            saveDataSetPlot=FALSE, useKenyaLocations=useKenya, urbanOverSamplefrac=urbanOverSamplefrac)
    }
  }
  coords = cbind(simulationData$xTrain[,1], simulationData$yTrain[,1])
  ys = simulationData$zTrain[,1]
  
  # generate lattice and simulate observations
  # coords = matrix(runif(2*n), ncol=2)
  if(useKenya) {
    xRangeDat = simulationData$xRange
    yRangeDat = simulationData$yRange
    # if(is.null(effRangeRange))
    #   effRangeRange=exp(c(-6, 7))
  } else {
    xRangeDat = c(-1, 1)
    yRangeDat = c(-1, 1)
  }
  if(is.null(latInfo))
    latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  
  AObs = makeA(coords, latInfo)
  # Q = makeQ(kappa=kappa, rho=rho, latInfo, alphas=alphas, normalized=normalize, fastNormalize=fastNormalize) 
  # L = as.matrix(t(chol(solve(Q))))
  # zsims = matrix(rnorm(nrow(Q)), ncol=1)
  # fieldSims = L %*% zsims
  # ys = as.numeric(AObs %*% fieldSims) + 1 # add a constant unit mean term to be estimated by INLA
  # # ys = 1 + as.numeric(AObs %*% fieldSims) + coords[,1] # x-valued mean term to be estimated by INLA
  # errs = rnorm(n, sd=sqrt(sigma2))
  # ys = ys + errs
  
  # plot the observations
  pdf(file=paste0("Figures/mixtureLKINLAObservations", plotNameRoot, ".pdf"), width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid, and add testing points
  mx = 100
  my = 100
  predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  if(useKenya) {
    # remove grid points outside of Kenya national boundaries
    load("../U5MR/adminMapData.RData")
    polys = adm0@polygons
    kenyaPoly = polys[[1]]@Polygons[[77]]@coords
    kenyaPolyProj = projKenya(kenyaPoly)
    inKenya = in.poly(predPts, kenyaPolyProj)
    predPts = predPts[inKenya,]
    
    # add other testing locations to matrix of prediction locations and remember which 
    predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
    predPts = rbind(predPts, cbind(simulationData$xTestRural[,1], simulationData$yTestRural[,1]))
    predPts = rbind(predPts, cbind(simulationData$xTestUrban[,1], simulationData$yTestUrban[,1]))
    plotGridI = 1:sum(inKenya)
    overallTestI = simulationData$overallTestI
    ruralTestI = simulationData$ruralTestI
    urbanTestI = simulationData$urbanTestI
    gridTestI = (max(urbanTestI) + 1):(max(urbanTestI) + length(simulationData$xGrid))
  } else {
    predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
  }
  predPts = rbind(predPts, cbind(simulationData$xGrid, simulationData$yGrid))
  ysTest = c(simulationData$zTest[,1], simulationData$zTestRural[,1], simulationData$zTestUrban[,1], simulationData$zGrid[,1])
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  # priorPar = getPrior(.1, .1, 10)
  # generate hyperparameters for pc priors
  # median effective range is .4 or 200 for kenya data (a fifth of the spatial domain diameter), median spatial variance is 1
  if(!useKenya)
    priorPar = getPCPrior(.4, .5, 1, nLayer=nLayer, separateRanges=separateRanges, latticeInfo=latInfo)
  else
    priorPar = getPCPrior(200, .5, 1, nLayer=nLayer, separateRanges=separateRanges, latticeInfo=latInfo)
  X = matrix(rep(1, nrow(coords)), ncol=1)
  # X = matrix(coords[,1], ncol=1)
  XPred = matrix(rep(1, nrow(predPts)), ncol=1)
  
  # add linear terms in lat/lon to covariate matrices if requested
  if(testCovs) {
    X = cbind(X, coords)
    XPred = cbind(XPred, predPts)
  }
  
  if(assumeMeanZero) {
    X = NULL
    XPred = NULL
  }
  
  # show priors on effective correlation, marginal variance, and error variance:
  if(!useKenya)
    xs1 = seq(.01, 2, l=500)
  else
    xs1 = seq(.01, 1000, l=500)
  if(!separateRanges) {
    pdf(file=paste0("Figures/mixtureLKINLAPriorEffRange", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
         xlab="Effective Correlation Range", main="Effective Correlation Prior", 
         ylab="Prior Density")
    abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
    dev.off()
  } else {
    for(i in 1:nLayer) {
      if(!useKenya)
        xs1 = seq(.01, 2, l=500)
      else
        xs1 = seq(.01, 1000, l=500)
      
      if(i == nLayer && useKenya)
        xs1 = seq(1, 200, l=500)
      else if(!useKenya)
        xs1 = seq(.001, .4, l=500)
      pdf(file=paste0("Figures/mixtureLKINLAPriorEff", i, "Range", plotNameRoot, ".pdf"), width=5, height=5)
      plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar[i]), type="l", col="blue", 
           xlab="Effective Correlation Range", main="Effective Correlation Prior", 
           ylab="Prior Density")
      abline(v=qinvexp(.5, rate=priorPar$corScalePar[i]), col="red")
      dev.off()
    }
  }
  
  if(priorPar$priorType == "orig") {
    xs2 = seq(.01, 10.5, l=500)
    pdf(file=paste0("Figures/mixtureLKINLAPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs2, invgamma::dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
         xlab="Marginal Variance", main="Marginal Variance Prior", 
         ylab="Prior Density")
    abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
    abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
    dev.off()
  } else if(priorPar$priorType == "pc") {
    xs2 = seq(.01, 11.5, l=500)
    pdf(file=paste0("Figures/mixtureLKINLAPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs2, dpcvar(xs2, alpha=priorPar$alpha, u=priorPar$u), type="l", col="blue", 
         xlab="Marginal Variance", main="Marginal Variance Prior", 
         ylab="Prior Density")
    abline(v=qpcvar(.1, alpha=priorPar$alpha, u=priorPar$u), col="red")
    abline(v=qpcvar(.9, alpha=priorPar$alpha, u=priorPar$u), col="red")
    abline(v=1, col="green")
    dev.off()
  }
  
  # xs2 = seq(.001, invgamma::qinvgamma(.905, shape=0.1, rate=0.1), l=500)
  # pdf(file="Figures/mixtureLKINLAPriorErrorVar.pdf", width=5, height=5)
  # plot(xs2, invgamma::dinvgamma(xs2, shape=0.1, rate=0.1), type="l", col="blue", 
  #      xlab="Error Variance", main="Error Variance Prior", 
  #      ylab="Prior Density")
  # abline(v=invgamma::qinvgamma(.1, shape=0.1, rate=0.1), col="red")
  # abline(v=invgamma::qinvgamma(.9, shape=0.1, rate=0.1), col="red")
  # dev.off()
  
  xs2 = seq(.01, 1, l=500)
  pdf(file=paste0("Figures/mixtureLKINLAPriorErrorVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(xs2, dpcvar(xs2, alpha=.05, u=1), type="l", col="blue", 
       xlab="Error Variance", main="Error Variance Prior", 
       ylab="Prior Density")
  abline(v=qpcvar(.1, alpha=.05, u=1), col="red")
  abline(v=qpcvar(.9, alpha=.05, u=1), col="red")
  abline(v=sqrt(.1), col="green")
  dev.off()
  
  # browser()
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/mixtureLKINLAPriorAlpha", l, plotNameRoot, ".pdf"), width=5, height=5)
    xs = seq(0, 1, l=500)
    tempYs = dbeta(xs, priorPar$alphaPar[l], sum(priorPar$alphaPar[-l]))
    plot(xs, tempYs, type="l", xlab=TeX(paste0("$\\alpha_", l, "$")), ylab="Density", 
         main=TeX(paste0("Marginal for $\\alpha_", l, "$")), xlim=c(0,1), ylim=c(0, max(tempYs[is.finite(tempYs)])))
    abline(v=qbeta(c(0.025, 0.975), priorPar$alphaPar[l], sum(priorPar$alphaPar[-l])), col="purple", lty=2)
    dev.off()
  }
  
  # prior on covariogram
  alphaVals = t(rdirichlet(100, alpha=priorPar$alphaPar))
  rhoVals = rpcvar(100, alpha=priorPar$alpha, u=priorPar$u)
  if(!separateRanges) {
    kappaVals = sqrt(8)/rinvexp(100, rate=priorPar$corScalePar) * latInfo[[1]]$latWidth
  } else {
    kappaVals = matrix(sqrt(8)/rinvexp(100*3, rate=priorPar$corScalePar) * sapply(latInfo, function(x) {x$latWidth}), nrow=nLayer)
  }
  nuggetVarVals = rpcvar(100, alpha=.05, u=1)
  if(loadPrecomputationResults)
    loadFilename = precomputationFileNameRoot
  else
    loadFilename = ""
  out = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaVals, 
                                     normalize=normalize, fastNormalize=fastNormalize, precomputationsFileNameRoot=loadFilename)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[sortI]
  lowerCov=out$lowerCov[sortI]
  covMat=out$covMat[sortI]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[sortI]
  lowerCor=out$lowerCor[sortI]
  corMat=out$corMat[sortI]
  
  # true correlation and covariance functions
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", distMat=x, smoothness=nu)}
  spatialCovFun = spatialCorFun
  mixtureCovFun = function(x) {
    out = spatialCorFun(x)
    out[x == 0] = 1 + sigma2
    out
  }
  mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sigma2)) }
  
  # plot the covariance an correlation priors
  yRange = range(c(covMean, lowerCov, upperCov, mixtureCovFun(d)))
  pdf(file=paste0("Figures/mixtureLKINLAPriorCov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  lines(d, mixtureCovFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKINLAPriorCor", plotNameRoot, ".pdf"), width=5, height=5)
  yRange = range(c(corMean, lowerCor, upperCor, mixtureCorFun(d)))
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Covariance", 
       ylim=c(0,1))
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  lines(d, mixtureCorFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  # fit the model
  time = system.time(fit <- fitLKINLAStandard2(coords, ys, predCoords=predPts, seed=seed, nLayer=nLayer, NC=NC,
                                               nBuffer=nBuffer, priorPar=priorPar, xObs=X, xPred=XPred, normalize=normalize, 
                                               intStrategy=intStrategy, strategy=strategy, fastNormalize=fastNormalize, 
                                               printVerboseTimings=printVerboseTimings, latInfo=latInfo, effRangeRange=effRangeRange, 
                                               separateRanges=separateRanges, clusterEffect=clusterEffect, 
                                               savePrecomputationResults=savePrecomputationResults, 
                                               loadPrecomputationResults=loadPrecomputationResults, 
                                               precomputationFileNameRoot=precomputationFileNameRoot))
  mod = fit$mod
  preds=fit$preds
  predSDs=fit$sigmas
  latInfo=fit$latInfo
  latWidth=fit$latWidth
  obsPreds=fit$obsPreds
  obsSDs=fit$obsSDs
  coefPreds = fit$coefPreds
  coefSDs = fit$coefSDs
  
  # print out the total time
  print(paste0("Total time: ", time[3]))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  if(!useKenya) {
    inRange = function(pts, rangeShrink=0) {
      inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
      inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
      inX & inY
    }
  } else {
    inRange = function(pts, rangeShrink=0) {
      rep(TRUE, nrow(pts))
    }
  }
  
  # show predictive surface, SD, and data
  
  if(nLayer==1) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=9, height=6)
    par(mfrow=c(2,3))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeSD = range(c(range(predSDs[inRange(predPts)]), coefSDs[[1]][inRange(gridPtsL1)], 
                         coefSDs[[2]][inRange(gridPtsL2)], coefSDs[[3]][inRange(gridPtsL3)]))
    gridPtsL1 = latInfo[[1]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    
    # quilt.plot(coords, obsPreds, main="Prediction mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat, 
    #            zlim=range(predSDs[inRange(predPts)]))
    plot.new()
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==2) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=12, height=6)
    par(mfrow=c(2,4))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", xlim=xRangeDat, ylim=yRangeDat)
    
    quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==3) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=15, height=6)
    par(mfrow=c(2,5), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==4) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=18, height=6)
    par(mfrow=c(2,6), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    gridPtsL4 = latInfo[[4]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefPreds[[4]], main="Basis Coefficient Mean (Layer 4)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefSDs[[4]], main="Basis Coefficient SD (Layer 4)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==5) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=21, height=6)
    par(mfrow=c(2,7), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    gridPtsL4 = latInfo[[4]]$latCoords
    gridPtsL5 = latInfo[[5]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefPreds[[4]], main="Basis Coefficient Mean (Layer 4)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL5[,1], gridPtsL5[,2], coefPreds[[5]], main="Basis Coefficient Mean (Layer 5)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefSDs[[4]], main="Basis Coefficient SD (Layer 4)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL5[,1], gridPtsL5[,2], coefSDs[[5]], main="Basis Coefficient SD (Layer 5)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKINLALeftOutResiduals", plotNameRoot, ".pdf"), width=5, height=5)
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  plot(preds[testIndices], ysTest-preds[testIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted", 
       ylab="Residuals", xlab="Fitted")
  abline(h=0, lty=2)
  dev.off()
  
  if(useKenya) {
    pdf(file=paste0("Figures/mixtureLKINLALeftOutResidualsLabeled", plotNameRoot, ".pdf"), width=5, height=5)
    testIndices = (length(preds) - length(ysTest) + 1):length(preds)
    gridIndices = 
    ylim = range(ysTest-preds[testIndices])
    xlim = range(preds[testIndices])
    plot(preds[testIndices][overallTestI], ysTest[overallTestI]-preds[testIndices][overallTestI], pch=19, cex=.1, col="black", main="Residuals versus fitted", 
         ylab="Residuals", xlab="Fitted", xlim=xlim, ylim=ylim)
    points(preds[testIndices][ruralTestI], ysTest[ruralTestI]-preds[testIndices][ruralTestI], pch=19, cex=.1, col="green")
    points(preds[testIndices][urbanTestI], ysTest[urbanTestI]-preds[testIndices][urbanTestI], pch=19, cex=.1, col="blue")
    points(preds[testIndices][gridTestI], ysTest[gridTestI]-preds[testIndices][gridTestI], pch=19, cex=.1, col="red")
    abline(h=0, lty=2)
    legend("topright", c("Overall", "Rural", "Urban", "Grid"), col=c("black", "green", "blue", "red"), pch=19)
    dev.off()
  }
  
  # calculate true effective range and marginal variance:
  latticeWidth = latInfo[[1]]$latWidth
  if(separateRanges)
    latticeWidth = sapply(latInfo, function(x) {x$latWidth})#
  # marginalVar = rho/(4*pi * kappa^2)
  # marginalVar = getMultiMargVar(kappa, rho, nLayer=nLayer, nu=nu, xRange=xRangeBasis, 
  #                               yRange=yRangeBasis, nx=nx, ny=ny)[1]
  marginalVar = rho
  
  # # plot marginals on interpretable scale (effective range, marginal variance)
  # effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
  # varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  # sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar$`Precision for the Gaussian observations`)
  # covNames = names(mod$marginals.fixed)
  # if(!assumeMeanZero) {
  #   XMarginals = list()
  #   for(i in 1:length(covNames)) {
  #     XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
  #     XMarginals = c(XMarginals, list(XMarginal))
  #   }
  # }
  # 
  # par(mfrow=c(1,1))
  # if(!separateRanges) {
  #   pdf(file=paste0("Figures/mixtureLKINLAEffRange", plotNameRoot, ".pdf"), width=5, height=5)
  #   plot(effRangeMarg, type="l", main="Marginal for effective range")
  #   abline(v=effectiveRange, col="green")
  #   abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  #   dev.off()
  # } else {
  #   for(l in 1:nLayer) {
  #     pdf(file=paste0("Figures/mixtureLKINLAEffRange", "Layer", l, plotNameRoot, ".pdf"), width=5, height=5)
  #     plot(effRangeMarg, type="l", main=paste0("Marginal for effective range (Layer ", l, ")"))
  #     abline(v=effectiveRange, col="green")
  #     abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  #     dev.off()
  #   }
  # }
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  if(!separateRanges) {
    effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
    varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  } else {
    numberThetas = nLayer + 1 + nLayer - 1
    allNames = paste0("Theta", 1:numberThetas, " for field")
    effRangeMargs = list()
    for(i in 1:nLayer) {
      effRangeMargs = c(effRangeMargs, list(inla.tmarginal(exp, mod$marginals.hyperpar[[allNames[i]]])))
    }
    varMarg = inla.tmarginal(exp, mod$marginals.hyperpar[[allNames[nLayer+1]]])
  }
  sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar[[1]])
  covNames = names(mod$marginals.fixed)
  XMarginals = list()
  if(length(covNames) != 0) {
    for(i in 1:length(covNames)) {
      XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
      XMarginals = c(XMarginals, list(XMarginal))
    }
  }
  
  par(mfrow=c(1,1))
  
  if(!separateRanges) {
    pdf(file=paste0("Figures/mixtureLKINLAEffRange", plotNameRoot, ".pdf"), width=5, height=5)
    plot(effRangeMarg, type="l", main="Marginal for effective range")
    abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
    dev.off()
  } else {
    for(i in 1:nLayer) {
      pdf(file=paste0("Figures/mixtureLKINLAEffRange", i, plotNameRoot, ".pdf"), width=5, height=5)
      plot(effRangeMargs[[i]], type="l", main="Marginal for effective range")
      abline(v=inla.qmarginal(c(.025, .975), effRangeMargs[[i]]), col="purple", lty=2)
      dev.off()
    }
  }
  
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file=paste0("Figures/mixtureLKINLAVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(varMarg, type="l", main="Marginal for spatial variance")
  abline(v=marginalVar, col="green")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file=paste0("Figures/mixtureLKINLASigma2", plotNameRoot, ".pdf"), width=5, height=5)
  plot(sigma2Marg, type="l", main="Marginal for error variance")
  abline(v=sigma2, col="green")
  abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
  dev.off()
  
  if(!assumeMeanZero) {
    for(i in 1:length(covNames)) {
      XMarginal = XMarginals[[i]]
      pdf(file=paste0("Figures/mixtureLKINLA", covNames[i], plotNameRoot, ".pdf"), width=5, height=5)
      plot(XMarginal, type="l", main="Marginal for fixed effect")
      abline(v=0, col="green")
      abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
      dev.off()
    }
  }
  
  # pdf(file="Figures/mixtureLKINLARho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  
  # do the same for kappa, rho
  # in order to get distribution for rho, must sample from joint hyperparameters
  kappaMarg = inla.tmarginal(function(x) {sqrt(8)/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  # thetasToRho = function(xs) {
  #   logCor = xs[2]
  #   logVar = xs[3]
  #   kappa = sqrt(8)/exp(logCor) * latticeWidth
  #   sigma2 = exp(logVar)
  #   sigma2 * 4*pi * kappa^2
  # }
  # samples = inla.hyperpar.sample(50000, mod, TRUE)
  # rhos = apply(samples, 1, thetasToRho)
  
  # pdf(file=paste0("Figures/mixtureLKINLAKappa", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  # abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  # dev.off()
  
  # pdf(file="Figures/mixtureLKINLARho.pdf", width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now generate marginals for the alpha parameters. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  out = inla.hyperpar.sample(nHyperSamples, mod, improve.marginals=TRUE)
  if(separateRanges)
    alphaI = (1 + nLayer+1 + 1):(1 + nLayer+1 + nLayer-1)
  else
    alphaI = 4:(3+nLayer-1)
  zSamples = matrix(out[,alphaI], ncol=length(alphaI))
  xSamples = t(matrix(apply(zSamples, 1, multivariateExpit), ncol=length(alphaI)))
  xSamples = rbind(xSamples, 1-colSums(xSamples))
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/mixtureLKINLAAlpha", l, plotNameRoot, ".pdf"), width=5, height=5)
    hist(xSamples[l,], xlab=TeX(paste0("$\\alpha_", l, "$")), main=TeX(paste0("Marginal for $\\alpha_", l, "$")), breaks=100, freq=F, xlim=c(0,1))
    abline(v=mean(xSamples[l,]), col="purple", lty=1)
    abline(v=quantile(probs=c(.025, .975), xSamples[l,]), col="purple", lty=2)
    dev.off()
  }
  
  ## plot covariance and correlation functions
  # first get the true covariance an correlation functions
  # spatialCovFun = function(x) {0.4 * Exp.cov(x, theta=0.1) + 0.6 * Exp.cov(x, theta=3)}
  # mixtureCovFun = function(x) {
  #   out = spatialCovFun(x)[1,]
  #   out[x == 0] = 1 + sqrt(.1)
  #   out
  # }
  # mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sqrt(.1))) }
  # spatialCorFun = function(x) {0.4 * Exp.cov(x, theta=thetas[1]) + 0.6 * Exp.cov(x, theta=thetas[2])}
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", distMat=x, smoothness=nu)}
  spatialCovFun = spatialCorFun
  mixtureCovFun = function(x) {
    out = spatialCorFun(x)
    out[x == 0] = 1 + sigma2
    out
  }
  mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sigma2)) }
  
  # first to transform all the hyperparameter samples to their relevant values
  # 1: error precision
  # 2: log effective range
  # 3: log spatial variance
  # 4-(3 + nLayer - 1): multivariateLogit alpha
  # nuggetVarVals = 1 / out[,1]
  nuggetVarVals = rep(0, ncol(xSamples))
  if(separateRanges) {
    kappaVals = t(sweep(sqrt(8)/exp(out[,2:(nLayer+1)]), 2, sapply(latInfo, function(x) {x$latWidth}), "*"))
    rhoVals = exp(out[,nLayer+2])
  } else {
    kappaVals = sqrt(8)/exp(out[,2]) * latticeWidth
    rhoVals = exp(out[,3])
  }
  alphaMat = xSamples
  
  # compute the covariance function for many different hyperparameter samples
  out = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaMat, precomputationsFileNameRoot=loadFilename)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[sortI]
  lowerCov=out$lowerCov[sortI]
  covMat=out$covMat[sortI,]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[sortI]
  lowerCor=out$lowerCor[sortI]
  corMat=out$corMat[sortI,]
  covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
                 corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  
  # plot the covariance function
  pdf(file=paste0("Figures/mixtureLKINLACov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance")
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  lines(d, mixtureCovFun(d), col="green")
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKINLACor", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Correlation")
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  lines(d, mixtureCorFun(d), col="green")
  dev.off()
  
  yRange = range(c(covMean, lowerCov, upperCov, mixtureCovFun(d)))
  pdf(file=paste0("Figures/mixtureLKINLACov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  lines(d, mixtureCovFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKINLACor", plotNameRoot, ".pdf"), width=5, height=5)
  yRange = range(c(corMean, lowerCor, upperCor, mixtureCorFun(d)))
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Correlation", 
       ylim=c(0,1))
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  lines(d, mixtureCorFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  # get scoring rules
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  leftOutIndices = (length(preds) - length(ysTest) + 1):(length(preds) - length(ysTest) + length(simulationData$zTest[,1]))
  gridIndices = (length(preds) - length(ysTest) + length(simulationData$zTest[,1]) + 1):length(preds)
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  # first calculate scoring rules at grid points
  truth = ysTest[gridIndicesTest]
  est = preds[gridIndices]
  vars = predSDs[gridIndices]^2
  lower = fit$lower[gridIndices]
  upper = fit$upper[gridIndices]
  
  # compute nearest neighbor distances and scores as a function of them
  gridPts = predPts[gridIndices,]
  distMat = rdist(coords, gridPts)
  nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  print("Binned grid scores:")
  gridScoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  print(gridScoringRules$pooledResults)
  print(gridScoringRules$binnedResults)
  
  # now calculate square rules at left out points
  truth = ysTest[leftOutIndicesTest]
  est = preds[leftOutIndices]
  vars = predSDs[leftOutIndices]^2
  lower = fit$lower[leftOutIndices]
  upper = fit$upper[leftOutIndices]
  leftOutScoringRules = getScores(truth, est, vars, lower, upper)
  leftOutScoringRules = data.frame(c(leftOutScoringRules, Time=time[3]))
  print("Binned left out scores:")
  print(leftOutScoringRules)
  
  scoringRules = list(gridScoringRules=gridScoringRules, leftOutScoringRules=leftOutScoringRules)
  # 
  # # get scoring rules
  # truth = ysTest
  # est = preds[testIndices]
  # vars = predSDs[testIndices]^2
  # lower = fit$lower[testIndices]
  # upper = fit$upper[testIndices]
  # 
  # # compute nearest neighbor distances and scores as a function of them
  # testPts = predPts[testIndices,]
  # distMat = rdist(coords, testPts)
  # nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  # print("Binned scores:")
  # scoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  # scoringRules$pooledResults = data.frame(c(scoringRules$pooledResults, Time=time[3]))
  # print(scoringRules$binnedResults)
  
  
  
  # print("Grid scores:")
  # print(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                 distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults)
  # 
  # # plot scores as a function of distance
  # distanceScores = getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                            distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults
  
  # pdf(file=paste0("Figures/mixtureSPDEScoreBias", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Bias, pch=19, col="blue", main="Bias", ylab="Bias", xlab="Nearest neighbor distance (km)")
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreVar", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Var, pch=19, col="blue", main="Variance", ylab="Variance", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Var)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$MSE, pch=19, col="blue", main="MSE", ylab="MSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$MSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreRMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$RMSE, pch=19, col="blue", main="RMSE", ylab="RMSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$RMSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCRPS", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$CRPS, pch=19, col="blue", main="CRPS", ylab="CRPS", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$CRPS)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCvg", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Coverage, pch=19, col="blue", main="80% Coverage", ylab="80% Coverage", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Coverage)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreWidth", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Width, pch=19, col="blue", main="Width", ylab="Width", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Width)))
  # abline(h=0, lty=2)
  # dev.off()
  
  if(!useKenya) {
    # print("Pooled scores:")
    # print(data.frame(c(getScores(truth, est, vars, lower, upper), Time=time[3])))
  } else {
    print("Pooled scores:")
    print(data.frame(c(getScores(truth, est, vars, lower, upper), Time=time[3])))
    print("Overall scores:")
    print(data.frame(c(getScores(truth[overallTestI], est[overallTestI], vars[overallTestI], lower[overallTestI], upper[overallTestI]), Time=time[3])))
    print("Rural scores:")
    print(data.frame(c(getScores(truth[ruralTestI], est[ruralTestI], vars[ruralTestI], lower[ruralTestI], upper[ruralTestI]), Time=time[3])))
    print("Urban scores:")
    print(data.frame(c(getScores(truth[urbanTestI], est[urbanTestI], vars[urbanTestI], lower[urbanTestI], upper[urbanTestI]), Time=time[3])))
    # print("Grid scores:")
    # print(data.frame(c(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI]), Time=time[3])))
  }
  
  # get aggregated predictions
  # A = t(sapply(1:(mx*my), getARow))
  # A = sweep(A, 1, rowSums(A), "/")
  # mx = 100
  # my = 100
  # predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  # testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  # gridIndices = (length(preds) - length(simulationData$xGrid) + 1):length(preds)
  # gridIndicesEst = (length(est) - length(simulationData$xGrid) + 1):length(est)
  
  A = makeNumericalIntegralMat(gridCoords, mx=3, my=3)
  aggregatedPreds = A %*% preds[gridIndices]
  
  truth = A %*% ysTest[gridIndicesTest]
  est = A %*% preds[gridIndices]
  aggregatedPredMat = A %*% fit$predMat[gridIndices,]
  vars = apply(aggregatedPredMat, 1, var)
  sds = apply(aggregatedPredMat, 1, sd)
  lower = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.1)})
  upper = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.9)})
  predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  print("Aggregated prediction summary statistics:")
  print(predictionMatrix)
  
  print("Pooled aggregated scores:")
  pooledAggregatedScores = getScores(truth, est, vars, lower, upper)
  print(pooledAggregatedScores)
  print("Left out region aggregated scores:")
  leftOutAggregatedScores = getScores(truth[5], est[5], vars[5], lower[5], upper[5])
  leftOutAggregatedScores$Var = sds[5]
  names(leftOutAggregatedScores)[2] = "Predictive.SD"
  print(leftOutAggregatedScores)
  print("Included regions aggregated scores:")
  includedAggregatedScores = getScores(truth[-5], est[-5], vars[-5], lower[-5], upper[-5])
  print(includedAggregatedScores)
  
  aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
                                includedAggregatedScores=includedAggregatedScores)
  
  fit$mod = NULL
  if(!gscratch)
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("savedOutput/simulations/mixtureLKINLA", plotNameRoot, ".RData"))
  else
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("/work/johnpai/mixtureLKINLA", plotNameRoot, ".RData"))
  
  invisible(list(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules))
}

# tests the fitLKINLAStandard function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
# thetas: originally was c(.1, 3), but 3 is too large
# testLKINLAModelMixtureMultiple = function(seed=1, nSamples=10, nLayer=3, nx=20, ny=nx, assumeMeanZero=TRUE, nu=1, 
#                                   nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=14, testCovs=FALSE, 
#                                   printVerboseTimings=FALSE, latInfo=NULL, n=900, thetas=NULL, 
#                                   testfrac=.1, plotNameRoot="", sigma2=.1^2, useKenya=FALSE, 
#                                   effRangeRange=NULL, urbanOverSamplefrac=0, 
#                                   intStrategy="ccd", strategy="gaussian", separateRanges=FALSE, 
#                                   leaveOutRegion=TRUE) {
testLKINLAModelMixtureMultiple = function(seed=1, nSamples=100, NC=14, nLayer=3, separateRanges=FALSE, n=900, nu=1, sigma2=0.1^2, 
                                          useKenya=FALSE, assumeMeanZero=TRUE, urbanOverSamplefrac=0, gscratch=TRUE, ...) {
  # set random seeds for each simulation
  set.seed(seed)
  allSeeds = sample(1:1000000, nSamples, replace = FALSE)
  
  # load in the results
  print("Loading in simulation results")
  # set plotNameRoot
  if(length(NC) == 1) {
    if(separateRanges)
      NC = c(14, 126) # by default, use two layers with the finest layer having resolution equal to 10km
  }
  if(separateRanges)
    nLayer = length(NC)
  
  # set plotNameRoot
  # > 2/.08 * 5
  # [1] 125
  # > 2/.8 * 5
  # [1] 12.5
  ncText = ""
  if(length(NC) == 1) {
    if(separateRanges)
      ncText = "_NC14_126"
    else {
      ncText = paste0("_NC", NC)
    }
  } else {
    tempText = do.call("paste0", as.list(c(NC[1], paste0("_", NC[-1]))))
    ncText = paste0("_NC", tempText)
  }
  plotNameRoot = paste0("_L", nLayer, ncText, "_sepRange", separateRanges, "_n", n, "_nu", nu, "_nugV", 
                        round(sigma2, 2), "_Kenya", useKenya, "_noInt", assumeMeanZero, "_urbOversamp", round(urbanOverSamplefrac, 4))
  
  
  # call testLKINLAModelMixture for each simulation requested
  precomputationFileNameRoot = paste0("precomputationMixLKINLA", plotNameRoot)
  temp = function(i) {
    print(paste0("Beginning simulation ", i, "/", nSamples))
    thisPlotNameRoot = paste0("sim", i)
    
    # make sure to save the precomputed results if this is the first run
    savePrecomputationResults = i == 1
    loadPrecomputationResults = i != 1
    do.call("testLKINLAModelMixture", c(list(seed = allSeeds[i], NC=NC, nLayer=nLayer, separateRanges=separateRanges, n=n, nu=nu, sigma2=sigma2, 
                                             useKenya=useKenya, urbanOverSamplefrac=urbanOverSamplefrac, assumeMeanZero=assumeMeanZero, 
                                             plotNameRoot=thisPlotNameRoot, gscratch=gscratch, 
                                             savePrecomputationResults=savePrecomputationResults, 
                                             loadPrecomputationResults=loadPrecomputationResults, 
                                             precomputationFileNameRoot=precomputationFileNameRoot), list(...)))
  }
  sapply(1:nSamples, temp)
  
  # save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("savedOutput/simulations/mixtureLKINLA", plotNameRoot, ".RData"))
  allScoringRulesGrid = list()
  allScoringRulesLeftOut = list()
  allFits = list()
  allCovInfo = list()
  allPredictionMatrices = list()
  allAggregatedScoringRules = list()
  for(i in 1:nSamples) {
    if(!gscratch)
      out = load(paste0("savedOutput/simulations/mixtureLKINLAsim", i, plotNameRoot, ".RData"))
    else
      out = load(paste0("/work/johnpai/mixtureLKINLAsim", i, plotNameRoot, ".RData"))
    allScoringRulesGrid = c(allScoringRulesGrid, list(scoringRules$gridScoringRules))
    allScoringRulesLeftOut = c(allScoringRulesLeftOut, list(scoringRules$leftOutScoringRules))
    allFits = c(allFits, list(fit))
    allCovInfo = c(allCovInfo, list(covInfo))
    allPredictionMatrices = c(allPredictionMatrices, list(predictionMatrix))
    allAggregatedScoringRules = c(allAggregatedScoringRules, list(aggregatedScoringRules))
  }
  
  ##### average results from each simulation
  # pointwise scoring rules
  allPooledScoringRulesGrid = do.call("rbind", lapply(allScoringRulesGrid, function(x) {x$pooledResults}))
  allBinnedScoringRulesGrid = lapply(allScoringRulesGrid, function(x) {x$binnedResults})
  binnedScoringRulesGrid = averageBinnedScores(allBinnedScoringRulesGrid)
  ns = binnedScoringRulesGrid[,2]
  pooledScoringRulesGrid = apply(binnedScoringRulesGrid, 2, function(x) {sum(x * (ns / sum(ns)))})
  pooledScoringRulesGrid = as.data.frame(matrix(pooledScoringRulesGrid, nrow=1))
  names(pooledScoringRulesGrid) = names(binnedScoringRulesGrid)
  
  fullPooledScoringRulesLeftOut = do.call("rbind", allScoringRulesLeftOut)
  pooledScoringRulesLeftOut = colMeans(fullPooledScoringRulesLeftOut)
  
  # covInfo
  # covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
  #                corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  d = allCovInfo[[1]]$d
  covMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$covMean})))
  upperCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCov})))
  lowerCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCov})))
  corMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$corMean})))
  upperCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCor})))
  lowerCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCor})))
  
  # aggregated scoring rules
  # aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
  #                               includedAggregatedScores=includedAggregatedScores)
  # predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  fullPredictionMatrix = do.call("rbind", allPredictionMatrices)
  leftOutIndices = seq(5, nrow(fullPredictionMatrix), by=9)
  leftOutPredictionMatrix = fullPredictionMatrix[leftOutIndices, ]
  leftInPredictionMatrix = fullPredictionMatrix[-leftOutIndices, ]
  leftOutScores = getScores(leftOutPredictionMatrix[,1], leftOutPredictionMatrix[,2], leftOutPredictionMatrix[,3]^2)
  leftInScores = getScores(leftInPredictionMatrix[,1], leftInPredictionMatrix[,2], leftInPredictionMatrix[,3]^2)
  aggregatedScores = getScores(fullPredictionMatrix[,1], fullPredictionMatrix[,2], fullPredictionMatrix[,3]^2)
  
  ##### Save results
  if(!gscratch) {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("savedOutput/simulations/mixtureLKINLAAll_nsim", nSamples, plotNameRoot, ".RData"))
  } else {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("/work/johnpai/mixtureLKINLAAll_nsim", nSamples, plotNameRoot, ".RData"))
  }
}

# tests the fitLKStandard function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKModelMixture = function(seed=548676, nLayer=3, nx=20, ny=nx, nu=1, assumeMeanZero=TRUE, 
                              nBuffer=5, normalize=TRUE, NC=14, testCovs=TRUE, 
                              printVerboseTimings=FALSE, n=900, separatea.wght=FALSE, 
                              plotNameRoot="", doMatern=FALSE, fixNu=FALSE, thetas=c(.08, .8) / sqrt(8), 
                              testfrac=1/9, leaveOutRegion=TRUE, sigma2 = 0.1^2, extraPlotName=plotNameRoot, 
                              gscratch=FALSE) {
  set.seed(seed)
  
  # if(useKenya)
  #   distanceBreaks = seq(0, 300, l=20)
  # else
  distanceBreaks = seq(0, 0.5, l=20)
  
  # set true parameter values
  rho = 1
  effectiveRange = thetas * sqrt(8)
  
  # load data set if necessary
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
  if(is.null(n)) {
    out = load("mixtureDataSet.RData")
  } else {
    mixtureCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
        0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
    nTest = round(testfrac * n)
    if(leaveOutRegion) {
      simulationData = getSimulationDataSetsGivenCovarianceTest(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                                nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                                saveDataSetPlot=FALSE, doPredGrid=TRUE)
    } else {
      simulationData = getSimulationDataSetsGivenCovariance(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                            nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                            saveDataSetPlot=FALSE, useKenyaLocations=useKenya, urbanOverSamplefrac=urbanOverSamplefrac)
    }
  }
  coords = cbind(simulationData$xTrain[,1], simulationData$yTrain[,1])
  ys = simulationData$zTrain[,1]
  
  # generate lattice and simulate observations
  # coords = matrix(runif(2*n), ncol=2)
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  
  # plot the observations
  pdf(file=paste0("Figures/mixtureObservations", extraPlotName, ".pdf"), width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid, and add testing points
  # urban and rural points are NULL, but they are left in the code in case kenya data used later
  xRange=c(-1,1)
  yRange=c(-1,1)
  mx = 100
  my = 100
  predPts = make.surface.grid(list(x=seq(xRange[1], xRange[2], l=mx), y=seq(yRange[1], yRange[2], l=my)))
  predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
  # ysTest = simulationData$zTest[,1]
  predPts = rbind(predPts, cbind(simulationData$xGrid, simulationData$yGrid))
  ysTest = c(simulationData$zTest[,1], simulationData$zTestRural[,1], simulationData$zTestUrban[,1], simulationData$zGrid[,1])
  
  # fit the model
  if(assumeMeanZero)
    m=0
  else
    m=1
  time = system.time(fit <- fitLKStandard(coords, ys, predCoords=predPts, nLayer=nLayer, NC=NC,
                                          nBuffer=nBuffer, normalize=normalize, fixedFunctionArgs=list(m=m), 
                                          xRangeDat=xRange, yRangeDat=yRange, separatea.wght=separatea.wght, 
                                          doMatern=doMatern, fixNu=fixNu))
  mod = fit$mod
  preds = fit$preds
  predSDs = fit$sigmas
  parameterSummaryTable = fit$parameterSummaryTable
  parSim = fit$parSim
  
  # print out the total time
  print(paste0("Total time: ", time[3]))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  inRange = function(pts, rangeShrink=0) {
    inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
    inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
    inX & inY
  }
  
  # show predictive surface, SD, and data
  
  pdf(file=paste0("Figures/mixtureLKPreds", extraPlotName, ".pdf"), width=8, height=8)
  par(mfrow=c(2,2))
  
  # obsInds = 1:n
  # predInds = (n+1):(n+mx*my)
  # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
  colRangeDat = range(c(ys, preds))
  colRangeSD = range(range(predSDs[inRange(predPts)]))
  
  quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  
  plot.new()
  quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
             xlim=xRangeDat, ylim=yRangeDat)
  dev.off()
  
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  leftOutIndices = (length(preds) - length(ysTest) + 1):(length(preds) - length(ysTest) + length(ys))
  gridIndices = (length(preds) - length(ysTest) + length(ys)):length(preds)
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  pdf(file=paste0("Figures/mixtureLKLeftOutResidualsGrid", extraPlotName, ".pdf"), width=5, height=5)
  plot(preds[gridIndices], ysTest[gridIndicesTest]-preds[gridIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted (Grid)", 
       ylab="Residuals", xlab="Fitted")
  abline(h=0, lty=2)
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKLeftOutResidualsLeftOut", extraPlotName, ".pdf"), width=5, height=5)
  plot(preds[leftOutIndices], ysTest[leftOutIndicesTest]-preds[leftOutIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted (Left out)", 
       ylab="Residuals", xlab="Fitted")
  abline(h=0, lty=2)
  dev.off()
  
  # calculate true effective range and marginal variance:
  # marginalVar = rho/(4*pi * kappa^2)
  # marginalVar = getMultiMargVar(kappa, rho, nLayer=nLayer, nu=nu, xRange=xRangeBasis, 
  #                               yRange=yRangeBasis, nx=nx, ny=ny)[1]
  
  # plot estimates on interpretable scale (effective range, marginal variance)
  
  # transform all the hyperparameter samples to their relevant values
  totalVariance = fit$totalVariance
  lambdaVals = fit$lambdaVals
  rhoVals = fit$rhoVals
  nuggetVarVals = fit$nuggetVarVals
  a.wghtVals = fit$a.wghtVals
  alphaVals = fit$alphaVals
  nuVals = fit$nuVals
  
  par(mfrow=c(1,1))
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file=paste0("Figures/mixtureLKVar", extraPlotName, ".pdf"), width=5, height=5)
  hist(rhoVals, main="Estimate of spatial variance", freq=FALSE, breaks=20)
  abline(v=1 + sqrt(.1), col="green")
  abline(v=quantile(probs=c(.025, .975), rhoVals), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file=paste0("Figures/mixtureLKSigma2", extraPlotName, ".pdf"), width=5, height=5)
  hist(nuggetVarVals, main="Estimate of error variance", freq=FALSE, breaks=20)
  abline(v=sigma2, col="green")
  abline(v=quantile(probs=c(.025, .975), nuggetVarVals), col="purple", lty=2)
  dev.off()
  # for(i in 1:length(covNames)) {
  #   XMarginal = XMarginals[[i]]
  #   pdf(file=paste0("Figures/mixtureLK", covNames[i], ".pdf"), width=5, height=5)
  #   plot(XMarginal, type="l", main="Marginal for fixed effect")
  #   abline(v=0, col="green")
  #   abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
  #   dev.off()
  # }
  
  pdf(file=paste0("Figures/mixtureLKRho", extraPlotName, ".pdf"), width=5, height=5)
  hist(rhoVals, main=TeX("Estimate of $\\rho$"), xlab=TeX("$\\rho$"), breaks=20, freq=FALSE)
  abline(v=quantile(probs=c(.025, .975), rhoVals), col="purple", lty=2)
  abline(v=rho, col="green")
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKa.wght", extraPlotName, ".pdf"), width=5, height=5)
  hist(a.wghtVals, xlab="a.wght", main="Marginal for a.wght", breaks=20, freq=FALSE)
  abline(v=quantile(probs=c(.025, .975), a.wghtVals), col="purple", lty=2)
  dev.off()
  
  # pdf(file=paste0("Figures/mixtureLKRho", extraPlotName, ".pdf"), width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now get distribution for the alpha parameters
  xSamples = alphaVals
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/mixtureLKAlpha", l, extraPlotName, ".pdf"), width=5, height=5)
    hist(xSamples[l,], xlab=TeX(paste0("$\\alpha_", l, "$")), main=TeX(paste0("Marginal for $\\alpha_", l, "$")), breaks=20, freq=F, xlim=c(0,1))
    abline(v=mean(xSamples[l,]), col="purple", lty=1)
    abline(v=quantile(probs=c(.025, .975), xSamples[l,]), col="purple", lty=2)
    dev.off()
  }
  
  ## plot covariance and correlation functions
  # first get the true covariance an correlation functions
  spatialCovFun = spatialCorFun
  mixtureCovFun = function(x) {
    out = spatialCorFun(x)[1,]
    out[x == 0] = 1 + sigma2
    out
  }
  mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sigma2)) }
  
  # compute the covariance function for many different hyperparameter samples
  out = covarianceDistributionLK(mod$LKinfo, alphaVals, lambdaVals, a.wghtVals, rhoVals)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[sortI]
  lowerCov=out$lowerCov[sortI]
  covMat=out$covMat[sortI,]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[sortI]
  lowerCor=out$lowerCor[sortI]
  corMat=out$corMat[sortI,]
  corMatNoNugget=out$corMatNoNugget[sortI,]
  covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
                 corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  
  getEffectiveRange = function(ds, cors) {
    firstI = match(TRUE, cors<.1)
    if(is.na(firstI)) {
      warning("effective range larger than spatial domain diameter. Returning spatial domain diameter as effective range")
      firstI = length(ds)
    } else if(firstI == 1)
      warning("effective range is extremely small and more resolution required to accurately determine it")
    ds[firstI]
  }
  effectiveRanges = apply(corMatNoNugget, 2, getEffectiveRange, ds=d)
  
  pdf(file=paste0("Figures/mixtureLKEffRange", extraPlotName, ".pdf"), width=5, height=5)
  hist(effectiveRanges, main="Estimated effective range", freq=FALSE, breaks=20)
  abline(v=effectiveRange, col="green")
  abline(v=quantile(probs=c(.025, .975), effectiveRanges), col="purple", lty=2)
  dev.off()
  
  # plot the covariance function
  yRange = range(c(covMean, lowerCov, upperCov, mixtureCovFun(d)))
  pdf(file=paste0("Figures/mixtureLKCov", extraPlotName, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Estimated covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  lines(d, mixtureCovFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKCor", extraPlotName, ".pdf"), width=5, height=5)
  plot(d, corMean, type="l", main="Estimated correlation function", xlab="Distance", ylab="Covariance", 
       ylim=c(0, 1))
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  lines(d, mixtureCorFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  # get scoring rules
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  leftOutIndices = (length(preds) - length(ysTest) + 1):(length(preds) - length(ysTest) + length(simulationData$zTest[,1]))
  gridIndices = (length(preds) - length(ysTest) + length(simulationData$zTest[,1]) + 1):length(preds)
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  # first calculate scoring rules at grid points
  truth = ysTest[gridIndicesTest]
  est = preds[gridIndices]
  vars = predSDs[gridIndices]^2
  lower = mod$lower[gridIndices]
  upper = mod$upper[gridIndices]
  
  # compute nearest neighbor distances and scores as a function of them
  gridPts = predPts[gridIndices,]
  distMat = rdist(coords, gridPts)
  nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  print("Binned grid scores:")
  gridScoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  print(gridScoringRules$pooledResults)
  print(gridScoringRules$binnedResults)
  
  # now calculate square rules at left out points
  truth = ysTest[leftOutIndicesTest]
  est = preds[leftOutIndices]
  vars = predSDs[leftOutIndices]^2
  lower = mod$lower[leftOutIndices]
  upper = mod$upper[leftOutIndices]
  leftOutScoringRules = getScores(truth, est, vars, lower, upper)
  leftOutScoringRules = data.frame(c(leftOutScoringRules, Time=time[3]))
  print("Binned left out scores:")
  print(leftOutScoringRules)
  
  scoringRules = list(gridScoringRules=gridScoringRules, leftOutScoringRules=leftOutScoringRules)
  # 
  # # get scoring rules
  # truth = ysTest
  # est = preds[testIndices]
  # vars = predSDs[testIndices]^2
  # lower = mod$lower[testIndices]
  # upper = mod$upper[testIndices]
  # 
  # # compute nearest neighbor distances and scores as a function of them
  # testPts = predPts[testIndices,]
  # distMat = rdist(coords, testPts)
  # nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  # print("Binned scores:")
  # scoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  # scoringRules$pooledResults = data.frame(c(scoringRules$pooledResults, Time=time[3]))
  # print(scoringRules$binnedResults)
  
  
  
  # print("Grid scores:")
  # print(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                 distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults)
  # 
  # # plot scores as a function of distance
  # distanceScores = getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                            distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults
  
  # pdf(file=paste0("Figures/mixtureSPDEScoreBias", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Bias, pch=19, col="blue", main="Bias", ylab="Bias", xlab="Nearest neighbor distance (km)")
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreVar", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Var, pch=19, col="blue", main="Variance", ylab="Variance", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Var)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$MSE, pch=19, col="blue", main="MSE", ylab="MSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$MSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreRMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$RMSE, pch=19, col="blue", main="RMSE", ylab="RMSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$RMSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCRPS", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$CRPS, pch=19, col="blue", main="CRPS", ylab="CRPS", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$CRPS)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCvg", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Coverage, pch=19, col="blue", main="80% Coverage", ylab="80% Coverage", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Coverage)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreWidth", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Width, pch=19, col="blue", main="Width", ylab="Width", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Width)))
  # abline(h=0, lty=2)
  # dev.off()
  
  # get aggregated predictions
  # A = t(sapply(1:(mx*my), getARow))
  # A = sweep(A, 1, rowSums(A), "/")
  # mx = 100
  # my = 100
  # predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  # testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  # gridIndices = (length(preds) - length(simulationData$xGrid) + 1):length(preds)
  # gridIndicesEst = (length(est) - length(simulationData$xGrid) + 1):length(est)
  
  A = makeNumericalIntegralMat(gridCoords, mx=3, my=3)
  aggregatedPreds = A %*% preds[gridIndices]
  
  truth = A %*% ysTest[gridIndicesTest]
  est = A %*% preds[gridIndices]
  aggregatedPredMat = A %*% fit$predMat[gridIndices,]
  vars = apply(aggregatedPredMat, 1, var)
  sds = apply(aggregatedPredMat, 1, sd)
  lower = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.1)})
  upper = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.9)})
  predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  print("Aggregated prediction summary statistics:")
  print(predictionMatrix)
  
  print("Pooled aggregated scores:")
  pooledAggregatedScores = getScores(truth, est, vars, lower, upper)
  print(pooledAggregatedScores)
  print("Left out region aggregated scores:")
  leftOutAggregatedScores = getScores(truth[5], est[5], vars[5], lower[5], upper[5])
  leftOutAggregatedScores$Var = sds[5]
  names(leftOutAggregatedScores)[2] = "Predictive.SD"
  print(leftOutAggregatedScores)
  print("Included regions aggregated scores:")
  includedAggregatedScores = getScores(truth[-5], est[-5], vars[-5], lower[-5], upper[-5])
  print(includedAggregatedScores)
  
  aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
                                includedAggregatedScores=includedAggregatedScores)
  
  fit$mod = NULL
  if(!gscratch)
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("savedOutput/simulations/mixtureLK", plotNameRoot, ".RData"))
  else
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("/work/johnpai/mixtureLK", plotNameRoot, ".RData"))
  
  invisible(list(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules))
}

# seed=1, nLayer=3, nx=20, ny=nx, nu=1, assumeMeanZero=TRUE, 
# nBuffer=5, normalize=TRUE, NC=14, testCovs=TRUE, 
# printVerboseTimings=FALSE, n=900, separatea.wght=FALSE, 
# plotNameRoot="", doMatern=FALSE, fixNu=FALSE, thetas=c(.08, .8) / sqrt(8), 
# testfrac=.1, leaveOutRegion=TRUE, sigma2 = 0.1^2, extraPlotName=plotNameRoot
testLKModelMixtureMultiple = function(seed=1, nSamples=100, gscratch=TRUE, loadResults=FALSE, startI=1, endI=nSamples, ...) {
  # set random seeds for each simulation
  set.seed(seed)
  allSeeds = sample(1:1000000, nSamples, replace = FALSE)
  
  # call testLKINLAModelMixture for each simulation requested
  temp = function(i) {
    print(paste0("Beginning simulation ", i, "/", nSamples))
    thisPlotNameRoot = paste0("sim", i)
    do.call("testLKModelMixture", c(list(seed = allSeeds[i], plotNameRoot=thisPlotNameRoot, gscratch=gscratch), list(...)))
  }
  if(!loadResults)
    sapply(startI:endI, temp)
  
  if(startI != 1 || endI != nSamples)
    return(invisible(NULL))
  
  # load in the results
  allScoringRulesGrid = list()
  allScoringRulesLeftOut = list()
  allFits = list()
  allCovInfo = list()
  allPredictionMatrices = list()
  allAggregatedScoringRules = list()
  for(i in 1:nSamples) {
    if(!gscratch)
      out = load(paste0("savedOutput/simulations/mixtureLKsim", i, ".RData"))
    else
      out = load(paste0("/work/johnpai/mixtureLKsim", i, ".RData"))
    allScoringRulesGrid = c(allScoringRulesGrid, list(scoringRules$gridScoringRules))
    allScoringRulesLeftOut = c(allScoringRulesLeftOut, list(scoringRules$leftOutScoringRules))
    allFits = c(allFits, list(fit))
    allCovInfo = c(allCovInfo, list(covInfo))
    allPredictionMatrices = c(allPredictionMatrices, list(predictionMatrix))
    allAggregatedScoringRules = c(allAggregatedScoringRules, list(aggregatedScoringRules))
  }
  
  ##### average results from each simulation
  # pointwise scoring rules
  allPooledScoringRulesGrid = do.call("rbind", lapply(allScoringRulesGrid, function(x) {x$pooledResults}))
  allBinnedScoringRulesGrid = lapply(allScoringRulesGrid, function(x) {x$binnedResults})
  binnedScoringRulesGrid = averageBinnedScores(allBinnedScoringRulesGrid)
  ns = binnedScoringRulesGrid[,2]
  pooledScoringRulesGrid = apply(binnedScoringRulesGrid, 2, function(x) {sum(x * (ns / sum(ns)))})
  pooledScoringRulesGrid = as.data.frame(matrix(pooledScoringRulesGrid, nrow=1))
  names(pooledScoringRulesGrid) = names(binnedScoringRulesGrid)
  
  fullPooledScoringRulesLeftOut = do.call("rbind", allScoringRulesLeftOut)
  pooledScoringRulesLeftOut = colMeans(fullPooledScoringRulesLeftOut)
  
  # covInfo
  # covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
  #                corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  d = allCovInfo[[1]]$d
  covMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$covMean})))
  upperCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCov})))
  lowerCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCov})))
  corMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$corMean})))
  upperCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCor})))
  lowerCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCor})))
  
  # aggregated scoring rules
  # aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
  #                               includedAggregatedScores=includedAggregatedScores)
  # predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  fullPredictionMatrix = do.call("rbind", allPredictionMatrices)
  leftOutIndices = seq(5, nrow(fullPredictionMatrix), by=9)
  leftOutPredictionMatrix = fullPredictionMatrix[leftOutIndices, ]
  leftInPredictionMatrix = fullPredictionMatrix[-leftOutIndices, ]
  leftOutScores = getScores(leftOutPredictionMatrix[,1], leftOutPredictionMatrix[,2], leftOutPredictionMatrix[,3]^2)
  leftInScores = getScores(leftInPredictionMatrix[,1], leftInPredictionMatrix[,2], leftInPredictionMatrix[,3]^2)
  aggregatedScores = getScores(fullPredictionMatrix[,1], fullPredictionMatrix[,2], fullPredictionMatrix[,3]^2)
  
  ##### Save results
  if(!gscratch) {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("savedOutput/simulations/mixtureLKAll_nsim", nSamples, ".RData"))
  } else {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("/work/johnpai/mixtureLKAll_nsim", nSamples, ".RData"))
  }
}

# recalculate covariances for the LatticeKrig model
fixLKModelMixtureCovariance = function(seed=1, nSamples=100, gscratch=TRUE, startI=1, endI=nSamples, ...) {
  # set random seeds for each simulation
  set.seed(seed)
  allSeeds = sample(1:1000000, nSamples, replace = FALSE)
  
  # load in the results
  for(i in startI:endI) {
    set.seed(allSeeds[i])
    print(paste0("Recalculating covariance for simulation ", i, "/", nSamples))
    if(!gscratch)
      out = load(paste0("savedOutput/simulations/mixtureLKsim", i, ".RData"))
    else
      out = load(paste0("/work/johnpai/mixtureLKsim", i, ".RData"))
    
    # calculate covariance
    out = covarianceDistributionLK(fit$LKinfo, fit$alphaVals, fit$lambdaVals, fit$a.wghtVals, fit$rhoVals)
    d = out$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    covMean = out$cov[sortI]
    upperCov=out$upperCov[sortI]
    lowerCov=out$lowerCov[sortI]
    covMat=out$covMat[sortI,]
    corMean = out$cor[sortI]
    upperCor=out$upperCor[sortI]
    lowerCor=out$lowerCor[sortI]
    corMat=out$corMat[sortI,]
    corMatNoNugget=out$corMatNoNugget[sortI,]
    covInfoNew = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
                   corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
    
    # save results
    covInfoOld = covInfo
    covInfo = covInfoNew
    if(!gscratch)
      save(scoringRules, fit, covInfo, covInfoOld, predictionMatrix, aggregatedScoringRules, file=paste0("savedOutput/simulations/mixtureLK", i, ".RData"))
    else
      save(scoringRules, fit, covInfo, covInfoOld, predictionMatrix, aggregatedScoringRules, file=paste0("/work/johnpai/mixtureLK", i, ".RData"))
  }
}

# tests the fitSPDE function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
# thetas: originally was c(.1, 3), but 3 is too large
testSPDEModelMixture = function(seed=1, nx=20, ny=nx, assumeMeanZero=TRUE, 
                                testCovs=FALSE, n=900, thetas=NULL, 
                                int.strategy="auto", strategy="gaussian", 
                                nPostSamples=1000, mesh=NULL, 
                                prior=NULL, testfrac=1/9, nu=1, nHyperSamples=1000, 
                                plotNameRoot="", sigma2 = 0.1^2, useKenya=FALSE, 
                                urbanOverSamplefrac=0, leaveOutRegion=TRUE, gscratch=FALSE) {
  set.seed(seed)
  
  if(useKenya)
    distanceBreaks = seq(0, 300, l=20)
  else
    distanceBreaks = seq(0, 0.5, l=20)
  
  # set plotNameRoot
  plotNameRoot = paste0(plotNameRoot, "_n", n, "_nu", nu, "_nugV", round(sigma2, 2), "_Kenya", useKenya, 
                        "_noInt", assumeMeanZero, "_urbOversamp", round(urbanOverSamplefrac, 4))
  
  # set true parameter values
  if(useKenya) {
    if(is.null(thetas))
      thetas=c(.08, .8) * (1000/2) / sqrt(8)
  } else {
    if(is.null(thetas))
      thetas=c(.08, .8) / sqrt(8)
  }
  rho = 1
  effectiveRange = thetas * sqrt(8)
  
  # set the SPDE mesh if necessary
  if(is.null(mesh)) {
    if(useKenya)
      mesh = getSPDEMeshKenya()
    else
      mesh = getSPDEMesh()
  }
  
  # load data set if necessary
  if(is.null(n)) {
    out = load("mixtureDataSet.RData")
  } else {
    mixtureCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
        0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
    nTest = round(testfrac * n)
    if(leaveOutRegion) {
      simulationData = getSimulationDataSetsGivenCovarianceTest(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                                nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                                saveDataSetPlot=FALSE, doPredGrid=TRUE)
    } else {
      simulationData = getSimulationDataSetsGivenCovariance(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                            nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                            saveDataSetPlot=FALSE, useKenyaLocations=useKenya, urbanOverSamplefrac=urbanOverSamplefrac)
    }
    
  }
  coords = cbind(simulationData$xTrain[,1], simulationData$yTrain[,1])
  ys = simulationData$zTrain[,1]
  
  # generate lattice and simulate observations
  # coords = matrix(runif(2*n), ncol=2)
  if(useKenya) {
    xRangeDat = simulationData$xRange
    yRangeDat = simulationData$yRange
  } else {
    xRangeDat = c(-1, 1)
    yRangeDat = c(-1, 1)
  }
  
  AObs = inla.spde.make.A(mesh, loc = coords)
  # Q = makeQ(kappa=kappa, rho=rho, latInfo, alphas=alphas, normalized=normalize, fastNormalize=fastNormalize) 
  # L = as.matrix(t(chol(solve(Q))))
  # zsims = matrix(rnorm(nrow(Q)), ncol=1)
  # fieldSims = L %*% zsims
  # ys = as.numeric(AObs %*% fieldSims) + 1 # add a constant unit mean term to be estimated by INLA
  # # ys = 1 + as.numeric(AObs %*% fieldSims) + coords[,1] # x-valued mean term to be estimated by INLA
  # errs = rnorm(n, sd=sqrt(sigma2))
  # ys = ys + errs
  
  # plot the observations
  pdf(file=paste0("Figures/mixtureSPDEObservations", plotNameRoot, ".pdf"), width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid, and add testing points
  mx = 100
  my = 100
  predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  if(useKenya) {
    # remove grid points outside of Kenya national boundaries
    load("../U5MR/adminMapData.RData")
    polys = adm0@polygons
    kenyaPoly = polys[[1]]@Polygons[[77]]@coords
    kenyaPolyProj = projKenya(kenyaPoly)
    inKenya = in.poly(predPts, kenyaPolyProj)
    predPts = predPts[inKenya,]
    
    # add other testing locations to matrix of prediction locations and remember which 
    predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
    predPts = rbind(predPts, cbind(simulationData$xTestRural[,1], simulationData$yTestRural[,1]))
    predPts = rbind(predPts, cbind(simulationData$xTestUrban[,1], simulationData$yTestUrban[,1]))
    plotGridI = 1:sum(inKenya)
    overallTestI = simulationData$overallTestI
    ruralTestI = simulationData$ruralTestI
    urbanTestI = simulationData$urbanTestI
    gridTestI = (max(urbanTestI) + 1):(max(urbanTestI) + length(simulationData$xGrid))
  } else {
    predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
  }
  predPts = rbind(predPts, cbind(simulationData$xGrid, simulationData$yGrid))
  ysTest = c(simulationData$zTest[,1], simulationData$zTestRural[,1], simulationData$zTestUrban[,1], simulationData$zGrid[,1])
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  # priorPar = getPrior(.1, .1, 10)
  # generate hyperparameters for pc priors
  # median effective range is .4 or 200 for Kenya data (a fifth of the spatial domain diameter), median spatial variance is 1
  # priorPar = getPCPrior(.4, .5, 1) 
  if(is.null(prior)) {
    if(!useKenya) {
      prior = getSPDEPrior(mesh, U=1, alpha=.5, medianRange=.4)
    } else {
      prior = getSPDEPrior(mesh, U=1, alpha=.5, medianRange=.4 * 1000 / 2)
    }
  }
  
  X = matrix(rep(1, nrow(coords)), ncol=1)
  # X = matrix(coords[,1], ncol=1)
  XPred = matrix(rep(1, nrow(predPts)), ncol=1)
  
  # add linear terms in lat/lon to covariate matrices if requested
  if(testCovs) {
    X = cbind(X, coords)
    XPred = cbind(XPred, predPts)
  }
  
  if(assumeMeanZero) {
    X = NULL
    XPred = NULL
  }
  
  # show priors on effective correlation, marginal variance, and error variance:
  if(!useKenya)
    xs1 = seq(.01, 1, l=500)
  else
    xs1 = seq(1, 1000, l=500)
  # pdf(file=paste0("Figures/mixtureSPDEPriorEffRange", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
  #      xlab="Effective Correlation Range", main="Effective Correlation Prior", 
  #      ylab="Prior Density")
  # abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  # dev.off()
  
  # if(priorPar$priorType == "orig") {
  #   xs2 = seq(.01, 10.5, l=500)
  #   pdf(file=paste0("Figures/mixtureSPDEPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
  #   plot(xs2, invgamma::dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
  #        xlab="Marginal Variance", main="Marginal Variance Prior", 
  #        ylab="Prior Density")
  #   abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  #   abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  #   dev.off()
  # } else if(priorPar$priorType == "pc") {
  #   xs2 = seq(.01, 11.5, l=500)
  #   pdf(file=paste0("Figures/mixtureSPDEPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
  #   plot(xs2, dpcvar(xs2, alpha=priorPar$alpha, u=priorPar$u), type="l", col="blue", 
  #        xlab="Marginal Variance", main="Marginal Variance Prior", 
  #        ylab="Prior Density")
  #   abline(v=qpcvar(.1, alpha=priorPar$alpha, u=priorPar$u), col="red")
  #   abline(v=qpcvar(.9, alpha=priorPar$alpha, u=priorPar$u), col="red")
  #   abline(v=1, col="green")
  #   dev.off()
  # }
  
  # xs2 = seq(.001, invgamma::qinvgamma(.905, shape=0.1, rate=0.1), l=500)
  # pdf(file="Figures/mixtureSPDEPriorErrorVar.pdf", width=5, height=5)
  # plot(xs2, invgamma::dinvgamma(xs2, shape=0.1, rate=0.1), type="l", col="blue", 
  #      xlab="Error Variance", main="Error Variance Prior", 
  #      ylab="Prior Density")
  # abline(v=invgamma::qinvgamma(.1, shape=0.1, rate=0.1), col="red")
  # abline(v=invgamma::qinvgamma(.9, shape=0.1, rate=0.1), col="red")
  # dev.off()
  
  xs2 = seq(.01, 1, l=500)
  pdf(file=paste0("Figures/mixtureSPDEPriorErrorVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(xs2, dpcvar(xs2, alpha=.05, u=1), type="l", col="blue", 
       xlab="Marginal Variance", main="Marginal Variance Prior", 
       ylab="Prior Density")
  abline(v=qpcvar(.1, alpha=.05, u=1), col="red")
  abline(v=qpcvar(.9, alpha=.05, u=1), col="red")
  abline(v=sqrt(.1), col="green")
  dev.off()
  
  # browser()
  
  # fit the model
  time = system.time(fit <- fitSPDE(coords, ys, predCoords=predPts, seed=seed, prior=prior, 
                                    xObs=X, xPred=XPred, int.strategy=int.strategy, strategy=strategy, 
                                    mesh=mesh, nPostSamples=nPostSamples))
  mod = fit$mod
  preds=fit$preds
  predSDs=fit$sigmas
  latInfo=fit$latInfo
  latWidth=fit$latWidth
  obsPreds=fit$obsPreds
  obsSDs=fit$obsSDs
  
  # print out the total time
  print(paste0("Total time: ", time[3]))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  if(!useKenya) {
    inRange = function(pts, rangeShrink=0) {
      inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
      inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
      inX & inY
    }
  } else {
    inRange = function(pts, rangeShrink=0) {
      rep(TRUE, nrow(pts))
    }
  }
  
  # show predictive surface, SD, and data
  
  pdf(file=paste0("Figures/mixtureSPDEPreds", plotNameRoot, ".pdf"), width=8, height=8)
  par(mfrow=c(2,2), mar=c(5.1, 4.1, 4.1, 6))
  
  # obsInds = 1:n
  # predInds = (n+1):(n+mx*my)
  # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
  # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
  colRangeDat = range(c(ys, obsPreds, preds))
  colRangeSD = range(c(predSDs, obsSDs))
  quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  
  quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
             xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureSPDELeftOutResiduals", plotNameRoot, ".pdf"), width=5, height=5)
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  plot(preds[testIndices], ysTest-preds[testIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted", 
       ylab="Residuals", xlab="Fitted")
  abline(h=0, lty=2)
  dev.off()
  
  if(useKenya) {
    pdf(file=paste0("Figures/mixtureSPDELeftOutResidualsLabeled", plotNameRoot, ".pdf"), width=5, height=5)
    testIndices = (length(preds) - length(ysTest) + 1):length(preds)
    gridIndices = 
      ylim = range(ysTest-preds[testIndices])
    xlim = range(preds[testIndices])
    plot(preds[testIndices][overallTestI], ysTest[overallTestI]-preds[testIndices][overallTestI], pch=19, cex=.1, col="black", main="Residuals versus fitted", 
         ylab="Residuals", xlab="Fitted", xlim=xlim, ylim=ylim)
    points(preds[testIndices][ruralTestI], ysTest[ruralTestI]-preds[testIndices][ruralTestI], pch=19, cex=.1, col="green")
    points(preds[testIndices][urbanTestI], ysTest[urbanTestI]-preds[testIndices][urbanTestI], pch=19, cex=.1, col="blue")
    points(preds[testIndices][gridTestI], ysTest[gridTestI]-preds[testIndices][gridTestI], pch=19, cex=.1, col="red")
    abline(h=0, lty=2)
    legend("topright", c("Overall", "Rural", "Urban", "Grid"), col=c("black", "green", "blue", "red"), pch=19)
    dev.off()
  }
  
  # calculate true effective range and marginal variance:
  marginalVar = rho
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  effRangeMarg = mod$marginals.hyperpar$`Range for field`
  varMarg = inla.tmarginal(function(x) {x^2}, mod$marginals.hyperpar$`Stdev for field`)
  sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar$`Precision for the Gaussian observations`)
  covNames = names(mod$marginals.fixed)
  
  if(!assumeMeanZero) {
    XMarginals = list()
    for(i in 1:length(covNames)) {
      XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
      XMarginals = c(XMarginals, list(XMarginal))
    }
  }
  
  par(mfrow=c(1,1))
  pdf(file=paste0("Figures/mixtureSPDEEffRange", plotNameRoot, ".pdf"), width=5, height=5)
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  abline(v=effectiveRange, col="green")
  abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file=paste0("Figures/mixtureSPDEVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(varMarg, type="l", main="Marginal for spatial variance")
  abline(v=marginalVar, col="green")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file=paste0("Figures/mixtureSPDESigma2", plotNameRoot, ".pdf"), width=5, height=5)
  plot(sigma2Marg, type="l", main="Marginal for error variance")
  abline(v=sigma2, col="green")
  abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
  dev.off()
  
  if(!assumeMeanZero) {
    for(i in 1:length(covNames)) {
      XMarginal = XMarginals[[i]]
      pdf(file=paste0("Figures/mixtureSPDE", covNames[i], plotNameRoot, ".pdf"), width=5, height=5)
      plot(XMarginal, type="l", main="Marginal for fixed effect")
      abline(v=0, col="green")
      abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
      dev.off()
    }
  }
  
  # pdf(file="Figures/mixtureSPDERho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now generate marginals for the alpha parameters. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  out = inla.hyperpar.sample(nHyperSamples, mod, improve.marginals=TRUE)
  
  ## plot covariance and correlation functions
  # first get the true covariance an correlation functions
  spatialCovFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu, distMat=x)}
  mixtureCovFun = function(x) {
    out = spatialCovFun(x)
    out[x == 0] = 1 + sigma2
    out
  }
  mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sigma2)) }
  
  # first to transform all the hyperparameter samples to their relevant values
  # 1: error precision
  # 2: log effective range
  # 3: log spatial variance
  # 4-(3 + nLayer - 1): multivariateLogit alpha
  nuggetVarVals = 1/out[,1]
  effectiveRangeVals = out[,2]
  varVals = out[,3]^2
  
  # compute the covariance function for many different hyperparameter samples
  out = covarianceDistributionSPDE(effectiveRangeVals, varVals, nuggetVarVals, mesh, xRangeDat=xRangeDat, yRangeDat=yRangeDat)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[,sortI]
  lowerCov=out$lowerCov[,sortI]
  covMat=out$covMat[sortI,]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[,sortI]
  lowerCor=out$lowerCor[,sortI]
  corMat=out$corMat[sortI,]
  covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
                 corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  
  # plot the covariance function
  yRange = range(c(covMean, lowerCov, upperCov, mixtureCovFun(d)))
  pdf(file=paste0("Figures/mixtureSPDECov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov[1,], lty=2)
  lines(d, upperCov[1,], lty=2)
  lines(d, lowerCov[2,], lty=3)
  lines(d, upperCov[2,], lty=3)
  lines(d, mixtureCovFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI", "95% CI"), lty=c(1, 1, 2, 3), col=c("green", "black", "black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureSPDECor", plotNameRoot, ".pdf"), width=5, height=5)
  yRange = range(c(corMean, lowerCor, upperCor, mixtureCorFun(d)))
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCor[1,], lty=2)
  lines(d, upperCor[1,], lty=2)
  lines(d, lowerCor[2,], lty=3)
  lines(d, upperCor[2,], lty=3)
  lines(d, mixtureCorFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI", "95% CI"), lty=c(1, 1, 2, 3), col=c("green", "black", "black", "black"))
  dev.off()
  
  # get scoring rules
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  leftOutIndices = (length(preds) - length(ysTest) + 1):(length(preds) - length(ysTest) + length(simulationData$zTest[,1]))
  gridIndices = (length(preds) - length(ysTest) + length(simulationData$zTest[,1]) + 1):length(preds)
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  # first calculate scoring rules at grid points
  truth = ysTest[gridIndicesTest]
  est = preds[gridIndices]
  vars = predSDs[gridIndices]^2
  lower = fit$lower[gridIndices]
  upper = fit$upper[gridIndices]
  
  # compute nearest neighbor distances and scores as a function of them
  gridPts = predPts[gridIndices,]
  distMat = rdist(coords, gridPts)
  nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  print("Binned grid scores:")
  gridScoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  print(gridScoringRules$pooledResults)
  print(gridScoringRules$binnedResults)
  
  # now calculate square rules at left out points
  truth = ysTest[leftOutIndicesTest]
  est = preds[leftOutIndices]
  vars = predSDs[leftOutIndices]^2
  lower = fit$lower[leftOutIndices]
  upper = fit$upper[leftOutIndices]
  leftOutScoringRules = getScores(truth, est, vars, lower, upper)
  leftOutScoringRules = data.frame(c(leftOutScoringRules, Time=time[3]))
  print("Binned left out scores:")
  print(leftOutScoringRules)
  
  scoringRules = list(gridScoringRules=gridScoringRules, leftOutScoringRules=leftOutScoringRules)
  # 
  # # get scoring rules
  # truth = ysTest
  # est = preds[testIndices]
  # vars = predSDs[testIndices]^2
  # lower = fit$lower[testIndices]
  # upper = fit$upper[testIndices]
  # 
  # # compute nearest neighbor distances and scores as a function of them
  # testPts = predPts[testIndices,]
  # distMat = rdist(coords, testPts)
  # nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  # print("Binned scores:")
  # scoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  # scoringRules$pooledResults = data.frame(c(scoringRules$pooledResults, Time=time[3]))
  # print(scoringRules$binnedResults)
  
  
  
  # print("Grid scores:")
  # print(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                 distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults)
  # 
  # # plot scores as a function of distance
  # distanceScores = getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                            distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults
  
  # pdf(file=paste0("Figures/mixtureSPDEScoreBias", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Bias, pch=19, col="blue", main="Bias", ylab="Bias", xlab="Nearest neighbor distance (km)")
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreVar", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Var, pch=19, col="blue", main="Variance", ylab="Variance", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Var)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$MSE, pch=19, col="blue", main="MSE", ylab="MSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$MSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreRMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$RMSE, pch=19, col="blue", main="RMSE", ylab="RMSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$RMSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCRPS", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$CRPS, pch=19, col="blue", main="CRPS", ylab="CRPS", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$CRPS)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCvg", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Coverage, pch=19, col="blue", main="80% Coverage", ylab="80% Coverage", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Coverage)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreWidth", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Width, pch=19, col="blue", main="Width", ylab="Width", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Width)))
  # abline(h=0, lty=2)
  # dev.off()
  
  if(!useKenya) {
    # print("Pooled scores:")
    # print(data.frame(c(getScores(truth, est, vars, lower, upper), Time=time[3])))
  } else {
    print("Pooled scores:")
    print(data.frame(c(getScores(truth, est, vars, lower, upper), Time=time[3])))
    print("Overall scores:")
    print(data.frame(c(getScores(truth[overallTestI], est[overallTestI], vars[overallTestI], lower[overallTestI], upper[overallTestI]), Time=time[3])))
    print("Rural scores:")
    print(data.frame(c(getScores(truth[ruralTestI], est[ruralTestI], vars[ruralTestI], lower[ruralTestI], upper[ruralTestI]), Time=time[3])))
    print("Urban scores:")
    print(data.frame(c(getScores(truth[urbanTestI], est[urbanTestI], vars[urbanTestI], lower[urbanTestI], upper[urbanTestI]), Time=time[3])))
    # print("Grid scores:")
    # print(data.frame(c(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI]), Time=time[3])))
  }
  
  # get aggregated predictions
  # A = t(sapply(1:(mx*my), getARow))
  # A = sweep(A, 1, rowSums(A), "/")
  # mx = 100
  # my = 100
  # predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  # testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  # gridIndices = (length(preds) - length(simulationData$xGrid) + 1):length(preds)
  # gridIndicesEst = (length(est) - length(simulationData$xGrid) + 1):length(est)
  
  A = makeNumericalIntegralMat(gridCoords, mx=3, my=3)
  aggregatedPreds = A %*% preds[gridIndices]
  
  truth = A %*% ysTest[gridIndicesTest]
  est = A %*% preds[gridIndices]
  aggregatedPredMat = A %*% fit$predMat[gridIndices,]
  vars = apply(aggregatedPredMat, 1, var)
  sds = apply(aggregatedPredMat, 1, sd)
  lower = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.1)})
  upper = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.9)})
  predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  print("Aggregated prediction summary statistics:")
  print(predictionMatrix)
  
  print("Pooled aggregated scores:")
  pooledAggregatedScores = getScores(truth, est, vars, lower, upper)
  print(pooledAggregatedScores)
  print("Left out region aggregated scores:")
  leftOutAggregatedScores = getScores(truth[5], est[5], vars[5], lower[5], upper[5])
  leftOutAggregatedScores$Var = sds[5]
  names(leftOutAggregatedScores)[2] = "Predictive.SD"
  print(leftOutAggregatedScores)
  print("Included regions aggregated scores:")
  includedAggregatedScores = getScores(truth[-5], est[-5], vars[-5], lower[-5], upper[-5])
  print(includedAggregatedScores)
  
  aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
                                includedAggregatedScores=includedAggregatedScores)
  
  fit$mod = NULL
  if(!gscratch)
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("savedOutput/simulations/mixtureSPDE", plotNameRoot, ".RData"))
  else
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("/work/johnpai/mixtureSPDE", plotNameRoot, ".RData"))
  
  invisible(list(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules))
}

# runs the testSPDEModelMixture function for multiple realizations, saves results
testSPDEModelMixtureMultiple = function(seed=1, nSamples=100, n=900, nu=1, sigma2=0.1^2, 
                                        useKenya=FALSE, assumeMeanZero=TRUE, urbanOverSamplefrac=0, 
                                        gscratch=TRUE, ...) {
  # set random seeds for each simulation
  set.seed(seed)
  allSeeds = sample(1:1000000, nSamples, replace = FALSE)
  
  # call testSPDEModelMixture for each simulation requested
  temp = function(i) {
    print(paste0("Beginning simulation ", i, "/", nSamples))
    thisPlotNameRoot = paste0("sim", i)
    do.call("testSPDEModelMixture", c(list(seed = allSeeds[i], n=n, nu=nu, sigma2=sigma2, gscratch=gscratch, 
                                             useKenya=useKenya, urbanOverSamplefrac=urbanOverSamplefrac, assumeMeanZero=assumeMeanZero, plotNameRoot=thisPlotNameRoot), list(...)))
  }
  sapply(1:nSamples, temp)
  
  # set plotNameRoot
  plotNameRoot = paste0("_n", n, "_nu", nu, "_nugV", round(sigma2, 2), "_Kenya", useKenya, 
                        "_noInt", assumeMeanZero, "_urbOversamp", round(urbanOverSamplefrac, 4))
  
  # save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("savedOutput/simulations/mixtureSPDE", plotNameRoot, ".RData"))
  allScoringRulesGrid = list()
  allScoringRulesLeftOut = list()
  allFits = list()
  allCovInfo = list()
  allPredictionMatrices = list()
  allAggregatedScoringRules = list()
  for(i in 1:nSamples) {
    if(!gscratch)
      out = load(paste0("savedOutput/simulations/mixtureSPDEsim", i, plotNameRoot, ".RData"))
    else
      out = load(paste0("/work/johnpai/mixtureSPDEsim", i, plotNameRoot, ".RData"))
    allScoringRulesGrid = c(allScoringRulesGrid, list(scoringRules$gridScoringRules))
    allScoringRulesLeftOut = c(allScoringRulesLeftOut, list(scoringRules$leftOutScoringRules))
    allFits = c(allFits, list(fit))
    allCovInfo = c(allCovInfo, list(covInfo))
    allPredictionMatrices = c(allPredictionMatrices, list(predictionMatrix))
    allAggregatedScoringRules = c(allAggregatedScoringRules, list(aggregatedScoringRules))
  }
  
  ##### average results from each simulation
  # pointwise scoring rules
  allPooledScoringRulesGrid = do.call("rbind", lapply(allScoringRulesGrid, function(x) {x$pooledResults}))
  allBinnedScoringRulesGrid = lapply(allScoringRulesGrid, function(x) {x$binnedResults})
  binnedScoringRulesGrid = averageBinnedScores(allBinnedScoringRulesGrid)
  ns = binnedScoringRulesGrid[,2]
  pooledScoringRulesGrid = apply(binnedScoringRulesGrid, 2, function(x) {sum(x * (ns / sum(ns)))})
  pooledScoringRulesGrid = as.data.frame(matrix(pooledScoringRulesGrid, nrow=1))
  names(pooledScoringRulesGrid) = names(binnedScoringRulesGrid)
  
  fullPooledScoringRulesLeftOut = do.call("rbind", allScoringRulesLeftOut)
  pooledScoringRulesLeftOut = colMeans(fullPooledScoringRulesLeftOut)
  
  # covInfo
  # covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
  #                corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  d = allCovInfo[[1]]$d
  covMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$covMean})))
  upperCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCov[1,]})))
  lowerCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCov[1,]})))
  corMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$corMean})))
  upperCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCor[1,]})))
  lowerCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCor[1,]})))
  
  # aggregated scoring rules
  # aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
  #                               includedAggregatedScores=includedAggregatedScores)
  # predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  fullPredictionMatrix = do.call("rbind", allPredictionMatrices)
  leftOutIndices = seq(5, nrow(fullPredictionMatrix), by=9)
  leftOutPredictionMatrix = fullPredictionMatrix[leftOutIndices, ]
  leftInPredictionMatrix = fullPredictionMatrix[-leftOutIndices, ]
  leftOutScores = getScores(leftOutPredictionMatrix[,1], leftOutPredictionMatrix[,2], leftOutPredictionMatrix[,3]^2)
  leftInScores = getScores(leftInPredictionMatrix[,1], leftInPredictionMatrix[,2], leftInPredictionMatrix[,3]^2)
  aggregatedScores = getScores(fullPredictionMatrix[,1], fullPredictionMatrix[,2], fullPredictionMatrix[,3]^2)
  
  ##### Save results
  if(!gscratch) {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("savedOutput/simulations/mixtureSPDEAll_nsim", nSamples, plotNameRoot, ".RData"))
  } else {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("/work/johnpai/mixtureSPDEAll_nsim", nSamples, plotNameRoot, ".RData"))
  }
}

# test how close we can get to the spatial correlation function:
testLKINLACorrelationApproximation = function(seed=1, nLayer=3, NP=200, 
                                        nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=13, 
                                        latInfo=NULL, thetas=c(.1, .4), 
                                        initialEffectiveRange=1, initialAlphas=rep(1/nLayer, nLayer-1), 
                                        nu=.5) {
  set.seed(seed)
  
  # set true parameter values
  rho = 1
  effectiveRange = thetas * sqrt(8)
  sigma2 = sqrt(.1)
  # spatialCorFun = function(x) {0.4 * Exp.cov(x, theta=thetas[1], distMat=x) + 0.6 * Exp.cov(x, theta=thetas[2], distMat=x)}
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", distMat=x, smoothness=nu)}
  
  # construct the lattice
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  if(is.null(latInfo))
    latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  
  # set initial parameter values
  initialKappa = (sqrt(8) * latInfo[[1]]$latWidth /initialEffectiveRange)
  
  xlim <- latInfo[[1]]$xRangeDat
  ux <- seq(xlim[1], xlim[2], , NP)
  ylim <- latInfo[[1]]$yRangeDat
  uy <- seq(ylim[1], ylim[2], , NP)
  center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  # precompute relevant matrices
  Qprecomputations = precomputationsQ2(latInfo)
  Acenter = makeA(center, latInfo)
  Ax = makeA(cbind(ux, uy[NP/2]), latInfo)
  Ay = makeA(cbind(ux[NP/2], uy), latInfo)
  
  testFun = function(parameters) {
    kappa = exp(parameters[1])
    logitAlphas = parameters[2:nLayer]
    alphas = multivariateExpit(logitAlphas)
    alphas = c(alphas, 1 - sum(alphas))
    test = getLKInlaCovarianceFun(kappa, 1, 0, alphas, latticeInfo = latInfo, 
                                  precomputedMatrices=Qprecomputations, precomputedAcenter=Acenter, precomputedAx=Ax, precomputedAy=Ay)
    
    ds = test[-1,1]
    mean((spatialCorFun(ds) - test[-1,2])^2 * (1/ds)) # inverse distance weighted mean squared error of the correlation function fit
    # mean((spatialCorFun(ds) - test[-1,2])^2)
  }
  
  print("Beginning optimization...")
  temp = optim(c(log(initialKappa), multivariateLogit(initialAlphas)), testFun)
  outPar = temp$par
  outkappa = exp(outPar[1])
  outlogitAlphas = outPar[2:nLayer]
  outalphas = multivariateExpit(outlogitAlphas)
  outalphas = c(outalphas, 1 - sum(outalphas))
  outEffectiveRange = (1/outkappa) * sqrt(8) * latInfo[[1]]$latWidth
  
  trueEffectiveRangeOverall = getTrueLKEffectiveRange(nLayer, NP, sigma2=0, rho=1, nBuffer, normalize, fastNormalize, 
                                                      NC, latInfo, outEffectiveRange, outalphas)
  trueEffectiveRanges = rep(0, nLayer)
  for(i in 1:nLayer) {
    theseAlphas = rep(0, nLayer)
    theseAlphas[i] = 1
    trueEffectiveRanges[i] = getTrueLKEffectiveRange(nLayer, NP, sigma2=0, rho=1, nBuffer, normalize, fastNormalize, 
                                                     NC, latInfo, outEffectiveRange, theseAlphas)
    
    print(paste0("Layer ", i, " weight: ", outalphas[i], ", true effective range: ", trueEffectiveRanges[i]))
  }
  
  test = getLKInlaCovarianceFun(outkappa, 1, 0, outalphas, latticeInfo = latInfo)
  ds = test[,1]
  pdf(paste0("Figures/approxCorrelation_nLayer", nLayer, "_NC", NC, "_nBuffer", nBuffer, "_nu", nu, ".pdf"), width=5, height=5)
  plot(ds, spatialCorFun(ds), type="l", col="green", ylim=c(0,1), ylab="Correlation", main="Correlation", xlab="Distance")
  lines(ds, test[,2], col="blue")
  legend("topright", c("Truth", "Approximate"), col=c("green", "blue"), lty=1)
  dev.off()
}

# get correlation function from the Matern family with smoothness given by nu that is best approximation to
# the exponential mixture correlation function
getCorrelationApproximation = function(thetas=c(.1, .4), nu=0.5) {
  NP = 200
  # set true parameter values
  rho = 1
  effectiveRange = thetas * sqrt(8)
  sigma2 = sqrt(.1)
  # spatialCorFun = function(x) {0.4 * Exp.cov(x, theta=thetas[1], distMat=x) + 0.6 * Exp.cov(x, theta=thetas[2], distMat=x)}
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", distMat=x, smoothness=nu)}
  
  # construct the lattice
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  ds = seq(0, 1, l=NP)[-1]
  
  testFun = function(parameters) {
    theta = exp(parameters[1])
    
    mean((spatialCorFun(ds) - Matern(ds, range=theta, smoothness=1))^2 / ds) # inverse distance weighted mean squared error of the correlation function fit
    # mean((spatialCorFun(ds) - test[-1,2])^2)
  }
  
  print("Beginning optimization...")
  temp = optim(log(.2), testFun)
  theta = exp(temp$par)
  
  pdf(paste0("Figures/approxCorrelation_MaternNu1.pdf"), width=5, height=5)
  plot(ds, spatialCorFun(ds), type="l", col="green", ylim=c(0,1), ylab="Correlation", main="Correlation", xlab="Distance")
  lines(ds, Matern(ds, range=theta, smoothness=1), col="blue")
  legend("topright", c("Truth", "Approximate"), col=c("green", "blue"), lty=1)
  dev.off()
  theta # 0.1256181
}

maternMixture = function(x, thetas=c(.1, .4), nu=0.5, alphas=c(0.5, 0.5)) {
  alphas[1] * Matern(x, range=thetas[1], smoothness=nu) + 
    alphas[2] * Matern(x, range=thetas[2], smoothness=nu)
}
maternMixtureCor = function(x1, x2=NULL, thetas=c(.1, .4), nu=0.5, alphas=c(0.5, 0.5), distMat=NA) {
  tempFun = function(d) {maternMixture(d, thetas, nu, alphas)}
  stationary.cov(x1, x2, Covariance=tempFun, distMat=distMat)
}
maternApproxCor = function(x1, x2=NULL, theta=0.1256181, smoothness=1, distMat=NA) {
  approxCorrelation=function(d) {Matern(d, range=theta, smoothness=smoothness)}
  stationary.cov(x1, x2, Covariance=approxCorrelation, distMat=distMat)
}

# nx, ny: resolution of prediction grid in the x and y directions
# mx, my: the number of aggregation regions (rectangles) in the x and y directions
testCorrelationApproximation = function(trueCorrelation = maternMixtureCor, 
                                        approxCorrelation = maternApproxCor, 
                                        nx=75, ny=75, seed=1, n=900, 
                                        mx=3, my=3) {
  set.seed(seed)
  
  # set true parameter values
  rho = 1
  # sigma2 = sqrt(.1)
  sigma2 = 0
  
  # generate data set
  nTest = 100
  simulationData = getSimulationDataSetsGivenCovarianceTest(trueCorrelation, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=0, 
                                                        nDataSets=2, plotNameRoot=paste0("(0.5*Exp(.1) + 0.5*Exp(.4))"), fileNameRoot="mix", 
                                                        saveDataSetPlot=FALSE, doPredGrid=TRUE)
  coords = cbind(simulationData$xTrain[,1], simulationData$yTrain[,1])
  testCoords = cbind(simulationData$xTest[,1], simulationData$yTest[,1])
  ys = simulationData$zTrain[,1]
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  gridValues = simulationData$zGrid[,1]
  
  # generate lattice and simulate observations
  # coords = matrix(runif(2*n), ncol=2)
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  
  # plot the observations
  pdf(file="Figures/mixtureMaternObservations.pdf", width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid, and add testing points
  xRange=c(-1,1)
  yRange=c(-1,1)
  predPts = rbind(gridCoords, testCoords)
  ysGrid = gridValues
  ysTest = simulationData$zTest[,1]
  
  # construct aggregation matrix for predictions by testing which prediction locations 
  # are in which aggregation regions
  xRegionGrid = seq(-1, 1, l=mx + 1)[-1]
  yRegionGrid = seq(-1, 1, l=my + 1)[-1]
  xRegion = function(x) {
    match(TRUE, x <= xRegionGrid)
  }
  yRegion = function(y) {
    match(TRUE, y <= yRegionGrid)
  }
  xRegionI = sapply(gridCoords[,1], xRegion)
  yRegionI = sapply(gridCoords[,2], yRegion)
  regionI = (yRegionI-1)*mx + xRegionI
  getARow = function(ai) {regionI == ai}
  
  # A = t(sapply(1:(mx*my), getARow))
  # A = sweep(A, 1, rowSums(A), "/")
  A = makeNumericalIntegralMat(gridCoords, mx=mx, my=my)
  
  # generate the true aggregated values of the field
  # NOTE: this includes any nugget effect. Is that a good idea??
  gridValuesAggregated = A %*% gridValues
  
  # fit the model true and approximate models
  time = system.time(fit <- GPpreds(coords, ys, predCoords=predPts, 
                                    xObs=NULL, xPred=NULL, cov.fun=trueCorrelation, A=A))
  predsGrid=fit$preds[1:ncol(A)]
  predSDsGrid=fit$sigmas[1:ncol(A)]
  predsTest=fit$preds[-(1:ncol(A))]
  predSDsTest=fit$sigmas[-(1:ncol(A))]
  predsAggregated=fit$predsAggregated
  predSDsAggregated=fit$sigmasAggregated
  
  time = system.time(fit <- GPpreds(coords, ys, predCoords=predPts, 
                                    xObs=NULL, xPred=NULL, cov.fun=approxCorrelation, A=A))
  predsApproxGrid=fit$preds[1:ncol(A)]
  predSDsApproxGrid=fit$sigmas[1:ncol(A)]
  predsApproxTest=fit$preds[-(1:ncol(A))]
  predSDsApproxTest=fit$sigmas[-(1:ncol(A))]
  predsApproxAggregated=fit$predsAggregated
  predSDsApproxAggregated=fit$sigmasAggregated
  
  # print results
  print(data.frame(list(truth=gridValuesAggregated, truePredsAggregated=predsAggregated, approxPredsAggregated=predsApproxAggregated, trueSDsAggregated=predSDsAggregated, approxSDsAggregated=predSDsApproxAggregated)))
  print(data.frame(list(trueMean=mean(ysTest), truePredsMean=mean(predsTest), approxPredsMean=mean(predsApproxTest), 
                        trueSDsMean=mean(predSDsTest), approxSDsMean=mean(predSDsApproxTest), 
                        trueMSE=mean((ysTest-predsTest)^2), approxMSE=mean((ysTest-predsApproxTest)^2))))
  
  # Calculate all scoring rules and print them
  aggregatedScoresTrue = getScores(gridValuesAggregated, predsAggregated, predSDsAggregated)
  aggregatedScoresApprox = getScores(gridValuesAggregated, predsApproxAggregated, predSDsApproxAggregated)
  scoresTrue = getScores(ysTest, predsTest, predSDsTest)
  scoresApprox = getScores(ysTest, predsApproxTest, predSDsApproxTest)
  
  print(aggregatedScoresTrue)
  print(aggregatedScoresApprox)
  print(scoresTrue)
  print(scoresApprox)
  
  # function for determining if points are in correct range
  inRange = function(pts, rangeShrink=0) {
    inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
    inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
    inX & inY
  }
  
  # show predictive surface, SD, and data
  browser()
  
  pdf(file="Figures/mixtureMaternPreds.pdf", width=12, height=8)
  par(mfrow=c(2,3), mar=c(5.1, 4.1, 4.1, 6))
  
  colRangeDat = range(c(ys, predsGrid, predsApproxGrid, gridValues))
  colRangeSD = range(c(predSDsGrid[inRange(gridCoords, rangeShrink=.01)], predSDsApproxGrid[inRange(gridCoords, rangeShrink=.01)]))
  quilt.plot(gridCoords[,1], gridCoords[,2], gridValues, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridCoords[,1], gridCoords[,2], predsGrid, main="Prediction Mean (True Covariance)", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridCoords[,1], gridCoords[,2], predsApproxGrid, main="Prediction Mean (Approx. Covariance)", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  
  quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridCoords[,1], gridCoords[,2], predSDsGrid, main="Prediction SD (True Covariance)",
             xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeSD)
  quilt.plot(gridCoords[,1], gridCoords[,2], predSDsApproxGrid, main="Prediction SD (Approx. Covariance)",
             xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeSD)
  dev.off()
  
  my_line <- function(x,y,...){
    abline(a = 0,b = 1,...)
    points(x,y,..., col="blue")
  }
  
  pdf(file="Figures/mixtureMaternLeftOutPairGrid.pdf", width=6, height=6)
  thisRange = range(c(cbind(Truth=gridValues, BLUP=predsGrid, ApproxCov=predsApproxGrid)))
  pairs(cbind(Truth=gridValues, BLUP=predsGrid, ApproxCov=predsApproxGrid), pch=19, cex=.1, xlim=thisRange, ylim=thisRange, 
        upper.panel=my_line, lower.panel=my_line, main="Full Domain Grid Pair Plot")
  dev.off()
  
  pdf(file="Figures/mixtureMaternLeftOutPairTest.pdf", width=6, height=6)
  thisRange = range(c(cbind(Truth=ysTest, BLUP=predsTest, ApproxCov=predsApproxTest)))
  pairs(cbind(Truth=ysTest, BLUP=predsTest, ApproxCov=predsApproxTest), pch=19, cex=.4, xlim=thisRange, ylim=thisRange, 
        upper.panel=my_line, lower.panel=my_line, main="Left Out Area Pair Plot")
  dev.off()
  
  # pdf(file="Figures/mixtureMaternLeftOutResidualsGrid.pdf", width=5, height=5)
  # testIndices = (length(predsGrid) - length(ysGrid) + 1):length(predsGrid)
  # plot(predsGrid[testIndices], ysGrid-predsGrid[testIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted", 
  #      ylab="Residuals", xlab="Fitted")
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file="Figures/mixtureMaternLeftOutResidualsTest.pdf", width=5, height=5)
  # testIndices = (length(predsGrid) - length(ysGrid) + 1):length(predsGrid)
  # plot(predsGrid[testIndices], ysGrid-predsGrid[testIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted", 
  #      ylab="Residuals", xlab="Fitted")
  # abline(h=0, lty=2)
  # dev.off()
}

# The same as getSimulationDataSetsGivenCovariance, but test observations all come from a left out region in the center
getSimulationDataSetsGivenCovarianceTest = function(corFun, nTotal=900, nTest=round(nTotal / 9), marginalVar=1, errorVar=sqrt(.1), nDataSets=100, 
                                                    printEvery=10, saveDataSetPlot=TRUE, fileNameRoot="", plotNameRoot="", 
                                                    doPredGrid=FALSE, nTestGrid=70^2) {
  
  # set the spatial domain
  xRange = c(-1, 1)
  yRange = c(-1, 1)
  
  # set the region grid
  mx=3
  my=3
  xRegionGrid = seq(-1, 1, l=mx + 1)[-1]
  yRegionGrid = seq(-1, 1, l=my + 1)[-1]
  centerxLeftI = floor((mx - 1) / 2)
  
  # if generating prediction grid, make that grid here
  if(doPredGrid) {
    nx = round(sqrt(nTestGrid))
    ny = nx
    if(nx * ny != nTestGrid)
      stop("nTestGrid must be a square number if we are constructing a prediction grid")
    xValuesGrid = seq(-1, 1, l=nx)
    yValuesGrid = xValuesGrid
    glist = make.surface.grid(list(x=xValuesGrid, y=yValuesGrid))
    xValuesGrid = glist[,1]
    yValuesGrid = glist[,2]
  } else {
    nx = 0
    ny = 0
    xValuesGrid = c()
    yValuesGrid = c()
  }
  
  # simulate observation spatial locations
  # xValues = matrix(runif(nTotal * nDataSets, xRange[1], xRange[2]), ncol=nDataSets)
  # yValues = matrix(runif(nTotal * nDataSets, yRange[1], yRange[2]), ncol=nDataSets)
  nTrain = nTotal - nTest
  coordsTrain = runifsqMissingRectangle(nTrain * nDataSets)
  coordsTest = runifsq(nTest * nDataSets, c(-1/3, 1/3), c(-1/3, 1/3))
  xValuesTrain = matrix(coordsTrain[,1], ncol=nDataSets)
  yValuesTrain = matrix(coordsTrain[,2], ncol=nDataSets)
  xValuesTest = matrix(coordsTest[,1], ncol=nDataSets)
  yValuesTest = matrix(coordsTest[,2], ncol=nDataSets)
  xValues = rbind(xValuesTrain, xValuesTest)
  yValues = rbind(yValuesTrain, yValuesTest)
  
  # preallocate observation matrix, and pregenerate standard normal draws
  observations = matrix(nrow=nTotal+nx*ny, ncol=nDataSets)
  zsims = matrix(rnorm((nTotal + nx*ny) * nDataSets), ncol=nDataSets)
  
  # generate spatial component of observation values
  for(i in 1:nDataSets) {
    if(i %% printEvery == 0 || i == 1)
      print(paste0("Simulating data set ", i, "/", nDataSets))
    
    thisx = xValues[,i]
    thisy = yValues[,i]
    L = t(chol(corFun(cbind(c(thisx, xValuesGrid), c(thisy, yValuesGrid)))))
    observations[,i] = L %*% zsims[,i]
  }
  
  # scale by marginal standard deviation and add in error variance
  observations = observations * sqrt(marginalVar) + matrix(rnorm((nTotal + nx*ny) * nDataSets, sd=sqrt(errorVar)), ncol=nDataSets)
  
  # separate out test and train results
  trainI = 1:(nTotal - nTest)
  testI = (nTotal - nTest + 1):nTotal
  gridI = (nTotal + 1):(nTotal + nx * ny)
  if(nTotal - nTest != 0) {
    xTrain = xValues[trainI,]
    yTrain = yValues[trainI,]
    zTrain = observations[trainI,]
  } else {
    xTrain = c()
    yTrain = c()
    zTrain = c()
  }
  if(nTest != 0) {
    xTest = xValues[testI,]
    yTest = yValues[testI,]
    zTest = observations[testI,]
  }
  else {
    xTest = c()
    yTest = c()
    zTest = c()
  }
  # get grid results if necessary
  if(doPredGrid) {
    xGrid = xValuesGrid
    yGrid = yValuesGrid
    zGrid = observations[gridI,]
  } else {
    xGrid = NULL
    yGrid = NULL
    zGrid = NULL
  }
  
  # put relevant values into a list
  out = list(xTrain=xTrain, yTrain=yTrain, zTrain=zTrain, xTest=xTest, yTest=yTest, zTest=zTest, 
             xGrid=xGrid, yGrid=yGrid, zGrid=zGrid, xValuesGrid=xValuesGrid, yValuesGrid=yValuesGrid, nx=nx, ny=ny, 
             corFun=corFun, marginalVar=marginalVar, errorVar=errorVar, xRange=xRange, yRange=yRange)
  
  # plot the results
  plotExampleDataSets(out, saveDataSetPlot=saveDataSetPlot, plotNameRoot=plotNameRoot, fileNameRoot=fileNameRoot)
  
  # return the results
  out
}

# try out several bases, and determine the number of bases elements and their resolution
testBasisResolution = function() {
  latticeInfo = makeLatGridsKenya(nLayer=3, NC=28, nBuffer=5)
  print(paste0("NC=", 28, ", nLayer=", 3))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
  
  latticeInfo = makeLatGridsKenya(nLayer=3, NC=14, nBuffer=5)
  print(paste0("NC=", 14, ", nLayer=", 3))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
  
  latticeInfo = makeLatGridsKenya(nLayer=2, NC=54, nBuffer=5)
  print(paste0("NC=", 54, ", nLayer=", 2))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
  
  latticeInfo = makeLatGridsKenya(nLayer=1, NC=107, nBuffer=5)
  print(paste0("NC=", 107, ", nLayer=", 1))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
  
  latticeInfo2 = makeLatGridsKenya(nLayer=1, NC=30, nBuffer=5)
  print(paste0("NC=", 30, ", nLayer=", 1))
  print(paste0("lattice width=", min(sapply(latticeInfo2, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo2, function(x) {x$nx * x$ny}))))
  
  latticeInfo2 = makeLatGridsKenya(nLayer=1, NC=14, nBuffer=5)
  print(paste0("NC=", 14, ", nLayer=", 1))
  print(paste0("lattice width=", min(sapply(latticeInfo2, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo2, function(x) {x$nx * x$ny}))))
  
  test = c(latticeInfo2, latticeInfo)
  
  # NC=c(30, 107), 13725
  
  invisible(NULL)
}

testBasisResolution2 = function(nLayer, NC, nBuffer=5) {
  latticeInfo = makeLatGridsKenya(nLayer=nLayer, NC=NC, nBuffer=nBuffer)
  print(paste0("NC=", NC, ", nLayer=", nLayer))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
}

testBasisResolution3 = function(nLayer, NC, nBuffer=5) {
  latticeInfo = makeLatGrids(c(-1,1), c(-1,1), nLayer=nLayer, NC=NC, nBuffer=nBuffer)
  print(paste0("NC=", NC, ", nLayer=", nLayer))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
}










