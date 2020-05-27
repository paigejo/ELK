# this script contains some miscellaneous, but useful functions

logit <- function(x) {
  log(x/(1-x))
}

expit <- function(x) {
  res = exp(x)/(1+exp(x))
  res[x > 100] = 1
  res[x < -100] = 0
  res
}

# Do precomputations for computing precision matrix for a single layer or a block diagonal sparse 
# precision matrix for multiple layers
# kappa: scale of Matern covariance with smoothness 1
# xRange, yRange: x and y intervals in space over which basis elements are placed
# nx, ny: number of basis elements when counting along the lattive in 
#         x and y directions respectively (for the first layer)
# rho: sill (marginal variance) for first layer
# nLayer: number of lattice layers
# thisLayer: user should always set to 1 to get full block diagonal precision matrix 
#            for all layers
# alphas: weights on the variances of each layer.  Scales rho for each layer.
# fastNormalize: simple normalization to make marginal variance = rho in center. Basis coefficients may have different variances
# assume there is a single kappa, rho is adjusted with alpha when there are multiple layers
# latticeInfo: an object returned by makeLatGrids
precomputationsQ = function(latticeInfo, thisLayer=1) {
  
  # save layer info quantities
  nLayer = length(latticeInfo)
  nx = latticeInfo[[thisLayer]]$nx
  ny = latticeInfo[[thisLayer]]$ny
  xRange=latticeInfo[[thisLayer]]$xRangeKnots
  yRange=latticeInfo[[thisLayer]]$yRangeKnots
  knotPts = latticeInfo[[thisLayer]]$latCoords
  
  # make (Laplacian) differential operators
  Dnx = bandSparse(nx, k=0:1, diag=list(rep(-2, nx), rep(1, nx-1)), symmetric=TRUE)
  Dny = bandSparse(ny, k=0:1, diag=list(rep(-2, ny), rep(1, ny-1)), symmetric=TRUE)
  
  # generate x and y (Laplacian) differential operators
  Inx = Diagonal(n=nx)
  Iny = Diagonal(n=ny)
  Bx = kronecker(Iny, Dnx)
  By = kronecker(Dny, Inx)
  Bxy = Bx + By
  
  ## now construct relevant row of A for value at midpoint at this layer, and Q matrix
  
  # make mid lattice point
  # xi = ceiling(nx/2)
  # yi = ceiling(ny/2)
  # midPt = matrix(c(seq(xRange[1], xRange[2], l=nx)[xi], seq(yRange[1], yRange[2], l=ny)[yi]), nrow=1)
  midPt = matrix(c((xRange[1] + xRange[2])/2, (yRange[1] + yRange[2])/2), nrow=1)
  
  Ai = makeA(midPt, latticeInfo, thisLayer=thisLayer, maxLayer=thisLayer)
  
  # return results
  # If multiple layers, return block diagonal sparse matrix
  if(thisLayer == nLayer) {
    return(list(list(Bxy=Bxy, Ai=Ai)))
  }
  else {
    return(c(list(list(Bxy=Bxy, Ai=Ai)), precomputationsQ(latticeInfo, thisLayer+1)))
  }
}

# Do precomputations for computing precision matrix for a single layer or a block diagonal sparse 
# precision matrix for multiple layers. Return pre- constructed block diagonal matrices instead of 
# individual matrices going into the block diagonal
# kappa: scale of Matern covariance with smoothness 1
# xRange, yRange: x and y intervals in space over which basis elements are placed
# nx, ny: number of basis elements when counting along the lattive in 
#         x and y directions respectively (for the first layer)
# rho: sill (marginal variance) for first layer
# nLayer: number of lattice layers
# thisLayer: user should always set to 1 to get full block diagonal precision matrix 
#            for all layers
# alphas: weights on the variances of each layer.  Scales rho for each layer.
# fastNormalize: simple normalization to make marginal variance = rho in center. Basis coefficients may have different variances
# assume there is a single kappa, rho is adjusted with alpha when there are multiple layers
# latticeInfo: an object returned by makeLatGrids
precomputationsQ2 = function(latticeInfo) {
  
  out = precomputationsQ(latticeInfo, 1)
  Bxys = lapply(out, function(x) {x$Bxy})
  Ais = lapply(out, function(x) {x$Ai})
  Bxy = bdiag(Bxys)
  mBxymBxyT = -Bxy - t(Bxy) # - Bxy - t(Bxy)
  BxyTBxy = t(Bxy) %*% Bxy # t(Bxy) %*% Bxy
  A = bdiag(Ais)
  At = t(A)
  ms = sapply(out, function(x) {length(x$Ai)})
  
  list(mBxymBxyT=mBxymBxyT, BxyTBxy=BxyTBxy, A=A, At=At, ms=ms)
}

# precompute normalization constants using natural smoothing spline on log log scale
precomputeNormalization = function(xRangeDat=c(-1,1), yRangeDat=c(-1,1), effRangeRange=NULL, nLayer=3, NC=13, 
                                   nBuffer=5, nKnots=NULL, saveResults=FALSE, doFinalTest=!saveResults, 
                                   latticeInfo=NULL, plotNormalizationSplines=TRUE) {
  # construct lattice info if necessary, precompute relevant matrices
  if(is.null(latticeInfo))
    latticeInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  else {
    xRangeDat = latticeInfo[[1]]$xRangeDat
    yRangeDat = latticeInfo[[1]]$yRangeDat
    nLayer = length(latticeInfo)
    NC = latticeInfo[[1]]$NC
    nBuffer = latticeInfo[[1]]$nBuffer
  }
  precomputationsQ = precomputationsQ2(latticeInfo)
  
  # if there are any missing layers, then we know that we are allowing for multiple kappas
  latticeWidths = sapply(latticeInfo, function(x) {x$latWidth})
  if(any(latticeWidths[2:length(latticeWidths)] != latticeWidths[1:(length(latticeWidths)-1)] / 2))
    singleKappa = FALSE
  else
    singleKappa = TRUE
  
  # set the range of effRange if necessary to be between half of the finest layer's lattice width and the data domain diameter
  if(is.null(effRangeRange))
    effRangeRange = c(latticeInfo[[nLayer]]$latWidth / 5, max(c(diff(xRangeDat), diff(yRangeDat))))
  
  # set the number of knot points so that the second to largest point is roughly 95% of the maximum knot point if necessary
  if(is.null(nKnots)) {
    width = abs(log(.95))
    nKnots = ceiling(diff(log(effRangeRange))/width) + 1
  }
  
  # set the values of effRange between which we want to interpolate
  effRangeKnots = exp(seq(log(effRangeRange[1]), log(effRangeRange[2]), l=nKnots))
  if(singleKappa) {
    kappaKnots = sqrt(8)/effRangeKnots * latticeInfo[[1]]$latWidth
  } else {
    effRangeKnots = matrix(effRangeKnots, nrow=1)
    kappaKnots = outer(sapply(latticeInfo, function(x) {x$latWidth}), c(sqrt(8)/effRangeKnots), "*")
  }
  
  # compute ctilde vector for each value of effective range. The ctilde value for one layer is independent of the 
  # kappa value in another layer. Hence, only univariate splines are necessary to precompute
  # NOTE: multiply these by alphas = c(1/nLayer, ..., 1/nLayer) now, then divide by alpha later depending on alpha
  if(singleKappa) {
    ctildes = sapply(kappaKnots, makeQPrecomputed, precomputedMatrices=precomputationsQ, latticeInfo=latticeInfo, 
                     alphas=rep(1/nLayer, nLayer), normalized=TRUE, fastNormalize=TRUE, returnctildes=TRUE) / nLayer
  } else {
    ctildes = apply(kappaKnots, 2, makeQPrecomputed, precomputedMatrices=precomputationsQ, latticeInfo=latticeInfo, 
                    alphas=rep(1/nLayer, nLayer), normalized=TRUE, fastNormalize=TRUE, returnctildes=TRUE) / nLayer
  }
  
  # estimate splines:
  # a^T %*% Q^(-1) %*% a /(rho * alpha) = ctildes
  getSplineFun = function(cts) {
    logFun = splinefun(log(effRangeKnots), log(cts), method = "hyman")
    function(x, alpha) {exp(logFun(log(x)))/alpha}
  }
  funs = apply(ctildes, 1, getSplineFun)
  fullFun = function(effRange, alphas) {
    # sapply(alphas, funs, x=effRange)
    res = numeric(nLayer)
    for(i in 1:nLayer) {
      if(length(effRange) == 1)
        res[i] = funs[[i]](effRange, alphas[i])
      else
        res[i] = funs[[i]](effRange[i], alphas[i])
    }
    res
  }
  
  # plot functions
  for(i in 1:nLayer) {
    effRanges = seq(effRangeRange[1], effRangeRange[2], l=500)
    thisctildes = funs[[i]](effRanges, alpha=1/nLayer)
    if(plotNormalizationSplines) {
      par(mfrow=c(1,1))
      plot(effRanges, thisctildes, type="l", col="blue", main=paste0("Layer ", i), xlab="Effective Range", ylab="ctilde")
      points(effRangeKnots, ctildes[i,]*nLayer, pch=19, cex=.3)
    }
  }
  
  for(i in 1:nLayer) {
    effRanges = seq(effRangeRange[1], effRangeRange[2], l=500)
    thisctildes = funs[[i]](effRanges, alpha=1/nLayer)
    if(plotNormalizationSplines) {
      par(mfrow=c(1,1))
      plot(log(effRanges), log(thisctildes), type="l", col="blue", main=paste0("Layer ", i), xlab="Log Effective Range", ylab="Log ctilde")
      points(log(effRangeKnots), log(ctildes[i,]*nLayer), pch=19, cex=.3)
    }
  }
  
  if(doFinalTest) {
    # the true value of alphas doesn't matter as long as it is different than what was used in the precomputation
    alphas = getAlphas(nLayer)
    
    if(singleKappa) {
      i = round(nKnots/2)
      thisKappa = kappaKnots[i]
      thisEffRange = effRangeKnots[i]
      thisctildes = fullFun(thisEffRange, alphas)
    } else {
      i = round(nKnots/2)
      thisKappa = kappaKnots[,i]
      thisEffRange = effRangeKnots[i]
      thisctildes = fullFun(thisEffRange, alphas)
    }
    Q = makeQPrecomputed(kappa=thisKappa, precomputedMatrices=precomputationsQ, latticeInfo=latticeInfo, 
                         alphas=alphas, normalized=TRUE, fastNormalize=TRUE)
    Q2 = makeQPrecomputed(kappa=thisKappa, precomputedMatrices=precomputationsQ, latticeInfo=latticeInfo, 
                          alphas=alphas, normalized=TRUE, fastNormalize=TRUE, ctildes=thisctildes)
    if(plotNormalizationSplines)
      print(mean(abs(Q - Q2)))
  }
  
  # save functions
  if(saveResults) {
    allNCs = sapply(latticeInfo, function(x) {x$NC})
    if(length(allNCs) == 1)
      ncText = paste0("_NC", allNCs)
    else {
      tempText = do.call("paste0", as.list(c(allNCs[1], paste0("_", allNCs[-1]))))
      ncText = paste0("_NC", tempText)
    }
    
    save(list(funs=funs, fullFun=fullFun, latInfo=latticeInfo), 
         file=paste0("ctildeSplines_nLayer", nLayer, ncText, 
                     "_xmin", round(xRangeDat[1], 1), "_xmax", round(xRangeDat[2], 1), 
                     "_ymin", round(yRangeDat[1], 1), "_ymax", round(yRangeDat[2], 1), 
                     ".RData"))
  }
  
  list(funs=funs, fullFun=fullFun)
}

# get the marginal variance for multi-resolution process
# tod: theta/delta, or theta/latticeWidth
# either nu or alphas must be non-null
getMultiMargVar = function(kappa=1, rho=1, tod=2.5, nLayer=3, nu=NULL, alphas=NULL, nx=NULL, ny=NULL, 
                           xRange=c(0,1), yRange=xRange, xRangeDat=c(-2,1), yRangeDat=xRangeDat, nBuffer=5) {
  # set alphas if nu has been set
  if(!is.null(nu)) {
    alphas = getAlphas(nLayer, nu)
  }
  
  # set nx and ny if necessary and add buffer to avoid edge effects
  if(is.null(nx) || is.null(ny)) {
    maxPt = ceiling(tod)*4 + 1
    nx = maxPt
    ny = maxPt
  }
  
  # generate knot lattice locations and filter out locations 
  # too far outside of the data domain
  origNX = nx
  origNY = ny
  knotXs = seq(xRange[1], xRange[2], l=nx)
  knotYs = seq(yRange[1], yRange[2], l=ny)
  if(sum(knotXs > xRangeDat[2]) > nBuffer)
    knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
  if(sum(knotXs < xRangeDat[1]) > nBuffer)
    knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
  if(sum(knotYs > yRangeDat[2]) > nBuffer)
    knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
  if(sum(knotYs < yRangeDat[1]) > nBuffer)
    knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
  nx = length(knotXs)
  ny = length(knotYs)
  
  # sum the variances of each layer weighted by alphas
  totalMargVar = c()
  for(l in 1:nLayer) {
    # get the layer marginal variances
    layerMargVar = as.numeric(getMargVar(kappa, rho, tod, origNX*2^(l-1), origNY*2^(l-1), xRange=xRange, yRange=yRange, 
                                         xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer))
    layerMargVar[1:2] = layerMargVar[1:2]*alphas[l]
    if(l == 1)
      totalMargVar = layerMargVar[1:2]
    else
      totalMargVar = totalMargVar + layerMargVar[1:2]
  }
  
  # add in a variance ratio column
  totalMargVar = c(totalMargVar, totalMargVar[1]/totalMargVar[2])
  names(totalMargVar) = c("actualVar", "theorVar", "inflation")
  
  totalMargVar
}

# compute the marginal variance for a given resolution layer
# tod: theta/delta, or theta/latticeWidth
getMargVar = function(kappa=1, rho=1, tod=2.5, nx=NULL, ny=NULL, xRange=c(-1,2), yRange=xRange, 
                      xRangeDat=c(0,1), yRangeDat=xRangeDat, nBuffer=5) {
  # set nx and ny if necessary and add buffer to avoid edge effects
  if(is.null(nx) || is.null(ny)) {
    maxPt = ceiling(tod)*4 + 1
    nx = maxPt
    ny = maxPt
  }
  
  # generate knot lattice locations and filter out locations 
  # too far outside of the data domain
  knotXs = seq(xRange[1], xRange[2], l=nx)
  knotYs = seq(yRange[1], yRange[2], l=ny)
  delta = knotXs[2]-knotXs[1]
  if(sum(knotXs > xRangeDat[2]) > nBuffer)
    knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
  if(sum(knotXs < xRangeDat[1]) > nBuffer)
    knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
  if(sum(knotYs > yRangeDat[2]) > nBuffer)
    knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
  if(sum(knotYs < yRangeDat[1]) > nBuffer)
    knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
  knotPts = make.surface.grid(list(x=knotXs, y=knotYs))
  
  # take the middle knot location
  midPt = matrix(c(knotXs[ceiling(length(knotXs)/2)], 
                   knotYs[ceiling(length(knotYs)/2)]), nrow=1)
  
  # compute variance of process at the middle knot
  A = as.matrix(makeA(midPt, xRange, nx, yRange, ny, tod*delta, 
                      xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer))
  Q = makeQ(kappa, rho, xRange, yRange, nx, ny, 
            xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer)
  # VarC = as.matrix(solve(Q))
  varMidPt = A %*% inla.qsolve(Q, t(A))
  
  # compare with theoretical marginal variance
  sigma2 = rho/(4*pi * kappa^2)
  inflation = varMidPt/sigma2
  
  # return results
  list(actualVar=varMidPt, theorVar=sigma2, inflation=inflation)
}


# estimates effective range for a given LatticeKrig model
getEffRange = function(predPts=NULL, xRangeKnot=c(0,1), xNKnot=10, yRangeKnot=c(0,1), yNKnot=10, theta=NULL, nLayer=1, thisLayer=1, 
                       xRangeDat=xRangeKnot, yRangeDat=yRangeKnot, nBuffer=5, mx=20, my=20) {
  
}

##### simulate from a latticeKrig model
# return a function with argument nsim for simulating a number of realizations from a latticeKrig model with no nugget.  
# coords: coordinates at which to simulate
# all other arguments: same as general latticeKrig arguments
LKSimulator = function(coords, NC=5, kappa=1, rho=1, nu=1.5, nBuffer=5, nLayer=3, normalize=TRUE) {
  warning("LKSimulator is not tested")
  # first make the grid on which to set the basis functions
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  knotGrid = makeLatGrid(xRange=xRangeDat, yRange=yRangeDat, NC=NC, nBuffer=nBuffer)
  xRangeKnots=knotGrid$xRangeKnots
  nx=knotGrid$nx
  yRangeKnots=knotGrid$yRangeKnots
  ny=knotGrid$ny
  
  # generate layer variance weights
  alphas = getAlphas(nLayer, nu)
  
  # generate basis function and precision matrices
  AObs = makeA(coords, xRangeKnots, nx, yRangeKnots, ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer)
  Q = makeQ(kappa=kappa, rho=rho, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, ny=ny, 
            nLayer=nLayer, alphas=alphas, xRangeDat=xRangeDat, yRangeDat=yRangeDat, 
            nBuffer=nBuffer, normalized = normalize)
  L = as.matrix(t(chol(solve(Q))))
  zsim = matrix(rnorm(nrow(Q)), ncol=1)
  fieldSim = L %*% zsim
  
  fieldSims
}

# return a function with argument nsim for simulating a number of realizations from a latticeKrig model with no nugget.  
# same as LKSimulator, but uses marginal variance, effective range parameterization instead of rho, kappa.
# coords: coordinates at which to simulate
# all other arguments: same as general latticeKrig arguments
LKSimulator2 = function(coords, nsim=1, NC=5, effRange=(max(coords[,1])-min(coords[,1]))/3, margVar=1, 
                        nu=1.5, nBuffer=5, nLayer=3, normalize=TRUE) {
  warning("LKSimulator2 is not tested")
  
  # first make the grid on which to set the basis functions
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  knotGrid = makeLatGrid(xRange=xRangeDat, yRange=yRangeDat, NC=NC, nBuffer=nBuffer)
  xRangeKnots=knotGrid$xRangeKnots
  nx=knotGrid$nx
  yRangeKnots=knotGrid$yRangeKnots
  ny=knotGrid$ny
  
  # convert from effRange, margVar to rho, kappa
  latticeWidth = (xRangeKnots[2] - xRangeKnots[1])/(nx-1)
  kappa = sqrt(8)/effRange * latticeWidth
  
  # since we are normalizing the process, rho is just sigmaSq
  rho = margVar
  
  # generate layer variance weights
  alphas = getAlphas(nLayer, nu)
  
  # generate basis function and precision matrices
  AObs = makeA(coords, xRangeKnots, nx, yRangeKnots, ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer)
  Q = makeQ(kappa=kappa, rho=rho, xRange=xRangeKnots, yRange=yRangeKnots, nx=nx, ny=ny, 
            nLayer=nLayer, alphas=alphas, xRangeDat=xRangeDat, yRangeDat=yRangeDat, 
            nBuffer=nBuffer, normalized = normalize)
  L = as.matrix(t(chol(solve(Q))))
  zsims = matrix(rnorm(nrow(Q)*nsim), ncol=nsim)
  fieldSims = matrix(as.numeric(AObs %*% L %*% zsims), ncol=nsim)
  
  fieldSims
}

simSPDE = function(coords, nsim=1, mesh=NULL, effRange=(max(coords[,1])-min(coords[,1]))/3, margVar=1) {
  # generate mesh grid if necessary
  if(is.null(mesh)) {
    mesh = getSPDEMeshGrid(coords, doPlot = FALSE)
  }
  
  # calculate SPDE model parameters based on Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
  meshSize <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  # it is easier to use theta and set sigma0 to 1 then to set sigma0 and the effective range directly
  # kappa0 <- sqrt(8)/effRange * meshSize # since nu = 1
  # kappa0 <- sqrt(8)/effRange # since nu = 1
  # kappa0 = sqrt(8) / 5
  # logKappa = log(kappa0)
  sigma0 = 1
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # logTau = log(tau0)
  
  # from page 5 of the paper listed above:
  logKappa = 0.5 * log(8)
  logTau = 0.5 * (lgamma(1) - (lgamma(2) + log(4*pi))) - logKappa
  theta = c(log(sqrt(margVar)), log(effRange))
  spde <- inla.spde2.matern(mesh, B.tau = cbind(logTau, -1, +1),
                            B.kappa = cbind(logKappa, 0, -1), theta.prior.mean = theta,
                            theta.prior.prec = c(0.1, 1))
  
  # generate A and Q precision matrix
  Q = inla.spde2.precision(spde, theta = theta)
  A = inla.spde.make.A(mesh, coords)
  
  # generate simulations
  simField = inla.qsample(nsim, Q)
  simDat = as.matrix(A %*% simField)
  
  simDat
}

makeQSPDE = function(mesh, effRange, margVar=1) {
  
  # calculate SPDE model parameters based on Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
  # it is easier to use theta and set sigma0 to 1 then to set sigma0 and the effective range directly
  # kappa0 <- sqrt(8)/effRange * meshSize # since nu = 1
  # kappa0 <- sqrt(8)/effRange # since nu = 1
  # kappa0 = sqrt(8) / 5
  # logKappa = log(kappa0)
  sigma0 = 1
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # logTau = log(tau0)
  
  # from page 5 of the paper listed above:
  logKappa = 0.5 * log(8)
  logTau = 0.5 * (lgamma(1) - (lgamma(2) + log(4*pi))) - logKappa
  theta = c(log(sqrt(margVar)), log(effRange))
  spde <- inla.spde2.matern(mesh, B.tau = cbind(logTau, -1, +1),
                            B.kappa = cbind(logKappa, 0, -1), theta.prior.mean = theta,
                            theta.prior.prec = c(0.1, 1))
  
  # generate Q precision matrix
  inla.spde2.precision(spde, theta = theta)
}

# computing precision matrix for a single layer or a block diagonal sparse 
# precision matrix for multiple layers
# kappa: scale of Matern covariance with smoothness 1
# xRange, yRange: x and y intervals in space over which basis elements are placed
# nx, ny: number of basis elements when counting along the lattive in 
#         x and y directions respectively (for the first layer)
# rho: sill (marginal variance) for first layer
# nLayer: number of lattice layers
# thisLayer: user should always set to 1 to get full block diagonal precision matrix 
#            for all layers
# alphas: weights on the variances of each layer.  Scales rho for each layer.
# fastNormalize: simple normalization to make marginal variance = rho in center. Basis coefficients may have different variances
# assume there is a single kappa, rho is adjusted with alpha when there are multiple layers
# latticeInfo: an object returned by makeLatGrids
makeQ = function(kappa=1, rho=1, latticeInfo, thisLayer=1, alphas=NULL, nu=NULL, 
                 normalized=FALSE, fastNormalize=FALSE, precomputedMatrices=NULL) {
  require(Matrix)
  require(spam)
  require(fields)
  
  # save base layer input quantities
  origRho = rho
  nLayer = length(latticeInfo)
  
  # make alphas according to nu relation if nu is set
  if(is.null(alphas) && !is.null(nu)) {
    alphas = getAlphas(nLayer, nu)
  }
  else if(is.null(nu) && nLayer != 1 && is.null(alphas)) {
    warning("Both alphas and nu are NULL. Defaulting to exponential covariance.")
    nu = 0.5
    alphas = getAlphas(nLayer, nu)
  }
  
  nx = latticeInfo[[thisLayer]]$nx
  ny = latticeInfo[[thisLayer]]$ny
  if(is.null(precomputedMatrices)) {
    # make (Laplacian) differential operators
    Dnx = bandSparse(nx, k=0:1, diag=list(rep(-2, nx), rep(1, nx-1)), symmetric=TRUE)
    Dny = bandSparse(ny, k=0:1, diag=list(rep(-2, ny), rep(1, ny-1)), symmetric=TRUE)
    
    # generate x and y (Laplacian) differential operators
    Inx = Diagonal(n=nx)
    Iny = Diagonal(n=ny)
    Bx = kronecker(Iny, Dnx)
    By = kronecker(Dny, Inx)
    Bxy = Bx + By
  }
  else
    Bxy = precomputedMatrices[[thisLayer]]$Bxy
  
  # make B, SAR regression matrix for Bc = e
  B = Diagonal(n=nx*ny, x=kappa^2) - Bxy
  
  # compute precision matrix
  Q = t(B) %*% B
  
  if(normalized) {
    
    # now construct relevant row of A for value at midpoint at this layer, and Q matrix
    if(is.null(precomputedMatrices)) {
      xRange=latticeInfo[[thisLayer]]$xRangeKnots
      yRange=latticeInfo[[thisLayer]]$yRangeKnots
      
      # make mid lattice point
      # xi = ceiling(nx/2)
      # yi = ceiling(ny/2)
      # midPt = matrix(c(seq(xRange[1], xRange[2], l=nx)[xi], seq(yRange[1], yRange[2], l=ny)[yi]), nrow=1)
      midPt = matrix(c((xRange[1] + xRange[2])/2, (yRange[1] + yRange[2])/2), nrow=1)
      
      Ai = makeA(midPt, latticeInfo, thisLayer=thisLayer, maxLayer=thisLayer)
    }
    else
      Ai = precomputedMatrices[[thisLayer]]$Ai
    
    # # test
    # sds2 = 1/diag(Q)
    # sdMat = Diagonal(x=sds2)
    # # Qnorm2 = sweep(sweep(Q, 1, sds2, "*"), 2, sds2, "*")
    # Qnorm2 = sdMat %*% Q %*% sdMat
    # 
    # # QnormInv = sweep(sweep(Qinv, 1, 1/sds), 2, 1/sds)
    # procVar2 = as.numeric(Ai %*% inla.qsolve(Qnorm2, t(Ai)))
    # # procVar = as.numeric(Ai %*% QnormInv %*% t(Ai))
    # Q2 = Qnorm2 * (procVar2 / rho)
    # # Q = Q2 # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE, newnormalize=TRUE)): ~5.8s
    
    #  test 2
    if(fastNormalize) {
      ctilde = as.numeric(Ai %*% inla.qsolve(Q, t(Ai)))
      Qtilde = ctilde * Q
      Q = Qtilde # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE, newnormalize=TRUE)): 2.72
    }
    else {
      # renormalize basis coefficients to have constant variance, and the process to have unit variance
      Qinv = inla.qsolve(Q, diag(nrow(Q)))
      sds = sqrt(diag(Qinv))
      sdMat = Diagonal(x=sds)
      # Qnorm = sweep(sweep(Q, 1, sds, "*"), 2, sds, "*")
      Qnorm = sdMat %*% Q %*% sdMat
      # QnormInv = sweep(sweep(Qinv, 1, 1/sds), 2, 1/sds)
      procVar = as.numeric(Ai %*% inla.qsolve(Qnorm, t(Ai)))
      # procVar = as.numeric(Ai %*% QnormInv %*% t(Ai))
      Q = Qnorm * procVar # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE) ~5.87
    }
    # # compute how good an approximation it was
    # hist(diag(solve(Q)))
    # hist(diag(solve(Q2)))
    # hist(diag(solve(Qtilde)))
    # image(Q)
    # image(Q2)
    # image(Qtilde)
  }
  if(nLayer == 1 && is.null(alphas))
    alphas = 1
  Q = Q * (1/alphas[thisLayer])
  
  # return results
  # If multiple layers, return block diagonal sparse matrix
  if(thisLayer == nLayer) {
    if(thisLayer == 1)
      Q = Q * (1/rho)
    return(Q)
  }
  else if((thisLayer == 1) && (nLayer != 1)) {
    Q = bdiag(c(list(Q), makeQ(kappa, origRho, latticeInfo, thisLayer+1, alphas, normalized=normalized, 
                               fastNormalize=fastNormalize)))
    Q = Q * (1/rho)
    return(Q)
  }
  else {
    return(c(list(Q), makeQ(kappa, origRho, latticeInfo, thisLayer+1, alphas, normalized=normalized, 
                            fastNormalize=fastNormalize)))
  }
}

meanSegmentLength = function(mesh, filterLargerThan=NULL) {
  t.sub = 1:nrow(mesh$graph$tv)
  idx = cbind(mesh$graph$tv[t.sub, c(1:3, 1), drop = FALSE], NA)
  x = mesh$loc[t(idx), 1]
  y = mesh$loc[t(idx), 2]
  indices = 1:4 + rep(seq(from=0, to=length(x)-5, by=5), each=4)
  segx = x[indices]
  segy = y[indices]
  coords = cbind(segx, segy)
  dists = rdist.vec(coords[1:(length(segx) - 1),], coords[2:length(segx),])
  dists = dists[-seq(from=4, to=length(dists), by=4)]
  
  if(!is.null(filterLargerThan))
    dists = dists[dists <= filterLargerThan]
  
  mean(dists)
}

LKINLA.cov = function(x1, x2, latticeInfo, kappa, alphas, rho=1, normalize=TRUE, 
                      fastNormalize=TRUE, 
                      precomputedMatrices=NULL, precomputedA1=NULL, precomputedA2=NULL, 
                      precomputationsFileNameRoot="") {
  ctildes = NULL
  
  if(precomputationsFileNameRoot != "") {
    # load in precomputations
    load(paste0("savedOutput/precomputations/", precomputationsFileNameRoot, ".RData"))
    
    if(!is.null(precomputedNormalizationFun)) {
      if(length(kappa) == 1)
        effectiveCor = sqrt(8) / kappa * latticeInfo[[1]]$latWidth
      else
        effectiveCor = sqrt(8) / kappa * sapply(latticeInfo, function(x) {x$latWidth})
      ctildes = precomputedNormalizationFun$fullFun(effectiveCor, alphas)
    }
    
  }
  if(!is.null(precomputedMatrices)) {
    Q = makeQPrecomputed(precomputedMatrices, kappa, rho, latticeInfo, alphas, normalized=normalize, 
                         fastNormalize=fastNormalize, ctildes=ctildes)
  } else {
    Q = makeQ(kappa, rho, latticeInfo, 1, alphas, normalized=normalize, fastNormalize=fastNormalize)
  }
  
  if(is.null(precomputedA1))
    A1 = makeA(x1, latticeInfo, 1)
  else
    A1 = precomputedA1
  if(is.null(precomputedA2))
    A2 = makeA(x2, latticeInfo, 1)
  else
    A2 = precomputedA2
  
  if(any(alphas == 0)) {
    # in this case, we must be careful to only consider the layers with nonzero weights
    includeI = rep(TRUE, nrow(Q))
    thisI = 1
    for(i in 1:length(latticeInfo)) {
      startI = thisI
      endI = startI - 1 + latticeInfo[[i]]$nx * latticeInfo[[i]]$ny
      if(alphas[i] == 0) {
        includeI[startI:endI] = FALSE
      }
      thisI = endI + 1
    }
    A1 = matrix(as.matrix(A1)[,includeI], nrow=nrow(A1))
    A2 = matrix(as.matrix(A2)[,includeI], nrow=nrow(A2))
    Q = Q[includeI,includeI]
  }
  # use A1 %*% inla.qsolve(Q, t(A2)) instead?
  A1 %*% inla.qsolve(Q, t(A2))
}

LK.cov = function(x1, x2, LKinfo, a.wght, alphas, lambda, sigma, rho) {
  # domainCoords = LKinfo$latticeInfo$grid.info$range
  # nBuffer = LKinfo$latticeInfo$NC.buffer
  # NC = max(LKinfo$latticeInfo$mxDomain[1,])
  # LKrigSetup(domainCoords, nlevel=length(alphas), a.wght=a.wght, normalize=normalize, 
  #            lambda=lambda, sigma=sigma, rho=rho)
  if(length(a.wght) > 1)
    a.wght = as.list(a.wght)
  LKinfo = LKinfoUpdate(LKinfo, a.wght=a.wght, alpha=alphas, lambda=lambda, sigma=sigma, rho=rho)
  LKrig.cov(x1, x2, LKinfo)
}

# make method for calculating individual covariance function
getLKInlaCovarianceFun = function(kappa, rho, nuggetVar, alphas, NP=200, latticeInfo, normalize=TRUE, fastNormalize=TRUE, 
                                  precomputedMatrices=NULL, precomputedAcenter=NULL, precomputedAx=NULL, precomputedAy=NULL) {
  # generate test locations based on code from LKrig.cov.plot
  xlim <- latticeInfo[[1]]$xRangeDat
  ux <- seq(xlim[1], xlim[2], , NP)
  ylim <- latticeInfo[[1]]$yRangeDat
  uy <- seq(ylim[1], ylim[2], , NP)
  center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  # calculate covariances
  x1 <- cbind(ux, rep(center[2], NP))
  x2 <- rbind(center)
  d <- c(rdist(x1, x2))
  y <- as.numeric(LKINLA.cov(x1, x2, latticeInfo, kappa, alphas, rho, normalize, fastNormalize, 
                             precomputedMatrices, precomputedAx, precomputedAcenter))
  y[NP/2] = y[NP/2] + nuggetVar
  x1 <- cbind(rep(center[1], NP), uy)
  d2 <- c(rdist(x1, x2))
  y2 <- as.numeric(LKINLA.cov(x1, x2, latticeInfo, kappa, alphas, rho, normalize, fastNormalize, 
                              precomputedMatrices, precomputedAy, precomputedAcenter))
  y2[NP/2] = y2[NP/2] + nuggetVar
  
  # average x and y covariances
  sortXI = sort(d, index.return=TRUE)$ix
  d = d[sortXI]
  y = y[sortXI]
  sortYI = sort(d2, index.return=TRUE)$ix
  d2 = d2[sortYI]
  y2 = y2[sortYI]
  d = rowMeans(cbind(d, d2))
  y = rowMeans(cbind(y, y2))
  return(cbind(d=d, cov=y, cor=y * (1 / max(y))))
}

covarianceDistributionLKINLA = function(latticeInfo, kappaVals, rhoVals=rep(1, length(kappaVals)), nuggetVarVals=rep(0, length(kappaVals)), 
                                        alphaMat, maxSamples=100, significanceCI=.8, normalize=TRUE, fastNormalize=TRUE, seed=NULL, NP = 200, 
                                        precomputationsFileNameRoot="", maxRadius=NULL) {
  if(!is.null(seed))
    set.seed(seed)
  nLayer = length(latticeInfo)
  
  # get hyperparameter samples
  sampleI = sample(1:length(rhoVals), min(maxSamples, length(rhoVals)))
  if(!is.null(dim(kappaVals))) {
    kappaVals = kappaVals[,sampleI]
    separateRanges = TRUE
    effectiveRanges = sweep(sqrt(8)/kappaVals, 1, sapply(latticeInfo, function(x){x$latWidth}), "*")
    minRange = min(apply(effectiveRanges, 1, min))
    maxRange = max(apply(cbind(5 * sapply(latticeInfo, function(x){x$latWidth}), effectiveRanges), 1, max))
  } else {
    kappaVals = kappaVals[sampleI]
    separateRanges = FALSE
    effectiveRanges = sqrt(8) * latticeInfo[[1]]$latWidth / kappaVals
    minRange = min(effectiveRanges) / 2^(length(latticeInfo) - 1)
    maxRange = max(c(effectiveRanges, latticeInfo[[1]]$latWidth * 5))
  }
  rhoVals = rhoVals[sampleI]
  nuggetVarVals = nuggetVarVals[sampleI]
  alphaMat = matrix(alphaMat[,sampleI], ncol=length(sampleI))
  
  # generate test locations based on code from LKrig.cov.plot (modify sampling points to be the correct resolution)
  # xlim <- latticeInfo[[1]]$xRangeDat
  # ux <- seq(xlim[1], xlim[2], , NP)
  # ylim <- latticeInfo[[1]]$yRangeDat
  # uy <- seq(ylim[1], ylim[2], , NP)
  # center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  if(is.null(maxRadius))
    maxRadius = maxRange * 2
  minStep = minRange / 10
  ThisNP = 2 * maxRadius / minStep
  xlim <- latticeInfo[[1]]$xRangeDat
  ylim <- latticeInfo[[1]]$yRangeDat
  # centerX = mean(xlim)
  # widthX = xlim[2] - centerX
  # deltaX = min(widthX, maxRadius)
  # xlim = c(centerX - deltaX, centerX + deltaX)
  # ux <- seq(xlim[1], xlim[2], , NP)
  # ylim <- latticeInfo[[1]]$yRangeDat
  # centerY = mean(ylim)
  # widthY = ylim[2] - centerY
  # deltaY = min(widthY, maxRadius)
  # ylim = c(centerY - deltaY, centerY + deltaY)
  # uy <- seq(ylim[1], ylim[2], , NP)
  # center <- rbind(c(ux[NP/2], uy[NP/2]))
  centerX = mean(xlim)
  widthX = xlim[2] - centerX
  deltaX = min(widthX, maxRadius)
  centerY = mean(ylim)
  widthY = ylim[2] - centerY
  deltaY = min(widthY, maxRadius)
  delta = min(deltaX, deltaY)
  xlim = c(centerX - delta, centerX + delta)
  # ux <- seq(xlim[1], xlim[2], , NP)
  ux <- c(seq(xlim[1], centerX-minRange/100, l=NP/2), centerX, seq(centerX+minRange/100, xlim[2], l=NP/2))
  ylim = c(centerY - delta, centerY + delta)
  # uy <- seq(ylim[1], ylim[2], , NP)
  uy <- c(seq(ylim[1], centerY-minRange/100, l=NP/2), centerY, seq(centerY+minRange/100, ylim[2], l=NP/2))
  center <- rbind(c(centerX, centerY))
  
  # precompute relevant matrices
  Qprecomputations = precomputationsQ2(latticeInfo)
  Acenter = makeA(center, latticeInfo)
  Ax = makeA(cbind(ux, center[2]), latticeInfo)
  Ay = makeA(cbind(center[1], uy), latticeInfo)
  
  # make method for calculating individual covariance function
  getOneCovariance = function(parameters) {
    # get relevant parameters
    if(!separateRanges) {
      kappa = parameters[1]
      rho = parameters[2]
      nuggetVar = parameters[3]
      alphas = parameters[-(1:3)]
    } else {
      kappa = parameters[1:nLayer]
      rho = parameters[1 + nLayer]
      nuggetVar = parameters[2 + nLayer]
      alphas = parameters[-(1:(2 + nLayer))]
    }
    
    # calculate covariances
    x1 <- cbind(ux, rep(center[2], NP+1))
    x2 <- rbind(center)
    d <- c(rdist(x1, x2))
    y <- as.numeric(LKINLA.cov(x1, x2, latticeInfo, kappa, alphas, rho, normalize, fastNormalize, 
                               Qprecomputations, Ax, Acenter, precomputationsFileNameRoot=precomputationsFileNameRoot))
    y[NP/2+1] = y[NP/2+1] + nuggetVar
    x1 <- cbind(rep(center[1], NP+1), uy)
    d2 <- c(rdist(x1, x2))
    y2 <- as.numeric(LKINLA.cov(x1, x2, latticeInfo, kappa, alphas, rho, normalize, fastNormalize, 
                                Qprecomputations, Ay, Acenter, precomputationsFileNameRoot=precomputationsFileNameRoot))
    y2[NP/2+1] = y2[NP/2+1] + nuggetVar
    
    # average x and y covariances
    # sortXI = sort(d, index.return=TRUE)$ix
    # d = d[sortXI]
    # y = y[sortXI]
    # sortYI = sort(d2, index.return=TRUE)$ix
    # d2 = d2[sortYI]
    # y2 = y2[sortYI]
    # d = rowMeans(cbind(d, d2))
    # y = rowMeans(cbind(y, y2))
    d = c(0, rowMeans(cbind(rev(d[1:(NP/2)]),  d[(NP/2+2):length(d)],  rev(d2[1:(NP/2)]),  d2[(NP/2+2):length(d2)])))
    y = c(mean(y[NP/2+1], y2[NP/2+1]), rowMeans(cbind(rev(y[1:(NP/2)]),  y[(NP/2+2):length(y)],  rev(y2[1:(NP/2)]),  y2[(NP/2+2):length(y2)])))
    sortXI = sort(d, index.return=TRUE)$ix
    d = d[sortXI]
    y = y[sortXI]
    return(cbind(d=d, cov=y, cor=y * (1 / max(y))))
  }
  
  # calculate covariances for each sample from the posterior
  if(!separateRanges)
    parameterMat = cbind(kappaVals, rhoVals, nuggetVarVals, t(alphaMat))
  else
    parameterMat = cbind(t(kappaVals), rhoVals, nuggetVarVals, t(alphaMat))
  # browser()
  out = apply(parameterMat, 1, getOneCovariance)
  d = out[1:(NP/2+1),1]
  covMat = out[(NP/2+2):(2*(NP/2+1)),]
  corMat = out[(2*(NP/2+1)+1):(3*(NP/2+1)),]
  
  # calculate summary statistics
  meanCov = rowMeans(covMat)
  lowerCov = apply(covMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCov = apply(covMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCor = rowMeans(corMat)
  lowerCor = apply(corMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCor = apply(corMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  
  # return results
  list(d=d, 
       cov=meanCov, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
       cor=meanCor, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
}

covarianceDistributionLK = function(latticeInfo, alphaVals, lambdaVals, a.wghtVals, rhoVals, 
                                    maxSamples=100, significanceCI=.8, normalize=TRUE, seed=NULL) {
  NP = 200
  if(!is.null(seed))
    set.seed(seed)
  
  # get number of layers
  nLayer = latticeInfo$nlevel
  a.wghtVals = matrix(a.wghtVals, ncol=length(lambdaVals))
  
  # get hyperparameter samples
  nuggetVarVals = rhoVals * lambdaVals
  sampleI = sample(1:length(rhoVals), maxSamples)
  alphaMat = alphaVals[,sampleI]
  lambdaVals = lambdaVals[sampleI]
  a.wghtVals = matrix(a.wghtVals[,sampleI], ncol=maxSamples)
  nuggetVarVals = nuggetVarVals[sampleI]
  rhoVals = rhoVals[sampleI]
  
  # generate test locations based on code from LKrig.cov.plot
  xlim <- latticeInfo$latticeInfo$rangeLocations[,1]
  ux <- seq(xlim[1], xlim[2], , NP)
  ylim <- latticeInfo$latticeInfo$rangeLocations[,2]
  uy <- seq(ylim[1], ylim[2], , NP)
  center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  # make method for calculating individual covariance function
  getOneCovariance = function(parameters) {
    # get relevant parameters
    a.wght = parameters[1:nrow(a.wghtVals)]
    rho = parameters[nrow(a.wghtVals) + 1]
    nuggetVar = parameters[nrow(a.wghtVals) + 2]
    lambda = parameters[nrow(a.wghtVals) + 3]
    alphas = parameters[-(1:(nrow(a.wghtVals) + 3))]
    
    # calculate covariances
    x1 <- cbind(ux, rep(center[2], NP))
    x2 <- rbind(center)
    d <- c(rdist(x1, x2))
    y <- as.numeric(LK.cov(x1, x2, latticeInfo, a.wght, alphas, lambda, sqrt(nuggetVar), rho))
    y[NP/2] = y[NP/2] + nuggetVar
    x1 <- cbind(rep(center[1], NP), uy)
    d2 <- c(rdist(x1, x2))
    y2 <- as.numeric(LK.cov(x1, x2, latticeInfo, a.wght, alphas, lambda, sqrt(nuggetVar), rho))
    y2[NP/2] = y2[NP/2] + nuggetVar
    
    # calculate covariances excluding nugget
    y3 = y
    y3[NP/2] = y[NP/2] - nuggetVar
    y4 = y2
    y4[NP/2] = y2[NP/2] - nuggetVar
    
    # average x and y covariances
    sortXI = sort(d, index.return=TRUE)$ix
    d = d[sortXI]
    y = y[sortXI]
    y3 = y3[sortXI]
    sortYI = sort(d2, index.return=TRUE)$ix
    d2 = d2[sortYI]
    y2 = y2[sortYI]
    y4 = y4[sortYI]
    d = rowMeans(cbind(d, d2))
    y = rowMeans(cbind(y, y2))
    yNoNugget = rowMeans(cbind(y3, y4))
    return(cbind(d=d, cov=y, cor=y * (1 / max(y)), covNoNugget=yNoNugget, norNoNugget=yNoNugget * (1 / max(yNoNugget))))
  }
  
  # calculate covariances for each sample from the posterior
  parameterMat = cbind(t(a.wghtVals), rhoVals, nuggetVarVals, lambdaVals, t(alphaMat))
  # browser()
  out = apply(parameterMat, 1, getOneCovariance)
  d = out[1:200,1]
  covMat = out[201:400,]
  corMat = out[401:600,]
  covMatNoNugget = out[601:800,]
  corMatNoNugget = out[801:1000,]
  
  # calculate summary statistics
  meanCov = rowMeans(covMat)
  lowerCov = apply(covMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCov = apply(covMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCor = rowMeans(corMat)
  lowerCor = apply(corMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCor = apply(corMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCovNoNugget = rowMeans(covMatNoNugget)
  lowerCovNoNugget = apply(covMatNoNugget, 1, quantile, probs=(1-significanceCI)/2)
  upperCovNoNugget = apply(covMatNoNugget, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCorNoNugget = rowMeans(corMatNoNugget)
  lowerCorNoNugget = apply(corMatNoNugget, 1, quantile, probs=(1-significanceCI)/2)
  upperCorNoNugget = apply(corMatNoNugget, 1, quantile, probs=1 - (1-significanceCI)/2)
  
  # return results
  list(d=d, 
       cov=meanCov, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
       cor=meanCor, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat, 
       covNoNugget=meanCovNoNugget, upperCovNoNugget=upperCovNoNugget, lowerCovNoNugget=lowerCovNoNugget, covMatNoNugget=covMatNoNugget, 
       corNoNugget=meanCorNoNugget, upperCorNoNugget=upperCorNoNugget, lowerCorNoNugget=lowerCorNoNugget, corMatNoNugget=corMatNoNugget)
}

covarianceDistributionSPDE = function(effectiveRangeVals, rhoVals=rep(1, length(effectiveRangeVals)), nuggetVarVals=rep(0, length(rhoVals)), 
                                      mesh, maxSamples=100, significanceCI=c(.8, .95), seed=NULL, xRangeDat=NULL, yRangeDat=NULL, NP = 200, 
                                      maxRadius=NULL) {
  if(!is.null(seed))
    set.seed(seed)
  
  # get hyperparameter samples
  sampleI = sample(1:length(effectiveRangeVals), maxSamples)
  effectiveRangeVals = effectiveRangeVals[sampleI]
  rhoVals = rhoVals[sampleI]
  nuggetVarVals = nuggetVarVals[sampleI]
  
  # generate test locations based on code from LKrig.cov.plot
  if(is.null(xRangeDat) || is.null(yRangeDat)) {
    idx = unique(c(mesh$segm$int$idx[,1], mesh$segm$int$idx[,2]))
    locs = mesh$loc[idx,]
    xlim = range(locs[,1])
    ylim = range(locs[,2])
  } else {
    xlim = xRangeDat
    ylim = yRangeDat
  }
  # ux <- seq(xlim[1], xlim[2], , NP)
  # uy <- seq(ylim[1], ylim[2], , NP)
  # center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  maxRange = max(effectiveRangeVals)
  minRange = min(effectiveRangeVals)
  if(is.null(maxRadius))
    maxRadius = maxRange * 2
  minStep = minRange / 10
  ThisNP = 2 * maxRadius / minStep
  centerX = mean(xlim)
  widthX = xlim[2] - centerX
  deltaX = min(widthX, maxRadius)
  centerY = mean(ylim)
  widthY = ylim[2] - centerY
  deltaY = min(widthY, maxRadius)
  delta = min(deltaX, deltaY)
  xlim = c(centerX - delta, centerX + delta)
  # ux <- seq(xlim[1], xlim[2], , NP)
  ux <- c(seq(xlim[1], centerX-minRange/100, l=NP/2), centerX, seq(centerX+minRange/100, xlim[2], l=NP/2))
  ylim = c(centerY - delta, centerY + delta)
  # uy <- seq(ylim[1], ylim[2], , NP)
  uy <- c(seq(ylim[1], centerY-minRange/100, l=NP/2), centerY, seq(centerY+minRange/100, ylim[2], l=NP/2))
  center <- rbind(c(centerX, centerY))
  
  # precompute relevant matrices
  Acenter = inla.spde.make.A(mesh, center)
  Ax = inla.spde.make.A(mesh, cbind(ux, center[2]))
  Ay = inla.spde.make.A(mesh, cbind(center[1], uy))
  
  # make method for calculating individual covariance function
  getOneCovariance = function(parameters) {
    # get relevant parameters
    effectiveRange = parameters[1]
    rho = parameters[2]
    nuggetVar = parameters[3]
    i = parameters[4]
    print(paste0("Calculating covariances for iteration ", i, "/", length(rhoVals)))
    
    # calculate covariances
    x1 <- cbind(ux, rep(center[2], NP+1))
    x2 <- rbind(center)
    d <- c(rdist(x1, x2))
    Q = makeQSPDE(mesh, effectiveRange, rho)
    y = as.numeric(Acenter %*% inla.qsolve(Q, t(Ax)))
    y[NP/2+1] = y[NP/2+1] + nuggetVar
    x1 <- cbind(rep(center[1], NP+1), uy)
    d2 <- c(rdist(x1, x2))
    y2 = as.numeric(Acenter %*% inla.qsolve(Q, t(Ay)))
    y2[NP/2+1] = y2[NP/2+1] + nuggetVar
    
    # average x and y covariances
    # sortXI = sort(d, index.return=TRUE)$ix
    # d = d[sortXI]
    # y = y[sortXI]
    # sortYI = sort(d2, index.return=TRUE)$ix
    # d2 = d2[sortYI]
    # y2 = y2[sortYI]
    # d = rowMeans(cbind(d, d2))
    # y = rowMeans(cbind(y, y2))
    # spatialCov = y[NP/2+1]
    # return(cbind(d=d, cov=y, cor=y * (1 / max(y))))
    d = c(0, rowMeans(cbind(rev(d[1:(NP/2)]),  d[(NP/2+2):length(d)],  rev(d2[1:(NP/2)]),  d2[(NP/2+2):length(d2)])))
    y = c(mean(y[NP/2+1], y2[NP/2+1]), rowMeans(cbind(rev(y[1:(NP/2)]),  y[(NP/2+2):length(y)],  rev(y2[1:(NP/2)]),  y2[(NP/2+2):length(y2)])))
    sortXI = sort(d, index.return=TRUE)$ix
    d = d[sortXI]
    y = y[sortXI]
    return(cbind(d=d, cov=y, cor=y * (1 / max(y))))
  }
  
  # calculate covariances for each sample from the posterior
  parameterMat = cbind(effectiveRangeVals, rhoVals, nuggetVarVals, 1:length(rhoVals))
  # browser()
  out = apply(parameterMat, 1, getOneCovariance)
  d = out[1:(NP/2+1),1]
  covMat = out[(NP/2+2):(2*(NP/2+1)),]
  corMat = out[(2*(NP/2+1)+1):(3*(NP/2+1)),]
  
  # calculate summary statistics
  meanCov = rowMeans(covMat)
  lowerCov = apply(covMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCov = apply(covMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCor = rowMeans(corMat)
  lowerCor = apply(corMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCor = apply(corMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  
  # return results
  list(d=d, 
       cov=meanCov, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
       cor=meanCor, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
}

# test how close we can get to the spatial correlation function:
getTrueLKEffectiveRange = function(nLayer=3, NP=200, sigma2 = 0, rho=1, 
                                   nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=13, 
                                   latInfo=NULL, effectiveRange=1, alphas=rep(1/nLayer, nLayer-1)) {
  
  
  # construct the lattice
  if(is.null(latInfo)) {
    xRangeDat = c(-1, 1)
    yRangeDat = c(-1, 1)
    latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  }
  # set true parameter values
  kappa = (sqrt(8) * latInfo[[1]]$latWidth /effectiveRange)
  
  xlim <- latInfo[[1]]$xRangeDat
  ux <- seq(xlim[1], xlim[2], , NP)
  ylim <- latInfo[[1]]$yRangeDat
  uy <- seq(ylim[1], ylim[2], , NP)
  center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  test = getLKInlaCovarianceFun(kappa, rho, sigma2, alphas, latticeInfo = latInfo, 
                                normalize=normalize, fastNormalize=fastNormalize)
  ds = test[,1]
  firstI = match(TRUE, test[,3] <= .1)
  if(is.na(firstI)) {
    warning(paste0("effective range larger than domain diameter, returned value is too small.  Correlation is ", 
                   test[nrow(test),3], " at distance ", ds[length(ds)]))
    firstI = length(ds)
  }
  ds[firstI]
}

# calculate the minimum value of NC in order to have at least a certain number of basis functions total 
# for a fixed value of nLayer and nBuffer
getMinimalLattice = function(nMin=5000, locs=cbind(c(-1, 1), c(-1, 1)), nLayer=3, nBuffer=5) {
  # set starting NC so that the minimal dimension has at least two basis functions
  xRange = range(locs[,1])
  yRange = range(locs[,2])
  minRange = min(c(diff(xRange), diff(yRange)))
  maxRange = max(c(diff(xRange), diff(yRange)))
  NC = ceiling(2 * maxRange / minRange)
  
  # keep increasing NC until there are enough basis functions
  N = 0
  while(N < nMin) {
    latticeInfo = makeLatGrids(xRange, yRange, NC, nBuffer, nLayer)
    N = sum(sapply(latticeInfo, function(x) {x$ny * x$nx}))
    print(paste0(N, " basis functions for NC=", NC))
    NC = NC + 1
  }
  NC = NC - 1
  
  # return results
  c(NC=NC, nbasis=N)
}

rpcvar = function(n, alpha=.01, u=1) {
  1 / inla.pc.rprec(n, alpha=alpha, u=u)
}

dpcvar = function(x, alpha=.01, u=1) {
  inla.pc.dprec(1 / x, alpha=alpha, u=u) / x^2
}

qpcvar = function(p, alpha=.01, u=1) {
  1 / inla.pc.qprec(1-p, alpha=alpha, u=u)
}

ppcvar = function(q, alpha=.01, u=1, tol = 1e-10) {
  fun = function(x) {dpcvar(x, alpha=alpha, u=u)}
  integrate(fun, lower = tol, upper=q)$value
}

# generate simulations on rectangular domain (e.g. unit square) from a 
# Matern Gaussian Process with zero mean
genSimsMatern = function(xRange=c(0,1), yRange=c(0,1), n=100, nsim=100, 
                         beta=(xRange[2]-xRange[1])/10, nu=1.5, sigmaSq=1, tauSq=sigmaSq/10) {
  # generate locations of data
  xs = matrix(runif(n*nsim)*(xRange[2] - xRange[1]) + xRange[1], ncol=nsim)
  ys = matrix(runif(n*nsim)*(yRange[2] - yRange[1]) + yRange[1], ncol=nsim)
  
  # simulate from standard normal
  zSims = matrix(rnorm(n*nsim), ncol=nsim)
  errs = matrix(rnorm(n*nsim, sd=sqrt(tauSq)), ncol=nsim)
  
  # generate a simulation
  genSim = function(i) {
    # print progress
    if(i %% 10 == 0)
      print(paste0("iteration ", i, "/", nsim))
    
    # use Cholesky decomposition of covariance matrix to simulate
    Sigma = stationary.cov(cbind(xs[,i], ys[,i]), Covariance="MaternLR", theta=beta, phi=sigmaSq, nu=nu)
    L = t(chol(Sigma))
    L %*% zSims[,i]
  }
  trueMat = sapply(1:nsim, genSim)
  obsMat = trueMat + errs
  
  list(xs=xs, ys=ys, trueMat=trueMat, obsMat=obsMat)
}

# modified version of fields packages Matern function to code LR2011
# parameterization of the Matern covariance
# range = beta (this is the effective range for 99% correlation, roughly)
# all other parameters are the same as in ?Matern from fields package
MaternLR = function (d, range = 1, beta=range, alpha = 1/beta, smoothness = 0.5, nu = smoothness, 
                     phi = 1) {
  Matern(sqrt(8*nu)*d*alpha, range, alpha, smoothness, nu, phi)
}

# construct numerical integration matrix, constructing mx*my aggregation regions by 
# diving xRange and yRange into mx x my grid of regions
# predPts: Ideally, predPoints should be a regular grid of locations, since they are 
# all weighed equally when being aggregated in each region
makeNumericalIntegralMat = function(predPts, xRange=c(-1, 1), yRange=c(-1, 1), mx=3, my=3) {
  # construct aggregation matrix for predictions by testing which prediction locations 
  # are in which aggregation regions
  xRegionGrid = seq(xRange[1], xRange[2], l=mx + 1)[-1]
  yRegionGrid = seq(yRange[1], yRange[2], l=my + 1)[-1]
  xRegion = function(x) {
    match(TRUE, x <= xRegionGrid)
  }
  yRegion = function(y) {
    match(TRUE, y <= yRegionGrid)
  }
  xRegionI = sapply(predPts[,1], xRegion)
  yRegionI = sapply(predPts[,2], yRegion)
  regionI = (yRegionI-1)*mx + xRegionI
  getARow = function(ai) {regionI == ai}
  
  A = t(sapply(1:(mx*my), getARow))
  A = sweep(A, 1, rowSums(A), "/")
  
  A
}

# Divide the domain into four parts: everything to the left of what's missing, everything to the right, and 
# the two rectangles above and below what's missing. Randomly draw how many points come from each, and sample 
# from each independently
runifsqMissingRectangle = function(n, xRange=c(-1, 1), yRange=c(-1, 1), xRangeMissing=c(-1/3, 1/3), yRangeMissing=c(-1/3, 1/3), randomPointOrdering=TRUE) {
  # Calculate the areas of the four rectangles. Probability of point being drawn from each is proportional to these areas
  areaLeft = diff(yRange) * (xRangeMissing[1] - xRange[1])
  areaRight = diff(yRange) * (xRange[2] - xRangeMissing[2])
  areaAbove = diff(xRangeMissing) * (yRange[2] - yRangeMissing[2])
  areaBelow = diff(xRangeMissing) * (yRangeMissing[1] - yRange[1])
  totalArea = areaLeft + areaRight + areaAbove + areaBelow
  probabilityLeft = areaLeft / totalArea
  probabilityRight = areaRight / totalArea
  probabilityAbove = areaAbove / totalArea
  probabilityBelow = areaBelow / totalArea
  
  # determine how many points are in each of the rectangles
  out = sample(1:4, n, replace=TRUE, prob=c(probabilityLeft, probabilityRight, probabilityAbove, probabilityBelow))
  nLeft = sum(out == 1)
  nRight = sum(out == 2)
  nAbove = sum(out == 3)
  nBelow = sum(out == 4)
  
  # sample the points
  out = rbind(runifsq(nLeft, c(xRange[1], xRangeMissing[1]), yRange), 
              runifsq(nRight, c(xRangeMissing[2], xRange[2]), yRange), 
              runifsq(nAbove, xRangeMissing, c(yRangeMissing[2], yRange[2])),
              runifsq(nBelow, xRangeMissing, c(yRange[1], yRangeMissing[1]))
  )
  
  # scramble the points if necessary
  if(randomPointOrdering) {
    reordering = sample(1:n, n, replace=FALSE)
    out = out[reordering,]
  }
  
  out
}

# generate uniform observations on a rectangle
runifsq = function(n, xRange=c(-1, 1), yRange=c(-1, 1)) {
  xs = runif(n, xRange[1], xRange[2])
  ys = runif(n, yRange[1], yRange[2])
  cbind(xs, ys)
}

# for plotting areal data
# mapDat: has attributes "polygons" and "data", where data has names either NAME_1 or name_1 giving area names
# plotVar: the variable to plot in each area
# scaleFun, scaleFunInverse: by default, the identity function. Can be logit/expit or log/exp
# n.ticks, min.n: arguemnts for the pretty() function to make tick marks
# ...: arguments to polygon function
plotMapDat = function(mapDat, plotVar, varAreas, zlim=NULL, cols=tim.colors(), 
                      legend.mar=7, new=FALSE, plotArgs=NULL, main=NULL, xlim=NULL, xlab=NULL, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                      ylim=NULL, ylab=NULL, n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, asp=1, legend.width=1.2, addColorBar=TRUE, 
                      legendArgs=list(), leaveRoomForLegend=TRUE, ...) {
  
  # do setup for ploting data by county if necessary
  if(!is.null(plotVar)) {
    if(is.null(zlim)) {
      zlim = range(plotVar)
    }
    
    # get region names from map data
    if(!is.null(mapDat@data$NAME_1)) {
      regionNames = mapDat@data$NAME_1
    } else if(!is.null(mapDat@data$name_1)) {
      regionNames = as.character(mapDat@data$name_1)
    } else {
      stop("mapDat has unrecognized region names")
    }
  }
  
  # generate new plot if necessary
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    if(is.null(main))
      main = ""
    
    if(is.null(plotArgs)) {
      plotArgs = list(main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, asp=asp)
    } else {
      plotArgs$main = main
      plotArgs$xlab = xlab
      plotArgs$ylab = ylab
      plotArgs$xlim = xlim
      plotArgs$ylim = ylim
      plotArgs$asp = asp
    }
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    do.call("plot", c(list(1, 2, type="n"), plotArgs))
  }
  
  # add polygons to plot
  polys = mapDat@polygons
  plotCounty = function(i) {
    countyPolys = polys[[i]]@Polygons
    
    if(is.null(plotVar)) {
      sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(countyPolys[[x]]@coords), list(...)))})
    }
    else {
      # get index of plotVar corresponding to this county
      thisI = which(varAreas == regionNames[i])
      
      # get color to plot
      vals = c(zlim, scaleFun(plotVar[thisI]))
      vals = vals-vals[1]
      vals = vals/(vals[2] - vals[1])
      col = cols[round(vals[3]*(length(cols)-1))+1]
      
      sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(countyPolys[[x]]@coords, col=col), list(...)))})
    }
    
  }
  sapply(1:length(polys), plotCounty)
  
  if(!is.null(plotVar) && addColorBar) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(ticks))
      ticks = scaleFun(pretty(scaleFunInverse(zlim), n=n.ticks, min.n=min.n))
    else
      ticks = scaleFun(ticks)
    if(is.null(tickLabels))
      tickLabels = scaleFunInverse(ticks)
    
    # par( oma=c( 0,0,0,3))
    
    # set list of arguments to image.plot
    legendArgs$zlim=zlim
    legendArgs$nlevel=length(cols)
    legendArgs$legend.only=TRUE
    legendArgs$horizontal=FALSE
    legendArgs$col=cols
    legendArgs$add = TRUE
    if(is.null(legendArgs$axis.args))
      legendArgs$axis.args=list(at=ticks, labels=tickLabels)
    else {
      legendArgs$axis.args$at=ticks
      legendArgs$axis.args$labels=tickLabels
    }
    legendArgs$legend.mar=legend.mar
    legendArgs$legend.width=legend.width
    do.call("image.plot", legendArgs)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}


makeRedBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3)
  else
    scale_colour_continuous_sequential(h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3, n_interp=n)
}

makeGreenBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
  else
    scale_colour_continuous_sequential(h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2, n_interp=n)
}

makePurpleYellowSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9)
  else
    scale_colour_continuous_sequential(h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9, n_interp=n)
}

makeRedBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev)
    else
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n, mid=center)
    }
  }
}

makeRedGrayBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    if(!ggplot)
      scale_colour_continuous_diverging(n_interp=n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedGrayBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(n_interp, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev, mid=center)
    }
  }
}

makeBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE, n_interp=n)
}

makeBlueYellowSequentialColors = function(n, ggplot=FALSE, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, rev=rev)
  else
    scale_colour_continuous_sequential(h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, n_interp=n, rev=rev)
}

makeYellowRedSequentialColors = function(n, ggplot=FALSE, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=15, h2=79, c1=100, c2=52, l1=55, l2=95, p1=1.2, rev=rev)
  else
    scale_colour_continuous_sequential(h1=15, h2=79, c1=100, c2=52, l1=55, l2=95, p1=1.2, n_interp=n, rev=rev)
}

makeRedGreenDivergingColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
  else
    scale_colour_continuous_sequential(h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5, n_interp=n)
}

# add in binomial variation to the probability sampling matrix
addBinomialVar = function(probMatrix, ns) {
  simulatedObservations = matrix(rbinom(n=length(probMatrix), size=rep(ns, ncol(probMatrix)), prob=c(as.matrix(probMatrix))), nrow=nrow(probMatrix))
  sweep(simulatedObservations, 1, 1/ns, "*")
}
