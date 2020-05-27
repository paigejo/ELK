##### This script is supposed to contain ONLY the functions and packages required by the 
##### latticeKrig rgeneric functions

# require the following packages
require(Matrix)
require(spam)
require(fields)
require(invgamma)
require(MCMCpack) # for ddirichlet

## code for calculating precision and basis matrices for coding LatticeKrig into INLA:

# profvis({
# library(fields)

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
makeQPrecomputed = function(precomputedMatrices, kappa=1, rho=1, latticeInfo, alphas=NULL, nu=NULL, 
                            normalized=FALSE, fastNormalize=FALSE, returnctildes=FALSE, ctildes=NULL) {
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
  if(nLayer == 1 && is.null(alphas))
    alphas = 1
  
  mBxymBxyT = precomputedMatrices$mBxymBxyT # - Bxy - t(Bxy)
  BxyTBxy = precomputedMatrices$BxyTBxy # t(Bxy) %*% Bxy
  ms = precomputedMatrices$ms
  
  # # make B, SAR regression matrix for Bc = e
  # B = Diagonal(n=nrow(Bxy), x=kappa^2) - Bxy
  # 
  # # compute precision matrix
  # Q = (1/rho) * t(B) %*% B
  if(length(kappa) == 1)
    Q = (1/rho) * ( Diagonal(n=nrow(mBxymBxyT), x=kappa^4) + kappa^2 * mBxymBxyT + BxyTBxy)
  else {
    # for multiple values of kappa, take advantage of block diagonal structure of Q
    nbasis = sapply(latticeInfo, function(x) {x$nx * x$ny})
    Q = (1/rho) * ( Diagonal(n=nrow(mBxymBxyT), x=rep(kappa^4, nbasis)) + 
                    Diagonal(n=nrow(mBxymBxyT), x=rep(kappa^2, nbasis)) %*% mBxymBxyT + 
                    BxyTBxy)
  }
  
  if(normalized) {
    
    # now get relevant row of A for value at midpoint of domain
    A = precomputedMatrices$A
    At = precomputedMatrices$At
    
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
      # ctilde = as.numeric(Ai %*% inla.qsolve(Q, t(Ai)))/rho
      # Qtilde = ctilde * Q
      # Q = Qtilde # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE, newnormalize=TRUE)): 2.72
      if(is.null(ctildes))
        ctildes = as.numeric(diag(A %*% inla.qsolve(Q, At)))/(rho * alphas)
      Qtilde = Diagonal(n=nrow(Q), x=rep(ctildes, ms)) %*% Q
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
      if(is.null(ctildes))
        ctildes = as.numeric(diag(A %*% inla.qsolve(Qnorm, t(A))))/(rho * alphas)
      Qtilde = Diagonal(n=nrow(Qnorm), x=rep(ctildes, ms)) %*% Qnorm
      Q = Qtilde
    }
    # # compute how good an approximation it was
    # hist(diag(solve(Q)))
    # hist(diag(solve(Q2)))
    # hist(diag(solve(Qtilde)))
    # image(Q)
    # image(Q2)
    # image(Qtilde)
  }
  
  if(!returnctildes)
    Q
  else
    ctildes
}

# computing basis matrix with Wendland covariance
# theta: 
# latticeInfo: an object returned by makeLatGrids
makeA = function(predPts=NULL, latticeInfo, thisLayer=1, theta=NULL, maxLayer=length(latticeInfo)) {
  require(Matrix)
  require(spam)
  require(fields)
  nLayer = length(latticeInfo)
  xRangeKnot = latticeInfo[[thisLayer]]$xRangeKnots
  yRangeKnot = latticeInfo[[thisLayer]]$yRangeKnots
  xNKnot = latticeInfo[[thisLayer]]$nx
  yNKnot = latticeInfo[[thisLayer]]$ny
  xRangeDat = latticeInfo[[thisLayer]]$xRangeDat
  yRangeDat = latticeInfo[[thisLayer]]$yRangeDat
  nBuffer = latticeInfo[[thisLayer]]$nBuffer
  NC = latticeInfo[[thisLayer]]$NC
  knotPts = latticeInfo[[thisLayer]]$latCoords
  
  # setup theta if null
  # if(is.null(theta))
    theta = 2.5*(xRangeKnot[2]-xRangeKnot[1])/(xNKnot-1)
  
  # set up prediction locations if necessary
  # have some sort of default behavior for setting predPts for ease of testing
  if(is.null(predPts)) {
    xRangePred = xRangeKnot
    yRangePred = yRangeKnot
    xNPred = xNKnot*3
    yNPred = yNKnot*3
    
    predPts = make.surface.grid(list(x=seq(xRangePred[1], xRangePred[2], l=xNPred), 
                                     y=seq(yRangePred[1], yRangePred[2], l=yNPred)))
  }
  
  # adjust basis elements to depend on layer
  origXNKnot = xNKnot
  origYNKnot = yNKnot
  origTheta = theta
  # xNKnot = xNKnot * 2^(thisLayer-1)
  # yNKnot = yNKnot * 2^(thisLayer-1)
  # theta = theta / 2^(thisLayer-1)
  
  # # generate knot lattice locations and filter out locations 
  # # too far outside of the data domain
  # knotXs = seq(xRangeKnot[1], xRangeKnot[2], l=xNKnot)
  # knotYs = seq(yRangeKnot[1], yRangeKnot[2], l=yNKnot)
  # if(sum(knotXs > xRangeDat[2]) > nBuffer)
  #   knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
  # if(sum(knotXs < xRangeDat[1]) > nBuffer)
  #   knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
  # if(sum(knotYs > yRangeDat[2]) > nBuffer)
  #   knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
  # if(sum(knotYs < yRangeDat[1]) > nBuffer)
  #   knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
  # knotPts = make.surface.grid(list(x=knotXs, y=knotYs))
  
  # NOTE: the wendland.cov function in `fields' does fancy things for only computing 
  # covariances for distances smaller than theta.  Returns a spam matrix
  thisA = as.dgCMatrix.spam(wendland.cov(predPts, knotPts, theta=theta, k=2))
  
  # return results recursively
  if(thisLayer == min(nLayer, maxLayer)) {
    return(thisA)
  }
  else {
    return(cbind(thisA, makeA(predPts, latticeInfo, thisLayer+1, theta=origTheta)))
  }
}
# out = makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE, newnormalize= TRUE)})

# makes alpha vector based on theoretical relationship in Nychka 2015 paper Cor 1 in 4.1
# NOTE: used to be based on ?LKrig, but it doesn't match the paper or 
#       what the package actually does (in LKrigSetupAlpha.default)
getAlphas = function(nLayer=3, nu=.5) {
  # 
  # alphas = exp(-2*(1:nLayer) * nu)
  # alphas = alphas/sum(alphas)
  # alphas
  if(nu <= exp(-100)) {
    alphas = rep(1 / nLayer, nLayer)
  } else if(nu>300) {
    alphas = c(1, rep(0, nLayer - 1))
  } else {
    thetaL=  2^(-(1:nLayer))
    alphas = thetaL^(2*nu)
    alphas = alphas/sum(alphas)
  }
  alphas
}

# make the graph of the SAR model for LatticeKrig
makeGraph = function(latticeInfo, thisLayer=1) {
  require(Matrix)
  require(spam)
  require(fields)
  nx = latticeInfo[[thisLayer]]$nx
  ny = latticeInfo[[thisLayer]]$ny
  
  # make (Laplacian) differential operators
  Dx = bandSparse(nx, k=0:1, diag=list(rep(1, nx), rep(1, nx-1)), symmetric=TRUE)
  Dy = bandSparse(ny, k=0:1, diag=list(rep(1, ny), rep(1, ny-1)), symmetric=TRUE)
  Dxx = bandSparse(nx, k=0:2, diag=list(rep(1, nx), rep(1, nx-1), rep(1, nx-2)), symmetric=TRUE)
  Dyy = bandSparse(ny, k=0:2, diag=list(rep(1, ny), rep(1, ny-1), rep(1, ny-2)), symmetric=TRUE)
  
  # generate x and y (Laplacian) differential operators
  Ix = Diagonal(n=nx)
  Iy = Diagonal(n=ny)
  Bx = kronecker(Iy, Dxx)
  By = kronecker(Dyy, Ix)
  Bxy = kronecker(Dy, Dx)
  
  # based on B, SAR regression matrix for Bc = e
  # G = (Bx + By + Bxy) > 0
  G = Bx + By + Bxy
  
  # if only need one layer, return what we have.  Otherwise, make sparse 
  # block diagonal matrix
  nLayer = length(latticeInfo)
  if(nLayer == thisLayer)
    return(G)
  else if((nLayer != 1) && (thisLayer == 1))
    return(bdiag(c(list(G), makeGraph(latticeInfo, thisLayer+1))))
  else
    return(c(list(G), makeGraph(latticeInfo, thisLayer+1)))
}

# based on:
# https://www.stat.washington.edu/~rje42/lca/html/logitTrans.html
# Used for optimization of values something to 1
multivariateLogit = function(x) {
  denominator = 1 - sum(x)
  log(x / denominator)
}
multivariateExpit = function(z) {
  denominator = 1 + sum(exp(z))
  out = exp(z) / denominator
  out
}

# this function calculates the density of a real valued vector, z, whose multivariateExpit is 
# Dirichlet with parameters alpha (alpha is the vector of shape parameters in ddirichlet). 
# Note that if length(z) == n, then length(alpha) = n+1, with the (n+1)st value corresponding to 
# 1 - sum(multivariateLogit(z)).
dexpitDirichlet = function(z, alpha, doLog=FALSE) {
  # first transform from z to probabilities summing to 1
  x = multivariateExpit(z)
  x = c(x, 1 - sum(x))
  
  # get the raw density (without the jacobian factor)
  rawDensity = ddirichlet(x, alpha, doLog=doLog)
  
  # get the jacobian matrix and its (log) determinant
  J = multivariateExpitJacobian(z)
  d = determinant(J, logarithm = TRUE)$modulus
  ifelse(doLog, rawDensity + d, rawDensity * exp(d))
}
# this function calculates the jacobian matrix of the multivariateExpit transformation
multivariateExpitJacobian = function(z) {
  # based on results from WolframAlpha
  normalizer = (1 + sum(exp(z)))^2
  mat = outer(z, z, FUN = function(zi, zj) {-exp(zi + zj)})
  diag(mat) = diag(mat) + exp(z)*(1+sum(exp(z)))
  mat * (1 / normalizer)
}
# Returns the belief that the probabilities of K rival events are x_i given that 
# each event has been observed alpha_i - 1 times.
ddirichlet = function(x, alpha, doLog=FALSE) {
  # ddirichlet from MCMCpack:
  # logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  # s <- sum((alpha - 1) * log(x))
  # exp(sum(s) - logD)
  
  # based on dirichlet.log.pdf from mmppr:
  if (length(x) == 1) {
    return(1)
  }
  logResult = sum((alpha - 1) * log(x + 1e-08)) - sum(lgamma(alpha)) + 
    lgamma(sum(alpha))
  
  ifelse(doLog, logResult, exp(logResult))
}