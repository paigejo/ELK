## create the rgeneric model for inla

## function for making prior on correlation scale and marginal sd (a PC prior)
# corScaleMed: median correlation scale
# sdAlpha, sdU: set PC prior so that P(sd > sdU) = sdAlpha
getPCPrior = function(corScaleMed, sdAlpha, sdU, dirichletConcentration=1.5, nLayer=3, separateRanges=FALSE, latticeInfo=NULL) {
  ## set the correlation scale hyperparameter
  if(!separateRanges) {
    corScalePar = corScaleMed*log(2)
  } else {
    # set the priors so that the coarsest layer scale parameter starts at corScaleMed*log(2), and 
    # the other scale parameters decay depending on the lattice widths:
    if(!is.null(latticeInfo)) {
      deltas = sapply(latticeInfo, function(x) {x$latWidth})
      corScalePar = corScaleMed*log(2) * (deltas / max(deltas))
    } else {
      deltas = 2^seq(0, -nLayer+1, by=-1)
      corScalePar = corScaleMed*log(2) * (deltas / max(deltas))
    }
    
  }
  
  ## set the layer weight hyperparameters
  alphaPar = rep(dirichletConcentration / nLayer, nLayer)
  
  list(corScalePar=corScalePar, alpha=sdAlpha, u=sdU, alphaPar=alphaPar, priorType="pc")
}

# generate lattice points for all layers
# NOTE: if any of the possibly NULL input variables are NULL, they are all set to the default
makeLatGrids = function(xRangeDat=c(0,1), yRangeDat=c(0,1), NC=5, nBuffer=5, nLayer=1) {
  
  grids = list()
  for(thisLayer in 1:nLayer) {
    i = thisLayer
    
    # set layer parameters
    if(length(NC) == 1)
      layerNC = (NC-1)*2^(i-1) + 1
    else
      layerNC = NC[i]
    
    if(length(nBuffer) == 1)
      thisBuffer = nBuffer
    else
      thisBuffer = nBuffer[i]
    
    grids = c(grids, list(makeLatGrid(xRangeDat, yRangeDat, layerNC, nBuffer)))
  }
  
  grids
}

# generate lattice points for all layers for Kenya datasets
# NOTE: if any of the possibly NULL input variables are NULL, they are all set to the default
makeLatGridsKenya = function(nLayer=1, NC=15, nBuffer=5) {
  out = load(paste0("dataPointsKenya.RData"))
  xRange = dataPointsKenya$xRange
  yRange = dataPointsKenya$yRange
  makeLatGrids(xRange, yRange, NC, nBuffer, nLayer)
}

# converts from LatticeKrig lattice parameters to LKinla's
# xRange: range of x-coords of data
# yRange: range of y-coords of data
# NC: number of basis functions to put on largest dimension of data
# nBuffer: number of lattice points outside of range of data to include as buffer
makeLatGrid = function(xRange, yRange, NC=5, nBuffer=5) {
  xLength = xRange[2] - xRange[1]
  yLength = yRange[2] - yRange[1]
  if(xLength >= yLength) {
    # lattice point xs are on edge of data domain, lattice point ys aren't
    latWidth = xLength/(NC-1)
    Buff = latWidth*nBuffer
    xRangeKnots = c(xRange[1] - Buff, xRange[2] + Buff)
    nx = NC + 2*nBuffer
    ny = floor(yLength/latWidth) + 2*nBuffer + 1
    yMidPt = mean(yRange)
    extraYLength = yLength - floor(yLength/latWidth)*latWidth
    yRangeKnots = c(yRange[1] + 0.5*extraYLength - Buff, yRange[2] - 0.5*extraYLength + Buff)
  }
  else {
    # lattice point ys are on edge of data domain, lattice point xs aren't
    latWidth = yLength/(NC-1)
    Buff = latWidth*nBuffer
    yRangeKnots = c(yRange[1] - Buff, yRange[2] + Buff)
    ny = NC + 2*nBuffer
    nx = floor(xLength/latWidth) + 2*nBuffer + 1
    xMidPt = mean(xRange)
    extraXLength = xLength - floor(xLength/latWidth)*latWidth
    xRangeKnots = c(xRange[1] + 0.5*extraXLength - Buff, xRange[2] - 0.5*extraXLength + Buff)
  }
  xGrid = seq(xRangeKnots[1], xRangeKnots[2], l=nx)
  yGrid = seq(yRangeKnots[1], yRangeKnots[2], l=ny)
  latCoords = make.surface.grid(list(x=xGrid, y=yGrid))
  
  list(xRangeKnots=xRangeKnots, nx=nx, yRangeKnots=yRangeKnots, ny=ny, 
       xRangeDat=xRange, yRangeDat=yRange, NC=NC, nBuffer=nBuffer, latWidth=latWidth, 
       xGrid=xGrid, yGrid=yGrid, latCoords=latCoords)
}

# compute number of basis elements per layer.  Returns a list with element at index i 
# equal to the number of basis elements in layer i.
getMs = function(xRangeDat=c(0,1), yRangeDat=c(0,1), NC=5, nBuffer=5, nLayer=1) {
  # makeLatGrids(xRangeDat=xRangeDat, yRangeDat=yRangeDat, NC=NC, nBuffer=nBuffer, nLayer=nLayer)$ms
  sapply(makeLatGrids(xRangeDat=xRangeDat, yRangeDat=yRangeDat, NC=NC, nBuffer=nBuffer, nLayer=nLayer), function(x) {x$nx*x$ny})
}

# compute number of basis elements per layer.  Returns a list with element at index i 
# equal to the number of basis elements in layer i.
getMsDomain = function(latInfo) {
  # makeLatGrids(xRangeDat=xRangeDat, yRangeDat=yRangeDat, NC=NC, nBuffer=nBuffer, nLayer=nLayer)$ms
  sapply(latInfo, function(x) {(x$nx-x$nBuffer)*(x$ny-x$nBuffer)})
}

# convert from explicit grid parameters to LatticeKrig grid parameters
rawGridToLK = function(xRangeKnot=c(0,1), xNKnot=10, yRangeKnot=c(0,1), yNKnot=10, nBuffer=5) {
  # first compute NC, the number of lattice points in the largest dimension
  xNKnotInRange = xNKnot - 2*nBuffer
  yNKnotInRange = yNKnot - 2*nBuffer
  NC = max(xNKnotInRange, yNKnotInRange)
  
  # subtract off the buffer to get the "range of the data" under the LatticeKrig parameterization
  latWidth = diff(xRangeKnot)/(xNKnot-1)
  xGrid = seq(xRangeKnot[1], xRangeKnot[2], l=xNKnot)
  yGrid = seq(yRangeKnot[1], yRangeKnot[2], l=yNKnot)
  xRangeDat = c(xGrid[nBuffer+1], xGrid[xNKnot - nBuffer])
  yRangeDat = c(yGrid[nBuffer+1], yGrid[yNKnot - nBuffer])
  list(NC=NC, nBuffer=nBuffer, xRangeDat=xRangeDat, yRangeDat=yRangeDat)
}

inla.rgeneric.lk.model.standard = function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL)
{
  # print(paste0("Starting rgeneric call with command ", cmd))
  startTime = proc.time()[3]
  envir = environment(sys.call()[[1]])
  
  # theta is of the form:
  # c(betas, effectiveCor, sigmaSq, kappa, rho, nu, alphas)
  interpret.theta = function()
  {
    ## internal helper-function to map the parameters from the internal-scale to the
    ## user-scale
    if(!is.null(theta)) {
      # get effective correlation range and marginal variance
      effectiveCor = exp(theta[1])
      sigmaSq = exp(theta[2])
      
      # compute layer weights, alpha_1, ..., alpha_L
      L = nLayer = length(latInfo)
      
      if(L != 1) {
        if(!exists("multivariateExpit")) {
          # load relevant external functions
          if(printVerboseTimings)
            print("sourcing LKinla.R...")
          inf = sessionInfo()
          if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
            source("~/git/LK-INLA/LKinla.R")
          else
            source("U:/git/LK-INLA/LKinla.R")
        }
        # alphas = getAlphas(L, nu)
        alphas = multivariateExpit(theta[3:(2 + L - 1)])
        alphas = c(alphas, 1 - sum(alphas))
      }
      else {
        alphas = NULL
      }
    }
    else {
      effectiveCor = NULL
      sigmaSq = NULL
      alphas = NULL
    }
    
    # precomputations: get lattice grid cell width, convert parameters from effective correlation 
    # and marginal variance to kappa and rho.  Use spline to convert from marginal variance to kappa
    latticeWidth = latInfo[[1]]$latWidth
    kap = sqrt(8)/effectiveCor * latticeWidth
    
    # since we are normalizing the process, rho is just sigmaSq
    rho = sigmaSq
    
    list(effectiveCor = effectiveCor, 
         sigmaSq = sigmaSq, 
         kappa = kap, 
         rho = rho, 
         alphas = alphas)
  }
  
  if(cmd != "initial") {
    pars = interpret.theta()
    effectiveCor = pars$effectiveCor
    sigmaSq = pars$sigmaSq
    kap = pars$kappa
    rho = pars$rho
    alphas = pars$alphas
  }
  
  # returns matrix just like Q except ony zero and 1 for nonzero elements
  graph = function(){
    if(!exists("makeGraph")) {
      # load relevant external functions
      if(printVerboseTimings)
        print("sourcing LKinla.R...")
      inf = sessionInfo()
      if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
        source("~/git/LK-INLA/LKinla.R")
      else
        source("U:/git/LK-INLA/LKinla.R")
    }
    
    makeGraph(latInfo)
  }
  
  # compute the precision matrix
  Q = function() {
    if(!exists("makeQPrecomputed")) {
      # load relevant external functions
      if(printVerboseTimings)
        print("sourcing LKinla.R...")
      inf = sessionInfo()
      if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
        source("~/git/LK-INLA/LKinla.R")
      else
        source("U:/git/LK-INLA/LKinla.R")
    }
    
    ctildes = precomputedNormalizationFun$fullFun(effectiveCor, alphas)
    makeQPrecomputed(precomputedMatrices=precomputedMatrices, kap, rho, latInfo, alphas=alphas, normalized=normalize, 
                     fastNormalize=fastNormalize, ctildes=ctildes)
  }
  
  # get mean of each latent coefficient
  mu = function() {
    # 
    # # get LatticeKrig grid parameters, number of knots in each layer
    # ms = getMs(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
    # 
    # # total number of basis functions should be this number
    # # return(rep(0, nx*ny*(4^nLayer - 1)/3))
    # rep(0, sum(ms))
    
    return(numeric(0))
  }
  
  log.norm.const = function()
  {
    ## let INLA compute it as -n/2 log(2pi) + 1/2 * log(|Q|)
    return (numeric(0))
  }
  
  log.prior = function() {
    require(invgamma)
    
    # get prior (note the Jacobian factors)
    # NOTE: unlike for the simple model, here "orig" versus "pc" refers to Dirichlet vs PC priors for alpha.
    # PC priors are always used for the marginal variance
    # corScalePar is the median of the prior belief of the spatial correlation scale parameter
    if(prior$priorType == "orig") {
      # rho = exp(theta)
      # P(theta < t) = P(rho < exp(t))
      # p_theta(t) = p_rho(exp(t)) |exp(t)|
      # log p_theta(theta) = log p_rho(rho) + log(rho)
      # 
      # sigmasq = exp(theta)
      # P(theta < t) = P(sigma < exp(0.5 t))
      # p_theta(t) = p_sigma(exp(t)) * |exp(t)|
      # log p_theta(theta) = log p_sigma(sigmasq) + log(sigmasq)
      out = dinvexp(effectiveCor, rate=prior$corScalePar, log=TRUE) + log(effectiveCor) + 
        invgamma::dinvgamma(sigmaSq, shape=prior$varPar1, rate=prior$varPar2, log=TRUE) + log(sigmaSq)
    } else if(prior$priorType == "pc") {
      # rho = exp(theta)
      # P(theta < t) = P(rho < exp(t))
      # p_theta(t) = p_rho(exp(t)) |exp(t)|
      # log p_theta(theta) = log p_rho(rho) + log(rho)
      # 
      # tau = exp(-theta) (here tau is 1/sigma^2, the precision)
      # P(theta < t) = P(tau < exp(-t))
      # p_theta(t) = p_tau(exp(-t)) * |-exp(-t)|
      # log p_theta(theta) = log p_tau(tau) + log(tau)
      out = dinvexp(effectiveCor, rate=prior$corScalePar, log=TRUE) + log(effectiveCor) + 
        inla.pc.dprec(1/sigmaSq, u=prior$u, alpha=prior$alpha, log=TRUE) - log(sigmaSq)
    }
    if(length(latInfo) == 1) {
      out
    } else {
      if(!exists("dexpitDirichlet")) {
        # load relevant external functions
        if(printVerboseTimings)
          print("sourcing LKinla.R...")
        inf = sessionInfo()
        if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
          source("~/git/LK-INLA/LKinla.R")
        else
          source("U:/git/LK-INLA/LKinla.R")
      }
      
      # the jacobian factor has already been included in dexpitDirichlet
      # Add on the density of a random vector whose multivariateExpit is Dirichlet
      out + dexpitDirichlet(multivariateLogit(alphas[1:(length(alphas)-1)]), prior$alphaPar, doLog=TRUE)
    }
  }
  
  initial = function() {
    # # initialize fixed effects and process variances
    # mod = lm(ys ~ X-1)
    # betas = coef(mod)
    # 
    # # remove from observations
    # ysCntr = ys - X %*% betas
    # 
    # # initialize covariance parameters
    # effectiveRangeInit = ((xRange[2] - xRange[1]) + (yRange[2] - yRange[1]))/2
    # xLength = xRange[2] - xRange[1]
    # squareWidth = xLength/(nx-1)
    # sigmaSq = sum((ysCntr - 0)^2)/n
    
    # get range of data
    xRangeDat = latInfo[[1]]$xRangeDat
    yRangeDat = latInfo[[1]]$yRangeDat
    
    # initialize process variance by estimating spatial process with OLS
    if(!exists("makeA")) {
      # load relevant external functions
      if(printVerboseTimings)
        print("sourcing LKinla.R...")
      inf = sessionInfo()
      if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
        source("~/git/LK-INLA/LKinla.R")
      else
        source("U:/git/LK-INLA/LKinla.R")
    }
    AMat = makeA(datCoords, latInfo)
    mod = lm(ys ~ cbind(X, as.matrix(AMat)) - 1)
    r2 = summary(mod)$r.squared
    s2 = summary(mod)$sigma^2
    # sigmaSq = (r2/(1-r2)) * s2
    sigmaSq = var(ys) * r2
    
    # initialize covariance parameters
    if(is.null(initialEffectiveRange)) {
      # We want the "middle" representable effective range to be an eighth of the domain diameter
      middleEffectiveRange = ((xRangeDat[2] - xRangeDat[1]) + (yRangeDat[2] - yRangeDat[1]))/8
      nLayer = length(latInfo)
      exponent = (nLayer-1)/2
      effectiveRangeInit = middleEffectiveRange * 2^exponent
    }
    else
      effectiveRangeInit = initialEffectiveRange
    
    # initialize the layer weights to be equal
    if(length(latInfo) == 1)
      c(log(effectiveRangeInit), log(sigmaSq))
    else {
      if(is.null(initialAlphas))
        initAlphas = rep(1/length(latInfo), length(latInfo)-1)
      else
        initAlphas = initialAlphas
      c(log(effectiveRangeInit), log(sigmaSq), multivariateLogit(initAlphas))
    }
  }
  
  quit = function() {
    return(invisible())
  }
  
  if(is.null(theta)) { theta = initial() }
  val = do.call(match.arg(cmd), args = list())
  
  totalTime = proc.time()[3] - startTime
  if(totalTime >= .1 && printVerboseTimings)
    print(paste0("rgeneric call with command ", cmd, " had the significant total time: ", totalTime))
  return (val)
}

inla.rgeneric.lk.model.full = function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL)
{
  # print(paste0("Starting rgeneric call with command ", cmd))
  startTime = proc.time()[3]
  envir = environment(sys.call()[[1]])
  
  # # convert to theoretical kappa estimate
  # kappaEst = sqrt(8)/effectiveRangeInit * latWidth
  # 
  # # now to the important part: get an initial estimate of marginal variance and a range (go with /100 and *100 of initial guess of kappa)
  # kappas <<- 10^(seq(log10(kappaEst)-2, log10(kappaEst)+2, l=100))
  # kappaWidths <<- log10(kappas[2]) - log10(kappas[1])
  # margVars <<- sapply(kappas, getMultiMargVar, rho=1, tod=2.5, nx=nx, ny=ny, nu=nu, nLayer=nLayer)[1,]
  # logMargVarSpline <<- splinefun(log(kappas), margVars, "natural")
  
  # theta is of the form:
  # c(betas, effectiveCor, sigmaSq, kappa, rho, nu, alphas)
  interpret.theta = function()
  {
    ## internal helper-function to map the parameters from the internal-scale to the
    ## user-scale
    if(!is.null(theta)) {
      # compute layer weights, alpha_1, ..., alpha_L
      L = nLayer = length(latInfo)
      
      # get effective correlation range and marginal variance
      effectiveCor = exp(theta[1:(0+L)])
      sigmaSq = exp(theta[1+L])
      
      if(L != 1) {
        if(!exists("multivariateExpit")) {
          # load relevant external functions
          if(printVerboseTimings)
            print("sourcing LKinla.R...")
          inf = sessionInfo()
          if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
            source("~/git/LK-INLA/LKinla.R")
          else
            source("U:/git/LK-INLA/LKinla.R")
        }
        # alphas = getAlphas(L, nu)
        alphas = multivariateExpit(theta[(2 + L):(2 + 2*L - 2)])
        alphas = c(alphas, 1 - sum(alphas))
      }
      else {
        alphas = NULL
      }
    }
    else {
      effectiveCor = NULL
      sigmaSq = NULL
      alphas = NULL
    }
    
    # precomputations: get lattice grid cell width, convert parameters from effective correlation 
    # and marginal variance to kappa and rho.  Use spline to convert from marginal variance to kappa
    latticeWidth = sapply(latInfo, function(x) {x$latWidth})
    kap = sqrt(8)/effectiveCor * latticeWidth
    
    # since we are normalizing the process, rho is just sigmaSq
    rho = sigmaSq
    
    # # kap = sqrt(8)/effectiveCor * latticeWidth
    # 
    # # If we're at a value of kappa outside out current reparameterization range, 
    # # adjust the range by refitting the spline function
    # if(kap > max(kappas)) {
    #   newKappas = 10^(seq(log10(kappas[length(kappas)] + kappaWidths), log10(kap)+.5, by=kappaWidths))
    #   kappas <<- c(kappas, newKappas)
    #   kappas = c(kappas, newKappas)
    #   margVars <<- c(margVars, sapply(newKappas, getMultiMargVar, rho=1, tod=2.5, nx=nx, ny=ny, nu=nu, nLayer=nLayer)[1,])
    #   logMargVarSpline <<- splinefun(log(kappas), margVars, "natural")
    # }
    # else if(kap < min(kappas)) {
    #   newKappas = 10^(seq(log10(kap) - .5, log10(kappas[1] - kappaWidths), by=kappaWidths))
    #   kappas <<- c(newKappas, kappas)
    #   margVars <<- c(sapply(newKappas, getMultiMargVar, rho=1, tod=2.5, nx=nx, ny=ny, nu=nu, nLayer=nLayer)[1,], margVars)
    #   logMargVarSpline <<- splinefun(log(kappas), margVars, "natural")
    # }
    
    # # now find rho using reparameterization:
    # # h(log kappa) = sigma^2/rho
    # # rho = sigma^2/h(log kappa)
    # 
    # # rho = sigmaSq * 4*pi * kappa^2
    # rho = sigmaSq / logMargVarSpline(log(kap))
    
    list(effectiveCor = effectiveCor, 
         sigmaSq = sigmaSq, 
         kappa = kap, 
         rho = rho, 
         alphas = alphas)
  }
  
  if(cmd != "initial") {
    pars = interpret.theta()
    effectiveCor = pars$effectiveCor
    sigmaSq = pars$sigmaSq
    kap = pars$kappa
    rho = pars$rho
    alphas = pars$alphas
  }
  
  # returns matrix just like Q except ony zero and 1 for nonzero elements
  graph = function(){
    if(!exists("makeGraph")) {
      # load relevant external functions
      if(printVerboseTimings)
        print("sourcing LKinla.R...")
      inf = sessionInfo()
      if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
        source("~/git/LK-INLA/LKinla.R")
      else
        source("U:/git/LK-INLA/LKinla.R")
    }
    
    makeGraph(latInfo)
  }
  
  # compute the precision matrix
  Q = function() {
    if(!exists("makeQPrecomputed")) {
      # load relevant external functions
      if(printVerboseTimings)
        print("sourcing LKinla.R...")
      inf = sessionInfo()
      if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
        source("~/git/LK-INLA/LKinla.R")
      else
        source("U:/git/LK-INLA/LKinla.R")
    }
    
    ctildes = precomputedNormalizationFun$fullFun(effectiveCor, alphas)
    makeQPrecomputed(precomputedMatrices=precomputedMatrices, kap, rho, latInfo, alphas=alphas, normalized=normalize, 
                     fastNormalize=fastNormalize, ctildes=ctildes)
  }
  
  # get mean of each latent coefficient
  mu = function() {
    # 
    # # get LatticeKrig grid parameters, number of knots in each layer
    # ms = getMs(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
    # 
    # # total number of basis functions should be this number
    # # return(rep(0, nx*ny*(4^nLayer - 1)/3))
    # rep(0, sum(ms))
    
    return(numeric(0))
  }
  
  log.norm.const = function()
  {
    ## let INLA compute it as -n/2 log(2pi) + 1/2 * log(|Q|)
    return (numeric(0))
  }
  
  log.prior = function() {
    require(invgamma)
    
    # get prior (note the Jacobian factors)
    # NOTE: unlike for the simple model, here "orig" versus "pc" refers to Dirichlet vs PC priors for alpha.
    # PC priors are always used for the marginal variance
    # corScalePar is the median of the prior belief of the spatial correlation scale parameter
    if(prior$priorType == "orig") {
      # rho = exp(theta)
      # P(theta < t) = P(rho < exp(t))
      # p_theta(t) = p_rho(exp(t)) |exp(t)|
      # log p_theta(theta) = log p_rho(rho) + log(rho)
      # 
      # sigmasq = exp(theta)
      # P(theta < t) = P(sigma < exp(0.5 t))
      # p_theta(t) = p_sigma(exp(t)) * |exp(t)|
      # log p_theta(theta) = log p_sigma(sigmasq) + log(sigmasq)
      out = sum(dinvexp(effectiveCor, rate=prior$corScalePar, log=TRUE)) + sum(log(effectiveCor)) + 
        invgamma::dinvgamma(sigmaSq, shape=prior$varPar1, rate=prior$varPar2, log=TRUE) + log(sigmaSq)
    } else if(prior$priorType == "pc") {
      # rho = exp(theta)
      # P(theta < t) = P(rho < exp(t))
      # p_theta(t) = p_rho(exp(t)) |exp(t)|
      # log p_theta(theta) = log p_rho(rho) + log(rho)
      # 
      # tau = exp(-theta) (here tau is 1/sigma^2, the precision)
      # P(theta < t) = P(tau < exp(-t))
      # p_theta(t) = p_tau(exp(-t)) * |-exp(-t)|
      # log p_theta(theta) = log p_tau(tau) + log(tau)
      out = sum(dinvexp(effectiveCor, rate=prior$corScalePar, log=TRUE)) + sum(log(effectiveCor)) + 
        inla.pc.dprec(1/sigmaSq, u=prior$u, alpha=prior$alpha, log=TRUE) - log(sigmaSq)
    }
    if(length(latInfo) == 1) {
      out
    } else {
      if(!exists("dexpitDirichlet")) {
        # load relevant external functions
        if(printVerboseTimings)
          print("sourcing LKinla.R...")
        inf = sessionInfo()
        if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
          source("~/git/LK-INLA/LKinla.R")
        else
          source("U:/git/LK-INLA/LKinla.R")
      }
      
      # the jacobian factor has already been included in dexpitDirichlet
      # Add on the density of a random vector whose multivariateExpit is Dirichlet
      out + dexpitDirichlet(multivariateLogit(alphas[1:(length(alphas)-1)]), prior$alphaPar, doLog=TRUE)
    }
  }
  
  initial = function() {
    # # initialize fixed effects and process variances
    # mod = lm(ys ~ X-1)
    # betas = coef(mod)
    # 
    # # remove from observations
    # ysCntr = ys - X %*% betas
    # 
    # # initialize covariance parameters
    # effectiveRangeInit = ((xRange[2] - xRange[1]) + (yRange[2] - yRange[1]))/2
    # xLength = xRange[2] - xRange[1]
    # squareWidth = xLength/(nx-1)
    # sigmaSq = sum((ysCntr - 0)^2)/n
    
    # get range of data
    xRangeDat = latInfo[[1]]$xRangeDat
    yRangeDat = latInfo[[1]]$yRangeDat
    
    # initialize process variance by estimating spatial process with OLS
    if(!exists("makeA")) {
      # load relevant external functions
      if(printVerboseTimings)
        print("sourcing LKinla.R...")
      inf = sessionInfo()
      if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)")
        source("~/git/LK-INLA/LKinla.R")
      else
        source("U:/git/LK-INLA/LKinla.R")
    }
    # AMat = makeA(datCoords, latInfo)
    theseYs = ys/ns
    # mod = lm(theseYs ~ cbind(X, as.matrix(AMat)) - 1)
    # r2 = summary(mod)$r.squared
    # s2 = summary(mod)$sigma^2
    # mod = lm.ridge(theseYs ~ cbind(X, AMat) - 1, lambda = seq(0,0.1,0.001))
    # mod = glmnet(cbind(X, as.matrix(AMat)), theseYs, alpha=0, intercept=FALSE)
    # r2 = mod$dev.ratio
    r2 = 1
    # sigmaSq = (r2/(1-r2)) * s2
    sigmaSq = var(theseYs) * r2
    
    # initialize covariance parameters
    if(is.null(initialEffectiveRange)) {
      # We want the "middle" representable effective range to be an eighth of the domain diameter
      middleEffectiveRange = ((xRangeDat[2] - xRangeDat[1]) + (yRangeDat[2] - yRangeDat[1]))/8
      
      # based on the middle effective range, set the effective ranges of each layer:
      nLayer = length(latInfo)
      exponents = seq(0.5*(nLayer-1), -0.5*(nLayer-1), by=-1)
      effectiveRangeInit = middleEffectiveRange * 2^exponents
    }
    else
      effectiveRangeInit = initialEffectiveRange
    
    # initialize the layer weights to be equal
    if(length(latInfo) == 1)
      c(log(effectiveRangeInit), log(sigmaSq))
    else {
      if(is.null(initialAlphas))
        initAlphas = rep(1/length(latInfo), length(latInfo)-1)
      else
        initAlphas = initialAlphas
      c(log(effectiveRangeInit), log(sigmaSq), multivariateLogit(initAlphas))
    }
  }
  
  quit = function() {
    return(invisible())
  }
  
  if(is.null(theta)) { theta = initial() }
  val = do.call(match.arg(cmd), args = list())
  
  totalTime = proc.time()[3] - startTime
  if(totalTime >= .1 && printVerboseTimings)
    print(paste0("rgeneric call with command ", cmd, " had the significant total time: ", totalTime))
  return (val)
}













