# given a dataset, fit LatticeKrig with package (alpha is fixed)
# obsCoords and obs: coordinates (nx2 matrix) and observed values (n-vector)
# nDat: resolution of coarsest lattice in largest direction over the domain 
# of the data (including the buffer there are more lattice points though)
# nLayer: number of layers to include
# NOTE: includes a first order linear polynomial in spatial coordinates, so do not include that in the covariates
fitLKSimple = function(obsCoords, obsValues, predCoords=obsCoords, xObs=NULL, xPred=NULL, NC=5, nLayer=3, normalize=FALSE, 
                       nBuffer=5, nu=1.5, verbose=TRUE, lambdaStart=.1, a.wghtStart=5, maxit=15, doSEs=TRUE, significanceCI=.8) {
  # if(!simpleMod)
  #   stop("non simple models not yet supported for LatticeKrig prediction function")
  
  # set up (the lattice, the arguments to LatticeKrig)
  LKinfo = LKrigSetup(obsCoords, nlevel=nLayer, nu=nu, NC=NC, normalize=normalize, 
                      lambda=lambdaStart, a.wght=a.wghtStart)
  
  # if are missing the predictive covariates, predict them using the observation covariates
  if(!is.null(xObs) && is.null(xPred)) {
    if(!is.matrix(xObs))
      xObs = as.matrix(xObs)
    xPred = matrix(nrow=nrow(predCoords), ncol=ncol(xObs))
    
    # make LatticeKrig predictive model for each covariate
    print(paste0("predicting covariates..."))
    for(i in 1:ncol(xObs)) {
      print(paste0("predicting covariate ", i, "/", ncol(xObs)))
      thisCov = xObs[,i]
      out = fitLKSimple(obsCoords, thisCov, NC=NC, nLayer=nLayer, normalize=normalize, 
                        nBuffer=nBuffer, nu=nu, verbose=verbose, lambdaStart=lambdaStart, 
                        a.wghtStart=a.wghtStart, maxit=maxit, doSEs=FALSE)
      xPred[,i] = out$preds
    }
  }
  
  # Maximum Likelihood
  if(verbose)
    print("Beginning LatticeKrig MLE fit...")
  # LKMLE = LKrig.MLE(obsCoords, obs, LKinfo=LKinfo)
  if(is.null(xObs) || is.null(xPred))
    LKMLE = LKrigFindLambdaAwght(obsCoords, obsValues, LKinfo=LKinfo, verbose=verbose, maxit=maxit)
  else
    LKMLE = LKrigFindLambdaAwght(obsCoords, obsValues, LKinfo=LKinfo, verbose=verbose, maxit=maxit, Z=xObs)
  if(verbose)
    print(LKMLE$summary)
  
  # final fit
  if(is.null(xObs) || is.null(xPred)) {
    out = LKrig(obsCoords, obsValues, LKinfo=LKMLE$LKinfo)
    preds = predict.LKrig(out, predCoords)
    if(doSEs)
      predSEs = predictSE.LKrig(out, predCoords)
    else
      predSEs = NULL
  }
  else {
    out = LKrig(obsCoords, obsValues, LKinfo=LKMLE$LKinfo, Z=xObs)
    preds = predict.LKrig(out, predCoords, Znew=xPred)
    if(doSEs)
      predSEs = predictSE.LKrig(out, predCoords, Znew=xPred)
    else
      predSEs=NULL
  }
  
  # calculate confidence intervals
  lower = preds + qnorm((1 - significanceCI) / 2, sd=predSEs)
  upper = preds + qnorm(1 - (1 - significanceCI) / 2, sd=predSEs)
  
  ## preds
  ## sigmas
  ## lower
  ## upper
  ## interceptSummary
  ## rangeSummary
  ## sdSummary
  ## varSummary
  return(list(mod=out, preds=preds, sigmas=predSEs, lower=lower, upper=upper, LKinfo=LKinfo))
}

# given a dataset, fit LatticeKrig optimizing over alpha and other parameters 
# obsCoords and obsValues: coordinates (nx2 matrix) and observed values (n-vector)
# nDat: resolution of coarsest lattice in largest direction over the domain 
# of the data (including the buffer there are more lattice points though)
# nLayer: number of layers to include
# NOTE: includes a first order linear polynomial in spatial coordinates by default, so do 
#       not include that in the covariates
# NOTE: lambda is (nugget variance) / rho
# NOTE: Due to profile likelihood method of LK, using the hessian of the likelihood for calculating 
#       parameter uncertainty will be anticonservative
fitLKStandard = function(obsCoords, obsValues, predCoords=obsCoords, xObs=NULL, xPred=NULL, NC=5, nLayer=3, normalize=TRUE, 
                         nBuffer=5, nu=1.5, verbose=TRUE, lambdaStart=.1, a.wghtStart=5, doSEs=TRUE, doCovSEs=doSEs, significanceCI=.8, 
                         lowerBoundLogLambda =-16,
                         upperBoundLogLambda = 4,
                         lowerBoundLogNu =-15,
                         upperBoundLogNu = 3,
                         lowerBoundLogitAlpha = rep(-10, nLayer-1),
                         upperBoundLogitAlpha= rep(10, nLayer-1),
                         lowerBoundOmega = -3,
                         upperBoundOmega = .75,
                         factr=1e7,
                         pgtol=1e-1,
                         maxit=15, 
                         nsimConditional=100, 
                         fixedFunctionArgs = list(m = 1), 
                         xRangeDat=NULL, yRangeDat=NULL, 
                         separatea.wght=FALSE, 
                         doMatern=FALSE, 
                         fixNu=FALSE) {
  
  if(separatea.wght) {
    a.wghtStart = as.list(rep(a.wghtStart, nLayer))
    lowerBoundOmega = rep(lowerBoundOmega, nLayer)
    upperBoundOmega = rep(upperBoundOmega, nLayer)
  }
  
  if(fixedFunctionArgs$m == 0) {
    # fixedFunction = function(x, Z=NULL, m=0, distance.type = "Euclidean") {
    #   # matrix(0, nrow=nrow(x), ncol=1)
    #   NULL
    # }
    fixedFunction = NULL
  } else {
    fixedFunction = "LKrigDefaultFixedFunction"
  }
  assign("fixedFunction", fixedFunction, envir=.GlobalEnv)
  
  # set spatial domain if not already set
  if(is.null(xRangeDat))
    xRangeDat = range(c(obsCoords[,1], predCoords[,1]))
  if(is.null(yRangeDat))
    yRangeDat = range(c(obsCoords[,2], predCoords[,2]))
  domainCoords = cbind(xRangeDat, yRangeDat)
  
  # if are missing the predictive covariates, predict them using the observation covariates
  if(!is.null(xObs) && is.null(xPred)) {
    if(!is.matrix(xObs))
      xObs = as.matrix(xObs)
    xPred = matrix(nrow=nrow(predCoords), ncol=ncol(xObs))
    
    # make LatticeKrig predictive model for each covariate
    print(paste0("predicting covariates..."))
    for(i in 1:ncol(xObs)) {
      print(paste0("predicting covariate ", i, "/", ncol(xObs)))
      thisCov = xObs[,i]
      out = fitLKStandard(obsCoords, thisCov, NC=NC, nLayer=nLayer, simpleMod=simpleMod, normalize=normalize, 
                          nBuffer=nBuffer, nu=nu, verbose=verbose, lambdaStart=lambdaStart, fixedFunctionArgs=fixedFunctionArgs, 
                          a.wghtStart=a.wghtStart, maxit=maxit, doSEs=FALSE, doMatern=doMatern, separatea.wght=separatea.wght, 
                          xRangeDat=xRangeDat, yRangeDat=yRangeDat)
      xPred[,i] = out$preds
    }
  }
  
  # do initial latticeKrig fit
  # set up the lattice, the arguments to LatticeKrig
  LKinfoStart = LKrigSetup(domainCoords, nlevel=nLayer, nu=nu, NC=NC, normalize=normalize, NC.buffer=nBuffer, 
                           lambda=lambdaStart, a.wght=a.wghtStart, fixedFunctionArgs=fixedFunctionArgs, alpha=rep(NA, nLayer))
  
  # make a function to convert from a vector of parameters to a corresponding set of different, named parameters
  getParameters = function(parameters) {
    # omega =  log( a.wght -4)/2
    # transform from optimized parameters to probabilities summing to 1 to get alphas
    if(nLayer != 1) {
      if(!doMatern) {
        thisnu = NULL
        alphas = multivariateExpit(parameters[1:(nLayer-1)])
        alphas = c(alphas, 1 - sum(alphas))
        log.lambda = parameters[nLayer-1 + 1]
        
        if(!separatea.wght)
          omega = parameters[nLayer-1 + 2]
        else
          omega = parameters[(nLayer-1 + 2):(2*nLayer)]
      } else {
        if(fixNu) {
          thisnu = nu
          alphas = getAlphas(nLayer, thisnu)
          log.lambda = parameters[1]
          
          if(!separatea.wght)
            omega = parameters[2]
          else
            omega = parameters[2:(1 + nLayer)]
        } else {
          thisnu = exp(parameters[1])
          alphas = getAlphas(nLayer, thisnu)
          log.lambda = parameters[2]
          
          if(!separatea.wght)
            omega = parameters[3]
          else
            omega = parameters[3:(2 + nLayer)]
        }
      }
      
      list(nu=thisnu, alphas=alphas, log.lambda=log.lambda, lambda=exp(log.lambda), omega=omega, a.wght=omega2Awght(omega, LKinfoStart))
    }
    else {
      log.lambda = parameters[1]
      omega = parameters[2]
      list(nu=NULL, alphas=1, log.lambda=log.lambda, lambda=exp(log.lambda), omega=omega, a.wght=omega2Awght(omega, LKinfoStart))
    }
  }
  
  # make a wrapper function around LKrig in order to optimize over all parameters including alpha
  outerFun = function(thesePar, thisVerbose=verbose) {
    parameterList = getParameters(thesePar)
    thisnu = parameterList$nu
    alphas = parameterList$alphas
    if(nLayer > 1)
      logitAlphas = multivariateLogit(alphas[1:(nLayer-1)])
    else
      logitAlphas = (lowerBoundLogitAlpha + upperBoundLogitAlpha)/2 # this value is just used to make sure we're in the optimization range
    log.lambda = parameterList$log.lambda
    lambda = parameterList$lambda
    omega = parameterList$omega
    a.wght = as.list(parameterList$a.wght)
    
    # make sure parameters are in the correct range
    if(log.lambda<lowerBoundLogLambda || log.lambda>upperBoundLogLambda  || 
       any(logitAlphas<lowerBoundLogitAlpha) || any(logitAlphas>upperBoundLogitAlpha) || 
       lowerBoundOmega<lowerBoundOmega || lowerBoundOmega>upperBoundOmega)
      return(-1000000)
    
    # set up the lattice, the arguments to LatticeKrig
    LKinfo = LKrigSetup(domainCoords, nlevel=nLayer, nu=thisnu, NC=NC, normalize=normalize, NC.buffer=nBuffer, 
                        lambda=lambda, a.wght=a.wght, alpha=alphas, 
                        fixedFunction=fixedFunction, fixedFunctionArgs=fixedFunctionArgs)
    
    # if(is.null(xObs) || is.null(xPred))
    #   LKMLE = LKrig(obsCoords, obsValues, LKinfo=LKinfo, verbose=verbose, maxit=maxit)
    # else
    #   LKMLE = LKrig(obsCoords, obsValues, LKinfo=LKinfo, verbose=verbose, maxit=maxit, Z=xObs)
    out = LKrig(obsCoords, obsValues, LKinfo=LKinfo, verbose=thisVerbose, Z=xObs)
    out$lnProfileLike.FULL
  }
  
  # get initial parameters and optimization bounds
  if(nLayer == 1) {
    init = c(log.lambda=log(lambdaStart), omega=Awght2Omega(unlist(a.wghtStart), LKinfoStart))
    lower = c(lowerBoundLogLambda, lowerBoundOmega)
    upper = c(upperBoundLogLambda, upperBoundOmega)
  } else {
    if(doMatern) {
      if(fixNu) {
        init = c(log.lambda=log(lambdaStart), omega=Awght2Omega(unlist(a.wghtStart), LKinfoStart))
        lower = c(lowerBoundLogLambda, lowerBoundOmega)
        upper = c(upperBoundLogLambda, upperBoundOmega)
      } else {
        init = c(log.nu=log(nu), log.lambda=log(lambdaStart), omega=Awght2Omega(unlist(a.wghtStart), LKinfoStart))
        lower = c(lowerBoundLogNu, lowerBoundLogLambda, lowerBoundOmega)
        upper = c(upperBoundLogNu, upperBoundLogLambda, upperBoundOmega)
      }
    } else {
      init = c(logit.alphas=multivariateLogit(rep(1 / nLayer, nLayer-1)), log.lambda=log(lambdaStart), omega=Awght2Omega(unlist(a.wghtStart), LKinfoStart))
      lower = c(lowerBoundLogitAlpha, lowerBoundLogLambda, lowerBoundOmega)
      upper = c(upperBoundLogitAlpha, upperBoundLogLambda, upperBoundOmega)
    }
  }
  
  # Maximum Likelihood
  if(verbose)
    print("Beginning LatticeKrig MLE fit...")
  
  result <- try(optim(init,
                      outerFun, 
                      lower=lower, 
                      upper=upper, 
                      # method="L-BFGS-B",
                                           method="BFGS",
                      control=list(fnscale = -1,factr=factr,
                                   pgtol=pgtol, maxit=maxit,
                                   ndeps = rep(.05,length(init)))
  ))
  if(verbose){
    cat("Results from optimize:", fill=TRUE)
    print( result )
  }
  
  # final fit
  parameterList = getParameters(result$par)
  nuMLE = parameterList$nu
  alphasMLE = parameterList$alphas
  log.lambdaMLE = parameterList$log.lambda
  lambdaMLE = parameterList$lambda
  omegaMLE = parameterList$omega
  a.wghtMLE = parameterList$a.wght
  
  # set up the lattice, the arguments to LatticeKrig, generate predictions and predictive uncertainties
  LKinfo = LKrigSetup(domainCoords, nlevel=nLayer, nu=nuMLE, NC=NC, normalize=normalize, NC.buffer=nBuffer, 
                      lambda=lambdaMLE, a.wght=as.list(a.wghtMLE), alpha=alphasMLE, 
                      fixedFunction=fixedFunction, fixedFunctionArgs=fixedFunctionArgs)
  if(is.null(xObs) || is.null(xPred)) {
    mod = LKrig(obsCoords, obsValues, LKinfo=LKinfo)
    preds = predict.LKrig(mod, predCoords)
    if(doSEs)
      predSimulations = LKrig.sim.conditional(mod, x.grid=predCoords, M=nsimConditional)
    else
      predSimulations = NULL
  } else {
    mod = LKrig(obsCoords, obsValues, LKinfo=LKinfo, Z=xObs)
    preds = predict.LKrig(mod, predCoords, Znew=xPred)
    if(doSEs)
      predSimulations = LKrig.sim.conditional(mod, x.grid=predCoords, Z.grid=xPred, M=nsimConditional)
    else
      predSimulations=NULL
  }
  
  # get predictive standard errors
  predSEsNoNugget = predSimulations$SE
  predSEs = sqrt(predSimulations$SE^2 + mod$sigma.MLE^2)
  
  # calculate prediction confidence intervals
  lower = preds + qnorm((1 - significanceCI) / 2, sd=predSEs)
  medians = preds
  upper = preds + qnorm(1 - (1 - significanceCI) / 2, sd=predSEs)
  lowerNoNugget = preds + qnorm((1 - significanceCI) / 2, sd=predSEsNoNugget)
  mediansNoNugget = preds
  upperNoNugget = preds + qnorm(1 - (1 - significanceCI) / 2, sd=predSEsNoNugget)
  
  ## now we calculate uncertainty intervals for all parameters
  # intercept
  if(fixedFunctionArgs$m >= 1) {
    interceptSummary = c(Est=mod$d.coef[1], SD=sd(predSimulations$d.coef.draw[1,]),
                         Qlower=quantile(probs=(1 - significanceCI) / 2, predSimulations$d.coef.draw[1,]),
                         Q50=quantile(probs=.5, predSimulations$d.coef.draw[1,]),
                         Qupper=quantile(probs=1 - (1 - significanceCI) / 2, predSimulations$d.coef.draw[1,]))
  } else
    interceptSummary = c(Est=0, SD=0, Qlower=0, Q50=0, Qupper=0)
  
  # to calculate summaries for the parameters, must calculate inverse of negative hessian
  # NOTE: Due to profile likelihood method of LK, using the hessian in this way for getting 
  #       uncertainty estimates will be anticonservative
  if(doCovSEs) {
    require(numDeriv)
    print("Calculating hessian...")
    hess = hessian(outerFun, result$par, thisVerbose=FALSE)
    parSigma = solve(-hess)
  }
  
  
  # make a function to transform from a vector of parameters to a final set of different, named parameters
  getParametersFinal = function(parameters) {
    # omega =  log( a.wght -4)/2
    # transform from optimized parameters to probabilities summing to 1 to get alphas
    if(nLayer != 1) {
      if(!doMatern) {
        thisnu = NULL
        alphas = multivariateExpit(parameters[1:(nLayer-1)])
        alphas = c(alphas, 1 - sum(alphas))
        log.lambda = parameters[nLayer-1 + 1]
        
        if(!separatea.wght)
          omega = parameters[nLayer-1 + 2]
        else
          omega = parameters[(nLayer-1 + 2):(2*nLayer)]
        
        c(alphas=alphas, lambda=exp(log.lambda), a.wght=omega2Awght(omega, LKinfoStart))
      } else {
        if(fixNu) {
          thisnu = nu
          alphas = getAlphas(nLayer, thisnu)
          log.lambda = parameters[1]
          
          if(!separatea.wght)
            omega = parameters[2]
          else
            omega = parameters[2:(1 + nLayer)]
          
          c(alphas=alphas, lambda=exp(log.lambda), a.wght=omega2Awght(omega, LKinfoStart))
        } else {
          thisnu = exp(parameters[1])
          alphas = getAlphas(nLayer, thisnu)
          log.lambda = parameters[2]
          
          if(!separatea.wght)
            omega = parameters[3]
          else
            omega = parameters[3:(2 + nLayer)]
          
          c(nu=thisnu, alphas=alphas, lambda=exp(log.lambda), a.wght=omega2Awght(omega, LKinfoStart))
        }
      }
    }
    else {
      log.lambda = parameters[1]
      omega = parameters[2]
      c(alphas=1, lambda=exp(log.lambda), a.wght=omega2Awght(omega, LKinfoStart))
    }
  }
  
  # simulate possible parameter values and do any necessary transformations for any parameters from the 
  # optimization scale
  if(doCovSEs)
    U = try(chol(parSigma))
  
  if(doCovSEs && class(U) != "try-error") {
    L = t(U)
    zSim = matrix(rnorm(nsimConditional * nrow(L)), nrow=nrow(L))
    parSim = L %*% zSim
    parSim = sweep(parSim, 1, result$par, "+")
    
    finalParSim = apply(parSim, 2, getParametersFinal)
    
    getSummaryStatistics = function(draws) {
      c(Est=mean(draws), SD=sd(draws), 
        Qlower=quantile(probs=(1 - significanceCI) / 2, draws), 
        Q50=quantile(probs=0.5, draws), 
        Qupper=quantile(probs=1 - (1 - significanceCI) / 2, draws))
    }
    
    parameterSummaryTable = t(apply(finalParSim, 1, getSummaryStatistics))
    summaryNames = c("Est", "SD", "Qlower", "Q50", "Qupper")
    colnames(parameterSummaryTable) = summaryNames
  } else {
    if(doCovSEs)
      warning("bad hessian: fixing singular hyperparameter distribution to the estimates")
    finalParSim = matrix(rep(getParametersFinal(result$par), nsimConditional), ncol=nsimConditional)
    parameterSummaryTable = cbind(getParametersFinal(result$par), NA, NA, NA, NA)
  }
  
  totalVariance = mod$rho.MLE + mod$sigma.MLE^2
  
  if(!separatea.wght) {
    if(!doMatern) {
      lambdaVals = finalParSim[nLayer + 1,]
      a.wghtVals = finalParSim[nrow(finalParSim),]
      alphaVals = finalParSim[1:nLayer,]
      nuVals = NULL
    } else {
      lambdaVals = finalParSim[nLayer + 2,]
      a.wghtVals = finalParSim[nrow(finalParSim),]
      alphaVals = finalParSim[2:(nLayer+1),]
      nuVals = finalParSim[1,]
    }
  } else {
    if(!doMatern) {
      lambdaVals = finalParSim[nLayer + 1,]
      alphaVals = finalParSim[1:nLayer,]
      nuVals = NULL
    } else {
      lambdaVals = finalParSim[nLayer + 2,]
      alphaVals = finalParSim[2:(nLayer+1),]
      nuVals = finalParSim[1,]
    }
    a.wghtVals = finalParSim[(nrow(finalParSim) - nLayer + 1):nrow(finalParSim),]
  }
  
  rhoVals = (1 / lambdaVals) * totalVariance / (1 + 1 / lambdaVals)
  nuggetVarVals = lambdaVals * totalVariance / (1 + lambdaVals)
  
  ## preds
  ## sigmas
  ## lower
  ## upper
  ## interceptSummary
  ## rangeSummary
  ## sdSummary
  ## varSummary
  return(list(mod=mod, preds=preds, sigmas=predSEs, lower=lower, medians=medians, upper=upper, 
              sigmasNoNugget=predSEsNoNugget, lowerNoNugget=lowerNoNugget, mediansNoNugget=mediansNoNugget, upperNoNugget=upperNoNugget, 
              parameterSummaryTable=parameterSummaryTable, LKinfo=LKinfo, 
              interceptSummary=interceptSummary, rangeSummary=c(), sdSummary=c(), varSummary=c(), parSim=finalParSim, fixHyperpar=is.na(parameterSummaryTable[1,5]), 
              rhoVals=rhoVals, nuggetVarVals=nuggetVarVals, lambdaVals=lambdaVals, alphaVals=alphaVals, nuVals=nuVals, a.wghtVals=a.wghtVals, totalVariance=totalVariance, 
              predMat=predSimulations$g.draw))
}

# modified version of LKrigMakewU from LatticeKrig to prevent bad indexing of fixedFunctionArgs instead of fixedFunction
# (R has indexing issues)
LKrigMakewU = function(object, verbose = FALSE) 
{
  LKinfo <- object$LKinfo
  if (!is.null(object$U)) {
    wU <- sqrt(object$weights) * object$U
  }
  else {
    if (!is.null(LKinfo[["fixedFunction"]])) {
      wU <- sqrt(object$weights) * do.call(LKinfo$fixedFunction, 
                                           c(list(x = object$x, Z = object$Z, distance.type = LKinfo$distance.type), 
                                             LKinfo$fixedFunctionArgs))
    }
    else {
      wU <- NULL
    }
  }
  if (verbose) {
    cat("dim wU:", dim(wU), fill = TRUE)
  }
  return(wU)
}
assign("LKrigMakewU", LKrigMakewU, .GlobalEnv)
assignInNamespace("LKrigMakewU", LKrigMakewU, "LatticeKrig")











