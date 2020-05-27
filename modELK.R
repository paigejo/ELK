# fit ELK model given data and relevant parameters
# obsCoords: matrix of spatial coordinates of the data
# obsValues: matrix of response at the observation coordinates
# predCoords: places at which to make predictions
# xObs: observation design matrix (intercept and covariates at obs locations)
# xPred: prediction design matrix (intercept and covariates at pred locations)
# NC: number of coarse lattice points over the longest dimension of the data
# nLayer: number of lattice layers
# seed: random seed set at beginning of function
# normalize: whether or not to normalize the process (highly recommended to set this to TRUE, 
#            since otherwise interpretation of correlation length and marginal variance priors 
#            will be highly suspect)
# intStrategy: see INLA:::control.inla()
# strategy: see INLA:::control.inla()
# predClusterI: whether or not to include cluster effects for which predictions
# priorPar default was getPrior(.1, .1, 10, nLayer=nLayer)
fitELK = function(obsCoords, obsValues, predCoords=obsCoords, seed=1, nLayer=3, NC=5,
                  nBuffer=5, priorPar=NULL, 
                  xObs=cbind(1, obsCoords), xPred=cbind(1, predCoords), normalize=TRUE, 
                  intStrategy="ccd", strategy="gaussian", fastNormalize=TRUE, 
                  predictionType=c("mean", "median"), significanceCI=0.8, 
                  printVerboseTimings=FALSE, nPostSamples=1000, family=c("normal", "binomial", "betabinomial"),
                  obsNs=rep(1, length(obsValues)), clusterEffect=TRUE, latInfo=NULL, 
                  initialEffectiveRange=NULL, initialAlphas=rep(1/nLayer, nLayer-1), 
                  effRangeRange=NULL, predClusterI=rep(TRUE, nrow(predCoords)), 
                  plotNormalizationSplines=FALSE, verbose=TRUE, separateRanges=FALSE, 
                  doValidation=FALSE, previousFit=NULL, precomputedNormalizationFun=NULL, 
                  savePrecomputationResults=FALSE, loadPrecomputationResults=FALSE, 
                  precomputationFileNameRoot="precomputationResults") {
  
  startTime = proc.time()[3]
  set.seed(seed)
  
  if(!normalize)
    warning("WARNING: normalize set to FALSE. Correlation length and marginal variance prior interpretation is suspect. It is highly recommended to set normalize to TRUE")
  
  # get the type of prediction the user wants
  predictionType = match.arg(predictionType)
  family = match.arg(family)
  
  # check to make sure inputs make sense
  if(!clusterEffect && family == "normal")
    stop("cluster effect must be included for the normal family")
  
  # get lattice points, prediction points
  
  if(is.null(latInfo)) {
    xRangeDat = range(obsCoords[,1])
    yRangeDat = range(obsCoords[,2])
    latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  } else {
    xRangeDat = latInfo[[1]]$xRangeDat
    yRangeDat = latInfo[[1]]$yRangeDat
  }
  nx = latInfo[[1]]$nx
  ny = latInfo[[1]]$ny
  
  if(is.null(priorPar))
    priorPar = getPCPrior(1444.772/5, .01, 1, nLayer=nLayer, separateRanges=separateRanges, latticeInfo=latInfo) # 1444.772/5
  
  # generate lattice basis matrix
  AObs = makeA(obsCoords, latInfo)
  
  ## run precomputations
  print("running precomputations...")
  startTimePrecomputations = proc.time()[3]
  
  # either run the precomputations or load them in
  if(!loadPrecomputationResults) {
    precomputedMatrices = precomputationsQ2(latInfo)
    if(is.null(precomputedNormalizationFun)) {
      precomputedNormalizationFun = precomputeNormalization(saveResults=FALSE, latticeInfo=latInfo, effRangeRange=effRangeRange, 
                                                            plotNormalizationSplines=plotNormalizationSplines)
    }
  } else {
    load(paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
  }
  
  # save the precomputations if necessary
  if(savePrecomputationResults) {
    save(precomputedMatrices, precomputedNormalizationFun, 
         file=paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
  }
  
  endTimePrecomputations = proc.time()[3]
  precomputationTime = endTimePrecomputations - startTimePrecomputations
  print(paste0("finished precomputations. Took ", round(precomputationTime / 60, 2), " minutes"))
  
  # define the LK model
  startTimeDefineModel = proc.time()[3]
  
  # set beta binomial prior if necessary
  if(family == "betabinomial") {
    if(clusterEffect)
      stop("cluster effect must not be set to TRUE for betabinomial model")
    
    # set median at .04 and upper 97.5th pctile at 0.2
    mu = logit(0.04)
    prec = 1/((logit(.2)-logit(.04))/qnorm(.975))^2
    control.family = list(hyper = list(rho = list(prior="logtnormal", param=c(mu, prec))))
  } else {
    control.family = list()
  }
  
  n = length(obsValues)
  if(separateRanges) {
    rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.full, latInfo=latInfo, ys=obsValues, 
                                prior=priorPar, normalize=normalize, precomputedMatrices=precomputedMatrices, 
                                X=xObs, datCoords=obsCoords, fastNormalize=fastNormalize, 
                                printVerboseTimings=printVerboseTimings, 
                                initialEffectiveRange=initialEffectiveRange, initialAlphas=initialAlphas, 
                                precomputedNormalizationFun=precomputedNormalizationFun, ns=obsNs)
  } else {
    rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.standard, latInfo=latInfo, ys=obsValues, 
                                prior=priorPar, normalize=normalize, precomputedMatrices=precomputedMatrices, 
                                X=xObs, datCoords=obsCoords, fastNormalize=fastNormalize, 
                                printVerboseTimings=printVerboseTimings, 
                                initialEffectiveRange=initialEffectiveRange, initialAlphas=initialAlphas, 
                                precomputedNormalizationFun=precomputedNormalizationFun, ns=obsNs)
  }
  
  # use these global variables for testing calls to inla.rgeneric.lk.model.simple
  # latInfo<<-latInfo; ys<<-obsValues; ns<<-obsNs;
  # prior<<-priorPar; normalize<<-normalize; precomputedMatrices<<-precomputedMatrices;
  # X<<-xObs; nu<<-nu; datCoords<<-obsCoords; fastNormalize<<-fastNormalize;
  # printVerboseTimings<<-printVerboseTimings; initialEffectiveRange<<-initialEffectiveRange;
  # initialAlphas<<-initialAlphas; precomputedNormalizationFun<<-precomputedNormalizationFun
  ## generate inla stack:
  # Stacked A matrix (A_s from notation of LR2015 Bayesian Spatial Modelling with R_INLA):
  # (AEst   0  )
  # ( 0   APred)
  # eta_s = (c^T c^T)^T
  # where c is the vector of lattice coefficients
  AEst = makeA(obsCoords, latInfo)
  APred = makeA(predCoords, latInfo)
  latticeInds = 1:ncol(AEst)
  clust = 1:length(obsValues)
  if(family == "normal") {
    if(!is.null(xObs)) {
      stack.est = inla.stack(A =list(AEst, 1), 
                             effects =list(field=latticeInds, X=xObs), 
                             data =list(y=obsValues, link=1), 
                             tag ="est", 
                             remove.unused=FALSE)
      stack.pred = inla.stack(A =list(APred, 1), 
                              effects =list(field=latticeInds, X=xPred), 
                              data =list(y=NA, link=1), 
                              tag ="pred", 
                              remove.unused=FALSE)
    } else {
      stack.est = inla.stack(A =list(AEst), 
                             effects =list(field=latticeInds), 
                             data =list(y=obsValues, link=1), 
                             tag ="est", 
                             remove.unused=FALSE)
      stack.pred = inla.stack(A =list(APred), 
                              effects =list(field=latticeInds), 
                              data =list(y=NA, link=1), 
                              tag ="pred", 
                              remove.unused=FALSE)
    }
  }
  else if(family == "binomial" || family == "betabinomial") {
    if(clusterEffect) {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1, 1), 
                               effects =list(field=latticeInds, clust=clust, X=xObs), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      } else {
        stack.est = inla.stack(A =list(AEst, 1), 
                               effects =list(field=latticeInds, clust=clust), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      }
    } else {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1), 
                               effects =list(field=latticeInds, X=xObs), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      } else {
        stack.est = inla.stack(A =list(AEst), 
                               effects =list(field=latticeInds), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      }
    }
    if(!is.null(xObs)) {
      stack.pred = inla.stack(A =list(APred, 1), 
                              effects =list(field=latticeInds, X=xPred), 
                              data =list(y=NA, link=1, Ntrials=rep(1, nrow(APred)), 
                                         tag ="pred"), 
                              remove.unused=FALSE)
    } else {
      stack.pred = inla.stack(A =list(APred), 
                              effects =list(field=latticeInds), 
                              data =list(y=NA, link=1, Ntrials=rep(1, nrow(APred)), 
                                         tag ="pred"), 
                              remove.unused=FALSE)
    }
  }
  stack.full = inla.stack(stack.est, stack.pred, 
                          remove.unused=FALSE)
  dat = inla.stack.data(stack.full, rgen=rgen, remove.unused=FALSE)
  
  # initialize the fitting process based on a previous optimum if necessary
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousFit)) {
    # modeControl$result = previousFit
    modeControl$theta = previousFit$mode$theta
    modeControl$x = previousFit$mode$x
    modeControl$restart = TRUE
  }
  endTimeDefineModel = proc.time()[3]
  totalTimeDefineModel = endTimeDefineModel - startTimeDefineModel
  
  # fit the model
  # control.inla = list(cmin = 0, int.strategy=int.strategy) 
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  controls = list(strategy=strategy, int.strategy=intStrategy) 
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  startTimeFitModel = proc.time()[3]
  if(family == "normal") {
    if(!is.null(xObs)) {
      mod = inla(y ~ - 1 + X + f(field, model=rgen), 
                 data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                 control.inla=controls, 
                 control.mode=modeControl, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                 control.fixed=list(quantiles=allQuantiles), 
                 control.family=list(hyper = list(prec = list(param=c(1, 0.05), prior="pc.prec"))))
    } else {
      mod = inla(y ~ - 1 + f(field, model=rgen), 
                 data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                 control.inla=controls, 
                 control.mode=modeControl, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                 control.fixed=list(quantiles=allQuantiles), 
                 control.family=list(hyper = list(prec = list(param=c(1, 0.05), prior="pc.prec"))))
    }
  } else if(family == "binomial" || family == "betabinomial") {
    if(clusterEffect) {
      # clusterList = list(param=c(.15, 0.01), prior="pc.prec")
      if(!is.null(xObs)) {
        mod = inla(y ~ - 1 + X + f(field, model=rgen) + 
                     f(clust, model="iid", hyper = list(prec = list(param=c(1, 0.01), prior="pc.prec"))), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      } else {
        mod = inla(y ~ - 1 + f(field, model=rgen) + 
                     f(clust, model="iid", hyper = list(prec = list(param=c(1, 0.01), prior="pc.prec"))), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      }
    } else {
      if(!is.null(xObs)) {
        mod = inla(y ~ - 1 + X + f(field, model=rgen), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      } else {
        mod = inla(y ~ - 1 + f(field, model=rgen), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE, quantiles=allQuantiles), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      }
    }
  }
  endTimeFitModel = proc.time()[3]
  totalTimeFitModel = endTimeFitModel - startTimeFitModel
  
  print(paste0("finished fitting model. Took ", round(totalTimeFitModel / 60, 2), " minutes"))
  
  # # improve the approximation of the posterior if requested by the user
  # if(improveHyperpar) {
  #   browser()
  #   mod = inla.hyperpar(mod)
  # }
  
  # browser()
  
  # get predictive surface, SD, and data
  index = inla.stack.index(stack.full, "pred")$data
  obsInds = 1:n
  predInds = (n+1):(n+nrow(predCoords))
  if(predictionType == "mean") {
    linpreds = mod[["summary.linear.predictor"]]$mean
    # preds = linpreds[predInds]
    # obsPreds = linpreds[obsInds]
  } else {
    linpreds = mod[["summary.linear.predictor"]]$`0.5quant`
    # preds = linpreds[predInds]
    # obsPreds = linpreds[obsInds]
  }
  linpred.sd = mod[["summary.linear.predictor"]]$sd
  # predSDs = linpred.sd[predInds]
  # obsSDs = linpred.sd[obsInds]
  
  # generate samples from posterior
  print("Sampling from posterior")
  startTimePosteriorSampling = proc.time()[3]
  postSamples = inla.posterior.sample(nPostSamples, mod)
  endTimePosteriorSampling = proc.time()[3]
  totalTimePosteriorSampling = endTimePosteriorSampling - startTimePosteriorSampling
  print(paste0("finished sampling from the posterior. Took ", round(totalTimePosteriorSampling / 60, 2), " minutes"))
  
  print("Processing posterior samples...")
  startTimeSampleProcessing = proc.time()[3]
  # get posterior hyperparameter samples and transform them as necessary
  # interpretation of hyperparameters: 
  # 1: error precision
  # 2: log effective range
  # 3: log spatial variance
  # 4-(3 + nLayer - 1): multivariateLogit alpha
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  if(family == "normal") {
    if(!separateRanges) {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[3])+1/x[1], spatialVar=exp(x[3]), errorVar=1/x[1], 
                                              totalSD=sqrt(exp(x[3])+1/x[1]), spatialSD=sqrt(exp(x[3])), errorSD=sqrt(1/x[1]), 
                                              spatialRange=exp(x[2]), alpha=multivariateExpit(x[4:(3 + nLayer - 1)]))})
      mat = rbind(mat, alpha=1-colSums(mat[8:(7+nLayer-1),]))
      hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", "spatialRange", 
                     paste0("alpha", 1:nLayer))
    } else {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[2+nLayer])+1/x[1], spatialVar=exp(x[2+nLayer]), errorVar=1/x[1], 
                                              totalSD=sqrt(exp(x[2+nLayer])+1/x[1]), spatialSD=sqrt(exp(x[2+nLayer])), errorSD=sqrt(1/x[1]), 
                                              spatialRange=exp(x[2:(1+nLayer)]), alpha=multivariateExpit(x[(3+nLayer):(3+2*nLayer-2)]))})
      mat = rbind(mat, alpha=1-colSums(matrix(mat[(6+nLayer+1):(6+nLayer+1 + nLayer-2),], nrow=nLayer-1)))
      hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", paste0("spatialRange", 1:nLayer), 
                     paste0("alpha", 1:nLayer))
    }
  } else if(family == "binomial") {
    if(clusterEffect) {
      if(!separateRanges) {
        logSpatialRangeI = 1
        logSpatialVarI = 2
        logitAlphaI = 3:(2 + nLayer - 1)
        precisionI = 3 + nLayer - 1
        mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[logSpatialVarI])+1/x[precisionI], spatialVar=exp(x[logSpatialVarI]), clusterVar=1/x[precisionI], 
                                                totalSD=sqrt(exp(x[logSpatialVarI])+1/x[precisionI]), spatialSD=sqrt(exp(x[logSpatialVarI])), clusterSD=sqrt(1/x[precisionI]), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        mat = rbind(mat, alpha=1-colSums(mat[8:(7+nLayer-1),]))
        hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", "spatialRange", 
                       paste0("alpha", 1:nLayer))
      } else {
        logSpatialRangeI = 1:nLayer
        logSpatialVarI = 1 + nLayer
        logitAlphaI = (logSpatialVarI+1):(logSpatialVarI+1 + nLayer-2)
        precisionI = max(logitAlphaI)+1
        mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[logSpatialVarI])+1/x[precisionI], spatialVar=exp(x[logSpatialVarI]), clusterVar=1/x[precisionI], 
                                                totalSD=sqrt(exp(x[logSpatialVarI])+1/x[precisionI]), spatialSD=sqrt(exp(x[logSpatialVarI])), clusterSD=sqrt(1/x[precisionI]), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        mat = rbind(mat, alpha=1-colSums(matrix(mat[(6+nLayer+1):(6+nLayer+1 + nLayer-2),], nrow=nLayer-1)))
        hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", paste0("spatialRange", 1:nLayer), 
                       paste0("alpha", 1:nLayer))
      }
    } else {
      if(!separateRanges) {
        logSpatialRangeI = 1
        logSpatialVarI = 2
        logitAlphaI = 3:(2 + nLayer - 1)
        mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        mat = rbind(mat, alpha=1-colSums(matrix(mat[4:(3+nLayer-1),], nrow=nLayer-1)))
        hyperNames = c("spatialVar", "spatialSD", "spatialRange", paste0("alpha", 1:nLayer))
      } else {
        logSpatialRangeI = 1:nLayer
        logSpatialVarI = 1 + nLayer
        logitAlphaI = (logSpatialVarI+1):(logSpatialVarI+1 + nLayer-2)
        mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        # mat = rbind(mat, alpha=1-colSums(mat[4:(3+nLayer-1),]))
        mat = rbind(mat, alpha=1-colSums(matrix(mat[(3+nLayer):(3+nLayer + nLayer-2),], nrow=nLayer-1)))
        hyperNames = c("spatialVar", "spatialSD", paste0("spatialRange", 1:nLayer), paste0("alpha", 1:nLayer))
      }
    }
  } else if(family == "betabinomial") {
    overdispersionI = 1
    if(!separateRanges) {
      logSpatialRangeI = 2
      logSpatialVarI = 3
      logitAlphaI = 4:(3 + nLayer - 1)
      mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                              spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]), 
                                              overdispersion=x[overdispersionI])})
      mat = rbind(mat[-nrow(mat),], alpha=1-colSums(matrix(mat[4:(3+nLayer-1),], nrow=nLayer-1)), mat[nrow(mat),])
      hyperNames = c("spatialVar", "spatialSD", "spatialRange", paste0("alpha", 1:nLayer), "overdispersion")
    } else {
      logSpatialRangeI = 2:(1+nLayer)
      logSpatialVarI = 2 + nLayer
      logitAlphaI = (logSpatialVarI+1):(logSpatialVarI+1 + nLayer-2)
      mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                              spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]), 
                                              overdispersion=x[overdispersionI])})
      # mat = rbind(mat, alpha=1-colSums(mat[4:(3+nLayer-1),]))
      mat = rbind(mat[-nrow(mat),], alpha=1-colSums(matrix(mat[(3+nLayer):(3+nLayer + nLayer-2),], nrow=nLayer-1)), mat[nrow(mat),])
      hyperNames = c("spatialVar", "spatialSD", paste0("spatialRange", 1:nLayer), paste0("alpha", 1:nLayer), "overdispersion")
    }
  }
  
  rownames(mat) = hyperNames
  
  getSummaryStatistics = function(draws) {
    c(Est=mean(draws), SD=sd(draws), 
      Qlower=quantile(probs=(1 - significanceCI) / 2, draws), 
      Q50=quantile(probs=0.5, draws), 
      Qupper=quantile(probs=1 - (1 - significanceCI) / 2, draws))
  }
  summaryNames = c("Est", "SD", "Qlower", "Q50", "Qupper")
  parameterSummaryTable = t(apply(mat, 1, getSummaryStatistics))
  colnames(parameterSummaryTable) = summaryNames
  
  # separate out default parameter summaries
  varI = which(grepl("clusterVar", hyperNames))
  rangeI = which(grepl("Range", hyperNames))
  alphaI = which(grepl("alpha", hyperNames))
  varSummary=parameterSummaryTable[varI,]
  rangeSummary=parameterSummaryTable[rangeI,]
  alphaSummary=parameterSummaryTable[alphaI,]
  if(family == "normal" || clusterEffect) {
    sdI = which(grepl("clusterSD", hyperNames))
    sdSummary=parameterSummaryTable[sdI,]
  } else {
    sdSummary = matrix(rep(0, 5), nrow=1)
  }
  overdispersionSummary = matrix(rep(0, 5), nrow=1)
  if(family == "betabinomial")
    overdispersionSummary=parameterSummaryTable[nrow(parameterSummaryTable),]
  
  # get samples of the latent field
  latentMat = sapply(postSamples, function(x) {x$latent})
  latentVarNames = rownames(postSamples[[1]]$latent)
  fieldIndices = which(grepl("field", latentVarNames))
  fixedIndices = which(grepl("X", latentVarNames))
  # if(clusterEffect)
  #   clustIndices = grepl("clust", latentVarNames)
  
  ## generate logit predictions (first without cluster effect then add the cluster effect in)
  # for prediction locations
  if(length(xPred) != 0)
    fixedPart = xPred  %*% latentMat[fixedIndices,]
  else
    fixedPart = 0
  predMat = fixedPart + APred %*% latentMat[fieldIndices,]
  
  # for observation locations
  if(length(xObs) != 0)
    fixedPart = xObs  %*% latentMat[fixedIndices,]
  else
    fixedPart = 0
  obsMat = fixedPart + AObs %*% latentMat[fieldIndices,]
  
  # get draws from basis function coefficients
  basisCoefMat = latentMat[fieldIndices,]
  
  # add in cluster effect if necessary
  if((family == "binomial" && clusterEffect) || family == "normal") {
    # get betabinomial overdispersion parameter
    clusterVarI = which(grepl("clusterVar", hyperNames))
    clusterVars = mat[clusterVarI,]
    predMatClustEffect = predMat + sweep(matrix(rnorm(length(predMat), sd=rep(sqrt(clusterVars), each=nrow(predMat))), nrow=nrow(predMat)), 1, predClusterI, "*")
    obsMatClustEffect = obsMat + matrix(rnorm(length(obsMat), sd=rep(sqrt(clusterVars), each=nrow(obsMat))), nrow=nrow(obsMat))
    rhos = NULL
  } else if(family == "betabinomial") {
    # get cluster induced overdispersion
    overdispersionI = which(grepl("overdispersion", hyperNames))
    rhos = mat[overdispersionI,]
    predMat = expit(predMat)
    as = sweep(predMat, 2, 1/rhos-1, "*")
    bs = sweep(1-predMat, 2, 1/rhos-1, "*")
    predMatClustEffect = matrix(rbeta(length(predMat), c(as.matrix(as)), c(as.matrix(bs))), nrow=nrow(predMat))
    obsMat = expit(obsMat)
    as = sweep(obsMat, 2, 1/rhos-1, "*")
    bs = sweep(1-obsMat, 2, 1/rhos-1, "*")
    obsMatClustEffect = matrix(rbeta(length(obsMat), c(as.matrix(as)), c(as.matrix(bs))), nrow=nrow(obsMat))
    clusterVars = NULL
  } else {
    clusterVars = NULL
    rhos=NULL
    predMatClustEffect = predMat
    obsMatClustEffect = obsMat
  }
  
  # transform predictions from logit to probability scale
  if(family == "binomial") {
    predMat = expit(predMat)
    predMatClustEffect = expit(predMatClustEffect)
    obsMat = expit(obsMat)
    obsMatClustEffect = expit(obsMatClustEffect)
  }
  
  # compute predictive credible intervals
  if(is.null(predMatClustEffect)) {
    preds = rowMeans(predMat)
    predSDs = apply(predMat, 1, sd)
    lowerPreds = apply(predMat, 1, quantile, probs=(1-significanceCI)/2)
    medianPreds = apply(predMat, 1, median)
    upperPreds = apply(predMat, 1, quantile, probs=1-(1-significanceCI)/2)
    obsPreds = rowMeans(obsMat)
    obsSDs = apply(obsMat, 1, sd)
    lowerObs = apply(obsMat, 1, quantile, probs=(1-significanceCI)/2)
    medianObs = apply(obsMat, 1, median)
    upperObs = apply(obsMat, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    if(family == "normal") {
      preds = rowMeans(predMat)
      medianPreds = apply(predMat, 1, median)
      obsPreds = rowMeans(obsMat)
      medianObs = apply(obsMat, 1, median)
    } else {
      preds = rowMeans(predMatClustEffect)
      medianPreds = apply(predMatClustEffect, 1, median)
      obsPreds = rowMeans(obsMatClustEffect)
      medianObs = apply(obsMatClustEffect, 1, median)
    }
    predSDs = apply(predMatClustEffect, 1, sd)
    lowerPreds = apply(predMatClustEffect, 1, quantile, probs=(1-significanceCI)/2)
    upperPreds = apply(predMatClustEffect, 1, quantile, probs=1-(1-significanceCI)/2)
    obsSDs = apply(obsMatClustEffect, 1, sd)
    lowerObs = apply(obsMatClustEffect, 1, quantile, probs=(1-significanceCI)/2)
    upperObs = apply(obsMatClustEffect, 1, quantile, probs=1-(1-significanceCI)/2)
  }
  
  if(!is.null(xObs) && all(xObs[,1]==1))
    interceptSummary=mod$summary.fixed[1,c(1, 2, 4, 3, 5)]
  else
    interceptSummary = matrix(rep(0, 5), nrow=1)
  
  if(!is.null(xObs))
    fixedEffectSummary = mod$summary.fixed[,c(1, 2, 4, 3, 5)]
  else
    fixedEffectSummary = mod$summary.fixed
  
  # compute basis function coefficient predictions and standard deviations
  # startI = n+nrow(predCoords)+1
  # endI = n+nrow(predCoords)+nx*ny
  startI = 1
  endI = nx*ny
  thesePreds = rowMeans(basisCoefMat)
  theseSDs = apply(basisCoefMat, 1, sd)
  coefPreds = list(layer1 = thesePreds[startI:endI])
  coefSDs = list(layer1 = theseSDs[startI:endI])
  if(nLayer >= 2) {
    for(i in 2:nLayer) {
      startI = endI + 1
      endI = startI + nrow(latInfo[[i]]$latCoords) - 1
      # coefPreds = c(coefPreds, list(linpreds[startI:endI]))
      # coefSDs = c(coefSDs, list(linpred.sd[startI:endI]))
      coefPreds = c(coefPreds, list(thesePreds[startI:endI]))
      coefSDs = c(coefSDs, list(theseSDs[startI:endI]))
    }
  }
  
  endTime = proc.time()[3]
  totalTimeSampleProcessing = endTime - startTimeSampleProcessing
  print(paste0("finished processing samples. Took ", round(totalTimeSampleProcessing / 60, 2), " minutes"))
  totalTime = endTime - startTime
  timings = data.frame(totalTime=totalTime, 
                       precomputationTime=precomputationTime, 
                       modelDefineTime=totalTimeDefineModel, 
                       modelFitTime=totalTimeFitModel, 
                       posteriorSamplingTime=totalTimePosteriorSampling, 
                       sampleProcessingTime=totalTimeSampleProcessing, 
                       otherTime=totalTime-(precomputationTime + totalTimeDefineModel + totalTimeFitModel + totalTimePosteriorSampling + totalTimeSampleProcessing))
  timings$precomputationTimePct = timings$precomputationTime / timings$totalTime
  timings$modelDefinePct = timings$modelDefineTime / timings$totalTime
  timings$modelFitTimePct = timings$modelFitTime / timings$totalTime
  timings$posteriorSamplingTimePct = timings$posteriorSamplingTime / timings$totalTime
  timings$sampleProcessingTimePct = timings$sampleProcessingTime / timings$totalTime
  timings$otherTimePct = timings$otherTime / timings$totalTime
  
  list(preds=preds, sigmas=predSDs, lower=lowerPreds, median=medianPreds, upper=upperPreds, 
       obsPreds=obsPreds, obsSDs=obsSDs, obsLower=lowerObs, obsMedian=medianObs, obsUpper=upperObs, 
       mod=mod, latInfo=latInfo, coefPreds=coefPreds, coefSDs=coefSDs, 
       interceptSummary=interceptSummary, fixedEffectSummary=fixedEffectSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, overdispersionSummary=overdispersionSummary, parameterSummaryTable=parameterSummaryTable, 
       alphaSummary=alphaSummary, timings=timings, priorPar=priorPar, precomputedNormalizationFun=precomputedNormalizationFun, 
       # the rest of the outputs are saved to be used for spatial aggregations later on
       predMat=predMatClustEffect, obsMat=obsMatClustEffect, hyperMat=hyperMat, clusterVars=clusterVars, rhos=rhos
  )
}











# fit ELK model given data and relevant parameters
# obsCoords: matrix of spatial coordinates of the data
# obsValues: matrix of response at the observation coordinates
# predCoords: places at which to make predictions
# xObs: observation design matrix (intercept and covariates at obs locations)
# xPred: prediction design matrix (intercept and covariates at pred locations)
# NC: number of coarse lattice points over the longest dimension of the data
# nLayer: number of lattice layers
# seed: random seed set at beginning of function
# normalize: whether or not to normalize the process (highly recommended to set this to TRUE, 
#            since otherwise interpretation of correlation length and marginal variance priors 
#            will be highly suspect)
# intStrategy: see INLA:::control.inla()
# strategy: see INLA:::control.inla()
# predClusterI: whether or not to include cluster effects for which predictions
# priorPar default was getPrior(.1, .1, 10, nLayer=nLayer)
fitELK2 = function(obsCoords, obsValues, predCoords=obsCoords, seed=1, nLayer=3, NC=5,
                   nBuffer=5, priorPar=NULL, 
                   xObs=cbind(1, obsCoords), xPred=cbind(1, predCoords), normalize=TRUE, 
                   intStrategy="ccd", strategy="gaussian", fastNormalize=TRUE, 
                   predictionType=c("mean", "median"), significanceCI=0.8, 
                   printVerboseTimings=FALSE, nPostSamples=1000, family=c("normal", "binomial", "betabinomial"),
                   obsNs=rep(1, length(obsValues)), clusterEffect=TRUE, latInfo=NULL, 
                   initialEffectiveRange=NULL, initialAlphas=rep(1/nLayer, nLayer-1), 
                   effRangeRange=NULL, predClusterI=rep(TRUE, nrow(predCoords)), 
                   plotNormalizationSplines=FALSE, verbose=TRUE, separateRanges=FALSE, 
                   doValidation=FALSE, previousFit=NULL, precomputedNormalizationFun=NULL, 
                   savePrecomputationResults=FALSE, loadPrecomputationResults=FALSE, 
                   precomputationFileNameRoot="precomputationResults") {
  
  startTime = proc.time()[3]
  set.seed(seed)
  
  if(length(predCoords)< 1)
    stop("must include prediction locations")
  
  if(!normalize)
    warning("WARNING: normalize set to FALSE. Correlation length and marginal variance prior interpretation is suspect. It is highly recommended to set normalize to TRUE")
  
  # get the type of prediction the user wants
  predictionType = match.arg(predictionType)
  family = match.arg(family)
  
  # check to make sure inputs make sense
  if(!clusterEffect && family == "normal")
    stop("cluster effect must be included for the normal family")
  
  # get lattice points, prediction points
  
  if(is.null(latInfo)) {
    xRangeDat = range(obsCoords[,1])
    yRangeDat = range(obsCoords[,2])
    latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  } else {
    xRangeDat = latInfo[[1]]$xRangeDat
    yRangeDat = latInfo[[1]]$yRangeDat
  }
  nx = latInfo[[1]]$nx
  ny = latInfo[[1]]$ny
  
  if(is.null(priorPar))
    priorPar = getPCPrior(1444.772/5, .01, 1, nLayer=nLayer, separateRanges=separateRanges, latticeInfo=latInfo) # 1444.772/5
  
  # generate lattice basis matrix
  AObs = makeA(obsCoords, latInfo)
  
  ## run precomputations
  print("running precomputations...")
  startTimePrecomputations = proc.time()[3]
  
  # either run the precomputations or load them in
  if(!loadPrecomputationResults) {
    precomputedMatrices = precomputationsQ2(latInfo)
    if(is.null(precomputedNormalizationFun)) {
      precomputedNormalizationFun = precomputeNormalization(saveResults=FALSE, latticeInfo=latInfo, effRangeRange=effRangeRange, 
                                                            plotNormalizationSplines=plotNormalizationSplines)
    }
  } else {
    load(paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
  }
  
  # save the precomputations if necessary
  if(savePrecomputationResults) {
    save(precomputedMatrices, precomputedNormalizationFun, 
         file=paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
  }
  
  endTimePrecomputations = proc.time()[3]
  precomputationTime = endTimePrecomputations - startTimePrecomputations
  print(paste0("finished precomputations. Took ", round(precomputationTime / 60, 2), " minutes"))
  
  # define the LK model
  startTimeDefineModel = proc.time()[3]
  
  # set beta binomial prior if necessary
  if(family == "betabinomial") {
    if(clusterEffect)
      stop("cluster effect must not be set to TRUE for betabinomial model")
    
    # set median at .04 and upper 97.5th pctile at 0.2
    mu = logit(0.04)
    prec = 1/((logit(.2)-logit(.04))/qnorm(.975))^2
    control.family = list(hyper = list(rho = list(prior="logtnormal", param=c(mu, prec))))
  } else {
    control.family = list()
  }
  
  n = length(obsValues)
  if(separateRanges) {
    rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.full, latInfo=latInfo, ys=obsValues, 
                                prior=priorPar, normalize=normalize, precomputedMatrices=precomputedMatrices, 
                                X=xObs, datCoords=obsCoords, fastNormalize=fastNormalize, 
                                printVerboseTimings=printVerboseTimings, 
                                initialEffectiveRange=initialEffectiveRange, initialAlphas=initialAlphas, 
                                precomputedNormalizationFun=precomputedNormalizationFun, ns=obsNs)
  } else {
    rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.standard, latInfo=latInfo, ys=obsValues, 
                                prior=priorPar, normalize=normalize, precomputedMatrices=precomputedMatrices, 
                                X=xObs, datCoords=obsCoords, fastNormalize=fastNormalize, 
                                printVerboseTimings=printVerboseTimings, 
                                initialEffectiveRange=initialEffectiveRange, initialAlphas=initialAlphas, 
                                precomputedNormalizationFun=precomputedNormalizationFun, ns=obsNs)
  }
  
  # use these global variables for testing calls to inla.rgeneric.lk.model.simple
  # latInfo<<-latInfo; ys<<-obsValues; ns<<-obsNs;
  # prior<<-priorPar; normalize<<-normalize; precomputedMatrices<<-precomputedMatrices;
  # X<<-xObs; nu<<-nu; datCoords<<-obsCoords; fastNormalize<<-fastNormalize;
  # printVerboseTimings<<-printVerboseTimings; initialEffectiveRange<<-initialEffectiveRange;
  # initialAlphas<<-initialAlphas; precomputedNormalizationFun<<-precomputedNormalizationFun
  ## generate inla stack:
  # Stacked A matrix (A_s from notation of LR2015 Bayesian Spatial Modelling with R_INLA):
  # (AEst   0  )
  # ( 0   APred)
  # eta_s = (c^T c^T)^T
  # where c is the vector of lattice coefficients
  AEst = makeA(obsCoords, latInfo)
  APred = makeA(predCoords, latInfo)
  latticeInds = 1:ncol(AEst)
  clust = 1:length(obsValues)
  if(family == "normal") {
    if(!is.null(xObs)) {
      stack.est = inla.stack(A =list(AEst, 1), 
                             effects =list(field=latticeInds, X=xObs), 
                             data =list(y=obsValues, link=1), 
                             tag ="est", 
                             remove.unused=FALSE)
      stack.pred = inla.stack(A =list(matrix(APred[1,], nrow=1), 1), 
                              effects =list(field=latticeInds, X=matrix(xPred[1,], nrow=1)), 
                              data =list(y=NA, link=1), 
                              tag ="pred", 
                              remove.unused=FALSE)
    } else {
      stack.est = inla.stack(A =list(AEst), 
                             effects =list(field=latticeInds), 
                             data =list(y=obsValues, link=1), 
                             tag ="est", 
                             remove.unused=FALSE)
      stack.pred = inla.stack(A =list(matrix(APred[1,], nrow=1)), 
                              effects =list(field=latticeInds), 
                              data =list(y=NA, link=1), 
                              tag ="pred", 
                              remove.unused=FALSE)
    }
  }
  else if(family == "binomial" || family == "betabinomial") {
    if(clusterEffect) {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1, 1), 
                               effects =list(field=latticeInds, clust=clust, X=xObs), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      } else {
        stack.est = inla.stack(A =list(AEst, 1), 
                               effects =list(field=latticeInds, clust=clust), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      }
    } else {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1), 
                               effects =list(field=latticeInds, X=xObs), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      } else {
        stack.est = inla.stack(A =list(AEst), 
                               effects =list(field=latticeInds), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      }
    }
    if(!is.null(xObs)) {
      stack.pred = inla.stack(A =list(matrix(APred[1,], nrow=1), 1), 
                              effects =list(field=latticeInds, X=matrix(xPred[1,], nrow=1)), 
                              data =list(y=NA, link=1, Ntrials=1, 
                                         tag ="pred"), 
                              remove.unused=FALSE)
    } else {
      stack.pred = inla.stack(A =list(matrix(APred[1,], nrow=1)), 
                              effects =list(field=latticeInds), 
                              data =list(y=NA, link=1, Ntrials=1, 
                                         tag ="pred"), 
                              remove.unused=FALSE)
    }
  }
  stack.full = inla.stack(stack.est, stack.pred, 
                          remove.unused=FALSE)
  dat = inla.stack.data(stack.full, rgen=rgen, remove.unused=FALSE)
  
  # initialize the fitting process based on a previous optimum if necessary
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousFit)) {
    # modeControl$result = previousFit
    modeControl$theta = previousFit$mode$theta
    modeControl$x = previousFit$mode$x
    modeControl$restart = TRUE
  }
  endTimeDefineModel = proc.time()[3]
  totalTimeDefineModel = endTimeDefineModel - startTimeDefineModel
  
  # fit the model
  # control.inla = list(cmin = 0, int.strategy=int.strategy) 
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  controls = list(strategy=strategy, int.strategy=intStrategy) 
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  startTimeFitModel = proc.time()[3]
  if(family == "normal") {
    if(!is.null(xObs)) {
      mod = inla(y ~ - 1 + X + f(field, model=rgen), 
                 data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                 control.inla=controls, 
                 control.mode=modeControl, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                 control.fixed=list(quantiles=allQuantiles), 
                 control.family=list(hyper = list(prec = list(param=c(1, 0.05), prior="pc.prec"))))
    } else {
      mod = inla(y ~ - 1 + f(field, model=rgen), 
                 data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                 control.inla=controls, 
                 control.mode=modeControl, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                 control.fixed=list(quantiles=allQuantiles), 
                 control.family=list(hyper = list(prec = list(param=c(1, 0.05), prior="pc.prec"))))
    }
  } else if(family == "binomial" || family == "betabinomial") {
    if(clusterEffect) {
      # clusterList = list(param=c(.15, 0.01), prior="pc.prec")
      if(!is.null(xObs)) {
        mod = inla(y ~ - 1 + X + f(field, model=rgen) + 
                     f(clust, model="iid", hyper = list(prec = list(param=c(1, 0.01), prior="pc.prec"))), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      } else {
        mod = inla(y ~ - 1 + f(field, model=rgen) + 
                     f(clust, model="iid", hyper = list(prec = list(param=c(1, 0.01), prior="pc.prec"))), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      }
    } else {
      if(!is.null(xObs)) {
        mod = inla(y ~ - 1 + X + f(field, model=rgen), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      } else {
        mod = inla(y ~ - 1 + f(field, model=rgen), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE, quantiles=allQuantiles), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      }
    }
  }
  endTimeFitModel = proc.time()[3]
  totalTimeFitModel = endTimeFitModel - startTimeFitModel
  
  print(paste0("finished fitting model. Took ", round(totalTimeFitModel / 60, 2), " minutes"))
  
  # # improve the approximation of the posterior if requested by the user
  # if(improveHyperpar) {
  #   browser()
  #   mod = inla.hyperpar(mod)
  # }
  
  # browser()
  
  # get predictive surface, SD, and data
  index = inla.stack.index(stack.full, "pred")$data
  obsInds = 1:n
  predInds = (n+1):(n+nrow(predCoords))
  if(predictionType == "mean") {
    linpreds = mod[["summary.linear.predictor"]]$mean
    # preds = linpreds[predInds]
    # obsPreds = linpreds[obsInds]
  } else {
    linpreds = mod[["summary.linear.predictor"]]$`0.5quant`
    # preds = linpreds[predInds]
    # obsPreds = linpreds[obsInds]
  }
  linpred.sd = mod[["summary.linear.predictor"]]$sd
  # predSDs = linpred.sd[predInds]
  # obsSDs = linpred.sd[obsInds]
  
  # generate samples from posterior
  print("Sampling from posterior")
  startTimePosteriorSampling = proc.time()[3]
  postSamples = inla.posterior.sample(nPostSamples, mod)
  endTimePosteriorSampling = proc.time()[3]
  totalTimePosteriorSampling = endTimePosteriorSampling - startTimePosteriorSampling
  print(paste0("finished sampling from the posterior. Took ", round(totalTimePosteriorSampling / 60, 2), " minutes"))
  
  print("Processing posterior samples...")
  startTimeSampleProcessing = proc.time()[3]
  # get posterior hyperparameter samples and transform them as necessary
  # interpretation of hyperparameters: 
  # 1: error precision
  # 2: log effective range
  # 3: log spatial variance
  # 4-(3 + nLayer - 1): multivariateLogit alpha
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  if(family == "normal") {
    if(!separateRanges) {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[3])+1/x[1], spatialVar=exp(x[3]), errorVar=1/x[1], 
                                              totalSD=sqrt(exp(x[3])+1/x[1]), spatialSD=sqrt(exp(x[3])), errorSD=sqrt(1/x[1]), 
                                              spatialRange=exp(x[2]), alpha=multivariateExpit(x[4:(3 + nLayer - 1)]))})
      mat = rbind(mat, alpha=1-colSums(mat[8:(7+nLayer-1),]))
      hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", "spatialRange", 
                     paste0("alpha", 1:nLayer))
    } else {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[2+nLayer])+1/x[1], spatialVar=exp(x[2+nLayer]), errorVar=1/x[1], 
                                              totalSD=sqrt(exp(x[2+nLayer])+1/x[1]), spatialSD=sqrt(exp(x[2+nLayer])), errorSD=sqrt(1/x[1]), 
                                              spatialRange=exp(x[2:(1+nLayer)]), alpha=multivariateExpit(x[(3+nLayer):(3+2*nLayer-2)]))})
      mat = rbind(mat, alpha=1-colSums(matrix(mat[(6+nLayer+1):(6+nLayer+1 + nLayer-2),], nrow=nLayer-1)))
      hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", paste0("spatialRange", 1:nLayer), 
                     paste0("alpha", 1:nLayer))
    }
  } else if(family == "binomial") {
    if(clusterEffect) {
      if(!separateRanges) {
        logSpatialRangeI = 1
        logSpatialVarI = 2
        logitAlphaI = 3:(2 + nLayer - 1)
        precisionI = 3 + nLayer - 1
        mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[logSpatialVarI])+1/x[precisionI], spatialVar=exp(x[logSpatialVarI]), clusterVar=1/x[precisionI], 
                                                totalSD=sqrt(exp(x[logSpatialVarI])+1/x[precisionI]), spatialSD=sqrt(exp(x[logSpatialVarI])), clusterSD=sqrt(1/x[precisionI]), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        mat = rbind(mat, alpha=1-colSums(mat[8:(7+nLayer-1),]))
        hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", "spatialRange", 
                       paste0("alpha", 1:nLayer))
      } else {
        logSpatialRangeI = 1:nLayer
        logSpatialVarI = 1 + nLayer
        logitAlphaI = (logSpatialVarI+1):(logSpatialVarI+1 + nLayer-2)
        precisionI = max(logitAlphaI)+1
        mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[logSpatialVarI])+1/x[precisionI], spatialVar=exp(x[logSpatialVarI]), clusterVar=1/x[precisionI], 
                                                totalSD=sqrt(exp(x[logSpatialVarI])+1/x[precisionI]), spatialSD=sqrt(exp(x[logSpatialVarI])), clusterSD=sqrt(1/x[precisionI]), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        mat = rbind(mat, alpha=1-colSums(matrix(mat[(6+nLayer+1):(6+nLayer+1 + nLayer-2),], nrow=nLayer-1)))
        hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", paste0("spatialRange", 1:nLayer), 
                       paste0("alpha", 1:nLayer))
      }
    } else {
      if(!separateRanges) {
        logSpatialRangeI = 1
        logSpatialVarI = 2
        logitAlphaI = 3:(2 + nLayer - 1)
        mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        mat = rbind(mat, alpha=1-colSums(matrix(mat[4:(3+nLayer-1),], nrow=nLayer-1)))
        hyperNames = c("spatialVar", "spatialSD", "spatialRange", paste0("alpha", 1:nLayer))
      } else {
        logSpatialRangeI = 1:nLayer
        logSpatialVarI = 1 + nLayer
        logitAlphaI = (logSpatialVarI+1):(logSpatialVarI+1 + nLayer-2)
        mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        # mat = rbind(mat, alpha=1-colSums(mat[4:(3+nLayer-1),]))
        mat = rbind(mat, alpha=1-colSums(matrix(mat[(3+nLayer):(3+nLayer + nLayer-2),], nrow=nLayer-1)))
        hyperNames = c("spatialVar", "spatialSD", paste0("spatialRange", 1:nLayer), paste0("alpha", 1:nLayer))
      }
    }
  } else if(family == "betabinomial") {
    overdispersionI = 1
    if(!separateRanges) {
      logSpatialRangeI = 2
      logSpatialVarI = 3
      logitAlphaI = 4:(3 + nLayer - 1)
      mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                              spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]), 
                                              overdispersion=x[overdispersionI])})
      mat = rbind(mat[-nrow(mat),], alpha=1-colSums(matrix(mat[4:(3+nLayer-1),], nrow=nLayer-1)), mat[nrow(mat),])
      hyperNames = c("spatialVar", "spatialSD", "spatialRange", paste0("alpha", 1:nLayer), "overdispersion")
    } else {
      logSpatialRangeI = 2:(1+nLayer)
      logSpatialVarI = 2 + nLayer
      logitAlphaI = (logSpatialVarI+1):(logSpatialVarI+1 + nLayer-2)
      mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                              spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]), 
                                              overdispersion=x[overdispersionI])})
      # mat = rbind(mat, alpha=1-colSums(mat[4:(3+nLayer-1),]))
      mat = rbind(mat[-nrow(mat),], alpha=1-colSums(matrix(mat[(3+nLayer):(3+nLayer + nLayer-2),], nrow=nLayer-1)), mat[nrow(mat),])
      hyperNames = c("spatialVar", "spatialSD", paste0("spatialRange", 1:nLayer), paste0("alpha", 1:nLayer), "overdispersion")
    }
  }
  
  rownames(mat) = hyperNames
  
  getSummaryStatistics = function(draws) {
    c(Est=mean(draws), SD=sd(draws), 
      Qlower=quantile(probs=(1 - significanceCI) / 2, draws), 
      Q50=quantile(probs=0.5, draws), 
      Qupper=quantile(probs=1 - (1 - significanceCI) / 2, draws))
  }
  summaryNames = c("Est", "SD", "Qlower", "Q50", "Qupper")
  parameterSummaryTable = t(apply(mat, 1, getSummaryStatistics))
  colnames(parameterSummaryTable) = summaryNames
  
  # separate out default parameter summaries
  varI = which(grepl("clusterVar", hyperNames))
  rangeI = which(grepl("Range", hyperNames))
  alphaI = which(grepl("alpha", hyperNames))
  varSummary=parameterSummaryTable[varI,]
  rangeSummary=parameterSummaryTable[rangeI,]
  alphaSummary=parameterSummaryTable[alphaI,]
  if(family == "normal" || clusterEffect) {
    sdI = which(grepl("clusterSD", hyperNames))
    sdSummary=parameterSummaryTable[sdI,]
  } else {
    sdSummary = matrix(rep(0, 5), nrow=1)
  }
  overdispersionSummary = matrix(rep(0, 5), nrow=1)
  if(family == "betabinomial")
    overdispersionSummary=parameterSummaryTable[nrow(parameterSummaryTable),]
  
  # get samples of the latent field
  latentMat = sapply(postSamples, function(x) {x$latent})
  latentVarNames = rownames(postSamples[[1]]$latent)
  fieldIndices = which(grepl("field", latentVarNames))
  fixedIndices = which(grepl("X", latentVarNames))
  # if(clusterEffect)
  #   clustIndices = grepl("clust", latentVarNames)
  
  ## generate logit predictions (first without cluster effect then add the cluster effect in)
  # for prediction locations
  if(length(xPred) != 0)
    fixedPart = xPred  %*% latentMat[fixedIndices,]
  else
    fixedPart = 0
  predMat = fixedPart + APred %*% latentMat[fieldIndices,]
  
  # for observation locations
  if(length(xObs) != 0)
    fixedPart = xObs  %*% latentMat[fixedIndices,]
  else
    fixedPart = 0
  obsMat = fixedPart + AObs %*% latentMat[fieldIndices,]
  
  # get draws from basis function coefficients
  basisCoefMat = latentMat[fieldIndices,]
  
  # add in cluster effect if necessary
  if((family == "binomial" && clusterEffect) || family == "normal") {
    # get betabinomial overdispersion parameter
    clusterVarI = which(grepl("clusterVar", hyperNames))
    clusterVars = mat[clusterVarI,]
    predMatClustEffect = predMat + sweep(matrix(rnorm(length(predMat), sd=rep(sqrt(clusterVars), each=nrow(predMat))), nrow=nrow(predMat)), 1, predClusterI, "*")
    obsMatClustEffect = obsMat + matrix(rnorm(length(obsMat), sd=rep(sqrt(clusterVars), each=nrow(obsMat))), nrow=nrow(obsMat))
    rhos = NULL
  } else if(family == "betabinomial") {
    # get cluster induced overdispersion
    overdispersionI = which(grepl("overdispersion", hyperNames))
    rhos = mat[overdispersionI,]
    predMat = expit(predMat)
    as = sweep(predMat, 2, 1/rhos-1, "*")
    bs = sweep(1-predMat, 2, 1/rhos-1, "*")
    predMatClustEffect = matrix(rbeta(length(predMat), c(as.matrix(as)), c(as.matrix(bs))), nrow=nrow(predMat))
    obsMat = expit(obsMat)
    as = sweep(obsMat, 2, 1/rhos-1, "*")
    bs = sweep(1-obsMat, 2, 1/rhos-1, "*")
    obsMatClustEffect = matrix(rbeta(length(obsMat), c(as.matrix(as)), c(as.matrix(bs))), nrow=nrow(obsMat))
    clusterVars = NULL
  } else {
    clusterVars = NULL
    rhos=NULL
    predMatClustEffect = predMat
    obsMatClustEffect = obsMat
  }
  
  # transform predictions from logit to probability scale
  if(family == "binomial") {
    predMat = expit(predMat)
    predMatClustEffect = expit(predMatClustEffect)
    obsMat = expit(obsMat)
    obsMatClustEffect = expit(obsMatClustEffect)
  }
  
  # compute predictive credible intervals
  if(is.null(predMatClustEffect)) {
    preds = rowMeans(predMat)
    predSDs = apply(predMat, 1, sd)
    lowerPreds = apply(predMat, 1, quantile, probs=(1-significanceCI)/2)
    medianPreds = apply(predMat, 1, median)
    upperPreds = apply(predMat, 1, quantile, probs=1-(1-significanceCI)/2)
    obsPreds = rowMeans(obsMat)
    obsSDs = apply(obsMat, 1, sd)
    lowerObs = apply(obsMat, 1, quantile, probs=(1-significanceCI)/2)
    medianObs = apply(obsMat, 1, median)
    upperObs = apply(obsMat, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    if(family == "normal") {
      preds = rowMeans(predMat)
      medianPreds = apply(predMat, 1, median)
      obsPreds = rowMeans(obsMat)
      medianObs = apply(obsMat, 1, median)
    } else {
      preds = rowMeans(predMatClustEffect)
      medianPreds = apply(predMatClustEffect, 1, median)
      obsPreds = rowMeans(obsMatClustEffect)
      medianObs = apply(obsMatClustEffect, 1, median)
    }
    predSDs = apply(predMatClustEffect, 1, sd)
    lowerPreds = apply(predMatClustEffect, 1, quantile, probs=(1-significanceCI)/2)
    upperPreds = apply(predMatClustEffect, 1, quantile, probs=1-(1-significanceCI)/2)
    obsSDs = apply(obsMatClustEffect, 1, sd)
    lowerObs = apply(obsMatClustEffect, 1, quantile, probs=(1-significanceCI)/2)
    upperObs = apply(obsMatClustEffect, 1, quantile, probs=1-(1-significanceCI)/2)
  }
  
  if(!is.null(xObs) && all(xObs[,1]==1))
    interceptSummary=mod$summary.fixed[1,c(1, 2, 4, 3, 5)]
  else
    interceptSummary = matrix(rep(0, 5), nrow=1)
  
  if(!is.null(xObs))
    fixedEffectSummary = mod$summary.fixed[,c(1, 2, 4, 3, 5)]
  else
    fixedEffectSummary = mod$summary.fixed
  
  # compute basis function coefficient predictions and standard deviations
  # startI = n+nrow(predCoords)+1
  # endI = n+nrow(predCoords)+nx*ny
  startI = 1
  endI = nx*ny
  thesePreds = rowMeans(basisCoefMat)
  theseSDs = apply(basisCoefMat, 1, sd)
  coefPreds = list(layer1 = thesePreds[startI:endI])
  coefSDs = list(layer1 = theseSDs[startI:endI])
  if(nLayer >= 2) {
    for(i in 2:nLayer) {
      startI = endI + 1
      endI = startI + nrow(latInfo[[i]]$latCoords) - 1
      # coefPreds = c(coefPreds, list(linpreds[startI:endI]))
      # coefSDs = c(coefSDs, list(linpred.sd[startI:endI]))
      coefPreds = c(coefPreds, list(thesePreds[startI:endI]))
      coefSDs = c(coefSDs, list(theseSDs[startI:endI]))
    }
  }
  
  endTime = proc.time()[3]
  totalTimeSampleProcessing = endTime - startTimeSampleProcessing
  print(paste0("finished processing samples. Took ", round(totalTimeSampleProcessing / 60, 2), " minutes"))
  totalTime = endTime - startTime
  timings = data.frame(totalTime=totalTime, 
                       precomputationTime=precomputationTime, 
                       modelDefineTime=totalTimeDefineModel, 
                       modelFitTime=totalTimeFitModel, 
                       posteriorSamplingTime=totalTimePosteriorSampling, 
                       sampleProcessingTime=totalTimeSampleProcessing, 
                       otherTime=totalTime-(precomputationTime + totalTimeDefineModel + totalTimeFitModel + totalTimePosteriorSampling + totalTimeSampleProcessing))
  timings$precomputationTimePct = timings$precomputationTime / timings$totalTime
  timings$modelDefinePct = timings$modelDefineTime / timings$totalTime
  timings$modelFitTimePct = timings$modelFitTime / timings$totalTime
  timings$posteriorSamplingTimePct = timings$posteriorSamplingTime / timings$totalTime
  timings$sampleProcessingTimePct = timings$sampleProcessingTime / timings$totalTime
  timings$otherTimePct = timings$otherTime / timings$totalTime
  
  list(preds=preds, sigmas=predSDs, lower=lowerPreds, median=medianPreds, upper=upperPreds, 
       obsPreds=obsPreds, obsSDs=obsSDs, obsLower=lowerObs, obsMedian=medianObs, obsUpper=upperObs, 
       mod=mod, latInfo=latInfo, coefPreds=coefPreds, coefSDs=coefSDs, 
       interceptSummary=interceptSummary, fixedEffectSummary=fixedEffectSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, overdispersionSummary=overdispersionSummary, parameterSummaryTable=parameterSummaryTable, 
       alphaSummary=alphaSummary, timings=timings, priorPar=priorPar, precomputedNormalizationFun=precomputedNormalizationFun, 
       # the rest of the outputs are saved to be used for spatial aggregations later on
       predMat=predMatClustEffect, obsMat=obsMatClustEffect, hyperMat=hyperMat, clusterVars=clusterVars, rhos=rhos
  )
}