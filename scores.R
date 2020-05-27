# this script contains all used scoring rules used to evaluate the models for a given data set

# this function computes all scoring rules
# truth: the true values
# est: the estimates, of the same length as truth. By default calculated from estMat
# var: the estimates, of the same length as truth. By default calculated from estMat
# lower: the lower end of the credible interval. By default calculated from estMat
# upper: the upper end of the credible interval. By default calculated from estMat
# estMat: a matrix of joint estimate draws, with number of rows equal to the length of truth, a number of 
#          columns equal to the number of draws. If not included, a gaussian distribution is assumed.
# significance: the significance level of the credible interval. By default 80%
# distances: the distances to the nearest observation if not NULL, scores are broken up 
#            as a function of nearest neighbor distances
# breaks: the number of equal spaced bins to break the scores into as a function of distance
# NOTE: Discrete, count level credible intervals are estimated based on the input estMat along with coverage and CRPS
getScores = function(truth, est=NULL, var=NULL, lower=NULL, upper=NULL, estMat=NULL, significance=.8, 
                     distances=NULL, breaks=30, doRandomReject=FALSE, doFuzzyReject=TRUE, getAverage=TRUE) {
  
  # if distances is included, must also break down scoring rules by distance bins
  if(!is.null(distances)) {
    # construct the distance bins with which to group the data and compute scores within
    if(length(breaks) == 1)
      breaks = seq(0, max(distances), l=breaks)
    binsI = cut(distances, breaks, labels=1:(length(breaks)-1), include.lowest=TRUE)
    centers = breaks[1:(length(breaks)-1)] + diff(breaks)/2
    uniqueBinsI = sort(unique(binsI))
    
    # determine the number of observations per bin
    nPerBin = as.numeric(table(binsI))
    
    # helper function to compute the scoring rules for a given bin
    getSubScores = function(uniqueBinI, truth, est, var, lower, upper, estMat, significance) {
      thisDatI = binsI == uniqueBinI
      
      newEstMat = NULL
      if(!is.null(estMat))
        newEstMat = matrix(estMat[thisDatI,], ncol=ncol(estMat))
      getScores(truth[thisDatI], est[thisDatI], var[thisDatI], lower[thisDatI], upper[thisDatI], 
                newEstMat, significance, doRandomReject=doRandomReject, doFuzzyReject=doFuzzyReject)
    }
    
    # calculate scores for each bin individually
    binnedScores = t(sapply(uniqueBinsI, getSubScores, truth=truth, est=est, var=var, lower=lower, upper=upper, 
                          estMat=estMat, significance=significance))
    
    # make sure each variable in binnedScores is a numeric, not a list...
    temp = matrix(unlist(binnedScores), nrow=length(uniqueBinsI))
    theseNames = colnames(binnedScores)
    binnedScores = data.frame(temp)
    names(binnedScores) = theseNames
    binnedScores = as.data.frame(cbind(NNDist=centers[uniqueBinsI], nPerBin=nPerBin[uniqueBinsI], binnedScores))
  }
  
  # compute central estimates if estMat is not null
  if(!is.null(estMat)) {
    if(is.null(est))
      est = rowMeans(estMat)
  }
  
  # first calculate bias, variance, and MSE
  out = mse(truth, est, getAverage=getAverage)
  thisMSE = out$MSE
  thisBias = out$bias
  thisVar = out$var
  
  # calculate coverage and credible interval width with and without binomial variation
  coverage = coverage(truth, est, var, lower, upper, estMat=estMat, 
                      significance=significance, returnIntervalWidth=TRUE, 
                      doRandomReject=doRandomReject, doFuzzyReject=doFuzzyReject, getAverage=getAverage)
  if(getAverage) {
    thisCoverage = coverage[1]
    thisWidth = coverage[2]
  } else {
    thisCoverage = coverage[,1]
    thisWidth = coverage[,2]
  }
  
  # calculate CRPS
  thisCRPS = crps(truth, est, var, estMat=estMat, getAverage=getAverage)
  
  # collect the results in a data frame
  results = matrix(c(thisBias, thisVar, thisMSE, sqrt(thisMSE), thisCRPS, thisCoverage, 
                     thisWidth), ncol=7)
  colnames(results) = c("Bias", "Var", "MSE", "RMSE", "CRPS", "Coverage", "Width")
  results = as.data.frame(results)
  
  # include both binned and pooled results in one final table
  if(!is.null(distances)) {
    results = list(pooledResults=results, binnedResults=binnedScores)
  }
  
  results
}

# calculate bias, variance, and MSE
mse <- function(truth, est, weights=NULL, getAverage=TRUE){
  if(!is.null(weights))
    weights = weights / sum(weights)
  
  res = est - truth
  
  if(!is.null(weights)) {
    thisVar = (res - sum(res*weights))^2
    if(getAverage) {
      MSE = sum(res^2 * weights)
      bias=sum(res * weights)
      thisVar = sum(thisVar * weights)
    } else {
      MSE = res^2
      bias = res
    }
    
    out = list(MSE=MSE, bias=bias, var=thisVar)
  }
  else {
    thisVar = (res - mean(res))^2
    if(getAverage) {
      MSE = mean(res^2)
      bias=mean(res)
      thisVar = mean(thisVar)
    } else {
      MSE = res^2
      bias=res
    }
    
    out = list(MSE=MSE, bias=bias, var=thisVar)
  }
  
  out
}

# either include both lower and upper, or include either: 
#    - the joint estimate draw matrix
#    - estimates and variances (assumes gaussian)
# truth: the true empirical proportions of mortality rates within the regions or enumeration areas of interest
# lower: the lower end of the credible interval
# upper: the upper end of the credible interval
# estMat: a matrix of joint draws of estimates, with number of rows equal to the length of truth, a number of 
#         columns equal to the number of draws. If not included, a lgaussian distribution is assumed. Can be 
#         Gaussian or discrete values such as empirical proportions
# significance: the significance level of the credible interval. By default 80%
# doRandomReject, doFuzzyReject: based on https://www.jstor.org/stable/pdf/20061193.pdf
# ns: a vector of maximum possible counts (denominators) for each observation. Used only for random/fuzzy reject. 
#     Can be left out, in which case it will be inferred from the minimum draw difference in each row of estMat.
coverage = function(truth, est=NULL, var=NULL, lower=NULL, upper=NULL, 
                    estMat=NULL, significance=.8, returnIntervalWidth=FALSE, 
                    doRandomReject=FALSE, doFuzzyReject=TRUE, getAverage=TRUE, ns=NULL){
  
  if(any(is.null(lower)) || any(is.null(upper))) {
    # if the user did not supply their own credible intervals, we must get them ourselves given the other information
    
    if(is.null(estMat) && (is.null(est) || is.null(var)))
    stop("either include both lower and upper, est and var, or estMat")
    
    if(!is.null(est) && !is.null(var) && is.null(estMat)) {
      # in this case, we must calculate lower and upper assuming gaussianity
      lower = qnorm((1 - significance) / 2, est, sqrt(var))
      upper = qnorm(1 - (1 - significance) / 2, est, sqrt(var))
    }
    else {
      # we don't have information about the predictive distribution, and don't assume normality. 
      # Instead, use the user supplied to probability matrix estMat
      
      # take the quantiles of the probability draws
      CIs = apply(estMat, 1, function(ps) {quantile(ps, probs=c((1 - significance) / 2, 1 - (1 - significance) / 2))})
      lower = CIs[1,]
      upper = CIs[2,]
    }
  }
  
  if(any(lower > upper)) {
    warning("lower > upper, reordering")
    tmp = lower
    wrongOrder = lower > upper
    lower[wrongOrder] = upper[wrongOrder]
    upper[wrongOrder] = tmp[wrongOrder]
  }
  
  res = lower <= truth & upper >= truth
  
  if(returnIntervalWidth)
    width = upper - lower
  
  if(doRandomReject) {
    # in this case, we sometimes randomly reject if the truth is at the edge of the coverage interval. First 
    # determine what values are at the edge of the intervals, then determine the probability of rejection 
    # for each, then randomly reject
    atLowerEdge = which(lower == truth)
    atUpperEdge = which(upper == truth)
    
    probRejectLower = sapply(atLowerEdge, function(i) {((1 - significance) / 2 - mean(estMat[i,] < lower[i])) / mean(estMat[i,] == lower[i])})
    probRejectUpper = sapply(atUpperEdge, function(i) {((1 - significance) / 2 - mean(estMat[i,] > upper[i])) / mean(estMat[i,] == upper[i])})
    
    if(doFuzzyReject) {
      rejectLower = probRejectLower
      rejectUpper = probRejectUpper
    } else {
      rejectLower = runif(length(atLowerEdge)) <= probRejectLower
      rejectUpper = runif(length(atUpperEdge)) <= probRejectUpper
    }
    
    # determine minimum differences between probabilities
    if(is.null(ns))
      deltas = apply(estMat, 1, function(x) {min(diff(sort(unique(x))))})
    else
      deltas = 1 / ns
    
    if(length(atLowerEdge) != 0) {
      res[atLowerEdge] = sapply(1:length(atLowerEdge), function(i) {min(res[atLowerEdge][i], (1-rejectLower[i]))})
      
      # if reject, reduce CI width
      width[atLowerEdge] = width[atLowerEdge] - deltas[atLowerEdge] * as.numeric(rejectLower)
    }
    if(length(atUpperEdge) != 0) {
      res[atUpperEdge] = sapply(1:length(atUpperEdge), function(i) {min(res[atUpperEdge][i], (1-rejectUpper[i]))})
      
      # if reject, reduce CI width
      width[atUpperEdge] = width[atUpperEdge] - deltas[atUpperEdge] * as.numeric(rejectUpper)
    }
    # res[atLowerEdge] = res[atLowerEdge] & (!rejectLower)
    # res[atUpperEdge] = res[atUpperEdge] & (!rejectUpper)
  }
  
  if(getAverage)
    allResults = c(coverage=mean(res))
  else
    allResults = c(coverage=res)
  
  if(returnIntervalWidth) {
    if(getAverage)
      allResults = c(allResults, width=mean(width))
    else
      allResults = cbind(allResults, width=width)
  }
  
  allResults
}

# truth: a vector of observations on the desired scale
# est: a vector of logit-scale predictions of the same length as truth 
# my.var: a vector of logit-scale predictive variances of the same length as truth
# estMat: if available, use these probability draws in the integration. Use this argument 
#         when a gaussian approximation to the (possibly transformed) posterior is unreasonable
crps <- function(truth, est=NULL, my.var=NULL, estMat=NULL, getAverage=TRUE){
  if(!is.null(est) && !is.null(my.var) && is.null(estMat)) {
    sig = sqrt(my.var)
    x0 <- (truth - est) / sig
    res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))
    
    ## sign as in Held (2008)
    res <- -res
  }
  else {
    # Integrate numerically using estMat
    if(is.null(estMat))
      stop("must include either or both est and my.var, or estMat")
    
    # the following code was commented out since it has been modified to be computationally efficient
    # # compute the crps for this row of truth
    # crpsRow = function(rowI) {
    #   thisTruth = truth[rowI]
    #   
    #   # either build the predictive cdf assuming normality on the logit scale or from the 
    #   # empirical distribution given by estMat if the user supplies it
    #   if(is.null(estMat)) {
    #     thisEst = est[rowI]
    #     thisVar = my.var[rowI]
    #     thisCdf = function(ws) {pnorm(logit(ws), thisEst, sqrt(thisVar))}
    #   } else {
    #     thisCdf = ecdf(estMat[rowI,])
    #   }
    #   
    #   intFun = function(ws) {
    #     (thisCdf(ws) - (ws >= thisTruth))^2
    #   }
    #   
    #   if(is.null(estMat)) {
    #     # when integrating we will set bounds on the integral to be reasonable to avoid 
    #     # faulty "divergent integral" error. The bounds will be 20 standard errors out 
    #     # of the estimate, making sure to include the truth.
    #     lowerBound = max(0, min(thisTruth - .01, expit(thisEst - 20 * sqrt(thisVar))))
    #     upperBound = min(1, max(thisTruth + .01, expit(thisEst + 20 * sqrt(thisVar))))
    #     integrate(intFun, lowerBound, upperBound)$value
    #   }
    #   else {
    #     # since we are using the empirical distribution, there is a closed form for the integral
    #     ps = estMat[rowI,]
    #     allPoints = sort(c(ps, thisTruth))
    #     deltas = diff(allPoints)
    #     sum(deltas * intFun(allPoints[1:length(ps)]))
    #   }
    # }
    # 
    # crpsRow2 = function(rowI) {
    #   thisTruth = truth[rowI]
    #   
    #   # either build the predictive cdf assuming normality on the logit scale or from the 
    #   # empirical distribution given by estMat if the user supplies it
    #   if(is.null(estMat)) {
    #     thisEst = est[rowI]
    #     thisVar = my.var[rowI]
    #     thisCdf = function(ws) {pnorm(logit(ws), thisEst, sqrt(thisVar))}
    #   } else {
    #     # thisCdf = ecdf(estMat[rowI,])
    #     sorted = sort(estMat[rowI,])
    #     thisCdf = approxfun(sorted, (1:length(sorted))/length(sorted), 
    #                         method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    #   }
    #   
    #   intFun = function(ws) {
    #     (thisCdf(ws) - (ws >= thisTruth))^2
    #   }
    #   
    #   if(is.null(estMat)) {
    #     # when integrating we will set bounds on the integral to be reasonable to avoid 
    #     # faulty "divergent integral" error. The bounds will be 20 standard errors out 
    #     # of the estimate, making sure to include the truth.
    #     lowerBound = max(0, min(thisTruth - .01, expit(thisEst - 20 * sqrt(thisVar))))
    #     upperBound = min(1, max(thisTruth + .01, expit(thisEst + 20 * sqrt(thisVar))))
    #     integrate(intFun, lowerBound, upperBound)$value
    #   }
    #   else {
    #     # since we are using the empirical distribution, there is a closed form for the integral
    #     allPoints = sort(c(sorted, thisTruth))
    #     firstGreater = match(TRUE, sorted >= thisTruth)
    #     if(is.na(firstGreater))
    #       allPoints = c(sorted, thisTruth)
    #     else if(firstGreater == 1)
    #       allPoints = c(thisTruth, sorted)
    #     else
    #       allPoints = c(sorted[1:(firstGreater - 1)], thisTruth, sorted[firstGreater:length(sorted)])
    #     
    #     deltas = diff(allPoints)
    #     sum(deltas * intFun(allPoints[1:length(sorted)]))
    #   }
    # }
    
    crpsRow = function(rowI) {
      thisTruth = truth[rowI]
      
      # build the predictive cdf assuming from the empirical distribution given by 
      # estMat
      
      # thisCdf = ecdf(estMat[rowI,])
      sorted = estMat[rowI,] # already sorted
      
      # since we are using the empirical distribution, there is a closed form for the integral
      allPoints = sort(c(sorted, thisTruth))
      deltas = diff(allPoints)
      firstGreater = match(TRUE, sorted >= thisTruth)
      vals = (1:length(sorted))/length(sorted)
      if(is.na(firstGreater))
        return(sum((vals)^2 * deltas))
      else if(firstGreater == 1)
        return(deltas[1] + sum((1-vals[1:(length(sorted)-1)])^2 * deltas[2:length(deltas)]))
      else {
        left = sum(vals[1:(firstGreater-1)]^2 * deltas[1:(firstGreater-1)])
        mid = sum((1 - vals[firstGreater-1])^2 * deltas[firstGreater])
        right = ifelse(firstGreater == length(vals), 0, sum((1 - vals[firstGreater:(length(vals)-1)])^2 * deltas[(firstGreater+1):length(deltas)]))
        return(left+mid+right)
      }
      
      # intFun = function(ws) {
      #   (thisCdf(ws) - (ws >= thisTruth))^2
      # }
    }
    
    if(!is.null(estMat))
      estMat = t(apply(estMat, 1, sort))
    res = sapply(1:length(truth), crpsRow)
  }
  
  if(getAverage)
    mean(res)
  else
    res
}

# averages a list of many tables, each returned from the getScores function with distanceBreaks set by user
averageBinnedScores = function(tableList) {
  
  if(length(tableList) == 1) {
    # base case
    
    return(as.data.frame(tableList[[1]]))
  } else if(is.null(tableList[[2]])) {
    # minor case (some of the tables might be NULL)
    
    return(averageBinnedScores(tableList[-2]))
  } else {
    # recursive case
    
    firstTable = tableList[[1]]
    secondTable = tableList[[2]]
    
    # make sure all tables are matrices
    firstTable = as.matrix(firstTable)
    secondTable = as.matrix(secondTable)
    
    # make sure tables have matched bins
    uniqueDists = sort(unique(c(firstTable[,1], secondTable[,1])))
    firstMatch = match(firstTable[,1], uniqueDists)
    secondMatch = match(secondTable[,1], uniqueDists)
    newFirstTable = matrix(0, nrow=length(uniqueDists), ncol=ncol(firstTable))
    newFirstTable[,1] = uniqueDists
    newFirstTable[firstMatch,] = firstTable
    colnames(newFirstTable) = colnames(firstTable)
    newSecondTable = matrix(0, nrow=length(uniqueDists), ncol=ncol(secondTable))
    newSecondTable[,1] = uniqueDists
    newSecondTable[secondMatch,] = secondTable
    colnames(newSecondTable) = colnames(secondTable)
    
    firstTable = newFirstTable
    secondTable = newSecondTable
    
    # calculate weights for averaging
    # ns1 = firstTable$nPerBin
    # ns2 = secondTable$nPerBin
    ns1 = firstTable[,2]
    ns2 = secondTable[,2]
    nsTotal = ns1 + ns2
    ws1 = ns1 / nsTotal
    ws2 = ns2 / nsTotal
    
    # perform weighted averaging
    newTable = firstTable
    newTable[,2] = ns1 + ns2
    newTable[,3:ncol(firstTable)] = sweep(firstTable[,3:ncol(firstTable)], 1, ws1, "*") + sweep(secondTable[,3:ncol(secondTable)], 1, ws2, "*")
    
    # return results recursively
    return(averageBinnedScores(c(list(newTable), tableList[-(1:2)])))
  }
}

# singleScores: data frame with area names or observation IDs and "NNDist" variable in the first two columns 
#               along the with other columns for other scoring rules
aggregateScoresByDistance = function(singleScores, breaks=30, dat, distanceBreaksType=c("quantiles", "even"), nPerBin=NULL, maxDist=Inf) {
  distanceBreaksType = match.arg(distanceBreaksType)
  
  # Now determine distance to which type of point we use
  distanceVar = "NNDist"
  distances = singleScores[[distanceVar]]
  
  # remove distances beyond the maximum of breaks
  if(length(breaks) != 1) {
    maxDist = min(maxDist, max(breaks))
  }
  badDistances = distances >= maxDist
  singleScores = singleScores[!badDistances,]
  distances = distances[!badDistances]
  
  # sort table by distances
  sortI = sort(distances, index.return=TRUE)$ix
  singleScores = singleScores[sortI,]
  distances = distances[sortI]
  
  # calculate default breaks for the bin limits if necessary
  if(length(breaks) == 1) {
    nBreaks = breaks
    if(distanceBreaksType == "even" && is.null(nPerBin)) {
      breaks = seq(0, max(distances), l=nBreaks)
    } else {
      if(is.null(nPerBin))
        nPerBin = ceiling(nrow(singleScores)/nBreaks)
      
      # get endpoints of the bins, average their values when calculating breaks
      startI = seq(nPerBin+1, nrow(singleScores), by=nPerBin)
      endI = startI - 1
      breaks = c(0, c(rowMeans(cbind(distances[startI], distances[endI]))), distances[length(distances)]+1e-6)
    }
  }
  
  # construct the distance bins with which to group the data and compute scores within
  binsI = cut(distances, breaks, labels=1:(length(breaks)-1), include.lowest=TRUE)
  centers = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  uniqueBinsI = sort(unique(binsI))
  
  # determine the number of observations per bin
  nPerBin = as.numeric(table(binsI))
  
  # helper function to compute the scoring rules for a given bin
  getSubScores = function(uniqueBinI) {
    thisDatI = binsI == uniqueBinI
    
    # thisSingleScoresBinomial = data.frame(c(list(Region=thisRegion, dataI=which(thisSampleI), NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(truth, est, vars, lower, upper, estMatBinomial, getAverage=FALSE), Time=time[3]))
    colMeans(singleScores[thisDatI,-c(1, 2)])
  }
  
  # calculate scores for each bin individually
  binnedScores = t(sapply(uniqueBinsI, getSubScores))
  
  # make sure each variable in binnedScores is a numeric, not a list...
  temp = matrix(unlist(binnedScores), nrow=length(uniqueBinsI))
  theseNames = colnames(binnedScores)
  binnedScores = data.frame(temp)
  names(binnedScores) = theseNames
  out = as.data.frame(cbind(nPerBin=nPerBin[uniqueBinsI], binnedScores))
  out[[distanceVar]] = centers[uniqueBinsI]
  out
}



















