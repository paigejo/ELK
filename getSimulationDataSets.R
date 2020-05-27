# this script is for generating the data sets for the simulation study

# generate the data sets that will be used in the simulation study
getSimulationDataSets = function(nTotal=1000, nTest=100, marginalVar=1, errorVar=sqrt(.1), nDataSets=100, seed=123) {
  if(is.null(seed))
    set.seed(seed)
  
  # set the correlation functions for the simulation study
  exponentialCorFuns = list(function(x) {Exp.cov(x, theta=0.1)}, 
                            function(x) {Exp.cov(x, theta=0.5)}, 
                            function(x) {Exp.cov(x, theta=1)})
  maternCorFuns = list(function(x) {stationary.cov(x, theta=0.1, Covariance="Matern", nu=1)}, 
                       function(x) {stationary.cov(x, theta=0.5, Covariance="Matern", nu=1)}, 
                       function(x) {stationary.cov(x, theta=1, Covariance="Matern", nu=1)})
  # mixtureCorFun = function(x) {0.4 * Exp.cov(x, theta=0.1) + 0.6 * Exp.cov(x, theta=3)}
  mixtureCorFun = function(x) {0.4 * Exp.cov(x, theta=0.1) + 0.6 * Exp.cov(x, theta=.4)}
  
  ##### simulate the corresponding data sets
  ## simulate the exponential correlation function data sets
  print("Simulating dataset with first exponential correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(exponentialCorFuns[[1]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Exp(0.1))", fileNameRoot="exp01")
  print("Saving dataset...")
  save(simulationData, file="exponential01DataSet.RData")
  
  print("Simulating data set with second exponential correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(exponentialCorFuns[[2]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Exp(0.5))", fileNameRoot="exp05")
  print("Saving dataset...")
  save(simulationData, file="exponential05DataSet.RData")
  
  print("Simulating data set with third exponential correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(exponentialCorFuns[[3]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Exp(1))", fileNameRoot="exp1")
  print("Saving dataset...")
  save(simulationData, file="exponential1DataSet.RData")
  
  ## simulate the Matern correlation function data sets
  print("Simulating dataset with first matern correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(maternCorFuns[[1]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Matern(1,0.1))", fileNameRoot="mat01")
  print("Saving dataset...")
  save(simulationData, file="matern01DataSet.RData")
  
  print("Simulating data set with second matern correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(maternCorFuns[[2]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Matern(1,.5))", fileNameRoot="mat05")
  print("Saving dataset...")
  save(simulationData, file="matern05DataSet.RData")
  
  print("Simulating data set with third matern correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(maternCorFuns[[3]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Matern(1,1))", fileNameRoot="mat1")
  print("Saving dataset...")
  save(simulationData, file="matern1DataSet.RData")
  
  # simulate the mixture correlation function data set
  print("Simulating data set with mixture correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(mixtureCorFun, nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(0.4*Exp(0.1) + 0.6*Exp(3))", fileNameRoot="mix")
  print("Saving dataset...")
  save(simulationData, file="mixtureDataSet.RData")
}

getSimulationDataSetsGivenCovariance = function(corFun, nTotal=1000, nTest=100, marginalVar=1, errorVar=sqrt(.1), nDataSets=100, 
                                                printEvery=10, saveDataSetPlot=TRUE, fileNameRoot="", plotNameRoot="", 
                                                doPredGrid=TRUE, nTestGrid=70^2, useKenyaLocations=FALSE, urbanOverSamplefrac=0) {
  
  # set the spatial domain
  xRange = c(-1, 1)
  yRange = c(-1, 1)
  
  # if generating prediction grid, make that grid here
  if(doPredGrid && !useKenyaLocations) {
    nx = round(sqrt(nTestGrid))
    ny = nx
    gridTotal = nx * ny
    if(nx * ny != nTestGrid)
      stop("nTestGrid must be a square number if we are constructing a prediction grid")
    xValuesGrid = seq(-1, 1, l=nx)
    yValuesGrid = xValuesGrid
    glist = make.surface.grid(list(x=xValuesGrid, y=yValuesGrid))
    xValuesGrid = glist[,1]
    yValuesGrid = glist[,2]
  } else if(doPredGrid && useKenyaLocations) {
    if(urbanOverSamplefrac != 0)
      oversampleText = as.character(round(urbanOverSamplefrac, 4))
    else
      oversampleText = ""
    out = load(paste0("dataPointsKenya", oversampleText, ".RData"))
    nx = length(dataPointsKenya$xGrid)
    ny = length(dataPointsKenya$xGrid)
    gridTotal = nx
    xValuesGrid = dataPointsKenya$xGrid
    yValuesGrid = dataPointsKenya$yGrid
  } else if(!doPredGrid) {
    nx = 0
    ny = 0
    gridTotal = 0
    xValuesGrid = c()
    yValuesGrid = c()
  }
  
  # simulate observation spatial locations
  if(useKenyaLocations) {
    xValues = rbind(dataPointsKenya$xTrain, dataPointsKenya$xTest, dataPointsKenya$xTestRural, dataPointsKenya$xTestUrban)[,1:nDataSets]
    yValues = rbind(dataPointsKenya$yTrain, dataPointsKenya$yTest, dataPointsKenya$yTestRural, dataPointsKenya$yTestUrban)[,1:nDataSets]
    nTotal = nrow(xValues)
    nTest = nrow(dataPointsKenya$xTest) + nrow(dataPointsKenya$xTestRural) + nrow(dataPointsKenya$xTestUrban)
    overallTestI = 1:nrow(dataPointsKenya$xTest)
    ruralTestI = (nrow(dataPointsKenya$xTest) + 1):(nrow(dataPointsKenya$xTest) + nrow(dataPointsKenya$xTestRural))
    urbanTestI = (nrow(dataPointsKenya$xTest) + nrow(dataPointsKenya$xTestRural) + 1):(nrow(dataPointsKenya$xTest) + nrow(dataPointsKenya$xTestRural) + nrow(dataPointsKenya$xTestUrban))
    xRange = dataPointsKenya$xRange
    yRange = dataPointsKenya$yRange
  } else {
    xValues = matrix(runif(nTotal * nDataSets, xRange[1], xRange[2]), ncol=nDataSets)
    yValues = matrix(runif(nTotal * nDataSets, yRange[1], yRange[2]), ncol=nDataSets)
    overallTestI = NULL
    ruralTestI = NULL
    urbanTestI = NULL
  }
  
  # preallocate observation matrix, and pregenerate standard normal draws
  observations = matrix(nrow=nTotal+gridTotal, ncol=nDataSets)
  zsims = matrix(rnorm((nTotal + gridTotal) * nDataSets), ncol=nDataSets)
  
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
  observations = observations * sqrt(marginalVar) + matrix(rnorm((nTotal + gridTotal) * nDataSets, sd=sqrt(errorVar)), ncol=nDataSets)
  
  # separate out test and train results
  trainI = 1:(nTotal - nTest)
  testI = (nTotal - nTest + 1):nTotal
  gridI = (nTotal + 1):(nTotal + gridTotal)
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
  } else if(!doPredGrid) {
    xGrid = NULL
    yGrid = NULL
    zGrid = NULL
  }
  
  # separate out test locations into overall, urban, and rural if necessary
  if(useKenyaLocations) {
    nPerCounty = nrow(dataPointsKenya$xTest) / 47
    overallI = 1:(47*nPerCounty)
    ruralI = (47*nPerCounty + 1):((47 + 45)*nPerCounty)
    urbanI = ((47 + 45)*nPerCounty + 1):((47 + 45 + 47)*nPerCounty)
    
    xTestRural = xTest[ruralI,]
    xTestUrban = xTest[urbanI,]
    xTest = xTest[overallI,]
    yTestRural = yTest[ruralI,]
    yTestUrban = yTest[urbanI,]
    yTest = yTest[overallI,]
    zTestRural = zTest[ruralI,]
    zTestUrban = zTest[urbanI,]
    zTest = zTest[overallI,]
    
    # put relevant values into a list (make sure to include urban and rural testing values separately)
    out = list(xTrain=xTrain, yTrain=yTrain, zTrain=zTrain, xTest=xTest, yTest=yTest, zTest=zTest, 
               xTestRural=xTestRural, yTestRural=yTestRural, zTestRural=zTestRural, 
               xTestUrban=xTestUrban, yTestUrban=yTestUrban, zTestUrban=zTestUrban, 
               xGrid=xGrid, yGrid=yGrid, zGrid=zGrid, xValuesGrid=xValuesGrid, yValuesGrid=yValuesGrid, nx=nx, ny=ny, 
               corFun=corFun, marginalVar=marginalVar, errorVar=errorVar, xRange=xRange, yRange=yRange, 
               overallTestI=overallTestI, ruralTestI=ruralTestI, urbanTestI=urbanTestI)
  } else {
    # put relevant values into a list
    out = list(xTrain=xTrain, yTrain=yTrain, zTrain=zTrain, xTest=xTest, yTest=yTest, zTest=zTest, 
               xGrid=xGrid, yGrid=yGrid, zGrid=zGrid, xValuesGrid=xValuesGrid, yValuesGrid=yValuesGrid, nx=nx, ny=ny, 
               corFun=corFun, marginalVar=marginalVar, errorVar=errorVar, xRange=xRange, yRange=yRange)
  }
  
  # plot the results
  plotExampleDataSets(out, saveDataSetPlot=saveDataSetPlot, plotNameRoot=plotNameRoot, fileNameRoot=fileNameRoot)
  
  # return the results
  out
}

plotExampleDataSets = function(simulatedDataSets, saveDataSetPlot=TRUE, plotNameRoot="", fileNameRoot="", nx=100, ny=100) {
  # load and variable names
  xTrain = simulatedDataSets$xTrain
  yTrain = simulatedDataSets$yTrain
  zTrain = simulatedDataSets$zTrain
  xRange = simulatedDataSets$xRange
  yRange = simulatedDataSets$yRange
  
  if(saveDataSetPlot)
    pdf(file=paste0("Figures/sampleDataSet_", fileNameRoot, ".pdf"))
  
  quilt.plot(xTrain[,1], yTrain[,1], zTrain[,1], xlim=xRange, ylim=yRange, 
             main=paste("Example simulated dataset", plotNameRoot), nx=nx, ny=ny)
  
  if(saveDataSetPlot)
    dev.off()
}




