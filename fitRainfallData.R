# script for fitting North American rainfall dataset

# fit North American rainfall data set
# LatticeKrig paper also uses NorthAmericanRainfall dataset
# 1720 stations precip in JJA=June,July,August based on data from 1950-2010 (61 years).  
# Also includes elevation.
# dataset in LK package only includes linear trends and intercepts for each station
# 1720 x 1 = 1720 observations (or, with full dataset, 1720 x 61 = 104,920 observations)
# int.strategy: integration strategy of INLA.
# doCovs: whether or not elevation is included as a covariate
testRainfallData = function(normalize=TRUE, NC=5, nLayer=3, nBuffer=5, seed=1, maxit=25, 
                            int.strategy="ccd", strategy="gaussian", fastNormalize=TRUE) {
  # load it
  data(NorthAmericanRainfall)
  x<- cbind(NorthAmericanRainfall$longitude,  NorthAmericanRainfall$latitude)
  y<- log(NorthAmericanRainfall$precip/10) # to mm units (take log first?)
  x = mercator(x)/1000 # m to km units
  X = matrix(NorthAmericanRainfall$elevation, ncol=1) # include elevation as covariate if desired
  
  # plot it
  png("Figures/Rainfall/rainfallObservations.png", width=500, height=500)
  par(mfrow=c(1,1))
  quilt.plot( x,y, main="Mean JJA Precipitation, 1950-2010 (log mm)", xlab="Longitude", ylab="Latitude")
  world( add=TRUE)
  dev.off()
  
  # set up prediction locations for comparison with LatticeKrig
  latLim = c(31.75, 47.25)
  lonLim = c(-115, -101)
  lims = mercator(cbind(lonLim, latLim))/1000
  eastLim = lims[,1]
  northLim = lims[,2]
  # predLon = seq(lonLim[1], lonLim[2], length = 50)
  # predLat = seq(latLim[1], latLim[2], length = 50)
  predEast = seq(eastLim[1], eastLim[2], length = 50)
  predNorth = seq(northLim[1], northLim[2], length = 50)
  predCoords = make.surface.grid(list(east = predEast, north = predNorth))
  predCoordsLonLat = mercator(predCoords*1000, inverse=TRUE)
  Xpred = matrix(getElevation(predCoordsLonLat[,1], predCoordsLonLat[,2]), ncol=1)
  
  # fit the models
  LKITime = system.time(LKIout <- fitLKINLAStandard2(x, y, predCoords, nu=NULL, seed, nLayer, NC, nBuffer, priorPar=getPCPrior(2500, .1, 700/2, nLayer=nLayer), 
                                                     cbind(1, x, X), cbind(1, predCoords, Xpred), normalize=normalize, intStrategy=int.strategy, 
                                                     fastNormalize=fastNormalize))
  
  LKTime = system.time(LKout <- fitLKStandard(x, y, predCoords, X, Xpred, NC, nLayer, normalize, nBuffer, nu=NULL, verbose=TRUE, 
                                              normalize=normalize, maxit=maxit))
  
  # print out computation times
  print(paste0("LatticeKrig computation time: ", LKTime))
  print(paste0("LatticeKrig-INLA computation time: ", LKITime))
  
  # transform log scale predictions and standard deviations to standard scale
  LKMeans = exp(LKout$preds + LKout$SEs^2/2)
  LKVars = (exp(LKout$SEs^2)-1) * LKMeans^2
  LKIMeans = exp(LKIout$preds + LKIout$sigmas^2/2)
  LKIVars = (exp(LKIout$sigmas^2)-1) * LKIMeans^2
  
  browser()
  
  # show predictions just for LatticeKrig
  zlimPreds = range(LKMeans)
  zlimSEs = range(sqrt(LKVars))
  pdf(file="Figures/Rainfall/rainfallPredsLK.pdf", width=7, height=7)
  par(mfrow=c(1,2))
  quilt.plot(predCoords, LKMeans, main="LatticeKrig predictions (in mm)", zlim=zlimPreds, 
             col=viridis(64))
  points(x, pch=19, cex=.4)
  quilt.plot(predCoords, sqrt(LKVars), main="LatticeKrig SEs", zlim=zlimSEs)
  points(x, pch=19, cex=.4)
  US(add=TRUE)
  dev.off()
  
  # show it all
  zlimPreds = range(LKMeans, LKIMeans)
  zlimSEs = range(sqrt(LKVars), sqrt(LKIVars))
  pdf(file="Figures/Rainfall/rainfallPreds.pdf", width=7, height=7)
  par(mfrow=c(2,2))
  quilt.plot(predCoords, LKMeans, main="LatticeKrig predictions", zlim=zlimPreds)
  quilt.plot(predCoords, sqrt(LKVars), main="LatticeKrig SEs", zlim=zlimSEs)
  quilt.plot(predCoords, LKIMeans, main="LatticeKrig-INLA predictions", zlim=zlimPreds)
  quilt.plot(predCoords, sqrt(LKIVars), main="LatticeKrig-INLA SDs", zlim=zlimSEs)
  US(add=TRUE)
  dev.off()
  
  pdf(file="Figures/Rainfall/rainfallPredsVs.pdf", width=10, height=5)
  par(mfrow=c(1,2))
  plot(LKMeans, LKIMeans, pch=".", main="LK-INLA vs LatticeKrig predictions", 
       xlab="LatticeKrig", ylab="LK-INLA")
  abline(0,1, col="green")
  plot(sqrt(LKVars), sqrt(LKIVars), pch=".", log="xy", xlab="LatticeKrig SEs", ylab="LK-INLA SDs", 
       main="LK-INLA vs LatticeKrig predictive uncertainty")
  abline(0,1, col="green")
  dev.off()
  
  kappaMarg = inla.tmarginal(function(x) {sqrt(8)/exp(x) * latticeWidth}, LKIout$mod$marginals.hyperpar$`Theta1 for field`)
  pdf(file=paste0("Figures/Rainfall/posteriorKappa.pdf"), width=5, height=5)
  plot(kappaMarg, type="l", xlab=TeX("$\\kappa$"), main=TeX("Marginal for $\\kappa$"))
  abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  abline(v=mean(sqrt(LKout$a.wghtVals - 4)), col="green")
  legend("topright", c("95% CI", "LK Est"), lty=c(2, 1), col=c("purple", "green"))
  dev.off()
  
  rhoMarg = inla.tmarginal(function(x) {exp(x)}, LKIout$mod$marginals.hyperpar$`Theta2 for field`)
  pdf(file="Figures/posteriorRho.pdf", width=5, height=5)
  plot(rhoMarg, type="l", xlab=TeX("$\\rho$"), main=TeX("Marginal for $\\rho$"))
  abline(v=inla.qmarginal(probs=c(.025, .975), rhoMarg), col="purple")
  abline(v=mean(LKout$rhoVals), col="green")
  legend("topright", c("95% CI", "LK Est"), lty=c(2, 1), col=c("purple", "green"))
  dev.off()
  
  ## Now generate marginals for the alpha parameters if they exist. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  out = inla.hyperpar.sample(20000, LKIout$mod, improve.marginals=TRUE)
  if(nLayer >= 2) {
    if(!separateRanges) {
      zSamples = out[,3:(2+nLayer-1)]
      xSamples = apply(zSamples, 1, multivariateExpit)
      xSamples = rbind(xSamples, 1-colSums(xSamples))
    } else {
      zSamples = matrix(out[,(nLayer+1+1):(nLayer + 1 + 1 + nLayer-2)], ncol=nLayer-1)
      xSamples = matrix(apply(zSamples, 1, multivariateExpit), nrow=nLayer-1)
      xSamples = rbind(xSamples, 1-colSums(xSamples))
    }
    
    for(l in 1:nLayer) {
      theseQuantiles = quantile(probs=c(.025, .975), xSamples[l,])
      thisQuantileRange = abs(diff(theseQuantiles))
      # if(thisQuantileRange <= 0.2) {
      #   plotRange = range(xSamples[l,])
      # } else {
      plotRange = c(0, 1)
      # }
      if(theseQuantiles[1]-0 <= 1-theseQuantiles[2])
        legendLocation = "topright"
      else
        legendLocation = "topleft"
      pdf(file=paste0("Figures/Rainfall/posteriorAlpha", l, ".pdf"), width=5, height=5)
      hist(xSamples[l,], xlab=TeX(paste0("$\\alpha_", l, "$")), main=TeX(paste0("Marginal for $\\alpha_", l, "$, LK estimate (green)")), breaks=100, freq=F, xlim=plotRange)
      abline(v=mean(xSamples[l,]), col="purple", lty=1)
      abline(v=theseQuantiles, col="purple", lty=2)
      abline(v=rowMeans(LKout$alphaVals)[l], col="green")
      legend(legendLocation, c("95% CI", "LK Est"), lty=c(2, 1), col=c("purple", "green"))
      dev.off()
    }
  }
  
  browser()
  
  # produce figures like in the LatticeKrig paper of central predictions and standard errors?
  
  list(LKout=LKout, LKIout=LKIout, LKTime=LKTime, LKITime=LKITime)
  # how good are the INLA results?  Compare with MCMC?
}

# get elevation at the given set of lat/lon coordinates
# uses GMTED elevation data available from https://www.usgs.gov/land-resources/eros/coastal-changes-and-impacts/gmted2010?qt-science_support_page_related_con=0#qt-science_support_page_related_con
getElevation = function(lon, lat) {
  # require(rgdal)
  # require(raster)
  # require(sp)
  # # ogrListLayers("gtopo30/gtopo30.shp")
  # shape = readOGR("gtopo30/")
  # shape <- shapefile("gtopo30/gtopo30.shp")
  # test = extract(shape, SpatialPoints(cbind(lon, lat)),method="bilinear")
  # datapol <- data.frame(shape)
  # pointtoplot <- data.frame(x=-20, y=40)
  # coordinates(pointtoplot) <- ~ x + y 
  # proj4string(pointtoplot) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  # test <- data.frame(xx=over(shape, pointtoplot))
  # combine <- cbind(test, datapol)
  # combine <- na.omit(combine)
  
  temp = raster("elevation.nc")
  extract(temp, SpatialPoints(cbind(lon, lat)),method="bilinear")
}

# uses GMTED elevation data available from https://www.usgs.gov/land-resources/eros/coastal-changes-and-impacts/gmted2010?qt-science_support_page_related_con=0#qt-science_support_page_related_con
# creates a single-layer netcdf file readable using the raster command
convertNC2raster = function() {
  # require(gdalUtils)
  # gdalinfo("GMTED2010_15n015_00625deg.hdf")
  # hdf_dataset <- system.file("GMTED2010_15n015_00625deg.hdf")
  # gdal_translate("GMTED2010_15n015_00625deg.hdf", sds=T, of="GTiff", dst_dataset = "test.tif")
  # sds
  # 
  # test <- raster("GMTED2010_15n015_00625deg.hdf")
  
  library(ncdf4)
  library(raster)
  library(rasterVis)
  library(maptools)
  library(maps)
  # tmpin <- raster("GMTED2010_15n015_00625deg.nc", values= TRUE)
  
  # open netcdf file
  ncObject = nc_open("GMTED2010_15n015_00625deg.nc")
  
  # get desired variables and their attributes
  latVar = ncvar_get(ncObject, "latitude")
  lonVar = ncvar_get(ncObject, "longitude")
  lat = ncatt_get(ncObject, "latitude")
  lon = ncatt_get(ncObject, "longitude")
  elevation = ncatt_get(ncObject, "elevation")
  elevationVar = ncvar_get(ncObject, "elevation")
  
  # reformat the variables
  londim <- ncdim_def("lon", "degrees_east", as.double(lonVar))
  latdim <- ncdim_def("lat", "degrees_north", as.double(latVar))
  tmp.def <- ncvar_def("elevation", "m", list(londim, latdim), -1e20, 
                       prec = "double")
  
  # create new netcdf file
  ncfname <- "elevation.nc"
  ncout <- nc_create(ncfname, list(tmp.def), force_v4 = TRUE)
  
  # put the array
  tmat = elevationVar
  ncvar_put(ncout, tmp.def, tmat)
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout, "lon", "axis", "X")  
  ncatt_put(ncout, "lat", "axis", "Y")
  
  # add global attributes
  title <- "Elevation netCDF file"
  ncatt_put(ncout, 0, "title", title)
  
  # close the file, writing data to disk
  nc_close(ncout)
  
  invisible(NULL)
}









