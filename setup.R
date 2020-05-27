# load required packages
library(Matrix)
library(spam)
library(fields)
library(LatticeKrig)
library(MCMCpack)
library(INLA)
library(invgamma)

# load other recommended packages
library(latex2exp)
library(xtable)
library(geosphere)
library(viridis)
library(colorspace)
library(sp)
library(raster)

# load required R scripts
setwd("~/git/ELK/")
source('ELK.R')
source('ELK_rgeneric.R')
source('modELK.R')
source('modLK.R')
source("scores.R")
source('utilityFuns.R')
source('getSimulationDataSets.R')
source('fitRainfallData.R')