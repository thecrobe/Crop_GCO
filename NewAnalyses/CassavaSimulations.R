library(ggplot2)
library(ranger)
library(Rcpp)
library(parallel)
library(dplyr)
library(spatialRF) 

cassava<-read.csv(file = "./NewAnalyses/cassava.simulationfinal.csv", header=T)

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(cassava$Latitude, cassava$Longitude)))

#coordinates of the cases
x<-cassava$Latitude
y<-cassava$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(1,5,10,50)
#random seed for reproducibility
random.seed <- 1

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 14,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

# 5x yield outside GCO simulation
Sim10x_Model <- spatialRF::rf_spatial(
  dependent.variable.name = "Sim5x",
  predictor.variable.names = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"),
  data = cassava,
  distance.matrix = distance.matrix,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
)

print(Sim5xSpatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc() #cleanup 
