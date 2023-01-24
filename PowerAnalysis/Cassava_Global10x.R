library(spatialRF) 
library(spdep)
library(sf)
library(spdplyr)
library(rgdal)
library(ggplot2)
library(ranger)
library(Rcpp)
library(parallel)
library(dplyr)

#Read In
#Fishnet- pixel = 100km^2
fishnet<- readOGR(dsn= "Global/GIS/", layer="Fishnet_yield_NoAntarctica")
summary(fishnet)
#Set Projection
proj4string(fishnet) <- CRS("+init=epsg:3786") 

#Some covariates
cassava.model<-na.omit(read.csv(file="Global/Models/Cassava_RF.csv", header=T))
mapping<-read.csv(file="Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, cassava.model, by='Fishnet_ID')
cassava<-sp::merge(m, mapping, by='Fishnet_ID')
cassava.p = subset(cassava, cassava_HgHa > 0) #selecting values > 0 

#Transform numerical covariates
cassava.num<-cassava.p@data %>% select(cassava_HgHa, AET_mean, cassava_Fertilizer, Pesticide,GDP_Mean) 
cassava.num$FISHNET_ID<-as.numeric(cassava.p@data$Fishnet_ID)
#Select categorical covariates
cassava.cat<-cassava.p %>% select(CassavaGCO, FISHNET_ID,Latitude,Longitude)

#Merge
cassava.final<-merge(cassava.catagory,cassava.num, by="FISHNET_ID")
cassava.final$cassavaGCO<-as.factor(cassava.final$CassavaGCO)
cassava.final$cassavaGCO
cassava.final$cassavaBinaryGCO<-(ifelse(cassava.final$CassavaGCO == 'Inside', 1, 0)) #set binary dummy

cassava.coords<-cassava.final %>% select(Latitude,Longitude) #select coordinates
summary(cassava.final)
dim(cassava.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(cassava.coords$Latitude, cassava.coords$Longitude)))
dim(cassava.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-cassava.coords$Latitude
y<-cassava.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
cassava.data<-cassava.final

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 14,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

# Spatial Model
Trial <- spatialRF::rf_spatial(
  dependent.variable.name = "Sim10x",
  predictor.variable.names = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"),
  data = cassava.data,
  distance.matrix = distance.matrix,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)

print(Sim10xSpatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()
