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
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
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
cassava.scaled<-data.frame(log10(cassava.num+1)) #psuedo count to avoid inf
cassava.scaled$FISHNET_ID<-as.numeric(cassava.p@data$Fishnet_ID)
#Select categorical covariates
cassava.cat<-cassava.p %>% select(COUNTRY.x, FISHNET_ID,Latitude,Longitude)

cassava.final<-cbind(cassava.cat,cassava.scaled)
cassava.final<-na.omit(cassava.final)
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
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
cassava.data<-data.frame(cassava.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 3,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = cassava.data, 
  dependent.variable.name = "cassava_HgHa",
  predictor.variable.names = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  cluster = local.cluster,
  verbose = FALSE)

print(model.non.spatial)

gc() #cleanup ram
# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

######

# Spatial Model
cassava.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster,
  max.spatial.predictors = 10
)
print(cassava.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

cassava.resid<-data.frame(cassava.spatial$residuals$values) #extract residuals
dim(cassava.resid)
#Plot residuals
ggplot(cassava.data, aes(x=Longitude, y=Latitude, color=cassava.resid$cassava.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
cassava.curves.df <- spatialRF::get_response_curves(cassava.spatial,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean"))
write.csv(cassava.spatial$predictions, file="NewAnalyses/SpatialRandomForests_NoGCO/cassava_nogco_predictions")
write.csv(cassava.spatial$residuals$values, file="NewAnalyses/SpatialRandomForests_NoGCO/cassava_nogco_resid")
write.csv(cassava.spatial$variable.importance, file="NewAnalyses/SpatialRandomForests_NoGCO/cassava_nogco_varImp")
write.csv(cassava.spatial$performance, file="NewAnalyses/SpatialRandomForests_NoGCO/cassava_nogco_performance")
write.csv(cassava.curves.df, file="NewAnalyses/SpatialRandomForests_NoGCO/cassava_nogco_curves")
