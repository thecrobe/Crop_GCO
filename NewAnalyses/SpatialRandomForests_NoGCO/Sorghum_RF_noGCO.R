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
fishnet<- readOGR(dsn= "GIS", layer="Fishnet_yield_NoAntarctica")
#Set Projection
proj4string(fishnet) <- CRS("+init=epsg:3786") 

#Some covariates
sorghum.model<-na.omit(read.csv(file="Global/Models/Sorghum_RF.csv", header=T))
mapping<-read.csv(file="Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, sorghum.model, by='Fishnet_ID')
sorghum<-sp::merge(m, mapping, by='Fishnet_ID')
sorghum.p = subset(sorghum, mean_sorgh > 0) #selecting values > 0 
#Transform numerical covariates

sorghum.num<-sorghum.p@data %>% select(mean_sorgh, Evapotranspiraton, sorghum_Fertilizer, Pesticide,GDP) 
sorghum.scaled<-data.frame(log10(sorghum.num+1)) #psuedo count to avoid inf
sorghum.scaled$FISHNET_ID<-as.numeric(sorghum.p@data$Fishnet_ID)
#Select categorical covariates
sorghum.cat<-sorghum.p %>% select(COUNTRY.x, FISHNET_ID,Latitude,Longitude)

sorghum.final<-cbind(sorghum.cat,sorghum.scaled)
sorghum.final<-na.omit(sorghum.final@data)
sorghum.coords<-sorghum.final %>% select(Latitude,Longitude) #select coordinates
summary(sorghum.final)
dim(sorghum.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(sorghum.coords$Latitude, sorghum.coords$Longitude)))
dim(sorghum.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-sorghum.coords$Latitude
y<-sorghum.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1,5,10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sorghum.data<-data.frame(sorghum.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = sorghum.data, 
  dependent.variable.name = "mean_sorgh",
  predictor.variable.names = c("Evapotranspiraton","sorghum_Fertilizer", "Pesticide", "GDP"),
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
sorghum.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster,
  max.spatial.predictors = 5,
)
print(sorghum.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

sorghum.resid<-data.frame(sorghum.spatial$residuals$values) #extract residuals
dim(sorghum.resid)
#Plot residuals
ggplot(sorghum.data, aes(x=Longitude, y=Latitude, color=sorghum.resid$sorghum.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
sorghum.curves.df <- spatialRF::get_response_curves(sorghum.spatial,variables = c("Evapotranspiraton","sorghum_Fertilizer", "Pesticide", "GDP"))
write.csv(sorghum.spatial$predictions, file="Global/SpatialRandomForests_NoGCO/sorghum_nogco_predictions")
write.csv(sorghum.spatial$residuals$values, file="Global/SpatialRandomForests_NoGCO/sorghum_nogco_resid")
write.csv(sorghum.spatial$variable.importance, file="Global/SpatialRandomForests_NoGCO/sorghum_nogco_varImp")
write.csv(sorghum.spatial$performance, file="Global/SpatialRandomForests_NoGCO/sorghum_nogco_performance")
write.csv(sorghum.curves.df, file="Global/SpatialRandomForests_NoGCO/sorghum_nogco_curves")
