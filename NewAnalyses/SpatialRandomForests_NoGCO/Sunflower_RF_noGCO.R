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
sunflower.model<-na.omit(read.csv(file="Global/Models/Sunflower_RF.csv", header=T))
mapping<-read.csv(file="Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, sunflower.model, by='Fishnet_ID')
sunflower<-sp::merge(m, mapping, by='Fishnet_ID')
sunflower.p = subset(sunflower, sunflower_HgHa > 0) #selecting values > 0 
#Transform numerical covariates

sunflower.num<-sunflower.p@data %>% select(sunflower_HgHa, AET_mean, sunflower_Fertilizer, Pesticide,GDP_Mean) 
sunflower.scaled<-data.frame(log10(sunflower.num+1)) #psuedo count to avoid inf
sunflower.scaled$FISHNET_ID<-as.numeric(sunflower.p@data$Fishnet_ID)
#Select categorical covariates
sunflower.cat<-sunflower.p %>% select(COUNTRY.x, FISHNET_ID,Latitude,Longitude)

sunflower.final<-cbind(sunflower.cat,sunflower.scaled)
sunflower.final<-na.omit(sunflower.final@data)
sunflower.coords<-sunflower.final %>% select(Latitude,Longitude) #select coordinates
summary(sunflower.final)
dim(sunflower.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(sunflower.coords$Latitude, sunflower.coords$Longitude)))
dim(sunflower.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-sunflower.coords$Latitude
y<-sunflower.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1,5,10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sunflower.data<-data.frame(sunflower.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = sunflower.data, 
  dependent.variable.name = "sunflower_HgHa",
  predictor.variable.names = c("AET_mean","sunflower_Fertilizer", "Pesticide", "GDP_Mean"),
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
sunflower.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster,
  max.spatial.predictors = 5,
)
print(sunflower.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

sunflower.resid<-data.frame(sunflower.spatial$residuals$values) #extract residuals
dim(sunflower.resid)
#Plot residuals
ggplot(sunflower.data, aes(x=Longitude, y=Latitude, color=sunflower.resid$sunflower.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
sunflower.curves.df <- spatialRF::get_response_curves(sunflower.spatial,variables = c("AET_mean","sunflower_Fertilizer", "Pesticide", "GDP_Mean"))
write.csv(sunflower.spatial$predictions, file="Global/SpatialRandomForests_NoGCO/sunflower_nogco_predictions")
write.csv(sunflower.spatial$residuals$values, file="Global/SpatialRandomForests_NoGCO/sunflower_nogco_resid")
write.csv(sunflower.spatial$variable.importance, file="Global/SpatialRandomForests_NoGCO/sunflower_nogco_varImp")
write.csv(sunflower.spatial$performance, file="Global/SpatialRandomForests_NoGCO/sunflower_nogco_performance")
write.csv(sunflower.curves.df, file="Global/SpatialRandomForests_NoGCO/sunflower_nogco_curves")
