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
groundnut.model<-na.omit(read.csv(file="Global/Models/Groundnut_RF.csv", header=T))
mapping<-read.csv(file="Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, groundnut.model, by='Fishnet_ID')
groundnut<-sp::merge(m, mapping, by='Fishnet_ID')
groundnut.p = subset(groundnut, groundnut_HgHa > 0) #selecting values > 0 
#Transform numerical covariates

groundnut.num<-groundnut.p@data %>% select(groundnut_HgHa, AET_mean, groundnut_Fertilizer, Pesticide,GDP_Mean) 
groundnut.scaled<-data.frame(log10(groundnut.num+1)) #psuedo count to avoid inf
groundnut.scaled$FISHNET_ID<-as.numeric(groundnut.p@data$Fishnet_ID)
#Select categorical covariates
groundnut.cat<-groundnut.p %>% select(COUNTRY.x, FISHNET_ID,Latitude,Longitude)

groundnut.final<-cbind(groundnut.cat,groundnut.scaled)
groundnut.final<-na.omit(groundnut.final)
groundnut.coords<-groundnut.final %>% select(Latitude,Longitude) #select coordinates
summary(groundnut.final)
dim(groundnut.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(groundnut.coords$Latitude, groundnut.coords$Longitude)))
dim(groundnut.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-groundnut.coords$Latitude
y<-groundnut.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
groundnut.data<-data.frame(groundnut.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 3,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = groundnut.data, 
  dependent.variable.name = "groundnut_HgHa",
  predictor.variable.names = c("AET_mean","groundnut_Fertilizer", "Pesticide", "GDP_Mean"),
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
groundnut.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster,
  max.spatial.predictors = 10
)
print(groundnut.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

groundnut.resid<-data.frame(groundnut.spatial$residuals$values) #extract residuals
dim(groundnut.resid)
#Plot residuals
ggplot(groundnut.data, aes(x=Longitude, y=Latitude, color=groundnut.resid$groundnut.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
groundnut.curves.df <- spatialRF::get_response_curves(groundnut.spatial,variables = c("AET_mean","groundnut_Fertilizer", "Pesticide", "GDP_Mean"))
write.csv(groundnut.spatial$predictions, file="NewAnalyses/SpatialRandomForests_NoGCO/groundnut_nogco_predictions")
write.csv(groundnut.spatial$residuals$values, file="NewAnalyses/SpatialRandomForests_NoGCO/groundnut_nogco_resid")
write.csv(groundnut.spatial$variable.importance, file="NewAnalyses/SpatialRandomForests_NoGCO/groundnut_nogco_varImp")
write.csv(groundnut.spatial$performance, file="NewAnalyses/SpatialRandomForests_NoGCO/groundnut_nogco_performance")
write.csv(groundnut.curves.df, file="NewAnalyses/SpatialRandomForests_NoGCO/groundnut_nogco_curves")

