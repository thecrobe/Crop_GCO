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
rye.model<-na.omit(read.csv(file="Global/Models/Rye_RF.csv", header=T))
mapping<-read.csv(file="Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, rye.model, by='Fishnet_ID')
rye<-sp::merge(m, mapping, by='Fishnet_ID')
rye.p = subset(rye, rye_HgHa > 0) #selecting values > 0 
#Transform numerical covariates

rye.num<-rye.p@data %>% select(rye_HgHa, AET_mean, rye_Fertilizer, Pesticide,GDP_Mean) 
rye.scaled<-data.frame(log10(rye.num+1)) #psuedo count to avoid inf
rye.scaled$FISHNET_ID<-as.numeric(rye.p@data$Fishnet_ID)
#Select categorical covariates
rye.cat<-rye.p %>% select(COUNTRY.x, FISHNET_ID,Latitude,Longitude)

rye.final<-cbind(rye.cat,rye.scaled)
rye.final<-na.omit(rye.final@data)
rye.coords<-rye.final %>% select(Latitude,Longitude) #select coordinates
summary(rye.final)
dim(rye.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rye.coords$Latitude, rye.coords$Longitude)))
dim(rye.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rye.coords$Latitude
y<-rye.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1,5,10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rye.data<-data.frame(rye.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = rye.data, 
  dependent.variable.name = "rye_HgHa",
  predictor.variable.names = c("AET_mean","rye_Fertilizer", "Pesticide", "GDP_Mean"),
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
rye.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster,
  max.spatial.predictors = 5,
)
print(rye.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

rye.resid<-data.frame(rye.spatial$residuals$values) #extract residuals
dim(rye.resid)
#Plot residuals
ggplot(rye.data, aes(x=Longitude, y=Latitude, color=rye.resid$rye.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
rye.curves.df <- spatialRF::get_response_curves(rye.spatial,variables = c("AET_mean","rye_Fertilizer", "Pesticide", "GDP_Mean"))
write.csv(rye.spatial$predictions, file="Global/SpatialRandomForests_NoGCO/rye_nogco_predictions")
write.csv(rye.spatial$residuals$values, file="Global/SpatialRandomForests_NoGCO/rye_nogco_resid")
write.csv(rye.spatial$variable.importance, file="Global/SpatialRandomForests_NoGCO/rye_nogco_varImp")
write.csv(rye.spatial$performance, file="Global/SpatialRandomForests_NoGCO/rye_nogco_performance")
write.csv(rye.curves.df, file="Global/SpatialRandomForests_NoGCO/rye_nogco_curves")

