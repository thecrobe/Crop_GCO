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
rice.model<-na.omit(read.csv(file="Global/Models/Rice_RF.csv", header=T))
mapping<-read.csv(file="Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, rice.model, by='Fishnet_ID')
rice<-sp::merge(m, mapping, by='Fishnet_ID')
rice.p = subset(rice, rice_HgHa > 0) #selecting values > 0 
#Transform numerical covariates

rice.num<-rice.p@data %>% select(rice_HgHa, AET_mean, rice_Fertilizer, Pesticide,GDP_Mean) 
rice.scaled<-data.frame(log10(rice.num+1)) #psuedo count to avoid inf
rice.scaled$FISHNET_ID<-as.numeric(rice.p@data$Fishnet_ID)
#Select categorical covariates
rice.cat<-rice.p %>% select(COUNTRY.x, FISHNET_ID,Latitude,Longitude)

rice.final<-cbind(rice.cat,rice.scaled)
rice.final<-na.omit(rice.final@data)
rice.coords<-rice.final %>% select(Latitude,Longitude) #select coordinates
summary(rice.final)
dim(rice.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rice.coords$Latitude, rice.coords$Longitude)))
dim(rice.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rice.coords$Latitude
y<-rice.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1,5,10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rice.data<-data.frame(rice.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = rice.data, 
  dependent.variable.name = "rice_HgHa",
  predictor.variable.names = c("AET_mean","rice_Fertilizer", "Pesticide", "GDP_Mean"),
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
rice.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster,
  max.spatial.predictors = 5,
)
print(rice.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

rice.resid<-data.frame(rice.spatial$residuals$values) #extract residuals
dim(rice.resid)
#Plot residuals
ggplot(rice.data, aes(x=Longitude, y=Latitude, color=rice.resid$rice.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
rice.curves.df <- spatialRF::get_response_curves(rice.spatial,variables = c("AET_mean","rice_Fertilizer", "Pesticide", "GDP_Mean"))
write.csv(rice.spatial$predictions, file="Global/SpatialRandomForests_NoGCO/rice_nogco_predictions")
write.csv(rice.spatial$residuals$values, file="Global/SpatialRandomForests_NoGCO/rice_nogco_resid")
write.csv(rice.spatial$variable.importance, file="Global/SpatialRandomForests_NoGCO/rice_nogco_varImp")
write.csv(rice.spatial$performance, file="Global/SpatialRandomForests_NoGCO/rice_nogco_performance")
write.csv(rice.curves.df, file="Global/SpatialRandomForests_NoGCO/rice_nogco_curves")
