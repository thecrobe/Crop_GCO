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
soybean.model<-na.omit(read.csv(file="Global/Models/Soybean_RF.csv", header=T))
mapping<-read.csv(file="Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, soybean.model, by='Fishnet_ID')
soybean<-sp::merge(m, mapping, by='Fishnet_ID')
soybean.p = subset(soybean, soybean_HgHa > 0) #selecting values > 0 
#Transform numerical covariates

soybean.num<-soybean.p@data %>% select(soybean_HgHa, AET_mean, soybean_Fertilizer, Pesticide,GDP_Mean) 
soybean.scaled<-data.frame(log10(soybean.num+1)) #psuedo count to avoid inf
soybean.scaled$FISHNET_ID<-as.numeric(soybean.p@data$Fishnet_ID)
#Select categorical covariates
soybean.cat<-soybean.p %>% select(COUNTRY.x, FISHNET_ID,Latitude,Longitude)

soybean.final<-cbind(soybean.cat,soybean.scaled)
soybean.final<-na.omit(soybean.final@data)
soybean.coords<-soybean.final %>% select(Latitude,Longitude) #select coordinates
summary(soybean.final)
dim(soybean.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(soybean.coords$Latitude, soybean.coords$Longitude)))
dim(soybean.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-soybean.coords$Latitude
y<-soybean.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1,5,10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
soybean.data<-data.frame(soybean.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = soybean.data, 
  dependent.variable.name = "soybean_HgHa",
  predictor.variable.names = c("AET_mean","soybean_Fertilizer", "Pesticide", "GDP_Mean"),
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
soybean.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster,
  max.spatial.predictors = 5,
)
print(soybean.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

soybean.resid<-data.frame(soybean.spatial$residuals$values) #extract residuals
dim(soybean.resid)
#Plot residuals
ggplot(soybean.data, aes(x=Longitude, y=Latitude, color=soybean.resid$soybean.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
soybean.curves.df <- spatialRF::get_response_curves(soybean.spatial,variables = c("AET_mean","soybean_Fertilizer", "Pesticide", "GDP_Mean"))
write.csv(soybean.spatial$predictions, file="Global/SpatialRandomForests_NoGCO/soybean_nogco_predictions")
write.csv(soybean.spatial$residuals$values, file="Global/SpatialRandomForests_NoGCO/soybean_nogco_resid")
write.csv(soybean.spatial$variable.importance, file="Global/SpatialRandomForests_NoGCO/soybean_nogco_varImp")
write.csv(soybean.spatial$performance, file="Global/SpatialRandomForests_NoGCO/soybean_nogco_performance")
write.csv(soybean.curves.df, file="Global/SpatialRandomForests_NoGCO/soybean_nogco_curves")
