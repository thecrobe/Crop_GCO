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
fishnet<- readOGR(dsn= "Crop_GCO/Global/GIS/", layer="Fishnet_yield_NoAntarctica")
summary(fishnet)
#Set Projection
proj4string(fishnet) <- CRS("+init=epsg:3786") 

#Some covariates
rapeseed.model<-na.omit(read.csv(file="Crop_GCO/Global/Models/Rapeseed_RF.csv", header=T))
mapping<-read.csv(file="Crop_GCO/Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, rapeseed.model, by='Fishnet_ID')
rapeseed<-sp::merge(m, mapping, by='Fishnet_ID')
rapeseed.p = subset(rapeseed, rapeseed_HgHa > 0) #selecting values > 0 

#Transform numerical covariates
rapeseed.num<-rapeseed.p@data %>% select(rapeseed_HgHa, AET_mean, rapeseed_Fertilizer, Pesticide,GDP_Mean) 
rapeseed.scaled<-data.frame(log10(rapeseed.num+1)) #psuedo count to avoid inf
rapeseed.scaled$FISHNET_ID<-as.numeric(rapeseed.p@data$Fishnet_ID)
#Select categorical covariates
rapeseed.cat<-rapeseed.p %>% select(COUNTRY.x,RapeseedGCO, FISHNET_ID,Latitude,Longitude)
rapeseed.final<-merge(rapeseed.catagory,rapeseed.scaled, by="FISHNET_ID")
rapeseed.final$rapeseedGCO<-as.factor(rapeseed.final$RapeseedGCO)
rapeseed.final$rapeseedBinaryGCO<- ifelse(rapeseed.final$rapeseedGCO == 'Inside', 1, 0) #set binary dummy
rapeseed.final<-na.omit(rapeseed.final)
rapeseed.coords<-rapeseed.final %>% select(Latitude,Longitude) #select coordinates
summary(rapeseed.final)
dim(rapeseed.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rapeseed.coords$Latitude, rapeseed.coords$Longitude)))
dim(rapeseed.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rapeseed.coords$Latitude
y<-rapeseed.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rapeseed.data<-rapeseed.final

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = rapeseed.data, 
  dependent.variable.name = "rapeseed_HgHa",
  predictor.variable.names = c("AET_mean","rapeseed_Fertilizer", "Pesticide", "GDP_Mean","rapeseedBinaryGCO"),
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
rapeseed.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster
)
print(rapeseed.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

rapeseed.resid<-data.frame(rapeseed.spatial$residuals$values) #extract residuals
dim(rapeseed.resid)
#Plot residuals
ggplot(rapeseed.data, aes(x=Longitude, y=Latitude, color=rapeseed.resid$rapeseed.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
rapeseed.curves.df <- spatialRF::get_response_curves(rapeseed.spatial,variables = c("AET_mean","rapeseed_Fertilizer", "Pesticide", "GDP_Mean","rapeseedBinaryGCO"))
write.csv(rapeseed.spatial$predictions, file="Crop_GCO/Global/ModelOutputs/rapeseed_pred.csv")
write.csv(rapeseed.spatial$residuals$values, file="Crop_GCO/Global/ModelOutputs/rapeseed_resid.csv")
write.csv(rapeseed.spatial$variable.importance, file="Crop_GCO/Global/ModelOutputs/rapeseed_VarImp.csv")
write.csv(rapeseed.spatial$performance, file="Crop_GCO/Global/ModelOutputs/rapeseed_preformance.csv")
write.csv(rapeseed.curves.df, file="Crop_GCO/Global/ModelOutputs/rapeseed_curves.csv")

