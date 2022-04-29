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
wheat.model<-na.omit(read.csv(file="Crop_GCO/Global/Models/Wheat_RF.csv", header=T))
mapping<-read.csv(file="Crop_GCO/Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, wheat.model, by='Fishnet_ID')
wheat<-sp::merge(m, mapping, by='Fishnet_ID')
wheat.p = subset(wheat, wheat_HgHa > 0) #selecting values > 0 
#Transform numerical covariates

wheat.num<-wheat.p@data %>% select(wheat_HgHa, AET_mean, wheat_Fertilizer, Pesticide,GDP_Mean) 
wheat.scaled<-data.frame(log10(wheat.num+1)) #psuedo count to avoid inf
wheat.scaled$FISHNET_ID<-as.numeric(wheat.p@data$Fishnet_ID)
#Select categorical covariates
wheat.cat<-wheat.p %>% select(COUNTRY.x,WheatGCO, FISHNET_ID,Latitude,Longitude)

wheat.final<-merge(wheat.cat,wheat.scaled, by="FISHNET_ID")
wheat.final$wheatGCO<-as.factor(wheat.final$WheatGCO)
wheat.final$wheatBinaryGCO<- ifelse(wheat.final$wheatGCO == 'Inside', 1, 0) #set binary dummy
wheat.final<-na.omit(wheat.final)
wheat.coords<-wheat.final %>% select(Latitude,Longitude) #select coordinates
summary(wheat.final)
dim(wheat.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(wheat.coords$Latitude, wheat.coords$Longitude)))
dim(wheat.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-wheat.coords$Latitude
y<-wheat.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
wheat.data<-data.frame(wheat.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 10,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = wheat.data, 
  dependent.variable.name = "wheat_HgHa",
  predictor.variable.names = c("AET_mean","wheat_Fertilizer", "Pesticide", "GDP_Mean","wheatBinaryGCO"),
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
wheat.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster
)
print(wheat.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

wheat.resid<-data.frame(wheat.spatial$residuals$values) #extract residuals
dim(wheat.resid)
#Plot residuals
ggplot(wheat.data, aes(x=Longitude, y=Latitude, color=wheat.resid$wheat.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
wheat.curves.df <- spatialRF::get_response_curves(wheat.spatial,variables = c("AET_mean","wheat_Fertilizer", "Pesticide", "GDP_Mean","wheatBinaryGCO"))
write.csv(wheat.spatial$predictions, file="Crop_GCO/Global/ModelOutputs/wheat_pred.csv")
write.csv(wheat.spatial$residuals$values, file="Crop_GCO/Global/ModelOutputs/wheat_resid.csv")
write.csv(wheat.spatial$variable.importance, file="Crop_GCO/Global/ModelOutputs/wheat_VarImp.csv")
write.csv(wheat.spatial$performance, file="Crop_GCO/Global/ModelOutputs/wheat_preformance.csv")
write.csv(wheat.curves.df, file="Crop_GCO/Global/ModelOutputs/wheat_curves.csv")

