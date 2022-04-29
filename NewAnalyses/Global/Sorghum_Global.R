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
sorghum.model<-na.omit(read.csv(file="Crop_GCO/Global/Models/Sorghum_RF.csv", header=T))
mapping<-read.csv(file="Crop_GCO/Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name
sorghum.model$Evapotranspiraton
#Merging fishnet with covariates
m <- sp::merge(fishnet, sorghum.model, by='Fishnet_ID')
sorghum<-sp::merge(m, mapping, by='Fishnet_ID')
sorghum.p = subset(sorghum, mean_sorgh > 0) #selecting values > 0 

#Transform numerical covariates
sorghum.num<-sorghum.p@data %>% select(mean_sorgh, Evapotranspiraton, sorghum_Fertilizer, Pesticide,GDP) 
sorghum.scaled<-data.frame(log10(sorghum.num+1)) #psuedo count to avoid inf
sorghum.scaled$FISHNET_ID<-as.numeric(sorghum.p@data$Fishnet_ID)
#Select categorical covariates
sorghum.cat<-sorghum.p %>% select(COUNTRY.x,SorghumGCO, FISHNET_ID,Latitude,Longitude)
sorghum.final<-merge(sorghum.catagory,sorghum.scaled, by="FISHNET_ID")
sorghum.final$sorghumGCO<-as.factor(sorghum.final$SorghumGCO)
sorghum.final$sorghumBinaryGCO<- ifelse(sorghum.final$sorghumGCO == 'Inside', 1, 0) #set binary dummy
sorghum.final<-na.omit(sorghum.final)
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
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sorghum.data<-sorghum.final

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 2,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = sorghum.data, 
  dependent.variable.name = "mean_sorgh",
  predictor.variable.names = c("Evapotranspiraton","sorghum_Fertilizer", "Pesticide", "GDP","sorghumBinaryGCO"),
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
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster
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
sorghum.curves.df <- spatialRF::get_response_curves(sorghum.spatial,variables = c("Evapotranspiraton","sorghum_Fertilizer", "Pesticide", "GDP","sorghumBinaryGCO"))
write.csv(sorghum.spatial$predictions, file="Crop_GCO/Global/ModelOutputs/sorghum_pred.csv")
write.csv(sorghum.spatial$residuals$values, file="Crop_GCO/Global/ModelOutputs/sorghum_resid.csv")
write.csv(sorghum.spatial$variable.importance, file="Crop_GCO/Global/ModelOutputs/sorghum_VarImp.csv")
write.csv(sorghum.spatial$performance, file="Crop_GCO/Global/ModelOutputs/sorghum_preformance.csv")
write.csv(sorghum.curves.df, file="Crop_GCO/Global/ModelOutputs/sorghum_curves.csv")

