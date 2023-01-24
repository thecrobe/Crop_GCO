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
potato.model<-na.omit(read.csv(file="Crop_GCO/Global/Models/Potato_RF.csv", header=T))
mapping<-read.csv(file="Crop_GCO/Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, potato.model, by='Fishnet_ID')
potato<-sp::merge(m, mapping, by='Fishnet_ID')
potato.p = subset(potato, mean_potat.x > 0) #selecting values > 0 
potato.p$po
#Transform numerical covariates
potato.p$fe
potato.num<-potato.p@data %>% select(mean_potat.x, AET_mean, Potato_Fe, Pesticide,GDP_Mean) 
potato.scaled<-data.frame(log10(potato.num+1)) #psuedo count to avoid inf
potato.scaled$FISHNET_ID<-as.numeric(potato.p@data$Fishnet_ID)
#Select categorical covariates
potato.cat<-potato.p %>% select(COUNTRY.x,PotatoGCO, FISHNET_ID,Latitude,Longitude)
potato.final<-merge(potato.catagory,potato.scaled, by="FISHNET_ID")
potato.final$PotatoGCO<-as.factor(potato.final$PotatoGCO)
potato.final$potatoBinaryGCO<- ifelse(potato.final$PotatoGCO == 'Inside', 1, 0) #set binary dummy
potato.final<-na.omit(potato.final)
potato.coords<-potato.final %>% select(Latitude,Longitude) #select coordinates
summary(potato.final)
dim(potato.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(potato.coords$Latitude, potato.coords$Longitude)))
dim(potato.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-potato.coords$Latitude
y<-potato.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
potato.data<-potato.final

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 3,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = potato.data, 
  dependent.variable.name = "mean_potat.x",
  predictor.variable.names = c("AET_mean","Potato_Fe", "Pesticide", "GDP_Mean","potatoBinaryGCO"),
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
potato.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster,
 max.spatial.predictors = 10
)
print(potato.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

potato.resid<-data.frame(potato.spatial$residuals$values) #extract residuals
dim(potato.resid)
#Plot residuals
ggplot(potato.data, aes(x=Longitude, y=Latitude, color=potato.resid$potato.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
potato.curves.df <- spatialRF::get_response_curves(potato.spatial,variables = c("AET_mean","Potato_Fe", "Pesticide", "GDP_Mean","potatoBinaryGCO"))


write.csv(potato.spatial$predictions, file="Crop_GCO/Global/ModelOutputs/potato_pred.csv")
write.csv(potato.spatial$residuals$values, file="Crop_GCO/Global/ModelOutputs/potato_resid.csv")
write.csv(potato.spatial$variable.importance, file="Crop_GCO/Global/ModelOutputs/potato_VarImp.csv")
write.csv(potato.spatial$performance, file="Crop_GCO/Global/ModelOutputs/potato_preformance.csv")
write.csv(potato.curves.df, file="Crop_GCO/Global/ModelOutputs/potato_curves.csv")

