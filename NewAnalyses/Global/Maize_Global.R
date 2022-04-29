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
maize.model<-na.omit(read.csv(file="Crop_GCO/Global/Models/Maize_RF.csv", header=T))
mapping<-read.csv(file="Crop_GCO/Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, maize.model, by='Fishnet_ID')
maize<-sp::merge(m, mapping, by='Fishnet_ID')
maize.p = subset(maize, maize_HgHa > 0) #selecting values > 0 
#Transform numerical covariates

maize.num<-maize.p@data %>% select(maize_HgHa, AET_mean, maize_Fertilizer, Pesticide,GDP_Mean) 
maize.scaled<-data.frame(log10(maize.num+1)) #psuedo count to avoid inf
maize.scaled$FISHNET_ID<-as.numeric(maize.p@data$Fishnet_ID)
#Select categorical covariates
maize.cat<-maize.p %>% select(COUNTRY.x,MaizeGCO, FISHNET_ID,Latitude,Longitude)

maize.final<-merge(maize.cat,maize.scaled, by="FISHNET_ID")
maize.final$maizeGCO<-as.factor(maize.final$MaizeGCO)
maize.final$maizeBinaryGCO<- ifelse(maize.final$maizeGCO == 'Inside', 1, 0) #set binary dummy
maize.final<-na.omit(maize.final)
maize.coords<-maize.final %>% select(Latitude,Longitude) #select coordinates
summary(maize.final)
dim(maize.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(maize.coords$Latitude, maize.coords$Longitude)))
dim(maize.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-maize.coords$Latitude
y<-maize.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
maize.data<-data.frame(maize.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 10,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = maize.data, 
  dependent.variable.name = "maize_HgHa",
  predictor.variable.names = c("AET_mean","maize_Fertilizer", "Pesticide", "GDP_Mean","maizeBinaryGCO"),
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
maize.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster
)
print(maize.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

maize.resid<-data.frame(maize.spatial$residuals$values) #extract residuals
dim(maize.resid)
#Plot residuals
ggplot(maize.data, aes(x=Longitude, y=Latitude, color=maize.resid$maize.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
maize.curves.df <- spatialRF::get_response_curves(maize.spatial,variables = c("AET_mean","maize_Fertilizer", "Pesticide", "GDP_Mean","maizeBinaryGCO"))
write.csv(maize.spatial$predictions, file="Crop_GCO/Global/ModelOutputs/maize_pred.csv")
write.csv(maize.spatial$residuals$values, file="Crop_GCO/Global/ModelOutputs/maize_resid.csv")
write.csv(maize.spatial$variable.importance, file="Crop_GCO/Global/ModelOutputs/maize_VarImp.csv")
write.csv(maize.spatial$performance, file="Crop_GCO/Global/ModelOutputs/maize_preformance.csv")
write.csv(maize.curves.df, file="Crop_GCO/Global/ModelOutputs/maize_curves.csv")

