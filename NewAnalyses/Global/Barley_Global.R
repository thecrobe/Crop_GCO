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
barley.model<-na.omit(read.csv(file="Crop_GCO/Global/Models/Barley_RF.csv", header=T))
mapping<-read.csv(file="Crop_GCO/Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, barley.model, by='Fishnet_ID')
barley<-sp::merge(m, mapping, by='Fishnet_ID')
barley.p = subset(barley, barley_HgHa > 0) #selecting values > 0 
#Transform numerical covariates

barley.num<-barley.p@data %>% select(barley_HgHa, AET_mean, Barley_Fertilizer, Pesticide,GDP_Mean) 
barley.scaled<-data.frame(log10(barley.num+1)) #psuedo count to avoid inf
barley.scaled$FISHNET_ID<-as.numeric(barley.p@data$Fishnet_ID)
#Select categorical covariates
barley.cat<-barley.p %>% select(COUNTRY.x,BarleyGCO, FISHNET_ID,Latitude,Longitude)

barley.final<-merge(barley.cat,barley.scaled, by="FISHNET_ID")
barley.final$barleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$barleyBinaryGCO<- ifelse(barley.final$barleyGCO == 'Inside', 1, 0) #set binary dummy
barley.final<-na.omit(barley.final)
barley.coords<-barley.final %>% select(Latitude,Longitude) #select coordinates
summary(barley.final)
dim(barley.final) #make sure dims are correct 

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(barley.coords$Latitude, barley.coords$Longitude)))
dim(barley.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-barley.coords$Latitude
y<-barley.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
barley.data<-data.frame(barley.final)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 3,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = barley.data, 
  dependent.variable.name = "barley_HgHa",
  predictor.variable.names = c("AET_mean","Barley_Fertilizer", "Pesticide", "GDP_Mean","barleyBinaryGCO"),
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
barley.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster
)
print(barley.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

barley.resid<-data.frame(barley.spatial$residuals$values) #extract residuals
dim(barley.resid)
#Plot residuals
ggplot(barley.data, aes(x=Longitude, y=Latitude, color=barley.resid$barley.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
barley.curves.df <- spatialRF::get_response_curves(barley.spatial,variables = c("AET_mean","Barley_Fertilizer", "Pesticide", "GDP_Mean","barleyBinaryGCO"))
write.csv(barley.spatial$predictions, file="Crop_GCO/Global/ModelOutputs/barley_pred.csv")
write.csv(barley.spatial$residuals$values, file="Crop_GCO/Global/ModelOutputs/barley_resid.csv")
write.csv(barley.spatial$variable.importance, file="Crop_GCO/Global/ModelOutputs/barley_VarImp.csv")
write.csv(barley.spatial$performance, file="Crop_GCO/Global/ModelOutputs/barley_preformance.csv")
write.csv(barley.curves.df, file="Crop_GCO/Global/ModelOutputs/barley_curves.csv")

