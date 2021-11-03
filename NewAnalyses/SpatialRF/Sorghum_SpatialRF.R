library(spatialRF)
library(spdep)
library(spdplyr)
library(rgdal)
library(ggplot2)
library(ranger)



#Read In
#Fishnet- pixel = 100km^2
fishnet<- readOGR(dsn= "GIS/", layer="Fishnetield_NoAntarctica")
summary(fishnet)
#Set Projection
proj4string(fishnet) <- CRS("+init=epsg:3786") 

#Some covariates
sorgh_model<-na.omit(read.csv(file="Models/Sorghum_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, sorgh_model, by='Fishnet_ID')
sorghum<-sp::merge(m, mapping, by='Fishnet_ID')
sorghum.p = subset(sorghum, mean_sorgh > 0) #selecting values > 0 

#Transform numerical covariates
sorghum.num<-sorghum.p@data %>% select(mean_sorgh, Evapotranspiraton, sorghum_Fertilizer, Pesticide,GDP) 
sorghum.scaled<-data.frame(log10(sorghum.num+1)) #psuedo count to avoid inf
sorghum.scaled$FISHNET_ID<-as.numeric(sorghum.p@data$Fishnet_ID)

#Select categorical covariates
sorghum.cat<-sorghum.p %>% select(COUNTRY.x,SorghumGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-sorghum.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
sorghum.catagory<-sp::merge(sorghum.cat,county.pivot, by="COUNTRY.x")
sorghum.catagory<-sorghum.catagory@data


########## Asia #########
#Select Continent
sorghum.comb<-merge(sorghum.catagory,sorghum.scaled, by="FISHNET_ID")
sorghum.final<-(filter((sorghum.comb),CONTINENT == "Asia"))
sorghum.final<-na.omit(sorghum.final)

summary(sorghum.final)
sorghum.final$SorghumGCO<-as.factor(sorghum.final$SorghumGCO)
sorghum.final$COUNTRY.x<-as.factor(sorghum.final$COUNTRY.x)
sorghum.coords<-sorghum.final %>% select(Latitude,Longitude) #select coordinates

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

model.non.spatial <- spatialRF::rf(
  data = sorghum.data,
  dependent.variable.name = "mean_sorgh",
  predictor.variable.names = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 12,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 12
)


write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sorghum_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sorghum_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sorghum.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")


#Variable Importance
spatialRF::plot_importance( model.spatial,verbose = FALSE) + ggtitle("Spatial model")
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sorghum_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "Evapotranspiraton", "sorghum_Fertilizer", "GDP"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sorghum_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
sorghum.comb<-merge(sorghum.catagory,sorghum.scaled, by="FISHNET_ID")
sorghum.final<-(filter((sorghum.comb),CONTINENT == "Africa"))
sorghum.final<-na.omit(sorghum.final)

summary(sorghum.final)
sorghum.final$SorghumGCO<-as.factor(sorghum.final$SorghumGCO)
sorghum.final$COUNTRY.x<-as.factor(sorghum.final$COUNTRY.x)
sorghum.coords<-sorghum.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sorghum.data<-sorghum.final

model.non.spatial <- spatialRF::rf(
  data = sorghum.data,
  dependent.variable.name = "mean_sorgh",
  predictor.variable.names = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 12,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 12
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sorghum_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sorghum_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sorghum.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")


#Variable Importance
spatialRF::plot_importance( model.spatial,verbose = FALSE) + ggtitle("Spatial model")
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sorghum_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "Evapotranspiraton", "sorghum_Fertilizer", "GDP"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sorghum_Africa_RespCurv.csv")

########## North America #########
#Select Continent
sorghum.comb<-merge(sorghum.catagory,sorghum.scaled, by="FISHNET_ID")
sorghum.final<-(filter((sorghum.comb),CONTINENT == "North America"))
sorghum.final<-na.omit(sorghum.final)

summary(sorghum.final)
sorghum.final$SorghumGCO<-as.factor(sorghum.final$SorghumGCO)
sorghum.final$COUNTRY.x<-as.factor(sorghum.final$COUNTRY.x)
sorghum.coords<-sorghum.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sorghum.data<-sorghum.final

model.non.spatial <- spatialRF::rf(
  data = sorghum.data,
  dependent.variable.name = "mean_sorgh",
  predictor.variable.names = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 12,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 12
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sorghum_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sorghum_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sorghum.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")


#Variable Importance
spatialRF::plot_importance(model.spatial,verbose = FALSE) + ggtitle("Spatial model")
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sorghum_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "Evapotranspiraton", "sorghum_Fertilizer", "GDP"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "Evapotranspiraton", "sorghum_Fertilizer", "GDP"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sorghum_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
sorghum.comb<-merge(sorghum.catagory,sorghum.scaled, by="FISHNET_ID")
sorghum.final<-(filter((sorghum.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
sorghum.final<-na.omit(sorghum.final)

summary(sorghum.final)
sorghum.final$SorghumGCO<-as.factor(sorghum.final$SorghumGCO)
sorghum.final$COUNTRY.x<-as.factor(sorghum.final$COUNTRY.x)
sorghum.coords<-sorghum.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sorghum.data<-sorghum.final

model.non.spatial <- spatialRF::rf(
  data = sorghum.data,
  dependent.variable.name = "mean_sorgh",
  predictor.variable.names = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 12,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 12
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sorghum_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sorghum_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sorghum.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")


#Variable Importance
spatialRF::plot_importance( model.spatial,verbose = FALSE) + ggtitle("Spatial model")
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sorghum_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "Evapotranspiraton", "sorghum_Fertilizer", "GDP"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sorghum_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
sorghum.comb<-merge(sorghum.catagory,sorghum.scaled, by="FISHNET_ID")
sorghum.final<-(filter((sorghum.comb),CONTINENT == "Europe"))
sorghum.final<-na.omit(sorghum.final)

summary(sorghum.final)
sorghum.final$SorghumGCO<-as.factor(sorghum.final$SorghumGCO)
sorghum.final$COUNTRY.x<-as.factor(sorghum.final$COUNTRY.x)
sorghum.coords<-sorghum.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sorghum.data<-sorghum.final

model.non.spatial <- spatialRF::rf(
  data = sorghum.data,
  dependent.variable.name = "mean_sorgh",
  predictor.variable.names = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 12,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 12
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sorghum_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sorghum_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sorghum.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")


#Variable Importance
spatialRF::plot_importance(model.spatial,verbose = FALSE) + ggtitle("Spatial model")
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sorghum_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "Evapotranspiraton", "sorghum_Fertilizer", "GDP"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sorghum_Europe_RespCurv.csv")

########## South America #########
#Select Continent
sorghum.comb<-merge(sorghum.catagory,sorghum.scaled, by="FISHNET_ID")
sorghum.final<-(filter((sorghum.comb),CONTINENT == "South America"))
sorghum.final<-na.omit(sorghum.final)

summary(sorghum.final)
sorghum.final$SorghumGCO<-as.factor(sorghum.final$SorghumGCO)
sorghum.final$COUNTRY.x<-as.factor(sorghum.final$COUNTRY.x)
sorghum.coords<-sorghum.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sorghum.data<-sorghum.final

model.non.spatial <- spatialRF::rf(
  data = sorghum.data,
  dependent.variable.name = "mean_sorgh",
  predictor.variable.names = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 12,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 12
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sorghum_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sorghum_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sorghum.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")


#Variable Importance
spatialRF::plot_importance( model.spatial,verbose = FALSE) + ggtitle("Spatial model")
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sorghum_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("Evapotranspiraton", "sorghum_Fertilizer", "Pesticide", "GDP"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "Evapotranspiraton", "sorghum_Fertilizer", "GDP"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sorghum_South America_RespCurv.csv")

