library(spatialRF)
library(spdep)
library(spdplyr)
library(rgdal)
library(ggplot2)
library(ranger)

#Read In
#Fishnet- pixel = 100km^2
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
summary(fishnet)
#Set Projection
proj4string(fishnet) <- CRS("+init=epsg:3786") 

#Some covariates
soybean_model<-na.omit(read.csv(file="Models/Soybean_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, soybean_model, by='Fishnet_ID')
soybean<-sp::merge(m, mapping, by='Fishnet_ID')
soybean.p = subset(soybean, mean_soybe > 0) #selecting values > 0 

#Transform numerical covariates
soybean.num<-soybean.p@data %>% select(mean_soybe, AET_mean, soybean_Fertilizer, Pesticide,GDP_Mean) 
soybean.scaled<-data.frame(log10(soybean.num+1)) #psuedo count to avoid inf
soybean.scaled$FISHNET_ID<-as.numeric(soybean.p@data$Fishnet_ID)

#Select categorical covariates
soybean.cat<-soybean.p %>% select(COUNTRY.x,SoybeanGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-soybean.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
soybean.catagory<-sp::merge(soybean.cat,county.pivot, by="COUNTRY.x")
soybean.catagory<-soybean.catagory@data


########## Asia #########
#Select Continent
soybean.comb<-merge(soybean.catagory,soybean.scaled, by="FISHNET_ID")
soybean.final<-(filter((soybean.comb),CONTINENT == "Asia"))
soybean.final<-na.omit(soybean.final)

summary(soybean.final)
soybean.final$SoybeanGCO<-as.factor(soybean.final$SoybeanGCO)
soybean.final$COUNTRY.x<-as.factor(soybean.final$COUNTRY.x)
soybean.coords<-soybean.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
soybean.data<-soybean.final

model.non.spatial <- spatialRF::rf(
  data = soybean.data,
  dependent.variable.name = "mean_soybe",
  predictor.variable.names = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"),
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


write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Soybean_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Soybean_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(soybean.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Soybean_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "soybean_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Soybean_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
soybean.comb<-merge(soybean.catagory,soybean.scaled, by="FISHNET_ID")
soybean.final<-(filter((soybean.comb),CONTINENT == "Africa"))
soybean.final<-na.omit(soybean.final)

summary(soybean.final)
soybean.final$SoybeanGCO<-as.factor(soybean.final$SoybeanGCO)
soybean.final$COUNTRY.x<-as.factor(soybean.final$COUNTRY.x)
soybean.coords<-soybean.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
soybean.data<-soybean.final

model.non.spatial <- spatialRF::rf(
  data = soybean.data,
  dependent.variable.name = "mean_soybe",
  predictor.variable.names = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Soybean_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Soybean_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(soybean.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Soybean_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "soybean_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Soybean_Africa_RespCurv.csv")

########## North America #########
#Select Continent
soybean.comb<-merge(soybean.catagory,soybean.scaled, by="FISHNET_ID")
soybean.final<-(filter((soybean.comb),CONTINENT == "North America"))
soybean.final<-na.omit(soybean.final)

summary(soybean.final)
soybean.final$SoybeanGCO<-as.factor(soybean.final$SoybeanGCO)
soybean.final$COUNTRY.x<-as.factor(soybean.final$COUNTRY.x)
soybean.coords<-soybean.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
soybean.data<-soybean.final

model.non.spatial <- spatialRF::rf(
  data = soybean.data,
  dependent.variable.name = "mean_soybe",
  predictor.variable.names = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Soybean_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Soybean_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(soybean.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Soybean_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "soybean_Fertilizer", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "soybean_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Soybean_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
soybean.comb<-merge(soybean.catagory,soybean.scaled, by="FISHNET_ID")
soybean.final<-(filter((soybean.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
soybean.final<-na.omit(soybean.final)

summary(soybean.final)
soybean.final$SoybeanGCO<-as.factor(soybean.final$SoybeanGCO)
soybean.final$COUNTRY.x<-as.factor(soybean.final$COUNTRY.x)
soybean.coords<-soybean.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
soybean.data<-soybean.final

model.non.spatial <- spatialRF::rf(
  data = soybean.data,
  dependent.variable.name = "mean_soybe",
  predictor.variable.names = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Soybean_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Soybean_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(soybean.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Soybean_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "soybean_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Soybean_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
soybean.comb<-merge(soybean.catagory,soybean.scaled, by="FISHNET_ID")
soybean.final<-(filter((soybean.comb),CONTINENT == "Europe"))
soybean.final<-na.omit(soybean.final)

summary(soybean.final)
soybean.final$SoybeanGCO<-as.factor(soybean.final$SoybeanGCO)
soybean.final$COUNTRY.x<-as.factor(soybean.final$COUNTRY.x)
soybean.coords<-soybean.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
soybean.data<-soybean.final

model.non.spatial <- spatialRF::rf(
  data = soybean.data,
  dependent.variable.name = "mean_soybe",
  predictor.variable.names = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Soybean_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Soybean_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(soybean.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Soybean_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "soybean_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Soybean_Europe_RespCurv.csv")

########## South America #########
#Select Continent
soybean.comb<-merge(soybean.catagory,soybean.scaled, by="FISHNET_ID")
soybean.final<-(filter((soybean.comb),CONTINENT == "South America"))
soybean.final<-na.omit(soybean.final)

summary(soybean.final)
soybean.final$SoybeanGCO<-as.factor(soybean.final$SoybeanGCO)
soybean.final$COUNTRY.x<-as.factor(soybean.final$COUNTRY.x)
soybean.coords<-soybean.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
soybean.data<-soybean.final

model.non.spatial <- spatialRF::rf(
  data = soybean.data,
  dependent.variable.name = "mean_soybe",
  predictor.variable.names = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Soybean_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Soybean_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(soybean.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Soybean_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "soybean_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "soybean_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Soybean_South America_RespCurv.csv")

