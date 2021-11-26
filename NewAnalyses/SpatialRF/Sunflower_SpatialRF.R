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
sunflower_model<-na.omit(read.csv(file="Models/Sunflower_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, sunflower_model, by='Fishnet_ID')
sunflower<-sp::merge(m, mapping, by='Fishnet_ID')
sunflower.p = subset(sunflower, mean_sunfl > 0) #selecting values > 0 

#Transform numerical covariates
sunflower.num<-sunflower.p@data %>% select(mean_sunfl, AET_mean, sunflower_Fertilizer, Pesticide,GDP_Mean) 
sunflower.scaled<-data.frame(log10(sunflower.num+1)) #psuedo count to avoid inf
sunflower.scaled$FISHNET_ID<-as.numeric(sunflower.p@data$Fishnet_ID)

#Select categorical covariates
sunflower.cat<-sunflower.p %>% select(COUNTRY.x,SunflowerGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-sunflower.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
sunflower.catagory<-sp::merge(sunflower.cat,county.pivot, by="COUNTRY.x")
sunflower.catagory<-sunflower.catagory@data


########## Asia #########
#Select Continent
sunflower.comb<-merge(sunflower.catagory,sunflower.scaled, by="FISHNET_ID")
sunflower.final<-(filter((sunflower.comb),CONTINENT == "Asia"))
sunflower.final<-na.omit(sunflower.final)

summary(sunflower.final)
sunflower.final$SunflowerGCO<-as.factor(sunflower.final$SunflowerGCO)
sunflower.final$COUNTRY.x<-as.factor(sunflower.final$COUNTRY.x)
sunflower.coords<-sunflower.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(sunflower.coords$Latitude, sunflower.coords$Longitude)))
dim(sunflower.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-sunflower.coords$Latitude
y<-sunflower.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sunflower.data<-sunflower.final

model.non.spatial <- spatialRF::rf(
  data = sunflower.data,
  dependent.variable.name = "mean_sunfl",
  predictor.variable.names = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"),
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


write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sunflower_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sunflower_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sunflower.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sunflower_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "sunflower_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sunflower_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
sunflower.comb<-merge(sunflower.catagory,sunflower.scaled, by="FISHNET_ID")
sunflower.final<-(filter((sunflower.comb),CONTINENT == "Africa"))
sunflower.final<-na.omit(sunflower.final)

summary(sunflower.final)
sunflower.final$SunflowerGCO<-as.factor(sunflower.final$SunflowerGCO)
sunflower.final$COUNTRY.x<-as.factor(sunflower.final$COUNTRY.x)
sunflower.coords<-sunflower.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(sunflower.coords$Latitude, sunflower.coords$Longitude)))
dim(sunflower.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-sunflower.coords$Latitude
y<-sunflower.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sunflower.data<-sunflower.final

model.non.spatial <- spatialRF::rf(
  data = sunflower.data,
  dependent.variable.name = "mean_sunfl",
  predictor.variable.names = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sunflower_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sunflower_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sunflower.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sunflower_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "sunflower_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sunflower_Africa_RespCurv.csv")

########## North America #########
#Select Continent
sunflower.comb<-merge(sunflower.catagory,sunflower.scaled, by="FISHNET_ID")
sunflower.final<-(filter((sunflower.comb),CONTINENT == "North America"))
sunflower.final<-na.omit(sunflower.final)

summary(sunflower.final)
sunflower.final$SunflowerGCO<-as.factor(sunflower.final$SunflowerGCO)
sunflower.final$COUNTRY.x<-as.factor(sunflower.final$COUNTRY.x)
sunflower.coords<-sunflower.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(sunflower.coords$Latitude, sunflower.coords$Longitude)))
dim(sunflower.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-sunflower.coords$Latitude
y<-sunflower.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sunflower.data<-sunflower.final

model.non.spatial <- spatialRF::rf(
  data = sunflower.data,
  dependent.variable.name = "mean_sunfl",
  predictor.variable.names = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sunflower_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sunflower_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sunflower.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sunflower_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "sunflower_Fertilizer", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "sunflower_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sunflower_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
sunflower.comb<-merge(sunflower.catagory,sunflower.scaled, by="FISHNET_ID")
sunflower.final<-(filter((sunflower.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
sunflower.final<-na.omit(sunflower.final)

summary(sunflower.final)
sunflower.final$SunflowerGCO<-as.factor(sunflower.final$SunflowerGCO)
sunflower.final$COUNTRY.x<-as.factor(sunflower.final$COUNTRY.x)
sunflower.coords<-sunflower.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(sunflower.coords$Latitude, sunflower.coords$Longitude)))
dim(sunflower.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-sunflower.coords$Latitude
y<-sunflower.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sunflower.data<-sunflower.final

model.non.spatial <- spatialRF::rf(
  data = sunflower.data,
  dependent.variable.name = "mean_sunfl",
  predictor.variable.names = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sunflower_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sunflower_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sunflower.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sunflower_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "sunflower_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sunflower_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
sunflower.comb<-merge(sunflower.catagory,sunflower.scaled, by="FISHNET_ID")
sunflower.final<-(filter((sunflower.comb),CONTINENT == "Europe"))
sunflower.final<-na.omit(sunflower.final)

summary(sunflower.final)
sunflower.final$SunflowerGCO<-as.factor(sunflower.final$SunflowerGCO)
sunflower.final$COUNTRY.x<-as.factor(sunflower.final$COUNTRY.x)
sunflower.coords<-sunflower.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(sunflower.coords$Latitude, sunflower.coords$Longitude)))
dim(sunflower.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-sunflower.coords$Latitude
y<-sunflower.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sunflower.data<-sunflower.final

model.non.spatial <- spatialRF::rf(
  data = sunflower.data,
  dependent.variable.name = "mean_sunfl",
  predictor.variable.names = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sunflower_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sunflower_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sunflower.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sunflower_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "sunflower_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sunflower_Europe_RespCurv.csv")

########## South America #########
#Select Continent
sunflower.comb<-merge(sunflower.catagory,sunflower.scaled, by="FISHNET_ID")
sunflower.final<-(filter((sunflower.comb),CONTINENT == "South America"))
sunflower.final<-na.omit(sunflower.final)

summary(sunflower.final)
sunflower.final$SunflowerGCO<-as.factor(sunflower.final$SunflowerGCO)
sunflower.final$COUNTRY.x<-as.factor(sunflower.final$COUNTRY.x)
sunflower.coords<-sunflower.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(sunflower.coords$Latitude, sunflower.coords$Longitude)))
dim(sunflower.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-sunflower.coords$Latitude
y<-sunflower.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
sunflower.data<-sunflower.final

model.non.spatial <- spatialRF::rf(
  data = sunflower.data,
  dependent.variable.name = "mean_sunfl",
  predictor.variable.names = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Sunflower_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Sunflower_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(sunflower.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Sunflower_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "sunflower_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "sunflower_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Sunflower_South America_RespCurv.csv")

