
library(spatialRF)
library(spdep)
library(spdplyr)
library(rgdal)
library(ggplot2)
library(ranger)

read

#Read In
#Fishnet- pixel = 100km^2
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
summary(fishnet)
#Set Projection
proj4string(fishnet) <- CRS("+init=epsg:3786") 

#Some covariates
ryemodel<-na.omit(read.csv(file="Models/Rye_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, ryemodel, by='Fishnet_ID')
rye<-sp::merge(m, mapping, by='Fishnet_ID')
rye.p = subset(rye, mean_rye_Y > 0) #selecting values > 0 

#Transform numerical covariates
rye.num<-rye.p@data %>% select(mean_rye_Y, AET_mean, rye_Fertilizer, Pesticide,GDP_Mean) 
rye.scaled<-data.frame(log10(rye.num+1)) #psuedo count to avoid inf
rye.scaled$FISHNET_ID<-as.numeric(rye.p@data$Fishnet_ID)

#Select categorical covariates
rye.cat<-rye.p %>% select(COUNTRY.x,RyeGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-rye.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
rye.catagory<-sp::merge(rye.cat,county.pivot, by="COUNTRY.x")
rye.catagory<-rye.catagory@data


########## Asia #########
#Select Continent
rye.comb<-merge(rye.catagory,rye.scaled, by="FISHNET_ID")
rye.final<-(filter((rye.comb),CONTINENT == "Asia"))
rye.final<-na.omit(rye.final)

summary(rye.final)
rye.final$RyeGCO<-as.factor(rye.final$RyeGCO)
rye.final$COUNTRY.x<-as.factor(rye.final$COUNTRY.x)
rye.coords<-rye.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rye.coords$Latitude, rye.coords$Longitude)))
dim(rye.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rye.coords$Latitude
y<-rye.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rye.data<-rye.final

model.non.spatial <- spatialRF::rf(
  data = rye.data,
  dependent.variable.name = "mean_rye_Y",
  predictor.variable.names = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"),
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


write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rye_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rye_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rye.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rye_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rye_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rye_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
rye.comb<-merge(rye.catagory,rye.scaled, by="FISHNET_ID")
rye.final<-(filter((rye.comb),CONTINENT == "Africa"))
rye.final<-na.omit(rye.final)

summary(rye.final)
rye.final$RyeGCO<-as.factor(rye.final$RyeGCO)
rye.final$COUNTRY.x<-as.factor(rye.final$COUNTRY.x)
rye.coords<-rye.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rye.coords$Latitude, rye.coords$Longitude)))
dim(rye.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rye.coords$Latitude
y<-rye.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rye.data<-rye.final

model.non.spatial <- spatialRF::rf(
  data = rye.data,
  dependent.variable.name = "mean_rye_Y",
  predictor.variable.names = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rye_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rye_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rye.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rye_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rye_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rye_Africa_RespCurv.csv")

########## North America #########
#Select Continent
rye.comb<-merge(rye.catagory,rye.scaled, by="FISHNET_ID")
rye.final<-(filter((rye.comb),CONTINENT == "North America"))
rye.final<-na.omit(rye.final)

summary(rye.final)
rye.final$RyeGCO<-as.factor(rye.final$RyeGCO)
rye.final$COUNTRY.x<-as.factor(rye.final$COUNTRY.x)
rye.coords<-rye.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rye.coords$Latitude, rye.coords$Longitude)))
dim(rye.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rye.coords$Latitude
y<-rye.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rye.data<-rye.final

model.non.spatial <- spatialRF::rf(
  data = rye.data,
  dependent.variable.name = "mean_rye_Y",
  predictor.variable.names = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rye_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rye_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rye.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rye_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "rye_Fertilizer", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rye_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rye_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
rye.comb<-merge(rye.catagory,rye.scaled, by="FISHNET_ID")
rye.final<-(filter((rye.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
rye.final<-na.omit(rye.final)

summary(rye.final)
rye.final$RyeGCO<-as.factor(rye.final$RyeGCO)
rye.final$COUNTRY.x<-as.factor(rye.final$COUNTRY.x)
rye.coords<-rye.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rye.coords$Latitude, rye.coords$Longitude)))
dim(rye.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rye.coords$Latitude
y<-rye.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rye.data<-rye.final

model.non.spatial <- spatialRF::rf(
  data = rye.data,
  dependent.variable.name = "mean_rye_Y",
  predictor.variable.names = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rye_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rye_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rye.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rye_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rye_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rye_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
rye.comb<-merge(rye.catagory,rye.scaled, by="FISHNET_ID")
rye.final<-(filter((rye.comb),CONTINENT == "Europe"))
rye.final<-na.omit(rye.final)

summary(rye.final)
rye.final$RyeGCO<-as.factor(rye.final$RyeGCO)
rye.final$COUNTRY.x<-as.factor(rye.final$COUNTRY.x)
rye.coords<-rye.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rye.coords$Latitude, rye.coords$Longitude)))
dim(rye.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rye.coords$Latitude
y<-rye.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rye.data<-rye.final

model.non.spatial <- spatialRF::rf(
  data = rye.data,
  dependent.variable.name = "mean_rye_Y",
  predictor.variable.names = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rye_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rye_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rye.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rye_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rye_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rye_Europe_RespCurv.csv")

########## South America #########
#Select Continent
rye.comb<-merge(rye.catagory,rye.scaled, by="FISHNET_ID")
rye.final<-(filter((rye.comb),CONTINENT == "South America"))
rye.final<-na.omit(rye.final)

summary(rye.final)
rye.final$RyeGCO<-as.factor(rye.final$RyeGCO)
rye.final$COUNTRY.x<-as.factor(rye.final$COUNTRY.x)
rye.coords<-rye.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rye.coords$Latitude, rye.coords$Longitude)))
dim(rye.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rye.coords$Latitude
y<-rye.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rye.data<-rye.final

model.non.spatial <- spatialRF::rf(
  data = rye.data,
  dependent.variable.name = "mean_rye_Y",
  predictor.variable.names = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rye_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rye_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rye.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rye_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rye_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rye_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rye_South America_RespCurv.csv")

