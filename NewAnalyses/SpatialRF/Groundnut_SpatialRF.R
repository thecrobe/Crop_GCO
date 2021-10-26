library(spatialRF)
library(spdep)
library(spdplyr)
library(rgdal)
library(ggplot2)
library(ranger)
library(dplyr)

#Read In
#Fishnet- pixel = 100km^2
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
summary(fishnet)
#Set Projection
proj4string(fishnet) <- CRS("+init=epsg:3786") 

#Some covariates
groundnutmodel<-na.omit(read.csv(file="Models/Groundnut_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, groundnutmodel, by='Fishnet_ID')
groundnut<-sp::merge(m, mapping, by='Fishnet_ID')
groundnut.p = subset(groundnut, groundnut_HgHa > 0) #selecting values > 0 

#Transform numerical covariates
groundnut.num<-groundnut.p@data %>% select(groundnut_HgHa, AET_mean, groundnut_Fertilizer, Pesticide,GDP_Mean) 
groundnut.scaled<-data.frame(log10(groundnut.num+1)) #psuedo count to avoid inf
groundnut.scaled$FISHNET_ID<-as.numeric(groundnut.p@data$Fishnet_ID)

#Select categorical covariates
groundnut.cat<-groundnut.p %>% select(COUNTRY.x,GroundnutGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-groundnut.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
groundnut.catagory<-sp::merge(groundnut.cat,county.pivot, by="COUNTRY.x")
groundnut.catagory<-groundnut.catagory@data



########## Asia #########
#Select Continent
groundnut.comb<-merge(groundnut.catagory,groundnut.scaled, by="FISHNET_ID")
groundnut.final<-(filter((groundnut.comb),CONTINENT == "Asia"))
groundnut.final<-na.omit(groundnut.final)

summary(groundnut.final)
groundnut.final$GroundnutGCO<-as.factor(groundnut.final$GroundnutGCO)
groundnut.final$COUNTRY.x<-as.factor(groundnut.final$COUNTRY.x)
groundnut.coords<-groundnut.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(groundnut.coords$Latitude, groundnut.coords$Longitude)))
dim(groundnut.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-groundnut.coords$Latitude
y<-groundnut.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
groundnut.data<-groundnut.final

model.non.spatial <- spatialRF::rf(
  data = groundnut.data,
  dependent.variable.name = "groundnut_HgHa",
  predictor.variable.names = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 4,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Groundnut_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Groundnut_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(groundnut.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Groundnut_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "groundnut_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Groundnut_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
groundnut.comb<-merge(groundnut.catagory,groundnut.scaled, by="FISHNET_ID")
groundnut.final<-(filter((groundnut.comb),CONTINENT == "Africa"))
groundnut.final<-na.omit(groundnut.final)

summary(groundnut.final)
groundnut.final$GroundnutGCO<-as.factor(groundnut.final$GroundnutGCO)
groundnut.final$COUNTRY.x<-as.factor(groundnut.final$COUNTRY.x)
groundnut.coords<-groundnut.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(groundnut.coords$Latitude, groundnut.coords$Longitude)))
dim(groundnut.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-groundnut.coords$Latitude
y<-groundnut.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
groundnut.data<-groundnut.final

model.non.spatial <- spatialRF::rf(
  data = groundnut.data,
  dependent.variable.name = "groundnut_HgHa",
  predictor.variable.names = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 4,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/SpatialRF/Groundnut_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Groundnut_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(groundnut.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Groundnut_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "groundnut_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Groundnut_Africa_RespCurv.csv")


########## North America #########
#Select Continent
groundnut.comb<-merge(groundnut.catagory,groundnut.scaled, by="FISHNET_ID")
groundnut.final<-(filter((groundnut.comb),CONTINENT == "North America"))
groundnut.final<-na.omit(groundnut.final)

summary(groundnut.final)
groundnut.final$GroundnutGCO<-as.factor(groundnut.final$GroundnutGCO)
groundnut.final$COUNTRY.x<-as.factor(groundnut.final$COUNTRY.x)
groundnut.coords<-groundnut.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(groundnut.coords$Latitude, groundnut.coords$Longitude)))
dim(groundnut.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-groundnut.coords$Latitude
y<-groundnut.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
groundnut.data<-groundnut.final

model.non.spatial <- spatialRF::rf(
  data = groundnut.data,
  dependent.variable.name = "groundnut_HgHa",
  predictor.variable.names = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 4,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Groundnut_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Groundnut_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(groundnut.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Groundnut_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "groundnut_Fertilizer", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "groundnut_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Groundnut_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
groundnut.comb<-merge(groundnut.catagory,groundnut.scaled, by="FISHNET_ID")
groundnut.final<-(filter((groundnut.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
groundnut.final<-na.omit(groundnut.final)

summary(groundnut.final)
groundnut.final$GroundnutGCO<-as.factor(groundnut.final$GroundnutGCO)
groundnut.final$COUNTRY.x<-as.factor(groundnut.final$COUNTRY.x)
groundnut.coords<-groundnut.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(groundnut.coords$Latitude, groundnut.coords$Longitude)))
dim(groundnut.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-groundnut.coords$Latitude
y<-groundnut.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
groundnut.data<-groundnut.final

model.non.spatial <- spatialRF::rf(
  data = groundnut.data,
  dependent.variable.name = "groundnut_HgHa",
  predictor.variable.names = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 4,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Groundnut_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Groundnut_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(groundnut.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Groundnut_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "groundnut_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Groundnut_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
groundnut.comb<-merge(groundnut.catagory,groundnut.scaled, by="FISHNET_ID")
groundnut.final<-(filter((groundnut.comb),CONTINENT == "Europe"))
groundnut.final<-na.omit(groundnut.final)

summary(groundnut.final)
groundnut.final$GroundnutGCO<-as.factor(groundnut.final$GroundnutGCO)
groundnut.final$COUNTRY.x<-as.factor(groundnut.final$COUNTRY.x)
groundnut.coords<-groundnut.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(groundnut.coords$Latitude, groundnut.coords$Longitude)))
dim(groundnut.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-groundnut.coords$Latitude
y<-groundnut.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
groundnut.data<-groundnut.final

model.non.spatial <- spatialRF::rf(
  data = groundnut.data,
  dependent.variable.name = "groundnut_HgHa",
  predictor.variable.names = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 4,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Groundnut_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Groundnut_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(groundnut.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Groundnut_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "groundnut_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Groundnut_Europe_RespCurv.csv")

########## South America #########
#Select Continent
groundnut.comb<-merge(groundnut.catagory,groundnut.scaled, by="FISHNET_ID")
groundnut.final<-(filter((groundnut.comb),CONTINENT == "South America"))
groundnut.final<-na.omit(groundnut.final)

summary(groundnut.final)
groundnut.final$GroundnutGCO<-as.factor(groundnut.final$GroundnutGCO)
groundnut.final$COUNTRY.x<-as.factor(groundnut.final$COUNTRY.x)
groundnut.coords<-groundnut.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(groundnut.coords$Latitude, groundnut.coords$Longitude)))
dim(groundnut.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-groundnut.coords$Latitude
y<-groundnut.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
groundnut.data<-groundnut.final

model.non.spatial <- spatialRF::rf(
  data = groundnut.data,
  dependent.variable.name = "groundnut_HgHa",
  predictor.variable.names = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 4,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Groundnut_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Groundnut_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(groundnut.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Groundnut_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "groundnut_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "groundnut_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Groundnut_South America_RespCurv.csv")
