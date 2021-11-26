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
cassavamodel<-na.omit(read.csv(file="Models/Cassava_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, cassavamodel, by='Fishnet_ID')
cassava<-sp::merge(m, mapping, by='Fishnet_ID')
cassava.p = subset(cassava, cassava_HgHa > 0) #selecting values > 0 

#Transform numerical covariates
cassava.num<-cassava.p@data %>% select(cassava_HgHa, AET_mean, cassava_Fertilizer, Pesticide,GDP_Mean) 
cassava.scaled<-data.frame(log10(cassava.num+1)) #psuedo count to avoid inf
cassava.scaled$FISHNET_ID<-as.numeric(cassava.p@data$Fishnet_ID)

#Select categorical covariates
cassava.cat<-cassava.p %>% select(COUNTRY.x,CassavaGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-cassava.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
cassava.catagory<-sp::merge(cassava.cat,county.pivot, by="COUNTRY.x")
cassava.catagory<-cassava.catagory@data



########## Asia #########
#Select Continent
cassava.comb<-merge(cassava.catagory,cassava.scaled, by="FISHNET_ID")
cassava.final<-(filter((cassava.comb),CONTINENT == "Asia"))
cassava.final<-na.omit(cassava.final)

summary(cassava.final)
cassava.final$cassavaGCO<-as.factor(cassava.final$cassavaGCO)
cassava.final$COUNTRY.x<-as.factor(cassava.final$COUNTRY.x)
cassava.coords<-cassava.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(cassava.coords$Latitude, cassava.coords$Longitude)))
dim(cassava.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-cassava.coords$Latitude
y<-cassava.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
cassava.data<-cassava.final

model.non.spatial <- spatialRF::rf(
  data = cassava.data,
  dependent.variable.name = "cassava_HgHa",
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/cassava_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/cassava_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(cassava.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/cassava_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "cassava_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/cassava_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
cassava.comb<-merge(cassava.catagory,cassava.scaled, by="FISHNET_ID")
cassava.final<-(filter((cassava.comb),CONTINENT == "Africa"))
cassava.final<-na.omit(cassava.final)

summary(cassava.final)
cassava.final$cassavaGCO<-as.factor(cassava.final$CassavaGCO)
cassava.final$COUNTRY.x<-as.factor(cassava.final$COUNTRY.x)
cassava.coords<-cassava.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(cassava.coords$Latitude, cassava.coords$Longitude)))
dim(cassava.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-cassava.coords$Latitude
y<-cassava.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
cassava.data<-cassava.final

model.non.spatial <- spatialRF::rf(
  data = cassava.data,
  dependent.variable.name = "cassava_HgHa",
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/SpatialRF/cassava_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/cassava_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(cassava.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/cassava_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "cassava_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/cassava_Africa_RespCurv.csv")


########## North America #########
#Select Continent
cassava.comb<-merge(cassava.catagory,cassava.scaled, by="FISHNET_ID")
cassava.final<-(filter((cassava.comb),CONTINENT == "North America"))
cassava.final<-na.omit(cassava.final)

summary(cassava.final)
cassava.final$cassavaGCO<-as.factor(cassava.final$CassavaGCO)
cassava.final$COUNTRY.x<-as.factor(cassava.final$COUNTRY.x)
cassava.coords<-cassava.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(cassava.coords$Latitude, cassava.coords$Longitude)))
dim(cassava.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-cassava.coords$Latitude
y<-cassava.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
cassava.data<-cassava.final

model.non.spatial <- spatialRF::rf(
  data = cassava.data,
  dependent.variable.name = "cassava_HgHa",
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/cassava_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/cassava_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(cassava.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/cassava_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "cassava_Fertilizer", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "cassava_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/cassava_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
cassava.comb<-merge(cassava.catagory,cassava.scaled, by="FISHNET_ID")
cassava.final<-(filter((cassava.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
cassava.final<-na.omit(cassava.final)

summary(cassava.final)
cassava.final$cassavaGCO<-as.factor(cassava.final$CassavaGCO)
cassava.final$COUNTRY.x<-as.factor(cassava.final$COUNTRY.x)
cassava.coords<-cassava.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(cassava.coords$Latitude, cassava.coords$Longitude)))
dim(cassava.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-cassava.coords$Latitude
y<-cassava.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
cassava.data<-cassava.final

model.non.spatial <- spatialRF::rf(
  data = cassava.data,
  dependent.variable.name = "cassava_HgHa",
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/cassava_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/cassava_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(cassava.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/cassava_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "cassava_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/cassava_AusOcean_RespCurv.csv")

########## Europe ######### NOT PRESENT
#Select Continent
cassava.comb<-merge(cassava.catagory,cassava.scaled, by="FISHNET_ID")
cassava.final<-(filter((cassava.comb),CONTINENT == "Europe"))
cassava.final<-na.omit(cassava.final)

summary(cassava.final)
cassava.final$cassavaGCO<-as.factor(cassava.final$CassavaGCO)
cassava.final$COUNTRY.x<-as.factor(cassava.final$COUNTRY.x)
cassava.coords<-cassava.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(cassava.coords$Latitude, cassava.coords$Longitude)))
dim(cassava.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-cassava.coords$Latitude
y<-cassava.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
cassava.data<-cassava.final

model.non.spatial <- spatialRF::rf(
  data = cassava.data,
  dependent.variable.name = "cassava_HgHa",
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/cassava_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/cassava_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(cassava.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/cassava_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "cassava_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/cassava_Europe_RespCurv.csv")

########## South America #########
#Select Continent
cassava.comb<-merge(cassava.catagory,cassava.scaled, by="FISHNET_ID")
cassava.final<-(filter((cassava.comb),CONTINENT == "South America"))
cassava.final<-na.omit(cassava.final)

summary(cassava.final)
cassava.final$cassavaGCO<-as.factor(cassava.final$CassavaGCO)
cassava.final$COUNTRY.x<-as.factor(cassava.final$COUNTRY.x)
cassava.coords<-cassava.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(cassava.coords$Latitude, cassava.coords$Longitude)))
dim(cassava.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-cassava.coords$Latitude
y<-cassava.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
cassava.data<-cassava.final

model.non.spatial <- spatialRF::rf(
  data = cassava.data,
  dependent.variable.name = "cassava_HgHa",
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,n.cores = 24,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

#Spatial Model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,n.cores = 24,
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/cassava_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/cassava_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(cassava.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/cassava_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "cassava_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "cassava_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/cassava_South America_RespCurv.csv")
