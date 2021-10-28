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
rapeseedmodel<-na.omit(read.csv(file="Models/Rapeseed_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, rapeseedmodel, by='Fishnet_ID')
rapeseed<-sp::merge(m, mapping, by='Fishnet_ID')
rapeseed.p = subset(rapeseed, mean_rapes > 0) #selecting values > 0 

#Transform numerical covariates
rapeseed.num<-rapeseed.p@data %>% select(mean_rapes, AET_mean, rapeseed_Fertilizer, Pesticide,GDP_Mean) 
rapeseed.scaled<-data.frame(log10(rapeseed.num+1)) #psuedo count to avoid inf
rapeseed.scaled$FISHNET_ID<-as.numeric(rapeseed.p@data$Fishnet_ID)

#Select categorical covariates
rapeseed.cat<-rapeseed.p %>% select(COUNTRY.x,RapeseedGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-rapeseed.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
rapeseed.catagory<-sp::merge(rapeseed.cat,county.pivot, by="COUNTRY.x")
rapeseed.catagory<-rapeseed.catagory@data


########## Asia #########
#Select Continent
rapeseed.comb<-merge(rapeseed.catagory,rapeseed.scaled, by="FISHNET_ID")
rapeseed.final<-(filter((rapeseed.comb),CONTINENT == "Asia"))
rapeseed.final<-na.omit(rapeseed.final)

summary(rapeseed.final)
rapeseed.final$RapeseedGCO<-as.factor(rapeseed.final$RapeseedGCO)
rapeseed.final$COUNTRY.x<-as.factor(rapeseed.final$COUNTRY.x)
rapeseed.coords<-rapeseed.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rapeseed.coords$Latitude, rapeseed.coords$Longitude)))
dim(rapeseed.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rapeseed.coords$Latitude
y<-rapeseed.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rapeseed.data<-rapeseed.final

model.non.spatial <- spatialRF::rf(
  data = rapeseed.data,
  dependent.variable.name = "mean_rapes",
  predictor.variable.names = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"),
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


write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rapeseed_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rapeseed_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rapeseed.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rapeseed_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rapeseed_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rapeseed_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
rapeseed.comb<-merge(rapeseed.catagory,rapeseed.scaled, by="FISHNET_ID")
rapeseed.final<-(filter((rapeseed.comb),CONTINENT == "Africa"))
rapeseed.final<-na.omit(rapeseed.final)

summary(rapeseed.final)
rapeseed.final$RapeseedGCO<-as.factor(rapeseed.final$RapeseedGCO)
rapeseed.final$COUNTRY.x<-as.factor(rapeseed.final$COUNTRY.x)
rapeseed.coords<-rapeseed.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rapeseed.coords$Latitude, rapeseed.coords$Longitude)))
dim(rapeseed.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rapeseed.coords$Latitude
y<-rapeseed.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rapeseed.data<-rapeseed.final

model.non.spatial <- spatialRF::rf(
  data = rapeseed.data,
  dependent.variable.name = "mean_rapes",
  predictor.variable.names = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rapeseed_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rapeseed_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rapeseed.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rapeseed_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rapeseed_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rapeseed_Africa_RespCurv.csv")


########## North America #########
#Select Continent
rapeseed.comb<-merge(rapeseed.catagory,rapeseed.scaled, by="FISHNET_ID")
rapeseed.final<-(filter((rapeseed.comb),CONTINENT == "North America"))
rapeseed.final<-na.omit(rapeseed.final)

summary(rapeseed.final)
rapeseed.final$RapeseedGCO<-as.factor(rapeseed.final$RapeseedGCO)
rapeseed.final$COUNTRY.x<-as.factor(rapeseed.final$COUNTRY.x)
rapeseed.coords<-rapeseed.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rapeseed.coords$Latitude, rapeseed.coords$Longitude)))
dim(rapeseed.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rapeseed.coords$Latitude
y<-rapeseed.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rapeseed.data<-rapeseed.final

model.non.spatial <- spatialRF::rf(
  data = rapeseed.data,
  dependent.variable.name = "mean_rapes",
  predictor.variable.names = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rapeseed_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rapeseed_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rapeseed.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rapeseed_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "rapeseed_Fertilizer", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rapeseed_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rapeseed_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
rapeseed.comb<-merge(rapeseed.catagory,rapeseed.scaled, by="FISHNET_ID")
rapeseed.final<-(filter((rapeseed.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
rapeseed.final<-na.omit(rapeseed.final)

summary(rapeseed.final)
rapeseed.final$RapeseedGCO<-as.factor(rapeseed.final$RapeseedGCO)
rapeseed.final$COUNTRY.x<-as.factor(rapeseed.final$COUNTRY.x)
rapeseed.coords<-rapeseed.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rapeseed.coords$Latitude, rapeseed.coords$Longitude)))
dim(rapeseed.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rapeseed.coords$Latitude
y<-rapeseed.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rapeseed.data<-rapeseed.final

model.non.spatial <- spatialRF::rf(
  data = rapeseed.data,
  dependent.variable.name = "mean_rapes",
  predictor.variable.names = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rapeseed_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rapeseed_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rapeseed.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rapeseed_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rapeseed_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rapeseed_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
rapeseed.comb<-merge(rapeseed.catagory,rapeseed.scaled, by="FISHNET_ID")
rapeseed.final<-(filter((rapeseed.comb),CONTINENT == "Europe"))
rapeseed.final<-na.omit(rapeseed.final)

summary(rapeseed.final)
rapeseed.final$RapeseedGCO<-as.factor(rapeseed.final$RapeseedGCO)
rapeseed.final$COUNTRY.x<-as.factor(rapeseed.final$COUNTRY.x)
rapeseed.coords<-rapeseed.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rapeseed.coords$Latitude, rapeseed.coords$Longitude)))
dim(rapeseed.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rapeseed.coords$Latitude
y<-rapeseed.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rapeseed.data<-rapeseed.final

model.non.spatial <- spatialRF::rf(
  data = rapeseed.data,
  dependent.variable.name = "mean_rapes",
  predictor.variable.names = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rapeseed_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rapeseed_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rapeseed.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rapeseed_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rapeseed_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rapeseed_Europe_RespCurv.csv")

########## South America #########
#Select Continent
rapeseed.comb<-merge(rapeseed.catagory,rapeseed.scaled, by="FISHNET_ID")
rapeseed.final<-(filter((rapeseed.comb),CONTINENT == "South America"))
rapeseed.final<-na.omit(rapeseed.final)

summary(rapeseed.final)
rapeseed.final$RapeseedGCO<-as.factor(rapeseed.final$RapeseedGCO)
rapeseed.final$COUNTRY.x<-as.factor(rapeseed.final$COUNTRY.x)
rapeseed.coords<-rapeseed.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rapeseed.coords$Latitude, rapeseed.coords$Longitude)))
dim(rapeseed.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rapeseed.coords$Latitude
y<-rapeseed.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rapeseed.data<-rapeseed.final

model.non.spatial <- spatialRF::rf(
  data = rapeseed.data,
  dependent.variable.name = "mean_rapes",
  predictor.variable.names = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rapeseed_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rapeseed_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rapeseed.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rapeseed_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rapeseed_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rapeseed_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rapeseed_South America_RespCurv.csv")

