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
Potato_model<-na.omit(read.csv(file="Models/Potato_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, Potato_model, by='Fishnet_ID')
Potato<-sp::merge(m, mapping, by='Fishnet_ID')
Potato.p = subset(Potato, mean_potat.x > 0) #selecting values > 0 

#Transform numerical covariates
Potato.num<-Potato.p@data %>% select(mean_potat.x, AET_mean, Potato_Fe, Pesticide,GDP_Mean) 
Potato.scaled<-data.frame(log10(Potato.num+1)) #psuedo count to avoid inf
Potato.scaled$FISHNET_ID<-as.numeric(Potato.p@data$Fishnet_ID)

#Select categorical covariates
Potato.cat<-Potato.p %>% select(COUNTRY.x,PotatoGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-Potato.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
Potato.catagory<-sp::merge(Potato.cat,county.pivot, by="COUNTRY.x")
Potato.catagory<-Potato.catagory@data


########## Asia #########
#Select Continent
Potato.comb<-merge(Potato.catagory,Potato.scaled, by="FISHNET_ID")
Potato.final<-(filter((Potato.comb),CONTINENT == "Asia"))
Potato.final<-na.omit(Potato.final)

summary(Potato.final)
Potato.final$PotatoGCO<-as.factor(Potato.final$PotatoGCO)
Potato.final$COUNTRY.x<-as.factor(Potato.final$COUNTRY.x)
Potato.coords<-Potato.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(Potato.coords$Latitude, Potato.coords$Longitude)))
dim(Potato.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-Potato.coords$Latitude
y<-Potato.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
Potato.data<-Potato.final

model.non.spatial <- spatialRF::rf(
  data = Potato.data,
  dependent.variable.name = "mean_potat.x",
  predictor.variable.names = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"),
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


write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Potato_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Potato_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(Potato.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Potato_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Potato_Fe", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Potato_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
Potato.comb<-merge(Potato.catagory,Potato.scaled, by="FISHNET_ID")
Potato.final<-(filter((Potato.comb),CONTINENT == "Africa"))
Potato.final<-na.omit(Potato.final)

summary(Potato.final)
Potato.final$PotatoGCO<-as.factor(Potato.final$PotatoGCO)
Potato.final$COUNTRY.x<-as.factor(Potato.final$COUNTRY.x)
Potato.coords<-Potato.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(Potato.coords$Latitude, Potato.coords$Longitude)))
dim(Potato.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-Potato.coords$Latitude
y<-Potato.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
Potato.data<-Potato.final

model.non.spatial <- spatialRF::rf(
  data = Potato.data,
  dependent.variable.name = "mean_potat.x",
  predictor.variable.names = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Potato_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Potato_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(Potato.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Potato_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Potato_Fe", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Potato_Africa_RespCurv.csv")

########## North America #########
#Select Continent
Potato.comb<-merge(Potato.catagory,Potato.scaled, by="FISHNET_ID")
Potato.final<-(filter((Potato.comb),CONTINENT == "North America"))
Potato.final<-na.omit(Potato.final)

summary(Potato.final)
Potato.final$PotatoGCO<-as.factor(Potato.final$PotatoGCO)
Potato.final$COUNTRY.x<-as.factor(Potato.final$COUNTRY.x)
Potato.coords<-Potato.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(Potato.coords$Latitude, Potato.coords$Longitude)))
dim(Potato.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-Potato.coords$Latitude
y<-Potato.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
Potato.data<-Potato.final

model.non.spatial <- spatialRF::rf(
  data = Potato.data,
  dependent.variable.name = "mean_potat.x",
  predictor.variable.names = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Potato_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Potato_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(Potato.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Potato_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "Potato_Fe", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Potato_Fe", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Potato_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
Potato.comb<-merge(Potato.catagory,Potato.scaled, by="FISHNET_ID")
Potato.final<-(filter((Potato.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
Potato.final<-na.omit(Potato.final)

summary(Potato.final)
Potato.final$PotatoGCO<-as.factor(Potato.final$PotatoGCO)
Potato.final$COUNTRY.x<-as.factor(Potato.final$COUNTRY.x)
Potato.coords<-Potato.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(Potato.coords$Latitude, Potato.coords$Longitude)))
dim(Potato.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-Potato.coords$Latitude
y<-Potato.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
Potato.data<-Potato.final

model.non.spatial <- spatialRF::rf(
  data = Potato.data,
  dependent.variable.name = "mean_potat.x",
  predictor.variable.names = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Potato_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Potato_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(Potato.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Potato_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Potato_Fe", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Potato_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
Potato.comb<-merge(Potato.catagory,Potato.scaled, by="FISHNET_ID")
Potato.final<-(filter((Potato.comb),CONTINENT == "Europe"))
Potato.final<-na.omit(Potato.final)

summary(Potato.final)
Potato.final$PotatoGCO<-as.factor(Potato.final$PotatoGCO)
Potato.final$COUNTRY.x<-as.factor(Potato.final$COUNTRY.x)
Potato.coords<-Potato.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(Potato.coords$Latitude, Potato.coords$Longitude)))
dim(Potato.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-Potato.coords$Latitude
y<-Potato.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
Potato.data<-Potato.final

model.non.spatial <- spatialRF::rf(
  data = Potato.data,
  dependent.variable.name = "mean_potat.x",
  predictor.variable.names = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Potato_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Potato_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(Potato.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Potato_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Potato_Fe", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Potato_Europe_RespCurv.csv")

########## South America #########
#Select Continent
Potato.comb<-merge(Potato.catagory,Potato.scaled, by="FISHNET_ID")
Potato.final<-(filter((Potato.comb),CONTINENT == "South America"))
Potato.final<-na.omit(Potato.final)

summary(Potato.final)
Potato.final$PotatoGCO<-as.factor(Potato.final$PotatoGCO)
Potato.final$COUNTRY.x<-as.factor(Potato.final$COUNTRY.x)
Potato.coords<-Potato.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(Potato.coords$Latitude, Potato.coords$Longitude)))
dim(Potato.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-Potato.coords$Latitude
y<-Potato.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
Potato.data<-Potato.final

model.non.spatial <- spatialRF::rf(
  data = Potato.data,
  dependent.variable.name = "mean_potat.x",
  predictor.variable.names = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Potato_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Potato_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(Potato.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Potato_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Potato_Fe", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Potato_Fe", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Potato_South America_RespCurv.csv")

