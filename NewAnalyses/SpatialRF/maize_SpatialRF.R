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
maizemodel<-na.omit(read.csv(file="Models/Maize_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, maizemodel, by='Fishnet_ID')
maize<-sp::merge(m, mapping, by='Fishnet_ID')
maize.p = subset(maize, mean_maize > 0) #selecting values > 0 

#Transform numerical covariates
maize.num<-maize.p@data %>% select(mean_maize, AET_mean, maize_Fertilizer, Pesticide,GDP_Mean) 
maize.scaled<-data.frame(log10(maize.num+1)) #psuedo count to avoid inf
maize.scaled$FISHNET_ID<-as.numeric(maize.p@data$Fishnet_ID)

#Select categorical covariates
maize.cat<-maize.p %>% select(COUNTRY.x,MaizeGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-maize.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
maize.catagory<-sp::merge(maize.cat,county.pivot, by="COUNTRY.x")
maize.catagory<-maize.catagory@data



########## Asia #########
#Select Continent
maize.comb<-merge(maize.catagory,maize.scaled, by="FISHNET_ID")
maize.final<-(filter((maize.comb),CONTINENT == "Asia"))
maize.final<-na.omit(maize.final)

summary(maize.final)
maize.final$MaizeGCO<-as.factor(maize.final$MaizeGCO)
maize.final$COUNTRY.x<-as.factor(maize.final$COUNTRY.x)
maize.coords<-maize.final %>% select(Latitude,Longitude) #select coordinates

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
maize.data<-maize.final

model.non.spatial <- spatialRF::rf(
  data = maize.data,
  dependent.variable.name = "mean_maize",
  predictor.variable.names = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Maize_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Maize_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(maize.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Maize_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "maize_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Maize_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
maize.comb<-merge(maize.catagory,maize.scaled, by="FISHNET_ID")
maize.final<-(filter((maize.comb),CONTINENT == "Africa"))
maize.final<-na.omit(maize.final)

summary(maize.final)
maize.final$MaizeGCO<-as.factor(maize.final$MaizeGCO)
maize.final$COUNTRY.x<-as.factor(maize.final$COUNTRY.x)
maize.coords<-maize.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
maize.data<-maize.final

model.non.spatial <- spatialRF::rf(
  data = maize.data,
  dependent.variable.name = "mean_maize",
  predictor.variable.names = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/SpatialRF/Maize_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Maize_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(maize.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Maize_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "maize_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Maize_Africa_RespCurv.csv")


########## North America #########
#Select Continent
maize.comb<-merge(maize.catagory,maize.scaled, by="FISHNET_ID")
maize.final<-(filter((maize.comb),CONTINENT == "North America"))
maize.final<-na.omit(maize.final)

summary(maize.final)
maize.final$MaizeGCO<-as.factor(maize.final$MaizeGCO)
maize.final$COUNTRY.x<-as.factor(maize.final$COUNTRY.x)
maize.coords<-maize.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
maize.data<-maize.final

model.non.spatial <- spatialRF::rf(
  data = maize.data,
  dependent.variable.name = "mean_maize",
  predictor.variable.names = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Maize_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Maize_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(maize.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Maize_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "maize_Fertilizer", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "maize_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Maize_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
maize.comb<-merge(maize.catagory,maize.scaled, by="FISHNET_ID")
maize.final<-(filter((maize.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
maize.final<-na.omit(maize.final)

summary(maize.final)
maize.final$MaizeGCO<-as.factor(maize.final$MaizeGCO)
maize.final$COUNTRY.x<-as.factor(maize.final$COUNTRY.x)
maize.coords<-maize.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
maize.data<-maize.final

model.non.spatial <- spatialRF::rf(
  data = maize.data,
  dependent.variable.name = "mean_maize",
  predictor.variable.names = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Maize_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Maize_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(maize.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Maize_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "maize_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Maize_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
maize.comb<-merge(maize.catagory,maize.scaled, by="FISHNET_ID")
maize.final<-(filter((maize.comb),CONTINENT == "Europe"))
maize.final<-na.omit(maize.final)

summary(maize.final)
maize.final$MaizeGCO<-as.factor(maize.final$MaizeGCO)
maize.final$COUNTRY.x<-as.factor(maize.final$COUNTRY.x)
maize.coords<-maize.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
maize.data<-maize.final

model.non.spatial <- spatialRF::rf(
  data = maize.data,
  dependent.variable.name = "mean_maize",
  predictor.variable.names = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Maize_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Maize_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(maize.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Maize_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "maize_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Maize_Europe_RespCurv.csv")

########## South America #########
#Select Continent
maize.comb<-merge(maize.catagory,maize.scaled, by="FISHNET_ID")
maize.final<-(filter((maize.comb),CONTINENT == "South America"))
maize.final<-na.omit(maize.final)

summary(maize.final)
maize.final$MaizeGCO<-as.factor(maize.final$MaizeGCO)
maize.final$COUNTRY.x<-as.factor(maize.final$COUNTRY.x)
maize.coords<-maize.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
maize.data<-maize.final

model.non.spatial <- spatialRF::rf(
  data = maize.data,
  dependent.variable.name = "mean_maize",
  predictor.variable.names = c("AET_mean", "maize_Fertilizer", "Pesticide", "GDP_Mean"),
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
  seed = random.seed,
  n.cores = 24
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Maize_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Maize_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(maize.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Maize_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "maize_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Maize_South America_RespCurv.csv")

