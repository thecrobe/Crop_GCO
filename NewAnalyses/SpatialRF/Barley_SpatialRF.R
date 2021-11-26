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
barleymodel<-na.omit(read.csv(file="Models/Barley_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, barleymodel, by='Fishnet_ID')
barley<-sp::merge(m, mapping, by='Fishnet_ID')
barley.p = subset(barley, mean_barle > 0) #selecting values > 0 

#Transform numerical covariates
barley.num<-barley.p@data %>% select(mean_barle, AET_mean, Barley_Fertilizer, Pesticide,GDP_Mean) 
barley.scaled<-data.frame(log10(barley.num+1)) #psuedo count to avoid inf
barley.scaled$FISHNET_ID<-as.numeric(barley.p@data$Fishnet_ID)

#Select categorical covariates
barley.cat<-barley.p %>% select(COUNTRY.x,BarleyGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-barley.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
barley.catagory<-sp::merge(barley.cat,county.pivot, by="COUNTRY.x")
barley.catagory<-barley.catagory@data



########## Asia #########
#Select Continent
barley.comb<-merge(barley.catagory,barley.scaled, by="FISHNET_ID")
barley.final<-(filter((barley.comb),CONTINENT == "Asia"))
barley.final<-na.omit(barley.final)

summary(barley.final)
barley.final$BarleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$COUNTRY.x<-as.factor(barley.final$COUNTRY.x)
barley.coords<-barley.final %>% select(Latitude,Longitude) #select coordinates

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
barley.data<-barley.final

model.non.spatial <- spatialRF::rf(
  data = barley.data,
  dependent.variable.name = "mean_barle",
  predictor.variable.names = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Barley_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Barley_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(barley.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Barley_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Barley_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Barley_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
barley.comb<-merge(barley.catagory,barley.scaled, by="FISHNET_ID")
barley.final<-(filter((barley.comb),CONTINENT == "Africa"))
barley.final<-na.omit(barley.final)

summary(barley.final)
barley.final$BarleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$COUNTRY.x<-as.factor(barley.final$COUNTRY.x)
barley.coords<-barley.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
barley.data<-barley.final

model.non.spatial <- spatialRF::rf(
  data = barley.data,
  dependent.variable.name = "mean_barle",
  predictor.variable.names = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/SpatialRF/Barley_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Barley_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(barley.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Barley_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Barley_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Barley_Africa_RespCurv.csv")


########## North America #########
#Select Continent
barley.comb<-merge(barley.catagory,barley.scaled, by="FISHNET_ID")
barley.final<-(filter((barley.comb),CONTINENT == "North America"))
barley.final<-na.omit(barley.final)

summary(barley.final)
barley.final$BarleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$COUNTRY.x<-as.factor(barley.final$COUNTRY.x)
barley.coords<-barley.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
barley.data<-barley.final

model.non.spatial <- spatialRF::rf(
  data = barley.data,
  dependent.variable.name = "mean_barle",
  predictor.variable.names = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Barley_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Barley_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(barley.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Barley_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "Barley_Fertilizer", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Barley_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Barley_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
barley.comb<-merge(barley.catagory,barley.scaled, by="FISHNET_ID")
barley.final<-(filter((barley.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
barley.final<-na.omit(barley.final)

summary(barley.final)
barley.final$BarleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$COUNTRY.x<-as.factor(barley.final$COUNTRY.x)
barley.coords<-barley.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
barley.data<-barley.final

model.non.spatial <- spatialRF::rf(
  data = barley.data,
  dependent.variable.name = "mean_barle",
  predictor.variable.names = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Barley_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Barley_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(barley.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Barley_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Barley_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Barley_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
barley.comb<-merge(barley.catagory,barley.scaled, by="FISHNET_ID")
barley.final<-(filter((barley.comb),CONTINENT == "Europe"))
barley.final<-na.omit(barley.final)

summary(barley.final)
barley.final$BarleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$COUNTRY.x<-as.factor(barley.final$COUNTRY.x)
barley.coords<-barley.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
barley.data<-barley.final

model.non.spatial <- spatialRF::rf(
  data = barley.data,
  dependent.variable.name = "mean_barle",
  predictor.variable.names = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Barley_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Barley_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(barley.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Barley_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Barley_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Barley_Europe_RespCurv.csv")

########## South America #########
#Select Continent
barley.comb<-merge(barley.catagory,barley.scaled, by="FISHNET_ID")
barley.final<-(filter((barley.comb),CONTINENT == "South America"))
barley.final<-na.omit(barley.final)

summary(barley.final)
barley.final$BarleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$COUNTRY.x<-as.factor(barley.final$COUNTRY.x)
barley.coords<-barley.final %>% select(Latitude,Longitude) #select coordinates

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
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
barley.data<-barley.final

model.non.spatial <- spatialRF::rf(
  data = barley.data,
  dependent.variable.name = "mean_barle",
  predictor.variable.names = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"),
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
  n.cores = 24,
  seed = random.seed
)

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Barley_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Barley_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(barley.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Barley_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "Barley_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Barley_South America_RespCurv.csv")
