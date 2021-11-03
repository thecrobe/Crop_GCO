
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
ricemodel<-na.omit(read.csv(file="Models/Rice_RF.csv", header=T))
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, ricemodel, by='Fishnet_ID')
rice<-sp::merge(m, mapping, by='Fishnet_ID')
rice.p = subset(rice, mean_rice_ > 0) #selecting values > 0 

#Transform numerical covariates
rice.num<-rice.p@data %>% select(mean_rice_, AET_mean, rice_Fertilizer, Pesticide,GDP_Mean) 
rice.scaled<-data.frame(log10(rice.num+1)) #psuedo count to avoid inf
rice.scaled$FISHNET_ID<-as.numeric(rice.p@data$Fishnet_ID)

#Select categorical covariates
rice.cat<-rice.p %>% select(COUNTRY.x,RiceGCO, FISHNET_ID,Latitude,Longitude, CONTINENT)

#Count number of pixels per country 
county.pivot<-rice.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
rice.catagory<-sp::merge(rice.cat,county.pivot, by="COUNTRY.x")
rice.catagory<-rice.catagory@data


########## Asia #########
#Select Continent
rice.comb<-merge(rice.catagory,rice.scaled, by="FISHNET_ID")
rice.final<-(filter((rice.comb),CONTINENT == "Asia"))
rice.final<-na.omit(rice.final)

summary(rice.final)
rice.final$RiceGCO<-as.factor(rice.final$RiceGCO)
rice.final$COUNTRY.x<-as.factor(rice.final$COUNTRY.x)
rice.coords<-rice.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rice.coords$Latitude, rice.coords$Longitude)))
dim(rice.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rice.coords$Latitude
y<-rice.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10, 50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rice.data<-rice.final

model.non.spatial <- spatialRF::rf(
  data = rice.data,
  dependent.variable.name = "mean_rice_",
  predictor.variable.names = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"),
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


write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rice_asia_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rice_asia_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rice.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rice_AsiaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"))
asia.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rice_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rice_Asia_RespCurv.csv")


########## Africa #########
#Select Continent
rice.comb<-merge(rice.catagory,rice.scaled, by="FISHNET_ID")
rice.final<-(filter((rice.comb),CONTINENT == "Africa"))
rice.final<-na.omit(rice.final)

summary(rice.final)
rice.final$RiceGCO<-as.factor(rice.final$RiceGCO)
rice.final$COUNTRY.x<-as.factor(rice.final$COUNTRY.x)
rice.coords<-rice.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rice.coords$Latitude, rice.coords$Longitude)))
dim(rice.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rice.coords$Latitude
y<-rice.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rice.data<-rice.final

model.non.spatial <- spatialRF::rf(
  data = rice.data,
  dependent.variable.name = "mean_rice_",
  predictor.variable.names = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rice_africa_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rice_africa_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rice.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rice_AfricaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"))
africa.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rice_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rice_Africa_RespCurv.csv")


########## North America #########
#Select Continent
rice.comb<-merge(rice.catagory,rice.scaled, by="FISHNET_ID")
rice.final<-(filter((rice.comb),CONTINENT == "North America"))
rice.final<-na.omit(rice.final)

summary(rice.final)
rice.final$RiceGCO<-as.factor(rice.final$RiceGCO)
rice.final$COUNTRY.x<-as.factor(rice.final$COUNTRY.x)
rice.coords<-rice.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rice.coords$Latitude, rice.coords$Longitude)))
dim(rice.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rice.coords$Latitude
y<-rice.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rice.data<-rice.final

model.non.spatial <- spatialRF::rf(
  data = rice.data,
  dependent.variable.name = "mean_rice_",
  predictor.variable.names = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rice_NorthAmerica_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rice_NorthAmerica_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rice.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rice_NorthAmerica_VarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables =c("Pesticide", "AET_mean", "rice_Fertilizer", "GDP_Mean"))
NorthAm.rc<-spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rice_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rice_NorthAmerica_RespCurv.csv")


########## Australia + Oceania #########
#Select Continent
rice.comb<-merge(rice.catagory,rice.scaled, by="FISHNET_ID")
rice.final<-(filter((rice.comb),CONTINENT == "Australia" | CONTINENT=="Oceania"))
rice.final<-na.omit(rice.final)

summary(rice.final)
rice.final$RiceGCO<-as.factor(rice.final$RiceGCO)
rice.final$COUNTRY.x<-as.factor(rice.final$COUNTRY.x)
rice.coords<-rice.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rice.coords$Latitude, rice.coords$Longitude)))
dim(rice.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rice.coords$Latitude
y<-rice.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rice.data<-rice.final

model.non.spatial <- spatialRF::rf(
  data = rice.data,
  dependent.variable.name = "mean_rice_",
  predictor.variable.names = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rice_ausocean_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rice_ausocean_pref.csv")

#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rice.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rice_AusOceanVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rice_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rice_AusOcean_RespCurv.csv")

########## Europe #########
#Select Continent
rice.comb<-merge(rice.catagory,rice.scaled, by="FISHNET_ID")
rice.final<-(filter((rice.comb),CONTINENT == "Europe"))
rice.final<-na.omit(rice.final)

summary(rice.final)
rice.final$RiceGCO<-as.factor(rice.final$RiceGCO)
rice.final$COUNTRY.x<-as.factor(rice.final$COUNTRY.x)
rice.coords<-rice.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rice.coords$Latitude, rice.coords$Longitude)))
dim(rice.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rice.coords$Latitude
y<-rice.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rice.data<-rice.final

model.non.spatial <- spatialRF::rf(
  data = rice.data,
  dependent.variable.name = "mean_rice_",
  predictor.variable.names = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rice_europe_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rice_europe_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rice.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rice_EuropeVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rice_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rice_Europe_RespCurv.csv")

########## South America #########
#Select Continent
rice.comb<-merge(rice.catagory,rice.scaled, by="FISHNET_ID")
rice.final<-(filter((rice.comb),CONTINENT == "South America"))
rice.final<-na.omit(rice.final)

summary(rice.final)
rice.final$RiceGCO<-as.factor(rice.final$RiceGCO)
rice.final$COUNTRY.x<-as.factor(rice.final$COUNTRY.x)
rice.coords<-rice.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(rice.coords$Latitude, rice.coords$Longitude)))
dim(rice.final) #check dims
dim(dist.mat) #check dims

#coordinates of the cases
x<-rice.coords$Latitude
y<-rice.coords$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 5, 10,50)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
rice.data<-rice.final

model.non.spatial <- spatialRF::rf(
  data = rice.data,
  dependent.variable.name = "mean_rice_",
  predictor.variable.names = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"),
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

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/SpatialRF/Rice_South America_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/SpatialRF/Rice_South America_pref.csv")
#Check residuals
spatialRF::plot_residuals_diagnostics(model.spatial,verbose = FALSE)

ggplot(rice.final, aes(x=Longitude, y=Latitude, color=model.spatial$residuals$values)) + 
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
write.csv(x=model.spatial$variable.importance, file="NewAnalyses/SpatialRF/Rice_South AmericaVarImp.csv")

#Response Curves
reponse.curves.df <- spatialRF::get_response_curves(model.spatial,variables = c("AET_mean", "rice_Fertilizer", "Pesticide", "GDP_Mean"))
spatialRF::plot_response_curves(model.spatial, quantiles = 0.5,ncol = 2,variables = c("Pesticide", "AET_mean", "rice_Fertilizer", "GDP_Mean"))
write.csv(x=reponse.curves.df,file="NewAnalyses/SpatialRF/Rice_South America_RespCurv.csv")

