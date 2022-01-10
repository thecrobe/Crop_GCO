library(spatialRF)
library(spdep)
library(spdplyr)
library(rgdal)
library(ggplot2)
library(ranger)
library(Rcpp)
library(parallel)

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


##### GLOBAL TRIAL 

#Select Continent
barley.final<-merge(barley.catagory,barley.scaled, by="FISHNET_ID")
barley.final$BarleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$BarleyGCO<- ifelse(barley.final$BarleyGCO == 'Inside', 1, 0)
barley.final$Country<-as.factor(barley.final$COUNTRY.x)
barley.final<-na.omit(barley.final)
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

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 12,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = barley.data,
  dependent.variable.name = "mean_barle",
  predictor.variable.names = c("AET_mean","Barley_Fertilizer", "Pesticide", "GDP_Mean","BarleyGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  cluster = local.cluster,
  verbose = FALSE)

# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

# Spatial Model
barley.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster
)

parallel::stopCluster(cl = local.cluster) #stop cluster
gc() #cleanup ram 

write.csv(x = model.spatial$residuals$values, file="NewAnalyses/GlobalSpatialRF/Barley_residuals.csv")
write.csv(x=model.spatial$performance, file="NewAnalyses/GlobalSpatialRF/Barley_pref.csv")




