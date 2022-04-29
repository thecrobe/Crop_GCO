library(spatialRF) 
library(spdep)
library(sf)
library(spdplyr)
library(rgdal)
library(ggplot2)
library(ranger)
library(Rcpp)
library(parallel)
library(dplyr)

#Read In
#Fishnet- pixel = 100km^2
fishnet<- readOGR(dsn= "Global/GIS/", layer="Fishnet_yield_NoAntarctica")
summary(fishnet)
#Set Projection
proj4string(fishnet) <- CRS("+init=epsg:3786") 

#Some covariates
cassava.model<-na.omit(read.csv(file="Global/Models/Cassava_RF.csv", header=T))
mapping<-read.csv(file="Global/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, cassava.model, by='Fishnet_ID')
cassava<-sp::merge(m, mapping, by='Fishnet_ID')
cassava.p = subset(cassava, cassava_HgHa > 0) #selecting values > 0 
cassava.p$cassava10x<-cassava.p$cassava_HgHa*10 # 10 times actual amount 

#Transform numerical covariates
cassava.num<-cassava.p@data %>% select(cassava_HgHa, AET_mean, cassava_Fertilizer, Pesticide,GDP_Mean,cassava10x) 
cassava.scaled<-data.frame(log10(cassava.num+1)) #psuedo count to avoid inf
cassava.scaled$FISHNET_ID<-as.numeric(cassava.p@data$Fishnet_ID)

#Select categorical covariates
cassava.cat<-cassava.p %>% select(COUNTRY.x,CassavaGCO, FISHNET_ID,Latitude,Longitude)

#Merge

cassava.final<-merge(cassava.catagory,cassava.scaled, by="FISHNET_ID")
cassava.final$cassavaGCO<-as.factor(cassava.final$CassavaGCO)
cassava.final$cassavaGCO
cassava.inside<-filter(cassava.final, CassavaGCO == "Inside")
cassava.inside$model10x
cassava.outside<-filter(cassava.final, CassavaGCO == "Outside")
cassava10x<-data.frame(cassava.p$cassava_HgHa*10) # 10 times actual amount 
cassava.inside$x10<-cassava10x$cassava.p.cassava_HgHa...10
dim(cassava.final)
inner_join(cassava.outside,cassava10x)

cassava.final$cassavaBinaryGCO<- ifelse(cassava.final$CassavaGCO == 'Inside', 1, 0) #set binary dummy



cassava.coords<-cassava.final %>% select(Latitude,Longitude) #select coordinates
summary(cassava.final)
dim(cassava.final) #make sure dims are correct 

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

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 20,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)

#NonSpatial Model
model.non.spatial <- spatialRF::rf(
  data = cassava.data, 
  dependent.variable.name = "cassava_HgHa",
  predictor.variable.names = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  cluster = local.cluster,
  verbose = FALSE)

print(model.non.spatial)
cassava.100x<-model.non.spatial

str(cassava.data)

gc() #cleanup ram
# Moran's Plot
spatialRF::plot_moran(model.non.spatial, verbose = FALSE)

######

# Spatial Model
cassava.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed,cluster = local.cluster
)
print(cassava.spatial)
parallel::stopCluster(cl = local.cluster) #stop cluster
gc()

cas.resid<-data.frame(cassava.spatial$residuals$values) #extract residuals
dim(cas.resid)
#Plot residuals
ggplot(cassava.data, aes(x=Longitude, y=Latitude, color=cas.resid$cassava.spatial.residuals.values)) + 
  geom_point() +
  scale_color_viridis_c(option = "F") +
  theme_bw() +
  labs(color = "Eigenvalue") +
  ggtitle("Variable: Residuals") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

#Get response curves
cassava.curves.df <- spatialRF::get_response_curves(cassava.spatial,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))


write.csv(cassava.spatial$predictions, file="./PowerAnalysis/Cassava_Global10xPred.csv")
write.csv(cassava.spatial$residuals$values, file="./PowerAnalysis/Cassava_Global10xresid.csv")
write.csv(cassava.spatial$variable.importance, file="./PowerAnalysis/Cassava_Global10xVarImp.csv")
write.csv(cassava.spatial$performance, file="./PowerAnalysis/Cassava_Global10xpreformance.csv")
write.csv(cassava.curves.df, file="./PowerAnalysis/Cassava_Global10xcurves.csv")
