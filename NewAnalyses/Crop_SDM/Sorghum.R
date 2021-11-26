library(dplyr)
library(maps)
library(mapdata)
library(dismo) 
library(rJava) 
library(maptools)
library(jsonlite)
library(ggplot2)


fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
fishnet.sorghum<-subset(fishnet, fishnet$mean_sorgh > 0 ) #yield > 0 
sorghum.xy<-data.frame(cbind(fishnet.sorghum$Longitude,fishnet.sorghum$Latitude))

grids <- list.files("wc2-5/" , pattern = "*.tif$")
#create a raster stack from the input raster files 
currentEnv <- raster::stack(paste0("wc2-5/", grids))



#Sorghum SDM 
occ.Sorghum<-sorghum.xy
model.extent<-extent(-180,180,-90,90)
modelEnv<-crop(currentEnv,model.extent)
fold <- kfold(occ.Sorghum, k=5) # add an index that makes five random groups of observations
occ.Sorghum.test <- occ.Sorghum[fold == 1, ] # hold out one fifth as test data
occ.Sorghum.train <- occ.Sorghum.test[fold != 1, ] # the other four fifths are training data
maxent()
Sorghum.me <- maxent(modelEnv, occ.Sorghum.train ) # just using the training data
print(Sorghum.me)
plot(Sorghum.me) #Variable Importance 
occ.Sorghum.pred <- predict(Sorghum.me, model
plot(occ.Sorghum.pred)
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Sorghum.me, p=occ.Sorghum.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
#writeRaster(Sorghum.pred, file = "SDMs/Sorghum.urvillei.tif")

