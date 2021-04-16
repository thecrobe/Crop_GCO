library(rgbif)
library(dplyr)
library(maps)
library(mapdata)
library(dismo) 
library(rJava) 
library(maptools)
library(jsonlite)
# Collect occurrences for pests

#Achatina.fulica
Achatina.fulica.gbf<-occ_search(scientificName = "Achatina fulica",hasCoordinate = TRUE, limit=100000)
Achatina.fulica<-data.frame(Achatina.fulica.gbf$data)
# Number of occurrences
dim(Achatina.fulica) #483 
# Coordinates
Achatina.fulica.xy<-Achatina.fulica  %>% select(scientificName,decimalLongitude, decimalLatitude)

#Claviceps.africana - TOO FEW
#Claviceps.africana.gbf<-occ_search(scientificName = "Claviceps africana",hasCoordinate = TRUE, limit=100000)
#Claviceps.africana<-data.frame(Claviceps.africana.gbf$data)
# Number of occurrences
#dim(Claviceps.africana) #13
# Coordinates
#Claviceps.africana.xy<-Claviceps.africana  %>% select(scientificName,decimalLongitude, decimalLatitude)

#Curvularia clavata - NOT PRESENT IN GBIF
#Curvularia.clavata.gbf<-occ_search(scientificName = "Claviceps clavata",hasCoordinate = TRUE)

#Eurystylus.oldi - TOO FEW
#Eurystylus.oldi.gbf<-occ_search(scientificName = "Eurystylus oldi",hasCoordinate = TRUE, limit=100000)
#Eurystylus.oldi<-data.frame(Eurystylus.oldi.gbf$data)
# Number of occurrences
#dim(Eurystylus.oldi) #4
# Coordinates
#Eurystylus.oldi.xy<-Eurystylus.oldi %>% select(scientificName,decimalLongitude, decimalLatitude)

#Maliarpha.separatella - TOO FEW
#Maliarpha.separatella.gbf<-occ_search(scientificName = " Maliarpha separatella ",hasCoordinate = TRUE, limit=100000)
#Maliarpha.separatella<-data.frame(Maliarpha.separatella.gbf$data)
# Number of occurrences
#dim(Maliarpha.separatella) #5
# Coordinates
#Maliarpha.separatella.xy<-Maliarpha.separatella %>% select(scientificName,decimalLongitude, decimalLatitude)

#Paspalum urvillei
Paspalum.urvillei.gbf<-occ_search(scientificName = "Paspalum urvillei",hasCoordinate = TRUE, limit=100000)
Paspalum.urvillei<-data.frame(Paspalum.urvillei.gbf$data)
# Number of occurrences
dim(Paspalum.urvillei) #5303
# Coordinates
Paspalum.urvillei.xy<-Paspalum.urvillei %>% select(scientificName,decimalLongitude, decimalLatitude)


#Peronosclerospora australiensis - NOT PRESENT 
#Peronosclerospora.australiensis.gbf<-occ_search(scientificName = "Peronosclerospora australiensis",hasCoordinate = TRUE)

#Peronosclerospora sargae - NOT PRESENT
#Peronosclerospora.sargae.gbf<-occ_search(scientificName = "Peronosclerospora sargae",hasCoordinate = TRUE)

#Spodoptera frugiperda
Spodoptera.frugiperda.gbf<-occ_search(scientificName = "Spodoptera frugiperda",hasCoordinate = TRUE, limit = 100000)
Spodoptera.frugiperda<-data.frame(Spodoptera.frugiperda.gbf$data)
# Number of occurrences
dim(Spodoptera.frugiperda) #2719
# Coordinates
Spodoptera.frugiperda.xy<-Spodoptera.frugiperda %>% select(scientificName,decimalLongitude, decimalLatitude)


#Verbesina encelioides
Verbesina.encelioides.gbf<-occ_search(scientificName = "Verbesina encelioides",limit = 100000)
Verbesina.encelioides<-data.frame(Verbesina.encelioides.gbf$data)
# Number of occurrences
dim(Verbesina.encelioides) #7795
# Coordinates
Verbesina.encelioides.xy<-Verbesina.encelioides %>% select(scientificName,decimalLongitude, decimalLatitude)

#Xanthium spinosum
Xanthium.spinosum.gbf<-occ_search(scientificName = "Xanthium spinosum",limit = 100000)
Xanthium.spinosum<-data.frame(Xanthium.spinosum.gbf$data)
Xanthium.spinosum
# Number of occurrences
dim(Xanthium.spinosum) #16833
# Coordinates
Xanthium.spinosum.xy<-Xanthium.spinosum %>% select(scientificName,decimalLongitude, decimalLatitude)

##### SDM MODELS #####
#Import Raster 
#generate a list of input rasters ("grids")
#pattern = "*.tif$" - filters for main raster files only and skips any associated files (e.g. world files)
grids <- list.files("wc2-5/" , pattern = "*.tif$")
#create a raster stack from the input raster files 
currentEnv <- raster::stack(paste0("wc2-5/", grids))

#Achatina.fulica SDM 
occ.Achatina.fulica<-Achatina.fulica.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)
fold <- kfold(occ.Achatina.fulica, k=5) # add an index that makes five random groups of observations
occ.Achatina.fulica.test <- occ.Achatina.fulica[fold == 1, ] # hold out one fifth as test data
occ.Achatina.fulica.train <- occ.Achatina.fulica.test[fold != 1, ] # the other four fifths are training data
maxent()
Achatina.fulica.me <- maxent(modelEnv, occ.Achatina.fulica.train ) # just using the training data
print(Achatina.fulica.me)
plot(Achatina.fulica.me) #Variable Importance 
occ.Achatina.fulica.pred <- predict(Achatina.fulica.me, modelEnv)
plot(occ.Achatina.fulica.pred)

#Paspalum.urvillei SDM 
occ.Paspalum.urvillei<-Paspalum.urvillei.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)

fold <- kfold(occ.Paspalum.urvillei, k=5) # add an index that makes five random groups of observations
occ.Paspalum.urvillei.test <- occ.Paspalum.urvillei[fold == 1, ] # hold out one fifth as test data
occ.Paspalum.urvillei.train <- occ.Paspalum.urvillei.test[fold != 1, ] # the other four fifths are training data
maxent()
Paspalum.urvillei.me <- maxent(modelEnv, occ.Paspalum.urvillei.train ) # just using the training data
print(Paspalum.urvillei.me) #model results
plot(Paspalum.urvillei.me) #Variable Importance 
Paspalum.urvillei.pred <- predict(Paspalum.urvillei.me, modelEnv)
plot(Paspalum.urvillei.pred) #check it out! 
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Paspalum.urvillei.me, p=occ.Paspalum.urvillei.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
#writeRaster(Paspalum.urvillei.pred, file = "SDMs/Sorghum_Paspalum.urvillei.tif")

#Spodoptera.frugiperda SDM 
occ.Spodoptera.frugiperda<-Spodoptera.frugiperda.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)

fold <- kfold(occ.Spodoptera.frugiperda, k=5) # add an index that makes five random groups of observations
occ.Spodoptera.frugiperda.test <- occ.Spodoptera.frugiperda[fold == 1, ] # hold out one fifth as test data
occ.Spodoptera.frugiperda.train <- occ.Spodoptera.frugiperda.test[fold != 1, ] # the other four fifths are training data
maxent()
Spodoptera.frugiperda.me <- maxent(modelEnv, occ.Spodoptera.frugiperda.train ) # just using the training data
print(Spodoptera.frugiperda.me) #model results
plot(Spodoptera.frugiperda.me) #Variable Importance 
Spodoptera.frugiperda.pred <- predict(Spodoptera.frugiperda.me, modelEnv)
plot(Spodoptera.frugiperda.pred) #check it out! 
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Spodoptera.frugiperda.me, p=occ.Spodoptera.frugiperda.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
#writeRaster(Spodoptera.frugiperda.pred, file = "SDMs/Sorghum_Spodoptera.frugiperda.tif")

#Verbesina.encelioides SDM 
occ.Verbesina.encelioides<-Verbesina.encelioides.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)

fold <- kfold(occ.Verbesina.encelioides, k=5) # add an index that makes five random groups of observations
occ.Verbesina.encelioides.test <- occ.Verbesina.encelioides[fold == 1, ] # hold out one fifth as test data
occ.Verbesina.encelioides.train <- occ.Verbesina.encelioides.test[fold != 1, ] # the other four fifths are training data
maxent()
Verbesina.encelioides.me <- maxent(modelEnv, occ.Verbesina.encelioides.train ) # just using the training data
print(Verbesina.encelioides.me) #model results
plot(Verbesina.encelioides.me) #Variable Importance 
Verbesina.encelioides.pred <- predict(Verbesina.encelioides.me, modelEnv)
plot(Verbesina.encelioides.pred) #check it out! 
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Verbesina.encelioides.me, p=occ.Verbesina.encelioides.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
#writeRaster(Verbesina.encelioides.pred, file = "SDMs/Sorghum_Verbesina.encelioides.tif")


#Xanthium.spinosum SDM 
occ.Xanthium.spinosum<-Xanthium.spinosum.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)

fold <- kfold(occ.Xanthium.spinosum, k=5) # add an index that makes five random groups of observations
occ.Xanthium.spinosum.test <- occ.Xanthium.spinosum[fold == 1, ] # hold out one fifth as test data
occ.Xanthium.spinosum.train <- occ.Xanthium.spinosum.test[fold != 1, ] # the other four fifths are training data
maxent()
Xanthium.spinosum.me <- maxent(modelEnv, occ.Xanthium.spinosum.train ) # just using the training data
print(Xanthium.spinosum.me) #model results
plot(Xanthium.spinosum.me) #Variable Importance 
Xanthium.spinosum.pred <- predict(Xanthium.spinosum.me, modelEnv)
plot(Xanthium.spinosum.pred) #check it out! 
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Xanthium.spinosum.me, p=occ.Xanthium.spinosum.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
#writeRaster(Xanthium.spinosum.pred, file = "SDMs/Sorghum_Xanthium.spinosum.tif")


#Zonal Statistics
sorghum.pest.rasters<-list(Xanthium.spinosum.pred,Verbesina.encelioides.pred,Spodoptera.frugiperda.pred,Paspalum.urvillei.pred,occ.Achatina.fulica.pred)
s <- stack(sorghum.pest.rasters)
fishnet <- shapefile("Shapefiles/Fishnet_Yield_Fertilizer_NoAnt.shp")
ex <- extract(s, fishnet, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
