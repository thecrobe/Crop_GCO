

library(rgbif)
library(dplyr)
library(maps)
library(mapdata)
library(dismo) 
library(rJava) 
library(maptools)
library(jsonlite)

#Import Raster 
#generate a list of input rasters ("grids")
#pattern = "*.tif$" - filters for main raster files only and skips any associated files (e.g. world files)
grids <- list.files("wc2-5/" , pattern = "*.tif$")
#create a raster stack from the input raster files 
currentEnv <- raster::stack(paste0("wc2-5/", grids))


#Achatina.fulica
Achatina.fulica.gbf<-occ_search(scientificName = "Achatina fulica",hasCoordinate = TRUE, limit=100000)
Achatina.fulica<-data.frame(Achatina.fulica.gbf$data)
# Number of occurrences
dim(Achatina.fulica) #483 
# Coordinates
Achatina.fulica.xy<-Achatina.fulica  %>% select(scientificName,decimalLongitude, decimalLatitude)

#Achatina.fulica SDM 
occ.Achatina.fulica<-Achatina.fulica.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)

fold <- kfold(occ.Achatina.fulica, k=5) # add an index that makes five random groups of observations
occ.Achatina.fulica.test <- occ.Achatina.fulica[fold == 1, ] # hold out one fifth as test data
occ.Achatina.fulica.train <- occ.Achatina.fulica.test[fold != 1, ] # the other four fifths are training data
maxent()
Achatina.fulica.me <- maxent(modelEnv, occ.Achatina.fulica.train ) # just using the training data
print(Achatina.fulica.me) #model results
plot(Achatina.fulica.me) #Variable Importance 
Achatina.fulica.pred <- predict(Achatina.fulica.me, modelEnv)
plot(Achatina.fulica.pred) #check it out! 
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Achatina.fulica.me, p=occ.Achatina.fulica.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
#writeRaster(Achatina.fulica.pred, file = "SDMs/Sorghum_Achatina.fulica.tif")


#Claviceps.africana - Too Few
#Claviceps.africana.gbf<-occ_search(scientificName = "Claviceps africana",hasCoordinate = TRUE, limit=100000)
#Claviceps.africana<-data.frame(Claviceps.africana.gbf$data)
# Number of occurrences
#dim(Claviceps.africana) #13
# Coordinates
#Claviceps.africana.xy<-Claviceps.africana  %>% select(scientificName,decimalLongitude, decimalLatitude)

#Curvularia clavata - NOT PRESENT IN GBIF
#Curvularia.clavata.gbf<-occ_search(scientificName = "Claviceps clavata",hasCoordinate = TRUE)

#Eurystylus.oldi - Too Few
#Eurystylus.oldi.gbf<-occ_search(scientificName = "Eurystylus oldi",hasCoordinate = TRUE, limit=100000)
#Eurystylus.oldi<-data.frame(Eurystylus.oldi.gbf$data)
# Number of occurrences
#dim(Eurystylus.oldi) #4
# Coordinates
#Eurystylus.oldi.xy<-Eurystylus.oldi %>% select(scientificName,decimalLongitude, decimalLatitude)

#Maliarpha.separatella - Too Few
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
writeRaster(Paspalum.urvillei.pred, file = "SDMs/Sorghum_Paspalum.urvillei.tif")


#Peronosclerospora australiensis - NOT PRESENT 
#Peronosclerospora.australiensis.gbf<-occ_search(scientificName = "Peronosclerospora australiensis",hasCoordinate = TRUE)

#Peronosclerospora sargae - Not Present
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
dim(Xanthium.spinosum) #16832
# Coordinates
Xanthium.spinosum.xy<-Xanthium.spinosum %>% select(scientificName,decimalLongitude, decimalLatitude)

