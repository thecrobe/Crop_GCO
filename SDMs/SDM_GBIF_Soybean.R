library(rgbif)
library(dplyr)
library(maps)
library(mapdata)
library(dismo) 
library(rJava) 
library(maptools)
library(jsonlite)
library(ggplot2)
library(RVAideMemoire)
# Collect occurrences for pests
# https://plantvillage.psu.edu/topics/soybean/infos/diseases_and_pests_description_uses_propagation

#Pseudomonas syringae
Pseudomonas.syringae.gbf<-occ_search(scientificName = "Pseudomonas syringae",hasCoordinate = TRUE, limit=100000)
Pseudomonas.syringae<-data.frame(Pseudomonas.syringae.gbf$data)
# Number of occurrences
dim(Pseudomonas.syringae) #1223
# Coordinates
Pseudomonas.syringae.xy<-Pseudomonas.syringae %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

#Xanthomonas campestris
Xanthomonas.campestris.gbf<-occ_search(scientificName = "Xanthomonas campestris",hasCoordinate = TRUE, limit=100000)
Xanthomonas.campestris<-data.frame(Xanthomonas.campestris.gbf$data)
# Number of occurrences
dim(Xanthomonas.campestris) #577
# Coordinates
Xanthomonas.campestris.xy<-Xanthomonas.campestris %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Rhizoctonia solani
Rhizoctonia.solani.gbf<-occ_search(scientificName = "Rhizoctonia solani",hasCoordinate = TRUE, limit=100000)
Rhizoctonia.solani<-data.frame(Rhizoctonia.solani.gbf$data)
# Number of occurrences
dim(Rhizoctonia.solani) #1270
# Coordinates
Rhizoctonia.solani.xy<-Rhizoctonia.solani %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Phakopsora pachyrhizi - TOO FEW
#Phakopsora.pachyrhizi.gbf<-occ_search(scientificName = "Phakopsora.pachyrhizi ",hasCoordinate = TRUE, limit=100000)
#Phakopsora.pachyrhizi<-data.frame(Phakopsora.pachyrhizi.gbf$data)
# Number of occurrences
#dim(Phakopsora.pachyrhizi) #268
# Coordinates
#Phakopsora.pachyrhizi.xy<-Phakopsora.pachyrhizi %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Sclerotinia clerotiorum - Not In Database
#Sclerotinia.clerotiorum.gbf<-occ_search(scientificName = "Sclerotinia clerotiorum",hasCoordinate = TRUE, limit=100000)
#Sclerotinia.clerotiorum<-data.frame(Sclerotinia.clerotiorum.gbf$data)
# Number of occurrences
#dim(Sclerotinia.clerotiorum) #483 
# Coordinates
#Sclerotinia.clerotiorum.xy<-Sclerotinia.clerotiorum %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Septoria glycines - TOO FEW
#Septoria.glycines.gbf<-occ_search(scientificName = "Septoria glycines",hasCoordinate = TRUE, limit=100000)
#Septoria.glycines<-data.frame(Septoria.glycines.gbf$data)
# Number of occurrences
#dim(Septoria.glycines) #5
# Coordinates
#Septoria.glycines.xy<-Septoria.glycines %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Spodoptera exigua
Spodoptera.exigua.gbf<-occ_search(scientificName = "Spodoptera exigua",hasCoordinate = TRUE, limit=100000)
Spodoptera.exigua<-data.frame(Spodoptera.exigua.gbf$data)
# Number of occurrences
dim(Spodoptera.exigua) #8338
# Coordinates
Spodoptera.exigua.xy<-Spodoptera.exigua %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Spodoptera praefica - TOO FEW
#Spodoptera.praefica.gbf<-occ_search(scientificName = "Spodoptera praefica",hasCoordinate = TRUE, limit=100000)
#Spodoptera.praefica<-data.frame(Spodoptera.praefica.gbf$data)
# Number of occurrences
#dim(Spodoptera.praefica) #4379
# Coordinates
#Spodoptera.praefica.xy<-Spodoptera.praefica %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Acalymma vittata - NOT IN DATABASE
#Acalymma.vittata.gbf<-occ_search(scientificName = "Acalymma vittata",hasCoordinate = TRUE, limit=100000)
#Acalymma.vittata<-data.frame(Acalymma.vittata.gbf$data)
# Number of occurrences
#dim(Acalymma.vittata) #4379
# Coordinates
#Acalymma.vittata.xy<-Acalymma.vittata %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Diabrotica undecimpunctata
Diabrotica.undecimpunctata.gbf<-occ_search(scientificName = "Diabrotica undecimpunctata",hasCoordinate = TRUE, limit=100000)
Diabrotica.undecimpunctata<-data.frame(Diabrotica.undecimpunctata.gbf$data)
# Number of occurrences
dim(Diabrotica.undecimpunctata) #11919
# Coordinates
Diabrotica.undecimpunctata.xy<-Diabrotica.undecimpunctata %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Epilachna svarivestis - TOO FEW
#Epilachna.varivestis.gbf<-occ_search(scientificName = "Epilachna varivestis",hasCoordinate = TRUE, limit=100000)
#Epilachna.varivestis<-data.frame(Epilachna.varivestis.gbf$data)
# Number of occurrences
#dim(Epilachna.varivestis) #459
# Coordinates
#Epilachna.varivestis.xy<-Epilachna.varivestis %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Septoria glycines - TOO FEW
#Septoria.glycines.gbf<-occ_search(scientificName = "Septoria glycines",hasCoordinate = TRUE, limit=100000)
#Septoria.glycines<-data.frame(Septoria.glycines.gbf$data)
# Number of occurrences
#dim(Septoria.glycines) #4379
# Coordinates
#Septoria.glycines.xy<-Septoria.glycines %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Cercospora kikuchii - NOT IN DATABASE
#Cercospora.kikuchii.gbf<-occ_search(scientificName = "Cercospora kikuchii",hasCoordinate = TRUE, limit=100000)
#Cercospora.kikuchii<-data.frame(Cercospora.kikuchii.gbf$data)
# Number of occurrences
#dim(Cercospora.kikuchii) #4379
# Coordinates
#Cercospora.kikuchii.xy<-Cercospora.kikuchii %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Macrophomina phaseolina - TOO FEW
#Macrophomina.phaseolina.gbf<-occ_search(scientificName = "Macrophomina phaseolina",hasCoordinate = TRUE, limit=100000)
#Macrophomina.phaseolina<-data.frame(Macrophomina.phaseolina.gbf$data)
# Number of occurrences
#dim(Macrophomina.phaseolina) #235 
# Coordinates
#Macrophomina.phaseolina.xy<-Macrophomina.phaseolina %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Cercospora sojina - TOO FEW
#Cercospora.sojina.gbf<-occ_search(scientificName = "Cercospora sojina",hasCoordinate = TRUE, limit=100000)
#Cercospora.sojina<-data.frame(Cercospora.sojina.gbf$data)
# Number of occurrences
#dim(Cercospora.sojina) #10
# Coordinates
#Cercospora.sojina.xy<-Cercospora.sojina %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)

# Phytophthora megasperma - TOO FEW
#Phytophthora.megasperma.gbf<-occ_search(scientificName = "Phytophthora megasperma",hasCoordinate = TRUE, limit=100000)
#Phytophthora.megasperma<-data.frame(Phytophthora.megasperma.gbf$data)
# Number of occurrences
#dim(Phytophthora.megasperma) #50
# Coordinates
#Phytophthora.megasperma.xy<-Phytophthora.megasperma %>% dplyr::select(scientificName,decimalLongitude, decimalLatitude)


##### SDM MODELS #####
#Import Raster 
#generate a list of input rasters ("grids")
#pattern = "*.tif$" - filters for main raster files only and skips any associated files (e.g. world files)
grids <- list.files("wc2-5/" , pattern = "*.tif$")
#create a raster stack from the input raster files 
currentEnv <- raster::stack(paste0("wc2-5/", grids))

#Pseudomonas.syringae SDM 
occ.Pseudomonas.syringae<-Pseudomonas.syringae.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)
fold <- kfold(occ.Pseudomonas.syringae, k=5) # add an index that makes five random groups of observations
occ.Pseudomonas.syringae.test <- occ.Pseudomonas.syringae[fold == 1, ] # hold out one fifth as test data
occ.Pseudomonas.syringae.train <- occ.Pseudomonas.syringae.test[fold != 1, ] # the other four fifths are training data
maxent()
Pseudomonas.syringae.me <- maxent(modelEnv, occ.Pseudomonas.syringae.train ) # just using the training data
print(Pseudomonas.syringae.me)
plot(Pseudomonas.syringae.me) #Variable Importance 
occ.Pseudomonas.syringae.pred <- predict(Pseudomonas.syringae.me, modelEnv)
plot(occ.Pseudomonas.syringae.pred)
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Pseudomonas.syringae.me, p=occ.Pseudomonas.syringae.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
Pseudomonas.syringae.pred <- predict(Pseudomonas.syringae.me, modelEnv)
#writeRaster(Pseudomonas.syringae.pred, file = "SDMs/Soybean_Pseudomonas.syringae.tif")

#Xanthomonas.campestris SDM 
occ.Xanthomonas.campestris<-Xanthomonas.campestris.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)
fold <- kfold(occ.Xanthomonas.campestris, k=5) # add an index that makes five random groups of observations
occ.Xanthomonas.campestris.test <- occ.Xanthomonas.campestris[fold == 1, ] # hold out one fifth as test data
occ.Xanthomonas.campestris.train <- occ.Xanthomonas.campestris.test[fold != 1, ] # the other four fifths are training data
maxent()
Xanthomonas.campestris.me <- maxent(modelEnv, occ.Xanthomonas.campestris.train ) # just using the training data
print(Xanthomonas.campestris.me)
plot(Xanthomonas.campestris.me) #Variable Importance 
occ.Xanthomonas.campestris.pred <- predict(Xanthomonas.campestris.me, modelEnv)
plot(occ.Xanthomonas.campestris.pred)
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Xanthomonas.campestris.me, p=occ.Xanthomonas.campestris.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
Xanthomonas.campestris.pred <- predict(Xanthomonas.campestris.me, modelEnv)
plot(Xanthomonas.campestris.pred)
#writeRaster(Xanthomonas.campestris.pred, file = "SDMs/Soybean_Xanthomonas.campestris.tif")

#Rhizoctonia.solani SDM 
occ.Rhizoctonia.solani<-Rhizoctonia.solani.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)
fold <- kfold(occ.Rhizoctonia.solani, k=5) # add an index that makes five random groups of observations
occ.Rhizoctonia.solani.test <- occ.Rhizoctonia.solani[fold == 1, ] # hold out one fifth as test data
occ.Rhizoctonia.solani.train <- occ.Rhizoctonia.solani.test[fold != 1, ] # the other four fifths are training data
maxent()
Rhizoctonia.solani.me <- maxent(modelEnv, occ.Rhizoctonia.solani.train ) # just using the training data
print(Rhizoctonia.solani.me)
plot(Rhizoctonia.solani.me) #Variable Importance 
occ.Rhizoctonia.solani.pred <- predict(Rhizoctonia.solani.me, modelEnv)
plot(occ.Rhizoctonia.solani.pred)
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Rhizoctonia.solani.me, p=occ.Rhizoctonia.solani.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
#writeRaster(occ.Rhizoctonia.solani.pred, file = "SDMs/Soybean_Rhizoctonia.solani.tif")

#Spodoptera.exigua SDM 
occ.Spodoptera.exigua<-Spodoptera.exigua.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)
fold <- kfold(occ.Spodoptera.exigua, k=5) # add an index that makes five random groups of observations
occ.Spodoptera.exigua.test <- occ.Spodoptera.exigua[fold == 1, ] # hold out one fifth as test data
occ.Spodoptera.exigua.train <- occ.Spodoptera.exigua.test[fold != 1, ] # the other four fifths are training data
maxent()
Spodoptera.exigua.me <- maxent(modelEnv, occ.Spodoptera.exigua.train ) # just using the training data
print(Spodoptera.exigua.me)
plot(Spodoptera.exigua.me) #Variable Importance 
occ.Spodoptera.exigua.pred <- predict(Spodoptera.exigua.me, modelEnv)
plot(occ.Spodoptera.exigua.pred)
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(occ.Spodoptera.exigua.me, p=occ.Spodoptera.exigua.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
#writeRaster(occ.Spodoptera.exigua.pred, file = "SDMs/Soybean_Spodoptera.exigua.tif")

#Diabrotica.undecimpunctata SDM 
occ.Diabrotica.undecimpunctata<-Diabrotica.undecimpunctata.xy[,-1]
model.extent<-extent(-180,180,-90,90)
modelEnv=crop(currentEnv,model.extent)
fold <- kfold(occ.Diabrotica.undecimpunctata, k=5) # add an index that makes five random groups of observations
occ.Diabrotica.undecimpunctata.test <- occ.Diabrotica.undecimpunctata[fold == 1, ] # hold out one fifth as test data
occ.Diabrotica.undecimpunctata.train <- occ.Diabrotica.undecimpunctata.test[fold != 1, ] # the other four fifths are training data
maxent()
Diabrotica.undecimpunctata.me <- maxent(modelEnv, occ.Diabrotica.undecimpunctata.train ) # just using the training data
print(Diabrotica.undecimpunctata.me)
plot(Diabrotica.undecimpunctata.me) #Variable Importance 
occ.Diabrotica.undecimpunctata.pred <- predict(Diabrotica.undecimpunctata.me, modelEnv)
plot(occ.Diabrotica.undecimpunctata.pred)
bg <- randomPoints(modelEnv, 1000) #make psuedorandom background points
e1 <- evaluate(Diabrotica.undecimpunctata.me, p=occ.Diabrotica.undecimpunctata.test, a=bg, x=modelEnv)
plot(e1, 'ROC')
#writeRaster(occ.Diabrotica.undecimpunctata.pred, file = "SDMs/Soybean_Diabrotica.undecimpunctata.tif")




