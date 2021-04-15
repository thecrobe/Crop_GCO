

library("rgbif")
# Collect occurrences for pests

#Achatina.fulica
Achatina.fulica.gbf<-occ_search(scientificName = "Achatina fulica",hasCoordinate = TRUE)
Achatina.fulica<-data.frame(Achatina.fulica.gbf$data)
plot(Achatina.fulica.map)
# Number of occurrences
dim(Achatina.fulica) #512 

#Claviceps.africana
Claviceps.africana.gbf<-occ_search(scientificName = "Claviceps africana",hasCoordinate = TRUE)
Claviceps.africana<-data.frame(Claviceps.africana.gbf$data)
# Number of occurrences
dim(Claviceps.africana) #13

#Curvularia clavata - NOT PRESENT IN GBIF
#Curvularia.clavata.gbf<-occ_search(scientificName = "Claviceps clavata",hasCoordinate = TRUE)

#Eurystylus.oldi
Eurystylus.oldi.gbf<-occ_search(scientificName = "Eurystylus oldi",hasCoordinate = TRUE)
Eurystylus.oldi<-data.frame(Eurystylus.oldi.gbf$data)
# Number of occurrences
dim(Eurystylus.oldi) #4

#Maliarpha.separatella
Maliarpha.separatella.gbf<-occ_search(scientificName = " Maliarpha separatella ",hasCoordinate = TRUE)
Maliarpha.separatella<-data.frame(Maliarpha.separatella.gbf$data)
# Number of occurrences
dim(Maliarpha.separatella) #5

#Paspalum urvillei
Paspalum.urvillei.gbf<-occ_search(scientificName = "Paspalum urvillei",hasCoordinate = TRUE)
Paspalum.urvillei<-data.frame(Paspalum.urvillei.gbf$data)
# Number of occurrences
dim(Paspalum.urvillei) #500

#Peronosclerospora australiensis - NOT PRESENT 
#Peronosclerospora.australiensis.gbf<-occ_search(scientificName = "Peronosclerospora australiensis",hasCoordinate = TRUE)

#Peronosclerospora sargae - Not Present
#Peronosclerospora.sargae.gbf<-occ_search(scientificName = "Peronosclerospora sargae",hasCoordinate = TRUE)

#Spodoptera frugiperda
Spodoptera.frugiperda.gbf<-occ_search(scientificName = "Spodoptera frugiperda",hasCoordinate = TRUE)
Spodoptera.frugiperda<-data.frame(Spodoptera.frugiperda.gbf$data)
# Number of occurrences
dim(Spodoptera.frugiperda) #500

#Verbesina encelioides
Verbesina.encelioides.gbf<-occ_search(scientificName = "Verbesina encelioides",limit = 10000)
Verbesina.encelioides<-data.frame(Verbesina.encelioides.gbf$data)
# Number of occurrences
dim(Verbesina.encelioides) #500

