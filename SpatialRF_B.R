SpatialRF_B
library(spatialRF)
library(spdep)
library(spdplyr)
library(rgdal)

#Read In
#Fishnet- pixel = 100km^2
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")

#Set Projection
proj4string(fishnet) <- CRS("+init=epsg:3786") 

#Some covariates
barleymodel<-read.csv(file="Models/Barley_RF.csv", header=T)
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID #Fixing name

#Merging fishnet with covariates
m <- sp::merge(fishnet, barleymodel, by='Fishnet_ID')
barley<-sp::merge(m, mapping, by='Fishnet_ID')
barley.p = subset(barley, mean_barle > 0) #selecting values > 0 

#Transform numerical covariates
barley.num<-barley.p@data %>% select(mean_barle, AET_mean, Barley_Fertilizer, Pesticide,GDP_Mean) 
row.names(barley.num)<-barley.p$FISHNET_ID
barley.scaled<-na.omit(data.frame(log10(barley.num+1)))
barley.scaled<- tibble::rownames_to_column(barley.scaled, "FISHNET_ID")

#Select categorical covariates
barley.cat<-barley.p %>% select(COUNTRY.x,BarleyGCO, FISHNET_ID,Latitude,Longitude)

#Count number of countries in case I want to remove tiny ones
county.pivot<-barley.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
barley.catagory<-sp::merge(barley.cat,county.pivot, by="COUNTRY.x")

#Combine data
barley.comb<-sp::merge(barley.catagory,barley.scaled, by="FISHNET_ID")

#Counties with at least 1000km  coverage
barley.final<-na.omit(filter(barley.comb,CountyPixelCount > 10))
barley.final<-na.omit(barley.final@data)
summary(barley.final)
barley.final$BarleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$COUNTRY.x<-as.factor(barley.final$COUNTRY.x)
barley.coords<-barley.final %>% select(Latitude,Longitude) #select coordinates

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(fishnet$Longitude,fishnet$Latitude)))
dim(barley.final) #check dims
dim(dist.mat) #check dims

#Parameters 
X <- c("COUNTRY.x", "BarleyGCO", "AET_mean", "Barley_Fertilizer", "Pesticide", "GDP_Mean")
Y <-"mean_barle"


#coordinates of the cases
xy <- barley.coords
#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1000, 2000, 4000, 8000)
#random seed for reproducibility
random.seed <- 1
# remove data from spatial dataframe
barley.data<-barley.final

model.non.spatial <- spatialRF::rf(
  data = barley.data,
  dependent.variable.name = "mean_barle",
  predictor.variable.names = c("Barley_Fertilizer", "Pesticide"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = FALSE)



