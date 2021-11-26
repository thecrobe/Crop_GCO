library(dplyr)
library(rgdal)
library(ggplot2)
library(wesanderson)
library(rgdal)
library(spdplyr)
library(sp)
library(rgeos)

#Read In
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_Null")
gco<-read.csv(file="Fishnets/GCO_Mapping.csv",header = T)
proj4string(fishnet) <- CRS("+init=epsg:3786")
data<-fishnet@data
data$FISHNET_ID<-data$Fishnet_ID


null<-inner_join(gco,data, by="FISHNET_ID")
null[null==0] <- NA
summary(null)

#Barley
ggplot(data, aes(x=null$Bar_null, y=null$BarleyGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Cassava
ggplot(data, aes(x=null$Cas_null, y=null$CassavaGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Groundnut
ggplot(data, aes(x=null$Groun_null, y=null$GroundnutGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Maize
ggplot(data, aes(x=null$Maiz_null, y=null$MaizeGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Potato
ggplot(data, aes(x=null$Potat_null, y=null$PotatoGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Rapeseed
ggplot(data, aes(x=null$Rap_null, y=null$RapeseedGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 


#Rice
ggplot(data, aes(x=null$Rice_null, y=null$RiceGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Rye
ggplot(data, aes(x=null$Rye_null, y=null$RyeGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Sorghum
ggplot(data, aes(x=null$Sorg_null, y=null$SorghumGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Soybean
ggplot(data, aes(x=null$Soy_null, y=null$SoybeanGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Sunflower
ggplot(data, aes(x=null$Sunf_null, y=null$SunflowerGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 

#Wheat
ggplot(data, aes(x=null$Wheat_null, y=null$WheatGCO)) + 
  geom_point() + 
  xlab("Slope") + 
  ylab("GCO") +
  theme_justin 


######### Distances to max pixel 
## Barley
max<-max(na.omit(data$Bar_null))
null.barley<-fishnet %>% filter(Bar_null == max ) #yield > 0 
xy <- coordinates(null.barley)
xy1<-data.frame(xy)
barleyGCO<-readOGR(dsn = "GIS/",layer="BarleyWheat_GC")
barleyGCOcentroid<-gCentroid(barleyGCO)
distance <- spDistsN1(xy,barleyGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 

## Cassava
max<-max(na.omit(data$Cas_null))
null.cas<-fishnet %>% filter(Cas_null == max ) #yield > 0 
xy <- coordinates(null.cas)
xy1<-data.frame(xy)
cassavaGCO<-readOGR(dsn = "GIS/",layer="Cassava_GC")
cassavaGCOcentroid<-gCentroid(cassavaGCO)
distance <- spDistsN1(xy,cassavaGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 

## Groundnut
max<-max(na.omit(data$Groun_null))
null.groun<-fishnet %>% filter(Groun_null == max ) #yield > 0 
xy <- coordinates(null.groun)
xy1<-data.frame(xy)
groundnutGCO<-readOGR(dsn = "GIS/",layer="Groundnut_GC")
groundnutGCOcentroid<-gCentroid(groundnutGCO)
distance <- spDistsN1(xy,groundnutGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 

## Maize
max<-max(na.omit(data$Maiz_null))
null.maize<-fishnet %>% filter(Maiz_null == max ) #yield > 0 
xy <- coordinates(null.maize)
xy1<-data.frame(xy)
maizeGCO<-readOGR(dsn = "GIS/",layer="Maize_GC")
maizeGCOcentroid<-gCentroid(maizeGCO)
distance <- spDistsN1(xy,maizeGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 

## Wheat
max<-max(na.omit(data$Wheat_null))
null.wheat<-fishnet %>% filter(Wheat_null == max ) #yield > 0 
xy <- coordinates(null.wheat)
xy1<-data.frame(xy)
wheatGCO<-readOGR(dsn = "GIS/",layer="BarleyWheat_GC")
wheatGCOcentroid<-gCentroid(wheatGCO)
distance <- spDistsN1(xy,wheatGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 


## Sunflower
max<-max(na.omit(data$Sunf_null))
null.sunf<-fishnet %>% filter(Sunf_null == max ) #yield > 0 
xy <- coordinates(null.sunf)
xy1<-data.frame(xy)
sunflowerGCO<-readOGR(dsn = "GIS/",layer="Sunflower_GC")
sunflowerGCOcentroid<-gCentroid(sunflowerGCO)
distance <- spDistsN1(xy,sunflowerGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 


## Soybean
max<-max(na.omit(data$Soy_null))
null.soy<-fishnet %>% filter(Soy_null == max ) #yield > 0 
xy <- coordinates(null.soy)
xy1<-data.frame(xy)
soyGCO<-readOGR(dsn = "GIS/",layer="Soy_GC")
soyGCOcentroid<-gCentroid(soyGCO)
distance <- spDistsN1(xy,soyGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 

## Sorghum
max<-max(na.omit(data$Sorg_null))
null.sorg<-fishnet %>% filter(Sorg_null == max ) #yield > 0 
xy <- coordinates(null.sorg)
xy1<-data.frame(xy)
sorgGCO<-readOGR(dsn = "GIS/",layer="Sorghum_GC")
sorgGCOcentroid<-gCentroid(sorgGCO)
distance <- spDistsN1(xy,sorgGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 

## Rye
max<-max(na.omit(data$Rye_null))
null.rye<-fishnet %>% filter(Rye_null == max ) #yield > 0 
xy <- coordinates(null.rye)
xy1<-data.frame(xy)
ryeGCO<-readOGR(dsn = "GIS/",layer="Rye")
ryeGCOcentroid<-gCentroid(ryeGCO)
distance <- spDistsN1(xy,ryeGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 

## Rapeseed
max<-max(na.omit(data$Rap_null))
null.rap<-fishnet %>% filter(Rap_null == max ) #yield > 0 
xy <- coordinates(null.rap)
xy1<-data.frame(xy)
rapGCO<-readOGR(dsn = "GIS/",layer="Rapeseed")
rapGCOcentroid<-gCentroid(rapGCO)
distance <- spDistsN1(xy,rapGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 

## Potato
max<-max(na.omit(data$Potat_null))
null.pot<-fishnet %>% filter(Potat_null == max ) #yield > 0 
xy <- coordinates(null.pot)
xy1<-data.frame(xy)
potGCO<-readOGR(dsn = "GIS/",layer="Potato")
potGCOcentroid<-gCentroid(potGCO)
distance <- spDistsN1(xy,potGCOcentroid, longlat = FALSE)
print(distance/1000) #KMs away form GCO 

