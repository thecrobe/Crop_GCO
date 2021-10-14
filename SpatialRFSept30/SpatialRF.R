#Read In
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
proj4string(fishnet) <- CRS("+init=epsg:3786")
barleymodel<-read.csv(file="Models/Barley_RF.csv", header=T)
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID
m <- sp::merge(fishnet, barleymodel, by='Fishnet_ID')
barley<-sp::merge(m, mapping, by='Fishnet_ID')
barley.p = subset(barley, mean_barle > 0)
dim(barley.p)
barley.num<-barley.p@data %>% select(mean_barle, AET_mean, Barley_Fertilizer, Pesticide,GDP_Mean)
barley.scaled<-data.frame(scale(barley.num,center = TRUE,scale = TRUE))
barley.cat<-barley.p@data %>% select(COUNTRY.x,BarleyGCO, FISHNET_ID,Latitude,Longitude)

barley.final %>%
  group_by(COUNTRY.x) %>% 
  summarize(CountyPixelCount = n())

lobsters %>%
  group_by(year) %>%
  summarize(count_by_year = n())

dim(barley.scaled)
dim(barley.cat)
barley.comb<-cbind(barley.scaled,barley.cat)
barley.final<-na.omit(barley.comb)
barley.coords<-barley.final %>% select(Latitude,Longitude)

dim(barley.final)
dim(barley.coords)

#Random Forest Model
barley.spatialrf<- grf(mean_barle~Barley_Fertilizer+GDP_Mean+Pesticide+AET_mean+BarleyGCO+COUNTRY.x, na.omit(barley.final), bw=10, kernel="adaptive", coords=barley.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)


print(barley.spatialrf)
warnings()
