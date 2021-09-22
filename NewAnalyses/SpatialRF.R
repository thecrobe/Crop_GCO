library(SpatialML)
library(ggplot2)
library(dplyr)
library(rgdal)

#Read In
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")

#Data Prep 
rf.barley<-read.csv(file="Models/Barley_RF.csv", header = T)   
fish.barley<-subset(fishnet, fishnet$mean_barle > 0 ) #yield > 0 
fish.barley.subset<-fish.barley@data %>% select(Fishnet_ID,mean_barle, Latitude, Longitude)
joined<-na.omit(inner_join(fish.barley.subset,rf.barley, by="Fishnet_ID"))
barley.data<-log10(joined %>% select(mean_barle,Barley_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
barley.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
barley.spatialrf<- grf(mean_barle~Barley_Fertilizer+GDP_Mean+Pesticide+AET_mean , barley.data, bw=100, kernel="adaptive", coords=barley.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(barley.spatialrf$LocalModelSummary)

#Prediction and Error
barley.predict<-predict(barley.spatialrf$Global.Model,barley.data)
pred.error<-data.frame(10^((barley.predict-barley.data$mean_barle)-1))
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Barley_Error.csv")


#Data Prep 
rf.cassava<-read.csv(file="Models/Cassava_RF.csv", header = T)   
fish.cassava<-subset(fishnet, fishnet$mean_cassa> 0 ) #yield > 0 
fish.cassava.subset<-fish.cassava@data %>% select(Fishnet_ID,mean_cassa, Latitude, Longitude)
joined<-na.omit(inner_join(fish.cassava.subset,rf.cassava, by="Fishnet_ID"))
cassava.data<-log10(joined %>% select(mean_cassa,cassava_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
cassava.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
cassava.spatialrf<- grf(mean_cassa~cassava_Fertilizer+GDP_Mean+Pesticide+AET_mean , cassava.data, bw=100, kernel="adaptive", coords=cassava.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(cassava.spatialrf$LocalModelSummary)

#Prediction and Error
cassava.predict<-predict(cassava.spatialrf$Global.Model,cassava.data)
pred.error<-data.frame(10^((cassava.predict-cassava.data$mean_cassa))-1)
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Cassava_Error.csv")


#Data Prep 
rf.groundnut<-read.csv(file="Models/Groundnut_RF.csv", header = T)  
fish.groundnut<-subset(fishnet, fishnet$mean_groun > 0 ) #yield > 0 
fish.groundnut.subset<-fish.groundnut@data %>% select(Fishnet_ID,mean_groun, Latitude, Longitude)
joined<-na.omit(inner_join(fish.groundnut.subset,rf.groundnut, by="Fishnet_ID"))
groundnut.data<-log10(joined %>% select(mean_groun,groundnut_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
groundnut.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
groundnut.spatialrf<- grf(mean_groun~groundnut_Fertilizer+GDP_Mean+Pesticide+AET_mean , groundnut.data, bw=100, kernel="adaptive", coords=groundnut.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(groundnut.spatialrf$LocalModelSummary)

#Prediction and Error
groundnut.predict<-predict(groundnut.spatialrf$Global.Model,groundnut.data)
pred.error<-data.frame(10^((groundnut.predict-groundnut.data$mean_groun))-1)
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Groundnut_Error.csv")



#Data Prep 
rf.maize<-read.csv(file="Models/Maize_RF.csv", header = T)  
fish.maize<-subset(fishnet, fishnet$mean_maize > 0 ) #yield > 0 
fish.maize.subset<-fish.maize@data %>% select(Fishnet_ID,mean_maize, Latitude, Longitude)
joined<-na.omit(inner_join(fish.maize.subset,rf.maize, by="Fishnet_ID"))
maize.data<-log10(joined %>% select(mean_maize,maize_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
maize.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
maize.spatialrf<- grf(mean_maize~maize_Fertilizer+GDP_Mean+Pesticide+AET_mean , maize.data, bw=100, kernel="adaptive", coords=maize.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)


print(maize.spatialrf$LocalModelSummary)

#Prediction and Error
maize.predict<-predict(maize.spatialrf$Global.Model,maize.data)
pred.error<-data.frame(10^((maize.predict-maize.data$mean_maize))-1)
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Maize_Error.csv")


#Data Prep 
rf.potato<-read.csv(file="Models/Potato_RF.csv", header = T)  
fish.potato<-subset(fishnet, fishnet$mean_potat > 0 ) #yield > 0 
fish.potato.subset<-fish.potato@data %>% select(Fishnet_ID,mean_potat, Latitude, Longitude)
joined<-na.omit(inner_join(fish.potato.subset,rf.potato, by="Fishnet_ID"))
potato.data<-log10(joined %>% select(mean_potat.x,Potato_Fe,GDP_Mean,Pesticide,AET_mean)+1)
potato.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
potato.spatialrf<- grf(mean_potat.x~Potato_Fe+GDP_Mean+Pesticide+AET_mean , potato.data, bw=100, kernel="adaptive", coords=potato.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(potato.spatialrf$LocalModelSummary)

#Prediction and Error
potato.predict<-predict(potato.spatialrf$Global.Model,potato.data)
pred.error<-data.frame(10^((potato.predict-potato.data$mean_potat.x))-1)
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Potato_Error.csv")



#Data Prep 
rf.rapes<-read.csv(file="Models/Rapeseed_RF.csv", header = T)  
fish.rapes<-subset(fishnet, fishnet$mean_rapes > 0 ) #yield > 0 
fish.rapes.subset<-fish.rapes@data %>% select(Fishnet_ID,mean_rapes, Latitude, Longitude)
joined<-na.omit(inner_join(fish.rapes.subset,rf.rapes, by="Fishnet_ID"))
rapes.data<-log10(joined %>% select(mean_rapes,rapeseed_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
rapes.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
rapes.spatialrf<- grf(mean_rapes~rapeseed_Fertilizer+GDP_Mean+Pesticide+AET_mean , rapes.data, bw=100, kernel="adaptive", coords=rapes.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(rapes.spatialrf$LocalModelSummary)

#Prediction and Error
rapes.predict<-predict(rapes.spatialrf$Global.Model,rapes.data)
pred.error<-data.frame(10^((rapes.predict-rapes.data$mean_rapes))-1)
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Rapes_Error.csv")

#Data Prep 
rf.rice<-read.csv(file="Models/Rice_RF.csv", header = T)  
fish.rice<-subset(fishnet, fishnet$mean_rice > 0 ) #yield > 0 
fish.rice<-subset(fishnet, fishnet$mean_rice_ > 0 ) #yield > 0 
fish.rice.subset<-fish.rice@data %>% select(Fishnet_ID,mean_rice_, Latitude, Longitude)
joined<-na.omit(inner_join(fish.rice.subset,rf.rice, by="Fishnet_ID"))
rice.data<-log10(joined %>% select(mean_rice_,Rice_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)

rice.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))
#Random Forest Model
rice.spatialrf<- grf(mean_rice_~rice_Fertilizer+GDP_Mean+Pesticide+AET_mean , rice.data, bw=100, kernel="adaptive", coords=rice.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(rice.spatialrf$LocalModelSummary)

#Prediction and Error
rice.predict<-predict(rice.spatialrf$Global.Model,rice.data)
pred.error<-data.frame(10^((rice.predict-rice.data$mean_rice))-1)
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Rice_Error.csv")


#Data Prep 
rf.rye<-read.csv(file="Models/Rye_RF.csv", header = T)  
fish.rye<-subset(fishnet, fishnet$mean_rye_Y > 0 ) #yield > 0 
fish.rye.subset<-fish.rye@data %>% select(Fishnet_ID,mean_rye_Y, Latitude, Longitude)
joined<-na.omit(inner_join(fish.rye.subset,rf.rye, by="Fishnet_ID"))
rye.data<-log10(joined %>% select(mean_rye_Y,rye_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
rye.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
rye.spatialrf<- grf(mean_rye_Y~rye_Fertilizer+GDP_Mean+Pesticide+AET_mean , rye.data, bw=100, kernel="adaptive", coords=rye.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(rye.spatialrf$LocalModelSummary)

#Prediction and Error
rye.predict<-predict(rye.spatialrf$Global.Model,rye.data)
pred.error<-data.frame(10^((rye.predict-rye.data$mean_rye))-1)
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Rye_Error.csv")


#Data Prep 
rf.sorgh<-read.csv(file="Models/Sorghum_RF.csv", header = T)  
fish.sorgh<-subset(fishnet, fishnet$mean_sorgh > 0 ) #yield > 0 
fish.sorgh.subset<-fish.sorgh@data %>% select(Fishnet_ID,mean_sorgh, Latitude, Longitude)
joined<-na.omit(inner_join(fish.sorgh.subset,rf.sorgh, by="Fishnet_ID"))
sorgh.data<-log10(joined %>% select(mean_sorgh,sorghum_Fertilizer,GDP,Pesticide,Evapotranspiraton)+1)
sorgh.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
sorgh.spatialrf<- grf(mean_sorgh~sorghum_Fertilizer+GDP+Pesticide+Evapotranspiration , sorgh.data, bw=100, kernel="adaptive", coords=sorgh.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(sorgh.spatialrf$LocalModelSummary)

#Prediction and Error
sorgh.predict<-predict(sorgh.spatialrf$Global.Model,sorgh.data)
pred.error<-data.frame(10^((sorgh.predict-sorgh.data$mean_sorgh)-1)
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Sorgh_Error.csv")


#Data Prep 
rf.soybe<-read.csv(file="Models/Soybean_RF.csv", header = T)  
fish.soybe<-subset(fishnet, fishnet$mean_soybe > 0 ) #yield > 0 
fish.soybe.subset<-fish.soybe@data %>% select(Fishnet_ID,mean_soybe, Latitude, Longitude)
joined<-na.omit(inner_join(fish.soybe.subset,rf.soybe, by="Fishnet_ID"))
soybe.data<-log10(joined %>% select(mean_soybe,soybean_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
soybe.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
soybe.spatialrf<- grf(mean_soybe~soybean_Fertilizer+GDP_Mean+Pesticide+AET_mean , soybe.data, bw=100, kernel="adaptive", coords=soybe.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(soybe.spatialrf$LocalModelSummary)

#Prediction and Error
soybe.predict<-predict(soybe.spatialrf$Global.Model,soybe.data)
pred.error<-data.frame(10^((soybe.predict-soybe.data$mean_soybe)-1)
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Soybe_Error.csv")


#Data Prep 
rf.sunfl<-read.csv(file="Models/Sunflower_RF.csv", header = T)  
fish.sunfl<-subset(fishnet, fishnet$mean_sunfl > 0 ) #yield > 0 
fish.sunfl.subset<-fish.sunfl@data %>% select(Fishnet_ID,mean_sunfl, Latitude, Longitude)
joined<-na.omit(inner_join(fish.sunfl.subset,rf.sunfl, by="Fishnet_ID"))
sunfl.data<-log10(joined %>% select(mean_sunfl,sunflower_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
sunfl.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
sunfl.spatialrf<- grf(mean_sunfl~sunflower_Fertilizer+GDP_Mean+Pesticide+AET_mean , sunfl.data, bw=100, kernel="adaptive", coords=sunfl.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(sunfl.spatialrf$LocalModelSummary)

#Prediction and Error
sunfl.predict<-predict(sunfl.spatialrf$Global.Model,sunfl.data)
pred.error<-data.frame(10^((sunfl.predict-sunfl.data$mean_sunfl)-1))
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Sunfl_Error.csv")

#Data Prep 
rf.wheat<-read.csv(file="Models/Wheat_RF.csv", header = T)  
fish.wheat<-subset(fishnet, fishnet$mean_wheat > 0 ) #yield > 0 
fish.wheat.subset<-fish.wheat@data %>% select(Fishnet_ID,mean_wheat, Latitude, Longitude)
joined<-na.omit(inner_join(fish.wheat.subset,rf.wheat, by="Fishnet_ID"))
wheat.data<-log10(joined %>% select(mean_wheat,wheat_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
wheat.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
wheat.spatialrf<- grf(mean_wheat~wheat_Fertilizer+GDP_Mean+Pesticide+AET_mean , wheat.data, bw=100, kernel="adaptive", coords=wheat.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(wheat.spatialrf$LocalModelSummary)

#Prediction and Error
wheat.predict<-predict(wheat.spatialrf$Global.Model,wheat.data)
pred.error<-data.frame(10^((wheat.predict-wheat.data$mean_wheat)-1))
summary(pred.error) 
write.csv(x = pred.error,file="SpatialRF/Wheat_Error.csv")

nb <- poly2nb(fish.wheat)
lw <- nb2listw(nb, zero.policy = TRUE)
moran.mc(pred.error$X10...wheat.predict...wheat.data.mean_wheat....1., lw, 999,zero.policy = TRUE) #Test for autocorrelation                       

wheat.fish<-fish.wheat
wheat.fish@data<- wheat.fish@data$%>% 
  subset(1,3) #keeps column 1 and column 3 in the spdf object.
dim(fish.wheat)
dim(pred.error)

mer

library(sp)
