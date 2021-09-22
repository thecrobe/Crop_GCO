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
pred.error<-data.frame(10^(barley.predict-barley.data$mean_barle))
summary(pred.error) 
#write.csv(x = pred.error,file="SpatialRF/Barley_Error.csv")


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
pred.error<-data.frame(10^(cassava.predict-cassava.data$mean_cassa))
summary(pred.error) 
#write.csv(x = pred.error,file="SpatialRF/Cassava_Error.csv")


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
pred.error<-data.frame(10^(groundnut.predict-groundnut.data$mean_groun))
summary(pred.error) 
#write.csv(x = pred.error,file="SpatialRF/Groundnut_Error.csv")



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
pred.error<-data.frame(10^(maize.predict-maize.data$mean_maize))
summary(pred.error) 
#write.csv(x = pred.error,file="SpatialRF/Maize_Error.csv")


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
pred.error<-data.frame(10^(potato.predict-potato.data$mean_potat.x))
summary(pred.error) 
#write.csv(x = pred.error,file="SpatialRF/Potato_Error.csv")



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
pred.error<-data.frame(10^(rapes.predict-rapes.data$mean_rapes))
summary(pred.error) 
#write.csv(x = pred.error,file="SpatialRF/Rapes_Error.csv")


#Data Prep 
rf.rice<-read.csv(file="Models/Rice_RF.csv", header = T)  
fish.rice<-subset(fishnet, fishnet$mean_rice_ > 0 ) #yield > 0 
fish.rice.subset<-fish.rice@data %>% select(Fishnet_ID,mean_rice_, Latitude, Longitude)
joined<-na.omit(inner_join(fish.rice.subset,rf.rice, by="Fishnet_ID"))
rice.data<-log10(joined %>% select(mean_rice_,rice_Fertilizer,GDP_Mean,Pesticide,AET_mean)+1)
rice.coords<-data.frame(cbind(joined$Latitude,joined$Longitude))

#Random Forest Model
rice.spatialrf<- grf(mean_rice_~rice_Fertilizer+GDP_Mean+Pesticide+AET_mean , rice.data, bw=100, kernel="adaptive", coords=rice.coords , ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)

print(rice.spatialrf$LocalModelSummary)

#Prediction and Error
rice.predict<-predict(rice.spatialrf$Global.Model,rice.data)
pred.error<-data.frame(10^(rice.predict-rice.data$mean_rice))
summary(pred.error) 
#write.csv(x = pred.error,file="SpatialRF/Rice_Error.csv")


