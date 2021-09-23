
#models
library(spatialEco)
library(rgdal)
library(rgeos)
library(spdep)
library(spatialreg)

#plots
library(RColorBrewer)
library(ggplot2)
library(wesanderson)
library(latticeExtra)
library(RColorBrewer)

#data wrangling
library(reshape2)
library(dplyr)


### Graphics

theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank())

fishnet<- readOGR(dsn= "GIS/", layer="FertNear")
proj4string(fishnet) <- CRS("+init=epsg:3786")

#Barley 
fert.barley<-subset(fishnet, fishnet$barley_Fer > 0 ) #yield > 0 
xy <- coordinates(fert.barley)
xy1<-data.frame(xy)
fert <- fert.barley@data
barleyGCO<-readOGR(dsn = "GIS/",layer="BarleyWheat_GC")
barleyGCOcentroid<-gCentroid(barleyGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,barleyGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
barley.fert<-fert %>% dplyr::select(barley_Fer,rescale_ND,COUNTRY)
barley.fert<-barley.fert %>% filter(barley_Fer > 0, na.rm=TRUE)
barley.fert<-barley.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
barley.fert$logBarleyfert<-log10(barley.fert$barley_Fer)


#Links
nb <- poly2nb(fert.barley)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.barley$logBarleyFert<-log10(fert.barley$barley_Fer)
#Spatial Error Model
barley.fert.sper <- errorsarlm(logBarleyFert ~ barley.fert$rescale_ND, data=fert.barley, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(barley.fert.sper) #Summary
fert.barley$residualsSpecError <- residuals(barley.fert.sper) #Residuals
moran.mc(fert.barley$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.barley$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Fertilizer Applications by Country Sum 
barley.pivot<-dcast(fert, COUNTRY ~., value.var="barley_Fer", fun.aggregate=sum)
barley.pivot<-arrange(barley.pivot,.)
tail(barley.pivot, 3) #ID Top 3 Fertilizer Countries
sum(barley.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$barley_Fer)/sum(barley.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"
summary(china)

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$barley_Fer)/sum(barley.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

canada<-filter(fert, COUNTRY == "Canada")
canada.Perc<-sum(canada$barley_Fer)/sum(barley.pivot$.) 
print(canada.Perc) # % of Global Fertilizer For Crop 
canada.color<-"#98C24B"

#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(fert$barley_Fer))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$barley_Fer)), color=china.color, alpha=0.1) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$barley_Fer)), color=usa.color, alpha=0.2) + 
  geom_point(data=canada, aes(x=rescale_ND, y=log10(canada$barley_Fer)), color=canada.color, alpha=0.1) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(barley.fert.sper)[2], slope=coef(barley.fert.sper)[3]), color="red", size=2)

#Cassava 
fert.cassava<-subset(fishnet, fishnet$cassava_Fe > 0 ) #yield > 0 
xy <- coordinates(fert.cassava)
xy1<-data.frame(xy)
fert <- fert.cassava@data
cassavaGCO<-readOGR(dsn = "GIS/",layer="Cassava_GC")
cassavaGCOcentroid<-gCentroid(cassavaGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,cassavaGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
cassava.fert<-fert %>% dplyr::select(cassava_Fe,rescale_ND,COUNTRY)
cassava.fert<-cassava.fert %>% filter(cassava_Fe > 0, na.rm=TRUE)
cassava.fert<-cassava.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
cassava.fert$logCassavafert<-log10(cassava.fert$cassava_Fe)


#Links
nb <- poly2nb(fert.cassava)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.cassava$logCassavaFert<-log10(fert.cassava$cassava_Fe)
#Spatial Error Model
cassava.fert.sper <- errorsarlm(logCassavaFert ~ cassava.fert$rescale_ND, data=fert.cassava, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(cassava.fert.sper) #Summary
fert.cassava$residualsSpecError <- residuals(cassava.fert.sper) #Residuals
moran.mc(fert.cassava$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.cassava$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Fertilizer Applications by Country Sum 
cassava.pivot<-dcast(fert, COUNTRY ~., value.var="cassava_Fe", fun.aggregate=sum)
cassava.pivot<-arrange(cassava.pivot,.)
tail(cassava.pivot, 3) #ID Top 3 Fertilizer Countries
sum(cassava.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$cassava_Fe)/sum(cassava.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"
summary(china)

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$cassava_Fe)/sum(cassava.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

brazil<-filter(fert, COUNTRY == "Brazil")
brazil.Perc<-sum(brazil$cassava_Fe)/sum(cassava.pivot$.) 
print(brazil.Perc) # % of Global Fertilizer For Crop 
brazil.color<-"#E28D15"


#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(fert$cassava_Fe))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$cassava_Fe)), color=china.color, alpha=0.1) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$cassava_Fe)), color=usa.color, alpha=0.2) + 
  geom_point(data=canada, aes(x=rescale_ND, y=log10(brazil$cassava_Fe)), color=brazil.color, alpha=0.1) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(cassava.fert.sper)[2], slope=coef(cassava.fert.sper)[3]), color="grey", size=2)

#Groundnut 
fert.groundnut<-subset(fishnet, fishnet$groundnut_ > 0 ) #yield > 0 
xy <- coordinates(fert.groundnut)
xy1<-data.frame(xy)
fert <- fert.groundnut@data
groundnutGCO<-readOGR(dsn = "GIS/",layer="Groundnut_GC")
groundnutGCOcentroid<-gCentroid(groundnutGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,groundnutGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
groundnut.fert<-fert %>% dplyr::select(groundnut_,rescale_ND,COUNTRY)
groundnut.fert<-groundnut.fert %>% filter(groundnut_ > 0, na.rm=TRUE)
groundnut.fert<-groundnut.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
groundnut.fert$logGroundnutfert<-log10(groundnut.fert$groundnut_)


#Links
nb <- poly2nb(fert.groundnut)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.groundnut$logGroundnutFert<-log10(fert.groundnut$groundnut_)
#Spatial Error Model
groundnut.fert.sper <- errorsarlm(logGroundnutFert ~ groundnut.fert$rescale_ND, data=fert.groundnut, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(groundnut.fert.sper) #Summary
fert.groundnut$residualsSpecError <- residuals(groundnut.fert.sper) #Residuals
moran.mc(fert.groundnut$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.groundnut$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Fertilizer Applications by Country Sum 
groundnut.pivot<-dcast(fert, COUNTRY ~., value.var="groundnut_", fun.aggregate=sum)
groundnut.pivot<-arrange(groundnut.pivot,.)
tail(groundnut.pivot, 3) #ID Top 3 Fertilizer Countries
sum(groundnut.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$groundnut_)/sum(groundnut.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"
summary(china)

india<-filter(fert, COUNTRY == "India")
india.Perc<-sum(india$groundnut_)/sum(groundnut.pivot$.) 
print(india.Perc) # % of Global Fertilizer For Crop 
india.color<-"#078D40"

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$groundnut_)/sum(groundnut.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"


#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(groundnut_))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$groundnut_)), color=china.color, alpha=0.1) + 
  geom_point(data=india, aes(x=rescale_ND, y=log10(india$groundnut_)), color=india.color, alpha=0.2) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$groundnut_)), color=usa.color, alpha=0.1) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(groundnut.fert.sper)[2], slope=coef(groundnut.fert.sper)[3]), color="red", size=2)


#Maize 
fert.maize<-subset(fishnet, fishnet$maize_Fert > 0 ) #yield > 0 
xy <- coordinates(fert.maize)
xy1<-data.frame(xy)
fert <- fert.maize@data
maizeGCO<-readOGR(dsn = "GIS/",layer="Maize_GC")
maizeGCOcentroid<-gCentroid(maizeGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,maizeGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
maize.fert<-fert %>% dplyr::select(maize_Fert,rescale_ND,COUNTRY)
maize.fert<-maize.fert %>% filter(maize_Fert > 0, na.rm=TRUE)
maize.fert<-maize.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
maize.fert$logMaizefert<-log10(maize.fert$maize_Fert)


#Links
nb <- poly2nb(fert.maize)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.maize$logMaizeFert<-log10(fert.maize$maize_Fert)
#Spatial Error Model
maize.fert.sper <- errorsarlm(logMaizeFert ~ maize.fert$rescale_ND, data=fert.maize, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(maize.fert.sper) #Summary
fert.maize$residualsSpecError <- residuals(maize.fert.sper) #Residuals
moran.mc(fert.maize$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.maize$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Fertilizer Applications by Country Sum 
maize.pivot<-dcast(fert, COUNTRY ~., value.var="maize_Fert", fun.aggregate=sum)
maize.pivot<-arrange(maize.pivot,.)
tail(maize.pivot, 3) #ID Top 3 Fertilizer Countries
sum(maize.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$maize_Fert)/sum(maize.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$maize_Fert)/sum(maize.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

brazil<-filter(fert, COUNTRY == "Brazil")
brazil.Perc<-sum(brazil$maize_Fert)/sum(maize.pivot$.) 
print(brazil.Perc) # % of Global Fertilizer For Crop 
brazil.color<-"#E28D15"

#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(maize_Fert))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$maize_Fert)), color=china.color, alpha=0.1) + 
  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$maize_Fert)), color=brazil.color, alpha=0.2) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$maize_Fert)), color=usa.color, alpha=0.1) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(maize.fert.sper)[2], slope=coef(maize.fert.sper)[3]), color="grey", size=2)

summary(maize.fert.sper)

#Potato 
fert.Potato<-subset(fishnet, fishnet$Potato_Fe > 0 ) #yield > 0 
xy <- coordinates(fert.Potato)
xy1<-data.frame(xy)
fert <- fert.Potato@data
PotatoGCO<-readOGR(dsn = "GIS/",layer="Potato")
PotatoGCOcentroid<-gCentroid(PotatoGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,PotatoGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
Potato.fert<-fert %>% dplyr::select(Potato_Fe,rescale_ND,COUNTRY)
Potato.fert<-Potato.fert %>% filter(Potato_Fe > 0, na.rm=TRUE)
Potato.fert<-Potato.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
Potato.fert$logPotatofert<-log10(Potato.fert$Potato_Fe)

#Links
nb <- poly2nb(fert.Potato)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.Potato$logPotatoFert<-log10(fert.Potato$Potato_Fe)
#Spatial Error Model
Potato.fert.sper <- errorsarlm(logPotatoFert ~ Potato.fert$rescale_ND, data=fert.Potato, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(Potato.fert.sper) #Summary
fert.Potato$residualsSpecError <- residuals(Potato.fert.sper) #Residuals
moran.mc(fert.Potato$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.Potato$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals

#Pivot Table To ID Top Fertilizer Applications by Country Sum 
potato.pivot<-dcast(fert, COUNTRY ~., value.var="Potato_Fe", fun.aggregate=sum)
potato.pivot<-arrange(potato.pivot,.)
tail(potato.pivot, 3) #ID Top 3 Fertilizer Countries
sum(potato.pivot$.) #Total Fertilizer
#Filter Countries
canada<-filter(fert, COUNTRY == "Canada") 
canada.Perc<-sum(canada$Potato_Fe)/sum(potato.pivot$.) 
print(canada.Perc) # % of Global Fertilizer For Crop 
canada.color<-"#98C24B"

china<-filter(fert, COUNTRY == "China")
china.Perc<-sum(china$Potato_Fe)/sum(potato.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$Potato_Fe)/sum(potato.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

sum(china.Perc,usa.Perc,india.Perc) 

#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(Potato_Fe))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$Potato_Fe)), color=usa.color, alpha=0.1) + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$Potato_Fe)), color=china.color, alpha=0.2) + 
  geom_point(data=canada, aes(x=rescale_ND, y=log10(canada$Potato_Fe)), color=canada.color, alpha=0.3) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Potato.fert.sper)[2], slope=coef(Potato.fert.sper)[3]), color="red", size=2)

#Rapeseed 
fert.Rapeseed<-subset(fishnet, fishnet$rapeseed_F > 0 ) #yield > 0 
xy <- coordinates(fert.Rapeseed)
xy1<-data.frame(xy)
fert <- fert.Rapeseed@data
RapeseedGCO<-readOGR(dsn = "GIS/",layer="Rapeseed")
RapeseedGCOcentroid<-gCentroid(RapeseedGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,RapeseedGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
Rapeseed.fert<-fert %>% dplyr::select(rapeseed_F,rescale_ND,COUNTRY)
Rapeseed.fert<-Rapeseed.fert %>% filter(rapeseed_F > 0, na.rm=TRUE)
Rapeseed.fert<-Rapeseed.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
Rapeseed.fert$logRapeseedfert<-log10(Rapeseed.fert$rapeseed_F)


#Links
nb <- poly2nb(fert.Rapeseed)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.Rapeseed$logRapeseedFert<-log10(fert.Rapeseed$rapeseed_F)
#Spatial Error Model
Rapeseed.fert.sper <- errorsarlm(logRapeseedFert ~ Rapeseed.fert$rescale_ND, data=fert.Rapeseed, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(Rapeseed.fert.sper) #Summary
fert.Rapeseed$residualsSpecError <- residuals(Rapeseed.fert.sper) #Residuals
moran.mc(fert.Rapeseed$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.Rapeseed$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals

#Pivot Table To ID Top Fertilizer Applications by Country Sum 
rapeseed.pivot<-dcast(fert, COUNTRY ~., value.var="rapeseed_F", fun.aggregate=sum)
rapeseed.pivot<-arrange(rapeseed.pivot,.)
tail(rapeseed.pivot, 3) #ID Top 3 Fertilizer Countries
sum(rapeseed.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$rapeseed_F)/sum(rapeseed.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"

brazil<-filter(fert, COUNTRY == "Brazil")
brazil.Perc<-sum(brazil$rapeseed_F)/sum(rapeseed.pivot$.) 
print(brazil.Perc) # % of Global Fertilizer For Crop 
brazil.color<-"#E28D15"

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$rapeseed_F)/sum(rapeseed.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

sum(china.Perc,brazil.Perc,usa.Perc)

#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(rapeseed_F))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$rapeseed_F)), color=china.color, alpha=0.1) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$rapeseed_F)), color=usa.color, alpha=0.2) + 
  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$rapeseed_F)), color=brazil.color, alpha=0.3) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Rapeseed.fert.sper)[2], slope=coef(Rapeseed.fert.sper)[3]), color="red", size=2)


#Rice 
fert.rice<-subset(fishnet, fishnet$rice_FertM > 0 ) #yield > 0 
xy <- coordinates(fert.rice)
xy1<-data.frame(xy)
fert <- fert.rice@data
riceGCO<-readOGR(dsn = "GIS/",layer="Rice_GC1")
riceGCOcentroid<-gCentroid(riceGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,riceGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
rice.fert<-fert %>% dplyr::select(rice_FertM,rescale_ND,COUNTRY)
rice.fert<-rice.fert %>% filter(rice_FertM > 0, na.rm=TRUE)
rice.fert<-rice.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
rice.fert$logRicefert<-log10(rice.fert$rice_FertM)


#Links
nb <- poly2nb(fert.rice)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.rice$logRiceFert<-log10(fert.rice$rice_FertM)
#Spatial Error Model
rice.fert.sper <- errorsarlm(logRiceFert ~ rice.fert$rescale_ND, data=fert.rice, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rice.fert.sper) #Summary
fert.rice$residualsSpecError <- residuals(rice.fert.sper) #Residuals
moran.mc(fert.rice$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.rice$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Fertilizer Applications by Country Sum 
rice.pivot<-dcast(fert, COUNTRY ~., value.var="rice_FertM", fun.aggregate=sum)
rice.pivot<-arrange(rice.pivot,.)
tail(rice.pivot, 3) #ID Top 3 Fertilizer Countries
sum(rice.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$rice_FertM)/sum(rice.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"

india<-filter(fert, COUNTRY == "India")
india.Perc<-sum(india$rice_FertM)/sum(rice.pivot$.) 
print(india.Perc) # % of Global Fertilizer For Crop 
india.color<-"#078D40"

brazil<-filter(fert, COUNTRY == "Brazil")
brazil.Perc<-sum(brazil$rice_FertM)/sum(rice.pivot$.) 
print(brazil.Perc) # % of Global Fertilizer For Crop 
brazil.color<-"#E28D15"



#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(rice_FertM))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$rice_FertM)), color=china.color, alpha=0.1) + 
  geom_point(data=india, aes(x=rescale_ND, y=log10(india$rice_FertM)), color=india.color, alpha=0.2) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$rice_FertM)), color=usa.color, alpha=0.1) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(rice.fert.sper)[2], slope=coef(rice.fert.sper)[3]), color="red", size=2)

#Rye 
fert.Rye<-subset(fishnet, fishnet$rye_FertMe > 0 ) #yield > 0 
xy <- coordinates(fert.Rye)
xy1<-data.frame(xy)
fert <- fert.Rye@data
RyeGCO<-readOGR(dsn = "GIS/",layer="Rye")
RyeGCOcentroid<-gCentroid(RyeGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,RyeGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
Rye.fert<-fert %>% dplyr::select(rye_FertMe,rescale_ND,COUNTRY)
Rye.fert<-Rye.fert %>% filter(rye_FertMe > 0, na.rm=TRUE)
Rye.fert<-Rye.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
Rye.fert$logRyefert<-log10(Rye.fert$rye_FertMe)


#Links
nb <- poly2nb(fert.Rye)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.Rye$logRyeFert<-log10(fert.Rye$rye_FertMe)
#Spatial Error Model
Rye.fert.sper <- errorsarlm(logRyeFert ~ Rye.fert$rescale_ND, data=fert.Rye, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(Rye.fert.sper) #Summary
fert.Rye$residualsSpecError <- residuals(Rye.fert.sper) #Residuals
moran.mc(fert.Rye$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.Rye$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Fertilizer Applications by Country Sum 
rye.pivot<-dcast(fert, COUNTRY ~., value.var="rye_FertMe", fun.aggregate=sum)
rye.pivot<-arrange(rye.pivot,.)
tail(rye.pivot, 3) #ID Top 3 Fertilizer Countries
sum(rye.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$rye_Fert)/sum(rye.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$rye_Fert)/sum(rye.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

russia<-filter(fert, COUNTRY == "Russian Federation")
russia.Perc<-sum(russia$rye_FertMe)/sum(rye.pivot$.) 
print(russia.Perc) # % of Global Fertilizer For Crop 
russia.color<-"#D3867A"


#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(rye_FertMe))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$rye_FertMe)), color=china.color, alpha=0.1) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$rye_FertMe)), color=usa.color, alpha=0.2) + 
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$rye_FertMe)), color=russia.color, alpha=0.3) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Rye.fert.sper)[2], slope=coef(Rye.fert.sper)[3]), color="red", size=2)


#Sorghum 
fert.sorghum<-subset(fishnet, fishnet$sorghum_Fe > 0 ) #yield > 0 
xy <- coordinates(fert.sorghum)
xy1<-data.frame(xy)
fert <- fert.sorghum@data
sorghumGCO<-readOGR(dsn = "GIS/",layer="Sorghum_GC")
sorghumGCOcentroid<-gCentroid(sorghumGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,sorghumGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
sorghum.fert<-fert %>% dplyr::select(sorghum_Fe,rescale_ND,COUNTRY)
sorghum.fert<-sorghum.fert %>% filter(sorghum_Fe > 0, na.rm=TRUE)
sorghum.fert<-sorghum.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
sorghum.fert$logSorghumfert<-log10(sorghum.fert$sorghum_Fe)


#Links
nb <- poly2nb(fert.sorghum)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.sorghum$logSorghumFert<-log10(fert.sorghum$sorghum_Fe)
#Spatial Error Model
sorghum.fert.sper <- errorsarlm(logSorghumFert ~ sorghum.fert$rescale_ND, data=fert.sorghum, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(sorghum.fert.sper) #Summary
fert.sorghum$residualsSpecError <- residuals(sorghum.fert.sper) #Residuals
moran.mc(fert.sorghum$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.sorghum$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Fertilizer Applications by Country Sum 
sorghum.pivot<-dcast(fert, COUNTRY ~., value.var="sorghum_Fe", fun.aggregate=sum)
sorghum.pivot<-arrange(sorghum.pivot,.)
tail(sorghum.pivot, 3) #ID Top 3 Fertilizer Countries
sum(sorghum.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$sorghum_Fe)/sum(sorghum.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"
summary(china)

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$sorghum_Fe)/sum(sorghum.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

mexico<-filter(fert, COUNTRY == "Mexico")
mexico.Perc<-sum(mexico$sorghum_Fe)/sum(sorghum.pivot$.) 
print(mexico.Perc) # % of Global Fertilizer For Crop 
mexico.color<-"#96AEEA"

#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(sorghum_Fe))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$sorghum_Fe)), color=china.color, alpha=0.1) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$sorghum_Fe)), color=usa.color, alpha=0.2) + 
  geom_point(data=mexico, aes(x=rescale_ND, y=log10(mexico$sorghum_Fe)), color=mexico.color, alpha=0.1) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(sorghum.fert.sper)[2], slope=coef(sorghum.fert.sper)[3]), color="red", size=2)

#Soybean
fert.Soybean<-subset(fishnet, fishnet$soybean_Fe > 0 ) #yield > 0 
xy <- coordinates(fert.Soybean)
xy1<-data.frame(xy)
fert <- fert.Soybean@data
SoybeanGCO<-readOGR(dsn = "GIS/",layer="Soy_GC")
SoybeanGCOcentroid<-gCentroid(SoybeanGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,SoybeanGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
Soybean.fert<-fert %>% dplyr::select(soybean_Fe,rescale_ND,COUNTRY)
Soybean.fert<-Soybean.fert %>% filter(soybean_Fe > 0, na.rm=TRUE)
Soybean.fert<-Soybean.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
Soybean.fert$logSoybeanfert<-log10(Soybean.fert$soybean_Fe)


#Links
nb <- poly2nb(fert.Soybean)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.Soybean$logSoybeanFert<-log10(fert.Soybean$soybean_Fe)
#Spatial Error Model
Soybean.fert.sper <- errorsarlm(logSoybeanFert ~ Soybean.fert$rescale_ND, data=fert.Soybean, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(Soybean.fert.sper) #Summary
fert.Soybean$residualsSpecError <- residuals(Soybean.fert.sper) #Residuals
moran.mc(fert.Soybean$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.Soybean$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Fertilizer Applications by Country Sum 
soybean.pivot<-dcast(fert, COUNTRY ~., value.var="soybean_Fe", fun.aggregate=sum)
soybean.pivot<-arrange(soybean.pivot,.)
tail(soybean.pivot, 3) #ID Top 3 Fertilizer Countries
sum(soybean.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$soybean_Fe)/sum(soybean.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"
summary(china)

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$soybean_Fe)/sum(soybean.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

brazil<-filter(fert, COUNTRY == "Brazil")
brazil.Perc<-sum(brazil$soybean_Fe)/sum(cassava.pivot$.) 
print(brazil.Perc) # % of Global Fertilizer For Crop 
brazil.color<-"#E28D15"

sum(china.Perc,usa.Perc,brazil.Perc)

#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(soybean_Fe))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$soybean_Fe)), color=china.color, alpha=0.1) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$soybean_Fe)), color=usa.color, alpha=0.2) + 
  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$soybean_Fe)), color=brazil.color, alpha=0.3) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Soybean.fert.sper)[2], slope=coef(Soybean.fert.sper)[3]), color="red", size=2)

#Sunflower 
fert.Sunflower<-subset(fishnet, fishnet$sunflower_ > 0 ) #yield > 0 
xy <- coordinates(fert.Sunflower)
xy1<-data.frame(xy)
fert <- fert.Sunflower@data
SunflowerGCO<-readOGR(dsn = "GIS/",layer="Sunflower_GC")
SunflowerGCOcentroid<-gCentroid(SunflowerGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,SunflowerGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
Sunflower.fert<-fert %>% dplyr::select(sunflower_,rescale_ND,COUNTRY)
Sunflower.fert<-Sunflower.fert %>% filter(sunflower_ > 0, na.rm=TRUE)
Sunflower.fert<-Sunflower.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
Sunflower.fert$logSunflowerfert<-log10(Sunflower.fert$sunflower_)


#Links
nb <- poly2nb(fert.Sunflower)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.Sunflower$logSunflowerFert<-log10(fert.Sunflower$sunflower_)
#Spatial Error Model
Sunflower.fert.sper <- errorsarlm(logSunflowerFert ~ Sunflower.fert$rescale_ND, data=fert.Sunflower, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(Sunflower.fert.sper) #Summary
fert.Sunflower$residualsSpecError <- residuals(Sunflower.fert.sper) #Residuals
moran.mc(fert.Sunflower$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.Sunflower$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Fertilizer Applications by Country Sum 
sunflower.pivot<-dcast(fert, COUNTRY ~., value.var="sunflower_", fun.aggregate=sum)
sunflower.pivot<-arrange(sunflower.pivot,.)
tail(sunflower.pivot, 3) #ID Top 3 Fertilizer Countries
sum(sunflower.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$sunflower_)/sum(sunflower.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"

france<-filter(fert, COUNTRY == "France")
france.Perc<-sum(usa$sunflower_)/sum(sunflower.pivot$.) 
print(france.Perc) # % of Global Fertilizer For Crop 
france.color<-"#CE465B"

india<-filter(fert, COUNTRY == "India")
india.Perc<-sum(india$sunflower_)/sum(rice.pivot$.) 
print(india.Perc) # % of Global Fertilizer For Crop 
india.color<-"#078D40"

sum(china.Perc,france.Perc,india.Perc)

#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(sunflower_))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$sunflower_)), color=china.color, alpha=0.1) + 
  geom_point(data=france, aes(x=rescale_ND, y=log10(france$sunflower_)), color=france.color, alpha=0.2) + 
  geom_point(data=india, aes(x=rescale_ND, y=log10(india$sunflower_)), color=india.color, alpha=0.3) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Sunflower.fert.sper)[2], slope=coef(Sunflower.fert.sper)[3]), color="red", size=2)

#Wheat 
fert.Wheat<-subset(fishnet, fishnet$wheat_Fert > 0 ) #yield > 0 
xy <- coordinates(fert.Wheat)
xy1<-data.frame(xy)
fert <- fert.Wheat@data
WheatGCO<-readOGR(dsn = "GIS/",layer="BarleyWheat_GC")
WheatGCOcentroid<-gCentroid(WheatGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,WheatGCOcentroid, longlat = FALSE)
fert$rescale_ND<-distances/(1000*1000)
Wheat.fert<-fert %>% dplyr::select(wheat_Fert,rescale_ND,COUNTRY)
Wheat.fert<-Wheat.fert %>% filter(wheat_Fert > 0, na.rm=TRUE)
Wheat.fert<-Wheat.fert %>% filter(rescale_ND > 0, na.rm=TRUE)
Wheat.fert$logWheatfert<-log10(Wheat.fert$wheat_Fert)


#Links
nb <- poly2nb(fert.Wheat)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fert.Wheat$logWheatFert<-log10(fert.Wheat$wheat_Fert)
#Spatial Error Model
Wheat.fert.sper <- errorsarlm(logWheatFert ~ Wheat.fert$rescale_ND, data=fert.Wheat, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(Wheat.fert.sper) #Summary
fert.Wheat$residualsSpecError <- residuals(Wheat.fert.sper) #Residuals
moran.mc(fert.Wheat$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

ggplot(xy1, aes(x=X1,y=X2, color=fert.Wheat$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals

#Pivot Table To ID Top Fertilizer Applications by Country Sum 
wheat.pivot<-dcast(fert, COUNTRY ~., value.var="wheat_Fert", fun.aggregate=sum)
wheat.pivot<-arrange(wheat.pivot,.)
tail(wheat.pivot, 3) #ID Top 3 Fertilizer Countries
sum(wheat.pivot$.) #Total Fertilizer

#Filter Countries
china<-filter(fert, COUNTRY == "China") 
china.Perc<-sum(china$wheat_Fert)/sum(wheat.pivot$.) 
print(china.Perc) # % of Global Fertilizer For Crop 
china.color<-"#1665AF"

usa<-filter(fert, COUNTRY == "United States")
usa.Perc<-sum(usa$wheat_Fert)/sum(wheat.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

india<-filter(fert, COUNTRY == "India")
india.Perc<-sum(india$wheat_Fert)/sum(wheat.pivot$.) 
print(india.Perc) # % of Global Fertilizer For Crop 
india.color<-"#078D40"



#Plot 
ggplot(fert, aes(x=rescale_ND,y=log10(wheat_Fert))) + 
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$wheat_Fert)), color=china.color, alpha=0.1) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$wheat_Fert)), color=usa.color, alpha=0.2) + 
  geom_point(data=india, aes(x=rescale_ND, y=log10(india$wheat_Fert)), color=india.color, alpha=0.3) + 
  ylab("Log10 Fertilizer") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Wheat.fert.sper)[2], slope=coef(Wheat.fert.sper)[3]), color="red", size=2)



?errorsarlm
