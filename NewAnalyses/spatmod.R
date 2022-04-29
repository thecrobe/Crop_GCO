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
library(Hmisc)

### Graphics

theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank())

fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
fishnet
#Barley 
fishnet.barley<-subset(fishnet, fishnet$mean_barle > 0 ) #yield > 0 
xy <- coordinates(fishnet.barley)
xy1<-data.frame(xy)
yield <- fishnet.barley@data
barleyGCO<-readOGR(dsn = "GIS/",layer="BarleyWheat_GC")
barleyGCOcentroid<-gCentroid(barleyGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,barleyGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
barley.near<-yield %>% dplyr::select(mean_barle,rescale_ND,COUNTRY)
barley.near<-barley.near %>% filter(mean_barle > 0, na.rm=TRUE)
barley.near<-barley.near %>% filter(rescale_ND > 0, na.rm=TRUE)
barley.near$logBarley<-log10(barley.near$mean_barle)


#Links
nb <- poly2nb(fishnet.barley)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.barley$logBarley<-log10(fishnet.barley$mean_barle)
#Spatial Error Model
barley.sper <- errorsarlm(logBarley ~ barley.near$rescale_ND, data=fishnet.barley, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(barley.sper) #Summary
fishnet.barley$residualsSpecError <- residuals(barley.sper) #Residuals
moran.mc(fishnet.barley$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.barley$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
barley.pivot<-dcast(barley.near, COUNTRY ~., value.var="mean_barle", fun.aggregate=sum)
barley.pivot<-arrange(barley.pivot,.)
tail(barley.pivot, 3) #ID Top 3 Producing Countries
barley.top<-tail(barley.pivot, 3)
sum(barley.pivot$.) #Total Production
sum(barley.top$.)/sum(barley.pivot$.) # % produced by top 3 


#Filter Countries
china<-filter(barley.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_barle)/sum(barley.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"
summary(china)

usa<-filter(barley.near, COUNTRY == "United States")
usa.Perc<-sum(usa$mean_barle)/sum(barley.pivot$.) 
print(usa.Perc) # % of Global Production For Crop 
usa.color<-"#8F4A33"

russia<-filter(barley.near, COUNTRY == "Russian Federation")
russia.Perc<-sum(russia$mean_barle)/sum(barley.pivot$.) 
print(russia.Perc) # % of Global Production For Crop 
russia.color<-"#D3867A"

#Plot
ggplot(barley.near, aes(x=rescale_ND, y=barley.near$logBarley)) +
  geom_point(alpha=0.3) +theme_justin + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_barle)), color=china.color, alpha=0.9) + 
#  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_barle)), color=usa.color, alpha=0.5) + 
#  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_barle)),color=russia.color, alpha=0.2) +
  geom_abline(aes(intercept=coef(barley.sper)[2], slope=coef(barley.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Cassava 
fishnet.cassava<-subset(fishnet, fishnet$mean_cassa > 0 ) #yield > 0 
xy <- coordinates(fishnet.cassava)
xy1<-data.frame(xy) 
yield <- fishnet.cassava@data
cassavaGCO<-readOGR(dsn = "GIS/",layer="Cassava_GC")
cassavaGCOcentroid<-gCentroid(cassavaGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,cassavaGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
cassava.near<-yield %>% dplyr::select(mean_cassa,rescale_ND,COUNTRY)
cassava.near<-cassava.near %>% filter(mean_cassa > 0, na.rm=TRUE)
cassava.near<-cassava.near %>% filter(rescale_ND > 0, na.rm=TRUE)
cassava.near$logCassava<-log10(cassava.near$mean_cassa)


#Links
nb <- poly2nb(fishnet.cassava)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.cassava$logCassava<-log10(fishnet.cassava$mean_cassa)
#Spatial Error Model
cassava.sper <- errorsarlm(logCassava ~ cassava.near$rescale_ND, data=fishnet.cassava, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(cassava.sper) #Summary
fishnet.cassava$residualsSpecError <- residuals(cassava.sper) #Residuals
moran.mc(fishnet.cassava$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.cassava$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals

#Pivot Table To ID Top Producers by Country Sum 
cassava.pivot<-dcast(cassava.near, COUNTRY ~., value.var="mean_cassa", fun.aggregate=sum)
cassava.pivot<-arrange(cassava.pivot,.)
tail(cassava.pivot, 3) #ID Top 3 Producing Countries
cassava.top<-tail(cassava.pivot, 3)
sum(cassava.pivot$.) #Total Production
sum(cassava.top$.)/sum(cassava.pivot$.) # % produced by top 3 


#Filter Countries
china<-filter(cassava.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_cassa)/sum(cassava.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"
summary(china)

#Brazil
brazil<-filter(cassava.near, COUNTRY == "Brazil") 
brazil.Perc<-sum(brazil$mean_cassa)/sum(cassava.pivot$.) 
print(brazil.Perc) # % of Global Production For Crop 
brazil.color<-"#E28D15"
summary(brazil)

#Indonesia
indo<-filter(cassava.near, COUNTRY == "Indonesia") 
indo.Perc<-sum(indo$mean_cassa)/sum(cassava.pivot$.) 
print(indo.Perc) # % of Global Production For Crop 
indo.color<-"#63377C"
summary(indo)

#Plot
ggplot(cassava.near, aes(x=rescale_ND, y=cassava.near$logCassava)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_cassa)), color=china.color, alpha=0.9) + 
  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_cassa)), color=brazil.color, alpha=0.5) + 
  geom_point(data=indo, aes(x=rescale_ND, y=log10(indo$mean_cassa)),color=indo.color, alpha=0.2) +
  geom_abline(aes(intercept=coef(cassava.sper)[2], slope=coef(cassava.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Groundnut 
fishnet.groundnut<-subset(fishnet, fishnet$mean_groun > 0 ) #yield > 0 
xy <- coordinates(fishnet.groundnut)
xy1<-data.frame(xy) 
yield <- fishnet.groundnut@data
groundnutGCO<-readOGR(dsn = "GIS/",layer="Groundnut_GC")
groundnutGCOcentroid<-gCentroid(groundnutGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,groundnutGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
groundnut.near<-yield %>% dplyr::select(mean_groun,rescale_ND,COUNTRY)
groundnut.near<-groundnut.near %>% filter(mean_groun > 0, na.rm=TRUE)
groundnut.near<-groundnut.near %>% filter(rescale_ND > 0, na.rm=TRUE)
groundnut.near$logGroundnut<-log10(groundnut.near$mean_groun)


#Links
nb <- poly2nb(fishnet.groundnut)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.groundnut$logGroundnut<-log10(fishnet.groundnut$mean_groun)
#Spatial Error Model
groundnut.sper <- errorsarlm(logGroundnut ~ groundnut.near$rescale_ND, data=fishnet.groundnut, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(groundnut.sper) #Summary
fishnet.groundnut$residualsSpecError <- residuals(groundnut.sper) #Residuals
moran.mc(fishnet.groundnut$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.groundnut$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
groundnut.pivot<-dcast(groundnut.near, COUNTRY ~., value.var="mean_groun", fun.aggregate=sum)
groundnut.pivot<-arrange(groundnut.pivot,.)
tail(groundnut.pivot, 3) #ID Top 3 Producing Countries
groundnut.top<-tail(groundnut.pivot, 3)
sum(groundnut.pivot$.) #Total Production
sum(groundnut.top$.)/sum(groundnut.pivot$.) # % produced by top 3 

#Filter Countries
china<-filter(groundnut.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_groun)/sum(groundnut.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"
summary(china)

indo<-filter(groundnut.near, COUNTRY == "Indonesia")
indo.Perc<-sum(indo$mean_groun)/sum(groundnut.pivot$.) 
print(indo.Perc) # % of Global Production For Crop 
indo.color<-"#63377C"

india<-filter(groundnut.near, COUNTRY == "India")
india.Perc<-sum(india$mean_groun)/sum(groundnut.pivot$.) 
print(india.Perc) # % of Global Production For Crop 
india.color<-"#078D40"

#Plot
ggplot(groundnut.near, aes(x=rescale_ND, y=logGroundnut)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_groun)), color=china.color, alpha=0.5) + 
  geom_point(data=indo, aes(x=rescale_ND, y=log10(indo$mean_groun)),color=indo.color, alpha=0.8) +
  geom_point(data=india, aes(x=rescale_ND, y=log10(india$mean_groun)),color=india.color, alpha=0.5) + 
  geom_abline(aes(intercept=coef(groundnut.sper)[2], slope=coef(groundnut.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Maize 
fishnet.maize<-subset(fishnet, fishnet$mean_maize > 0 ) #yield > 0 
xy <- coordinates(fishnet.maize)
xy1<-data.frame(xy) 
yield <- fishnet.maize@data
maizeGCO<-readOGR(dsn = "GIS/",layer="Maize_GC")

maizeGCOcentroid<-gCentroid(maizeGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,maizeGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
maize.near<-yield %>% dplyr::select(mean_maize,rescale_ND,COUNTRY)
maize.near<-maize.near %>% filter(mean_maize > 0, na.rm=TRUE)
maize.near<-maize.near %>% filter(rescale_ND > 0, na.rm=TRUE)
maize.near$logMaize<-log10(maize.near$mean_maize)


#Links
nb <- poly2nb(fishnet.maize)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.maize$logMaize<-log10(fishnet.maize$mean_maize)
#Spatial Error Model
maize.sper <- errorsarlm(logMaize ~ maize.near$rescale_ND, data=fishnet.maize, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(maize.sper) #Summary
fishnet.maize$residualsSpecError <- residuals(maize.sper) #Residuals
moran.mc(fishnet.maize$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.maize$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
maize.pivot<-dcast(maize.near, COUNTRY ~., value.var="mean_maize", fun.aggregate=sum)
maize.pivot<-arrange(maize.pivot,.)
tail(maize.pivot, 3) #ID Top 3 Producing Countries
maize.top<-tail(maize.pivot, 3)
sum(maize.pivot$.) #Total Production
sum(maize.top$.)/sum(maize.pivot$.) # % produced by top 3 

#Filter Countries
china<-filter(maize.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_maize)/sum(maize.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"
summary(china)

usa<-filter(maize.near, COUNTRY == "United States")
usa.Perc<-sum(usa$mean_maize)/sum(maize.pivot$.) 
print(usa.Perc) # % of Global Production For Crop 
usa.color<-"#8F4A33"

russia<-filter(maize.near, COUNTRY == "Russian Federation")
russia.Perc<-sum(russia$mean_maize)/sum(maize.pivot$.) 
print(russia.Perc) # % of Global Production For Crop 
russia.color<-"#D3867A"


#Plot
ggplot(maize.near, aes(x=rescale_ND, y=logMaize)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_maize)), color=usa.color, alpha=0.8) + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_maize)),color=china.color, alpha=0.5) +
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_maize)),color=russia.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(maize.sper)[2], slope=coef(maize.sper)[2]), color="grey", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")



#Potato 
fishnet.potato<-subset(fishnet, fishnet$mean_potat > 0 ) #yield > 0 
xy <- coordinates(fishnet.potato)
xy1<-data.frame(xy) 
yield <- fishnet.potato@data
potatoGCO<-readOGR(dsn = "GIS/",layer="Potato")

potatoGCOcentroid<-gCentroid(potatoGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,potatoGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
potat.near<-yield %>% dplyr::select(mean_potat,rescale_ND,COUNTRY)
potat.near<-potat.near %>% filter(mean_potat > 0, na.rm=TRUE)
potat.near<-potat.near %>% filter(rescale_ND > 0, na.rm=TRUE)
potat.near$logPotato<-log10(potat.near$mean_potat)


#Links
nb <- poly2nb(fishnet.potato)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.potato$logPotato<-log10(fishnet.potato$mean_potat)
#Spatial Error Model
potato.sper <- errorsarlm(logPotato ~ potato.near$rescale_ND, data=fishnet.potato, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(potato.sper) #Summary
fishnet.potato$residualsSpecError <- residuals(potato.sper) #Residuals
moran.mc(fishnet.potato$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.potato$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Filter Countries

#Pivot Table To ID Top Producers by Country Sum 
potat.pivot<-dcast(potat.near, COUNTRY ~., value.var="mean_potat", fun.aggregate=sum)
potat.pivot<-arrange(potat.pivot,.)
tail(potat.pivot, 3) #ID Top 3 Producing Countries
potat.top<-tail(potat.pivot, 3)
sum(potat.pivot$.) #Total Production
sum(potat.top$.)/sum(potat.pivot$.) # % produced by top 3 

usa<-filter(potat.near, COUNTRY == "United States")
usa.Perc<-sum(usa$mean_potat)/sum(potat.pivot$.) 
print(usa.Perc) # % of Global Production For Crop 
usa.color<-"#8F4A33"

china<-filter(potat.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_potat)/sum(potat.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"

canada<-filter(potat.near, COUNTRY == "Canada") 
canada.Perc<-sum(canada$mean_potat)/sum(potat.pivot$.) 
print(canada.Perc) # % of Global Production For Crop 
canada.color<-"#98C24B"


#Plot
ggplot(potat.near, aes(x=rescale_ND, y=logPotato)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_potat)),color=china.color, alpha=0.8) +
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_potat)),color=usa.color, alpha=0.4) + 
  geom_point(data=canada, aes(x=rescale_ND, y=log10(canada$mean_potat)), color=canada.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(potato.sper)[2], slope=coef(potato.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Rapeseed 
fishnet.rapeseed<-subset(fishnet, fishnet$mean_rapes > 0 ) #yield > 0 
xy <- coordinates(fishnet.rapeseed)
xy1<-data.frame(xy) 
yield <- fishnet.rapeseed@data
rapeseedGCO<-readOGR(dsn = "GIS/",layer="Rapeseed")

rapeseedGCOcentroid<-gCentroid(rapeseedGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,rapeseedGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
rapes.near<-yield %>% dplyr::select(mean_rapes,rescale_ND,COUNTRY)
rapes.near<-rapes.near %>% filter(mean_rapes > 0, na.rm=TRUE)
rapes.near<-rapes.near %>% filter(rescale_ND > 0, na.rm=TRUE)
rapes.near$logRapeseed<-log10(rapes.near$mean_rapes)


#Links
nb <- poly2nb(fishnet.rapeseed)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.rapeseed$logRapeseed<-log10(fishnet.rapeseed$mean_rapes)
#Spatial Error Model
rapeseed.sper <- errorsarlm(logRapeseed ~ rapeseed.near$rescale_ND, data=fishnet.rapeseed, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rapeseed.sper) #Summary
fishnet.rapeseed$residualsSpecError <- residuals(rapeseed.sper) #Residuals
moran.mc(fishnet.rapeseed$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.rapeseed$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
rapes.pivot<-dcast(rapes.near, COUNTRY ~., value.var="mean_rapes", fun.aggregate=sum)
rapes.pivot<-arrange(rapes.pivot,.)
tail(rapes.pivot, 3) #ID Top 3 Producing Countries
rapes.top<-tail(rapes.pivot, 3)
sum(rapes.pivot$.) #Total Production
sum(rapes.top$.)/sum(rapes.pivot$.) # % produced by top 3 

#Filter Countries
china<-filter(rapes.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_rapes)/sum(rapes.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"

brazil<-filter(rapes.near, COUNTRY == "Brazil")
brazil.Perc<-sum(brazil$mean_rapes)/sum(rapes.pivot$.) 
print(brazil.Perc) # % of Global Production For Crop 
brazil.color<-"#E28D15"

russia<-filter(rapes.near, COUNTRY == "Russian Federation")
russia.Perc<-sum(russia$mean_rapes)/sum(rapes.pivot$.) 
print(russia.Perc) # % of Global Production For Crop 
russia.color<-"#D3867A"


#Plot
ggplot(rapes.near, aes(x=rescale_ND, y=logRapeseed)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_rapes)),color=china.color, alpha=0.5) +
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_rapes)),color=russia.color, alpha=0.3) + 
  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_rapes)), color=brazil.color, alpha=0.8) + 
  geom_abline(aes(intercept=coef(rapeseed.sper)[2], slope=coef(rapeseed.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Rice 
fishnet.rice<-subset(fishnet, fishnet$mean_rice_ > 0 ) #yield > 0 
xy <- coordinates(fishnet.rice)
xy1<-data.frame(xy) 
yield <- fishnet.rice@data
riceGCO<-readOGR(dsn = "GIS/",layer="Rice_GC1")

riceGCOcentroid<-gCentroid(riceGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,riceGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
Rice.near<-yield %>% dplyr::select(mean_rice_,rescale_ND,COUNTRY)
Rice.near<-Rice.near %>% filter(mean_rice_ > 0, na.rm=TRUE)
Rice.near<-Rice.near %>% filter(rescale_ND > 0, na.rm=TRUE)
Rice.near$logRice<-log10(Rice.near$mean_rice_)


#Links
nb <- poly2nb(fishnet.rice)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.rice$logRice<-log10(fishnet.rice$mean_rice_)
#Spatial Error Model
rice.sper <- errorsarlm(logRice ~ rice.near$rescale_ND, data=fishnet.rice, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rice.sper) #Summary
fishnet.rice$residualsSpecError <- residuals(rice.sper) #Residuals
moran.mc(fishnet.rice$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.rice$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
Rice.pivot<-dcast(Rice.near, COUNTRY ~., value.var="mean_rice_", fun.aggregate=sum)
Rice.pivot<-arrange(Rice.pivot,.)
tail(Rice.pivot, 3) #ID Top 3 Producing Countries
Rice.top<-tail(Rice.pivot, 3)
sum(Rice.pivot$.) #Total Production
sum(Rice.top$.)/sum(Rice.pivot$.) # % produced by top 3 

#Filter Countries
china<-filter(Rice.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_Rice)/sum(Rice.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"

brazil<-filter(Rice.near, COUNTRY == "Brazil")
brazil.Perc<-sum(brazil$mean_Rice)/sum(Rice.pivot$.) 
print(brazil.Perc) # % of Global Production For Crop 
brazil.color<-"#E28D15"

india<-filter(Rice.near, COUNTRY == "India")
india.Perc<-sum(india$mean_rice_)/sum(Rice.pivot$.) 
print(india.Perc) # % of Global Production For Crop 
india.color<-"#078D40"


#Plot
ggplot(Rice.near, aes(x=rescale_ND, y=logRice)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_rice_)),color=china.color, alpha=0.5) +
  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_rice_)), color=brazil.color, alpha=0.8) + 
  geom_point(data=india, aes(x=rescale_ND, y=log10(india$mean_rice_)),color=russia.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(rice.sper)[2], slope=coef(rice.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Rye 
fishnet.rye<-subset(fishnet, fishnet$mean_rye_Y > 0 ) #yield > 0 
xy <- coordinates(fishnet.rye)
xy1<-data.frame(xy) 
yield <- fishnet.rye@data
ryeGCO<-readOGR(dsn = "GIS/",layer="Rye")

ryeGCOcentroid<-gCentroid(ryeGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,ryeGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
Rye.near<-yield %>% dplyr::select(mean_rye_Y,rescale_ND,COUNTRY)
Rye.near<-Rye.near %>% filter(mean_rye_Y > 0, na.rm=TRUE)
Rye.near<-Rye.near %>% filter(rescale_ND > 0, na.rm=TRUE)
Rye.near$logRye<-log10(rye.near$mean_rye_Y)


#Links
nb <- poly2nb(fishnet.rye)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.rye$logRye<-log10(fishnet.rye$mean_rye_Y)
#Spatial Error Model
rye.sper <- errorsarlm(logRye ~ rye.near$rescale_ND, data=fishnet.rye, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rye.sper) #Summary
fishnet.rye$residualsSpecError <- residuals(rye.sper) #Residuals
moran.mc(fishnet.rye$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.rye$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
rye.pivot<-dcast(Rye.near, COUNTRY ~., value.var="mean_rye_Y", fun.aggregate=sum)
rye.pivot<-arrange(rye.pivot,.)
tail(rye.pivot, 3) #ID Top 3 Producing Countries
rye.top<-tail(rye.pivot, 3)
sum(rye.pivot$.) #Total Production
sum(rye.top$.)/sum(rye.pivot$.) # % produced by top 3 

#Filter Countries

russia<-filter(Rye.near, COUNTRY == "Russian Federation")
russia.Perc<-sum(russia$mean_rye_Y)/sum(rye.pivot$.) 
print(russia.Perc) # % of Global Production For Crop 
russia.color<-"#D3867A"

china<-filter(Rye.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_rye)/sum(rye.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"

usa<-filter(Rye.near, COUNTRY == "United States")
usa.Perc<-sum(usa$mean_rye_Y)/sum(rye.pivot$.) 
print(usa.Perc) # % of Global Fertilizer For Crop 
usa.color<-"#8F4A33"

#Plot
ggplot(Rye.near, aes(x=rescale_ND, y=logRye)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_rye_Y)), color=russia.color, alpha=0.8) + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_rye_Y)),color=china.color, alpha=0.8) +
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_rye_Y)),color=usa.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(rye.sper)[2], slope=coef(rye.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Sorghum 
fishnet.sorghum<-subset(fishnet, fishnet$mean_sorgh> 0 ) #yield > 0 
xy <- coordinates(fishnet.sorghum)
xy1<-data.frame(xy) 
yield <- fishnet.sorghum@data
sorghumGCO<-readOGR(dsn = "GIS/",layer="Sorghum_GC")
sorghumGCOcentroid<-gCentroid(sorghumGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,sorghumGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
sorghum.near<-yield %>% dplyr::select(mean_sorgh,rescale_ND,COUNTRY)
sorghum.near<-sorghum.near %>% filter(mean_sorgh > 0, na.rm=TRUE)
sorghum.near<-sorghum.near %>% filter(rescale_ND > 0, na.rm=TRUE)
sorghum.near$logSorghum<-log10(sorghum.near$mean_sorgh)


#Links
nb <- poly2nb(fishnet.sorghum)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.sorghum$logSorghum<-log10(fishnet.sorghum$mean_sorgh)
#Spatial Error Model
sorghum.sper <- errorsarlm(logSorghum ~ sorghum.near$rescale_ND, data=fishnet.sorghum, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(sorghum.sper) #Summary
fishnet.sorghum$residualsSpecError <- residuals(sorghum.sper) #Residuals
moran.mc(fishnet.sorghum$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.sorghum$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
sorghum.pivot<-dcast(sorghum.near, COUNTRY ~., value.var="mean_sorgh", fun.aggregate=sum)
sorghum.pivot<-arrange(sorghum.pivot,.)
tail(sorghum.pivot, 3) #ID Top 3 Producing Countries
sorghum.top<-tail(sorghum.pivot, 3)
sum(sorghum.pivot$.) #Total Production
sum(sorghum.top$.)/sum(sorghum.pivot$.) # % produced by top 3 

#Filter Countries

china<-filter(sorghum.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_sorgh)/sum(sorghum.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"

usa<-filter(sorghum.near, COUNTRY == "United States")
usa.Perc<-sum(usa$mean_sorgh)/sum(sorghum.pivot$.) 
print(usa.Perc) # % of Global Production For Crop 
usa.color<-"#8F4A33"

russia<-filter(sorghum.near, COUNTRY == "Russian Federation")
russia.Perc<-sum(russia$mean_sorgh)/sum(sorghum.pivot$.) 
print(russia.Perc) # % of Global Production For Crop 
russia.color<-"#D3867A"


#Plot
ggplot(sorghum.near, aes(x=rescale_ND, y=logSorghum)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_sorgh)),color=china.color, alpha=0.8) +
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_sorgh)),color=usa.color, alpha=0.5) + 
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_sorgh)), color=russia.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(sorghum.sper)[2], slope=coef(sorghum.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Soybean 
fishnet.soybean<-subset(fishnet, fishnet$mean_soybe > 0 ) #yield > 0 
xy <- coordinates(fishnet.soybean)
xy1<-data.frame(xy) 
yield <- fishnet.soybean@data
soybeanGCO<-readOGR(dsn = "GIS/",layer="Soy_GC")
soybeanGCOcentroid<-gCentroid(soybeanGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,soybeanGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
soybean.near<-yield %>% dplyr::select(mean_soybe,rescale_ND,COUNTRY)
soybean.near<-soybean.near %>% filter(mean_soybe > 0, na.rm=TRUE)
soybean.near<-soybean.near %>% filter(rescale_ND > 0, na.rm=TRUE)
soybean.near$logSoybean<-log10(soybean.near$mean_soybe)


#Links
nb <- poly2nb(fishnet.soybean)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.soybean$logSoybean<-log10(fishnet.soybean$mean_soybe)
#Spatial Error Model
soybean.sper <- errorsarlm(logSoybean ~ soybean.near$rescale_ND, data=fishnet.soybean, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(soybean.sper) #Summary
fishnet.soybean$residualsSpecError <- residuals(soybean.sper) #Residuals
moran.mc(fishnet.soybean$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.soybean$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
soybe.pivot<-dcast(soybean.near, COUNTRY ~., value.var="mean_soybe", fun.aggregate=sum)
soybe.pivot<-arrange(soybe.pivot,.)
tail(soybe.pivot, 3) #ID Top 3 Producing Countries
soybe.top<-tail(soybe.pivot, 3)
sum(soybe.pivot$.) #Total Production
sum(soybe.top$.)/sum(soybe.pivot$.) # % produced by top 3 

#Filter Countries
china<-filter(soybean.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_soybe)/sum(soybe.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"

usa<-filter(soybean.near, COUNTRY == "United States")
usa.Perc<-sum(usa$mean_soybe)/sum(soybe.pivot$.) 
print(usa.Perc) # % of Global Production For Crop 
usa.color<-"#8F4A33"

russia<-filter(soybean.near, COUNTRY == "Russian Federation")
russia.Perc<-sum(russia$mean_soybe)/sum(soybe.pivot$.) 
print(russia.Perc) # % of Global Production For Crop 
russia.color<-"#D3867A"


#Plot
ggplot(soybean.near, aes(x=rescale_ND, y=logSoybean)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_soybe)),color=china.color, alpha=0.8) +
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_soybe)),color=usa.color, alpha=0.3) + 
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_soybe)), color=russia.color, alpha=0.1) + 
  geom_abline(aes(intercept=coef(soybean.sper)[2], slope=coef(soybean.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Sunflower 
fishnet.sunflower<-subset(fishnet, fishnet$mean_sunfl > 0 ) #yield > 0 
xy <- coordinates(fishnet.sunflower)
xy1<-data.frame(xy) 
yield <- fishnet.sunflower@data
sunflowerGCO<-readOGR(dsn = "GIS/",layer="Sunflower_GC")
sunflowerGCOcentroid<-gCentroid(sunflowerGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,sunflowerGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
sunflower.near<-yield %>% dplyr::select(mean_sunfl,rescale_ND,COUNTRY)
sunflower.near<-sunflower.near %>% filter(mean_sunfl > 0, na.rm=TRUE)
sunflower.near<-sunflower.near %>% filter(rescale_ND > 0, na.rm=TRUE)
sunflower.near$logSunflower<-log10(sunflower.near$mean_sunfl)


#Links
nb <- poly2nb(fishnet.sunflower)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.sunflower$logSunflower<-log10(fishnet.sunflower$mean_sunfl)
#Spatial Error Model
sunflower.sper <- errorsarlm(logSunflower ~ sunflower.near$rescale_ND, data=fishnet.sunflower, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(sunflower.sper) #Summary
fishnet.sunflower$residualsSpecError <- residuals(sunflower.sper) #Residuals
moran.mc(fishnet.sunflower$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.sunflower$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
sunflower.pivot<-dcast(sunflower.near, COUNTRY ~., value.var="mean_sunfl", fun.aggregate=sum)
sunflower.pivot<-arrange(sunflower.pivot,.)
tail(sunflower.pivot, 3) #ID Top 3 Producing Countries
sunflower.top<-tail(sunflower.pivot, 3)
sum(sunflower.pivot$.) #Total Production
sum(sunflower.top$.)/sum(sunflower.pivot$.) # % produced by top 3 

#Filter Countries
usa<-filter(sunflower.near, COUNTRY == "United States")
usa.Perc<-sum(usa$mean_sunfl)/sum(sunflower.pivot$.) 
print(usa.Perc) # % of Global Production For Crop 
usa.color<-"#8F4A33"

china<-filter(sunflower.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_sunfl)/sum(sunflower.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"

brazil<-filter(sunflower.near, COUNTRY == "Brazil")
brazil.Perc<-sum(brazil$mean_sunfle)/sum(cassava.pivot$.) 
print(brazil.Perc) # % of Global Production For Crop 
brazil.color<-"#E28D15"

#Plot
ggplot(sunflower.near, aes(x=rescale_ND, y=logSunflower)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_sunfl)),color=china.color, alpha=0.8) +
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_sunfl)), color=usa.color, alpha=0.4) + 
  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_sunfl)), color=brazil.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(sunflower.sper)[2], slope=coef(sunflower.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Wheat 
fishnet.wheat<-subset(fishnet, fishnet$mean_wheat > 0 ) #yield > 0 
xy <- coordinates(fishnet.wheat)
xy1<-data.frame(xy) 
yield <- fishnet.wheat@data
wheatGCO<-readOGR(dsn = "GIS/",layer="BarleyWheat_GC")
wheatGCOcentroid<-gCentroid(wheatGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,wheatGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
wheat.near<-yield %>% dplyr::select(mean_wheat,rescale_ND,COUNTRY)
wheat.near<-wheat.near %>% filter(mean_wheat > 0, na.rm=TRUE)
wheat.near<-wheat.near %>% filter(rescale_ND > 0, na.rm=TRUE)
wheat.near$logWheat<-log10(wheat.near$mean_wheat)

#Links
nb <- poly2nb(fishnet.wheat)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.wheat$logWheat<-log10(fishnet.wheat$mean_wheat)
#Spatial Error Model
wheat.sper <- errorsarlm(logWheat ~ wheat.near$rescale_ND, data=fishnet.wheat, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(wheat.sper) #Summary
fishnet.wheat$residualsSpecError <- residuals(wheat.sper) #Residuals
moran.mc(fishnet.wheat$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.wheat$residualsSpecError, )) + 
  geom_point(size=1.0) #plot residuals


#Pivot Table To ID Top Producers by Country Sum 
wheat.pivot<-dcast(wheat.near, COUNTRY ~., value.var="mean_wheat", fun.aggregate=sum)
wheat.pivot<-arrange(wheat.pivot,.)
tail(wheat.pivot, 3) #ID Top 3 Producing Countries
wheat.top<-tail(wheat.pivot, 3)
sum(wheat.pivot$.) #Total Production
sum(wheat.top$.)/sum(wheat.pivot$.) # % produced by top 3 

#Filter Countries

usa<-filter(wheat.near, COUNTRY == "United States")
usa.Perc<-sum(usa$mean_wheat)/sum(wheat.pivot$.) 
print(usa.Perc) # % of Global Production For Crop 
usa.color<-"#8F4A33"

china<-filter(wheat.near, COUNTRY == "China") 
china.Perc<-sum(china$mean_wheat)/sum(wheat.pivot$.) 
print(china.Perc) # % of Global Production For Crop 
china.color<-"#1665AF"

russia<-filter(wheat.near, COUNTRY == "Russian Federation")
russia.Perc<-sum(russia$mean_wheat)/sum(wheat.pivot$.) 
print(russia.Perc) # % of Global Production For Crop 
russia.color<-"#D3867A"


#Plot
ggplot(wheat.near, aes(x=rescale_ND, y=logWheat)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_wheat)),color=china.color, alpha=0.8) +
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_wheat)),color=usa.color, alpha=0.4) + 
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_wheat)), color=russia.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(wheat.sper)[2], slope=coef(wheat.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")




