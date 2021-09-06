
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
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_barle)), color=china.color, alpha=0.9) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_barle)), color=usa.color, alpha=0.5) + 
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_barle)),color=russia.color, alpha=0.2) +
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
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.groundnut$residualsSpecError, )) + 
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
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_groun)), color        =china.color, alpha=0.5) + 
  geom_point(data=indo, aes(x=rescale_ND, y=log10(indo$mean_groun)),                color=indo.color, alpha=0.8) +
  geom_point(data=india, aes(x=rescale_ND, y=log10(india$mean_groun)),           color=india.color, alpha=0.5) + 
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
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.maize$residualsSpecError, )) + 
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
potat.near<-potato.near %>% filter(mean_potat > 0, na.rm=TRUE)
potat.near<-potato.near %>% filter(rescale_ND > 0, na.rm=TRUE)
potat.near$logPotato<-log10(potato.near$mean_potat)


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
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.potato$residualsSpecError, )) + 
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
fishnet.rapeseed<-subset(fishnet, fishnet$mean_rapeseed > 0 ) #yield > 0 
xy <- coordinates(fishnet.rapeseed)
xy1<-data.frame(xy) 
yield <- fishnet.rapeseed@data
rapeseedGCO<-readOGR(dsn = "GIS/",layer="Rapeseed_GC")

rapeseedGCOcentroid<-gCentroid(rapeseedGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,rapeseedGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
rapeseed.near<-yield %>% dplyr::select(mean_rapes,rescale_ND,COUNTRY)
rapeseed.near<-rapeseed.near %>% filter(mean_rapes > 0, na.rm=TRUE)
rapeseed.near<-rapeseed.near %>% filter(rescale_ND > 0, na.rm=TRUE)
rapeseed.near$logRapeseed<-log10(rapeseed.near$mean_rapes)


#Links
nb <- poly2nb(fishnet.rapeseed)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.rapeseed$logRapeseed<-log10(fishnet.rapeseed$mean_rapeseed)
#Spatial Error Model
rapeseed.sper <- errorsarlm(logRapeseed ~ rapeseed.near$rescale_ND, data=fishnet.rapeseed, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rapeseed.sper) #Summary
fishnet.rapeseed$residualsSpecError <- residuals(rapeseed.sper) #Residuals
moran.mc(fishnet.rapeseed$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.rapeseed$residualsSpecError, )) + 
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
ggplot(rapes.near, aes(x=rescale_ND, y=lograpes)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_rapes)),color=china.color, alpha=0.5) +
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_rapes)),color=russia.color, alpha=0.3) + 
  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_rapes)), color=brazil.color, alpha=0.8) + 
  geom_abline(aes(intercept=coef(rapes.sper)[2], slope=coef(rapes.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


#Rice 
fishnet.rice<-subset(fishnet, fishnet$mean_rice_ > 0 ) #yield > 0 
xy <- coordinates(fishnet.rice)
xy1<-data.frame(xy) 
yield <- fishnet.rice@data
riceGCO<-readOGR(dsn = "GIS/",layer="Rice_GC")

riceGCOcentroid<-gCentroid(riceGCO)

#Distances + Yield Transformation
distances <- spDistsN1(xy,riceGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
rice.near<-yield %>% dplyr::select(mean_rice_,rescale_ND,COUNTRY)
rice.near<-rice.near %>% filter(mean_rice_ > 0, na.rm=TRUE)
rice.near<-rice.near %>% filter(rescale_ND > 0, na.rm=TRUE)
rice.near$logRice<-log10(rice.near$mean_rice_)


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
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.rice$residualsSpecError, )) + 
  geom_point(size=1.0) #plot residuals


#Plot
ggplot(Rice.near, aes(x=rescale_ND, y=logRice)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_rice_)),color=china.color, alpha=0.5) +
  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_rice_)), color=brazil.color, alpha=0.8) + 
  geom_point(data=india, aes(x=rescale_ND, y=log10(india$mean_rice_)),color=russia.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(Rice.power)[1], slope=coef(Rice.power)[2]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")

