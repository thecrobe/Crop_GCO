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
library(gridExtra)

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

#Plot
barley.plot<-ggplot(barley.near, aes(x=rescale_ND, y=barley.near$logBarley)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
  geom_abline(aes(intercept=coef(barley.sper)[2], slope=coef(barley.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title="Barley",subtitle="Slope=0.015, p < 0.0001")

=======
  geom_point(alpha=0.3) +theme_justin + 
  geom_abline(aes(intercept=coef(barley.sper)[2], slope=coef(barley.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")
>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8


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

#Plot
cassava.plot<-ggplot(cassava.near, aes(x=rescale_ND, y=cassava.near$logCassava)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_cassa)), color=china.color, alpha=0.9) + 
#  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_cassa)), color=brazil.color, alpha=0.5) + 
#  geom_point(data=indo, aes(x=rescale_ND, y=log10(indo$mean_cassa)),color=indo.color, alpha=0.2) +
  geom_abline(aes(intercept=coef(cassava.sper)[2], slope=coef(cassava.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") + labs(title="Cassava",subtitle="Slope=0.039, p < 0.0001")
=======
  geom_point(alpha=0.3) +theme_justin +
  geom_abline(aes(intercept=coef(cassava.sper)[2], slope=coef(cassava.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")
>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8


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


#Plot
groundnut.plot<-ggplot(groundnut.near, aes(x=rescale_ND, y=logGroundnut)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_groun)), color=china.color, alpha=0.5) + 
#  geom_point(data=indo, aes(x=rescale_ND, y=log10(indo$mean_groun)),color=indo.color, alpha=0.8) +
#  geom_point(data=india, aes(x=rescale_ND, y=log10(india$mean_groun)),color=india.color, alpha=0.5) + 
  geom_abline(aes(intercept=coef(groundnut.sper)[2], slope=coef(groundnut.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title="Groundnut", subtitle="Slope=0.02, p < 0.0001")
=======
  geom_point(alpha=0.3) +theme_justin + 
  geom_abline(aes(intercept=coef(groundnut.sper)[2], slope=coef(groundnut.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")
>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8


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

#Plot
maize.plot<-ggplot(maize.near, aes(x=rescale_ND, y=logMaize)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_maize)), color=usa.color, alpha=0.8) + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_maize)),color=china.color, alpha=0.5) +
#  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_maize)),color=russia.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(maize.sper)[2], slope=coef(maize.sper)[3]), color="grey", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title="Maize", subtitle="Slope=0.002, p=0.437")
=======
  geom_point(alpha=0.3) +theme_justin +
  geom_abline(aes(intercept=coef(maize.sper)[2], slope=coef(maize.sper)[2]), color="grey", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")


>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8

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
potato.sper <- errorsarlm(logPotato ~ potat.near$rescale_ND, data=fishnet.potato, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(potato.sper) #Summary
fishnet.potato$residualsSpecError <- residuals(potato.sper) #Residuals
moran.mc(fishnet.potato$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.potato$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals

#Plot
potato.plot<-ggplot(potat.near, aes(x=rescale_ND, y=logPotato)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_potat)),color=china.color, alpha=0.8) +
#  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_potat)),color=usa.color, alpha=0.4) + 
#  geom_point(data=canada, aes(x=rescale_ND, y=log10(canada$mean_potat)), color=canada.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(potato.sper)[2], slope=coef(potato.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title="Potato", subtitle="Slope=0.02, p < 0.0001")
=======
  geom_point(alpha=0.3) +theme_justin  + 
  geom_abline(aes(intercept=coef(potato.sper)[2], slope=coef(potato.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")
>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8


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
rapeseed.sper <- errorsarlm(logRapeseed ~ rapes.near$rescale_ND, data=fishnet.rapeseed, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rapeseed.sper) #Summary
fishnet.rapeseed$residualsSpecError <- residuals(rapeseed.sper) #Residuals
moran.mc(fishnet.rapeseed$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.rapeseed$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Plot
rapeseed.plot<-ggplot(rapes.near, aes(x=rescale_ND, y=logRapeseed)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_rapes)),color=china.color, alpha=0.5) +
#  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_rapes)),color=russia.color, alpha=0.3) + 
#  geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_rapes)), color=brazil.color, alpha=0.8) + 
  geom_abline(aes(intercept=coef(rapeseed.sper)[2], slope=coef(rapeseed.sper)[3]), color="grey", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title="Rapeseed", subtitle="Slope=-0.001, p =0.0501")
=======
  geom_point(alpha=0.3) +theme_justin  + 
  geom_abline(aes(intercept=coef(rapeseed.sper)[2], slope=coef(rapeseed.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")
>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8


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
rice.sper <- errorsarlm(logRice ~ Rice.near$rescale_ND, data=fishnet.rice, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rice.sper) #Summary
fishnet.rice$residualsSpecError <- residuals(rice.sper) #Residuals
moran.mc(fishnet.rice$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.rice$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals

#Plot
rice.plot<-ggplot(Rice.near, aes(x=rescale_ND, y=logRice)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
  #geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_rice_)),color=china.color, alpha=0.5) +
  #geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_rice_)), color=brazil.color, alpha=0.8) + 
  #geom_point(data=india, aes(x=rescale_ND, y=log10(india$mean_rice_)),color=russia.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(rice.sper)[2], slope=coef(rice.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title="Rice", subtitle="Slope = -0.017, p < 0.0001")
=======
  geom_point(alpha=0.3) +theme_justin  + 
  geom_abline(aes(intercept=coef(rice.sper)[2], slope=coef(rice.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")
>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8


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
Rye.near$logRye<-log10(Rye.near$mean_rye_Y)


#Links
nb <- poly2nb(fishnet.rye)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
fishnet.rye$logRye<-log10(fishnet.rye$mean_rye_Y)
#Spatial Error Model
rye.sper <- errorsarlm(logRye ~ Rye.near$rescale_ND, data=fishnet.rye, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rye.sper) #Summary
fishnet.rye$residualsSpecError <- residuals(rye.sper) #Residuals
moran.mc(fishnet.rye$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
#plot residuals
ggplot(xy1, aes(x=X1,y=X2, color=fishnet.rye$residualsSpecError)) + 
  geom_point(size=1.0) #plot residuals


#Plot
rye.plot<-ggplot(Rye.near, aes(x=rescale_ND, y=logRye)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_rye_Y)), color=russia.color, alpha=0.8) + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_rye_Y)),color=china.color, alpha=0.8) +
#  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_rye_Y)),color=usa.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(rye.sper)[2], slope=coef(rye.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title="Rye", subtitle = "Slope=-0.017, p=0.0004")
=======
  geom_point(alpha=0.3) +theme_justin  + 
  geom_abline(aes(intercept=coef(rye.sper)[2], slope=coef(rye.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")

>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8

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

#Plot
sorghum.plot<-ggplot(sorghum.near, aes(x=rescale_ND, y=logSorghum)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_sorgh)),color=china.color, alpha=0.8) +
#  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_sorgh)),color=usa.color, alpha=0.5) + 
#  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_sorgh)), color=russia.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(sorghum.sper)[2], slope=coef(sorghum.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title="Sorghum", subtitle="Slope=0.046, p < 0.0001")
=======
  geom_point(alpha=0.3) +theme_justin + 
  geom_abline(aes(intercept=coef(sorghum.sper)[2], slope=coef(sorghum.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")
>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8


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

#Plot
<<<<<<< HEAD
soybean.plot<-ggplot(soybean.near, aes(x=rescale_ND, y=logSoybean)) +
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_soybe)),color=china.color, alpha=0.8) +
#  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_soybe)),color=usa.color, alpha=0.3) + 
#  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_soybe)), color=russia.color, alpha=0.1) + 
  geom_abline(aes(intercept=coef(soybean.sper)[2], slope=coef(soybean.sper)[3]), color="grey", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title="Soybean", subtitle="Slope=0.123, p=0.259")
=======
soybean.plot<-ggplot(soybean.near, aes(x=rescale_ND, y=logSoybean))  + 
  geom_abline(aes(intercept=coef(soybean.sper)[2], slope=coef(soybean.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")
>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8


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


#Plot
sunflower.plot<-ggplot(sunflower.near, aes(x=rescale_ND, y=logSunflower)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_sunfl)),color=china.color, alpha=0.8) +
#  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_sunfl)), color=usa.color, alpha=0.4) + 
# geom_point(data=brazil, aes(x=rescale_ND, y=log10(brazil$mean_sunfl)), color=brazil.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(sunflower.sper)[2], slope=coef(sunflower.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") +labs(title ="Sunflower", subtitle="Slope=-0.007, p=0.0476")
=======
  geom_point(alpha=0.3) +theme_justin + 
  geom_abline(aes(intercept=coef(sunflower.sper)[2], slope=coef(sunflower.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")
>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8


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

#Plot
wheat.plot<-ggplot(wheat.near, aes(x=rescale_ND, y=logWheat)) +
<<<<<<< HEAD
  geom_point(alpha=0.1) +theme_justin + 
#  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_wheat)),color=china.color, alpha=0.8) +
#  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_wheat)),color=usa.color, alpha=0.4) + 
#  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_wheat)), color=russia.color, alpha=0.3) + 
  geom_abline(aes(intercept=coef(wheat.sper)[2], slope=coef(wheat.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield") + labs(title="Wheat", subtitle="Slope=0.133, p=0.0036")

#Arrange final plot
grid.arrange(barley.plot,cassava.plot,groundnut.plot,maize.plot,potato.plot,rapeseed.plot,rice.plot,rye.plot,sorghum.plot,soybean.plot,sunflower.plot,wheat.plot, nrow = 3)
=======
  geom_point(alpha=0.3) +theme_justin  + 
  geom_abline(aes(intercept=coef(wheat.sper)[2], slope=coef(wheat.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")

>>>>>>> 69b2d5123fa5200145e09d95d88506efbd2f53f8



