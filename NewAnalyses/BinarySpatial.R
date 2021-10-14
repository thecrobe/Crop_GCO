
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

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


#colors
sighigh<-"#00ADB6"
siglow<-"#FFA800"
insig<-"#7DA3B2"

##### Models
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
proj4string(fishnet) <- CRS("+init=epsg:3786") #world equidistant cylindrical (sphere)
fish<-fishnet@data
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID
joined<-na.omit(inner_join(fish,mapping, by="Fishnet_ID"))

#Barley 
fishnet$BarleyGCO<-joined$BarleyGCO
mean.barley<-subset(fishnet, fishnet$mean_barle > 0 ) #yield > 0 
mean.barley.data<-mean.barley@data

#Links
nb <- poly2nb(mean.barley)
lw <- nb2listw(nb, zero.policy = TRUE)
#Spatial Error Model
barley.binary.sper <- errorsarlm(log10(mean_barle) ~ BarleyGCO, data=mean.barley.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(barley.binary.sper) #Summary
mean.barley$residualsSpecError <- residuals(barley.binary.sper) #Residuals
moran.barley<-moran.mc(mean.barley$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.barley)

#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
barley<-na.omit(yield.mapping %>% dplyr::select(barley_HgHa,BarleyGCO))

#Inside vs Outside GCO 
barley.inside<-filter(barley, BarleyGCO == "Inside")
barley.outside<-filter(barley, BarleyGCO == "Outside")
barley.m<-melt(barley)
summary(barley.inside)
summary(barley.outside)


#Plot
ggplot(barley.m, aes(x=barley.m$BarleyGCO, y=barley.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 

#Cassava 
fishnet$CassavaGCO<-joined$CassavaGCO
mean.cassava<-subset(fishnet, fishnet$mean_cassa > 0 ) #yield > 0 
mean.cassava.data<-mean.cassava@data

#Links
nb <- poly2nb(mean.cassava)
lw <- nb2listw(nb, zero.policy = TRUE)
#Spatial Error Model
cassava.binary.sper <- errorsarlm(log10(mean_cassa) ~ CassavaGCO, data=mean.cassava.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(cassava.binary.sper) #Summary
mean.cassava$residualsSpecError <- residuals(cassava.binary.sper) #Residuals
moran.cassava<-moran.mc(mean.cassava$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.cassava)

#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
cassava<-na.omit(yield.mapping %>% dplyr::select(cassava_HgHa,CassavaGCO))

#Inside vs Outside GCO 
cassava.inside<-filter(cassava, CassavaGCO == "Inside")
cassava.outside<-filter(cassava, CassavaGCO == "Outside")
cassava.m<-melt(cassava)
summary(cassava.inside)
summary(cassava.outside)


#Plot
ggplot(cassava.m, aes(x=cassava.m$CassavaGCO, y=cassava.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 

#Groundnut 
fishnet$GroundnutGCO<-joined$GroundnutGCO
mean.groundnut<-subset(fishnet, fishnet$mean_groun > 0 ) #yield > 0 
mean.groundnut.data<-mean.groundnut@data

#Links
nb <- poly2nb(mean.groundnut)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
groundnut.binary.sper <- errorsarlm(log10(mean_groun) ~ GroundnutGCO, data=mean.groundnut.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(groundnut.binary.sper) #Summary
mean.groundnut$residualsSpecError <- residuals(groundnut.binary.sper) #Residuals
moran.groundnut<-moran.mc(mean.groundnut$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.groundnut)

#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
groundnut<-na.omit(yield.mapping %>% dplyr::select(groundnut_HgHa,GroundnutGCO))

#Inside vs Outside GCO 
groundnut.inside<-filter(groundnut, GroundnutGCO == "Inside")
groundnut.outside<-filter(groundnut, GroundnutGCO == "Outside")
groundnut.m<-melt(groundnut)
summary(groundnut.inside)
summary(groundnut.outside)


#Plot
ggplot(groundnut.m, aes(x=groundnut.m$GroundnutGCO, y=groundnut.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


#Maize 
fishnet$MaizeGCO<-joined$MaizeGCO
mean.maize<-subset(fishnet, fishnet$mean_maize > 0 ) #yield > 0 
mean.maize.data<-mean.maize@data

#Links
nb <- poly2nb(mean.maize)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
maize.binary.sper <- errorsarlm(log10(mean_maize) ~ MaizeGCO, data=mean.maize.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(maize.binary.sper) #Summary
mean.maize$residualsSpecError <- residuals(maize.binary.sper) #Residuals
moran.maize<-moran.mc(mean.maize$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.maize)

#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
maize<-na.omit(yield.mapping %>% dplyr::select(maize_HgHa,MaizeGCO))

#Inside vs Outside GCO 
maize.inside<-filter(maize, MaizeGCO == "Inside")
maize.outside<-filter(maize, MaizeGCO == "Outside")
maize.m<-melt(maize)
summary(maize.inside)
summary(maize.outside)

#Plot
ggplot(maize.m, aes(x=maize.m$MaizeGCO, y=maize.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 

#Potato 
fishnet$PotatoGCO<-joined$PotatoGCO
mean.potato<-subset(fishnet, fishnet$mean_potat > 0 ) #yield > 0 
mean.potato.data<-mean.potato@data

#Links
nb <- poly2nb(mean.potato)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
potato.binary.sper <- errorsarlm(log10(mean_potat) ~ PotatoGCO, data=mean.potato.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(potato.binary.sper) #Summary
mean.potato$residualsSpecError <- residuals(potato.binary.sper) #Residuals
moran.potato<-moran.mc(mean.potato$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.potato)

#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
potato<-na.omit(yield.mapping %>% dplyr::select(mean_potato,PotatoGCO))

#Inside vs Outside GCO 
potato.inside<-filter(potato, PotatoGCO == "Inside")
potato.outside<-filter(potato, PotatoGCO == "Outside")
potato.m<-melt(potato)
summary(potato.inside)
summary(potato.outside)

#Plot
ggplot(potato.m, aes(x=potato.m$PotatoGCO, y=potato.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


#Rapeseed 
fishnet$RapeseedGCO<-joined$RapeseedGCO
mean.rapeseed<-subset(fishnet, fishnet$mean_rapes > 0 ) #yield > 0 
mean.rapeseed.data<-mean.rapeseed@data

#Links
nb <- poly2nb(mean.rapeseed)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
rapeseed.binary.sper <- errorsarlm(log10(mean_rapes) ~ RapeseedGCO, data=mean.rapeseed.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rapeseed.binary.sper) #Summary
mean.rapeseed$residualsSpecError <- residuals(rapeseed.binary.sper) #Residuals
moran.rapeseed<-moran.mc(mean.rapeseed$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.rapeseed)

#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
rapeseed<-na.omit(yield.mapping %>% dplyr::select(rapeseed_HgHa,RapeseedGCO))

#Inside vs Outside GCO 
rapeseed.inside<-filter(rapeseed, RapeseedGCO == "Inside")
rapeseed.outside<-filter(rapeseed, RapeseedGCO == "Outside")
rapeseed.m<-melt(rapeseed)
summary(rapeseed.inside)
summary(rapeseed.outside)

#Plot
ggplot(rapeseed.m, aes(x=rapeseed.m$RapeseedGCO, y=rapeseed.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO')


#Rice 
fishnet$RiceGCO<-joined$RiceGCO
mean.rice<-subset(fishnet, fishnet$mean_rice_ > 0 ) #yield > 0 
mean.rice.data<-mean.rice@data

#Links
nb <- poly2nb(mean.rice)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
rice.binary.sper <- errorsarlm(log10(mean_rice_) ~ RiceGCO, data=mean.rice.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rice.binary.sper) #Summary
mean.rice$residualsSpecError <- residuals(rice.binary.sper) #Residuals
moran.rice<-moran.mc(mean.rice$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.rice)


#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
rice<-na.omit(yield.mapping %>% dplyr::select(rice_HgHa,RiceGCO))

#Inside vs Outside GCO 
rice.inside<-filter(rice, RiceGCO == "Inside")
rice.outside<-filter(rice, RiceGCO == "Outside")
rice.m<-melt(rice)
summary(rice.inside)
summary(rice.outside)

#Plot
ggplot(rice.m, aes(x=rice.m$RiceGCO, y=rice.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO')


#Rye 
fishnet$RyeGCO<-joined$RyeGCO
mean.rye<-subset(fishnet, fishnet$mean_rye_Y > 0 ) #yield > 0 
mean.rye.data<-mean.rye@data

#Links
nb <- poly2nb(mean.rye)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
rye.binary.sper <- errorsarlm(log10(mean_rye_Y) ~ RyeGCO, data=mean.rye.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(rye.binary.sper) #Summary
mean.rye$residualsSpecError <- residuals(rye.binary.sper) #Residuals
moran.rye<-moran.mc(mean.rye$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.rye)

#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
rye<-na.omit(yield.mapping %>% dplyr::select(rye_HgHa,RyeGCO))

#Inside vs Outside GCO 
rye.inside<-filter(rye, RyeGCO == "Inside")
rye.outside<-filter(rye, RyeGCO == "Outside")
rye.m<-melt(rye)
summary(rye.inside)
summary(rye.outside)

#Plot
ggplot(rye.m, aes(x=rye.m$RyeGCO, y=rye.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO')



#Sorghum 
fishnet$SorghumGCO<-joined$SorghumGCO
mean.sorghum<-subset(fishnet, fishnet$mean_sorgh > 0 ) #yield > 0 
mean.sorghum.data<-mean.sorghum@data

#Links
nb <- poly2nb(mean.sorghum)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
sorghum.binary.sper <- errorsarlm(log10(mean_sorgh) ~ SorghumGCO, data=mean.sorghum.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(sorghum.binary.sper) #Summary
mean.sorghum$residualsSpecError <- residuals(sorghum.binary.sper) #Residuals
moran.sorghum<-moran.mc(mean.sorghum$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.sorghum)

#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
sorghum<-na.omit(yield.mapping %>% dplyr::select(sorghum_HgHa,SorghumGCO))

#Inside vs Outside GCO 
sorghum.inside<-filter(sorghum, SorghumGCO == "Inside")
sorghum.outside<-filter(sorghum, SorghumGCO == "Outside")
sorghum.m<-melt(sorghum)
summary(sorghum.inside)
summary(sorghum.outside)

#Plot
ggplot(sorghum.m, aes(x=sorghum.m$SorghumGCO, y=sorghum.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO')


#Soybean 
fishnet$SoybeanGCO<-joined$SoybeanGCO
mean.soybean<-subset(fishnet, fishnet$mean_soybe > 0 ) #yield > 0 
mean.soybean.data<-mean.soybean@data

#Links
nb <- poly2nb(mean.soybean)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
soybean.binary.sper <- errorsarlm(log10(mean_soybe) ~ SoybeanGCO, data=mean.soybean.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(soybean.binary.sper) #Summary
mean.soybean$residualsSpecError <- residuals(soybean.binary.sper) #Residuals
moran.soybean<-moran.mc(mean.soybean$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.soybean)

#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
soybean<-na.omit(yield.mapping %>% dplyr::select(soybean_HgHa,SoybeanGCO))

#Inside vs Outside GCO 
soybean.inside<-filter(soybean, SoybeanGCO == "Inside")
soybean.outside<-filter(soybean, SoybeanGCO == "Outside")
soybean.m<-melt(soybean)
summary(soybean.inside)
summary(soybean.outside)

#Plot
ggplot(soybean.m, aes(x=soybean.m$SoybeanGCO, y=soybean.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO')


#Sunflower 
fishnet$SunflowerGCO<-joined$SunflowerGCO
mean.sunflower<-subset(fishnet, fishnet$mean_sunfl > 0 ) #yield > 0 
mean.sunflower.data<-mean.sunflower@data

#Links
nb <- poly2nb(mean.sunflower)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
sunflower.binary.sper <- errorsarlm(log10(mean_sunfl) ~ SunflowerGCO, data=mean.sunflower.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(sunflower.binary.sper) #Summary
mean.sunflower$residualsSpecError <- residuals(sunflower.binary.sper) #Residuals
moran.sunflower<-moran.mc(mean.sunflower$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.sunflower)


#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
sunflower<-na.omit(yield.mapping %>% dplyr::select(sunflower_HgHa,SunflowerGCO))

#Inside vs Outside GCO 
sunflower.inside<-filter(sunflower, SunflowerGCO == "Inside")
sunflower.outside<-filter(sunflower, SunflowerGCO == "Outside")
sunflower.m<-melt(sunflower)
summary(sunflower.inside)
summary(sunflower.outside)

#Plot
ggplot(sunflower.m, aes(x=sunflower.m$SunflowerGCO, y=sunflower.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO')



#Wheat 
fishnet$WheatGCO<-joined$WheatGCO
mean.wheat<-subset(fishnet, fishnet$mean_wheat > 0 ) #yield > 0 
mean.wheat.data<-mean.wheat@data

#Links
nb <- poly2nb(mean.wheat)
lw <- nb2listw(nb, zero.policy = TRUE)

#Spatial Error Model
wheat.binary.sper <- errorsarlm(log10(mean_wheat) ~ WheatGCO, data=mean.wheat.data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(wheat.binary.sper) #Summary
mean.wheat$residualsSpecError <- residuals(wheat.binary.sper) #Residuals
moran.wheat<-moran.mc(mean.wheat$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation
print(moran.wheat)


#merge
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
wheat<-na.omit(yield.mapping %>% dplyr::select(wheat_HgHa,WheatGCO))

#Inside vs Outside GCO 
wheat.inside<-filter(wheat, WheatGCO == "Inside")
wheat.outside<-filter(wheat, WheatGCO == "Outside")
wheat.m<-melt(wheat)
summary(wheat.inside)
summary(wheat.outside)

#Plot
ggplot(wheat.m, aes(x=wheat.m$WheatGCO, y=wheat.m$value)) + 
  geom_jitter(color=insig)+
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="black") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO')






