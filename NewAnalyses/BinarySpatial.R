
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

#colors
sighigh<-"#00ADB6"
siglow<-"#FFA800"
insig<-"#7DA3B2"

##### Models
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
fish<-fishnet@data
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID
joined<-na.omit(inner_join(fish,mapping, by="Fishnet_ID"))

#Barley 
fishnet$BarleyGCO<-joined$BarleyGCO
mean.barley<-subset(fishnet@data, fishnet$mean_barle > 0 ) #yield > 0 

#Links
nb <- poly2nb(mean.barley)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
mean.barley<-log10(mean.barley$mean_barle)
#Spatial Error Model
barley.binary.sper <- errorsarlm(mean_barle ~ BarleyGCO, data=mean.barley, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(barley.fert.sper) #Summary
mean.barley$residualsSpecError <- residuals(barley.binary.sper) #Residuals
moran.mc(mean.barley$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)


#merge
yield.mapping<-merge(mapping,yield)
# Barley 
barley<-na.omit(yield.mapping %>% dplyr::select(barley_HgHa,BarleyGCO))

#Inside vs Outside GCO 
barley.inside<-filter(barley, BarleyGCO == "Inside")
barley.outside<-filter(barley, BarleyGCO == "Outside")
barley.m<-melt(barley)
summary(barley.inside)
summary(barley.outside)


#Plot
ggplot(barley.m, aes(x=barley.m$BarleyGCO, y=barley.m$value)) + 
  geom_boxplot(fill=sighigh, fill="transparent", width=0.4)+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Hg/Ha)") +xlab("") + labs(color='GCO') 


#Cassava 
fishnet$CassavaGCO<-joined$CassavaGCO
mean.cassava<-subset(fishnet@data, fishnet$mean_barle > 0 ) #yield > 0 

#Links
nb <- poly2nb(mean.cassava)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
mean.cassava<-log10(mean.cassava$mean_barle)
#Spatial Error Model
cassava.binary.sper <- errorsarlm(mean_barle ~ CassavaGCO, data=mean.cassava, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(cassava.fert.sper) #Summary
mean.cassava$residualsSpecError <- residuals(cassava.binary.sper) #Residuals
moran.mc(mean.cassava$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)


#merge
yield.mapping<-merge(mapping,yield)
# Cassava 
cassava<-na.omit(yield.mapping %>% dplyr::select(cassava_HgHa,CassavaGCO))

#Inside vs Outside GCO 
cassava.inside<-filter(cassava, CassavaGCO == "Inside")
cassava.outside<-filter(cassava, CassavaGCO == "Outside")
cassava.m<-melt(cassava)
summary(cassava.inside)
summary(cassava.outside)


#Plot
ggplot(cassava.m, aes(x=cassava.m$CassavaGCO, y=cassava.m$value)) + 
  geom_boxplot(fill=sighigh, fill="transparent", width=0.4)+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Hg/Ha)") +xlab("") + labs(color='GCO')


#Groundnut 
fishnet$GroundnutGCO<-joined$GroundnutGCO
mean.groundnut<-subset(fishnet@data, fishnet$mean_barle > 0 ) #yield > 0 

#Links
nb <- poly2nb(mean.groundnut)
lw <- nb2listw(nb, zero.policy = TRUE)
#Log Transform yields
mean.groundnut<-log10(mean.groundnut$mean_barle)
#Spatial Error Model
groundnut.binary.sper <- errorsarlm(mean_barle ~ GroundnutGCO, data=mean.groundnut, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
summary(groundnut.fert.sper) #Summary
mean.groundnut$residualsSpecError <- residuals(groundnut.binary.sper) #Residuals
moran.mc(mean.groundnut$residualsSpecError, lw, 999,zero.policy = TRUE) #Test for autocorrelation

mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)


#merge
yield.mapping<-merge(mapping,yield)
# Groundnut 
groundnut<-na.omit(yield.mapping %>% dplyr::select(groundnut_HgHa,GroundnutGCO))

#Inside vs Outside GCO 
groundnut.inside<-filter(groundnut, GroundnutGCO == "Inside")
groundnut.outside<-filter(groundnut, GroundnutGCO == "Outside")
groundnut.m<-melt(groundnut)
summary(groundnut.inside)
summary(groundnut.outside)


#Plot
ggplot(groundnut.m, aes(x=groundnut.m$GroundnutGCO, y=groundnut.m$value)) + 
  geom_boxplot(fill=sighigh, fill="transparent", width=0.4)+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Hg/Ha)") +xlab("") + labs(color='GCO')


