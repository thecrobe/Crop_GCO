library(rgdal)
library(dplyr)
library(brms)
library(spdep)
library(spdplyr)

#Read In
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
proj4string(fishnet) <- CRS("+init=epsg:3786")
barleymodel<-read.csv(file="Models/Barley_RF.csv", header=T)
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID

m <- sp::merge(fishnet, barleymodel, by='Fishnet_ID')
barley<-sp::merge(m, mapping, by='Fishnet_ID')
barley.p = subset(barley, barley_HgHa > 0)
county.pivot<-barley.p@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())
barley<-sp::merge(barley.p,county.pivot, by="COUNTRY.x")
barley.final<-na.omit(filter(barley,CountyPixelCount > 500))

W.nb <- poly2nb(barley.final, row.names = barley.final@data$OBJECTID)
W <- nb2mat(W.nb, style="B",zero.policy = TRUE)


fit<- brm(log10(barley_HgHa) ~ Pesticide +car(gr=NA, type="icar",M = W), data=barley.final@data, data2=list(W=W), family=gaussian, iter=500, thin=5, cores=2, chains=1, seed=10)

summary(fit)
