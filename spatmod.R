
#models
library(wPerm)
library(caret)
library(wPerm)
library(spatialEco)
library(rgdal)
library(RVAideMemoire)
library(nlme)
library(rgeos)

#plots
library(RColorBrewer)
library(ggplot2)
library(wesanderson)

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

fishnet<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/Fishnet/", layer="Fishnet_yield_NoAntarctica")
fishnet.barley<-subset(fishnet, fishnet$mean_barle > 0 ) #yield > 0 
xy <- coordinates(fishnet.barley)
yield <- fishnet.barley@data

barleyGCO<-readOGR(dsn = "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/GCO_Shapefiles/",layer="BarleyWheat_GC")

barleyGCOcentroid<-gCentroid(barleyGCO)

library(dplyr)
distances <- spDistsN1(xy,barleyGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
barley.near<-yield %>% dplyr::select(mean_barle,rescale_ND,COUNTRY)
barley.near<-barley.near %>% filter(mean_barle > 0, na.rm=TRUE)
barley.near<-barley.near %>% filter(rescale_ND > 0, na.rm=TRUE)
barley.near$logBarley<-log10(barley.near$mean_barle)


library(latticeExtra)
## Loading required package: lattice
library(RColorBrewer)
grps <- 10
brks <- quantile(fishnet.barley$mean_barle, 0:(grps-1)/(grps-1), na.rm=TRUE)
p <- spplot(fishnet.barley, "mean_barle", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")), col="transparent" )
## Warning in wkt(obj): CRS object has no comment
## Warning in wkt(obj): CRS object has no comment
p + layer(sp.polygons(fishnet.barley))

f1 <- mean_barle ~ barley.near$rescale_ND 
m1 <- lm(f1, data=fishnet.barley)
summary(m1)

fishnet.barley$residuals <- residuals(m1)
brks <- quantile(fishnet.barley$residuals, 0:(grps-1)/(grps-1), na.rm=TRUE)
spplot(fishnet.barley, "residuals", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")), col="transparent")
## Warning in wkt(obj): CRS object has no comment
## Warning in wkt(obj): CRS object has no comment

nb <- poly2nb(fishnet.barley)
resnb <- sapply(nb, function(x) mean(fishnet.barley$residuals[x]))
cor(fishnet.barley$residuals, resnb)
## [1] 0.6311218
plot(fishnet.barley$residuals, resnb, xlab='Residuals', ylab='Mean adjacent residuals')
lw <- nb2listw(nb, zero.policy = TRUE)
moran.mc(fishnet.barley$residuals, lw, 999,zero.policy = TRUE)
m1e <- errorsarlm(f1, data=fishnet.barley, lw, tol.solve=1.0e-30,zero.policy  = TRUE)

