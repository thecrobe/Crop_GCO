
#models
library(wPerm)
library(caret)
library(wPerm)
library(spatialEco)
library(rgdal)
library(RVAideMemoire)
library(nlme)
library(rgeos)
library(spdep)

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

fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
fishnet.barley<-subset(fishnet, fishnet$mean_barle > 0 ) #yield > 0 
xy <- coordinates(fishnet.barley)
yield <- fishnet.barley@data

barleyGCO<-readOGR(dsn = "GIS/",layer="BarleyWheat_GC")

barleyGCOcentroid<-gCentroid(barleyGCO)

library(dplyr)
distances <- spDistsN1(xy,barleyGCOcentroid, longlat = FALSE)
yield$rescale_ND<-distances/(1000*1000)
barley.near<-yield %>% dplyr::select(mean_barle,rescale_ND,COUNTRY)
barley.near<-barley.near %>% filter(mean_barle > 0, na.rm=TRUE)
barley.near<-barley.near %>% filter(rescale_ND > 0, na.rm=TRUE)
barley.near$logBarley<-log10(barley.near$mean_barle)


library(latticeExtra)
library(RColorBrewer)

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
resnb.se <- sapply(nb, function(x) mean(fishnet.barley$residualsSpecError[x])) #avg value for polygon neighbors
plot(fishnet.barley$residualsSpecError, resnb.se, xlab='Residuals', ylab='Mean adjacent residuals') #check plot
brks <- quantile(fishnet.barley$residualsSpecError, 0:(grps-1)/(grps-1), na.rm=TRUE) #apply breaks
p <- spplot(fishnet.barley, "residualsSpecError", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")),
            col="transparent") 
print( p + layer(sp.polygons(fishnet)) ) #plot residuals


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
ggplot(barley.near, aes(x=rescale_ND, y=barley.near$mean_barle)) +
  geom_point(alpha=0.3) +theme_justin + 
  geom_point(data=china, aes(x=rescale_ND, y=log10(china$mean_barle)), color=china.color, alpha=0.9) + 
  geom_point(data=usa, aes(x=rescale_ND, y=log10(usa$mean_barle)), color=usa.color, alpha=0.5) + 
  geom_point(data=russia, aes(x=rescale_ND, y=log10(russia$mean_barle)),color=russia.color, alpha=0.2) +
  geom_abline(aes(intercept=coef(barley.sper)[1], slope=coef(barley.sper)[3]), color="red", size=2) +labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")

