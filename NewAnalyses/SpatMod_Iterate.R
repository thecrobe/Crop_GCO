library(spdep)
library(spatialreg)
library(rgdal)
library(ggplot2)
library(foreach)
library(doParallel)
library(rgeos)
library(dplyr)

#ggplot graphic parameters
theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

fishnet<- readOGR(dsn= "/GIS/Fishnet/", layer="Fishnet_yield_NoAntarctica")

#Barley
fishnet.barley<-subset(fishnet, fishnet$mean_barle > 0 ) #yield > 0 
xy <- coordinates(fishnet.barley)
xy1<-data.frame(xy)
yield <- fishnet.barley@data
barleyGCO<-readOGR(dsn = "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/GCO_Shapefiles/",layer="BarleyWheat_GC")
barleyGCOcentroid<-gCentroid(barleyGCO)

#Register cluster
cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)
N=100 #sample size
output <- data.frame(slope=numeric(N), x = numeric(N), y = numeric(N))
iteration<-foreach(i = 1:N) %dopar% {
  #load package 
  library(dplyr)
  library(sp)
  # pick a point
  random.point <- sample(1:nrow(xy), 1)
  # calculate distances from that point
  distances <- spDistsN1(xy,xy[random.point,], longlat = FALSE)
  yield$rescale_ND<-distances/(1000*1000)
  barley.near<-yield %>% select(mean_barle,rescale_ND)
  barley.near<-barley.near %>% filter(mean_barle > 0, na.rm=TRUE)
  barley.near<-barley.near %>% filter(rescale_ND > 0, na.rm=TRUE)
  barley.near$logBarley<-log10(barley.near$mean_barle)
  #fit model 
  barley.sper<-lm(logBarley ~ rescale_ND, data=barley.near)
  # associate it to the point
  output[i,] <- c(coef(barley.sper)[2], xy[random.point,][1],xy[random.point,][2])
}
x<-data.frame(do.call(rbind,iteration))

ggplot(x, aes(x=x$V2,y=x$V3,color=x$rescale_ND)) +geom_point()
