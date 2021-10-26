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

fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")

#Barley data prep
fishnet.barley<-subset(fishnet, fishnet$mean_barle > 0 ) #yield > 0 
fishnet.barley$logBarley<-log10(fishnet.barley$mean_barle)
xy <- coordinates(fishnet.barley)
xy1<-data.frame(xy)
barleyGCO<-readOGR(dsn = "GIS/",layer="BarleyWheat_GC")
barleyGCOcentroid<-gCentroid(barleyGCO)

#Neighborhoods
nb <- poly2nb(fishnet.barley)
lw <- nb2listw(nb,zero.policy = T)
print(lw)

#Register cluster
cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)
N=100 #sample size
output <- data.frame(slope=numeric(N), x = numeric(N), y = numeric(N))
iteration<-foreach(i = 1:N) %dopar% {
  #load package 
  library(dplyr)
  library(sp) #try :: 
  library(spdep)
  # pick a point
  random.point <- sample(1:nrow(xy), 1)
  # calculate distances from that point
  distances <- spDistsN1(xy,xy[random.point,], longlat = FALSE)
  #fit model 
  barley.sper<-errorsarlm(logBarley ~ distances, data=fishnet.barley@data, lw, tol.solve=1.0e-30,zero.policy  = TRUE)
  # associate it to the point
  output[i,] <- c(coef(barley.sper)[3], xy[random.point,][1],xy[random.point,][2])

}
stopCluster(cl)
x<-data.frame(do.call(rbind,iteration))

ggplot(x, aes(x=x$V2,y=x$V3,color=x$rescale_ND)) +geom_point()
