library(spdep)
library(spatialreg)
library(rgdal)
library(ggplot2)
library(rgeos)
library(dplyr)

#ggplot graphic parameters
theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")

#Barley data prep
fishnet.barley<-subset(fishnet, fishnet$mean_barle > 0 ) #yield > 0 
fishnet.barley$logBarley<-log10(fishnet.barley$mean_barle)
xy <- coordinates(fishnet.barley)

#Neighborhoods
nb <- poly2nb(fishnet.barley,queen = TRUE)
lw <- nb2listw(nb,zero.policy = TRUE,style = "W")

#Register cluster
library(RhpcBLASctl)
blas_set_num_threads(8)
N=100 #sample size
output <- data.frame(slope=numeric(N), pval=numeric(N), x=numeric(N), y=numeric(N), moran_i=numeric(N), moran_p=numeric(N))
iteration<-for(i in 1:N){
  # pick a point
  random.point <- sample(1:nrow(xy), 1)
  # calculate distances from that point
  dist <- spDistsN1(xy,xy[random.point,], longlat = FALSE)
  fishnet.barley$distances<-dist/(1000*1000)
  #fit model 
  barley.sper<-errorsarlm(logBarley ~ distances, listw = lw, data = fishnet.barley, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(barley.sper)
  sa<-moran.mc(barley.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(barley.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/barley_null.csv")

#Cassava data prep
fishnet.cassava<-subset(fishnet, fishnet$mean_cassa > 0 ) #yield > 0 
fishnet.cassava$logCassava<-log10(fishnet.cassava$mean_cassa)
xy <- coordinates(fishnet.cassava)

#Neighborhoods
nb <- poly2nb(fishnet.cassava,queen = TRUE)
lw <- nb2listw(nb,zero.policy = TRUE,style = "W")

#Register cluster
library(RhpcBLASctl)
blas_set_num_threads(8)
N=100 #sample size
output <- data.frame(slope=numeric(N), pval=numeric(N), x=numeric(N), y=numeric(N), moran_i=numeric(N), moran_p=numeric(N))
iteration<-for(i in 1:N){
  # pick a point
  random.point <- sample(1:nrow(xy), 1)
  # calculate distances from that point
  dist <- spDistsN1(xy,xy[random.point,], longlat = FALSE)
  fishnet.cassava$distances<-dist/(1000*1000)
  #fit model 
  cassava.sper<-errorsarlm(logCassava ~ distances, listw = lw, data = fishnet.cassava, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(cassava.sper)
  sa<-moran.mc(cassava.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(cassava.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/cassava_null.csv")


#Groundnut data prep
fishnet.groundnut<-subset(fishnet, fishnet$mean_groun > 0 ) #yield > 0 
fishnet.groundnut$logGroundnut<-log10(fishnet.groundnut$mean_groun)
xy <- coordinates(fishnet.groundnut)

#Neighborhoods
nb <- poly2nb(fishnet.groundnut,queen = TRUE)
lw <- nb2listw(nb,zero.policy = TRUE,style = "W")

#Register cluster
library(RhpcBLASctl)
blas_set_num_threads(8)
N=100 #sample size
output <- data.frame(slope=numeric(N), pval=numeric(N), x=numeric(N), y=numeric(N), moran_i=numeric(N), moran_p=numeric(N))
iteration<-for(i in 1:N){
  # pick a point
  random.point <- sample(1:nrow(xy), 1)
  # calculate distances from that point
  dist <- spDistsN1(xy,xy[random.point,], longlat = FALSE)
  fishnet.groundnut$distances<-dist/(1000*1000)
  #fit model 
  groundnut.sper<-errorsarlm(logGroundnut ~ distances, listw = lw, data = fishnet.groundnut, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(groundnut.sper)
  sa<-moran.mc(groundnut.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(groundnut.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], sa$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/groundnut_null.csv")

#Maize data prep
fishnet.maize<-subset(fishnet, fishnet$mean_maize > 0 ) #yield > 0 
fishnet.maize$logMaize<-log10(fishnet.maize$mean_maize)
xy <- coordinates(fishnet.maize)

#Neighborhoods
nb <- poly2nb(fishnet.maize,queen = TRUE)
lw <- nb2listw(nb,zero.policy = TRUE,style = "W")

#Register cluster
library(RhpcBLASctl)
blas_set_num_threads(8)
N=100 #sample size
output <- data.frame(slope=numeric(N), pval=numeric(N), x=numeric(N), y=numeric(N), moran_i=numeric(N), moran_p=numeric(N))
iteration<-for(i in 1:N){
  # pick a point
  random.point <- sample(1:nrow(xy), 1)
  # calculate distances from that point
  dist <- spDistsN1(xy,xy[random.point,], longlat = FALSE)
  fishnet.maize$distances<-dist/(1000*1000)
  #fit model 
  maize.sper<-errorsarlm(logMaize ~ distances, listw = lw, data = fishnet.maize, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(maize.sper)
  sa<-moran.mc(maize.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(maize.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/maize_null.csv")


#Rapeseed data prep
fishnet.Rapeseed<-subset(fishnet, fishnet$mean_rapes > 0 ) #yield > 0 
fishnet.Rapeseed$logRapeseed<-log10(fishnet.Rapeseed$mean_rapes)
xy <- coordinates(fishnet.Rapeseed)

#Neighborhoods
nb <- poly2nb(fishnet.Rapeseed,queen = TRUE)
lw <- nb2listw(nb,zero.policy = TRUE,style = "W")

#Register cluster
library(RhpcBLASctl)
blas_set_num_threads(8)
N=100 #sample size
output <- data.frame(slope=numeric(N), pval=numeric(N), x=numeric(N), y=numeric(N), moran_i=numeric(N), moran_p=numeric(N))
iteration<-for(i in 1:N){
  # pick a point
  random.point <- sample(1:nrow(xy), 1)
  # calculate distances from that point
  dist <- spDistsN1(xy,xy[random.point,], longlat = FALSE)
  fishnet.Rapeseed$distances<-dist/(1000*1000)
  #fit model 
  Rapeseed.sper<-errorsarlm(logRapeseed ~ distances, listw = lw, data = fishnet.Rapeseed, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(Rapeseed.sper)
  moran<-moran.mc(Rapeseed.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(Rapeseed.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/Rapeseed_null.csv")

