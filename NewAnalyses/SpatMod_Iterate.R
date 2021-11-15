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
  moran<-moran.mc(maize.sper$residuals, lw, 1000, zero.policy = T)
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


#Rice data prep
fishnet.Rice<-subset(fishnet, fishnet$mean_rice_ > 0 ) #yield > 0 
fishnet.Rice$logRice<-log10(fishnet.Rice$mean_rice_)
xy <- coordinates(fishnet.Rice)

#Neighborhoods
nb <- poly2nb(fishnet.Rice,queen = TRUE)
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
  fishnet.Rice$distances<-dist/(1000*1000)
  #fit model 
  Rice.sper<-errorsarlm(logRice ~ distances, listw = lw, data = fishnet.Rice, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(Rice.sper)
  moran<-moran.mc(Rice.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(Rice.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/Rice_null.csv")


#Rye data prep
fishnet.Rye<-subset(fishnet, fishnet$mean_rye_Y > 0 ) #yield > 0 
fishnet.Rye$logRye<-log10(fishnet.Rye$mean_rye_Y)
xy <- coordinates(fishnet.Rye)

#Neighborhoods
nb <- poly2nb(fishnet.Rye,queen = TRUE)
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
  fishnet.Rye$distances<-dist/(1000*1000)
  #fit model 
  Rye.sper<-errorsarlm(logRye ~ distances, listw = lw, data = fishnet.Rye, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(Rye.sper)
  moran<-moran.mc(Rye.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(Rye.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/Rye_null.csv")

#Sorghum data prep
fishnet.Sorghum<-subset(fishnet, fishnet$mean_sorgh > 0 ) #yield > 0 
fishnet.Sorghum$logSorghum<-log10(fishnet.Sorghum$mean_sorgh)
xy <- coordinates(fishnet.Sorghum)

#Neighborhoods
nb <- poly2nb(fishnet.Sorghum,queen = TRUE)
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
  fishnet.Sorghum$distances<-dist/(1000*1000)
  #fit model 
  Sorghum.sper<-errorsarlm(logSorghum ~ distances, listw = lw, data = fishnet.Sorghum, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(Sorghum.sper)
  moran<-moran.mc(Sorghum.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(Sorghum.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/Sorghum_null.csv")

#Soybean data prep
fishnet.Soybean<-subset(fishnet, fishnet$mean_soybe > 0 ) #yield > 0 
fishnet.Soybean$logSoybean<-log10(fishnet.Soybean$mean_soybe)
xy <- coordinates(fishnet.Soybean)

#Neighborhoods
nb <- poly2nb(fishnet.Soybean,queen = TRUE)
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
  fishnet.Soybean$distances<-dist/(1000*1000)
  #fit model 
  Soybean.sper<-errorsarlm(logSoybean ~ distances, listw = lw, data = fishnet.Soybean, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(Soybean.sper)
  moran<-moran.mc(Soybean.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(Soybean.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/Soybean_null.csv")


#Sunflower data prep
fishnet.Sunflower<-subset(fishnet, fishnet$mean_sunfl > 0 ) #yield > 0 
fishnet.Sunflower$logSunflower<-log10(fishnet.Sunflower$mean_sunfl)
xy <- coordinates(fishnet.Sunflower)
fishnet$
#Neighborhoods
nb <- poly2nb(fishnet.Sunflower,queen = TRUE)
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
  fishnet.Sunflower$distances<-dist/(1000*1000)
  #fit model 
  Sunflower.sper<-errorsarlm(logSunflower ~ distances, listw = lw, data = fishnet.Sunflower, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(Sunflower.sper)
  moran<-moran.mc(Sunflower.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(Sunflower.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/Sunflower_null.csv")


#Wheat data prep
fishnet.Wheat<-subset(fishnet, fishnet$mean_wheat > 0 ) #yield > 0 
fishnet.Wheat$logWheat<-log10(fishnet.Wheat$mean_wheat)
xy <- coordinates(fishnet.Wheat)

  #Neighborhoods
nb <- poly2nb(fishnet.Wheat,queen = TRUE)
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
  fishnet.Wheat$distances<-dist/(1000*1000)
  #fit model 
  Wheat.sper<-errorsarlm(logWheat ~ distances, listw = lw, data = fishnet.Wheat, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(Wheat.sper)
  moran<-moran.mc(Wheat.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(Wheat.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/Wheat_null.csv")


#Potato data prep
fishnet.Potato<-subset(fishnet, fishnet$mean_potat > 0 ) #yield > 0 
fishnet.Potato$logPotato<-log10(fishnet.Potato$mean_potat)
xy <- coordinates(fishnet.Potato)

#Neighborhoods
nb <- poly2nb(fishnet.Potato,queen = TRUE)
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
  fishnet.Potato$distances<-dist/(1000*1000)
  #fit model 
  Potato.sper<-errorsarlm(logPotato ~ distances, listw = lw, data = fishnet.Potato, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(Potato.sper)
  moran<-moran.mc(Potato.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(Potato.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/Potato_null.csv")


#Soybean data prep
fishnet.Soybean<-subset(fishnet, fishnet$mean_soybe > 0 ) #yield > 0 
fishnet.Soybean$logSoybean<-log10(fishnet.Soybean$mean_soybe)
xy <- coordinates(fishnet.Soybean)

#Neighborhoods
nb <- poly2nb(fishnet.Soybean,queen = TRUE)
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
  fishnet.Soybean$distances<-dist/(1000*1000)
  #fit model 
  Soybean.sper<-errorsarlm(logSoybean ~ distances, listw = lw, data = fishnet.Soybean, quiet=TRUE,zero.policy = TRUE)
  sumr<-summary(Soybean.sper)
  moran<-moran.mc(Soybean.sper$residuals, lw, 1000, zero.policy = T)
  # associate it to the point
  output[i,] <- c(coef(Soybean.sper)[3],sumr$Wald1[3], xy[random.point,][1],xy[random.point,][2], moran$p.value)
}
gc()
ggplot(output, aes(x=output$x,y=output$y,color=log10(output$slope))) +geom_point()
write.csv(x = output, file="NewAnalyses/Null/Soybean_null.csv")























