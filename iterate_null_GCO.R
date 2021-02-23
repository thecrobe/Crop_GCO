library(nlme)
library(rgdal)
library(dplyr)

#### readins and preprocess ######################
fishnet<- readOGR(dsn= "Fishnets", layer="Fishnet_yield_world")

### Get Lat-Long Baby Girl 
xy <- coordinates(fishnet)
metadata <- fishnet@data
head(data.frame(metadata))
###Generate Random Point GCO
random.point <- sample(1:nrow(xy), 1)
###Calculate Distances
distances <- spDistsN1(xy,xy[random.point,], longlat = FALSE)

metadata.df<-data.frame(metadata)
yield<-cbind(distances,metadata.df)
yield$rescale_ND<-yield$distances/(1000*1000)

########################################### END READINS


# the loop

N <-1000
output <- data.frame(slope=numeric(N), x = numeric(N), y = numeric(N), P=numeric(N))
for(i in 1:N){
  # pick a point
  random.point <- sample(1:nrow(xy), 1)
  # calculate distances from that point
  distances <- spDistsN1(xy,xy[random.point,], longlat = FALSE)
  yield$rescale_ND<-distances/(1000*1000)
  sorghum.near<-yield %>% select(mean_sorgh,rescale_ND)
  sorghum.near<-sorghum.near %>% filter(mean_sorgh > 0, na.rm=TRUE)
  sorghum.near<-sorghum.near %>% filter(rescale_ND > 0, na.rm=TRUE)
  sorghum.near$logSorghum<-log10(sorghum.near$mean_sorgh)
  
  # fit yield ~ distance model
  sorghum.ols <- gls(logSorghum ~ rescale_ND, data = sorghum.near)
  
  # associate it to the point
  output[i,] <- c(coef(sorghum.ols)[2], xy[random.point,], summary(sorghum.ols)$tTable[2,4] )
  
}


# plot it!
plot(output[,2:3], type="n")
text(output$x, output$y, labels = round(output$slope,2), data = output, cex = 0.8)





