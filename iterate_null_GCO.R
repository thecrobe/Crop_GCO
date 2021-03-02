library(nlme)
library(rgdal)
library(dplyr)
library(ggplot2)

#ggplot graphic parameters
theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank())


#### readins and preprocess ######################
fishnet<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/Fishnet/", layer="Fishnet_yield_NoAntarctica")


### Get Lat-Long Baby Girl 
xy <- coordinates(fishnet)
?coordinates
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

N <-10000
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

#Plot it!
ggplot(output, aes(x=x, y=y, color=slope)) +geom_point() +theme_justin 

#write.csv(output, file="Near_Null/Sorghum_Null_10000.csv")

#Sorghum GCO Distance Sensitivity
sorghum<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/GCO_Shapefiles/", layer="Sorghum_Centroid")



plot(sorghum)
coords<-data.matrix(cbind(output$x,output$y))
sorg.distances <- spDistsN1(coords, sorghum,longlat = FALSE)
data<-data.frame(cbind(output, sorg.distances))
data$x_rescale<-data$x/(1000*1000)

ggplot(data, aes(x=data$x_rescale, y=data$slope, color=data$slope)) +theme_justin +
  geom_point() +geom_hline(yintercept = 0.0, linetype="dashed", color="Black",size=1.5) +xlab("Longitudinal Distance From GCO (1000's KM)") +ylab("Intercept Slope") +scale_color_gradient2(low = "#582999", high = "#147A0B")
