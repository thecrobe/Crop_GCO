library(dplyr)
library(rgdal)
library(ggplot2)

##### Soybean
#data read-ins
data1<-read.csv(file="YieldGap/SoyCountry.csv", header=T)
x<-data1
x[x==0] <- NA
summary(x)
data<-na.omit(x)
data<-x

#Import Fishnet
fishnet<- readOGR(dsn= "Shapefiles/", layer="SoybeanGap")
### Get Lat-Long Baby Girl and apply metadata
xy.soy <- coordinates(fishnet)
metadata <- fishnet@data
sapply(metadata, "class")

soy<- readOGR(dsn= "Shapefiles/", layer="SoybeanCentroid") #Import GCO Centroid
plot(soy)
plot(fishnet)
soy.c<-coordinates(soy) #Extract Lat-Long
distances <- spDistsN1(xy.soy,soy.c, longlat = FALSE) #Calculate Near Distances
data<-cbind(data.frame(metadata),data.frame(distances)) #Merge into a DF
data$rescale_ND<-data$distances/(1000*1000) #Rescale distances 
data<-data %>% filter(data$MEAN> 0, na.rm=TRUE) #select  > 0 

plot(density(data$MEAN))
qqnorm(data$MEAN)

plot(density(data$rescale_ND))
qqnorm(data$rescale_ND)

ggplot(data, aes(x=data$rescale_ND, y=data$MEAN)) + 
  geom_point(color="#8F5E77",alpha=0.3) +theme_justin +
  geom_abline(aes(intercept=coef(soy.gap.power)[1], slope=coef(soy.gap.power)[2]), color="black", size=2) +xlab("Distance From GCO (1000's km)") +ylab("Mean Yield Gap (tons/ha)")


##  Model - Barley Fertilizer Near
library(nlme)
# regular OLS no variance structure
soy.gap.ols <- gls(MEAN ~ rescale_ND, data = data)
# varFixed (variance changes linearly with X)
soy.gap.fixed <- update(soy.gap.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
soy.gap.power <- update(soy.gap.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
soy.gap.exp <- update(soy.gap.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
soy.gap.ConstPower <- update(soy.gap.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(soy.gap.ols, soy.gap.fixed, soy.gap.power, soy.gap.ConstPower, soy.gap.exp)
summary(soy.gap.power)
coef(soy.gap.exp) #for every 1000km escaped, yield gaps decrease 0.023 tons/hc
summary(soy.gap.exp)

##### Sorghum

#data read-ins
data1<-read.csv(file="YieldGap/SorghumCountry.csv", header=T)
x<-data1
x[x==0] <- NA
summary(x)
data<-na.omit(x)
data<-x

#Import Fishnet
fishnet<- readOGR(dsn= "Shapefiles/", layer="SorghumGap")
### Get Lat-Long Baby Girl and apply metadata
xy.sorghum <- coordinates(fishnet)
metadata <- fishnet@data
sapply(metadata, "class")

sorghum<- readOGR(dsn= "Shapefiles/", layer="SorghumCentroid") #Import GCO Centroid
plot(sorghum)
plot(fishnet)
sorghum.c<-coordinates(sorghum) #Extract Lat-Long
distances <- spDistsN1(xy.sorghum,sorghum.c, longlat = FALSE) #Calculate Near Distances
data<-cbind(data.frame(metadata),data.frame(distances)) #Merge into a DF
data$rescale_ND<-data$distances/(1000*1000) #Rescale distances 
data<-data %>% filter(data$MEAN> 0, na.rm=TRUE) #select  > 0 

plot(density(data$MEAN))
qqnorm(data$MEAN)

plot(density(data$rescale_ND))
qqnorm(data$rescale_ND)

ggplot(data, aes(x=data$rescale_ND, y=data$MEAN)) + 
  geom_point(color="#D9C5A9",alpha=0.3) +theme_justin +
  geom_abline(aes(intercept=coef(sorghum.gap.power)[1], slope=coef(sorghum.gap.power)[2]), color="black", size=2) +xlab("Distance From GCO (1000's km)") +ylab("Mean Yield Gap (tons/ha)")


##  Model - Barley Fertilizer Near
library(nlme)
# regular OLS no variance structure
sorghum.gap.ols <- gls(MEAN ~ log10(rescale_ND), data = data)
# varFixed (variance changes linearly with X)
sorghum.gap.fixed <- update(sorghum.gap.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
sorghum.gap.power <- update(sorghum.gap.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
sorghum.gap.exp <- update(sorghum.gap.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
sorghum.gap.ConstPower <- update(sorghum.gap.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(sorghum.gap.ols, sorghum.gap.fixed, sorghum.gap.power, sorghum.gap.ConstPower, sorghum.gap.exp)
summary(sorghum.gap.power)
coef(sorghum.gap.exp) #for every 1000km escaped, yield gaps decrease 0.023 tons/hc
# NOT SIGNIFICANT




