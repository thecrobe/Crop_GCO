##### Pest Analysis Models

```{r,  message=TRUE, warning=FALSE, cache=TRUE, message=FALSE}
library(rgdal)
library(maptools)
library(rgeos)
library(ggplot2)
library(dplyr)
library(reshape2)
library(sf)

#ggplot graphic parameters
theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank())
### Sorghum
#Import Fishnet
sorghum<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/Pests/", layer="Sorghum_Pests")
### Get Lat-Long Baby Girl and apply metadata
sorghum.xy <- coordinates(sorghum)
sorghum.metadata <- sorghum@data
sapply(sorghum.metadata, "class")
plot(sorghum) #Fishnet only contains pixels with fertilizer input
```

```{r,  message=TRUE, warning=FALSE, cache=TRUE, message=FALSE}
sorghum.gco<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/GCO_Shapefiles/", layer="Sorghum_GC") #Import GCO Centroid
plot(sorghum.gco)
sorghum.c<-gCentroid(sorghum.gco) #Extract Lat-Long
distances <- spDistsN1(sorghum.xy,sorghum.c,longlat = FALSE) #Calculate Near Distances
data<-cbind(data.frame(sorghum.metadata),data.frame(distances)) #Merge into a DF
data$rescale_ND<-data$distances/(1000*1000) #Rescale distances 


s.achatina<-data %>% filter(data$Achatina > .7, na.rm=TRUE) #select greater than 70
dim(s.achatina) # data points

library(nlme)
# regular OLS no variance structure
Achatina.ols <- gls(Achatina ~ rescale_ND, data = s.achatina)
# varFixed (variance changes linearly with X)
Achatina.fixed <- update(Achatina.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Achatina.power <- update(Achatina.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Achatina.exp <- update(Achatina.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Achatina.ConstPower <- update(Achatina.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Achatina.ols, Achatina.fixed, Achatina.power, Achatina.ConstPower, Achatina.exp)
summary(barley.power)
summary(Achatina.ols)

ggplot(s.achatina, aes(x=rescale_ND, y=Achatina)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Achatina.ols)[1], slope=coef(Achatina.ols)[2]), color="blue", size=2) +
  theme_justin 
sorghum$Paspalum


s.Paspalum<-data %>% filter(data$Paspalum > .7, na.rm=TRUE) #select greater than 70
dim(s.Paspalum) # data points

library(nlme)
# regular OLS no variance structure
Paspalum.ols <- gls(Paspalum ~ rescale_ND, data = s.Paspalum)
# varFixed (variance changes linearly with X)
Paspalum.fixed <- update(Paspalum.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Paspalum.power <- update(Paspalum.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Paspalum.exp <- update(Paspalum.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Paspalum.ConstPower <- update(Paspalum.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Paspalum.ols, Paspalum.fixed, Paspalum.power, Paspalum.ConstPower, Paspalum.exp)
summary(Paspalum.ols)

ggplot(s.Paspalum, aes(x=rescale_ND, y=Paspalum)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Paspalum.ols)[1], slope=coef(Paspalum.ols)[2]), color="grey", size=2) +
  theme_justin

s.Spodoptera<-data %>% filter(data$Spodoptera > .7, na.rm=TRUE) #select greater than 70
dim(s.Spodoptera) # data points

library(nlme)
# regular OLS no variance structure
Spodoptera.ols <- gls(Spodoptera ~ rescale_ND, data = s.Spodoptera)
# varFixed (variance changes linearly with X)
Spodoptera.fixed <- update(Spodoptera.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Spodoptera.power <- update(Spodoptera.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Spodoptera.exp <- update(Spodoptera.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Spodoptera.ConstPower <- update(Spodoptera.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Spodoptera.ols, Spodoptera.fixed, Spodoptera.power, Spodoptera.ConstPower, Spodoptera.exp)
summary(Spodoptera.ols)

ggplot(s.Spodoptera, aes(x=rescale_ND, y=Spodoptera)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Spodoptera.ols)[1], slope=coef(Spodoptera.ols)[2]), color="grey", size=2) +
  theme_justin


s.Verbesina<-data %>% filter(data$Verbesina > .7, na.rm=TRUE) #select greater than 70
dim(s.Verbesina) # data points

library(nlme)
# regular OLS no variance structure
Verbesina.ols <- gls(Verbesina ~ rescale_ND, data = s.Verbesina)
# varFixed (variance changes linearly with X)
Verbesina.fixed <- update(Verbesina.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Verbesina.power <- update(Verbesina.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Verbesina.exp <- update(Verbesina.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Verbesina.ConstPower <- update(Verbesina.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Verbesina.ols, Verbesina.fixed, Verbesina.power, Verbesina.ConstPower, Verbesina.exp)
summary(Verbesina.ols)

ggplot(s.Verbesina, aes(x=rescale_ND, y=Verbesina)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Verbesina.ols)[1], slope=coef(Verbesina.ols)[2]), color="blue", size=2) +
  theme_justin

s.Xanthium<-data %>% filter(data$Xanthium > .7, na.rm=TRUE) #select greater than 70
dim(s.Xanthium) # data points

library(nlme)
# regular OLS no variance structure
Xanthium.ols <- gls(Xanthium ~ rescale_ND, data = s.Xanthium)
# varFixed (variance changes linearly with X)
Xanthium.fixed <- update(Xanthium.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Xanthium.power <- update(Xanthium.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Xanthium.exp <- update(Xanthium.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Xanthium.ConstPower <- update(Xanthium.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Xanthium.ols, Xanthium.fixed, Xanthium.power, Xanthium.ConstPower, Xanthium.exp)
summary(Xanthium.ols)

ggplot(s.Xanthium, aes(x=rescale_ND, y=Xanthium)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Xanthium.ols)[1], slope=coef(Xanthium.ols)[2]), color="blue", size=2) +
  theme_justin

```






