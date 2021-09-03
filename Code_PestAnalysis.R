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

########## Sorghum ##########
#Import Fishnet
sorghum<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/Pests/", layer="Sorghum_Pests")
### Get Lat-Long and apply metadata
sorghum.xy <- coordinates(sorghum)
sorghum.metadata <- sorghum@data
sapply(sorghum.metadata, "class")


sorghum.gco<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/GCO_Shapefiles/", layer="Sorghum_GC") #Import GCO Centroid
plot(sorghum.gco)
sorghum.c<-gCentroid(sorghum.gco) #Extract Lat-Long
distances <- spDistsN1(sorghum.xy,sorghum.c,longlat = FALSE) #Calculate Near Distances
data<-cbind(data.frame(sorghum.metadata),data.frame(distances)) #Merge into a DF
data$rescale_ND<-data$distances/(1000*1000) #Rescale distances 

s.achatina<-data %>% filter(data$Achatina > 0.5, na.rm=TRUE) #select greater than 70
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
summary(Achatina.ols)

ggplot(s.achatina, aes(x=rescale_ND, y=Achatina)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Achatina.ConstPower)[1], slope=coef(Achatina.ConstPower)[2]), color="blue", size=2) +
  theme_justin 

s.Paspalum<- data %>% filter(Paspalum > .5, na.rm=TRUE) #select greater than 70
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
  geom_abline(aes(intercept=coef(Paspalum.ols)[1], slope=coef(Paspalum.ols)[2]), color="blue", size=2) +
  theme_justin

s.Spodoptera<-data %>% filter(data$Spodoptera > .5, na.rm=TRUE) #select greater than 70
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
  geom_abline(aes(intercept=coef(Spodoptera.ols)[1], slope=coef(Spodoptera.ols)[2]), color="blue", size=2) +
  theme_justin


s.Verbesina<-data %>% filter(data$Verbesina > .5, na.rm=TRUE) #select greater than 70
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
summary(Verbesina.exp)

ggplot(s.Verbesina, aes(x=rescale_ND, y=Verbesina)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Verbesina.exp)[1], slope=coef(Verbesina.exp)[2]), color="blue", size=2) +
  theme_justin

s.Xanthium<-data %>% filter(data$Xanthium > .5, na.rm=TRUE) #select greater than 70
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
summary(Xanthium.exp)

ggplot(s.Xanthium, aes(x=rescale_ND, y=Xanthium)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Xanthium.exp)[1], slope=coef(Xanthium.exp)[2]), color="blue", size=2) +
  theme_justin

########## Soybean ##########

### Sorghum
#Import Fishnet
soybean<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/Pests/", layer="Soybean_Pests")
### Get Lat-Long and apply metadata
soybean.xy <- coordinates(soybean)
soybean.metadata <- soybean@data
sapply(soybean.metadata, "class")

soybean.gco<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/GCO_Shapefiles/", layer="Soy_GC") #Import GCO Centroid
plot(soybean.gco)
soybean.c<-gCentroid(soybean.gco) #Extract Lat-Long
distances <- spDistsN1(soybean.xy,soybean.c,longlat = FALSE) #Calculate Near Distances
data<-cbind(data.frame(soybean.metadata),data.frame(distances)) #Merge into a DF
data$rescale_ND<-data$distances/(1000*1000) #Rescale distances 


s.Diabrotica<-data %>% filter(data$Diabrotica > 0.5, na.rm=TRUE) #select greater than 70
dim(s.Diabrotica) # data points

library(nlme)
# regular OLS no variance structure
Diabrotica.ols <- gls(Diabrotica ~ rescale_ND, data = s.Diabrotica)
# varFixed (variance changes linearly with X)
Diabrotica.fixed <- update(Diabrotica.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Diabrotica.power <- update(Diabrotica.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Diabrotica.exp <- update(Diabrotica.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Diabrotica.ConstPower <- update(Diabrotica.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Diabrotica.ols, Diabrotica.fixed, Diabrotica.power, Diabrotica.ConstPower, Diabrotica.exp)
summary(Diabrotica.ols)

ggplot(s.Diabrotica, aes(x=rescale_ND, y=Diabrotica)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Diabrotica.ols)[1], slope=coef(Diabrotica.ols)[2]), color="grey", size=2) +
  theme_justin 

s.Pseudo<-data %>% filter(data$Pseudo > 0.5, na.rm=TRUE) #select greater than 70
dim(s.Pseudo) # data points

library(nlme)
# regular OLS no variance structure
Pseudo.ols <- gls(Pseudo ~ rescale_ND, data = s.Pseudo)
# varFixed (variance changes linearly with X)
Pseudo.fixed <- update(Pseudo.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Pseudo.power <- update(Pseudo.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Pseudo.exp <- update(Pseudo.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Pseudo.ConstPower <- update(Pseudo.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Pseudo.ols, Pseudo.fixed, Pseudo.power, Pseudo.ConstPower, Pseudo.exp)
summary(Pseudo.power)

ggplot(s.Pseudo, aes(x=rescale_ND, y=Pseudo)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Pseudo.power)[1], slope=coef(Pseudo.power)[2]), color="blue", size=2) +
  theme_justin

s.Rhizoc<-data %>% filter(data$Rhizoc > 0.5, na.rm=TRUE) #select greater than 70
dim(s.Rhizoc) # data points

library(nlme)
# regular OLS no variance structure
Rhizoc.ols <- gls(Rhizoc ~ rescale_ND, data = s.Rhizoc)
# varFixed (variance changes linearly with X)
Rhizoc.fixed <- update(Rhizoc.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Rhizoc.power <- update(Rhizoc.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Rhizoc.exp <- update(Rhizoc.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Rhizoc.ConstPower <- update(Rhizoc.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Rhizoc.ols, Rhizoc.fixed, Rhizoc.power, Rhizoc.ConstPower, Rhizoc.exp)
summary(Rhizoc.ols)

ggplot(s.Rhizoc, aes(x=rescale_ND, y=Rhizoc)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Rhizoc.ols)[1], slope=coef(Rhizoc.ols)[2]), color="grey", size=2) +
  theme_justin

s.Spodop<-data %>% filter(data$Spodop > 0.5, na.rm=TRUE) #select greater than 70
dim(s.Spodop) # data points

library(nlme)
# regular OLS no variance structure
Spodop.ols <- gls(Spodop ~ rescale_ND, data = s.Spodop)
# varFixed (variance changes linearly with X)
Spodop.fixed <- update(Spodop.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Spodop.power <- update(Spodop.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Spodop.exp <- update(Spodop.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Spodop.ConstPower <- update(Spodop.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Spodop.ols, Spodop.fixed, Spodop.power, Spodop.ConstPower, Spodop.exp)
summary(Spodop.ols)

ggplot(s.Spodop, aes(x=rescale_ND, y=Spodop)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Spodop.ols)[1], slope=coef(Spodop.ols)[2]), color="grey", size=2) +
  theme_justin

s.Xantho<-data %>% filter(data$Xantho > 0.5, na.rm=TRUE) #select greater than 70
dim(s.Xantho) # data points

library(nlme)
# regular OLS no variance structure
Xantho.ols <- gls(Xantho ~ rescale_ND, data = s.Xantho)
# varFixed (variance changes linearly with X)
Xantho.fixed <- update(Xantho.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
Xantho.power <- update(Xantho.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
Xantho.exp <- update(Xantho.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
Xantho.ConstPower <- update(Xantho.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(Xantho.ols, Xantho.fixed, Xantho.power, Xantho.ConstPower, Xantho.exp)
summary(Xantho.ols)

ggplot(s.Xantho, aes(x=rescale_ND, y=Xantho)) +
  geom_point() +
  ylab("Presence Probability") + 
  xlab("Distance From GCO (1000's km)") +
  geom_abline(aes(intercept=coef(Xantho.ols)[1], slope=coef(Xantho.ols)[2]), color="grey", size=2) +
  theme_justin


