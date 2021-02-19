library(rgdal)

###Import Fishnet

fishnet<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/Fishnet/", layer="Fishnet_yield_NoAntarctica")

### Get Lat-Long Baby Girl 
xy <- coordinates(fishnet)
metadata <- fishnet@data
head(data.frame(metadata))
###Generate Random Point GCO
random.point <- sample(1:nrow(xy), 1)
###Calculate Distances
distances <- spDistsN1(xy,xy[random.point,], longlat = FALSE)

###Plot###
xy.df<-data.frame(xy) #Set XY into a DF
#Basemap
plot(X2 ~ X1, xy.df, cex=0.01) 
#Points with color gradient
fun_color_range <- colorRampPalette(c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529'))
my_colors <- fun_color_range(100)  
points(X2 ~ X1, data= xy.df, pch= 16, col = my_colors[as.numeric(cut(distances, breaks = 100))])
#Add Random Point
points(X2 ~ X1, data= xy.df[random.point,], pch= 16, cex=1)

#Combine DF to keep FishnetID
metadata.df<-data.frame(metadata)
yield<-cbind(distances,metadata.df)
yield$rescale_ND<-yield$distances/(1000*1000)

#Sorghum Random GCO - 1
sorghum.near<-yield %>% select(mean_sorgh,rescale_ND)
sorghum.near<-sorghum.near %>% filter(mean_sorgh > 0, na.rm=TRUE)
sorghum.near<-sorghum.near %>% filter(rescale_ND > 0, na.rm=TRUE)
sorghum.near$logSorghum<-log10(sorghum.near$mean_sorgh)
summary(sorghum.near)
#Plot
ggplot(sorghum.near, aes(x=rescale_ND, y=logSorghum)) +
  geom_point(alpha=0.1, color="#00A9B4") + geom_smooth(method = "lm",color="#DB3500") + labs(color='GCO') +theme_justin +xlab("Distance From GCO (1000's km)") +ylab ("Log10 Yield")

###  Model - Sorghum Near Null - 1 
library(nlme)
# regular OLS no variance structure
sorghum.ols <- gls(logSorghum ~ rescale_ND, data = sorghum.near)
# varFixed (variance changes linearly with X)
sorghum.fixed <- update(sorghum.ols, .~., weights = varFixed(~rescale_ND))
# varPower (variance changes as a power function with X)
sorghum.power <- update(sorghum.ols, . ~ ., weights = varPower(form = ~rescale_ND))
# varExp (variance changes as an exponential function of x)
sorghum.exp <- update(sorghum.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error
# varConstPower (constant plus a power function of X (useful if X includes 0))
sorghum.ConstPower <- update(sorghum.ols, . ~., weights = varConstPower(form = ~ rescale_ND))
# compare all models by AIC
AIC(sorghum.ols, sorghum.fixed, sorghum.power, sorghum.ConstPower, sorghum.exp)
#Exp
summary(sorghum.exp)


######### IDEAS!! 
# Map the slopes for each crop's randomized GCO by Color (+/-) and Hue (Magnitude)
# Graph of Slopes vs Distance from each crop's GCO 
