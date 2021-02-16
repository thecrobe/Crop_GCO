library(rgdal)

###Import Fishnet

fishnet<- readOGR(dsn= "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/GIS/Fishnet/", layer="Fishnet_NoAntarctica")

### Get Lat-Long Baby Girl 
xy <- coordinates(fishnet)
metadata <- fishnet@data
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
df<-cbind(xy.df,distances)
#write.csv(x = df, file="distance_australia.csv") #Export
