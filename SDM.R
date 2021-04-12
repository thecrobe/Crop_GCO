library("sp")
library("raster")
library("maptools")
library("rgdal")
library("dismo")

#Soybean Observations
obs.data<-read.csv(file = "SDM/Sorghum_Pests.csv", header=T)
summary(obs.data)

# Determine geographic extent of our data
max.lat <- ceiling(max(obs.data$Latitude))
min.lat <- floor(min(obs.data$Latitude))
max.lon <- ceiling(max(obs.data$Longitude))
min.lon <- floor(min(obs.data$Longitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

# Only pull out those columns of interest and in the order wanted
obs.data <- obs.data[, c("Longitude", "Latitude")]

#Import Raster 
#generate a list of input rasters ("grids")
#pattern = "*.tif$" - filters for main raster files only and skips any associated files (e.g. world files)
grids <- list.files("SDM/" , pattern = "*.tif$")
#create a raster stack from the input raster files 
s <- raster::stack(paste0("SDM/", grids))

# Crop the bioclim data to geographic extent of pests
bioclim.data <- crop(x = s, y = geographic.extent)

# Create pseudo-absence, or background, points
# Use the bioclim data files for sampling resolution
tif.files <- list.files(path = "SDM/", 
                        pattern = "*.tif$", 
                        full.names = TRUE)
# We only need one file, so use the first one in the list of .tif files
mask <- raster(tif.files[1])

# Randomly sample points (same number as our observed points)
background <- randomPoints(mask = mask,     # Provides resolution of sampling points
                           n = nrow(obs.data),      # Number of random points
                           ext = geographic.extent, # Spatially restricts sampling
                           extf = 1.25)             # Expands sampling a little bit

# Arbitrarily assign group 3 as the testing data group
testing.group <- 3

# Create vector of group memberships
group.presence <- kfold(x = obs.data, k = 5) # kfold is in dismo package

# Separate observations into training and testing groups
presence.train <- obs.data[group.presence != testing.group, ]
presence.test <- obs.data[group.presence == testing.group, ]

# Repeat the process for pseudo-absence points
group.background <- kfold(x = background, k = 5)
background.train <- background[group.background != testing.group, ]
background.test <- background[group.background == testing.group, ]

# Build a model using training data
bc.model <- bioclim(x = bioclim.data, p = presence.train)

# Predict presence from model
predict.presence <- dismo::predict(object = bc.model, 
                                   x = bioclim.data, 
                                   ext = geographic.extent)

plot(predict.presence) #check it out 

# Use testing data for model evaluation
bc.eval <- evaluate(p = presence.test,   # The presence testing data
                    a = background.test, # The absence testing data
                    model = bc.model,    # The model we are evaluating
                    x = bioclim.data)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc.threshold <- threshold(x = bc.eval, stat = "spec_sens")

# Load map data for plotting
data(wrld_simpl)

# Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict.presence > bc.threshold, 
     add = TRUE, 
     legend = TRUE)

x<-predict.presence>bc.threshold


# And add those observations
points(x = obs.data$Longitude, 
       y = obs.data$Latitude, 
       col = "black",
       pch = "+", 
       cex = 0.6)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

#Export for mapping in ArcGis (PRETTIER) 
