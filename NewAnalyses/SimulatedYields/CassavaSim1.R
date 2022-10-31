library(spatialRF) 
library(ggplot2)
library(ranger)
library(Rcpp)
library(parallel)
library(dplyr)
library(cowplot)

#import artificial yields
cassava<-read.csv(file="NewAnalyses/SimulatedYields/ArtificialCassava.csv", header=T)

#coordinates of the cases
x<-cassava$Latitude
y<-cassava$Longitude
xy <- data.frame(cbind(x,y))

#Create distance matrix
dist.mat <- data.matrix(dist(cbind(cassava$Latitude, cassava$Longitude)))

#coordinates of the cases
x<-cassava$Latitude
y<-cassava$Longitude
xy <- data.frame(cbind(x,y))

#distance matrix
distance.matrix <- dist.mat
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(1,5,10,50)
#random seed for reproducibility
random.seed <- 1
#log transform
cassava.data<-log10(cassava[,4:22]+1)

#Set up cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 4,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster)


# Spatial Model
Sim0.1x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim0.1x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim0.1x)

#Get response curves
Sim0.1x.df <- spatialRF::get_response_curves(Sim0.1x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim0.1x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim0.1xcassava_pred.csv")
write.csv(Sim0.1x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim0.1xcassava_resid.csv")
write.csv(Sim0.1x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim0.1xcassava_VarImp.csv")
write.csv(Sim0.1x$performance, file="Crop_GCO/Global/ArtificialYields//Sim0.1xcassava_preformance.csv")
write.csv(Sim0.1x.df, file="Crop_GCO/Global/ArtificialYields//Sim0.1xcassava_curves.csv")

Sim0.1x.repeat <- spatialRF::rf_repeat(
  model = Sim0.1x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)
print(Sim0.1x.repeat)
Sim0.1x.repeat$importance

# Spatial Model
Sim0.2x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim0.2x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim0.2x)

#Get response curves
Sim0.2x.df <- spatialRF::get_response_curves(Sim0.2x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim0.2x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim0.2xcassava_pred.csv")
write.csv(Sim0.2x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim0.2xcassava_resid.csv")
write.csv(Sim0.2x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim0.2xcassava_VarImp.csv")
write.csv(Sim0.2x$performance, file="Crop_GCO/Global/ArtificialYields//Sim0.2xcassava_preformance.csv")
write.csv(Sim0.2x.df, file="Crop_GCO/Global/ArtificialYields//Sim0.2xcassava_curves.csv")

Sim0.2x.repeat <- spatialRF::rf_repeat(
  model = Sim0.2x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)
print(Sim0.2x.repeat)
Sim0.2x.repeat$importance


# Spatial Model
Sim0.5x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim0.5x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim0.5x)

#Get response curves
Sim0.5x.df <- spatialRF::get_response_curves(Sim0.5x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim0.5x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim0.5xcassava_pred.csv")
write.csv(Sim0.5x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim0.5xcassava_resid.csv")
write.csv(Sim0.5x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim0.5xcassava_VarImp.csv")
write.csv(Sim0.5x$performance, file="Crop_GCO/Global/ArtificialYields//Sim0.5xcassava_preformance.csv")
write.csv(Sim0.5x.df, file="Crop_GCO/Global/ArtificialYields//Sim0.5xcassava_curves.csv")

Sim0.5x.repeat <- spatialRF::rf_repeat(
  model = Sim0.5x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)
print(Sim0.5x.repeat)
Sim0.5x.repeat$importance

# Spatial Model
Sim0.8x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim0.8x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim0.8x)

#Get response curves
Sim0.8x.df <- spatialRF::get_response_curves(Sim0.8x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim0.8x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim0.8xcassava_pred.csv")
write.csv(Sim0.8x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim0.8xcassava_resid.csv")
write.csv(Sim0.8x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim0.8xcassava_VarImp.csv")
write.csv(Sim0.8x$performance, file="Crop_GCO/Global/ArtificialYields//Sim0.8xcassava_preformance.csv")
write.csv(Sim0.8x.df, file="Crop_GCO/Global/ArtificialYields//Sim0.8xcassava_curves.csv")

Sim0.8x.repeat <- spatialRF::rf_repeat(
  model = Sim0.8x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)
print(Sim0.8x.repeat)
Sim0.8x.repeat$importance



# Spatial Model
Sim1.0x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim1.0x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim1.0x)

Sim1.0x.repeat <- spatialRF::rf_repeat(
  model = Sim1.0x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)

print(Sim1.0x.repeat)
Sim1.0x.repeat$importance

#Get response curves
Sim1.0x.df <- spatialRF::get_response_curves(Sim1.0x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim1.0x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim1.0xcassava_pred.csv")
write.csv(Sim1.0x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim1.0xcassava_resid.csv")
write.csv(Sim1.0x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim1.0xcassava_VarImp.csv")
write.csv(Sim1.0x$performance, file="Crop_GCO/Global/ArtificialYields//Sim1.0xcassava_preformance.csv")
write.csv(Sim1.0x.df, file="Crop_GCO/Global/ArtificialYields//Sim1.0xcassava_curves.csv")

# Spatial Model
Sim1.1x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim1.1x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim1.1x)

#Get response curves
Sim1.1x.df <- spatialRF::get_response_curves(Sim1.1x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim1.1x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim1.1xcassava_pred.csv")
write.csv(Sim1.1x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim1.1xcassava_resid.csv")
write.csv(Sim1.1x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim1.1xcassava_VarImp.csv")
write.csv(Sim1.1x$performance, file="Crop_GCO/Global/ArtificialYields//Sim1.1xcassava_preformance.csv")
write.csv(Sim1.1x.df, file="Crop_GCO/Global/ArtificialYields//Sim1.1xcassava_curves.csv")

Sim1.1x.repeat <- spatialRF::rf_repeat(
  model = Sim1.1x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)
print(Sim1.1x.repeat)
Sim1.1x.repeat$importance

# Spatial Model
Sim1.2x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim1.2x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim1.2x)

#Get response curves
Sim1.2x.df <- spatialRF::get_response_curves(Sim1.2x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim1.2x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim1.2xcassava_pred.csv")
write.csv(Sim1.2x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim1.2xcassava_resid.csv")
write.csv(Sim1.2x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim1.2xcassava_VarImp.csv")
write.csv(Sim1.2x$performance, file="Crop_GCO/Global/ArtificialYields//Sim1.2xcassava_preformance.csv")
write.csv(Sim1.2x.df, file="Crop_GCO/Global/ArtificialYields//Sim1.2xcassava_curves.csv")

Sim1.2x.repeat <- spatialRF::rf_repeat(
  model = Sim1.2x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)
print(Sim1.2x.repeat)
Sim1.2x.repeat$importance


# Spatial Model
Sim1.5x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim1.5x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 30,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim1.5x)

#Get response curves
Sim1.5x.df <- spatialRF::get_response_curves(Sim1.5x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim1.5x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim1.5xcassava_pred.csv")
write.csv(Sim1.5x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim1.5xcassava_resid.csv")
write.csv(Sim1.5x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim1.5xcassava_VarImp.csv")
write.csv(Sim1.5x$performance, file="Crop_GCO/Global/ArtificialYields//Sim1.5xcassava_preformance.csv")
write.csv(Sim1.5x.df, file="Crop_GCO/Global/ArtificialYields//Sim1.5xcassava_curves.csv")

Sim1.5x.repeat <- spatialRF::rf_repeat(
  model = Sim1.5x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)

print(Sim1.5x.repeat)
Sim1.5x.repeat$importance

# Spatial Model
Sim1.8x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim1.8x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim1.8x)

#Get response curves
Sim1.8x.df <- spatialRF::get_response_curves(Sim1.8x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim1.8x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim1.8xcassava_pred.csv")
write.csv(Sim1.8x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim1.8xcassava_resid.csv")
write.csv(Sim1.8x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim1.8xcassava_VarImp.csv")
write.csv(Sim1.8x$performance, file="Crop_GCO/Global/ArtificialYields//Sim1.8xcassava_preformance.csv")
write.csv(Sim1.8x.df, file="Crop_GCO/Global/ArtificialYields//Sim1.8xcassava_curves.csv")

Sim1.8x.repeat <- spatialRF::rf_repeat(
  model = Sim1.8x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)
print(Sim1.8x.repeat)
Sim1.8x.repeat$importance


# Spatial Model
Sim2.0x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim2.0x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim2.0x)

#Get response curves
Sim2.0x.df <- spatialRF::get_response_curves(Sim2.0x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim2.0x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim2.0xcassava_pred.csv")
write.csv(Sim2.0x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim2.0xcassava_resid.csv")
write.csv(Sim2.0x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim2.0xcassava_VarImp.csv")
write.csv(Sim2.0x$performance, file="Crop_GCO/Global/ArtificialYields//Sim2.0xcassava_preformance.csv")
write.csv(Sim2.0x.df, file="Crop_GCO/Global/ArtificialYields//Sim2.0xcassava_curves.csv")

Sim2.0x.repeat <- spatialRF::rf_repeat(
  model = Sim2.0x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)


# Spatial Model
Sim2.5x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim2.5x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim2.5x)

#Get response curves
Sim2.5x.df <- spatialRF::get_response_curves(Sim2.5x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim2.5x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim2.5xcassava_pred.csv")
write.csv(Sim2.5x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim2.5xcassava_resid.csv")
write.csv(Sim2.5x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim2.5xcassava_VarImp.csv")
write.csv(Sim2.5x$performance, file="Crop_GCO/Global/ArtificialYields//Sim2.5xcassava_preformance.csv")
write.csv(Sim2.5x.df, file="Crop_GCO/Global/ArtificialYields//Sim2.5xcassava_curves.csv")


Sim2.5x.repeat <- spatialRF::rf_repeat(
  model = Sim2.5x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)


# Spatial Model
Sim3x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim3x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster)

gc()
print(Sim3.0x)

#Get response curves
Sim3.0x.df <- spatialRF::get_response_curves(Sim3.0x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim3.0x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim3.0xcassava_pred.csv")
write.csv(Sim3.0x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim3.0xcassava_resid.csv")
write.csv(Sim3.0x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim3.0xcassava_VarImp.csv")
write.csv(Sim3.0x$performance, file="Crop_GCO/Global/ArtificialYields//Sim3.0xcassava_preformance.csv")
write.csv(Sim3.0x.df, file="Crop_GCO/Global/ArtificialYields//Sim3.0xcassava_curves.csv")

Sim3.0x.repeat <- spatialRF::rf_repeat(
  model = Sim3x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)

# Spatial Model
Sim4.0x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim4.0x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim4.0x)

#Get response curves
Sim4.0x.df <- spatialRF::get_response_curves(Sim4.0x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim4.0x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim4.0xcassava_pred.csv")
write.csv(Sim4.0x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim4.0xcassava_resid.csv")
write.csv(Sim4.0x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim4.0xcassava_VarImp.csv")
write.csv(Sim4.0x$performance, file="Crop_GCO/Global/ArtificialYields//Sim4.0xcassava_preformance.csv")
write.csv(Sim4.0x.df, file="Crop_GCO/Global/ArtificialYields//Sim4.0xcassava_curves.csv")

Sim4.0x.repeat <- spatialRF::rf_repeat(
  model = Sim4.0x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)

# Spatial Model
Sim5.0x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim5.0x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava.data,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
gc()
print(Sim5.0x)

Sim5.0x.repeat <- spatialRF::rf_repeat(
  model = Sim5.0x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)

#Get response curves
Sim5.0x.df <- spatialRF::get_response_curves(Sim5.0x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim5.0x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim5.0xcassava_pred.csv")
write.csv(Sim5.0x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim5.0xcassava_resid.csv")
write.csv(Sim5.0x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim5.0xcassava_VarImp.csv")
write.csv(Sim5.0x$performance, file="Crop_GCO/Global/ArtificialYields//Sim5.0xcassava_preformance.csv")
write.csv(Sim5.0x.df, file="Crop_GCO/Global/ArtificialYields//Sim5.0xcassava_curves.csv")




