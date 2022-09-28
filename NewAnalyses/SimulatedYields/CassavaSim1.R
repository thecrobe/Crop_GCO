library(spatialRF) 
library(ggplot2)
library(ranger)
library(Rcpp)
library(parallel)
library(dplyr)
library(cowplot)

inflations <- c(0.1,0.5,1.0,1.5,2.0,2.5,3,3.5,4.0,4.5,5.0,10.0,15,25) # fill in the desired levels of outside GCO inflation

#import artificial yields
cassava<-read.csv(file="NewAnalyses/SimulatedYields/ArtificialCassava.csv", header=T)

## cut out the previously inflated

cassava <- cassava[,1:9]

# this was really hard to figure out! multiply by inflations only if cassaveBinaryGCO==0
mask_mult <- function(x,y) {((y*x)+1)-(1*y)} # helper function

inflation_array <- t(outer(inflations, cassava$cassavaBinaryGCO==0, FUN = mask_mult))
inflated_values <- inflation_array * cassava$cassava_HgHa

inflated_values_df <- data.frame(inflated_values)
names(inflated_values_df) <- paste0("Sim",inflations,"x")

cassava <- cbind(cassava, inflated_values_df)

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



### I would suggest rewriting this with apply or similar to avoid repeating the same code block for each level of inflation. I will show how next time.

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
Sim1.0x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim1x"), 
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
Sim2.0x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim2x"), 
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
Sim3.5x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim3.5x"), 
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
print(Sim3.5x)

#Get response curves
Sim3.5x.df <- spatialRF::get_response_curves(Sim3.5x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim3.5x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim3.5xcassava_pred.csv")
write.csv(Sim3.5x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim3.5xcassava_resid.csv")
write.csv(Sim3.5x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim3.5xcassava_VarImp.csv")
write.csv(Sim3.5x$performance, file="Crop_GCO/Global/ArtificialYields//Sim3.5xcassava_preformance.csv")
write.csv(Sim3.5x.df, file="Crop_GCO/Global/ArtificialYields//Sim3.5xcassava_curves.csv")

Sim3.5x.repeat <- spatialRF::rf_repeat(
  model = Sim3.5x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)

# Spatial Model
Sim4.0x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim4x"), 
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
Sim4.5x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim4.5x"), 
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
print(Sim4.5x)

#Get response curves
Sim4.5x.df <- spatialRF::get_response_curves(Sim4.5x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim4.5x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim4.5xcassava_pred.csv")
write.csv(Sim4.5x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim4.5xcassava_resid.csv")
write.csv(Sim4.5x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim4.5xcassava_VarImp.csv")
write.csv(Sim4.5x$performance, file="Crop_GCO/Global/ArtificialYields//Sim4.5xcassava_preformance.csv")
write.csv(Sim4.5x.df, file="Crop_GCO/Global/ArtificialYields//Sim4.5xcassava_curves.csv")

Sim4.5x.repeat <- spatialRF::rf_repeat(
  model = Sim4.5x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)


# Spatial Model
Sim5.0x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim5x"), 
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


# Spatial Model

Sim10x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim10x"), 
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
print(Sim10x)

#Get response curves
Sim10x.df <- spatialRF::get_response_curves(Sim10x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim10x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim10xcassava_pred.csv")
write.csv(Sim10x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim10xcassava_resid.csv")
write.csv(Sim10x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim10xcassava_VarImp.csv")
write.csv(Sim10x$performance, file="Crop_GCO/Global/ArtificialYields//Sim10xcassava_preformance.csv")
write.csv(Sim10x.df, file="Crop_GCO/Global/ArtificialYields//Sim10xcassava_curves.csv")

Sim10.0x.repeat <- spatialRF::rf_repeat(
  model = Sim10x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)

# Spatial Model
Sim50x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim25x"), 
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
print(Sim50x)

#Get response curves
Sim50x.df <- spatialRF::get_response_curves(Sim50x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim50x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim50xcassava_pred.csv")
write.csv(Sim50x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim50xcassava_resid.csv")
write.csv(Sim50x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim50xcassava_VarImp.csv")
write.csv(Sim50x$performance, file="Crop_GCO/Global/ArtificialYields//Sim50xcassava_preformance.csv")
write.csv(Sim50x.df, file="Crop_GCO/Global/ArtificialYields//Sim50xcassava_curves.csv")

Sim50x.repeat <- spatialRF::rf_repeat(
  model = Sim50x, 
  repetitions = 5,
  seed = random.seed,
  verbose = FALSE,
  cluster = local.cluster
)


#Stop cluster
parallel::stopCluster(cl = local.cluster) #stop cluster

# Summarize results 
summary<-read.csv(file="NewAnalyses/SimulatedYields/Cassava_SimRepeat.csv")

gco<-dplyr::filter(summary,Variable == "cassavaBinaryGCO")

ggplot(summary, aes(x=Simulation,y=Fert_Standard)) + 
  geom_point() + 
  stat_smooth(method="lm") + 
  xlab("Yield Inflation Factor") + 
  ylab("GCO Variable Importance") + 
  theme_cowplot(12)

mod<-lm(Importance~Simulation, data=gco)
plot(mod) #check fit
summary(mod)

