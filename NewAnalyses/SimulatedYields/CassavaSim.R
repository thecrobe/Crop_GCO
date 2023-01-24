library(spatialRF) 
library(ggplot2)
library(ranger)
library(Rcpp)
library(parallel)
library(dplyr)
library(cowplot)
#import artificial yields
cassava<-read.csv(file="NewAnalyses/SimulatedYields/ArtificialCassava.csv", header=T)

## cut out the previously inflated
cassava <- cassava[,1:9]
#pick inflation factors
inflations <- c(0.1,0.2,0.5,1.0,1.2,1.5,2.0,2.5,5,10) # fill in the desired levels of outside GCO inflation
# this was really hard to figure out! multiply by inflations only if cassaveBinaryGCO==0
mask_mult <- function(x,y) {((y*x)+1)-(1*y)} # helper function
inflation_array <- t(outer(inflations, cassava$cassavaBinaryGCO==0, FUN = mask_mult))
inflated_values <- log10(inflation_array * cassava$cassava_HgHa+1)
inflated_values_df <- data.frame(inflated_values)
names(inflated_values_df) <- paste0("Sim",inflations,"x")
cassava <- na.omit(cbind(cassava, inflated_values_df))

#coordinates of the cases
x<-cassava$Latitude
y<-cassava$Longitude
xy <- data.frame(cbind(x,y))

#Create distance matrix
distance.matrix <- data.matrix(dist(cbind(cassava$Latitude, cassava$Longitude)))
#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(1,5,10,50)
#random seed for reproducibility
random.seed <- 3828


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
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
)
print(Sim0.1x)

#Get response curves
Sim0.1x.df <- spatialRF::get_response_curves(Sim0.1x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim0.1x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim0.1xcassava_pred.csv")
write.csv(Sim0.1x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim0.1xcassava_resid.csv")
write.csv(Sim0.1x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim0.1xcassava_VarImp.csv")
write.csv(Sim0.1x$performance, file="Crop_GCO/Global/ArtificialYields//Sim0.1xcassava_preformance.csv")
write.csv(Sim0.1x.df, file="Crop_GCO/Global/ArtificialYields//Sim0.1xcassava_curves.csv")


# Spatial Model
Sim0.2x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim0.2x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
)
print(Sim0.2x)

#Get response curves
Sim0.2x.df <- spatialRF::get_response_curves(Sim0.2x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim0.2x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim0.2xcassava_pred.csv")
write.csv(Sim0.2x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim0.2xcassava_resid.csv")
write.csv(Sim0.2x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim0.2xcassava_VarImp.csv")
write.csv(Sim0.2x$performance, file="Crop_GCO/Global/ArtificialYields//Sim0.2xcassava_preformance.csv")
write.csv(Sim0.2x.df, file="Crop_GCO/Global/ArtificialYields//Sim0.2xcassava_curves.csv")


# Spatial Model
Sim0.5x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim0.5x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
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

# Spatial Model
Sim1x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim1x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
)
gc()
print(Sim1x)

#Get response curves
Sim1x.df <- spatialRF::get_response_curves(Sim1x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim1x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim1xcassava_pred.csv")
write.csv(Sim1x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim1xcassava_resid.csv")
write.csv(Sim1x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim1xcassava_VarImp.csv")
write.csv(Sim1x$performance, file="Crop_GCO/Global/ArtificialYields//Sim1xcassava_preformance.csv")
write.csv(Sim1x.df, file="Crop_GCO/Global/ArtificialYields//Sim1xcassava_curves.csv")

# Spatial Model
Sim1.2x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim1.2x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,cluster = local.cluster
)
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
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
)

print(Sim1.5x)

#Get response curves
Sim1.5x.df <- spatialRF::get_response_curves(Sim1.5x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim1.5x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim1.5xcassava_pred.csv")
write.csv(Sim1.5x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim1.5xcassava_resid.csv")
write.csv(Sim1.5x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim1.5xcassava_VarImp.csv")
write.csv(Sim1.5x$performance, file="Crop_GCO/Global/ArtificialYields//Sim1.5xcassava_preformance.csv")
write.csv(Sim1.5x.df, file="Crop_GCO/Global/ArtificialYields//Sim1.5xcassava_curves.csv")


# Spatial Model
Sim2.0x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim2x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
)

print(Sim2.0x)

#Get response curves
Sim2.0x.df <- spatialRF::get_response_curves(Sim2.0x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim2.0x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim2.0xcassava_pred.csv")
write.csv(Sim2.0x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim2.0xcassava_resid.csv")
write.csv(Sim2.0x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim2.0xcassava_VarImp.csv")
write.csv(Sim2.0x$performance, file="Crop_GCO/Global/ArtificialYields//Sim2.0xcassava_preformance.csv")
write.csv(Sim2.0x.df, file="Crop_GCO/Global/ArtificialYields//Sim2.0xcassava_curves.csv")


# Spatial Model
Sim2.5x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim2.5x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
)
print(Sim2.5x)

#Get response curves
Sim2.5x.df <- spatialRF::get_response_curves(Sim2.5x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim2.5x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim2.5xcassava_pred.csv")
write.csv(Sim2.5x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim2.5xcassava_resid.csv")
write.csv(Sim2.5x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim2.5xcassava_VarImp.csv")
write.csv(Sim2.5x$performance, file="Crop_GCO/Global/ArtificialYields//Sim2.5xcassava_preformance.csv")
write.csv(Sim2.5x.df, file="Crop_GCO/Global/ArtificialYields//Sim2.5xcassava_curves.csv")


# Spatial Model
Sim5x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim5x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
)
print(Sim5x)

#Get response curves
Sim5x.df <- spatialRF::get_response_curves(Sim5x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim5x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim5xcassava_pred.csv")
write.csv(Sim5x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim5xcassava_resid.csv")
write.csv(Sim5x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim5xcassava_VarImp.csv")
write.csv(Sim5x$performance, file="Crop_GCO/Global/ArtificialYields//Sim5xcassava_preformance.csv")
write.csv(Sim5x.df, file="Crop_GCO/Global/ArtificialYields//Sim5xcassava_curves.csv")



# Spatial Model
Sim10x <- spatialRF::rf_spatial(
  dependent.variable.name = ("Sim10x"), 
  predictor.variable.names = c("AET_mean", "cassava_Fertilizer","Pesticide","GDP_Mean","cassavaBinaryGCO"),
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10,
  xy = xy,
  data = cassava,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  seed = random.seed,
  cluster = local.cluster
)
print(Sim10x)

#Get response curves
Sim10x.df <- spatialRF::get_response_curves(Sim10x,variables = c("AET_mean","cassava_Fertilizer", "Pesticide", "GDP_Mean","cassavaBinaryGCO"))
write.csv(Sim10x$predictions, file="Crop_GCO/Global/ArtificialYields/Sim10xcassava_pred.csv")
write.csv(Sim10x$residuals$values, file="Crop_GCO/Global/ArtificialYields//Sim10xcassava_resid.csv")
write.csv(Sim10x$variable.importance, file="Crop_GCO/Global/ArtificialYields//Sim10xcassava_VarImp.csv")
write.csv(Sim10x$performance, file="Crop_GCO/Global/ArtificialYields//Sim10xcassava_preformance.csv")
write.csv(Sim10x.df, file="Crop_GCO/Global/ArtificialYields//Sim10xcassava_curves.csv")


#cleanup and shut down cluster
gc()
stopCluster(local.cluster)

#### Import summarized results 
data<-read.csv(file="SimulatedYields/ArtificialYields/Cassava_SimSummary_data.csv")
gco<-dplyr::filter(data,Variable == "GCO") #subset GCO variable
ggplot(gco, aes(x=log10(Inflation.Factor),y=(Standardized.Importance))) + 
  geom_point(size=2) + 
  geom_hline(yintercept = 1) +
  stat_smooth(method="loess",se = FALSE) + 
  xlab("Log10 Yield Inflation Factor") + 
  ylab("GCO Variable Importance") + 
  theme_cowplot(12)

