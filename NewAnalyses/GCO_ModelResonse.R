library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
setwd(dir = "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/Crop_GCO/NewAnalyses/")

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#Graphics
theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + theme(legend.position = "none")

mapping<-read.csv(file="../Fishnets/GCO_Mapping.csv", header=T)

#barley
barley.real<-na.omit(read.csv(file="../Global/Models/Barley_RF.csv"))
barley.real$FISHNET_ID<-barley.real$Fishnet_ID
   
barley.pred<-read.csv(file="Global/ModelOutputs/barley_pred.csv",header=T)
barley.real$Predictions<-barley.pred$values
barley.real$PredOverReal<-barley.real$Predictions/barley.real$barley_HgHa
barley<-left_join(barley.real,mapping)

barley.plot<-ggplot(barley, aes(x=log10(PredOverReal),color=BarleyGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Barley") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Barley") + 
  labs(color = "GCO")


#cassava
cassava.real<-na.omit(read.csv(file="../Global/Models/Cassava_RF.csv"))
cassava.real$FISHNET_ID<-cassava.real$Fishnet_ID

cassava.pred<-read.csv(file="Global/ModelOutputs/cassava_pred.csv",header=T)
cassava.real$Predictions<-cassava.pred$values
cassava.real$PredOverReal<-cassava.real$Predictions/cassava.real$cassava_HgHa
cassava<-left_join(cassava.real,mapping)

cassava.plot<-ggplot(cassava, aes(x=log10(PredOverReal),color=CassavaGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Cassava") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Cassava")

#groundnut
groundnut.real<-na.omit(read.csv(file="../Global/Models/Groundnut_RF.csv"))
groundnut.real$FISHNET_ID<-groundnut.real$Fishnet_ID

groundnut.pred<-read.csv(file="Global/ModelOutputs/groundnut_pred.csv",header=T)
groundnut.real$Predictions<-groundnut.pred$values
groundnut.real$PredOverReal<-groundnut.real$Predictions/groundnut.real$groundnut_HgHa
groundnut<-left_join(groundnut.real,mapping)

groundnut.plot<-ggplot(groundnut, aes(x=log10(PredOverReal),color=GroundnutGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Groundnut") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Groundnut")


#maize
maize.real<-na.omit(read.csv(file="../Global/Models/Maize_RF.csv"))
maize.real$FISHNET_ID<-maize.real$Fishnet_ID

maize.pred<-read.csv(file="Global/ModelOutputs/maize_pred.csv",header=T)
maize.real$Predictions<-maize.pred$values
maize.real$PredOverReal<-maize.real$Predictions/maize.real$maize_HgHa
maize<-left_join(maize.real,mapping)

maize.plot<-ggplot(maize, aes(x=log10(PredOverReal),color=MaizeGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Maize") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Maize")

#potato
potato.real<-na.omit(read.csv(file="../Global/Models/Potato_RF.csv"))
potato.real$FISHNET_ID<-potato.real$Fishnet_ID

potato.pred<-read.csv(file="Global/ModelOutputs/potato_pred.csv",header=T)
potato.real$Predictions<-potato.pred$values
potato.real$PredOverReal<-potato.real$Predictions/potato.real$mean_potat
potato<-left_join(potato.real,mapping)

potato.plot<-ggplot(potato, aes(x=log10(PredOverReal),color=PotatoGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Potato") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Potato")


#rapeseed
rapeseed.real<-na.omit(read.csv(file="../Global/Models/Rapeseed_RF.csv"))
rapeseed.real$FISHNET_ID<-rapeseed.real$Fishnet_ID

rapeseed.pred<-read.csv(file="Global/ModelOutputs/rapeseed_pred.csv",header=T)
rapeseed.real$Predictions<-rapeseed.pred$values
rapeseed.real$PredOverReal<-rapeseed.real$Predictions/rapeseed.real$rapeseed_HgHa
rapeseed<-left_join(rapeseed.real,mapping)

rapeseed.plot<-ggplot(rapeseed, aes(x=log10(PredOverReal),color=RapeseedGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Rapeseed") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Rapeseed")


#Rice
Rice.real<-na.omit(read.csv(file="../Global/Models/Rice_RF.csv"))
Rice.real$FISHNET_ID<-Rice.real$Fishnet_ID

Rice.pred<-read.csv(file="Global/ModelOutputs/Rice_pred.csv",header=T)
Rice.real$Predictions<-Rice.pred$values
Rice.real$PredOverReal<-Rice.real$Predictions/Rice.real$rice_HgHa
Rice<-left_join(Rice.real,mapping)

rice.plot<-ggplot(Rice, aes(x=log10(PredOverReal),color=RiceGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Rice") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Rice")


#Rye
Rye.real<-na.omit(read.csv(file="../Global/Models/Rye_RF.csv"))
Rye.real$FISHNET_ID<-Rye.real$Fishnet_ID

Rye.pred<-read.csv(file="Global/ModelOutputs/Rye_pred.csv",header=T)
Rye.real$Predictions<-Rye.pred$values
Rye.real$PredOverReal<-Rye.real$Predictions/Rye.real$rye_HgHa
Rye<-left_join(Rye.real,mapping)

rye.plot<-ggplot(Rye, aes(x=log10(PredOverReal),color=RyeGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Rye") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Rye")


#Sorghum
Sorghum.real<-na.omit(read.csv(file="../Global/Models/Sorghum_RF.csv"))
Sorghum.real$FISHNET_ID<-Sorghum.real$Fishnet_ID

Sorghum.pred<-read.csv(file="Global/ModelOutputs/Sorghum_pred.csv",header=T)
Sorghum.real$Predictions<-Sorghum.pred$values
Sorghum.real$PredOverReal<-Sorghum.real$Predictions/Sorghum.real$Yield
Sorghum<-left_join(Sorghum.real,mapping)

sorghum.plot<-ggplot(Sorghum, aes(x=log10(PredOverReal),color=SorghumGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Sorghum") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Sorghum")


#Soybean
Soybean.real<-na.omit(read.csv(file="../Global/Models/Soybean_RF.csv"))
Soybean.real$FISHNET_ID<-Soybean.real$Fishnet_ID

Soybean.pred<-read.csv(file="Global/ModelOutputs/Soybean_pred.csv",header=T)
Soybean.real$Predictions<-Soybean.pred$values
Soybean.real$PredOverReal<-Soybean.real$Predictions/Soybean.real$soybean_HgHa
Soybean<-left_join(Soybean.real,mapping)

soybean.plot<-ggplot(Soybean, aes(x=log10(PredOverReal),color=SoybeanGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Soybean") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Soybean")


#Sunflower
Sunflower.real<-na.omit(read.csv(file="../Global/Models//Sunflower_RF.csv"))
Sunflower.real$FISHNET_ID<-Sunflower.real$Fishnet_ID

Sunflower.pred<-read.csv(file="Global/ModelOutputs/Sunflower_pred.csv",header=T)
Sunflower.real$Predictions<-Sunflower.pred$values
Sunflower.real$PredOverReal<-Sunflower.real$Predictions/Sunflower.real$sunflower_HgHa
Sunflower<-left_join(Sunflower.real,mapping)

sunflower.plot<-ggplot(Sunflower, aes(x=log10(PredOverReal),color=SunflowerGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Sunflower") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Sunflower")


#Wheat
Wheat.real<-na.omit(read.csv(file="../Global/Models//Wheat_RF.csv"))
Wheat.real$FISHNET_ID<-Wheat.real$Fishnet_ID

Wheat.pred<-read.csv(file="Global/ModelOutputs/Wheat_pred.csv",header=T)
Wheat.real$Predictions<-Wheat.pred$values
Wheat.real$PredOverReal<-Wheat.real$Predictions/Wheat.real$wheat_HgHa
Wheat<-left_join(Wheat.real,mapping)

wheat.plot<-ggplot(Wheat, aes(x=(PredOverReal),color=WheatGCO)) + 
  geom_density(size=2) + 
  geom_vline(xintercept = 0,linetype="dotdash", size=1) +
  theme_justin + 
  xlab("Wheat") + xlab("Log10 Prediction Error") +
  scale_color_npg() + ggtitle("Wheat")


multipanel<-ggarrange(barley.plot,cassava.plot,groundnut.plot,maize.plot,potato.plot,rapeseed.plot,rice.plot,rye.plot,sorghum.plot,soybean.plot,sunflower.plot,wheat.plot,ncol = 4,nrow = 3,common.legend = TRUE)

ggsave(plot = multipanel,filename = "PredictionError.png",path = "../Figures/",dpi = 300,height = 8,width = 10)















