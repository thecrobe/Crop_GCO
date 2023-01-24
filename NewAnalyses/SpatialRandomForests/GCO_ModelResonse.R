library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

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
barley<-na.omit(read.csv(file="Global/Models/Barley_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/barley_resid.csv",header=T)
barley$residuals<-10^(residuals$x)-1
barley<-left_join(barley,mapping)

barley.plot<-ggplot(barley, aes(x=(residuals),color=BarleyGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Barley") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Barley") + 
  labs(color = "GCO")

#cassava
cassava<-na.omit(read.csv(file="Global/Models/Cassava_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/cassava_resid.csv",header=T)
cassava$residuals<-10^(residuals$x)-1
cassava<-left_join(cassava,mapping)

cassava.plot<-ggplot(cassava, aes(x=(residuals),color=CassavaGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Cassava") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Cassava") + 
  labs(color = "GCO")

#groundnut
groundnut<-na.omit(read.csv(file="Global/Models/Groundnut_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/groundnut_resid.csv",header=T)
groundnut$residuals<-10^(residuals$x)-1
groundnut<-left_join(groundnut,mapping)

groundnut.plot<-ggplot(groundnut, aes(x=(residuals),color=GroundnutGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Groundnut") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Groundnut") + 
  labs(color = "GCO")

#maize
maize<-na.omit(read.csv(file="Global/Models/Maize_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/maize_resid.csv",header=T)
maize$residuals<-10^(residuals$x)-1
maize<-left_join(maize,mapping)

maize.plot<-ggplot(maize, aes(x=(residuals),color=MaizeGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Maize") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Maize") + 
  labs(color = "GCO")

#potato
potato<-na.omit(read.csv(file="Global/Models/Potato_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/potato_resid.csv",header=T)
potato$residuals<-10^(residuals$x)-1
potato<-left_join(potato,mapping)

potato.plot<-ggplot(potato, aes(x=(residuals),color=PotatoGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Potato") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Potato") + 
  labs(color = "GCO")

#rapeseed
rapeseed<-na.omit(read.csv(file="Global/Models/Rapeseed_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/rapeseed_resid.csv",header=T)
rapeseed$residuals<-10^(residuals$x)-1
rapeseed<-left_join(rapeseed,mapping)

rapeseed.plot<-ggplot(rapeseed, aes(x=(residuals),color=RapeseedGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Rapeseed") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Rapeseed") + 
  labs(color = "GCO")

#rice
rice<-na.omit(read.csv(file="Global/Models/Rice_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/rice_resid.csv",header=T)
rice$residuals<-10^(residuals$x)-1
rice<-left_join(rice,mapping)

rice.plot<-ggplot(rice, aes(x=(residuals),color=RiceGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Rice") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Rice") + 
  labs(color = "GCO")

#rye
rye<-na.omit(read.csv(file="Global/Models/Rye_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/rye_resid.csv",header=T)
rye$residuals<-10^(residuals$x)-1
rye<-left_join(rye,mapping)

rye.plot<-ggplot(rye, aes(x=(residuals),color=RyeGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Rye") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Rye") + 
  labs(color = "GCO")

#sorghum
sorghum<-na.omit(read.csv(file="Global/Models/Sorghum_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/sorghum_resid.csv",header=T)
sorghum$residuals<-10^(residuals$x)-1
sorghum<-left_join(sorghum,mapping)

sorghum.plot<-ggplot(sorghum, aes(x=(residuals),color=SorghumGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Sorghum") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Sorghum") + 
  labs(color = "GCO")

#soybean
soybean<-na.omit(read.csv(file="Global/Models/Soybean_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/soybean_resid.csv",header=T)
soybean$residuals<-10^(residuals$x)-1
soybean<-left_join(soybean,mapping)

soybean.plot<-ggplot(soybean, aes(x=(residuals),color=SoybeanGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Soybean") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Soybean") + 
  labs(color = "GCO")

#sunflower
sunflower<-na.omit(read.csv(file="Global/Models/Sunflower_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/sunflower_resid.csv",header=T)
sunflower$residuals<-10^(residuals$x)-1
sunflower<-left_join(sunflower,mapping)

sunflower.plot<-ggplot(sunflower, aes(x=(residuals),color=SunflowerGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Sunflower") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Sunflower") + 
  labs(color = "GCO")

#wheat
wheat<-na.omit(read.csv(file="Global/Models/Wheat_RF.csv"))
residuals<-read.csv(file="NewAnalyses/Global/ModelOutputs/wheat_resid.csv",header=T)
wheat$residuals<-10^(residuals$x)-1
wheat<-left_join(wheat,mapping)

wheat.plot<-ggplot(wheat, aes(x=(residuals),color=WheatGCO)) + 
  geom_density(size=2) + 
  theme_justin + 
  xlab("Wheat") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Wheat") + 
  labs(color = "GCO")

multipanel<-ggarrange(barley.plot,cassava.plot,groundnut.plot,maize.plot,potato.plot,rapeseed.plot,rice.plot,rye.plot,sorghum.plot,soybean.plot,sunflower.plot,wheat.plot,ncol = 4,nrow = 3,common.legend = TRUE)

ggsave(plot = multipanel,filename = "PredictionError.png",path = "../Figures/",dpi = 300,height = 8,width = 10)















