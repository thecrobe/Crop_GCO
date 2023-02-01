library(dplyr)
library(ggplot2)
library(ggpubr)

#Graphics
theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + theme(legend.position = "none")

mapping<-read.csv(file="NewAnalyses/SpatialRandomForests/GCO_Mapping.csv", header=T)

#Barley
barley<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/barley_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/barley_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
barley$residuals<-10^(barley$x)-1 #backtransform

barley.plot<-ggplot(barley, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Barley") + 
  labs(color = "GCO")


#Cassava
cassava<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/cassava_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/cassava_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
cassava$residuals<-10^(cassava$x)-1 #backtransform

cassava.plot<-ggplot(cassava, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Cassava") + 
  labs(color = "GCO")

#Groundnut
groundnut<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/groundnut_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/groundnut_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
groundnut$residuals<-10^(groundnut$x)-1 #backtransform

groundnut.plot<-ggplot(groundnut, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Groundnut") + 
  labs(color = "GCO")

#Maize
maize<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/maize_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/maize_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
maize$residuals<-10^(maize$x)-1 #backtransform

maize.plot<-ggplot(maize, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Maize") + 
  labs(color = "GCO")

#Rapeseed
rapeseed<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/rapeseed_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/rapeseed_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
rapeseed$residuals<-10^(rapeseed$x)-1 #backtransform

rapeseed.plot<-ggplot(rapeseed, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Rapeseed") + 
  labs(color = "GCO")

#Rice
rice<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/rice_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/rice_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
rice$residuals<-10^(rice$x)-1 #backtransform

rice.plot<-ggplot(rice, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Rice") + 
  labs(color = "GCO")

#Rye
rye<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/rye_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/rye_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
rye$residuals<-10^(rye$x)-1 #backtransform

rye.plot<-ggplot(rye, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Rye") + 
  labs(color = "GCO")

#Sorghum
sorghum<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/sorghum_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/sorghum_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
sorghum$residuals<-10^(sorghum$x)-1 #backtransform

sorghum.plot<-ggplot(sorghum, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Sorghum") + 
  labs(color = "GCO")

#Soybean
soybean<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/soybean_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/soybean_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
soybean$residuals<-10^(soybean$x)-1 #backtransform

soybean.plot<-ggplot(soybean, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Soybean") + 
  labs(color = "GCO")

#Sunflower
sunflower<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/sunflower_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/sunflower_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
sunflower$residuals<-10^(sunflower$x)-1 #backtransform

sunflower.plot<-ggplot(sunflower, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
  xlab("Residuals") +
  scale_color_npg() + ggtitle("Sunflower") + 
  labs(color = "GCO")

#Wheat
wheat<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/wheat_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/wheat_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
wheat$residuals<-10^(wheat$x)-1 #backtransform

wheat.plot<-ggplot(wheat, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
   xlab("Residuals") +
  scale_color_npg() + ggtitle("Wheat") + 
  labs(color = "GCO")

#Potato
potato<-read.csv(file="NewAnalyses/SpatialRandomForests/ModelOutputs/potato_resid.csv",header=T)
noGCO<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/potato_nogco_resid")
noGCO$residuals<-10^(noGCO$x)-1
potato$residuals<-10^(potato$x)-1 #backtransform

potato.plot<-ggplot(potato, aes(x=(residuals))) + 
  geom_density(size=1, color="blue") + 
  geom_density(size=1,data = noGCO, aes(x=residuals), color="red") +
  theme_justin + 
  xlab("Residuals") +
  scale_color_npg() + ggtitle("Potato") + 
  labs(color = "GCO")


### Plot
multipanel<-ggarrange(barley.plot,cassava.plot,groundnut.plot,maize.plot,potato.plot,rapeseed.plot,rice.plot,rye.plot,sorghum.plot,soybean.plot,sunflower.plot,wheat.plot,ncol = 4,nrow = 3,common.legend = TRUE)

ggsave(plot = multipanel,filename = "Residuals_GCO_NoGCO.tiff",path = "../Figures/",dpi = 300,height = 8,width = 10)















