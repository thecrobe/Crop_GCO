library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

#Graphics
theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + theme(legend.position = "none")

mapping<-read.csv(file="Global/GCO_Mapping.csv", header=T)

#barley
barley<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/barley_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/barley_nogco_resid",header=T)
barley$residuals<-10^(residuals$barley.spatial.residuals.values)-1
barley<-left_join(barley,mapping)

barley.plot<-ggplot(barley, aes(x=(residuals),color=BarleyGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Barley") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Barley") + 
  labs(color = "GCO")

#cassava
cassava<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/cassava_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/cassava_nogco_resid",header=T)
cassava$residuals<-10^(residuals$x)-1
cassava<-cbind(cassava,mapping)
join<-inner_join(cassava, mapping,by="FISHNET_ID")

cassava.plot<-ggplot(join, aes(x=(residuals),color=CassavaGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Cassava") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Cassava") + 
  labs(color = "GCO")

#groundnut
groundnut<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/groundnut_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/groundnut_nogco_resid",header=T)
groundnut$residuals<-10^(residuals$x)-1
groundnut<-left_join(groundnut,mapping)

groundnut.plot<-ggplot(groundnut, aes(x=(residuals),color=GroundnutGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Groundnut") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Groundnut") + 
  labs(color = "GCO")

#maize
maize<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/maize_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/maize_nogco_resid",header=T)
maize$residuals<-10^(residuals$x)-1
maize<-left_join(maize,mapping)

maize.plot<-ggplot(maize, aes(x=(residuals),color=MaizeGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Maize") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Maize") + 
  labs(color = "GCO") 

#potato
potato<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/potato_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/potato_nogco_resid",header=T)
potato$residuals<-10^(residuals$x)-1
potato<-left_join(potato,mapping)

potato.plot<-ggplot(potato, aes(x=(residuals),color=PotatoGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Potato") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Potato") + 
  labs(color = "GCO")

#rapeseed
rapeseed<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/rapeseed_df.csv"))
residuals<-read.csv(file="NewAnalyses//SpatialRandomForests_NoGCO/rapeseed_nogco_resid",header=T)
rapeseed$residuals<-10^(residuals$x)-1
rapeseed<-left_join(rapeseed,mapping)

rapeseed.plot<-ggplot(rapeseed, aes(x=(residuals),color=RapeseedGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Rapeseed") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Rapeseed") + 
  labs(color = "GCO")

#rice
rice<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/rice_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/rice_nogco_resid",header=T)
rice$residuals<-10^(residuals$x)-1
rice<-left_join(rice,mapping)

rice.plot<-ggplot(rice, aes(x=(residuals),color=RiceGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Rice") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Rice") + 
  labs(color = "GCO")

#rye
rye<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/rye_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/rye_nogco_resid.csv",header=T)
rye$residuals<-10^(residuals$x)-1
rye<-left_join(rye,mapping,by="FISHNET_ID")

rye.plot<-ggplot(rye, aes(x=(residuals),color=RyeGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Rye") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Rye") + 
  labs(color = "GCO")

#sorghum
sorghum<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/sorghum_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/sorghum_nogco_resid",header=T)
sorghum$residuals<-10^(residuals$x)-1
sorghum<-left_join(sorghum,mapping)

sorghum.plot<-ggplot(sorghum, aes(x=(residuals),color=SorghumGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Sorghum") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Sorghum") + 
  labs(color = "GCO")

#soybean
soybean<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/soybean_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/soybean_nogco_resid",header=T)
soybean$residuals<-10^(residuals$x)-1
soybean<-left_join(soybean,mapping)

soybean.plot<-ggplot(soybean, aes(x=(residuals),color=SoybeanGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Soybean") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Soybean") + 
  labs(color = "GCO")

#sunflower
sunflower<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/sunflower_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/sunflower_nogco_resid",header=T)
sunflower$residuals<-10^(residuals$x)-1
sunflower<-left_join(sunflower,mapping)

sunflower.plot<-ggplot(sunflower, aes(x=(residuals),color=SunflowerGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Sunflower") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Sunflower") + 
  labs(color = "GCO")

#wheat
wheat<-na.omit(read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/wheat_df.csv"))
residuals<-read.csv(file="NewAnalyses/SpatialRandomForests_NoGCO/wheat_nogco_resid",header=T)
wheat$residuals<-10^(residuals$x)-1
wheat<-left_join(wheat,mapping)

wheat.plot<-ggplot(wheat, aes(x=(residuals),color=WheatGCO)) + 
  geom_density(size=1, linetype="dashed") + 
  theme_justin + 
  xlab("Wheat") + xlab("Prediction Error") +
  scale_color_npg() + ggtitle("Wheat") + 
  labs(color = "GCO")

multipanel<-ggarrange(barley.plot,cassava.plot,groundnut.plot,maize.plot,potato.plot,rapeseed.plot,rice.plot,rye.plot,sorghum.plot,soybean.plot,sunflower.plot,wheat.plot,ncol = 4,nrow = 3,common.legend = TRUE)

ggsave(plot = multipanel,filename = "PredictionError.png",path = "../Figures/",dpi = 300,height = 8,width = 10)















