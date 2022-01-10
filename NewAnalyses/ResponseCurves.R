library(ggplot2)
library(dplyr)
library(ggsci)
library(viridis)


theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank())
#Read in Response Curves
data<-read.csv(file="NewAnalyses/ResponseCurves.csv", header=T)

fertilizer<- data %>% filter(predictor.name == "Fertilizer" & quantile==0.5)
summary(fertilizer)
ggplot(fertilizer, aes(x=predictor,y=response)) +
  facet_wrap(~Crop,scales = "free")  +
  geom_point(aes(color=Region),alpha=0.3) +
  theme_justin +
  xlab("Log10 Fertilizer") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  scale_color_npg() 

pesticide<- data %>% filter(predictor.name == "Pesticide" & quantile==0.5)

ggplot(pesticide, aes(x=predictor,y=response)) +
  facet_wrap(~Crop,scales = "free")  +
  geom_point(aes(color=Region),alpha=0.1) +
  theme_justin +
  xlab("Log10 Pesticide") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  scale_color_npg() 

gdp<- data %>% filter(predictor.name == "GDP" & quantile==0.5)

ggplot(gdp, aes(x=predictor,y=response)) +
  facet_wrap(~Crop,scales = "free")  +
  geom_point(aes(color=Region),alpha=0.1) +
  theme_justin +
  xlab("Log10 GDP") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  scale_color_npg() 

aet<- data %>% filter(predictor.name == "Evapotranspiration" & quantile==0.5)

ggplot(aet, aes(x=predictor,y=response)) +
  facet_wrap(~Crop,scales = "free")  +
  geom_point(aes(color=Region),alpha=0.1) +
  theme_justin +
  xlab("Log 10 Evapotranspiration") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  scale_color_npg() 


######### 

varimp<-read.csv(file="NewAnalyses/VarImp.csv",header=T)

ggplot(varimp, aes(x=Variable, y=Region,fill=Fert.Standard)) + 
  geom_tile() + 
  scale_fill_gradient2(midpoint = 1.0, low="#3B9AB2", high = "#F21A00") +
  facet_wrap(~Crop,nrow = 2,ncol = 6) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  xlab("Driver") + ylab("Global Region") 


######### Model Pref 

pref<-read.csv(file="NewAnalyses/Model_Preformance.csv", header=T)

ggplot(pref, aes(x=Crop,y=Region, fill=pref$r.squared.oob)) + 
  geom_tile() +
  geom_point(aes(size=pref$rmse.oob, alpha=0.5, color="#111111")) +
  theme_bw() +
  
  scale_fill_viridis(option="mako")





