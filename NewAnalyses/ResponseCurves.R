library(ggplot2)
library(dplyr)
library(ggsci)
library(viridis)
library(ggpubr)

#Read in Response Curves
data<-read.csv(file="NewAnalyses/ResponseCurves.csv", header=T)

fertilizer<- data %>% filter(predictor.name == "Fertilizer")

ggplot(fertilizer, aes(x=predictor,y=response)) +
  facet_wrap(~Crop,scales = "free",ncol=6)  +
  geom_point(aes(color=Region),alpha=0.3) +
  geom_smooth(method="lm",span=0, color="black", size=2,se = TRUE) +
  theme_justin +
  xlab("Log10 Fertilizer") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  scale_color_npg() 

pesticide<- data %>% filter(predictor.name == "Pesticide")

p_plot<-ggplot(pesticide, aes(x=predictor,y=response)) +
  facet_grid(~Crop,scales = "free_y")  +
  geom_point(aes(color=Region),alpha=0.1) +
  geom_smooth(method="lm", color="black", size=2,se = TRUE) +
  theme_justin +
  xlab("Log10 Pesticide") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  scale_color_npg() 

gdp<- data %>% filter(predictor.name == "GDP")

g_plot<-ggplot(gdp, aes(x=predictor,y=response)) +
  facet_grid(~Crop,scales = "free_y")  +
  geom_point(aes(color=Region),alpha=0.1) +
  geom_smooth(method="lm", color="black", size=2,se = TRUE) +
  theme_justin +
  xlab("Log10 GDP") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  scale_color_npg() 

aet<- data %>% filter(predictor.name == "Evapotranspiration")

a_plot<-ggplot(aet, aes(x=predictor,y=response)) +
  facet_grid(~Crop,scales = "free_y")  +
  geom_point(aes(color=Region),alpha=0.1) +
  geom_smooth(method="lm", color="black", size=2,se = TRUE) +
  theme_justin +
  xlab("Log 10 Evapotranspiration") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  scale_color_npg() 

#Combine Figure
figure <- ggarrange(p_plot,a_plot,g_plot,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

######### 

varimp<-read.csv(file="NewAnalyses/VarImp.csv",header=T)

ggplot(varimp, aes(x=Variable, y=Region,fill=Fert.Standard)) + 
  geom_tile() + 
  scale_fill_gradient2(midpoint = 1.0, low="#3B9AB2", high = "#F21A00") +
  facet_grid(~Crop) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

######### Model Pref 

pref<-read.csv(file="NewAnalyses/Model_Preformance.csv", header=T)

ggplot(pref, aes(x=Crop,y=Region, fill=pref$r.squared.oob)) + 
  geom_tile() +
  geom_point(aes(size=pref$rmse.oob, alpha=0.5, color="#111111")) +
  theme_bw() +
  scale_fill_viridis(option="mako")





