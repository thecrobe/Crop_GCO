library(ggplot2) #plots 
library(ggpubr) #plots
library(dplyr) #data wrangling
library(viridis) #colors

#Set graphic theme
theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank())


# Model variable importance
varimp<-read.csv(file="NewAnalyses/Global_VarImp.csv",header=T)
varimp$Variable <- factor(varimp$Variable, levels=c("Fertilizer", "Evapotranspiration", "Pesticide","GDP","GCO"))
varimp$Crop<-factor(varimp$Crop, levels=unique(varimp$Crop))

ggplot(varimp, aes(x=Variable, y=Crop,fill=FertStandardized)) +
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_fill_viridis_c()+
  xlab("Driver") + ylab("Crop") + 
  labs(fill="Variable Importance")




#Read in Response Curves
data<-read.csv(file="NewAnalyses/Global_ResponseCurves.csv", header=T)

fertilizer<- data %>% filter(predictor.name == "Fertilizer" & quantile==0.5)
fertilizer.plot<-ggplot(fertilizer, aes(x=predictor,y=response)) +
  facet_wrap(~Crop,scales = "free")  +
  geom_point() +
  theme_justin +
  xlab("Log10 Fertilizer") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  labs(title = "Fertilizer Response Curves")

pesticide<- data %>% filter(predictor.name == "Pesticide" & quantile==0.5)
pesticide.plot<-ggplot(pesticide, aes(x=predictor,y=response)) +
  facet_wrap(~Crop,scales = "free")  +
  geom_point() +
  theme_justin +
  xlab("Log10 Pesticide") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  labs(title = "Pesticide Response Curves")

gdp<- data %>% filter(predictor.name == "GDP" & quantile==0.5)
gdp.plot<-ggplot(gdp, aes(x=predictor,y=response)) +
  facet_wrap(~Crop,scales = "free")  +
  geom_point() +
  theme_justin +
  xlab("Log10 GDP") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  labs(title = "GDP Response Curves")


aet<- data %>% filter(predictor.name == "Evapotranspiration" & quantile==0.5)

aet.plot<-ggplot(aet, aes(x=predictor,y=response)) +
  facet_wrap(~Crop,scales = "free")  +
  geom_point() +
  theme_justin +
  xlab("Log10 Evapotranspiration") +
  ylab ("Log10 Yield") +
  labs(color="Predictor") +
  labs(title = "Evapotranspiration Response Curves")


ggsave(plot = gdp.plot,
       dpi = 300,
       path = "../Figures/",
       filename = "GDP_ResponseCurve.tiff"
       ,width = 8,height = 5)
