library(tidyverse)
library(cowplot)
data<-read.csv(file="/Users/justinstewart/Downloads/VarImp_Simulation.csv")

gco<-dplyr::filter(data,Variable == "GCO")


ggplot(gco, aes(x=Simulation,y=(FertStandard), color=Variable)) + 
  geom_point(size=2) + 
  stat_smooth(method="loess",se = FALSE) + 
  xlab("Yield Inflation Factor") + 
  ylab("GCO Variable Importance") + 
  theme_cowplot(12)
