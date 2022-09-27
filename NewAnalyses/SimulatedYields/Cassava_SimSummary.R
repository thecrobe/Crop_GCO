library(tidyverse)
library(cowplot)
data<-read.csv(file="/Users/justinstewart/Dropbox/Collaborations/TobyKiers/Crop_Productivity_GCO/Analysis/Crop_GCO/Cassava_SImRepeat.csv")

gco<-dplyr::filter(data,Variable == "cassavaBinaryGCO")

ggplot(gco, aes(x=gco$Simulation,y=gco$Fert_Standard)) + 
  geom_point() + 
  stat_smooth(method="lm") + 
  xlab("Yield Inflation Factor") + 
  ylab("GCO Variable Importance") + 
  theme_cowplot(12)

mod<-lm(Importance~Simulation, data=gco)
summary(mod)
