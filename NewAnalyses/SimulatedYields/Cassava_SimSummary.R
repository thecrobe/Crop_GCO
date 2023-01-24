library(tidyverse)
library(cowplot)
library(ggpubr)

#Import data
data<-read.csv(file="NewAnalyses/SimulatedYields/ArtificialYields/Cassava_SimSummary_data.csv")
data<-read.csv(file="NewAnalyses/SimulatedYields/ArtificialYields/Cassava_SimSummary_data.csv")

#Subset GCO 
gco<-dplyr::filter(data,Variable == "GCO")

#Variable Importance
ggplot(gco, aes(x=log10(Inflation.Factor),y=log10((Standardized.Importance)))) + # Plot standard to fertilizer
  stat_smooth(method="loess", color="grey", size=1.5) + 
  #geom_vline(xintercept = log10(1.08), color="#E47719", size=1.5) + # GDP
  #geom_vline(xintercept = log10(1.02), color="#0090D1", size=1.5) + # AET
  #geom_vline(xintercept = log10(3), color="#0090D1", size=1.5) + # Pesticide
  geom_point(size=3) +
  annotate(xmin = 0, xmax = log10(1.05), 
           ymin = -Inf, ymax = Inf, 
           geom = "rect", alpha = 0.5, fill="red") +
  xlab("Log10 Yield Inflation Factor") + 
  ylab("Log10 GCO Importance Relative To Fertilizer")+ 
  scale_color_jco() +
  theme_cowplot(12)

#Variable Importance
ggplot(gco, aes(x=log10(Inflation.Factor),y=log10((Standardized.Importance)))) + # Plot standard to fertilizer
  annotate(xmin = 0, xmax = log10(1.5), 
           ymin = -Inf, ymax = Inf, 
           geom = "rect", alpha = 0.4, fill="#0090D1") + 
  geom_vline(xintercept = log10(1.08), color="black", size=1.5) + # GDP
  annotate("text", x = .12, y = .5, label = "GDP", color="black",fontface="bold") +
  geom_vline(xintercept = log10(1.02), color="black", size=1.5) + # AET
  annotate("text", x = -.09, y = .3 , label = "AET",fontface="bold") +
  geom_vline(xintercept = log10(3), color="black", size=1.5) + # Pesticide
  annotate("text", x = .63, y = .7 , label = "Pesticide",fontface="bold") +
  xlab("Log10 Yield Inflation Factor") + 
  ylab("Log10 GCO Importance Relative To Fertilizer")+ 
  stat_smooth(method="gam", color="grey", size=1.5) + 
  geom_point(size=3) +
  scale_color_jco() +
  theme_cowplot(12)
