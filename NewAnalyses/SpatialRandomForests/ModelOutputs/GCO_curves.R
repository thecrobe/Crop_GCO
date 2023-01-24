#data wrangling
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)

#Custom Function to plot
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


#####Barley 
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
barley<-na.omit(yield.mapping %>% dplyr::select(barley_HgHa,BarleyGCO))

#Inside vs Outside GCO 
barley.inside<-filter(barley, BarleyGCO == "Inside")
barley.outside<-filter(barley, BarleyGCO == "Outside")
barley.m<-melt(barley)
summary(barley.inside)
summary(barley.outside)

#Plot
Barley<-ggplot(barley.m, aes(x=barley.m$BarleyGCO, y=log10(barley.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Barley") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


#####Cassava 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
cassava<-na.omit(yield.mapping %>% dplyr::select(cassava_HgHa,CassavaGCO))

#Inside vs Outside GCO 
cassava.inside<-filter(cassava, CassavaGCO == "Inside")
cassava.outside<-filter(cassava, CassavaGCO == "Outside")
cassava.m<-melt(cassava)
summary(cassava.inside)
summary(cassava.outside)

#Plot
Cassava<-ggplot(cassava.m, aes(x=cassava.m$CassavaGCO, y=log10(cassava.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Cassava") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


#####Maize 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
Maize<-na.omit(yield.mapping %>% dplyr::select(maize_HgHa,MaizeGCO))

#Inside vs Outside GCO 
Maize.inside<-filter(Maize, MaizeGCO == "Inside")
Maize.outside<-filter(Maize, MaizeGCO == "Outside")
Maize.m<-melt(Maize)

#Plot
Maize<-ggplot(Maize.m, aes(x=Maize.m$MaizeGCO, y=log10(Maize.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Maize") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 

#####Potato 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
Potato<-na.omit(yield.mapping %>% dplyr::select(mean_potato,PotatoGCO))

#Inside vs Outside GCO 
Potato.inside<-filter(Potato, PotatoGCO == "Inside")
Potato.outside<-filter(Potato, PotatoGCO == "Outside")
Potato.m<-melt(Potato)

#Plot
Potato<-ggplot(Potato.m, aes(x=Potato.m$PotatoGCO, y=log10(Potato.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Potato") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


#####Rapeseed 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
Rapeseed<-na.omit(yield.mapping %>% dplyr::select(rapeseed_HgHa,RapeseedGCO))

#Inside vs Outside GCO 
Rapeseed.inside<-filter(Rapeseed, RapeseedGCO == "Inside")
Rapeseed.outside<-filter(Rapeseed, RapeseedGCO == "Outside")
Rapeseed.m<-melt(Rapeseed)

#Plot
Rapeseed<-ggplot(Rapeseed.m, aes(x=Rapeseed.m$RapeseedGCO, y=log10(Rapeseed.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Rapeseed") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


#####Rice 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
Rice<-na.omit(yield.mapping %>% dplyr::select(rice_HgHa,RiceGCO))

#Inside vs Outside GCO 
Rice.inside<-filter(Rice, RiceGCO == "Inside")
Rice.outside<-filter(Rice, RiceGCO == "Outside")
Rice.m<-melt(Rice)

#Plot
Rice<-ggplot(Rice.m, aes(x=Rice.m$RiceGCO, y=log10(Rice.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Rice") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


#####Rye 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
Rye<-na.omit(yield.mapping %>% dplyr::select(rye_HgHa,RyeGCO))

#Inside vs Outside GCO 
Rye.inside<-filter(Rye, RyeGCO == "Inside")
Rye.outside<-filter(Rye, RyeGCO == "Outside")
Rye.m<-melt(Rye)

#Plot
Rye<-ggplot(Rye.m, aes(x=Rye.m$RyeGCO, y=log10(Rye.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Rye") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


#####Sorghum 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
Sorghum<-na.omit(yield.mapping %>% dplyr::select(sorghum_HgHa,SorghumGCO))

#Inside vs Outside GCO 
Sorghum.inside<-filter(Sorghum, SorghumGCO == "Inside")
Sorghum.outside<-filter(Sorghum, SorghumGCO == "Outside")
Sorghum.m<-melt(Sorghum)

#Plot
Sorghum<-ggplot(Sorghum.m, aes(x=Sorghum.m$SorghumGCO, y=log10(Sorghum.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Sorghum") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 



#####Soybean 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
Soybean<-na.omit(yield.mapping %>% dplyr::select(soybean_HgHa,SoybeanGCO))

#Inside vs Outside GCO 
Soybean.inside<-filter(Soybean, SoybeanGCO == "Inside")
Soybean.outside<-filter(Soybean, SoybeanGCO == "Outside")
Soybean.m<-melt(Soybean)

#Plot
Soybean<-ggplot(Soybean.m, aes(x=Soybean.m$SoybeanGCO, y=log10(Soybean.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Soybean") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


#####Sunflower 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
Sunflower<-na.omit(yield.mapping %>% dplyr::select(sunflower_HgHa,SunflowerGCO))

#Inside vs Outside GCO 
Sunflower.inside<-filter(Sunflower, SunflowerGCO == "Inside")
Sunflower.outside<-filter(Sunflower, SunflowerGCO == "Outside")
Sunflower.m<-melt(Sunflower)

#Plot
Sunflower<-ggplot(Sunflower.m, aes(x=Sunflower.m$SunflowerGCO, y=log10(Sunflower.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Sunflower") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 



#####Wheat 
mapping<-read.csv("Fishnets/GCO_Mapping.csv", header=T, fileEncoding = "UTF-8-BOM")
yield<-read.csv("Fishnets/yield_fishnet.csv", header=T)
yield.mapping<-merge(mapping,yield)
Wheat<-na.omit(yield.mapping %>% dplyr::select(wheat_HgHa,WheatGCO))

#Inside vs Outside GCO 
Wheat.inside<-filter(Wheat, WheatGCO == "Inside")
Wheat.outside<-filter(Wheat, WheatGCO == "Outside")
Wheat.m<-melt(Wheat)

#Plot
Wheat<-ggplot(Wheat.m, aes(x=Wheat.m$WheatGCO, y=log10(Wheat.m$value))) + 
  geom_jitter(color="black", alpha=0.2)+
  ggtitle("Wheat") + 
  stat_summary(fun.data=data_summary, mult=1, geom="pointrange", color="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  ylab("Yield (Tons/Ha)") +xlab("") + labs(color='GCO') 


plots=list(Barley,Cassava,Groundnut,Maize,Potato,Rapeseed,Rice,Rye,Sorghum,Soybean,Sunflower,Wheat)
ggarrange(plotlist =plots,ncol = 4,nrow = 3)


