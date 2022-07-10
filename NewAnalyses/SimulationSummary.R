library(ggplot2)
library(dplyr)
library(brms)
library(bayesplot)
library(tidybayes)

#Graphic setting
theme_justin<-theme_bw() +theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) 

#read in data
simulation<-read.csv(file="NewAnalyses/Simulations/YieldSimulation/ArtificialSummary.csv", header=T)
str(simulation)

#subset GCO and plot
GCO<-simulation %>% filter (Variable == "GCO" )

gcoplot<-ggplot(GCO, aes(x=log10(YieldFactor),y=log10(FertilizerStandard))) +
  geom_point(size=3) + 
  geom_smooth(method="lm") +
  theme_justin + 
  xlab("Log10 Artifical Yield Factor") + 
  ylab("Log10 Variable Importance")
print(gcoplot)

#Model relationship
bmod<-brm(formula = log10(YieldFactor) ~ log10(FertilizerStandard) ,family = "gaussian",chains = 4,iter = 10000,thin = 4,cores=4,data=simulation)

plot(bmod) #check fit
summary(bmod) #check fit
pp_check(bmod,ndraws = 1000) # posterior predictive check
bayes_R2(bmod) #bayesian R2

modplot<-  bmod %>%
  spread_draws(b_log10FertilizerStandard,sigma) %>%
  mutate(GCO = b_log10FertilizerStandard) %>%
  ggplot(aes(y = 1, x = GCO, ,fill = stat(x <0))) +
  stat_halfeye() + 
  theme_justin +
  xlab("GCO Effect") + 
  ylab("Density") + 
  scale_fill_manual(values = c("#C1C1C1","#353535"),name = "GCO Effect", labels = c("Yes", "No")) 

sigmaplot<-  bmod %>%
  spread_draws(b_log10FertilizerStandard,sigma) %>%
  mutate(GCO = b_log10FertilizerStandard) %>%
  ggplot(aes(y = 1, x = sigma, ,fill = stat(x <0))) +
  stat_halfeye() + 
  theme_justin +
  xlab("Residual SD") + 
  ylab("Density") + 
  scale_fill_manual(values = c("#C1C1C1","#353535"),name = "GCO Effect", labels = c("Yes", "No")) 

multipanel<-ggarrange(gcoplot,modplot,labels = "auto")
print(multipanel)

