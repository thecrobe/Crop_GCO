#Readins
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
barleymodel<-read.csv(file="Models/Barley_RF.csv", header=T)
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID
m <- sp::merge(fishnet, barleymodel, by='Fishnet_ID')
barley<-sp::merge(m, mapping, by='Fishnet_ID')


barley.num<-barley@data %>% select(mean_barle, AET_mean, Barley_Fertilizer, Pesticide,GDP_Mean)
barley.scaled<-data.frame(scale(barley.num,center = TRUE,scale = TRUE))
barley.cat<-barley@data %>% select(COUNTRY.x,BarleyGCO, FISHNET_ID)

dim(barley.scaled)
dim(barley.cat)
barley.final<-cbind(barley.scaled,barley.cat)

barley.adj <- poly2nb(fishnet)

W.barley <- nb2mat(barley.adj, style = "B",zero.policy = TRUE) 
W.barley.rs <- nb2mat(barley.adj, style = "W",zero.policy=TRUE) 

barley.form  <- mean_barle ~ AET_mean + GDP_Mean + GDP_Mean + Pesticide + Barley_Fertilizer + BarleyGCO 

barley.iid <- inla(update(barley.form, . ~. + f(FISHNET_ID, model = "iid")),
                   data = as.data.frame(barley.final),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.predictor = list(compute = TRUE))

summary(barley.iid)

#Besag's improper
barley.besag <- inla(update(barley.form, . ~. +
                              f(ID, model = "besag", graph = W.barley)), 
                     data = as.data.frame(barley.final),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.predictor = list(compute = TRUE)
)
boston.tr$BESAG <- tmarg(boston.besag$marginals.fitted.values)


                      