library(rgdal)
library(dplyr)
library(brms)
library(spdep)
library(spdplyr)

#Read In
fishnet<- readOGR(dsn= "GIS/", layer="Fishnet_yield_NoAntarctica")
proj4string(fishnet) <- CRS("+init=epsg:3786")
barleymodel<-read.csv(file="Models/Barley_RF.csv", header=T)
mapping<-read.csv(file="Fishnets/GCO_Mapping.csv", header=T)
mapping$Fishnet_ID<-mapping$FISHNET_ID

m <- sp::merge(fishnet, barleymodel, by='Fishnet_ID')
barley<-sp::merge(m, mapping, by='Fishnet_ID')
barley.p = subset(barley, mean_barle > 0)

barley.num<-barley.p@data %>% select(mean_barle, AET_mean, Barley_Fertilizer, Pesticide,GDP_Mean)
barley.scaled<-data.frame(scale(barley.num,center = TRUE,scale = TRUE))
barley.scaled$FISHNET_ID<-barley.p$FISHNET_ID
barley.cat<-barley.p %>% select(COUNTRY.x,BarleyGCO, FISHNET_ID,Latitude,Longitude)

county.pivot<-barley.cat@data %>% group_by(COUNTRY.x) %>% dplyr::summarize(CountyPixelCount = dplyr::n())

barley.catagory<-sp::merge(barley.cat,county.pivot, by="COUNTRY.x")
barley.comb<-sp::merge(barley.catagory,barley.scaled, by="FISHNET_ID")

barley.final<-na.omit(filter(barley.comb,CountyPixelCount > 500))
barley.final$BarleyGCO<-as.factor(barley.final$BarleyGCO)
barley.final$COUNTRY.x<-as.factor(barley.final$COUNTRY.x)
barley.coords<-barley.final %>% select(Latitude,Longitude)
dim(barley.final)

barley.adj <- nb2mat(poly2nb(barley.final), style = "B",zero.policy = TRUE)
#rownames(barley.adj)<-barley.final$FISHNET_ID
rownames(barley.final@data) = barley.final$FISHNET_ID

dist.mat <- as.matrix(dist(cbind(barley.final$Latitude, barley.final$Longitude))) 
W <- matrix(dist.mat, nrow = nrow(dist.mat), ncol = ncol(dist.mat))
W[W < 100] <- 1 
W[W > 100] <- 0 
rownames(W) <- barley.final$FISHNET_ID
rownames(barley.final@data)<-barley.final$FISHNET_ID


barley.icar<-brm(mean_barle ~ Barley_Fertilizer +car(type="icar", M=W,gr = COUNTRY.x ),data = barley.final@data, data2 =list(W=distance),family="gaussian", warmup = 1000, iter = 5000,thin = 5,cores = 4, chains = 4, seed = 10,control = list(max_treedepth = 12))

str(W)
str(barley.final)

saveRDS(barley.car, file = "BRMS/barley_car_model.rds")
summary(barley.car)
bayes_R2(barley.car)
library(brms)


barley.car<-readRDS("BRMS/barley_car_model.rds")
barley.car$fit
z<-summary(barley.car)




## Not run: 
# generate some spatial data
east <- north <- 1:10
Grid <- expand.grid(east, north)
K <- nrow(Grid)

# set up distance and neighbourhood matrices
distance <- as.matrix(dist(Grid))
W <- array(0, c(K, K))
W[distance == 1] <- 1 	

# generate the covariates and response data
x1 <- rnorm(K)
x2 <- rnorm(K)
theta <- rnorm(K, sd = 0.05)
phi <- rmulti_normal(
  1, mu = rep(0, K), Sigma = 0.4 * exp(-0.1 * distance)
)
eta <- x1 + x2 + phi
prob <- exp(eta) / (1 + exp(eta))
size <- rep(50, K)
y <- rbinom(n = K, size = size, prob = prob)
dat <- data.frame(y, size, x1, x2)

# fit a CAR model
fit <- brm(y | trials(size) ~ x1 + x2 + car(W), 
           data = dat, data2 = list(W = W),
           family = binomial()) 
summary(fit)

## End(Not run)