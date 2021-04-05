
library(nlme)

cassava.near<-read.csv(file = "Near_Coded/Cassava_Near.csv", header=T)
cassava.near$logYield<-log10(cassava.near$mean_cassa)

cassava.near$rescale_ND <- cassava.near$NEAR_DIST/(1000*1000) # rescaled the distance variables to 1000s of km (assuming it was originally in metres). Avoids tiny parameters downstream


# regular OLS no variance structure
cassava.ols <- gls(logYield ~ rescale_ND, data = cassava.near)

#When adding the variance structure you can specify the relevant covariate.
#I think if you leave it out it just assumes that you want to model the variance as a funciton of fitted(),
#which, in the case of a single predictor probably works out to be the same thing
#but better to be explicit I think


# varFixed (variance changes linearly with X)
cassava.fixed <- update(cassava.ols, .~., weights = varFixed(~rescale_ND))

# varPower (variance changes as a power function with X)
cassava.power <- update(cassava.ols, . ~ ., weights = varPower(form = ~rescale_ND))

# varExp (variance changes as an exponential function of x)
cassava.exp <- update(cassava.ols, . ~., weights = varExp(form = ~ rescale_ND)) # error

# varConstPower (constant plus a power function of X (useful if X includes 0))

cassava.ConstPower <- update(cassava.ols, . ~., weights = varConstPower(form = ~ rescale_ND))

# compare all models by AIC
AIC(cassava.ols, cassava.fixed, cassava.power, cassava.ConstPower, cassava.exp)

#varExp model wins in this case
summary(cassava.exp)

# some old school predict and plot to visualize the fits
newx <- seq(1,16, length.out = 15)
newy <- predict(cassava.exp, newdata = data.frame(rescale_ND =newx))


# very lame hack to get grey background  <embarassed emoji>
plot.new()
polygon(c(-min(cassava.near[,1])^2,-min(cassava.near[,1])^2,max(cassava.near[,1])^2,max(cassava.near[,1])^2),c(-min(cassava.near[,2])^2,max(cassava.near[,2])^2,max(cassava.near[,2])^2,-min(cassava.near[,2])^2), col="gray90")
par(new=T)
#####

plot(logYield ~ rescale_ND, cassava.near, ylim = c(-2,3), pch = 16, col = "#0000FF20",
     main = "Cassava GLS",
     xlab = "Distance from GCO (x 1000 km)")
lines(newy ~ newx, col = "red", lwd = 2, lty = 2)

# predicted sigmas
# the formula here will be different depending on which varStruct you use
# see Chapter 3 of Zuur Mixed Effects Models and Extensions for Ecology in R for explanation
sds <- cassava.exp$sigma * exp(2*newx*attr(cassava.exp$apVar, "Pars")[1])

# plotting 1.98*estimated sd of variation, to give a kind of 95% psuedo-prediction interval
lines(newy + 1.98*sds ~ newx, col = "red", lty =2)
lines(newy - 1.98*sds ~ newx, col = "red", lty =2)

# add the homogeneous variance model for comparison
newy.ols <- predict(cassava.ols, newdata = data.frame(rescale_ND =newx))
sds.ols <- cassava.ols$sigma

lines(newy.ols ~ newx, col = "green4", lwd = 2, lty = 2)
lines(newy.ols + 1.98*sds.ols ~ newx, col = "green4", lty =2)
lines(newy.ols - 1.98*sds.ols ~ newx, col = "green4", lty =2)

legend("topleft", legend = c("OLS", "varExp"), lty = 2, lwd = 2, col=c("red", "green4"), bty ="n")

## on the original scale
plot(mean_cassa ~ rescale_ND, data = cassava.near, ylab = "Yield Original Scale",xlab = "Distance from GCO (x 1000 km)",
     ylim = c(0,60),
     main = "Cassava GLS Original Scale",
     col = "#0000FF20")
lines(10^(newy) ~ newx, col = "red", lwd = 2, lty = 2)
lines(10^(newy +1.98*sds) ~ newx, col = "red", lty =2)
lines(10^(newy - 1.98*sds) ~ newx, col = "red", lty =2)

lines(10^(newy.ols) ~ newx, col = "green4", lwd = 2, lty = 2)
lines(10^(newy.ols + 1.98*sds.ols) ~ newx, col = "green4", lty =2)
lines(10^(newy.ols - 1.98*sds.ols) ~ newx, col = "green4", lty =2)
legend("topleft", legend = c("OLS", "varExp"), lty = 2, lwd = 2, col=c("red", "green4"), bty ="n")

########################################
########################################
# Some other ideas I had while looking at this data

## any pattern in residuals by country?


plot(resid(cassava.exp) ~ factor(cassava.near$COUNTRY), las = 2, cex.axis = 0.7)
abline(h= 0, col = "red", lwd =2)
abline(v = 1:100, col = "gray80", lty = 2)

## Can we identify some countries on the original plot?

# make a quick function to speed things up
country_hull <- function(countryname){
  sub <- cassava.near[cassava.near$COUNTRY == countryname,]
  pts <- with(sub, chull(rescale_ND, logYield))
  pts <- c(pts,pts[1])
  lines(sub[,c("rescale_ND","logYield")][pts,])
  text(x = mean(sub$rescale_ND), y = max(sub$logYield)+0.2,labels = countryname) 
}

#redraw the plot from above
plot(logYield ~ rescale_ND, cassava.near, ylim = c(-2,3), pch = 16, col = "#0000FF20",
     main = "Cassava GLS",
     xlab = "Distance from GCO (x 1000 km)")

# add some countries!
country_hull("Panama")
country_hull("China")
country_hull("Gabon")
country_hull("Eritrea")

# OK. this is maybe not so useful, but it helped me think about/investigate
# whether the unexplained variation is related  to some kind of country-level factors

## is there a "China effect"?

nochina <- cassava.near[cassava.near$COUNTRY!="China",]

nochina.ols <-  gls(logYield ~ rescale_ND, data = nochina)
summary(nochina.ols)
summary(cassava.ols)

nochina.exp <- update(nochina.ols, .~.,weights = varExp(form = ~rescale_ND))
