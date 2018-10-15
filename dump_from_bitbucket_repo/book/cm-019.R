oopar <- par(mfrow=c(1,2))
plot(dist ~ speed, data=cars, main="numerical: scatterplot")
plot(dist ~ qspeed, data=cars, main="factor: boxplots")
par(oopar)


