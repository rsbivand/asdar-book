print(xyplot(log(zinc)~sqrt(dist), as.data.frame(meuse), asp = 1), split = c(1,1,2,1), more = TRUE)
zn.lm <- lm(log(zinc)~sqrt(dist), meuse)
meuse$fitted.s <- predict(zn.lm, meuse) - mean(predict(zn.lm, meuse))
meuse$residuals <- residuals(zn.lm)
print(spplot(meuse, c("fitted.s", "residuals"), cuts=4, col.regions=grey.colors(5, 0.95, 0.55, 2.2)), split = c(2,1,2,1))


