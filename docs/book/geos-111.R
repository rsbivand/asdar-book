source("m1m2.R")
layout(matrix(1:2, 1, 2))
omar <- par("mar")
par(mar = rep(0,4))
plot(meuse)
points(meuse[m1 < quantile(m1,.1),], pch=1)
points(meuse[m1 > quantile(m1,.9),], pch=16)
plot(meuse)
points(meuse[m2 < quantile(m2,.1),], pch=16)
points(meuse[m2 > quantile(m2,.9),], pch=1)
par(mar = omar)


