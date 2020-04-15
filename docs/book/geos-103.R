meuse.grid$part <- meuse.grid$part.a + meuse.grid$part.b
layout(matrix(1:2, 1, 2))
omar <- par("mar")
par(mar = rep(0,4))
image(meuse.grid["part"], col = 'gray') #$
lines(area)
par(mar = c(2,2,0.5,0)+.1)
hist(res, main = NULL, xlab = NULL, ylab = NULL)
par(mar = omar)


