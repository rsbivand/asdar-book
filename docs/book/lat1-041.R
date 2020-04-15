library(RColorBrewer)
pal <- grey.colors(9, 1, 0.5, 2.2)
oopar <- par(mfrow=c(1,3), mar=c(1,1,3,1)+0.1)
z <- t(listw2mat(Sy0_lw_W))
brks <- c(0,0.1,0.143,0.167,0.2,0.5,1)
nbr3 <- length(brks)-3
image(1:63, 1:63, z[,ncol(z):1], breaks=brks, col=pal[c(1,(9-nbr3):9)],
 main="W style", axes=FALSE)
box()
z <- t(listw2mat(Sy0_lw_B))
image(1:63, 1:63, z[,ncol(z):1], col=pal[c(1,9)], main="B style", axes=FALSE)
box()
z <- t(listw2mat(Sy0_lw_idwB))
brks <- c(0,0.35,0.73,0.93,1.2,2.6)
nbr3 <- length(brks)-3
image(1:63, 1:63, z[,ncol(z):1], breaks=brks, col=pal[c(1,(9-nbr3):9)],
 main="IDW B style", axes=FALSE)
box()
par(oopar)


