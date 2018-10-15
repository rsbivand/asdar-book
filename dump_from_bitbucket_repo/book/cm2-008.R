oopar <- par(mar=c(1,1,1,1)+0.1, mfrow=c(1,2))
x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
grys <- grey.colors(11, 0.9, 0.45, 2.2)
image(y, x, t(volcano)[ncol(volcano):1,], breaks=seq(90,200,10), col=grys, asp=1, axes=FALSE)
contour(y, x, t(volcano)[ncol(volcano):1,], levels=seq(90,200,10), asp=1, axes=FALSE)
par(oopar)


