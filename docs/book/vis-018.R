def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2),1,2))
plot(meuse, axes=TRUE, cex = 0.6)
plot(meuse.sr, add=TRUE)
title("Sample locations")

par(mar=rep(0,4)+.1)
plot(meuse, axes=FALSE, cex = 0.6)
plot(meuse.sr, add=TRUE)
box()
par(def.par)
cat("\n")


