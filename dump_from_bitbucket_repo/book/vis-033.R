def.par <- par(no.readonly = TRUE)
par(mar=c(0.1,0.1,2,0.1))
cols <- grey.colors(4, 0.95, 0.55, 2.2)
image(zn.idw, col = cols, breaks=log(c(100,200,400,800,1800)))
plot(meuse.sr, add = TRUE)
plot(meuse, pch = 1, cex = sqrt(meuse$zinc)/20, add = TRUE)
legVals <- c(100, 200, 500, 1000, 2000)
legend("left", legend=legVals, pch = 1, pt.cex = sqrt(legVals)/20, bty = "n", title="measured, ppm", cex=0.8, y.inter=0.9)
legend("topleft", fill = cols, legend=c("100-200","200-400","400-800","800-1800"), bty = "n", title = "interpolated, ppm", cex=0.8, y.inter=0.9)
title("measured and interpolated zinc")
par(def.par)
cat("\n")


