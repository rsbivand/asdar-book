def.par <- par(no.readonly = TRUE)
opar <- par(mar=c(0,0,0,0))
image(meuse.grid, col = "lightgrey")
plot(meuse.sr, col = "grey", add = TRUE)
plot(meuse, add = TRUE)
par(def.par)
cat("\n")


