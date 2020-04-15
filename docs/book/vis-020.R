def.par <- par(no.readonly = TRUE)
par(mar = c(1,1,1,1))
plot(meuse, cex = 0.6)
plot(meuse.sr, col = "lightgrey", add = TRUE)
plot(meuse, cex = 0.6, add = TRUE)
SpatialPolygonsRescale(layout.scale.bar(), offset = c(180200,329600),
	scale = 1000, fill=c("transparent","black"), plot.grid = FALSE)
text(x = c(180200,181200), y = rep(329750, 2), c("0", "1 km"))
SpatialPolygonsRescale(layout.north.arrow(), offset = c(178750,332500), 
	scale = 400, plot.grid = FALSE)
box()
par(def.par)
cat("\n")


