grys <- grey.colors(11, 0.90, 0.50, 2.2)
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")
cl = ContourLines2SLDF(contourLines(as.image.SpatialGridDataFrame(meuse.grid["dist"])))
print(spplot(cl, colorkey=list(height=0.8, width=0.6), col.regions=grys), split = c(1,1,2,2), more=TRUE)
cuts = (0:10)/10
print(spplot(meuse.grid, "dist"
, colorkey=list(labels=list(at=cuts), at=cuts), col.regions=grys
, pretty = TRUE
), split = c(2,1,2,2), more = TRUE)
meuse.grid$f = factor(meuse.grid$ffreq, labels = c("annual", "2-5 yrs", "> 5 yrs"))
print(spplot(meuse.grid, "f", colorkey=list(height=0.4, width=0.6), col.regions=grey.colors(3, 0.90, 0.50, 2.2)), split = c(1,2,2,2), more=FALSE)
cat("\n")


