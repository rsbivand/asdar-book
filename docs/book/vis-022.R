def.par <- par(no.readonly = TRUE)
nc <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], proj4string=CRS("+proj=longlat +datum=NAD27"))
rrt <- nc$SID74/nc$BIR74
brks <- quantile(rrt, seq(0,1,1/5))
cols <- grey.colors(length(brks)-1, 0.95, 0.55, 2.2)
plot(nc, col=cols[findInterval(rrt, brks, all.inside=TRUE)], axes = TRUE)
par(def.par)
cat("\n")


