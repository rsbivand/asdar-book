oopar <- par(mar=c(1,1,1,1)+0.1)
grys <- grey.colors(8, 0.55, 0.95, 2.2)
image(auck_el1, "band1", col=grys)
plot(auck_gshhs, add=TRUE, pbg="white")
transect_sp <- SpatialPoints(coords=cbind(seq(174.458,175.3,0.000833333),
  c(-37.03625)), proj4string=CRS("+proj=longlat +ellps=WGS84"))
plot(transect_sp, add=TRUE, pch="-", cex=2)
legend_image(c(174.2,174.25), c(-37.5,-37.2), auck_el1$band1, vertical=TRUE, offset.leg=0.8, col=grys)
par(oopar)


