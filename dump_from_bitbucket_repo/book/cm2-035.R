oopar <- par(mar=c(3,2,1,1)+0.1)
plot(nc_sc_va90, border="grey", axes=TRUE)
plot(ncscva_MA, lwd=2, add=TRUE)
text(coordinates(ncscva_MA), labels=sapply(slot(ncscva_MA, "polygons"), function(x) slot(x, "ID")), cex=0.6)
par(oopar)


