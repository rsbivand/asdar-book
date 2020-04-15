def.par <- par(no.readonly = TRUE)
"region" <-
structure(list(x = c(180335.075864401, 180788.782724475, 181041.811550285, 
181416.992223039, 181443.167618812, 181111.612605681, 180718.981669079, 
180291.450204778, 180291.450204778), y = c(332617.465051335, 
332124.858644764, 332081.647556468, 332677.960574949, 333412.549075975, 
333723.668911704, 333377.980205339, 332833.520492813, 332677.960574949
)), .Names = c("x", "y"))
"meuse.id" <-
structure(list(ind = as.integer(c(82, 106, 109, 118, 155)), pos = as.integer(c(1, 
4, 4, 4, 4))), .Names = c("ind", "pos"))
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
par(mar = c(0,0,1,0))
plot(meuse)
box()
text(coordinates(meuse)[meuse.id$ind,], labels=meuse.id$ind, pos=meuse.id$pos)
title("identify")

plot(meuse)
box()
lines(region, type="o")
n = length(region$x)
p = Polygon(cbind(region$x,region$y)[c(1:n,1),], hole=FALSE)
ps = Polygons(list(p), ID = "region")
sps = SpatialPolygons(list(ps))
plot(meuse[!is.na(overlay(meuse,sps)),],pch=16,cex=.5,add=TRUE)
title("locator")
"pt" <- structure(list(x = -78.6948396159729, y = 35.8043970058349), 
	.Names = c("x", "y"))
par(def.par)
cat("\n")


