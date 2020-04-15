def.par <- par(no.readonly = TRUE)
opar = par(no.readonly = TRUE)
par(mar = rep(1,4))
library(sp)
data(volcano)
grys <- grey.colors(8, 0.55, 0.95, grey_gamma)
layout(matrix(c(1,2,1,3,1,4),3,2,byrow=TRUE), c(3,1))
image(volcano, axes=F, col=grys, asp=1, main="a")
contour(volcano, add=T)
box()
image(volcano, axes=F, col='white', asp=1, main="b")
library(maptools)
x2 = ContourLines2SLDF(contourLines(volcano))
plot(x2, add=T)
box()
image(volcano, axes=F, col='white', asp=1, main="c")
plot(SpatialPolygons(list(Polygons(list(Polygon(coordinates(x2[x2$level == 140,]))), ID="x"))),
	add = T)
box()
image(volcano, axes=F, col=grys, asp=1, main="d")
x3l1 = coordinates(x2[x2$level == 160,])[[1]][[1]]
x3l2 = coordinates(x2[x2$level == 160,])[[1]][[2]]
x3 = SpatialPolygons(list(Polygons(list(Polygon(x3l1,hole=F), Polygon(x3l2,hole=T)), ID=c("x"))))

SP2TRI = function(x, debug = TRUE){
    p = x@polygons[[1]] # object of class Polygons
    p1 = p@Polygons[[1]] # outer Polygon
    p2 = p@Polygons[[2]] # inner Polygon
    stopifnot(!p1@hole)
    stopifnot(p2@hole)
    # find nearest point
    allcoords = rbind(p1@coords, p2@coords)
    n1 = nrow(p1@coords)
    n2 = nrow(p2@coords)
    dists = as.matrix(dist(allcoords))[((n1+1):(n1+n2)),1:n1]
    wm = which.min(dists)[1]
    ind1 = (wm %/% n2) + 1
    ind2 = wm %% n2
    if (debug)
        print(c(ind1,ind2))
    #plot polygon points:
    p1c = p1@coords
    p2c = p2@coords
    #plot shortest distance:
    if (debug)
        lines(rbind(p1c[ind1,], p2c[ind2,]))
    if (debug)
        points(rbind(p1c, p2c))
    p = rbind(p1c[c(ind1:n1,1:ind1),], p2c[c(ind2:n2,1:ind2),], p1c[ind1,])
    #polygon(p, col = 'red', border = NULL)
    polygon(p, angle=45, border = NA, density = 12)
}
plot(x3, col = 'transparent', add = T)
SP2TRI(x3, F)

box()
par(opar)
par(def.par)
cat("\n")
layout(matrix(1))
grey_gamma <- 2.2


