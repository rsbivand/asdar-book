###################################################
# vis_mod.R
# packages: sp, maptools, gstat, RColorBrewer, classInt
# datasets: 


###################################################
### chunk number 1: 
###################################################
rm(list=ls())


###################################################
### chunk number 8:  
###################################################
library(sp)
data(meuse)
coordinates(meuse) <- c("x", "y")
plot(meuse)
title("points")


###################################################
### chunk number 9:  
###################################################
cc <- coordinates(meuse)
m.sl <- SpatialLines(list(Lines(list(Line(cc)), "1")))
plot(m.sl)
title("lines")


###################################################
### chunk number 10:  
###################################################
data(meuse.riv)
meuse.lst <- list(Polygons(list(Polygon(meuse.riv)), "meuse.riv")) 
meuse.sr <- SpatialPolygons(meuse.lst)
plot(meuse.sr, col = "grey")
title("polygons")


###################################################
### chunk number 11:  
###################################################
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")
image(meuse.grid, col = "grey")
title("grid")


###################################################
### chunk number 13:  
###################################################
image(meuse.grid, col = "lightgrey")
plot(meuse.sr, col = "grey", add = TRUE)
plot(meuse, add = TRUE)


###################################################
### chunk number 15:  
###################################################
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2),1,2))
plot(meuse.sr, axes = TRUE)
plot(meuse.sr, axes = FALSE)
axis(1, at = c(178000 + 0:2 * 2000), cex.axis = .7)
axis(2, at = c(326000 + 0:3 * 4000), cex.axis = .7)
box()
par(def.par)


###################################################
### chunk number 17:  
###################################################
oldpar = par(no.readonly = TRUE)
layout(matrix(c(1,2),1,2))
plot(meuse, axes = TRUE, cex = 0.6)
plot(meuse.sr, add = TRUE)
title("Sample locations")
 
par(mar=c(0,0,0,0)+.1)
plot(meuse, axes = FALSE, cex = 0.6)
plot(meuse.sr, add = TRUE)
box()
par(oldpar)


###################################################
### chunk number 19:  
###################################################
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


###################################################
### chunk number 21: 
###################################################
library(maptools) 


###################################################
### chunk number 22: 
###################################################
def.par <- par(no.readonly = TRUE)
nc <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], proj4string=CRS("+proj=longlat +datum=NAD27"))
rrt <- nc$SID74/nc$BIR74
brks <- quantile(rrt, seq(0,1,1/5))
cols <- grey.colors(length(brks)-1, 0.95, 0.55, 2.2)
plot(nc, col=cols[findInterval(rrt, brks, all.inside=TRUE)], axes = TRUE)
par(def.par)


###################################################
### chunk number 27:  
###################################################
def.par <- par(no.readonly = TRUE)
pin <- par("pin")
dxy <- apply(bbox(meuse), 1, diff)
ratio <- dxy[1]/dxy[2]
par(pin=c(ratio * pin[2], pin[2]), xaxs="i", yaxs="i")
plot(meuse, pch = 1)
box()
par(def.par)


###################################################
### chunk number 31: 
###################################################
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) <- TRUE
library(gstat)
zn.idw <- krige(log(zinc) ~ 1, meuse, meuse.grid)


###################################################
### chunk number 33: 
###################################################
def.par <- par(no.readonly = TRUE)
par(mar=c(0.1,0.1,2,0.1))
cols <- grey.colors(4, 0.95, 0.55, 2.2)
image(zn.idw, col = cols, breaks=log(c(100,200,400,800,1800)))
plot(meuse.sr, add = TRUE)
plot(meuse, pch = 1, cex = sqrt(meuse$zinc)/20, add = TRUE)
legVals <- c(100, 200, 500, 1000, 2000)
legend("left", legend=legVals, pch = 1, pt.cex = sqrt(legVals)/20, bty = "n",
 title="measured, ppm", cex=0.8, y.inter=0.9)
legend("topleft", fill = cols, legend=c("100-200","200-400","400-800",
 "800-1800"), bty = "n", title = "interpolated, ppm", cex=0.8, y.inter=0.9)
title("measured and interpolated zinc")
par(def.par)


###################################################
### chunk number 34: 
###################################################
library(gstat)
library(sp)
data(meuse)
coordinates(meuse) <- ~x+y
data(meuse.grid) 
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- T
zn <- krige(zinc~1,meuse,meuse.grid)
zn$direct <- zn$var1.pred
zn$log <- exp(krige(log(zinc)~1,meuse,meuse.grid)$var1.pred)


###################################################
### chunk number 36: 
###################################################
library(lattice)
print(levelplot(z~x+y|name, spmap.to.lev(zn[c("direct", "log")]), asp = "iso",
        cuts=4, col.regions=grey.colors(5, 0.90, 0.50, 2.2)),
	split = c(1,1,1,2), more = TRUE)
print(spplot(zn[c("direct", "log")], cuts=4,
        col.regions=grey.colors(5, 0.90, 0.50, 2.2)), split = c(1,2,1,2))


###################################################
### chunk number 37: 
###################################################
cuts=c(0,20,50,200,500,2000)
grys <- grey.colors(7, 0.90, 0.50, 2.2)
print(spplot(meuse[1:4], main = "ppm", cuts=cuts, cex=.5, col.regions=grys),
      split=c(1,1,2,1),more=T)
meuse$lead.st = as.vector(scale(meuse$lead))
meuse$zinc.st = as.vector(scale(meuse$zinc))
meuse$copper.st = as.vector(scale(meuse$copper))
meuse$cadmium.st = as.vector(scale(meuse$cadmium))
cuts=c(-1.2,0,1,2,3,5)
print(spplot(meuse, c("lead.st", "zinc.st", "cadmium.st", "copper.st"),
	main = "standardised", cex = .5, cuts = cuts, col.regions=grys),
        split=c(2,1,2,1))


###################################################
### chunk number 39: 
###################################################
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")
cl = ContourLines2SLDF(contourLines(as.image.SpatialGridDataFrame(meuse.grid["dist"])))
## modified to match sp_0.9-29, RSB, 090105
grys <- grey.colors(nlevels(cl$level), 0.90, 0.50, 2.2)
print(spplot(cl, colorkey=list(height=0.8, width=0.6), col.regions=grys), split = c(1,1,2,2), more=TRUE)
cuts = (0:10)/10
grys <- grey.colors(length(cuts), 0.90, 0.50, 2.2)
print(spplot(meuse.grid, "dist"
, colorkey=list(labels=list(at=cuts), at=cuts), col.regions=grys
, pretty = TRUE
), split = c(2,1,2,2), more = TRUE)
meuse.grid$f = factor(meuse.grid$ffreq, labels = c("annual", "2-5 yrs", "> 5 yrs"))
print(spplot(meuse.grid, "f", colorkey=list(height=0.4, width=0.6), col.regions=grey.colors(3, 0.90, 0.50, 2.2)), split = c(1,2,2,2), more=FALSE)


###################################################
### chunk number 40: 
###################################################
river <- list("sp.polygons", meuse.sr)
north <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(178750,332500),
    scale = 400)
scale <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(180200,329800), scale = 1000, fill=c("transparent","black"))
txt1 <- list("sp.text", c(180200, 329950), "0")
txt2 <- list("sp.text", c(181200, 329950), "1 km")
pts <- list("sp.points", meuse, pch = 3, col = "black")
meuse.layout <- list(river, north, scale, txt1, txt2, pts)


###################################################
### chunk number 42: 
###################################################
grys <- grey.colors(7, 0.90, 0.50, 2.2)
print(spplot(zn["log"], sp.layout = meuse.layout, cuts=6, col.regions=grys))


###################################################
### chunk number 45: 
###################################################
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
if (packageVersion("sp") < "1.1.0") {
  plot(meuse[!is.na(overlay(meuse,sps)),],pch=16,cex=.5,add=TRUE)
} else {
  plot(meuse[sps,],pch=16,cex=.5,add=TRUE)
}
title("locator")
par(def.par)


###################################################
### chunk number 46: 
###################################################
library(maptools)
prj <- CRS("+proj=longlat +datum=NAD27")
nc_shp <- system.file("shapes/sids.shp", package="maptools")[1]
nc <- readShapePoly(nc_shp, proj4string=prj)
"pt" <- structure(list(x = -78.6948396159729, y = 35.8043970058349), 
	.Names = c("x", "y"))


###################################################
### chunk number 48: 
###################################################
print(pt)
if (packageVersion("sp") < "1.1.0") {
  print(overlay(nc, SpatialPoints(cbind(pt$x,pt$y),proj4string=prj)))
} else {
  print(slot(nc[SpatialPoints(cbind(pt$x,pt$y),proj4string=prj),], "data"))
}


###################################################
### chunk number 54: 
###################################################
library(RColorBrewer)
library(classInt)
# RSB quietening greys
pal <- grey.colors(4, 0.95, 0.55, 2.2)
#pal <- brewer.pal(4, "Greys")[2:4]
q5 <- classIntervals(meuse$zinc, n=5, style="quantile")
q5
diff(q5$brks)


###################################################
### chunk number 56: 
###################################################
fj5 <- classIntervals(meuse$zinc, n=5, style="fisher")
fj5
diff(fj5$brks)


###################################################
### chunk number 58: 
###################################################
def.par <- par(no.readonly = TRUE)
par(mar=c(3,3,3,1)+0.1, mfrow=c(1,2))
plot(q5, pal=pal, main="Quantile", xlab="", ylab="")
plot(fj5, pal=pal, main="Fisher-Jenks", xlab="", ylab="")
par(def.par)


###################################################
### chunk number 59: 
###################################################
def.par <- par(no.readonly = TRUE)
par(mar=c(1,1,3,1)+0.1, mfrow=c(1,2))
q5Colours <- findColours(q5, pal)
plot(meuse, col=q5Colours, pch=19)
points(meuse, pch=1)
box()
title(main="Quantile")
legend("topleft", fill=attr(q5Colours, "palette"), legend=names(attr(q5Colours, "table")), bty="n", cex=0.8, y.intersp=0.8)
fj5Colours <- findColours(fj5, pal)
plot(meuse, col=fj5Colours, pch=19)
points(meuse, pch=1)
box()
title(main="Fisher-Jenks")
legend("topleft", fill=attr(fj5Colours, "palette"), legend=names(attr(fj5Colours, "table")), bty="n", cex=0.8, y.intersp=0.8)
par(def.par)



