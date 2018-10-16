### R code from vignette source 'vis.Rnw'

###################################################
### code chunk number 1: vis.Rnw:6-16
###################################################
if (!exists("book_R_dont_trash")) rm(list=ls())
require(lattice)
op <- options()
options("width"=70, warn=1, str = strOptions(strict.width="wrap", vec.len=2), useFancyQuotes="TeX")
.epsNo <- 0
#library(digest)
#if (!exists("online")) online, 4, nchar(file)) <- TRUE

#if (!exists("chkDigest")) chkDigest <- TRUE
#intamap <- "http://intamap.geo.uu.nl/~roger/ASDAR/data"


###################################################
### code chunk number 2: figreset (eval = FALSE)
###################################################
## .iwidth <- 5
## .iheight <- 6
## .ipointsize <- 12


###################################################
### code chunk number 3: vis.Rnw:24-25
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 4: afig (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
## lattice.options(default.theme = col.whitebg())


###################################################
### code chunk number 5: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 6: afig_l (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
## .iheight, pointsize = .ipointsize, horizontal=FALSE)
## lattice.options(default.theme = col.whitebg())


###################################################
### code chunk number 7: zfig_l (eval = FALSE)
###################################################
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
## cat("\n")


###################################################
### code chunk number 8: vis.Rnw:138-171
###################################################
.iwidth <- 5
.iheight <- 2.5
.ipointsize <- 10
#.pwd <- 0.7
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
oldpar = par(mar = c(0,0,1,0))
library(sp)
data(meuse)
coordinates(meuse) = c("x", "y")
#layout(matrix(1:4, 2, 2, byrow = TRUE))
layout(matrix(1:4, 1, 4, byrow = TRUE))
par(mar = c(0,0,1,0))
plot(meuse, cex = 0.6)
title("points")

cc = coordinates(meuse)
m.sl = SpatialLines(list(Lines(list(Line(cc)), "mess")))
plot(m.sl)
title("lines")

data(meuse.riv)
meuse.lst = list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))
meuse.pol = SpatialPolygons(meuse.lst)
plot(meuse.pol, col = "grey")
title("polygons")

data(meuse.grid)
coordinates(meuse.grid) = c("x", "y")
meuse.grid <- as(meuse.grid, "SpatialPixels")
image(meuse.grid, col = "grey")
title("grid")
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 9: vis.Rnw:180-185 (eval = FALSE)
###################################################
## library(sp)
## data(meuse)
## coordinates(meuse) <- c("x", "y")
## plot(meuse)
## title("points")


###################################################
### code chunk number 10: vis.Rnw:197-201 (eval = FALSE)
###################################################
## cc <- coordinates(meuse)
## m.sl <- SpatialLines(list(Lines(list(Line(cc)), "line1")))
## plot(m.sl)
## title("lines")


###################################################
### code chunk number 11: vis.Rnw:213-218 (eval = FALSE)
###################################################
## data(meuse.riv)
## meuse.lst <- list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))
## meuse.pol <- SpatialPolygons(meuse.lst)
## plot(meuse.pol, col = "grey")
## title("polygons")


###################################################
### code chunk number 12: vis.Rnw:230-235 (eval = FALSE)
###################################################
## data(meuse.grid)
## coordinates(meuse.grid) <- c("x", "y")
## meuse.grid <- as(meuse.grid, "SpatialPixels")
## image(meuse.grid, col = "grey")
## title("grid")


###################################################
### code chunk number 13: vis.Rnw:266-278
###################################################
.iwidth <- 4
.iheight <- 5
.ipointsize <- 10
.pwd <- 0.5
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
opar <- par(mar=c(0,0,0,0))
# "lightsteelblue2", col="khaki2"
image(meuse.grid, col = "khaki2")
plot(meuse.pol, col = "lightsteelblue2", add = TRUE)
plot(meuse, add = TRUE, col = "brown", cex = .5)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 14: vis.Rnw:289-292 (eval = FALSE)
###################################################
## image(meuse.grid, col = "lightgrey")
## plot(meuse.pol, col = "grey", add = TRUE)
## plot(meuse, add = TRUE)


###################################################
### code chunk number 15: vis.Rnw:317-326
###################################################
.pwd <- 0.5
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
layout(matrix(c(1,2),1,2))
plot(meuse.pol, axes = TRUE); title("axes = TRUE")
plot(meuse.pol, axes = FALSE); title("axes added")
axis(1, at = c(178000 + 0:2 * 2000), cex.axis = .7)
axis(2, at = c(326000 + 0:3 * 4000), cex.axis = .7)
box()
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")


###################################################
### code chunk number 16: vis.Rnw:353-359 (eval = FALSE)
###################################################
## layout(matrix(c(1,2),1,2))
## plot(meuse.pol, axes = TRUE)
## plot(meuse.pol, axes = FALSE)
## axis(1, at = c(178000 + 0:2 * 2000), cex.axis = .7)
## axis(2, at = c(326000 + 0:3 * 4000), cex.axis = .7)
## box()


###################################################
### code chunk number 17: vis.Rnw:379-395
###################################################
.ipointsize <- 9
.iwidth <- 4
.iheight <- 3
.pwd <- .5
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
layout(matrix(c(1,2),1,2))
plot(meuse, axes=TRUE, cex = 0.6)
plot(meuse.pol, add=TRUE)
title("Sample locations")

par(mar=rep(0,4)+.1)
plot(meuse, axes=FALSE, cex = 0.6)
plot(meuse.pol, add=TRUE)
box()
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 18: vis.Rnw:426-437 (eval = FALSE)
###################################################
## oldpar = par(no.readonly = TRUE)
## layout(matrix(c(1,2),1,2))
## plot(meuse, axes = TRUE, cex = 0.6)
## plot(meuse.pol, add = TRUE)
## title("Sample locations")
## 
## par(mar=c(0,0,0,0)+.1)
## plot(meuse, axes = FALSE, cex = 0.6)
## plot(meuse.pol, add = TRUE)
## box()
## par(oldpar)


###################################################
### code chunk number 19: vis.Rnw:467-481
###################################################
.pwd <- 0.45
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
par(mar = c(1,1,1,1))
plot(meuse, cex = 0.7)
plot(meuse.pol, col = "lightgrey", add = TRUE)
plot(meuse, cex = 0.6, add = TRUE)
SpatialPolygonsRescale(layout.scale.bar(), offset = c(180200,329600),
    scale = 1000, fill=c("transparent","black"), plot.grid = FALSE)
text(x = c(180200,181200), y = rep(329750, 2), c("0", "1 km"))
SpatialPolygonsRescale(layout.north.arrow(), offset = c(178750,332500),
    scale = 400, plot.grid = FALSE)
box()
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")


###################################################
### code chunk number 20: vis.Rnw:498-507 (eval = FALSE)
###################################################
## plot(meuse)
## plot(meuse.pol, add=TRUE)
## plot(meuse)
## SpatialPolygonsRescale(layout.scale.bar(), offset = locator(1),
##     scale = 1000, fill=c("transparent","black"), plot.grid = FALSE)
## text(locator(1), "0")
## text(locator(1), "1 km")
## SpatialPolygonsRescale(layout.north.arrow(), offset = locator(1),
##     scale = 400, plot.grid = FALSE)


###################################################
### code chunk number 21: vis.Rnw:537-538
###################################################
library(maptools) # requires sp


###################################################
### code chunk number 22: vis.Rnw:543-567
###################################################
.pwd <- 0.7
.iwidth <- 5
.iheight <- 3.25
.ipointsize <- 8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
#par(pin = c(5,3))
nc <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1],
  proj4string=CRS("+proj=longlat +datum=NAD27"))
# is.projected(nc)
rrt <- nc$SID74/nc$BIR74
brks <- quantile(rrt, seq(0,1,1/5))
#cols <- grey((length(brks):2)/length(brks))
# RSB quieten greys
library(RColorBrewer)
cols <- brewer.pal(5, "Reds")
#cols <- grey.colors(length(brks)-1, 0.95, 0.55, 2.2)
#dens <- (2:length(brks))*3
#par(pin = c(4, 2.5), mar = c(2,2,1,1)+0.1)
plot(nc, col=cols[findInterval(rrt, brks, all.inside=TRUE)], axes = FALSE)
box()
degAxis(1)
degAxis(2, at=34:37)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 23: vis.Rnw:598-617 (eval = FALSE)
###################################################
## library(maptools)
## library(maps)
## wrld <- map("world", interior=FALSE, xlim=c(-179,179), 
##    ylim=c(-89,89), plot=FALSE)
## wrld_p <- pruneMap(wrld, xlim=c(-179,179))
## llCRS <- CRS("+proj=longlat +ellps=WGS84")
## wrld_sp <- map2SpatialLines(wrld_p, proj4string=llCRS)
## prj_new <- CRS("+proj=moll")
## library(rgdal)
## wrld_proj <- spTransform(wrld_sp, prj_new)
## wrld_grd <- gridlines(wrld_sp, easts=c(-179,seq(-150,150,50), 179.5),              
##   norths=seq(-75,75,15), ndiscr=100)
## wrld_grd_proj <- spTransform(wrld_grd, prj_new)
## at_sp <- gridat(wrld_sp, easts=0, norths=seq(-75,75,15), offset=0.3)
## at_proj <- spTransform(at_sp, prj_new)
## plot(wrld_proj, col="grey60")
## plot(wrld_grd_proj, add=TRUE, lty=3, col="grey70")
## text(coordinates(at_proj), pos=at_proj$pos, offset=at_proj$offset,                 
##   labels=parse(text=as.character(at_proj$labels)), cex=0.6)


###################################################
### code chunk number 24: vis.Rnw:636-638
###################################################
par("pin")
par(pin = c(4,4))


###################################################
### code chunk number 25: vis.Rnw:649-651 (eval = FALSE)
###################################################
## dev.off()
## X11(width = 10, height = 10)


###################################################
### code chunk number 26: vis.Rnw:659-660 (eval = FALSE)
###################################################
## pdf("file.pdf", width=5, height=7)


###################################################
### code chunk number 27: vis.Rnw:676-682 (eval = FALSE)
###################################################
## pin <- par("pin")
## dxy <- apply(bbox(meuse), 1, diff)
## ratio <- dxy[1]/dxy[2]
## par(pin=c(ratio * pin[2], pin[2]), xaxs="i", yaxs="i")
## plot(meuse, pch = 1)
## box()


###################################################
### code chunk number 28: vis.Rnw:688-700
###################################################
.pwd <- 0.45
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
pin <- par("pin")
dxy <- apply(bbox(meuse),1,diff)
ratio <- dxy[1]/dxy[2]
par(pin = c(ratio * pin[2], pin[2]))
par(xaxs="i")
par(yaxs="i")
plot(meuse, pch = 1)
box()
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 29: vis.Rnw:723-724 (eval = FALSE)
###################################################
## par(mfrow = c(2,3))


###################################################
### code chunk number 30: vis.Rnw:731-732 (eval = FALSE)
###################################################
## layout(matrix(1:6, 2, 3, byrow = TRUE))


###################################################
### code chunk number 31: vis.Rnw:780-785
###################################################
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) <- TRUE
library(gstat)
zn.idw <- krige(log(zinc) ~ 1, meuse, meuse.grid)


###################################################
### code chunk number 32: vis.Rnw:790-814
###################################################
.iwidth <- 3.5
.iheight <- 4
.ipointsize <- 9
.pwd <- 0.5
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
par(mar=c(0.1,0.1,2,0.1))
# RSB quietening greys
library(RColorBrewer)
cols <- brewer.pal(4, "Reds")
#cols <- grey.colors(4, 0.95, 0.55, 2.2)
#image(zn.idw, col = gray(5:8/10), breaks=log(c(100,200,400,800,1800)))
image(zn.idw, col = cols, breaks=log(c(100,200,400,800,1800)))
plot(meuse.pol, add = TRUE)
plot(meuse, pch = 1, cex = sqrt(meuse$zinc)/20, add = TRUE)
legVals <- c(100, 200, 500, 1000, 2000)
legend("left", legend=legVals, pch = 1, pt.cex = sqrt(legVals)/20, bty = "n",
  title="measured, ppm", cex=0.8, y.inter=0.9)
#legend("topleft", fill = gray(5:8/10), legend=c("100-200","200-400","400-800",
#"800-1800"), bty = "n", title = "interpolated, ppm", cex=0.8, y.inter=0.9)
legend("topleft", fill = cols, legend=c("100-200","200-400","400-800",
  "800-1800"), bty = "n", title = "interpolated, ppm", cex=0.8, y.inter=0.9)
title("measured and interpolated zinc")
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 33: vis.Rnw:822-832
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
### code chunk number 34: vis.Rnw:837-853
###################################################
.iwidth <- 8
.iheight <- 6
.pwd <- .8
#.ipointsize <- 9
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
# RSB quietening greys
print(levelplot(z~x+y|name, spmap.to.lev(zn[c("direct", "log")]), asp = "iso",
        cuts=4, col.regions= #grey.colors(5, 0.90, 0.50, 2.2)),
			brewer.pal(5, "Reds")),
    split = c(1,1,1,2), more = TRUE)
print(spplot(zn[c("direct", "log")], cuts=4,
        col.regions= #grey.colors(5, 0.90, 0.50, 2.2)), 
			brewer.pal(5, "Reds")),
		split = c(1,2,1,2))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")


###################################################
### code chunk number 35: vis.Rnw:896-905 (eval = FALSE)
###################################################
## grays = gray.colors(4, 0.55, 0.95)
## image(zn.idw, col = grays, breaks=log(c(100,200,400,800,1800)))
## plot(meuse.pol, add = TRUE)
## plot(meuse, pch = 1, cex = sqrt(meuse$zinc)/20, add = TRUE)
## legVals <- c(100, 200, 500, 1000, 2000)
## legend("left", legend=legVals, pch = 1, pt.cex = sqrt(legVals)/20, bty = "n",
##  title = "measured")
## legend("topleft", legend=c("100-200","200-400","400-800","800-1800"),
##  fill = grays, bty = "n", title = "interpolated")


###################################################
### code chunk number 36: vis.Rnw:946-948 (eval = FALSE)
###################################################
## library(lattice)
## levelplot(z~x+y|name, spmap.to.lev(zn[c("direct", "log")]), asp = "iso")


###################################################
### code chunk number 37: vis.Rnw:962-963 (eval = FALSE)
###################################################
## spplot(zn[c("direct", "log")])


###################################################
### code chunk number 38: vis.Rnw:997-1014
###################################################
.iwidth <- 6
.iheight <- 4.5
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
cuts=c(0,20,50,200,500,2000)
# RSB quietening greys
#grys <- grey.colors(7, 0.90, 0.50, 2.2)
grys <- brewer.pal(7, "Reds")
print(spplot(meuse[1:4], main = "ppm", cuts=cuts, cex=.5, col.regions=grys, key.space="right"),
      split=c(1,1,1,2),more=TRUE)
meuse$lead.st = as.vector(scale(meuse$lead))
meuse$zinc.st = as.vector(scale(meuse$zinc))
meuse$copper.st = as.vector(scale(meuse$copper))
meuse$cadmium.st = as.vector(scale(meuse$cadmium))
cuts=c(-1.2,0,1,2,3,5)
print(spplot(meuse, c("cadmium.st", "copper.st", "lead.st", "zinc.st"), key.space="right", main = "standardised", cex = .5, cuts = cuts, col.regions=grys),
        split=c(1,2,1,2))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 39: vis.Rnw:1025-1051
###################################################
.iwidth <- 6
.iheight <- 4
.pwd=.95
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
# RSB quietening greys
grys <- grey.colors(9, 0.90, 0.50, 2.2)
grys <- brewer.pal(9, "Reds")
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")
cl = ContourLines2SLDF(contourLines(as.image.SpatialGridDataFrame(
  meuse.grid["dist"])))
print(spplot(cl, colorkey=list(height=0.8, width=0.6), col.regions=grys), 
  split = c(1,1,3,1), more=TRUE)
grys <- grey.colors(11, 0.90, 0.50, 2.2)
grys <- brewer.pal(6, "Reds")
cuts = (0:5)/5
print(spplot(meuse.grid, "dist", at=cuts
, colorkey=list(labels=list(at=cuts), at=cuts), col.regions=grys)
, split = c(2,1,3,1), more = TRUE)
meuse.grid$f = factor(meuse.grid$ffreq, labels = c("annual", "2-5 yrs", 
  "> 5 yrs"))
print(spplot(meuse.grid, "f", colorkey=list(height=0.4, width=0.6), 
  col.regions= brewer.pal(3, "Reds") #grey.colors(3, 0.90, 0.50, 2.2)
  ), split = c(3,1,3,1), more=FALSE)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")


###################################################
### code chunk number 40: vis.Rnw:1098-1105 (eval = FALSE)
###################################################
## library(maptools)
## data(meuse.grid)
## coordinates(meuse.grid) <- c("x", "y")
## meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")
## im <- as.image.SpatialGridDataFrame(meuse.grid["dist"])
## cl <- ContourLines2SLDF(contourLines(im))
## spplot(cl)


###################################################
### code chunk number 41: vis.Rnw:1151-1161
###################################################
river <- list("sp.polygons", meuse.pol)
north <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = 
  c(178750,332500),
    scale = 400)
scale <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = 
  c(180200, 329800), scale = 1000, fill=c("transparent","black"))
txt1 <- list("sp.text", c(180200, 329950), "0")
txt2 <- list("sp.text", c(181200, 329950), "1 km")
pts <- list("sp.points", meuse, pch = 3, col = "black")
meuse.layout <- list(river, north, scale, txt1, txt2, pts)


###################################################
### code chunk number 42: vis.Rnw:1163-1164 (eval = FALSE)
###################################################
## spplot(zn["log"], sp.layout = meuse.layout)


###################################################
### code chunk number 43: vis.Rnw:1170-1179
###################################################
.iwidth <- 4
.iheight <- 5
.pwd=.45
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
# RSB quietening greys
#grys <- grey.colors(7, 0.90, 0.50, 2.2)
grys <- brewer.pal(7, "Reds")
print(spplot(zn["log"], sp.layout = meuse.layout, cuts=5, col.regions=grys))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")


###################################################
### code chunk number 44: vis.Rnw:1248-1256
###################################################
.iwidth <- 4
.iheight <- 5
.pwd=.45
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
m = as(meuse, "data.frame")
library(ggplot2)
ggplot(m, aes(x,y)) + geom_point() + coord_equal()
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")


###################################################
### code chunk number 45: vis.Rnw:1272-1274 (eval = FALSE)
###################################################
## library(ggplot2)
## methods(fortify)


###################################################
### code chunk number 46: vis.Rnw:1287-1289 (eval = FALSE)
###################################################
## m = as(meuse, "data.frame")
## ggplot(m, aes(x,y)) + geom_point() + coord_equal()


###################################################
### code chunk number 47: vis.Rnw:1324-1330 (eval = FALSE)
###################################################
## library(latticeExtra)
## p = spplot(meuse["zinc"])
## m = SpatialPolygonsDataFrame(meuse.pol, data.frame(col=1), match.ID = FALSE)
## l = spplot(m)
## l + p
## p + l


###################################################
### code chunk number 48: vis.Rnw:1344-1356
###################################################
.iwidth <- 6
.iheight <- 5
.pwd=.65
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(latticeExtra)
reds <- brewer.pal(7, "Reds")
p = spplot(meuse["zinc"], cuts = 5, col.regions = reds)
m = SpatialPolygonsDataFrame(meuse.pol, data.frame(col=1), match.ID = FALSE)
l = spplot(m, cuts = 5, col.regions = reds)
print(l + p, more = TRUE, split = c(1,1,2,1))
print(p + l, split = c(2,1,2,1))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")


###################################################
### code chunk number 49: vis.Rnw:1391-1393 (eval = FALSE)
###################################################
## plot(meuse)
## meuse.id <- identify(coordinates(meuse))


###################################################
### code chunk number 50: vis.Rnw:1403-1410 (eval = FALSE)
###################################################
## plot(meuse)
## region <- locator(type="o")
## n <- length(region$x)
## p <- Polygon(cbind(region$x,region$y)[c(1:n,1),], hole=FALSE)
## ps <- Polygons(list(p), ID = "region")
## sps <- SpatialPolygons(list(ps))
## plot(meuse[sps,], pch=16, cex=.5, add=TRUE)


###################################################
### code chunk number 51: vis.Rnw:1421-1454
###################################################
.iwidth <- 8
.iheight <- 5.1
.pwd <- 0.7
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
"region" <-
structure(list(x = c(180335.075864401, 180788.782724475, 181041.811550285,
181416.992223039, 181443.167618812, 181111.612605681, 180718.981669079,
180291.450204778, 180291.450204778), y = c(332617.465051335,
332124.858644764, 332081.647556468, 332677.960574949, 333412.549075975,
333723.668911704, 333377.980205339, 332833.520492813, 332677.960574949
)), .Names = c("x", "y"))
"meuse.id" <-
structure(list(ind = as.integer(c(82, 106, 109, 118, 155)), pos = 
  as.integer(c(1, 4, 4, 4, 4))), .Names = c("ind", "pos"))
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
par(mar = c(1,1,1,1))
plot(meuse)
box()
text(coordinates(meuse)[meuse.id$ind,], labels=meuse.id$ind, pos=meuse.id$pos, col = 'red')
title("identify")

plot(meuse)
box()
lines(region, type="o", col = 'red')
n = length(region$x)
p = Polygon(cbind(region$x,region$y)[c(1:n,1),], hole=FALSE)
ps = Polygons(list(p), ID = "region")
sps = SpatialPolygons(list(ps))
points(meuse[sps,], pch=16, cex=.5)
title("locator")
"pt" <- structure(list(x = -78.6948396159729, y = 35.8043970058349),
    .Names = c("x", "y"))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")


###################################################
### code chunk number 52: vis.Rnw:1471-1475
###################################################
library(maptools)
prj <- CRS("+proj=longlat +datum=NAD27")
nc_shp <- system.file("shapes/sids.shp", package="maptools")[1]
nc <- readShapePoly(nc_shp, proj4string=prj)


###################################################
### code chunk number 53: vis.Rnw:1477-1479 (eval = FALSE)
###################################################
## plot(nc)
## pt <- locator(type="p")


###################################################
### code chunk number 54: vis.Rnw:1481-1484
###################################################
print(pt)
pt.sp = SpatialPoints(cbind(pt$x,pt$y),proj4string=prj)
over(pt.sp, nc)


###################################################
### code chunk number 55: vis.Rnw:1500-1501 (eval = FALSE)
###################################################
## ids <- spplot(meuse, "zinc", identify=TRUE)


###################################################
### code chunk number 56: vis.Rnw:1515-1519 (eval = FALSE)
###################################################
## library(lattice)
## trellis.focus("panel", column = 1, row = 1)
## ids <- panel.identify()
## trellis.unfocus()


###################################################
### code chunk number 57: vis.Rnw:1529-1533 (eval = FALSE)
###################################################
## library(grid)
## trellis.focus("panel", column = 1, row = 1)
## as.numeric(grid.locator())
## trellis.unfocus()


###################################################
### code chunk number 58: vis.Rnw:1575-1577 (eval = FALSE)
###################################################
## rw.colors <- colorRampPalette(c("red", "white"))
## image(meuse.grid["dist"], col=rw.colors(10))


###################################################
### code chunk number 59: vis.Rnw:1594-1596 (eval = FALSE)
###################################################
## library(RColorBrewer)
## example(brewer.pal)


###################################################
### code chunk number 60: vis.Rnw:1648-1658
###################################################
library(RColorBrewer)
library(classInt)
# RSB quietening greys
#pal <- grey.colors(4, 0.95, 0.55, 2.2)
#pal <- brewer.pal(4, "Greys")[2:4]
# RSB colours!
pal <- brewer.pal(5, "Reds")
q5 <- classIntervals(meuse$zinc, n=5, style="quantile")
q5
diff(q5$brks)


###################################################
### code chunk number 61: vis.Rnw:1660-1661 (eval = FALSE)
###################################################
## plot(q5, pal=pal)


###################################################
### code chunk number 62: vis.Rnw:1688-1691
###################################################
fj5 <- classIntervals(meuse$zinc, n=5, style="fisher")
fj5
diff(fj5$brks)


###################################################
### code chunk number 63: vis.Rnw:1692-1693 (eval = FALSE)
###################################################
## plot(fj5, pal=pal)


###################################################
### code chunk number 64: vis.Rnw:1700-1711
###################################################
.iwidth <- 4
.iheight <- 2.5
.ipointsize <- 9
.pwd <- 0.85
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
oopar <- par(mar=c(3,3,3,1)+0.1, mfrow=c(1,2))
plot(q5, pal=pal, main="Quantile", xlab="", ylab="")
plot(fj5, pal=pal, main="Fisher-Jenks", xlab="", ylab="")
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 65: vis.Rnw:1722-1745
###################################################
.iwidth <- 7
.iheight <- 5
.ipointsize <- 12
.pwd <- 0.95
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-vis-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
oopar <- par(mar=c(1,1,3,1)+0.1, mfrow=c(1,2))
q5Colours <- findColours(q5, pal)
plot(meuse, col=q5Colours, pch=19)
points(meuse, pch=1)
box()
title(main="Quantile")
legend("topleft", fill=attr(q5Colours, "palette"), legend=names(attr(q5Colours,
  "table")), bty="n", cex=0.8, y.intersp=0.8)
fj5Colours <- findColours(fj5, pal)
plot(meuse, col=fj5Colours, pch=19)
points(meuse, pch=1)
box()
title(main="Fisher-Jenks")
legend("topleft", fill=attr(fj5Colours, "palette"),
 legend=names(attr(fj5Colours, "table")), bty="n", cex=0.8, y.intersp=0.8)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 66: vis.Rnw:1765-1769 (eval = FALSE)
###################################################
## q5Colours <- findColours(q5, pal)
## plot(meuse, col=q5Colours, pch=19)
## legend("topleft", fill=attr(q5Colours, "palette"),
##  legend=names(attr(q5Colours, "table")), bty="n")


###################################################
### code chunk number 67: vis.Rnw:1816-1818 (eval = FALSE)
###################################################
## cuts = (0:10)/10
## spplot(meuse.grid, "dist", colorkey=list(labels=list(at=cuts)), at=cuts)


###################################################
### code chunk number 68: vis.Rnw:1826-1834
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
cat("\n")
ver <- system("svnversion", intern=TRUE)
cat("%SVN version", ver, "\n")
sT <- capture.output(print(Sys.time()))
cat("\n")
cat(paste("%", sT, sep=" "), sep="\n")


