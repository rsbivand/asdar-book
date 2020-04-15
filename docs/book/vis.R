###################################################
### chunk number 1: 
###################################################
rm(list=ls())
.owidth <- getOption("width")
options("width"=70)
owarn <- options("warn")$warn
options(warn=1)
.epsNo <- 0


###################################################
### chunk number 2: figreset eval=FALSE
###################################################
## .iwidth <- 5.5
## .iheight <- 5.5
## .ipointsize <- 10
## .pwd <- 0.95


###################################################
### chunk number 3: 
###################################################
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95


###################################################
### chunk number 4: afig eval=FALSE
###################################################
## .epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
## def.par <- par(no.readonly = TRUE)


###################################################
### chunk number 5: zfig eval=FALSE
###################################################
## par(def.par)
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
## cat("\n")


###################################################
### chunk number 6: afig_l eval=FALSE
###################################################
## .epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)


###################################################
### chunk number 7: zfig_l eval=FALSE
###################################################
## dev.null <- dev.off()
## #system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
## cat("\n")


###################################################
### chunk number 8:  eval=FALSE
###################################################
## library(sp)
## data(meuse)
## coordinates(meuse) <- c("x", "y")
## plot(meuse)
## title("points")


###################################################
### chunk number 9:  eval=FALSE
###################################################
## cc <- coordinates(meuse)
## m.sl <- SpatialLines(list(Lines(list(Line(cc)))))
## plot(m.sl)
## title("lines")


###################################################
### chunk number 10:  eval=FALSE
###################################################
## data(meuse.riv)
## meuse.lst <- list(Polygons(list(Polygon(meuse.riv)), "meuse.riv")) 
## meuse.sr <- SpatialPolygons(meuse.lst)
## plot(meuse.sr, col = "grey")
## title("polygons")


###################################################
### chunk number 11:  eval=FALSE
###################################################
## data(meuse.grid)
## coordinates(meuse.grid) <- c("x", "y")
## meuse.grid <- as(meuse.grid, "SpatialPixels")
## image(meuse.grid, col = "grey")
## title("grid")


###################################################
### chunk number 12: 
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
oldpar = par(mar = c(0,0,1,0))
library(sp)
data(meuse)
coordinates(meuse) = c("x", "y")
layout(matrix(1:4, 2, 2, byrow = TRUE))
par(mar = c(0,0,1,0))
plot(meuse, cex = 0.6)
title("points")

cc = coordinates(meuse)
m.sl = SpatialLines(list(Lines(list(Line(cc)), "mess")))
plot(m.sl)
title("lines")

data(meuse.riv)
meuse.lst = list(Polygons(list(Polygon(meuse.riv)), "meuse.riv")) 
meuse.sr = SpatialPolygons(meuse.lst)
plot(meuse.sr, col = "grey")
title("polygons")

data(meuse.grid)
coordinates(meuse.grid) = c("x", "y")
meuse.grid <- as(meuse.grid, "SpatialPixels")
image(meuse.grid, col = "grey")
title("grid")
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95


###################################################
### chunk number 13:  eval=FALSE
###################################################
## image(meuse.grid, col = "lightgrey")
## plot(meuse.sr, col = "grey", add = TRUE)
## plot(meuse, add = TRUE)


###################################################
### chunk number 14: 
###################################################
.iwidth <- 3
.iheight <- 3
.ipointsize <- 10
.pwd <- 0.5
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
opar <- par(mar=c(0,0,0,0))
image(meuse.grid, col = "lightgrey")
plot(meuse.sr, col = "grey", add = TRUE)
plot(meuse, add = TRUE)
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95


###################################################
### chunk number 15:  eval=FALSE
###################################################
## layout(matrix(c(1,2),1,2))
## plot(meuse.sr, axes = TRUE)
## plot(meuse.sr, axes = FALSE)
## axis(1, at = c(178000 + 0:2 * 2000), cex.axis = .7)
## axis(2, at = c(326000 + 0:3 * 4000), cex.axis = .7)
## box()


###################################################
### chunk number 16: 
###################################################
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2),1,2))
plot(meuse.sr, axes = TRUE); title("axes = TRUE")
plot(meuse.sr, axes = FALSE); title("axes added")
axis(1, at = c(178000 + 0:2 * 2000), cex.axis = .7)
axis(2, at = c(326000 + 0:3 * 4000), cex.axis = .7)
box()
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")


###################################################
### chunk number 17:  eval=FALSE
###################################################
## oldpar = par(no.readonly = TRUE)
## layout(matrix(c(1,2),1,2))
## plot(meuse, axes = TRUE, cex = 0.6)
## plot(meuse.sr, add = TRUE)
## title("Sample locations")
## 
## par(mar=c(0,0,0,0)+.1)
## plot(meuse, axes = FALSE, cex = 0.6)
## plot(meuse.sr, add = TRUE)
## box()
## par(oldpar)


###################################################
### chunk number 18: 
###################################################
.ipointsize <- 9
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2),1,2))
plot(meuse, axes=TRUE, cex = 0.6)
plot(meuse.sr, add=TRUE)
title("Sample locations")

par(mar=rep(0,4)+.1)
plot(meuse, axes=FALSE, cex = 0.6)
plot(meuse.sr, add=TRUE)
box()
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95


###################################################
### chunk number 19:  eval=FALSE
###################################################
## plot(meuse)
## plot(meuse.sr, add=TRUE)
## plot(meuse)
## SpatialPolygonsRescale(layout.scale.bar(), offset = locator(1),
## 	scale = 1000, fill=c("transparent","black"), plot.grid = FALSE)
## text(locator(1), "0")
## text(locator(1), "1 km")
## SpatialPolygonsRescale(layout.north.arrow(), offset = locator(1), 
## 	scale = 400, plot.grid = FALSE)


###################################################
### chunk number 20: 
###################################################
.pwd <- 0.75
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
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
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")


###################################################
### chunk number 21: 
###################################################
library(maptools) # requires sp


###################################################
### chunk number 22: 
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 9
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
#par(pin = c(5,3))
nc <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], proj4string=CRS("+proj=longlat +datum=NAD27"))
# is.projected(nc)
rrt <- nc$SID74/nc$BIR74
brks <- quantile(rrt, seq(0,1,1/5))
#cols <- grey((length(brks):2)/length(brks))
# RSB quieten greys
cols <- grey.colors(length(brks)-1, 0.95, 0.55, 2.2)
#dens <- (2:length(brks))*3
#par(pin = c(4, 2.5), mar = c(2,2,1,1)+0.1)
plot(nc, col=cols[findInterval(rrt, brks, all.inside=TRUE)], axes = TRUE)
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95


###################################################
### chunk number 23:  eval=FALSE
###################################################
## library(maptools)
## library(maps)
## wrld <- map("world", interior=FALSE, xlim=c(-179,179), ylim=c(-89,89), plot=FALSE)
## wrld_p <- pruneMap(wrld, xlim=c(-179,179))
## llCRS <- CRS("+proj=longlat +ellps=WGS84")
## wrld_sp <- map2SpatialLines(wrld_p, proj4string=llCRS)
## prj_new <- CRS("+proj=moll")
## library(rgdal)
## wrld_proj <- spTransform(wrld_sp, prj_new)
## wrld_grd <- gridlines(wrld_sp, easts=c(-179,seq(-150,150,50),179.5), norths=seq(-75,75,15), ndiscr=100)
## wrld_grd_proj <- spTransform(wrld_grd, prj_new)
## at_sp <- gridat(wrld_sp, easts=0, norths=seq(-75,75,15), offset=0.3)
## at_proj <- spTransform(at_sp, prj_new)
## plot(wrld_proj, col="grey60")
## plot(wrld_grd_proj, add=TRUE, lty=3, col="grey70")
## text(coordinates(at_proj), pos=at_proj$pos, offset=at_proj$offset, labels=parse(text=as.character(at_proj$labels)), cex=0.6)


###################################################
### chunk number 24:  eval=FALSE
###################################################
## par("pin")
## par(pin = c(4,4))


###################################################
### chunk number 25:  eval=FALSE
###################################################
## dev.off()
## X11(width = 10, height = 10)


###################################################
### chunk number 26:  eval=FALSE
###################################################
## postscript("file.ps", width=10, height=10)


###################################################
### chunk number 27:  eval=FALSE
###################################################
## pin <- par("pin")
## dxy <- apply(bbox(meuse), 1, diff)
## ratio <- dxy[1]/dxy[2]
## par(pin=c(ratio * pin[2], pin[2]), xaxs="i", yaxs="i")
## plot(meuse, pch = 1)
## box()


###################################################
### chunk number 28: 
###################################################
.pwd <- 0.5
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
pin <- par("pin")
dxy <- apply(bbox(meuse),1,diff)
ratio <- dxy[1]/dxy[2]
par(pin = c(ratio * pin[2], pin[2]))
par(xaxs="i")
par(yaxs="i")
plot(meuse, pch = 1)
box()
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95


###################################################
### chunk number 29:  eval=FALSE
###################################################
## par(mfrow = c(2,3))


###################################################
### chunk number 30:  eval=FALSE
###################################################
## layout(matrix(1:6, 2, 3, byrow = TRUE))


###################################################
### chunk number 31: 
###################################################
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) <- TRUE
library(gstat)
zn.idw <- krige(log(zinc) ~ 1, meuse, meuse.grid)


###################################################
### chunk number 32:  eval=FALSE
###################################################
## grays = gray.colors(4, 0.55, 0.95)
## image(zn.idw, col = grays, breaks=log(c(100,200,400,800,1800)))
## plot(meuse.sr, add = TRUE)
## plot(meuse, pch = 1, cex = sqrt(meuse$zinc)/20, add = TRUE)
## legVals <- c(100, 200, 500, 1000, 2000)
## legend("left", legend=legVals, pch = 1, pt.cex = sqrt(legVals)/20, bty = "n",
##  title = "measured")
## legend("topleft", legend=c("100-200","200-400","400-800","800-1800"),
##  fill = grays, bty = "n", title = "interpolated")


###################################################
### chunk number 33: 
###################################################
.pwd <- 0.75
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
par(mar=c(0.1,0.1,2,0.1))
# RSB quietening greys
cols <- grey.colors(4, 0.95, 0.55, 2.2)
#image(zn.idw, col = gray(5:8/10), breaks=log(c(100,200,400,800,1800)))
image(zn.idw, col = cols, breaks=log(c(100,200,400,800,1800)))
plot(meuse.sr, add = TRUE)
plot(meuse, pch = 1, cex = sqrt(meuse$zinc)/20, add = TRUE)
legVals <- c(100, 200, 500, 1000, 2000)
legend("left", legend=legVals, pch = 1, pt.cex = sqrt(legVals)/20, bty = "n", title="measured, ppm", cex=0.8, y.inter=0.9)
#legend("topleft", fill = gray(5:8/10), legend=c("100-200","200-400","400-800","800-1800"), bty = "n", title = "interpolated, ppm", cex=0.8, y.inter=0.9)
legend("topleft", fill = cols, legend=c("100-200","200-400","400-800","800-1800"), bty = "n", title = "interpolated, ppm", cex=0.8, y.inter=0.9)
title("measured and interpolated zinc")
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95


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
### chunk number 35:  eval=FALSE
###################################################
## library(lattice)
## levelplot(z~x+y|name, spmap.to.lev(zn[c("direct", "log")]), asp = "iso")
## spplot(zn[c("direct", "log")])


###################################################
### chunk number 36: 
###################################################
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(lattice)
# RSB quietening greys
print(levelplot(z~x+y|name, spmap.to.lev(zn[c("direct", "log")]), asp = "iso",
        cuts=4, col.regions=grey.colors(5, 0.90, 0.50, 2.2)),
	split = c(1,1,1,2), more = TRUE)
print(spplot(zn[c("direct", "log")], cuts=4,
        col.regions=grey.colors(5, 0.90, 0.50, 2.2)), split = c(1,2,1,2))
dev.null <- dev.off()
#system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")


###################################################
### chunk number 37: 
###################################################
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
cuts=c(0,20,50,200,500,2000)
# RSB quietening greys
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
dev.null <- dev.off()
#system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")


###################################################
### chunk number 38:  eval=FALSE
###################################################
## library(maptools)
## data(meuse.grid)
## coordinates(meuse.grid) <- c("x", "y")
## meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")
## im <- as.image.SpatialGridDataFrame(meuse.grid["dist"])
## cl <- ContourLines2SLDF(contourLines(im))
## spplot(cl)


###################################################
### chunk number 39: 
###################################################
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
# RSB quietening greys
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
dev.null <- dev.off()
#system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")


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
### chunk number 41:  eval=FALSE
###################################################
## spplot(zn["log"], sp.layout = meuse.layout)


###################################################
### chunk number 42: 
###################################################
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
# RSB quietening greys
grys <- grey.colors(7, 0.90, 0.50, 2.2)
print(spplot(zn["log"], sp.layout = meuse.layout, cuts=6, col.regions=grys))
dev.null <- dev.off()
#system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")


###################################################
### chunk number 43:  eval=FALSE
###################################################
## plot(meuse)
## meuse.id <- identify(coordinates(meuse))


###################################################
### chunk number 44:  eval=FALSE
###################################################
## plot(meuse)
## region <- locator(type="o")
## n <- length(region$x)
## p <- Polygon(cbind(region$x,region$y)[c(1:n,1),], hole=FALSE)
## ps <- Polygons(list(p), ID = "region")
## sps <- SpatialPolygons(list(ps))
## plot(meuse[!is.na(overlay(meuse,sps)),],pch=16,cex=.5,add=TRUE)


###################################################
### chunk number 45: 
###################################################
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
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
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")


###################################################
### chunk number 46: 
###################################################
library(maptools)
prj <- CRS("+proj=longlat +datum=NAD27")
nc_shp <- system.file("shapes/sids.shp", package="maptools")[1]
nc <- readShapePoly(nc_shp, proj4string=prj)


###################################################
### chunk number 47:  eval=FALSE
###################################################
## plot(nc)
## pt <- locator(type="p")


###################################################
### chunk number 48: 
###################################################
print(pt)
overlay(nc, SpatialPoints(cbind(pt$x,pt$y),proj4string=prj))


###################################################
### chunk number 49:  eval=FALSE
###################################################
## ids <- spplot(meuse, "zinc", identify=TRUE)


###################################################
### chunk number 50:  eval=FALSE
###################################################
## library(lattice)
## trellis.focus("panel", column = 1, row = 1)
## ids <- panel.identify()
## trellis.unfocus()


###################################################
### chunk number 51:  eval=FALSE
###################################################
## library(grid)
## trellis.focus("panel", column = 1, row = 1)
## as.numeric(grid.locator())
## trellis.unfocus()


###################################################
### chunk number 52:  eval=FALSE
###################################################
## rw.colors <- colorRampPalette(c("red", "white"))
## image(meuse.grid["dist"], col=rw.colors(10))


###################################################
### chunk number 53:  eval=FALSE
###################################################
## library(RColorBrewer)
## example(brewer.pal)


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
### chunk number 55:  eval=FALSE
###################################################
## plot(q5, pal=pal)


###################################################
### chunk number 56: 
###################################################
fj5 <- classIntervals(meuse$zinc, n=5, style="fisher")
fj5
diff(fj5$brks)


###################################################
### chunk number 57:  eval=FALSE
###################################################
## plot(fj5, pal=pal)


###################################################
### chunk number 58: 
###################################################
.iwidth <- 4
.iheight <- 2.5
.ipointsize <- 9
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
oopar <- par(mar=c(3,3,3,1)+0.1, mfrow=c(1,2))
plot(q5, pal=pal, main="Quantile", xlab="", ylab="")
plot(fj5, pal=pal, main="Fisher-Jenks", xlab="", ylab="")
par(oopar)
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95


###################################################
### chunk number 59: 
###################################################
.iwidth <- 4
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-VIS-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
oopar <- par(mar=c(1,1,3,1)+0.1, mfrow=c(1,2))
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
par(oopar)
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
.iwidth <- 5.5
.iheight <- 5.5
.ipointsize <- 10
.pwd <- 0.95


###################################################
### chunk number 60:  eval=FALSE
###################################################
## q5Colours <- findColours(q5, pal)
## plot(meuse, col=q5Colours, pch=19)
## legend("topleft", fill=attr(q5Colours, "palette"),
##  legend=names(attr(q5Colours, "table")), bty="n")


###################################################
### chunk number 61:  eval=FALSE
###################################################
## cuts = (0:10)/10
## spplot(meuse.grid, "dist", colorkey=list(labels=list(at=cuts)), at=cuts)


###################################################
### chunk number 62: 
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
sT <- capture.output(print(Sys.time()))
cat("\n")
cat(paste("%", sT, sep=" "), sep="\n")


###################################################
### chunk number 63: 
###################################################
options(warn=owarn)


