###################################################
### chunk number 1: 
###################################################
rm(list=ls())
library(digest)
op <- options()
options("width"=70, warn=1, str = strOptions(strict.width="wrap", vec.len=2))
.epsNo <- 0
library(digest)
if (!exists("online")) online <- TRUE
if (!exists("chkDigest")) chkDigest <- TRUE
intamap <- "http://intamap.geo.uu.nl/~roger/ASDAR/data"


###################################################
### chunk number 2: figreset eval=FALSE
###################################################
## .iwidth <- 5
## .iheight <- 6
## .ipointsize <- 12


###################################################
### chunk number 3: 
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 4: afig eval=FALSE
###################################################
## .epsNo <- .epsNo + 1; file <- paste("Fig-cm-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)


###################################################
### chunk number 5: zfig eval=FALSE
###################################################
## dev.null <- dev.off()
## system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 6: afig_png eval=FALSE
###################################################
## .epsNo <- .epsNo + 1; file <- paste("Fig-cm-", .epsNo, ".png", sep="")
## png(filename=file, width = .iwidth, height = .iheight, pointsize = .ipointsize)


###################################################
### chunk number 7: zfig_png eval=FALSE
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 8: 
###################################################
pi * 10^2


###################################################
### chunk number 9: 
###################################################
"*"(pi, "^"(10, 2))


###################################################
### chunk number 10: 
###################################################
pi * (1:10)^2


###################################################
### chunk number 11: 
###################################################
x <- pi * 10^2
x
print(x)
print(x, digits=12)


###################################################
### chunk number 12: 
###################################################
class(x)
typeof(x)


###################################################
### chunk number 13: 
###################################################
class(cars)
typeof(cars)
names(cars)
summary(cars)


###################################################
### chunk number 14: 
###################################################
str(cars)


###################################################
### chunk number 15: 
###################################################
class(dist ~ speed)


###################################################
### chunk number 16: 
###################################################
lm(dist ~ speed, data=cars)


###################################################
### chunk number 17: 
###################################################
cars$qspeed <- cut(cars$speed, breaks=quantile(cars$speed), include.lowest=TRUE)
is.factor(cars$qspeed)


###################################################
### chunk number 18:  eval=FALSE
###################################################
## plot(dist ~ speed, data=cars)
## plot(dist ~ qspeed, data=cars)


###################################################
### chunk number 19: 
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1,2))
plot(dist ~ speed, data=cars, main="numerical: scatterplot")
plot(dist ~ qspeed, data=cars, main="factor: boxplots")
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 20: 
###################################################
lm(dist ~ qspeed, data=cars)


###################################################
### chunk number 21: 
###################################################
library(sp)


###################################################
### chunk number 22:  eval=FALSE
###################################################
## getClass("Spatial")


###################################################
### chunk number 23: 
###################################################
x <- getClass("Spatial")
xx <- capture.output(methods:::print.classRepresentation(x))
cat(xx[1:5], sep="\n")
cat(strwrap(xx[6:length(xx)], width=getOption("width"), exdent=5), sep="\n")


###################################################
### chunk number 24: 
###################################################
getClass("CRS")


###################################################
### chunk number 25: 
###################################################
m <- matrix(c(0,0,1,1), ncol=2, dimnames=list(NULL, c("min", "max")))
crs <- CRS(projargs=as.character(NA))
crs
S <- Spatial(bbox=m, proj4string=crs)
S


###################################################
### chunk number 26:  eval=FALSE
###################################################
## Spatial(matrix(c(350, 85, 370, 95), ncol=2, dimnames=list(NULL, c("min", "max"))), proj4string=CRS("+longlat"))


###################################################
### chunk number 27: 
###################################################
res <- try(Spatial(matrix(c(350, 85, 370, 95), ncol=2, dimnames=list(NULL, c("min", "max"))), proj4string=CRS("+longlat")), silent=TRUE)
cat(strwrap(res), sep="\n")


###################################################
### chunk number 28: 
###################################################
setwd("hsd_data")


###################################################
### chunk number 29: 
###################################################
CRAN_df <- read.table("CRAN051001a.txt", header=TRUE)
CRAN_mat <- cbind(CRAN_df$long, CRAN_df$lat)
row.names(CRAN_mat) <- 1:nrow(CRAN_mat)
str(CRAN_mat)


###################################################
### chunk number 30: 
###################################################
setwd("..")


###################################################
### chunk number 31:  eval=FALSE
###################################################
## getClass("SpatialPoints")


###################################################
### chunk number 32: 
###################################################
x <- getClass("SpatialPoints")
xx <- capture.output(methods:::print.classRepresentation(x))
cat(xx[1:5], sep="\n")
cat(strwrap(xx[6:length(xx)], width=getOption("width"), exdent=5), sep="\n")


###################################################
### chunk number 33: 
###################################################
llCRS <- CRS("+proj=longlat +ellps=WGS84")
CRAN_sp <- SpatialPoints(CRAN_mat, proj4string=llCRS)
summary(CRAN_sp)


###################################################
### chunk number 34: 
###################################################
bbox(CRAN_sp)


###################################################
### chunk number 35: 
###################################################
proj4string(CRAN_sp)
proj4string(CRAN_sp) <- CRS(as.character(NA))
proj4string(CRAN_sp)
proj4string(CRAN_sp) <- llCRS


###################################################
### chunk number 36: 
###################################################
brazil <- which(CRAN_df$loc == "Brazil")
brazil
coordinates(CRAN_sp)[brazil,]


###################################################
### chunk number 37: 
###################################################
summary(CRAN_sp[brazil,])


###################################################
### chunk number 38: 
###################################################
south_of_equator <- which(coordinates(CRAN_sp)[,2] < 0)
summary(CRAN_sp[-south_of_equator,])


###################################################
### chunk number 39: 
###################################################
str(row.names(CRAN_df))


###################################################
### chunk number 40: 
###################################################
CRAN_spdf1 <- SpatialPointsDataFrame(CRAN_mat, CRAN_df, proj4string=llCRS, match.ID=TRUE)
CRAN_spdf1[4,]
str(CRAN_spdf1$loc)
str(CRAN_spdf1[["loc"]])


###################################################
### chunk number 41: 
###################################################
s <- sample(nrow(CRAN_df))
CRAN_spdf2 <- SpatialPointsDataFrame(CRAN_mat, CRAN_df[s,], proj4string=llCRS, match.ID=TRUE)
all.equal(CRAN_spdf2, CRAN_spdf1)
CRAN_spdf2[4,]


###################################################
### chunk number 42: 
###################################################
CRAN_df1 <- CRAN_df
row.names(CRAN_df1) <- sample(c(outer(letters, letters, paste, sep="")), nrow(CRAN_df1))


###################################################
### chunk number 43:  eval=FALSE
###################################################
## CRAN_spdf3 <- SpatialPointsDataFrame(CRAN_mat, CRAN_df1, proj4string=llCRS, match.ID=TRUE)


###################################################
### chunk number 44: 
###################################################
res <- try(CRAN_spdf3 <- SpatialPointsDataFrame(CRAN_mat, CRAN_df1, proj4string=llCRS, match.ID=TRUE), silent=TRUE)
cat(strwrap(res), sep="\n")


###################################################
### chunk number 45:  eval=FALSE
###################################################
## getClass("SpatialPointsDataFrame")


###################################################
### chunk number 46: 
###################################################
x <- getClass("SpatialPointsDataFrame")
xx <- capture.output(methods:::print.classRepresentation(x))
cat(xx[1:5], sep="\n")
cat(strwrap(xx[6:length(xx)], width=getOption("width"), exdent=5), sep="\n")


###################################################
### chunk number 47: 
###################################################
names(CRAN_spdf1)
str(model.frame(lat ~ long, data=CRAN_spdf1), give.attr=FALSE)


###################################################
### chunk number 48: 
###################################################
CRAN_spdf4 <- SpatialPointsDataFrame(CRAN_sp, CRAN_df)
all.equal(CRAN_spdf4, CRAN_spdf2)


###################################################
### chunk number 49: 
###################################################
CRAN_df0 <- CRAN_df
coordinates(CRAN_df0) <- CRAN_mat
proj4string(CRAN_df0) <- llCRS
all.equal(CRAN_df0, CRAN_spdf2)
str(CRAN_df0, max.level=2)


###################################################
### chunk number 50: 
###################################################
CRAN_df1 <- CRAN_df
names(CRAN_df1)
coordinates(CRAN_df1) <- c("long", "lat")
proj4string(CRAN_df1) <- llCRS
str(CRAN_df1, max.level=2)


###################################################
### chunk number 51: 
###################################################
setwd("hsd_data")


###################################################
### chunk number 52: 
###################################################
turtle_df <- read.csv("seamap105_mod.csv")
summary(turtle_df)


###################################################
### chunk number 53: 
###################################################
setwd("..")


###################################################
### chunk number 54: 
###################################################
timestamp <- as.POSIXlt(strptime(as.character(turtle_df$obs_date), "%m/%d/%Y %H:%M:%S"), "GMT")
turtle_df1 <- data.frame(turtle_df, timestamp=timestamp)
turtle_df1$lon <- ifelse(turtle_df1$lon < 0, turtle_df1$lon+360, turtle_df1$lon)
turtle_sp <- turtle_df1[order(turtle_df1$timestamp),]
coordinates(turtle_sp) <- c("lon", "lat")
proj4string(turtle_sp) <- CRS("+proj=longlat +ellps=WGS84")


###################################################
### chunk number 55: 
###################################################
library(maptools)
gshhs.c.b <- system.file("share/gshhs_c.b", package="maptools")
pac <- Rgshhs(gshhs.c.b, level=1, xlim=c(130,250), ylim=c(15,60), verbose=FALSE)


###################################################
### chunk number 56: 
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
plot(pac$SP, axes=TRUE, col="grey85", xaxs="i", yaxs="i")
plot(turtle_sp, add=TRUE)
m_rle <- rle(months(turtle_sp$timestamp))
clen <- cumsum(m_rle$lengths[-length(m_rle$lengths)])-1
crds <- coordinates(turtle_sp)
text(crds[clen,], labels=m_rle$values[-1], pos=3, offset=1.5, srt=45)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 57: 
###################################################
getClass("Line")
getClass("Lines")


###################################################
### chunk number 58:  eval=FALSE
###################################################
## getClass("SpatialLines")


###################################################
### chunk number 59: 
###################################################
x <- getClass("SpatialLines")
xx <- capture.output(methods:::print.classRepresentation(x))
cat(xx[1:5], sep="\n")
cat(strwrap(xx[6:length(xx)], width=getOption("width"), exdent=5), sep="\n")


###################################################
### chunk number 60: 
###################################################
library(maps)
japan <- map("world", "japan", plot=FALSE)
p4s <- CRS("+proj=longlat +ellps=WGS84")
SLjapan <- map2SpatialLines(japan, proj4string=p4s)
str(SLjapan, max.level=2)


###################################################
### chunk number 61: 
###################################################
Lines_len <- sapply(slot(SLjapan, "lines"), function(x) length(slot(x, "Lines")))
table(Lines_len)


###################################################
### chunk number 62: 
###################################################
volcano_sl <- ContourLines2SLDF(contourLines(volcano))
t(slot(volcano_sl, "data"))


###################################################
### chunk number 63: 
###################################################
setwd("hsd_data")


###################################################
### chunk number 64: 
###################################################
llCRS <- CRS("+proj=longlat +ellps=WGS84")
auck_shore <- MapGen2SL("auckland_mapgen.dat", llCRS)
summary(auck_shore)


###################################################
### chunk number 65: 
###################################################
setwd("..")


###################################################
### chunk number 66: 
###################################################
lns <- slot(auck_shore, "lines")
table(sapply(lns, function(x) length(slot(x, "Lines"))))
islands_auck <- sapply(lns, function(x) {
  crds <- slot(slot(x, "Lines")[[1]], "coords")
  identical(crds[1,], crds[nrow(crds),])
})
table(islands_auck)


###################################################
### chunk number 67: 
###################################################
getClass("Polygon")


###################################################
### chunk number 68: 
###################################################
getClass("Polygons")


###################################################
### chunk number 69:  eval=FALSE
###################################################
## getClass("SpatialPolygons")


###################################################
### chunk number 70: 
###################################################
x <- getClass("SpatialPolygons")
xx <- capture.output(methods:::print.classRepresentation(x))
cat(xx[1:5], sep="\n")
cat(strwrap(xx[6:length(xx)], width=getOption("width"), exdent=5), sep="\n")


###################################################
### chunk number 71: 
###################################################
islands_sl <- auck_shore[islands_auck]
list_of_Lines <- slot(islands_sl, "lines")
islands_sp <- SpatialPolygons(lapply(list_of_Lines, function(x) {
    Polygons(list(Polygon(slot(slot(x, "Lines")[[1]], "coords"))),
      ID=slot(x, "ID"))
  }),
  proj4string=CRS("+proj=longlat +ellps=WGS84"))
summary(islands_sp)
slot(islands_sp, "plotOrder")
order(sapply(slot(islands_sp, "polygons"),
  function(x) slot(x, "area")), decreasing=TRUE)


###################################################
### chunk number 72: 
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1,2), mar=c(3,3,1,1)+0.1)
plot(auck_shore)
legend("bottomleft", legend="a)", bty="n")
plot(auck_shore)
plot(islands_sp, add=TRUE, col="grey")
legend("bottomleft", legend="b)", bty="n")
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 73: 
###################################################
library(maps)
library(maptools)
state.map <- map("state", plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(state.map$names, ":"), function(x) x[1])
state.sp <- map2SpatialPolygons(state.map, IDs=IDs,
  proj4string=CRS("+proj=longlat +ellps=WGS84"))


###################################################
### chunk number 74: 
###################################################
setwd("hsd_data")


###################################################
### chunk number 75: 
###################################################
sat <- read.table("state.sat.data_mod.txt", row.names=5, header=TRUE)
str(sat)
id <- match(row.names(sat), sapply(slot(state.sp, "polygons"),
  function(x) slot(x, "ID")))
row.names(sat)[is.na(id)]
state.spdf <- SpatialPolygonsDataFrame(state.sp, sat)
str(slot(state.spdf, "data"))
str(state.spdf, max.level=2)


###################################################
### chunk number 76: 
###################################################
setwd("..")


###################################################
### chunk number 77:  eval=FALSE
###################################################
## rownames(sat)[3] <- "Arizona"
## SpatialPolygonsDataFrame(state.sp, sat)


###################################################
### chunk number 78: 
###################################################
rownames(sat)[3] <- "Arizona"
res <- try(SpatialPolygonsDataFrame(state.sp, sat), silent=TRUE)
cat(strwrap(res), sep="\n")


###################################################
### chunk number 79: 
###################################################
DC <- "district of columbia"
not_dc <- !(row.names(slot(state.spdf, "data")) == DC)
state.spdf1 <- state.spdf[not_dc,]
length(slot(state.spdf1, "polygons"))
summary(state.spdf1)


###################################################
### chunk number 80: 
###################################################
#high <- Rgshhs("/home/rsb/tmp/GSHHS/gshhs_h.b", xlim=c(277,278), ylim=c(45.7,46.2))
if (online) {
  con <- url(paste(intamap, "high.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("high.RData")
if (chkDigest && !identical(digest(high), "7f3ee8487d4e53a8a3ddd4c8e91358f1"))
  stop("high.RData digest error")
#load(url("http://intamap.geo.uu.nl/~roger/ASDAR/data/", open="rb"))
#if (!identical(digest(high), "7f3ee8487d4e53a8a3ddd4c8e91358f1"))
#  stop("high.RData digest error")
manitoulin_sp <- high$SP


###################################################
### chunk number 81: 
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(1,1,1,1)+0.1)
plot(manitoulin_sp, pbg="grey75", col="grey95")
text(t(sapply(slot(slot(manitoulin_sp, "polygons")[[1]], "Polygons"), function(x) slot(x, "labpt")))[-c(1,2),], label=high$polydata$level[-c(1,2)], col="black", font=2)
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 82: 
###################################################
length(slot(manitoulin_sp, "polygons"))
sapply(slot(slot(manitoulin_sp, "polygons")[[1]], "Polygons"), function(x) slot(x, "hole"))
sapply(slot(slot(manitoulin_sp, "polygons")[[1]], "Polygons"), function(x) slot(x, "ringDir"))


###################################################
### chunk number 83: 
###################################################
getClass("GridTopology")


###################################################
### chunk number 84: 
###################################################
bb <- bbox(manitoulin_sp)
bb
cs <- c(0.01, 0.01)
cc <- bb[,1]+(cs/2)
cd <- ceiling(diff(t(bb))/cs)
manitoulin_grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
manitoulin_grd


###################################################
### chunk number 85:  eval=FALSE
###################################################
## getClass("SpatialGrid")


###################################################
### chunk number 86: 
###################################################
x <- getClass("SpatialGrid")
xx <- capture.output(methods:::print.classRepresentation(x))
cat(xx[1:5], sep="\n")
cat(strwrap(xx[6:length(xx)], width=getOption("width"), exdent=5), sep="\n")


###################################################
### chunk number 87: 
###################################################
p4s <- CRS(proj4string(manitoulin_sp))
manitoulin_SG <- SpatialGrid(manitoulin_grd, proj4string=p4s)
summary(manitoulin_SG)


###################################################
### chunk number 88: 
###################################################
#library(rgdal)
#auck_el1 <- readGDAL("/home/rsb/tmp/70042108/70042108.tif")
if (online) {
  con <- url(paste(intamap, "auck_el1.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("high.RData")
if (chkDigest && !identical(digest(auck_el1), "d603173d5bf2d3e7bca177d0c673d569"))
  stop("auck_el1.RData digest error")
#load(url("http://intamap.geo.uu.nl/~roger/ASDAR/data/auck_el1.RData", open="rb"))
#if (!identical(digest(auck_el1), "d603173d5bf2d3e7bca177d0c673d569"))
#  stop("auck_el1.RData digest error")
#library(maptools)
#auck2 <- Rgshhs("/home/rsb/tmp/GSHHS/gshhs_f.b", xlim=c(174.2,175.3), ylim=c(-37.5,-36.5), level=2)
#auck_gshhs <- auck2$SP
load(url("http://intamap.geo.uu.nl/~roger/ASDAR/data/auck_gshhs.RData", open="rb"))
if (!identical(digest(auck_gshhs), "b5992b6958d92c5928127ec734a60bde"))
  stop("auck_gshhs.RData digest error")
source("legend_image.R")


###################################################
### chunk number 89: 
###################################################
class(auck_el1)
slot(auck_el1, "grid")
slot(auck_el1, "grid.index")
slot(auck_el1, "coords")
slot(auck_el1, "bbox")
object.size(auck_el1)
object.size(slot(auck_el1, "data"))


###################################################
### chunk number 90: 
###################################################
is.na(auck_el1$band1) <- auck_el1$band1 <= 0 | auck_el1$band1 > 1e+4
summary(auck_el1$band1)


###################################################
### chunk number 91: 
###################################################
.iwidth <- 5
.iheight <- 4.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(1,1,1,1)+0.1)
grys <- grey.colors(8, 0.55, 0.95, 2.2)
image(auck_el1, "band1", col=grys)
plot(auck_gshhs, add=TRUE, pbg="white")
transect_sp <- SpatialPoints(coords=cbind(seq(174.458,175.3,0.000833333),
  c(-37.03625)), proj4string=CRS("+proj=longlat +ellps=WGS84"))
plot(transect_sp, add=TRUE, pch="-", cex=2)
legend_image(c(174.2,174.25), c(-37.5,-37.2), auck_el1$band1, vertical=TRUE, offset.leg=0.8, col=grys)
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 92: 
###################################################
#load("/home/rsb/tmp/70042108/auck_el2.RData")
load(url("http://intamap.geo.uu.nl/~roger/ASDAR/data/auck_el2.RData", open="rb"))
if (!identical(digest(auck_el2), "f842372f0b2660c53f7a072015cfbce7"))
  stop("auck_el2.RData digest error")


###################################################
### chunk number 93:  eval=FALSE
###################################################
## auck_el2 <- as(auck_el1, "SpatialPixelsDataFrame")


###################################################
### chunk number 94: 
###################################################
object.size(auck_el2)
object.size(slot(auck_el2, "grid.index"))
object.size(slot(auck_el2, "coords"))
sum(is.na(auck_el1$band1)) + nrow(slot(auck_el2, "coords"))
prod(slot(slot(auck_el2, "grid"), "cells.dim"))


###################################################
### chunk number 95: 
###################################################
auck_el_500 <- auck_el2[auck_el2$band1 > 500,]


###################################################
### chunk number 96: 
###################################################
cat(paste("Warning messages:", "1: grid has empty column/rows in dimension 1 in:", "     points2grid(points, tolerance)", "2: grid has empty column/rows in dimension 2 in:", "     points2grid(points, tolerance)", sep="\n"), "\n")


###################################################
### chunk number 97: 
###################################################
summary(auck_el_500)
object.size(auck_el_500)


###################################################
### chunk number 98: 
###################################################
data(meuse.grid)
mg_SP <- SpatialPoints(cbind(meuse.grid$x, meuse.grid$y))
summary(mg_SP)
mg_SPix0 <- SpatialPixels(mg_SP)
summary(mg_SPix0)
prod(slot(slot(mg_SPix0, "grid"), "cells.dim"))


###################################################
### chunk number 99: 
###################################################
mg_SPix1 <- as(mg_SP, "SpatialPixels")
summary(mg_SPix1)


###################################################
### chunk number 100: 
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
sT <- capture.output(print(Sys.time()))
cat("\n")
cat(paste("%", sT, sep=" "), sep="\n")


###################################################
### chunk number 101: 
###################################################
options(op)


