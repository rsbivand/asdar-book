###################################################
# csdacm_mod.R
# packages: sp, gstat, rgdal
# datasets: 70042108.zip, seamap105_mod.csv


###################################################
### chunk number 1: 
###################################################
rm(list=ls())
if ((site <- Sys.getenv("ASDAR_DOWNLOAD")) != "") {
  download.file(paste(site, "seamap105_mod.csv", sep="/"), "seamap105_mod.csv")
  download.file(paste(site, "70042108.zip", sep="/"), "70042108.zip")
}


###################################################
### chunk number 6: 
###################################################
myfun <- function(x) {
	x + 2
}


###################################################
### chunk number 7: 
###################################################
myfun(1:3)


###################################################
### chunk number 8: 
###################################################
myfun(x=1:3)


###################################################
### chunk number 9: 
###################################################
plotXplus2Yminus3 <- function(x, y, ...) {
	plot(x = x + 2, y = y - 3, ...)
}


###################################################
### chunk number 10: 
###################################################
methods("plot")


###################################################
### chunk number 11: 
###################################################
library(sp)
showMethods("plot")


###################################################
### chunk number 12: 
###################################################
x <- rnorm(10)
class(x) <- "foo"
x


###################################################
### chunk number 13: 
###################################################
plot.foo <- function(x, y, ...) {
	plot.default(x, type = 'l', ...)
}


###################################################
### chunk number 14: 
###################################################
class(x) <- c("foo", "bar")


###################################################
### chunk number 16: 
###################################################
data(meuse)
class(meuse)
class(lm(log(zinc)~sqrt(dist), meuse))


###################################################
### chunk number 20: 
###################################################
isGeneric("show")


###################################################
### chunk number 22: 
###################################################
library(sp)
setClass("trip", representation("SpatialPointsDataFrame", TOR.columns = "character"),
    validity <- function(object) {
        if (length(object@TOR.columns) != 2)
            stop("Time/id column names must be of length 2")
		if (!all(object@TOR.columns %in% names(object@data)))
			stop("Time/id columns must be present in attribute table")
        TRUE
    }
)
showClass("trip")


###################################################
### chunk number 23: 
###################################################
trip.default <- function(obj, TORnames) {
    if (!is(obj, "SpatialPointsDataFrame"))
        stop("trip only supports SpatialPointsDataFrame") 
	if (is.numeric(TORnames))
		TORnames <- names(obj)[TORnames]
    new("trip", obj, TOR.columns = TORnames)
}

if (!isGeneric("trip"))
    setGeneric("trip", function(obj, TORnames)
        standardGeneric("trip"))

setMethod("trip", signature(obj = "SpatialPointsDataFrame", TORnames = "ANY"), trip.default)


###################################################
### chunk number 25: 
###################################################
turtle <- read.csv("seamap105_mod.csv")


###################################################
### chunk number 27: 
###################################################
timestamp <- as.POSIXlt(strptime(as.character(turtle$obs_date), "%m/%d/%Y %H:%M:%S"), "GMT")
turtle <- data.frame(turtle, timestamp = timestamp)
turtle$lon <- ifelse(turtle$lon < 0, turtle$lon+360, turtle$lon)
turtle <- turtle[order(turtle$timestamp),]
coordinates(turtle) <- c("lon", "lat")
proj4string(turtle) <- CRS("+proj=longlat +ellps=WGS84")
turtle$id <- c(rep(1, 200), rep(2, nrow(coordinates(turtle)) - 200))
turtle_trip <- trip(turtle, c("timestamp", "id"))
summary(turtle_trip)


###################################################
### chunk number 28: 
###################################################
summary.trip <- function(object, ...) {
	cat("Object of class \"trip\"\nTime column: ")
	print(object@TOR.columns[1])
	cat("Identifier column: ")
	print(object@TOR.columns[2])
	print(summary(as(object, "Spatial")))
	print(summary(object@data))
}
setMethod("summary", "trip", summary.trip)
summary(turtle_trip)


###################################################
### chunk number 29: 
###################################################
setGeneric("lines", function(x, ...) standardGeneric("lines"))
setMethod("lines", signature(x = "trip"),
    function(x, ..., col = NULL) {
	tor <- x@TOR.columns
	if (is.null(col)) {
	  l <- length(unique(x[[tor[2]]]))
          col <- hsv(seq(0, 0.5, length = l))
	}
        coords <- coordinates(x)
        lx <- split(1:nrow(coords), x[[tor[2]]])
        for (i in 1:length(lx))
        	lines(coords[lx[[i]], ], col = col[i], ...)
    }
)


###################################################
### chunk number 31: 
###################################################
setClass("SpatialMultiPoints", representation("SpatialLines"), 
	validity <- function(object) {
		if (any(unlist(lapply(object@lines, function(x) length(x@Lines))) != 1))
# NOT TOO WIDE
			stop("Only Lines objects with one Line element")
		TRUE
	}
)
SpatialMultiPoints <- function(object) new("SpatialMultiPoints", object)


###################################################
### chunk number 33: 
###################################################
n <- 5
set.seed(1)
x1 <- cbind(rnorm(n),rnorm(n, 0, 0.25))
x2 <- cbind(rnorm(n),rnorm(n, 0, 0.25))
x3 <- cbind(rnorm(n),rnorm(n, 0, 0.25))
L1 <- Lines(list(Line(x1)), ID="mp1")
L2 <- Lines(list(Line(x2)), ID="mp2")
L3 <- Lines(list(Line(x3)), ID="mp3")
s <- SpatialLines(list(L1,L2,L3))
smp <- SpatialMultiPoints(s)


###################################################
### chunk number 34: 
###################################################
plot.SpatialMultiPoints <- function(x, ..., pch = 1:length(x@lines), col = 1, cex = 1) {
	n <- length(x@lines)
	if (length(pch) < n)
		pch <- rep(pch, length.out = n)
	if (length(col) < n)
		col <- rep(col, length.out = n)
	if (length(cex) < n)
		cex <- rep(cex, length.out = n)
	plot(as(x, "Spatial"),  ...)
	for (i in 1:n)
		points(x@lines[[i]]@Lines[[1]]@coords, pch = pch[i], col = col[i], cex = cex[i])
}
setMethod("plot", signature(x = "SpatialMultiPoints", y = "missing"),
    function(x, y, ...) plot.SpatialMultiPoints(x, ...))


###################################################
### chunk number 35: 
###################################################
cName <- "SpatialMultiPointsDataFrame"
setClass(cName, representation("SpatialLinesDataFrame"), 
	validity <- function(object) {
                lst <- lapply(object@lines, function(x) length(x@Lines))
		if (any(unlist(lst) != 1))
			stop("Only Lines objects with single Line")
		TRUE
	}
)
SpatialMultiPointsDataFrame <- function(object) {
   new("SpatialMultiPointsDataFrame", object)
}


###################################################
### chunk number 36: 
###################################################
df <- data.frame(x1 = 1:3, x2 = c(1,4,2), row.names = c("mp1", "mp2", "mp3"))
smp_df <- SpatialMultiPointsDataFrame(SpatialLinesDataFrame(smp, df))
setMethod("plot", signature(x = "SpatialMultiPointsDataFrame", y = "missing"),
    function(x, y, ...) plot.SpatialMultiPoints(x, ...))
grys <- c("grey10", "grey40", "grey80")


###################################################
### chunk number 38: 
###################################################
plot(smp_df, col = grys[smp_df[["x1"]]], pch = smp_df[["x2"]], cex = 2,
 axes = TRUE) 


###################################################
### chunk number 39: 
###################################################
data(meuse.grid)
gridded(meuse.grid)=~x+y
xx <- spsample(meuse.grid,  type="hexagonal", cellsize=200)
class(xx)


###################################################
### chunk number 40: 
###################################################
HexPts <- spsample(meuse.grid,  type="hexagonal", cellsize=200)


###################################################
### chunk number 42: 
###################################################
HexPols <- HexPoints2SpatialPolygons(HexPts)
if (packageVersion("sp") < "1.1.0") {
  df <- as.data.frame(meuse.grid)[overlay(meuse.grid, HexPts),]
} else {
   df <- as.data.frame(meuse.grid)[over(HexPts, as(meuse.grid, "SpatialGrid")),]
}
HexPolsDf <- SpatialPolygonsDataFrame(HexPols, df, match.ID = FALSE)


###################################################
### chunk number 44: 
###################################################
grys <- grey.colors(11, 0.95, 0.55, 2.2)
print(spplot(meuse.grid["dist"], cuts=10, col.regions=grys, sp.layout = list("sp.points", HexPts, col = 1)),
	split = c(1, 1, 2, 1), more = TRUE)
print(spplot(HexPolsDf["dist"], cuts=10, col.regions=grys),
	split = c(2, 1, 2, 1), more = FALSE)


###################################################
### chunk number 45: 
###################################################
setClass("SpatialHexGrid", representation("SpatialPoints", dx = "numeric"),
	validity <- function(object) {
		if (object@dx <= 0)
			stop("dx should be positive")
		TRUE
	}
)


###################################################
### chunk number 47: 
###################################################
setClass("SpatialHexGridDataFrame", representation("SpatialPointsDataFrame", dx = "numeric"),
# NOT TOO WIDE
	validity <- function(object) {
		if (object@dx <= 0)
			stop("dx should be positive")
		TRUE
	}
)


###################################################
### chunk number 49: 
###################################################
HexPts <- spsample(meuse.grid,  type="hexagonal", cellsize=200)
Hex <- new("SpatialHexGrid", HexPts, dx = 200)
if (packageVersion("sp") < "1.1.0") {
  df <- as.data.frame(meuse.grid)[overlay(meuse.grid, Hex),]
spdf <- SpatialPointsDataFrame(HexPts, df)
HexDf <- new("SpatialHexGridDataFrame", spdf, dx = 200)
}

###################################################
### chunk number 50: 
###################################################
if (packageVersion("sp") < "1.1.0") {
is(HexDf, "SpatialHexGrid")
setIs("SpatialHexGridDataFrame", "SpatialHexGrid")
is(HexDf, "SpatialHexGrid")
}

###################################################
### chunk number 52: 
###################################################
setAs("SpatialHexGrid", "SpatialPolygons", 
	function(from) 
		HexPoints2SpatialPolygons(from, from@dx)
)
setAs("SpatialHexGridDataFrame", "SpatialPolygonsDataFrame", 
	function(from)
		SpatialPolygonsDataFrame(as(obj, "SpatialPolygons"), obj@data, match.ID = FALSE)
)


###################################################
### chunk number 54: 
###################################################
setMethod("plot", signature(x = "SpatialHexGrid", y = "missing"),
    function(x, y, ...) plot(as(x, "SpatialPolygons"), ...)
)
setMethod("spplot", signature(obj = "SpatialHexGridDataFrame"),
    function(obj, ...)
		spplot(SpatialPolygonsDataFrame( as(obj, "SpatialPolygons"), obj@data, match.ID = FALSE), ...)
)
setMethod("spsample", "SpatialHexGrid", function(x, n, type, ...) 
	spsample(as(x, "SpatialPolygons"), n = n, type = type, ...)
)
if (packageVersion("sp") < "1.1.0") {
setMethod("overlay", c("SpatialHexGrid", "SpatialPoints"), function(x, y, ...) 
	overlay(as(x, "SpatialPolygons"), y)
)
}

###################################################
### chunk number 57: 
###################################################
bbox(Hex)
bbox(as(Hex, "SpatialPolygons"))


###################################################
### chunk number 58: 
###################################################
n <- 10
x <- data.frame(expand.grid(x1 = 1:n, x2 = 1:n, x3 = 1:n), z = rnorm(n^3))
coordinates(x) <- ~x1+x2+x3
gridded(x) <- TRUE
fullgrid(x) <- TRUE
summary(x)


###################################################
### chunk number 60: 
###################################################
setClass("SpatialTimeGrid", "SpatialGrid",
	validity <- function(object) {
		stopifnot(dimensions(object) == 3)
		TRUE
	}
)


###################################################
### chunk number 62: 
###################################################
setClass("SpatialTimeGridDataFrame", "SpatialGridDataFrame",
	validity <- function(object) {
		stopifnot(dimensions(object) == 3)
		TRUE
	}
)
#setIs("SpatialTimeGridDataFrame", "SpatialTimeGrid")
x <- new("SpatialTimeGridDataFrame", x)


###################################################
### chunk number 63: 
###################################################
summary.SpatialTimeGridDataFrame <- function(object, ...) {
	cat("Object of class SpatialTimeGridDataFrame\n")
	x <- gridparameters(object)
	t0 <- ISOdate(1970,1,1,0,0,0)
	t1 <- t0 + x[3,1]
	cat(paste("first time step:", t1, "\n"))
	t2 <- t0 + x[3,1] + (x[3,3] - 1) * x[3,2]
	cat(paste("last time step: ", t2, "\n"))
	cat(paste("time step:      ", x[3,2], "\n"))
	summary(as(object, "SpatialGridDataFrame"))
}


###################################################
### chunk number 65: 
###################################################
setMethod("summary", "SpatialTimeGridDataFrame", summary.SpatialTimeGridDataFrame)
summary(x)


###################################################
### chunk number 67: 
###################################################
subs.SpatialTimeGridDataFrame <- function(x, i, j, ..., drop=FALSE) {
	t <- coordinates(x)[,3] + ISOdate(1970,1,1,0,0,0)
	if (missing(j))
		j <- TRUE
	sel <- t %in% i
	if (! any(sel))
		stop("selection results in empty set")
	fullgrid(x) <- FALSE
	if (length(i) > 1) {
		x <- x[i = sel, j = j,...]
		fullgrid(x) <- TRUE
		as(x, "SpatialTimeGridDataFrame")
	} else {
		gridded(x) <- FALSE
		x <- x[i = sel, j = j,...]
		cc <- coordinates(x)[,1:2]
		p4s <- CRS(proj4string(x))
		SpatialPixelsDataFrame(cc, x@data, proj4string = p4s)
	}
}
setMethod("[", c("SpatialTimeGridDataFrame", "POSIXct", "ANY"), 
	subs.SpatialTimeGridDataFrame)
t1 <- as.POSIXct("1970-01-01 0:00:03", tz = "GMT")
t2 <- as.POSIXct("1970-01-01 0:00:05", tz = "GMT")
summary(x[c(t1,t2)])
summary(x[t1])


###################################################
### chunk number 68: 
###################################################
spplot.stgdf <- function(obj, zcol = 1, ..., format = NULL) {
	if (length(zcol) != 1)
		stop("can only plot a single attribute")
	if (is.null(format)) format <- "%Y-%m-%d %H:%M:%S"
	cc <- coordinates(obj)
	df <- unstack(data.frame(obj[[zcol]], cc[,3]))
	ns <- as.character(coordinatevalues(getGridTopology(obj))[[3]] + ISOdate(1970,1,1,0,0,0), format = format)
	cc2d <- cc[cc[,3] == min(cc[,3]), 1:2]
	obj <- SpatialPixelsDataFrame(cc2d, df)
	spplot(obj, names.attr = ns,...)
}
setMethod("spplot", "SpatialTimeGridDataFrame", spplot.stgdf)


###################################################
### chunk number 69: 
###################################################
print(spplot(x, format = "%H:%M:%S", as.table=TRUE))


###################################################
### chunk number 72: 
###################################################
library(gstat)
data(meuse)
coordinates(meuse) <- ~x+y
v <- vgm(.5, "Sph", 800, .05)
sim <- krige(log(zinc)~1, meuse, meuse.grid, v, nsim=100, nmax=30)
sim@data <- exp(sim@data)


###################################################
### chunk number 73: 
###################################################
quantile.Spatial <- function(x, ..., byLayer = FALSE) {
	stopifnot("data" %in% slotNames(x))
	apply(x@data, ifelse(byLayer, 2, 1), quantile, ...)
}


###################################################
### chunk number 74: 
###################################################
sim$lower <- quantile.Spatial(sim[1:100], probs = 0.025)
sim$upper <- quantile.Spatial(sim[1:100], probs = 0.975)


###################################################
### chunk number 75: 
###################################################
medians <- quantile.Spatial(sim[1:100], probs = 0.5, byLayer = TRUE)


###################################################
### chunk number 78: 
###################################################
fractionBelow <- function(x, q, byLayer = FALSE) {
	stopifnot(is(x, "Spatial") || !("data" %in% slotNames(x)))
	apply(x@data < q, ifelse(byLayer, 2, 1), function(r) sum(r)/length(r))
}


###################################################
### chunk number 80: 
###################################################
over500 <- 1 - fractionBelow(sim[1:100], 200, byLayer = TRUE)
summary(over500)
quantile(over500, c(0.025, 0.975))


###################################################
### chunk number 81:  
###################################################
fname <- unzip(zipfile = "70042108.zip", files="70042108.tif")
#fname <- zip.file.extract(file="70042108.tif", zipname = "70042108.zip")
#file.copy(fname, "./70042108.tif", overwrite=TRUE)
library(rgdal)
#x <- readGDAL("70042108.tif", output.dim = c(120, 132))
x <- readGDAL(fname, output.dim = c(120, 132))
x$band1[x$band1 <= 0] <- NA
spplot(x, col.regions=bpy.colors())


###################################################
### chunk number 82: 
###################################################
x <- GDAL.open("70042108.tif")
class(x)
x.subs <- x[1:100, 1:100, 1]
class(x.subs)
gridparameters(x.subs)


###################################################
### chunk number 84: 
###################################################
setClass("SpatialGDAL",
    representation("Spatial", grid = "GridTopology", grod = "GDALReadOnlyDataset", 
# NOT TOO WIDE
		name = "character"))
setClass("SpatialGDALWrite", "SpatialGDAL")


###################################################
### chunk number 86:  !!NOTE!! Takes much time with bla <- 20 !
###################################################
x <- open.SpatialGDAL("70042108.tif")
nrows <- GDALinfo("70042108.tif")["rows"]
ncols <- GDALinfo("70042108.tif")["columns"]
xout <- copy.SpatialGDAL(x, "70042108out.tif")
bls <- 200 
for (i in 1:(nrows/bls - 1)) {
	r <- 1+(i-1)*bls
	for (j in 1:(ncols/bls - 1)) {
		c <- 1+(j-1)*bls
		x.in <- x[r:(r+bls),c:(c+bls)]
		xout[r:(r+bls),c:(c+bls)] <- x.in$band1 + 10 #$
	}
	cat(paste("row-block", i, "\n"))
}
close(x)
close(xout)



