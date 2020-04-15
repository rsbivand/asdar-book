###################################################
### chunk number 1: 
###################################################
rm(list=ls())
library(digest)
.owidth <- getOption("width")
options("width"=70)
owarn <- options("warn")$warn
options(warn=1)
.epsNo <- 0


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
## .epsNo <- .epsNo + 1; file <- paste("Fig-die-", .epsNo, ".eps", sep="")
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
## .epsNo <- .epsNo + 1; file <- paste("Fig-die-", .epsNo, ".png", sep="")
## png(filename=file, width = .iwidth, height = .iheight, pointsize = .ipointsize)


###################################################
### chunk number 7: zfig_png eval=FALSE
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 8: 
###################################################
library(rgdal)


###################################################
### chunk number 9: 
###################################################
gdl <- Sys.getenv("GDAL_DATA")
pl <- Sys.getenv("PROJ_LIB")
Smess <- paste('Geospatial Data Abstraction Library ',
    'extensions to R successfully loaded\n',
    'Loaded GDAL runtime: ', getGDALVersionInfo(), '\n',
    paste(names(gdl[1]), ": ", gdl[1], sep=""), "\n",
    'Loaded PROJ.4 runtime: ', getPROJ4VersionInfo(), '\n',
    paste(names(pl[1]), ": ", pl[1], sep=""), "\n", sep="")
cat(Smess)


###################################################
### chunk number 10: 
###################################################
EPSG <- make_EPSG()
EPSG[grep("^# ED50$", EPSG$note),]


###################################################
### chunk number 11: 
###################################################
CRS("+init=epsg:4230")


###################################################
### chunk number 12: 
###################################################
ED50 <- CRS("+init=epsg:4230 +towgs84=-87,-96,-120,0,0,0,0")
ED50


###################################################
### chunk number 13: 
###################################################
IJ.east <- as(char2dms("4d31\'00\"E"), "numeric")
IJ.north <- as(char2dms("52d28\'00\"N"), "numeric")
IJ.ED50 <- SpatialPoints(cbind(x=IJ.east, y=IJ.north), ED50)
res <- spTransform(IJ.ED50, CRS("+proj=longlat +datum=WGS84"))
x <- as(dd2dms(coordinates(res)[1]), "character")
y <- as(dd2dms(coordinates(res)[2], TRUE), "character")
cat(x, y, "\n")
spDistsN1(coordinates(IJ.ED50), coordinates(res), longlat=TRUE)*1000
library(maptools)
gzAzimuth(coordinates(IJ.ED50), coordinates(res))


###################################################
### chunk number 14: 
###################################################
EPSG[grep("Atlas", EPSG$note), 1:2]


###################################################
### chunk number 15:  eval=FALSE
###################################################
## CRS("+init=epsg:2163")


###################################################
### chunk number 16: 
###################################################
cat(strwrap(CRSargs(CRS("+init=epsg:2163")), exdent=4), sep="\n")


###################################################
### chunk number 17: 
###################################################
proj <- projInfo("proj")
proj[proj$name == "laea",]
ellps <- projInfo("ellps")
ellps[grep("a=6370997", ellps$major),]


###################################################
### chunk number 18: 
###################################################
IJ.dms.E <- "4d31\'00\"E"
IJ.dms.N <- "52d28\'00\"N"


###################################################
### chunk number 19: 
###################################################
IJ_east <- char2dms(IJ.dms.E)
IJ_north <- char2dms(IJ.dms.N)
IJ_east
IJ_north
getSlots("DMS")


###################################################
### chunk number 20: 
###################################################
c(as(IJ_east, "numeric"), as(IJ_north, "numeric"))


###################################################
### chunk number 21: 
###################################################
download.file("http://www.sph.emory.edu/~lwaller/book/ch9/scot.shp", "scot.shp", mode="wb")
download.file("http://www.sph.emory.edu/~lwaller/book/ch9/scot.dbf", "scot.dbf", mode="wb")
download.file("http://www.sph.emory.edu/~lwaller/book/ch9/scot.shx", "scot.shx", mode="wb")
download.file("http://www.sph.emory.edu/~lwaller/book/ch2/scotland.dat", "scotland.dat", mode="w")


###################################################
### chunk number 22: 
###################################################
scot_LL <- readOGR(".", "scot")
proj4string(scot_LL) <- CRS("+proj=longlat ellps=WGS84")
scot_LL$ID


###################################################
### chunk number 23: 
###################################################
scot_dat <- read.table("scotland.dat", skip=1)
names(scot_dat) <- c("District", "Observed", "Expected", "PcAFF", "Latitude", "Longitude")
scot_dat$District
library(maptools)
scot_dat1 <- scot_dat[match(scot_LL$ID, scot_dat$District),]
row.names(scot_dat1) <- sapply(slot(scot_LL, "polygons"), function(x) slot(x, "ID"))
scot_LLa <- spCbind(scot_LL, scot_dat1)
all.equal(scot_LLa$ID, scot_LLa$District)
names(scot_LLa)


###################################################
### chunk number 24: 
###################################################
library(spdep)
O <- scot_LLa$Observed
E <- scot_LLa$Expected
scot_LLa$SMR <- probmap(O, E)$relRisk/100
library(DCluster)
scot_LLa$smth <- empbaysmooth(O, E)$smthrr


###################################################
### chunk number 25: 
###################################################
scot_BNG <- spTransform(scot_LLa, CRS("+init=epsg:27700"))


###################################################
### chunk number 26: 
###################################################
library(lattice)
trellis.par.set(canonical.theme(color = FALSE))
library(RColorBrewer)


###################################################
### chunk number 27: 
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
spplot(scot_BNG, c("SMR", "smth"),
 at=c(0, 0.5, 1, 1.5, 2.5, 7),
 col.regions=grey.colors(5, 0.95, 0.55, 2.2))
# RSB quietening greys
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 28: 
###################################################
writeOGR(scot_LLa["ID"], dsn="scot_district.kml", layer="borders",
 driver="KML")
llCRS <- CRS("+proj=longlat ellps=WGS84")
scot_SP_LL <- SpatialPointsDataFrame(coordinates(scot_LLa),
 proj4string=llCRS, data=as(scot_LLa, "data.frame")[c("NAME",
 "Observed", "Expected", "SMR", "smth")])
writeOGR(scot_SP_LL, dsn="scot_rates.kml", layer="rates", driver="KML")


###################################################
### chunk number 29: 
###################################################
drv <- "ESRI Shapefile"
writeOGR(scot_BNG, dsn=".", layer="scot_BNG", driver=drv)


###################################################
### chunk number 30: 
###################################################
list.files(pattern="^scot_BNG")


###################################################
### chunk number 31: 
###################################################
getinfo.shape("scot_BNG.shp")


###################################################
### chunk number 32: 
###################################################
download.file("http://intamap.geo.uu.nl/~roger/ASDAR/data/70042108.zip", "70042108.zip", mode="wb")
fname <- zip.file.extract(file="70042108.tif", zipname = "70042108.zip")
file.copy(fname, "./70042108.tif", overwrite=TRUE)


###################################################
### chunk number 33: 
###################################################
auck_el1 <- readGDAL("70042108.tif")
summary(auck_el1)
is.na(auck_el1$band1) <- auck_el1$band1 <= 0 | auck_el1$band1 > 1e+4


###################################################
### chunk number 34: 
###################################################
x <- GDAL.open("70042108.tif")
xx <- getDriver(x)
xx
getDriverLongName(xx)
x
dim(x)
GDAL.close(x)


###################################################
### chunk number 35: 
###################################################
GDALinfo("70042108.tif")


###################################################
### chunk number 36: 
###################################################
library(maptools)
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) = TRUE
proj4string(meuse.grid) <- CRS(paste("+init=epsg:28992", "+towgs84=565.237,50.0087,465.658,-0.406857,0.350733,-1.87035,4.0812"))
data(meuse)
coordinates(meuse) <- c("x", "y")
proj4string(meuse) <- CRS(proj4string(meuse.grid))


###################################################
### chunk number 37: 
###################################################
library(gstat)
log_zinc <- krige(log(zinc)~1, meuse, meuse.grid)["var1.pred"]
proj4string(log_zinc) <- CRS(proj4string(meuse.grid))
summary(log_zinc)


###################################################
### chunk number 38: 
###################################################
writeGDAL(log_zinc, fname="log_zinc.tif", driver="GTiff", type="Float32", options="INTERLEAVE=PIXEL")


###################################################
### chunk number 39:  eval=FALSE
###################################################
## GDALinfo("log_zinc.tif")


###################################################
### chunk number 40: 
###################################################
options(width=55)
GDALinfo("log_zinc.tif")
options(width=70)


###################################################
### chunk number 41: 
###################################################
library(maptools)
grd <- as(meuse.grid, "SpatialPolygons")
proj4string(grd) <- CRS(proj4string(meuse))
grd.union <- unionSpatialPolygons(grd, rep("x", length(slot(grd, "polygons"))))
ll <- CRS("+proj=longlat +datum=WGS84")
grd.union.ll <- spTransform(grd.union, ll)


###################################################
### chunk number 42: 
###################################################
llGRD <- GE_SpatialGrid(grd.union.ll)
llGRD_in <- overlay(llGRD$SG, grd.union.ll)
llSGDF <- SpatialGridDataFrame(grid=slot(llGRD$SG, "grid"), proj4string=CRS(proj4string(llGRD$SG)), data=data.frame(in0=llGRD_in))
llSPix <- as(llSGDF, "SpatialPixelsDataFrame")


###################################################
### chunk number 43: 
###################################################
meuse_ll <- spTransform(meuse, CRS("+proj=longlat +datum=WGS84"))
llSPix$pred <- idw(log(zinc)~1, meuse_ll, llSPix)$var1.pred


###################################################
### chunk number 44: 
###################################################
png(file="zinc_IDW.png", width=llGRD$width, height=llGRD$height, bg="transparent")
par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
image(llSPix, "pred", col=bpy.colors(20))
dev.off()
kmlOverlay(llGRD, "zinc_IDW.kml", "zinc_IDW.png")


###################################################
### chunk number 45: 
###################################################
system("g.version", intern=TRUE)
library(spgrass6)
gmeta6()


###################################################
### chunk number 46: 
###################################################
spear <- readRAST6(c("elevation.dem", "geology"), cat=c(FALSE, TRUE))


###################################################
### chunk number 47: 
###################################################
summary(spear)


###################################################
### chunk number 48: 
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(4,4,2,2)+0.1)
plot(ecdf(spear$elevation.dem), pch=".", main="")
par(oopar)
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 49: 
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(4,7,2,2)+0.1, las=1)
tg <- table(spear$geology)
boxplot(elevation.dem ~ geology, spear, width=tg, horizontal=TRUE)
par(oopar)
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 50: 
###################################################
table(spear$geology)
system("r.stats --q -cl geology", intern=TRUE)


###################################################
### chunk number 51: 
###################################################
bugsDF <- readVECT6("bugsites")


###################################################
### chunk number 52: 
###################################################
vInfo("streams")


###################################################
### chunk number 53: 
###################################################
streams <- readVECT6("streams", type="line,boundary", remove.duplicates=FALSE)


###################################################
### chunk number 54: 
###################################################
summary(bugsDF)


###################################################
### chunk number 55: 
###################################################
system("g.gisenv set='LOCATION_NAME=snow2'")
system("g.region rast=snowcost_broad")


###################################################
### chunk number 56: 
###################################################
buildings <- readVECT6("vsnow4")
sohoSG <- readRAST6(c("snowcost_broad", "snowcost_not_broad"))


###################################################
### chunk number 57: 
###################################################
deaths <- readVECT6("deaths3")
o <- overlay(sohoSG, deaths)
deaths <- spCbind(deaths, as(o, "data.frame"))
deaths$b_nearer <- deaths$snowcost_broad < deaths$snowcost_not_broad


###################################################
### chunk number 58: 
###################################################
by(deaths$Num_Cases, deaths$b_nearer, sum)


###################################################
### chunk number 59: 
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1,2), mar=c(5,3,1,1)+0.1)
b_wid <- table(deaths$b_nearer)
boxplot(snowcost_broad ~ b_nearer, deaths, width=b_wid, ylim=c(0,450), ylab="distance", xlab="Broad Street", col=grey.colors(1, 0.8, 0.8, 2.2))
boxplot(snowcost_not_broad ~ b_nearer, deaths, width=b_wid, ylim=c(0,450), xlab="Other pump", col=grey.colors(1, 0.8, 0.8, 2.2))
par(oopar)
# RSB quietening greys
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 60: 
###################################################
nb_pump <- readVECT6("vpump_not_broad")
b_pump <- readVECT6("vpump_broad")


###################################################
### chunk number 61: 
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(1,1,1,1)+0.1)
#image(sohoSG, "snowcost_broad", breaks=seq(0,750,50),
# col=rev(heat.colors(15)))
library(RColorBrewer)
#gcols <- colorRampPalette(brewer.pal(8, "Greys")[-1])(15)
# RSB quietening greys
gcols <- grey.colors(15, 0.95, 0.55, 2.2)
image(sohoSG, "snowcost_broad", breaks=seq(0,750,50),
 col=gcols)
plot(buildings, col="white", add=TRUE)
plot(buildings, angle=45, density=10, col="grey70", add=TRUE)
#symbols(coordinates(deaths), circles=4*sqrt(deaths$Num_Cases),
# inches=FALSE, add=TRUE, bg=c("brown2","grey40")[deaths$b_nearer+1])
symbols(coordinates(deaths), circles=4*sqrt(deaths$Num_Cases),
 inches=FALSE, add=TRUE, bg=c("grey75","grey50")[deaths$b_nearer+1])
source("legend_image.R") #from geoR
rect(528900, 180550, 529040, 180990, border=NA, col="white")
text(528970, 180950, "metres from\nBroad Street\npump", cex=0.6)
#legend_image(c(528930, 528960), c(180600, 180900),
# sohoSG$snowcost_broad, vertical=TRUE, breaks=seq(0,750,50),
# col=rev(heat.colors(15)))
legend_image(c(528930, 528960), c(180600, 180900),
 sohoSG$snowcost_broad, vertical=TRUE, breaks=seq(0,750,50),
 col=gcols)
plot(nb_pump, add=TRUE, pch=8, cex=1.3, lwd=2)
plot(b_pump, add=TRUE, pch=4, cex=1.5, lwd=8, col="white")
plot(b_pump, add=TRUE, pch=4, cex=1.5, lwd=6)
#plot(b_pump, add=TRUE, pch=4, cex=1.3, lwd=2, col="plum1")
rect(528900, 181330, 529140, 181380, border=NA, col="white")
legend(c(528910, 529100), c(181350, 181380),
 legend=c("Broad Street pump","other pumps"), pch=c(4,8), bty="n",
 cex=0.6, y.inter=0.7)
rect(528900, 181270, 529180, 181335, border=NA, col="white")
#legend(c(528910, 529100), c(181275, 181325),
# legend=c("nearer Broad Street pump","nearer other pump"),
# fill=c("grey40","brown2"), bty="n", cex=0.6, y.inter=0.7)
legend(c(528910, 529100), c(181275, 181325),
 legend=c("nearer Broad Street pump","nearer other pump"),
 fill=c("grey50","grey75"), bty="n", cex=0.6, y.inter=0.7)
box()
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 62: 
###################################################
system("g.gisenv set='LOCATION_NAME=spearfish60'")


###################################################
### chunk number 63: 
###################################################
sp2Mondrian(scot_BNG, "scot_BNG.txt")


###################################################
### chunk number 64: 
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
sT <- capture.output(print(Sys.time()))
cat("\n")
cat("%", sT, "\n")
cat("\n% ")
cat(system("g.version", intern=TRUE), "\n")


###################################################
### chunk number 65: 
###################################################
options(warn=owarn)


