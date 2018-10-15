### R code from vignette source 'die.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: die.Rnw:2-7
###################################################
if (!exists("book_R_dont_trash")) rm(list=ls())
op <- options()
options("width"=70, warn=1, str = strOptions(strict.width="wrap", vec.len=2), useFancyQuotes="TeX")
.epsNo <- 0
library(lattice)


###################################################
### code chunk number 2: figreset (eval = FALSE)
###################################################
## .iwidth <- 5
## .iheight <- 6
## .ipointsize <- 12


###################################################
### code chunk number 3: die.Rnw:15-16
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 4: afig (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)


###################################################
### code chunk number 5: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 6: afig_png (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".png", sep="")
## png(filename=file, width = .iwidth, height = .iheight, pointsize = .ipointsize)


###################################################
### code chunk number 7: zfig_png (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 8: die.Rnw:38-40
###################################################
ret <- system("g.version")
stopifnot(ret == 0L)


###################################################
### code chunk number 9: die.Rnw:95-96
###################################################
library(rgdal)


###################################################
### code chunk number 10: die.Rnw:98-121
###################################################
#gdl <- Sys.getenv("GDAL_DATA")
#pl <- Sys.getenv("PROJ_LIB")
    gdl <- getGDAL_DATA_Path()
    pl <- getPROJ4libPath()
    if (nchar(pl) == 0) 
        pl <- "(autodetected)"

    fn <- system.file("SVN_VERSION", package = "rgdal")
    if (file.exists(fn)) {
        svn_version <- scan(fn, what = character(1), sep = "\n", 
            quiet = TRUE)
    } else {
        svn_version <- "(unknown)"
    }
Smess <- paste("rgdal: version: ", utils::packageDescription("rgdal")$Version, 
        ", (SVN revision ", svn_version, ")\n", 
    'Geospatial Data Abstraction Library ',
    'extensions to R successfully loaded\n',
    'Loaded GDAL runtime: ', getGDALVersionInfo(), '\n',
    paste("Path to GDAL shared files: ", gdl[1], sep=""), "\n",
    'Loaded PROJ.4 runtime: ', getPROJ4VersionInfo(), '\n',
    paste("Path to PROJ.4 shared files: ", pl[1], sep=""), "\n", sep="")
cat(Smess)


###################################################
### code chunk number 11: die.Rnw:209-211 (eval = FALSE)
###################################################
## NEWS <- "http://svn.osgeo.org/metacrs/proj/trunk/proj/NEWS"
## PROJ4_NEWS <- readLines(url(NEWS))


###################################################
### code chunk number 12: die.Rnw:213-215
###################################################
#save(PROJ4_NEWS, file="PROJ4_NEWS.RData")
load("PROJ4_NEWS.RData")


###################################################
### code chunk number 13: die.Rnw:217-219
###################################################
lns <- grep("Release Notes|EPSG", PROJ4_NEWS)
head(PROJ4_NEWS[lns])


###################################################
### code chunk number 14: die.Rnw:264-266
###################################################
EPSG <- make_EPSG()
EPSG[grep("^# ED50$", EPSG$note),]


###################################################
### code chunk number 15: die.Rnw:306-307
###################################################
CRS("+init=epsg:4230")


###################################################
### code chunk number 16: die.Rnw:330-332
###################################################
ED50 <- CRS("+init=epsg:4230 +towgs84=-87,-96,-120,0,0,0,0")
ED50


###################################################
### code chunk number 17: die.Rnw:358-368
###################################################
IJ.east <- as(char2dms("4d31\'00\"E"), "numeric")
IJ.north <- as(char2dms("52d28\'00\"N"), "numeric")
IJ.ED50 <- SpatialPoints(cbind(x=IJ.east, y=IJ.north), proj4string=ED50)
res <- spTransform(IJ.ED50, CRS("+proj=longlat +datum=WGS84"))
x <- as(dd2dms(coordinates(res)[1]), "character")
y <- as(dd2dms(coordinates(res)[2], TRUE), "character")
cat(x, y, "\n")
spDistsN1(coordinates(IJ.ED50), coordinates(res), longlat=TRUE)*1000
library(maptools)
gzAzimuth(coordinates(IJ.ED50), coordinates(res))


###################################################
### code chunk number 18: die.Rnw:427-428
###################################################
owarn <- options(warn=-1)


###################################################
### code chunk number 19: die.Rnw:430-434
###################################################
proj4string(IJ.ED50) <- CRS("+init=epsg:4230")
res <- spTransform(IJ.ED50, CRS("+proj=longlat +datum=WGS84"))
spDistsN1(coordinates(IJ.ED50), coordinates(res), longlat=TRUE)*1000
gzAzimuth(coordinates(IJ.ED50), coordinates(res))


###################################################
### code chunk number 20: die.Rnw:436-437
###################################################
options(warn=owarn$warn)


###################################################
### code chunk number 21: die.Rnw:449-450
###################################################
EPSG[grep("Atlas", EPSG$note), 1:2]


###################################################
### code chunk number 22: die.Rnw:452-453 (eval = FALSE)
###################################################
## CRS("+init=epsg:2163")


###################################################
### code chunk number 23: die.Rnw:455-456
###################################################
cat(strwrap(CRSargs(CRS("+init=epsg:2163")), exdent=4), sep="\n")


###################################################
### code chunk number 24: die.Rnw:475-479
###################################################
proj <- projInfo("proj")
proj[proj$name == "laea",]
ellps <- projInfo("ellps")
ellps[grep("a=6370997", ellps$major),]


###################################################
### code chunk number 25: die.Rnw:521-523
###################################################
IJ.dms.E <- "4d31\'00\"E"
IJ.dms.N <- "52d28\'00\"N"


###################################################
### code chunk number 26: die.Rnw:533-538
###################################################
IJ_east <- char2dms(IJ.dms.E)
IJ_north <- char2dms(IJ.dms.N)
IJ_east
IJ_north
getSlots("DMS")


###################################################
### code chunk number 27: die.Rnw:552-553
###################################################
c(as(IJ_east, "numeric"), as(IJ_north, "numeric"))


###################################################
### code chunk number 28: die.Rnw:672-673
###################################################
head(ogrDrivers(), n=10)


###################################################
### code chunk number 29: die.Rnw:696-697 (eval = FALSE)
###################################################
## vignette("OGR_shape_encoding", package="rgdal")


###################################################
### code chunk number 30: die.Rnw:730-731
###################################################
setwd("../Data")


###################################################
### code chunk number 31: die.Rnw:733-735
###################################################
scot_dat <- read.table("scotland.dat", skip=1)
names(scot_dat) <- c("District", "Observed", "Expected", "PcAFF", "Latitude", "Longitude")


###################################################
### code chunk number 32: die.Rnw:737-738
###################################################
setwd("../die")


###################################################
### code chunk number 33: die.Rnw:757-758
###################################################
setwd("../Data")


###################################################
### code chunk number 34: die.Rnw:760-761
###################################################
ogrInfo(".", "scot")


###################################################
### code chunk number 35: die.Rnw:763-764
###################################################
setwd("../die")


###################################################
### code chunk number 36: die.Rnw:788-789
###################################################
setwd("../Data")


###################################################
### code chunk number 37: die.Rnw:791-794
###################################################
scot_LL <- readOGR(dsn=".", layer="scot")
proj4string(scot_LL)
proj4string(scot_LL) <- CRS("+proj=longlat +ellps=WGS84")


###################################################
### code chunk number 38: die.Rnw:796-797
###################################################
setwd("../die")


###################################################
### code chunk number 39: die.Rnw:803-805
###################################################
sapply(slot(scot_LL, "data"), class)
scot_LL$ID


###################################################
### code chunk number 40: die.Rnw:817-825
###################################################
scot_dat$District
ID_D <- match(scot_LL$ID, scot_dat$District)
scot_dat1 <- scot_dat[ID_D,]
row.names(scot_dat1) <- row.names(scot_LL)
library(maptools)
scot_LLa <- spCbind(scot_LL, scot_dat1)
all.equal(scot_LLa$ID, scot_LLa$District)
names(scot_LLa)


###################################################
### code chunk number 41: die.Rnw:846-852
###################################################
library(spdep)
O <- scot_LLa$Observed
E <- scot_LLa$Expected
scot_LLa$SMR <- probmap(O, E)$relRisk/100
library(DCluster)
scot_LLa$smth <- empbaysmooth(O, E)$smthrr


###################################################
### code chunk number 42: die.Rnw:863-864
###################################################
scot_BNG <- spTransform(scot_LLa, CRS("+init=epsg:27700"))


###################################################
### code chunk number 43: die.Rnw:868-869
###################################################
library(RColorBrewer)


###################################################
### code chunk number 44: die.Rnw:874-887
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
trellis.par.set(canonical.theme("postscript", color=TRUE))
spplot(scot_BNG, c("SMR", "smth"),
 at=c(0, 0.25, 0.5, 0.8, 1, 1.5, 2.5, 4.5, 7),
 col.regions=rev(brewer.pal(8, "RdBu")))
# at=c(0, 0.5, 1, 1.5, 2.5, 7),
# col.regions=grey.colors(5, 0.95, 0.55, 2.2))
# RSB quietening greys
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 45: die.Rnw:907-909
###################################################
drv <- "ESRI Shapefile"
writeOGR(scot_BNG, dsn=".", layer="scot_BNG", driver=drv)


###################################################
### code chunk number 46: die.Rnw:911-912
###################################################
list.files(pattern="^scot_BNG")


###################################################
### code chunk number 47: die.Rnw:941-942
###################################################
load("geohub.RData")


###################################################
### code chunk number 48: die.Rnw:944-946 (eval = FALSE)
###################################################
## dsn <- "WFS:http://geohub.jrc.ec.europa.eu/effis/ows"
## ogrListLayers(dsn)


###################################################
### code chunk number 49: die.Rnw:949-950
###################################################
print(layers)


###################################################
### code chunk number 50: die.Rnw:952-953 (eval = FALSE)
###################################################
## Fires <- readOGR(dsn, "EFFIS:FiresAll")


###################################################
### code chunk number 51: die.Rnw:955-956
###################################################
cat(paste(strwrap(geohub, exdent=5), collapse="\n"), "\n")


###################################################
### code chunk number 52: die.Rnw:958-959
###################################################
setwd("../Data")


###################################################
### code chunk number 53: die.Rnw:961-962
###################################################
Fires <- readOGR(".", "fires_120104")


###################################################
### code chunk number 54: die.Rnw:964-965
###################################################
setwd("../die")


###################################################
### code chunk number 55: die.Rnw:967-968
###################################################
names(Fires)


###################################################
### code chunk number 56: die.Rnw:986-996
###################################################
x <- c(-15, -15, 38, 38, -15)
y <- c(28, 62, 62, 28, 28)
crds <- cbind(x=x, y=y)
bb <- SpatialPolygons(list(Polygons(list(Polygon(coords=crds)), "1")))
library(maptools)
data(wrld_simpl)
proj4string(bb) <- CRS(proj4string(wrld_simpl))
library(rgeos)
slbb <- gIntersection(bb, as(wrld_simpl, "SpatialLines"))
spl <- list("sp.lines", slbb, lwd=0.7, col="khaki4")


###################################################
### code chunk number 57: die.Rnw:1013-1018
###################################################
Fires$dt <- as.Date(as.character(Fires$FireDate), format="%d-%m-%Y")
Fires0 <- Fires[-which(coordinates(Fires)[,2] < 0),]
Fires1 <- Fires0[order(Fires0$dt),]
library(spacetime)
Fires2 <- STIDF(as(Fires1, "SpatialPoints"), Fires1$dt, as(Fires1, "data.frame"))


###################################################
### code chunk number 58: die.Rnw:1020-1021 (eval = FALSE)
###################################################
## stplot(Fires2, number=3, sp.layout=spl, cex=0.5)


###################################################
### code chunk number 59: die.Rnw:1028-1040
###################################################
.iwidth <- 5
.iheight <- 4.5
.ipointsize <- 8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
th <- canonical.theme("postscript", color=TRUE)
th$plot.symbol$col <- "darkorange3"
th$plot.symbol$pch <- 16L
th$plot.symbol$cex <- 0.5
trellis.par.set(th)
stplot(Fires2, number=3, sp.layout=spl)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 60: die.Rnw:1057-1060
###################################################
names(Fires1)[1] <- "name"
GR_Fires <- Fires1[Fires1$Country == "GR",]
writeOGR(GR_Fires, "EFFIS.gpx", "waypoints", driver="GPX", dataset_options="GPX_USE_EXTENSIONS=YES")


###################################################
### code chunk number 61: die.Rnw:1071-1073
###################################################
GR <- readOGR("EFFIS.gpx", "waypoints")
GR[1,c(5,24:28)]


###################################################
### code chunk number 62: die.Rnw:1075-1076
###################################################
unlink("EFFIS.gpx")


###################################################
### code chunk number 63: die.Rnw:1107-1108
###################################################
getinfo.shape("scot_BNG.shp")


###################################################
### code chunk number 64: die.Rnw:1110-1111
###################################################
unlink("scot_BNG.*")


###################################################
### code chunk number 65: die.Rnw:1175-1176
###################################################
setwd("../Data")


###################################################
### code chunk number 66: die.Rnw:1178-1181
###################################################
auck_el1 <- readGDAL("70042108.tif")
summary(auck_el1)
is.na(auck_el1$band1) <- auck_el1$band1 <= 0 | auck_el1$band1 > 1e+4


###################################################
### code chunk number 67: die.Rnw:1183-1184
###################################################
setwd("../die")


###################################################
### code chunk number 68: die.Rnw:1201-1202
###################################################
setwd("../Data")


###################################################
### code chunk number 69: die.Rnw:1204-1211
###################################################
x <- GDAL.open("70042108.tif")
xx <- getDriver(x)
xx
getDriverLongName(xx)
x
dim(x)
GDAL.close(x)


###################################################
### code chunk number 70: die.Rnw:1213-1214
###################################################
setwd("../die")


###################################################
### code chunk number 71: die.Rnw:1233-1234
###################################################
setwd("../Data")


###################################################
### code chunk number 72: die.Rnw:1236-1237
###################################################
GDALinfo("70042108.tif")


###################################################
### code chunk number 73: die.Rnw:1239-1240
###################################################
setwd("../die")


###################################################
### code chunk number 74: die.Rnw:1263-1272
###################################################
brks <- c(0,10,20,50,100,150,200,300,400,500,600,700)
pal <- terrain.colors(11)
pal
length(pal) == length(brks)-1
auck_el1$band1 <- findInterval(auck_el1$band1, vec=brks, all.inside=TRUE)-1
writeGDAL(auck_el1, "demIndex.tif", drivername="GTiff", type="Byte", colorTable=list(pal), mvFlag=length(brks)-1)
Gi <- GDALinfo("demIndex.tif", returnColorTable=TRUE)
CT <- attr(Gi, "ColorTable")[[1]]
CT[CT > "#000000"]


###################################################
### code chunk number 75: die.Rnw:1274-1275
###################################################
unlink("demIndex.tif")


###################################################
### code chunk number 76: die.Rnw:1296-1303
###################################################
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) <- TRUE
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
data(meuse)
coordinates(meuse) <- c("x", "y")
proj4string(meuse) <- CRS(proj4string(meuse.grid))


###################################################
### code chunk number 77: die.Rnw:1308-1310
###################################################
library(gstat)
log_zinc <- idw(log(zinc)~1, meuse, meuse.grid)["var1.pred"]


###################################################
### code chunk number 78: die.Rnw:1312-1314
###################################################
summary(log_zinc)
writeGDAL(log_zinc, fname="log_zinc.tif", drivername="GTiff", type="Float32", options="INTERLEAVE=PIXEL")


###################################################
### code chunk number 79: die.Rnw:1316-1317 (eval = FALSE)
###################################################
## GDALinfo("log_zinc.tif")


###################################################
### code chunk number 80: die.Rnw:1319-1322
###################################################
options(width=55)
GDALinfo("log_zinc.tif")
options(width=70)


###################################################
### code chunk number 81: die.Rnw:1324-1325
###################################################
unlink("log_zinc.tif")


###################################################
### code chunk number 82: die.Rnw:1346-1354
###################################################
Soil <- meuse.grid["soil"]
table(Soil$soil)
Soil$soil <- as.integer(Soil$soil)-1
Cn <- c("Rd10A", "Rd90C/VII", "Bkd26/VII")
writeGDAL(Soil, "Soil.tif", drivername="GTiff", type="Byte", catNames=list(Cn), mvFlag=length(Cn))
Gi <- GDALinfo("Soil.tif", returnCategoryNames=TRUE)
attr(Gi, "CATlist")[[1]]
summary(readGDAL("Soil.tif"))


###################################################
### code chunk number 83: die.Rnw:1356-1357
###################################################
unlink("demIndex.tif")


###################################################
### code chunk number 84: die.Rnw:1368-1369
###################################################
head(gdalDrivers(), n=10)


###################################################
### code chunk number 85: die.Rnw:1378-1379
###################################################
writeGDAL(log_zinc, fname="log_zinc.rda", drivername="R")


###################################################
### code chunk number 86: die.Rnw:1381-1382 (eval = FALSE)
###################################################
## GDALinfo("log_zinc.rda")


###################################################
### code chunk number 87: die.Rnw:1384-1387
###################################################
options(width=55)
GDALinfo("log_zinc.rda")
options(width=70)


###################################################
### code chunk number 88: die.Rnw:1389-1390
###################################################
unlink("log_zinc.rda*")


###################################################
### code chunk number 89: die.Rnw:1407-1408
###################################################
setwd("../Data")


###################################################
### code chunk number 90: die.Rnw:1410-1413 (eval = FALSE)
###################################################
## service_xml <- "frmt_wms_openstreetmap_tms.xml"
## offset <- c(19339000, 34546000)
## osm <- readGDAL(service_xml, offset=offset, region.dim=c(2000, 2000), output.dim=c(1000, 1000))


###################################################
### code chunk number 91: die.Rnw:1415-1416
###################################################
setwd("../die")


###################################################
### code chunk number 92: die.Rnw:1418-1420
###################################################
load("getosm.RData")
cat(paste(strwrap(getosm, exdent=5), collapse="\n"), "\n")


###################################################
### code chunk number 93: die.Rnw:1422-1423
###################################################
setwd("../Data")


###################################################
### code chunk number 94: die.Rnw:1425-1426
###################################################
osm <- readGDAL("osm_bergen_120105.tif")


###################################################
### code chunk number 95: die.Rnw:1428-1429
###################################################
setwd("../die")


###################################################
### code chunk number 96: die.Rnw:1431-1432
###################################################
summary(osm)


###################################################
### code chunk number 97: die.Rnw:1441-1450
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(1,1,1,1))
image(osm, red=1, green=2, blue=2)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 98: die.Rnw:1497-1498
###################################################
load("RgoogleMaps.RData")


###################################################
### code chunk number 99: die.Rnw:1501-1503 (eval = FALSE)
###################################################
## library(RgoogleMaps)
## myMap <- GetMap(center=c(60.395, 5.322), zoom =16, destfile = "MyTile2.png", maptype = "mobile")


###################################################
### code chunk number 100: die.Rnw:1505-1513
###################################################
BB <- do.call("rbind", myMap$BBOX)
dBB <- rev(diff(BB))
DIM12 <- dim(myMap$myTile)[1:2]
cs <- dBB/DIM12
cc <- c(BB[1,2] + cs[1]/2, BB[1,1] + cs[2]/2)
GT <- GridTopology(cc, cs, DIM12)
p4s <- CRS("+proj=longlat +datum=WGS84")
SG_myMap <- SpatialGridDataFrame(GT, proj4string=p4s, data=data.frame(r=c(t(myMap$myTile[,,1]))*255, g=c(t(myMap$myTile[,,2]))*255, b=c(t(myMap$myTile[,,3]))*255))


###################################################
### code chunk number 101: die.Rnw:1515-1516 (eval = FALSE)
###################################################
## myMap1 <- GetMap.OSM(lonR = c(5.3190, 5.3280), latR = c(60.392, 60.398), scale=4000, destfile = "MyTile.png")


###################################################
### code chunk number 102: die.Rnw:1535-1537
###################################################
library(osmar)
load("osmar.RData")


###################################################
### code chunk number 103: die.Rnw:1540-1544 (eval = FALSE)
###################################################
## library(osmar)
## api <- osmsource_api()
## box <- corner_bbox(5.3190, 60.392, 5.3280, 60.398)
## torget <- get_osm(box, source = api)


###################################################
### code chunk number 104: die.Rnw:1546-1547
###################################################
torget1 <- as_sp(torget, "lines")


###################################################
### code chunk number 105: die.Rnw:1548-1549 (eval = FALSE)
###################################################
## sort(table(torget1$user), decreasing=TRUE)[1:3]


###################################################
### code chunk number 106: die.Rnw:1551-1554
###################################################
mappers3 <- sort(table(torget1$user), decreasing=TRUE)[1:3]
names(mappers3) <- iconv(names(mappers3), "UTF-8", "UTF-8")
mappers3


###################################################
### code chunk number 107: die.Rnw:1572-1576
###################################################
bybane <- find(torget, way(tags(k == "light_rail")))
bybane <- find_down(torget, way(bybane))
bybane <- subset(torget, ids=bybane)
bybane <- as_sp(bybane, "lines")


###################################################
### code chunk number 108: die.Rnw:1583-1596
###################################################
.iwidth <- 5
.iheight <- 2.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(1,1,1,1), mfrow=c(1,2))
image(SG_myMap, red=1, green=2, blue=3)
plot(torget1, add=TRUE)
plot(bybane, add=TRUE, lwd=5, col="orange2")
plot(0:1, 0:1, type = "n", axes = FALSE, asp=1)
rasterImage(myMap1[[4]], 0, 0, 1, 1)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 109: die.Rnw:1611-1612 (eval = FALSE)
###################################################
## writeOGR(Fires[,c("gml_id", "FireDate", "Area_HA")], dsn="fires.kml", layer="fires", driver="KML")


###################################################
### code chunk number 110: die.Rnw:1614-1615
###################################################
writeOGR(Fires[,c("gml_id", "FireDate", "Area_HA")], dsn="fires.kml", layer="fires", driver="KML", overwrite_layer=TRUE)


###################################################
### code chunk number 111: die.Rnw:1673-1679
###################################################
library(maptools)
grd <- as(meuse.grid, "SpatialPolygons")
proj4string(grd) <- CRS(proj4string(meuse))
grd.union <- unionSpatialPolygons(grd, rep("x", length(slot(grd, "polygons"))))
ll <- CRS("+proj=longlat +datum=WGS84")
grd.union.ll <- spTransform(grd.union, ll)


###################################################
### code chunk number 112: die.Rnw:1703-1707
###################################################
llGRD <- GE_SpatialGrid(grd.union.ll)
llGRD_in <- over(llGRD$SG, grd.union.ll)
llSGDF <- SpatialGridDataFrame(grid=slot(llGRD$SG, "grid"), proj4string=CRS(proj4string(llGRD$SG)), data=data.frame(in0=llGRD_in))
llSPix <- as(llSGDF, "SpatialPixelsDataFrame")


###################################################
### code chunk number 113: die.Rnw:1722-1724
###################################################
meuse_ll <- spTransform(meuse, CRS("+proj=longlat +datum=WGS84"))
llSPix$pred <- gstat::idw(log(zinc)~1, meuse_ll, llSPix)$var1.pred


###################################################
### code chunk number 114: die.Rnw:1744-1749
###################################################
png(file="zinc_IDW.png", width=llGRD$width, height=llGRD$height, bg="transparent")
par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
image(llSPix, "pred", col=bpy.colors(20))
dev.off()
kmlOverlay(llGRD, "zinc_IDW.kml", "zinc_IDW.png")


###################################################
### code chunk number 115: die.Rnw:1880-1881
###################################################
library(spgrass6)


###################################################
### code chunk number 116: die.Rnw:1883-1884 (eval = FALSE)
###################################################
## execGRASS("g.region", flags="p")


###################################################
### code chunk number 117: die.Rnw:1886-1887
###################################################
cat(execGRASS("g.region", flags="p", intern=TRUE), sep="\n")


###################################################
### code chunk number 118: die.Rnw:1906-1907
###################################################
spear <- readRAST6(c("elevation.dem", "geology"), cat=c(FALSE, TRUE), close_OK=FALSE)


###################################################
### code chunk number 119: die.Rnw:1909-1910 (eval = FALSE)
###################################################
## spear <- readRAST6(c("elevation.dem", "geology"), cat=c(FALSE, TRUE))


###################################################
### code chunk number 120: die.Rnw:1912-1913
###################################################
summary(spear)


###################################################
### code chunk number 121: die.Rnw:1924-1933
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(4,4,2,2)+0.1)
plot(ecdf(spear$elevation.dem), pch=".", main="")
par(oopar)
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 122: die.Rnw:1942-1952
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(4,7,2,2)+0.1, las=1)
tg <- table(spear$geology)
boxplot(elevation.dem ~ geology, spear, width=tg, horizontal=TRUE)
par(oopar)
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 123: die.Rnw:1961-1962
###################################################
table(spear$geology)


###################################################
### code chunk number 124: die.Rnw:1963-1964 (eval = FALSE)
###################################################
## execGRASS("r.stats", input="geology", flags=c("quiet", "c", "l"))


###################################################
### code chunk number 125: die.Rnw:1966-1967
###################################################
cat(execGRASS("r.stats", input="geology", flags=c("quiet", "c", "l"), intern=TRUE), sep="\n")


###################################################
### code chunk number 126: die.Rnw:1998-1999
###################################################
bugsDF <- readVECT6("bugsites")


###################################################
### code chunk number 127: die.Rnw:2001-2002
###################################################
vInfo("streams")


###################################################
### code chunk number 128: die.Rnw:2004-2005
###################################################
streams <- readVECT6("streams", type="line,boundary", remove.duplicates=FALSE)


###################################################
### code chunk number 129: die.Rnw:2029-2030
###################################################
summary(bugsDF)


###################################################
### code chunk number 130: die.Rnw:2089-2091
###################################################
execGRASS("g.gisenv", set="LOCATION_NAME=snow2")
execGRASS("g.region", rast="snowcost_broad")


###################################################
### code chunk number 131: die.Rnw:2144-2145
###################################################
sohoSG <- readRAST6(c("snowcost_broad", "snowcost_not_broad"), close_OK=FALSE)


###################################################
### code chunk number 132: die.Rnw:2147-2148 (eval = FALSE)
###################################################
## sohoSG <- readRAST6(c("snowcost_broad", "snowcost_not_broad"))


###################################################
### code chunk number 133: die.Rnw:2149-2151
###################################################
buildings <- readVECT6("vsnow4")
proj4string(sohoSG) <- CRS(proj4string(buildings))


###################################################
### code chunk number 134: die.Rnw:2168-2173
###################################################
deaths <- readVECT6("deaths3")
o <- over(deaths, sohoSG)
library(maptools)
deaths <- spCbind(deaths, o)
deaths$b_nearer <- deaths$snowcost_broad < deaths$snowcost_not_broad


###################################################
### code chunk number 135: die.Rnw:2175-2176
###################################################
by(deaths$Num_Cases, deaths$b_nearer, sum)


###################################################
### code chunk number 136: die.Rnw:2182-2194
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1,2), mar=c(5,3,1,1)+0.1)
b_wid <- table(deaths$b_nearer)
boxplot(snowcost_broad ~ b_nearer, deaths, width=b_wid, ylim=c(0,450), ylab="distance", xlab="Broad Street", col=grey.colors(1, 0.8, 0.8, 2.2))
boxplot(snowcost_not_broad ~ b_nearer, deaths, width=b_wid, ylim=c(0,450), xlab="Other pump", col=grey.colors(1, 0.8, 0.8, 2.2))
par(oopar)
# RSB quietening greys
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 137: die.Rnw:2203-2205
###################################################
nb_pump <- readVECT6("vpump_not_broad")
b_pump <- readVECT6("vpump_broad")


###################################################
### code chunk number 138: die.Rnw:2229-2277
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-die-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(1,1,1,1)+0.1)
library(RColorBrewer)
image(sohoSG, "snowcost_broad", breaks=seq(0,750,50),
 col=colorRampPalette(brewer.pal(7, "Reds"))(15))
#gcols <- colorRampPalette(brewer.pal(8, "Greys")[-1])(15)
# RSB quietening greys
#gcols <- grey.colors(15, 0.95, 0.55, 2.2)
#image(sohoSG, "snowcost_broad", breaks=seq(0,750,50),
# col=gcols)
plot(buildings, col="white", add=TRUE)
plot(buildings, angle=45, density=10, col="grey70", add=TRUE)
symbols(coordinates(deaths), circles=4*sqrt(deaths$Num_Cases),
 inches=FALSE, add=TRUE, bg=c("brown2","grey40")[deaths$b_nearer+1])
#symbols(coordinates(deaths), circles=4*sqrt(deaths$Num_Cases),
# inches=FALSE, add=TRUE, bg=c("grey75","grey50")[deaths$b_nearer+1])
#source("legend_image.R") #from geoR
rect(528900, 180550, 529040, 180990, border=NA, col="white")
text(528970, 180950, "metres from\nBroad Street\npump", cex=0.6)
library(geoR)
legend.krige(c(528930, 528960), c(180600, 180900),
 sohoSG$snowcost_broad, vertical=TRUE, breaks=seq(0,750,50),
 col=colorRampPalette(brewer.pal(7, "Reds"))(15))
#legend.krige(c(528930, 528960), c(180600, 180900),
# sohoSG$snowcost_broad, vertical=TRUE, breaks=seq(0,750,50),
# col=gcols)
plot(nb_pump, add=TRUE, pch=8, cex=1.3, lwd=2)
plot(b_pump, add=TRUE, pch=4, cex=1.5, lwd=8, col="white")
plot(b_pump, add=TRUE, pch=4, cex=1.5, lwd=6)
#plot(b_pump, add=TRUE, pch=4, cex=1.3, lwd=2, col="plum1")
rect(528900, 181330, 529140, 181380, border=NA, col="white")
legend(c(528910, 529100), c(181350, 181380),
 legend=c("Broad Street pump","other pumps"), pch=c(4,8), bty="n",
 cex=0.6, y.inter=0.7)
rect(528900, 181270, 529180, 181335, border=NA, col="white")
legend(c(528910, 529100), c(181275, 181325),
 legend=c("nearer Broad Street pump","nearer other pump"),
 fill=c("grey40","brown2"), bty="n", cex=0.6, y.inter=0.7)
#legend(c(528910, 529100), c(181275, 181325),
# legend=c("nearer Broad Street pump","nearer other pump"),
# fill=c("grey50","grey75"), bty="n", cex=0.6, y.inter=0.7)
box()
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 139: die.Rnw:2292-2301 (eval = FALSE)
###################################################
## library(rgeos)
## vsnow4buf <- gBuffer(buildings, width=-4)
## GRD <- gmeta2grd()
## SG <- SpatialGrid(GRD, proj4string=CRS(proj4string(vsnow4buf)))
## o <- over(SG, vsnow4buf)
## crs <- CRS(proj4string(vsnow4buf)
## SGDF <- SpatialGridDataFrame(GRD, proj4string=crs), data=data.frame(o=o))
## SGDF$o[is.na(SGDF$o)] <- 2.5
## SGDF$o[SGDF$o == 1] <- NA


###################################################
### code chunk number 140: die.Rnw:2307-2312 (eval = FALSE)
###################################################
## library(gdistance)
## r <- as(SGDF, "RasterLayer")
## tr <- transition(r, mean, 8)
## d_b_pump <- rSPDistance(tr, deaths, b_pump, theta=1e-12)
## d_nb_pump <- rSPDistance(tr, deaths, nb_pump, theta=1e-12)


###################################################
### code chunk number 141: die.Rnw:2314-2315
###################################################
load("gdistance_snow.RData")


###################################################
### code chunk number 142: die.Rnw:2317-2321
###################################################
deaths$g_snowcost_broad <- d_b_pump[,1]
deaths$g_snowcost_not_broad <- apply(d_nb_pump, 1, min)
deaths$g_b_nearer <- deaths$g_snowcost_broad < deaths$g_snowcost_not_broad
by(deaths$Num_Cases, deaths$g_b_nearer, sum)


###################################################
### code chunk number 143: die.Rnw:2332-2333
###################################################
execGRASS("g.gisenv", set="LOCATION_NAME=spearfish60")


###################################################
### code chunk number 144: die.Rnw:2597-2605
###################################################
sI <- toLatex(sessionInfo())
cat(paste("%", sI), sep="\n")
cat("\n")
ver <- system("svnversion", intern=TRUE)
cat("%SVN version", ver, "\n")
cat("\n")
sT <- capture.output(print(Sys.time()))
cat(paste("%", sT, sep=" "), sep="\n")


###################################################
### code chunk number 145: die.Rnw:2608-2609
###################################################
options(op)


