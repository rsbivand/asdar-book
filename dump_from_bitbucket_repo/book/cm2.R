###################################################
### chunk number 1: 
###################################################
rm(list=ls())
library(digest)
if (!exists("chkDigest")) chkDigest <- TRUE
if (!exists("online")) online <- TRUE
intamap <- "http://intamap.geo.uu.nl/~roger/ASDAR/data"
op <- options()
options("width"=70, warn=1, str = strOptions(strict.width="wrap", vec.len=2))
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
## .epsNo <- .epsNo + 1; file <- paste("Fig-cm2-", .epsNo, ".eps", sep="")
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
## .epsNo <- .epsNo + 1; file <- paste("Fig-cm2-", .epsNo, ".png", sep="")
## png(filename=file, width = .iwidth, height = .iheight, pointsize = .ipointsize)


###################################################
### chunk number 7: zfig_png eval=FALSE
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 8: 
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(1,1,1,1)+0.1, mfrow=c(1,2))
x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
grys <- grey.colors(11, 0.9, 0.45, 2.2)
# RSB quietening greys
image(y, x, t(volcano)[ncol(volcano):1,], breaks=seq(90,200,10), col=grys, asp=1, axes=FALSE)
contour(y, x, t(volcano)[ncol(volcano):1,], levels=seq(90,200,10), asp=1, axes=FALSE)
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 9: 
###################################################
#library(rgdal)
#auck_el1 <- readGDAL("/home/rsb/tmp/70042108/70042108.tif")
library(sp)
if (online) {
  con <- url(paste(intamap, "auck_el1.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("auck_el1.RData")
#if (chkDigest && !identical(digest(auck_el1), "d603173d5bf2d3e7bca177d0c673d569"))
  #stop("auck_el1.RData digest error")
transect_sp <- SpatialPoints(coords=cbind(seq(174.458,175.3,0.000833333),
  c(-37.03625)), proj4string=CRS("+proj=longlat +ellps=WGS84"))


###################################################
### chunk number 10: 
###################################################
summary(transect_sp)
transect_el1 <- overlay(auck_el1, transect_sp)
summary(transect_el1)


###################################################
### chunk number 11: 
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(4,4,1,1)+0.1, mfrow=c(2,1))
plot(coordinates(transect_el1)[,1], transect_el1$band1, type="l", main="", xlab="", ylab="elevation, m", axes=FALSE)
abline(h=0)
box()
axis(2)
axis(1, at=axTicks(1), labels=parse(text=sp:::degreeLabelsEW(axTicks(1))))
plot(ecdf(transect_el1$band1), verticals= TRUE, do.p = FALSE, main="", xlab="elevation, m", ylab="")
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 12: 
###################################################
#library(Rgshhs)
#auck2 <- Rgshhs("/home/rsb/tmp/GSHHS/gshhs_f.b", xlim=c(174.2,175.3), ylim=c(-37.5,-36.5), level=2)
#auck_gshhs <- auck2$SP
if (online) {
  con <- url(paste(intamap, "auck_gshhs.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("auck_gshhs.RData")
if (chkDigest && !identical(digest(auck_gshhs), "b5992b6958d92c5928127ec734a60bde"))
  stop("auck_gshhs.RData digest error")
#load("/home/rsb/tmp/70042108/auck_el2.RData")
if (online) {
  con <- url(paste(intamap, "auck_el2.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("auck_el2.RData")
#if (chkDigest && !identical(digest(summary(auck_el2)), "afe7a5812a09a763ba976908d02615de"))
  #stop("auck_el2.RData digest error")


###################################################
### chunk number 13: 
###################################################
set.seed(9876)
polygon_random <- spsample(auck_gshhs, 1000, type="random")
polygon_random_el1 <- overlay(auck_el1, polygon_random)
grid_random <- spsample(auck_el2, 1000, type="random")
grid_random_el1 <- overlay(auck_el1, grid_random)
grid_regular <- spsample(auck_el2, 1000, type="regular")
grid_regular_el1 <- overlay(auck_el1, grid_regular)


###################################################
### chunk number 14: 
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(2,1,2,1)+0.1, mfrow=c(2,2))
plot(auck_gshhs)
plot(polygon_random, add=TRUE, cex=0.15)
title(main="polygon_random")
plot(grid_random, cex=0.15)
title(main="grid_random")
plot(grid_regular, cex=0.15)
title(main="grid_regular")
plot(ecdf(transect_el1$band1), verticals=TRUE, do.p=FALSE, ann=FALSE, col.hor="grey", col.vert="grey")
title(main="ECDF")
plot(ecdf(polygon_random_el1$band1), verticals=TRUE, do.p=FALSE, add=TRUE)
plot(ecdf(grid_random_el1$band1), verticals=TRUE, do.p=FALSE, add=TRUE)
plot(ecdf(grid_regular_el1$band1), verticals=TRUE, do.p=FALSE, add=TRUE)
abline(h=c(0.25,0.5,0.75), lty=2, col="grey")
legend(c(350,650), c(0,0.2), c("transect", "samples"), lty=c(1,1), col=c("grey", "black"), bty="n")
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 15: 
###################################################
tab <- rbind(transect=c(fivenum(transect_el1$band1), nrow(coordinates(transect_el1))), polygon_random=c(fivenum(polygon_random_el1$band1), nrow(coordinates(polygon_random_el1))), grid_random=c(fivenum(grid_random_el1$band1), nrow(coordinates(grid_random_el1))), grid_regular=c(fivenum(grid_regular_el1$band1), nrow(coordinates(grid_regular_el1))))
colnames(tab) <- c("minimum", "lower-hinge", "median", "upper-hinge", "maximum", "n")
tab


###################################################
### chunk number 16:  eval=FALSE
###################################################
## library(rgdal)
## nc90 <- readOGR(".", "co37_d90")
## proj4string(nc90) <- CRS("+proj=longlat +datum=NAD27")
## sc90 <- readOGR(".", "co45_d90")
## proj4string(sc90) <- CRS("+proj=longlat +datum=NAD27")
## va90 <- readOGR(".", "co51_d90")
## proj4string(va90) <- CRS("+proj=longlat +datum=NAD27")


###################################################
### chunk number 17: 
###################################################
library(maptools)
if (online) {
  con <- url(paste(intamap, "nc90.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("nc90.RData")
if (chkDigest && !identical(digest(nc90), "f70a058491549b7c3c2e153e0007e72c"))
  stop("nc90.RData digest error")
if (online) {
  con <- url(paste(intamap, "sc90.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("sc90.RData")
if (chkDigest && !identical(digest(sc90), "aadccadd2950fcd524aee77ea93bd1e3"))
  stop("sc90.RData digest error")
if (online) {
  con <- url(paste(intamap, "va90.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("va90.RData")
if (chkDigest && !identical(digest(va90), "1f304208f606372ab4b0d9be22cbeb65"))
  stop("va90.RData digest error")


###################################################
### chunk number 18: 
###################################################
library(maptools)


###################################################
### chunk number 19: 
###################################################
names(sc90)
sc90a <- spChFIDs(sc90, paste(sc90$ST, sc90$CO, sep=""))
sc90a <- sc90a[,-(1:4)]
names(sc90a)


###################################################
### chunk number 20: 
###################################################
proj4string(sc90a) <- CRS(proj4string(sc90a))


###################################################
### chunk number 21: 
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(3,2,1,1)+0.1)
plot(va90, xlim=c(-85,-75), ylim=c(32,40), axes=TRUE, border="grey10")
plot(nc90, add=TRUE, border="grey40")
plot(sc90, add=TRUE, border="grey70")
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 22: 
###################################################
names(nc90)


###################################################
### chunk number 23:  eval=FALSE
###################################################
## nc90a <- spChFIDs(nc90, paste(nc90$ST, nc90$CO, sep=""))


###################################################
### chunk number 24: 
###################################################
try1 <- try(spChFIDs(nc90, paste(nc90$ST, nc90$CO, sep="")))


###################################################
### chunk number 25: 
###################################################
cat(try1)


###################################################
### chunk number 26: 
###################################################
table(table(paste(nc90$ST, nc90$CO, sep="")))


###################################################
### chunk number 27: 
###################################################
nc90a <- unionSpatialPolygons(nc90, IDs=paste(nc90$ST, nc90$CO, sep=""))


###################################################
### chunk number 28: 
###################################################
nc90_df <- as(nc90, "data.frame")[!duplicated(nc90$CO),-(1:4)]
row.names(nc90_df) <- paste(nc90_df$ST, nc90_df$CO, sep="")
nc90b <- SpatialPolygonsDataFrame(nc90a, nc90_df)


###################################################
### chunk number 29: 
###################################################
va90a <- spChFIDs(va90, paste(va90$ST, va90$CO, sep=""))
va90a <- va90a[,-(1:4)]
va90_pl <- slot(va90a, "polygons")
va90_pla <- lapply(va90_pl, checkPolygonsHoles)
p4sva <- CRS(proj4string(va90a))
vaSP <- SpatialPolygons(va90_pla, proj4string=p4sva)
va90b <- SpatialPolygonsDataFrame(vaSP, data=as(va90a, "data.frame"))


###################################################
### chunk number 30: 
###################################################
nc_sc_va90 <- spRbind(spRbind(nc90b, sc90a), va90b)
FIPS <- sapply(slot(nc_sc_va90, "polygons"), function(x) slot(x, "ID"))
str(FIPS)
length(slot(nc_sc_va90, "polygons"))


###################################################
### chunk number 31: 
###################################################
if (online) {
  con <- url(paste(intamap, "90mfips.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("90mfips.RData")
if (chkDigest && !identical(digest(t1), "27be5c453b102366f4b814077c5a2f65"))
  stop("90mfips.RData digest error")


###################################################
### chunk number 32:  eval=FALSE
###################################################
## t1 <- read.fwf("90mfips.txt", skip=21,
##  widths=c(4,4,4,4,2,6,2,3,3,1,7,5,3,51), colClasses = "character")


###################################################
### chunk number 33: 
###################################################
t2 <- t1[1:2004,c(1,7,8,14)]
t3 <- t2[complete.cases(t2),]
cnty1 <- t3[t3$V7 != "  ",]
ma1 <- t3[t3$V7 == "  ",c(1,4)]
cnty2 <- cnty1[which(!is.na(match(cnty1$V7, c("37", "45", "51")))),]
cnty2$FIPS <- paste(cnty2$V7, cnty2$V8, sep="")


###################################################
### chunk number 34: 
###################################################
MA_FIPS <- cnty2$V1[match(FIPS, cnty2$FIPS)]
MA <- ma1$V14[match(MA_FIPS, ma1$V1)]
MA_df <- data.frame(MA_FIPS=MA_FIPS, MA=MA, row.names=FIPS)
nc_sc_va90a <- spCbind(nc_sc_va90, MA_df)
ncscva_MA <- unionSpatialPolygons(nc_sc_va90a, nc_sc_va90a$MA_FIPS)


###################################################
### chunk number 35: 
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-cm2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(3,2,1,1)+0.1)
plot(nc_sc_va90, border="grey", axes=TRUE)
plot(ncscva_MA, lwd=2, add=TRUE)
text(coordinates(ncscva_MA), labels=sapply(slot(ncscva_MA, "polygons"), function(x) slot(x, "ID")), cex=0.6)
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 36: 
###################################################
np <- sapply(slot(ncscva_MA, "polygons"), function(x) length(slot(x, "Polygons")))
table(np)
MA_fips <- sapply(slot(ncscva_MA, "polygons"), function(x) slot(x, "ID"))
MA_name <- ma1$V14[match(MA_fips, ma1$V1)]
data.frame(MA_fips, MA_name)[np > 1,]


###################################################
### chunk number 37: 
###################################################
hels <- matrix(c(24.97, 60.17), nrow=1)
p4s <- CRS("+proj=longlat +datum=WGS84")
Hels <- SpatialPoints(hels, proj4string=p4s)
d041224 <- as.POSIXct("2004-12-24", tz="EET")
sunriset(Hels, d041224, direction="sunrise", POSIXct.out=TRUE)


###################################################
### chunk number 38: 
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
sT <- capture.output(print(Sys.time()))
cat("\n")
cat(paste("%", sT, sep=" "), sep="\n")


###################################################
### chunk number 39: 
###################################################
options(op)


