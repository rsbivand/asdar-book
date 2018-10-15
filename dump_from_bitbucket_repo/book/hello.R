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
## .iwidth <- 4
## .iheight <- 4
## .ipointsize <- 10
## .pwd <- 0.95
## grey_gamma <- 2.2


###################################################
### chunk number 3: 
###################################################
.iwidth <- 4
.iheight <- 4
.ipointsize <- 10
.pwd <- 0.95
grey_gamma <- 2.2


###################################################
### chunk number 4: afig eval=FALSE
###################################################
## .epsNo <- .epsNo + 1; file <- paste("Fig-hlo-", .epsNo, ".eps", sep="")
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
### chunk number 6: zfigr eval=FALSE
###################################################
## par(def.par)
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth,angle=90]{", file, "}", sep="")
## cat("\n")


###################################################
### chunk number 7: afig_l eval=FALSE
###################################################
## .epsNo <- .epsNo + 1; file <- paste("Fig-hlo-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)


###################################################
### chunk number 8: zfig_l eval=FALSE
###################################################
## dev.null <- dev.off()
## #system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
## cat("\n")


###################################################
### chunk number 9: 
###################################################
library("pkgDepTools")
# library("Biobase")
library("Rgraphviz")
allDeps <- makeDepGraph("http://cran.ii.uib.no", type="source", keep.builtin=FALSE, dosize=FALSE)
allDepsOnMe <- reverseEdgeDirections(allDeps)
categoryNodes <- c("sp", names(acc(allDepsOnMe, "sp")[[1]]))
categoryGraph <- subGraph(categoryNodes, allDepsOnMe)
nn <- makeNodeAttrs(categoryGraph, shape="ellipse")
fc <- rep("#ffffff", numNodes(categoryGraph))
names(fc) <- names(nn$label)
ours <- match(c("sp", "DCluster", "gstat", "maptools", "rgdal", "spdep", "spgrass6", "spgwr", "splancs"), nn$label)
fc[ours] <- "#e0e0e0"
nn$fillcolor <- fc


###################################################
### chunk number 10: 
###################################################
.iwidth <- 7
.iheight <- 4
.ipointsize <- 10
.pwd <- 1.4
.epsNo <- .epsNo + 1; file <- paste("Fig-hlo-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
plot(categoryGraph, nodeAttrs=nn)
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth,angle=90]{", file, "}", sep="")
cat("\n")
.iwidth <- 4
.iheight <- 4
.ipointsize <- 10
.pwd <- 0.95
grey_gamma <- 2.2


###################################################
### chunk number 11: 
###################################################
library(maptools)
library(maps)
library(rgdal)


###################################################
### chunk number 12: 
###################################################
.iwidth <- 6
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-hlo-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
def.par <- par(no.readonly = TRUE)
# volc.tab = read.table("data.xy")
volc.tab = read.table("hsd_data/data1964al.xy")
volc = SpatialPoints(volc.tab[c(2,1)])
llCRS <- CRS("+proj=longlat +ellps=WGS84")
proj4string(volc) <- llCRS
prj_new = CRS("+proj=moll")
volc_proj = spTransform(volc, prj_new)
 wrld <- map("world", interior=FALSE, xlim=c(-179,179), ylim=c(-89,89), plot=FALSE)
 wrld_p <- pruneMap(wrld, xlim=c(-179,179))
 wrld_sp <- map2SpatialLines(wrld_p, proj4string=llCRS)
 wrld_proj <- spTransform(wrld_sp, prj_new)
 #save(c("wrld_proj", "wrld_sp"), file = "hsd_data/wrld.RData")
 #load("hsd_data/wrld.RData")
wrld_grd <- gridlines(wrld_sp, easts=c(-179,seq(-150,150,50),179.5), norths=seq(-75,75,15), ndiscr=100)
wrld_grd_proj <- spTransform(wrld_grd, prj_new)
at_sp <- gridat(wrld_sp, easts=0, norths=seq(-75,75,15), offset=0.3)
at_proj <- spTransform(at_sp, prj_new)
opar = par(no.readonly = TRUE)
par(mar=c(1,1,1,1)+0.1, xpd=NA)
plot(wrld_proj, col="grey50")
plot(wrld_grd_proj, add=TRUE, lty=3, col="grey50")
points(volc_proj, cex = .8, pch = 3)
text(coordinates(at_proj), pos=at_proj$pos, offset=at_proj$offset, labels=parse(text=as.character(at_proj$lab)), cex=0.6)
#legend(c(-18000000,-11000000), c(-6000000,-2500000), pch=3, cex = .2, legend=c("volcano"), bty="n") 
par(opar)
# $
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
.iwidth <- 4
.iheight <- 4
.ipointsize <- 10
.pwd <- 0.95
grey_gamma <- 2.2


###################################################
### chunk number 13: 
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 9
.epsNo <- .epsNo + 1; file <- paste("Fig-hlo-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
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
#image(volcano, axes=F, col='white', asp=1, main="d")
image(volcano, axes=F, col=grys, asp=1, main="d")
x3l1 = coordinates(x2[x2$level == 160,])[[1]][[1]]
x3l2 = coordinates(x2[x2$level == 160,])[[1]][[2]]
#plot(SpatialPolygons(list(Polygons(list(Polygon(x3l1,hole=F), Polygon(x3l2,hole=T)), ID=c("x")))),
#	pbg = 'white', col = grey(.5), add = T)
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
# $
par(def.par)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", file, "}", sep="")
cat("\n")
layout(matrix(1))
.iwidth <- 4
.iheight <- 4
.ipointsize <- 10
.pwd <- 0.95
grey_gamma <- 2.2


###################################################
### chunk number 14: 
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
sT <- capture.output(print(Sys.time()))
cat("%", sT, "\n")


