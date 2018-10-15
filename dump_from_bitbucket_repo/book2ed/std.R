### R code from vignette source 'std.Rnw'

###################################################
### code chunk number 1: std.Rnw:9-14
###################################################
if (!exists("book_R_dont_trash")) rm(list=ls())
op <- options()
options("width"=70, warn=1, str = strOptions(strict.width="wrap", vec.len=2), 
	useFancyQuotes="TeX")
.epsNo <- 0


###################################################
### code chunk number 2: figreset (eval = FALSE)
###################################################
## .iwidth <- 5
## .iheight <- 6
## .ipointsize <- 12


###################################################
### code chunk number 3: std.Rnw:22-23
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 4: afig (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-std-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
## lattice.options(default.theme = col.whitebg())


###################################################
### code chunk number 5: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 6: afig_png (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-std-", .epsNo, ".png", sep="")
## png(filename=file, width = .iwidth, height = .iheight, pointsize = .ipointsize)
## lattice.options(default.theme = col.whitebg())


###################################################
### code chunk number 7: zfig_png (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 8: afig_l (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-std-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
## .iheight, pointsize = .ipointsize, horizontal=FALSE)
## lattice.options(default.theme = col.whitebg())


###################################################
### code chunk number 9: zfig_l (eval = FALSE)
###################################################
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
## cat("\n")


###################################################
### code chunk number 10: std.Rnw:63-67
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
sT <- capture.output(print(Sys.time()))
cat("%", sT, "\n")


###################################################
### code chunk number 11: std.Rnw:408-412
###################################################
ecd.ll <- as.matrix(read.table("ECDovelatlon.dat", header = FALSE))
library(sp)
ecd.ll <- SpatialPoints(ecd.ll[,c(2,1)])
proj4string(ecd.ll) <- CRS("+proj=longlat +datum=WGS84")


###################################################
### code chunk number 12: std.Rnw:421-425
###################################################
library(xts)
library(spacetime)
ecd.years <- 1986:2003
ecd.y <- as.Date(paste(ecd.years, "-01-01", sep=""), "%Y-%m-%d")


###################################################
### code chunk number 13: std.Rnw:434-438
###################################################
ecd <- read.table("ECDoveBBS1986_2003.dat", header=FALSE)
ecd[ecd == -1] <- NA
ecd.st <- STFDF(ecd.ll, ecd.y, data.frame(counts = as.vector(as.matrix(ecd))))
dim(ecd.st)


###################################################
### code chunk number 14: std.Rnw:453-456
###################################################
ecd.st2 <- stConstruct(ecd, ecd.ll, list(counts = names(ecd)),
    TimeObj = ecd.y, interval = TRUE)
all.equal(ecd.st2, ecd.st)


###################################################
### code chunk number 15: std.Rnw:503-508
###################################################
library(maps)
m <- map("state", "florida", fill = TRUE, plot = FALSE)
library(maptools)
FL <- map2SpatialPolygons(m, "FL")
proj4string(FL) <- proj4string(ecd.st)


###################################################
### code chunk number 16: std.Rnw:518-522
###################################################
dim(ecd.st[FL,])
dim(ecd.st[, "1998::2003"])
dim(ecd.st[,,"counts"])
dim(ecd.st[FL, "1998::2003", "counts"])


###################################################
### code chunk number 17: std.Rnw:531-535
###################################################
mode(ecd.st[[1]])
length(ecd.st[[1]])
length(ecd.st[["counts"]])
length(ecd.st$counts)


###################################################
### code chunk number 18: std.Rnw:543-544
###################################################
ecd.st$sqrtcounts <- sqrt(ecd.st$counts)


###################################################
### code chunk number 19: std.Rnw:562-563 (eval = FALSE)
###################################################
## over(x, y)


###################################################
### code chunk number 20: std.Rnw:578-579
###################################################
bb <- STF(FL, ecd.y[c(4,6,8,10,12)])


###################################################
### code chunk number 21: std.Rnw:587-588
###################################################
over(bb, ecd.st) #, fn=sum, na.rm=TRUE)


###################################################
### code chunk number 22: std.Rnw:600-601
###################################################
over(bb, ecd.st, fn=sum, na.rm=TRUE)


###################################################
### code chunk number 23: std.Rnw:609-610
###################################################
bb.counts <- new("STFDF", bb, data = over(bb, ecd.st, fn=sum, na.rm=TRUE))


###################################################
### code chunk number 24: std.Rnw:619-620
###################################################
aggregate(ecd.st, bb, sum, na.rm = TRUE)


###################################################
### code chunk number 25: std.Rnw:655-656
###################################################
ecd.5y <- aggregate(ecd.st, "5 years", mean, na.rm = TRUE)


###################################################
### code chunk number 26: std.Rnw:665-666 (eval = FALSE)
###################################################
## vignette("sto")


###################################################
### code chunk number 27: std.Rnw:713-725
###################################################
.iwidth <- 8
.iheight <- 6
.pwd <- .7
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-std-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = 
.iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(RColorBrewer)
print(stplot(ecd.5y[FL,], c("1986-1990", "1991-1995", "1996-2000", "2001-2003"),
	col.regions = brewer.pal(6, "Reds"),
	cex=1, 
	cuts=c(0,5,10,20,40,80,131), sp.layout = list("sp.polygons", FL, col = "grey"),
	ylim = bbox(FL)[2,], scales = list(draw=FALSE),colorkey=TRUE))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
cat("\n")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 28: std.Rnw:775-779
###################################################
ecd.FL <- ecd.st[FL, , "sqrtcounts"]
x <- as(ecd.FL, "xts")
x[is.na(x)] <- 0
o <- order(as.vector(1:18 %*% x) / apply(x,2,sum))


###################################################
### code chunk number 29: std.Rnw:784-790 (eval = FALSE)
###################################################
## library(RColorBrewer)
## pal <- brewer.pal(6, "Reds")
## cuts <- c(0,2,4,6,8,10,12)
## ck <- list(at = cuts, labels = as.character(cuts^2))
## stplot(ecd.FL[o,], mode = "xt", col.regions = pal, cuts = 6, asp = .5, 
## 	xlab = "Sites, ordered by time", colorkey = ck)


###################################################
### code chunk number 30: std.Rnw:802-819
###################################################
.iwidth <- 8
.iheight <- 4
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-std-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
ecd.FL <- ecd.st[FL,,"sqrtcounts"]
x <- as(ecd.FL, "xts")
x[is.na(x)] <- 0
o <- order(as.vector(1:18 %*% x) / apply(x,2,sum))
library(RColorBrewer)
pal <- brewer.pal(6, "Reds")
cuts <- c(0,2,4,6,8,10,12)
ck <- list(at = cuts, labels = as.character(cuts^2))
stplot(ecd.FL[o,], mode = "xt", col.regions = pal, cuts = 6, asp = .5, 
	xlab = "Sites, ordered by time", 
	scales=list(x=list(rot = 90, cex=.5)),
	colorkey = ck)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 31: std.Rnw:829-848
###################################################
.iwidth <- 8
.iheight <- 4
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-std-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(maps)
states.m = map("state", plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(states.m$names, ":"), function(x) x[1])
library(maptools)
states = map2SpatialPolygons(states.m, IDs=IDs)
yrs = 1970:1986
time = as.POSIXct(paste(yrs, "-01-01", sep=""), tz = "GMT")
library(plm)
data("Produc")
Produc.st = STFDF(states[-8], time, Produc[order(Produc[2], Produc[1]),])
print(stplot(Produc.st[1:2,,5:8], mode = "tp", key.space = "bottom"),
	more = TRUE, split = c(1,1,2,1))
print(stplot(Produc.st[c(1:3,5),,5:6], mode = "ts", key.space = "bottom"),
	more = FALSE, split = c(2,1,2,1))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 32: std.Rnw:969-977
###################################################
sI <- toLatex(sessionInfo())
cat(paste("%", sI), sep="\n")
cat("\n")
ver <- system("svnversion", intern=TRUE)
cat("%SVN version", ver, "\n")
sT <- capture.output(print(Sys.time()))
cat("\n")
cat(paste("%", sT, sep=" "), sep="\n")


###################################################
### code chunk number 33: std.Rnw:980-981
###################################################
options(op)


