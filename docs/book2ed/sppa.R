### R code from vignette source 'sppa.Rnw'

###################################################
### code chunk number 1: sppa.Rnw:6-10
###################################################
if (!exists("book_R_dont_trash")) rm(list=ls())
op <- options()
options("width"=70, warn=1, str = strOptions(strict.width="wrap", vec.len=2), useFancyQuotes="TeX")
.epsNo <- 0


###################################################
### code chunk number 2: figreset (eval = FALSE)
###################################################
## .iwidth <- 5
## .iheight <- 6
## .ipointsize <- 12
## .pwd <- 0.95


###################################################
### code chunk number 3: sppa.Rnw:19-20
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 4: afig (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
## lattice.options(default.theme = col.whitebg())


###################################################
### code chunk number 5: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 6: sppa.Rnw:230-232
###################################################
library(spatstat)
data(japanesepines)


###################################################
### code chunk number 7: sppa.Rnw:234-235
###################################################
summary(japanesepines)


###################################################
### code chunk number 8: sppa.Rnw:258-259
###################################################
library(maptools)


###################################################
### code chunk number 9: sppa.Rnw:261-263
###################################################
spjpines <- as(japanesepines, "SpatialPoints")
summary(spjpines)


###################################################
### code chunk number 10: sppa.Rnw:275-277
###################################################
spjpines1 <- elide(spjpines, scale=TRUE, unitsq=TRUE)
summary(spjpines1)


###################################################
### code chunk number 11: sppa.Rnw:294-296
###################################################
pppjap <- as(spjpines1, "ppp")
summary(pppjap)


###################################################
### code chunk number 12: sppa.Rnw:300-311
###################################################
data(redwoodfull)
spred <- as(redwoodfull, "SpatialPoints")
data(cells)
spcells <- as(cells, "SpatialPoints")
dpp<-data.frame(rbind(coordinates(spjpines1), coordinates(spred), 
   coordinates(spcells)))
njap<-nrow(coordinates(spjpines1))
nred<-nrow(coordinates(spred))
ncells<-nrow(coordinates(spcells))
dpp<-cbind(dpp,c(rep("JAPANESE",njap), rep("REDWOOD", nred), rep("CELLS", ncells))) 
names(dpp)<-c("x", "y", "DATASET")


###################################################
### code chunk number 13: sppa.Rnw:319-327
###################################################
.iwidth <- 6
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
print(xyplot(y~x|DATASET, data=dpp, pch=19, aspect=1))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 14: sppa.Rnw:398-399
###################################################
setwd("../Data")


###################################################
### code chunk number 15: sppa.Rnw:401-406
###################################################
library(rgdal)
spasthma <- readOGR(".", "spasthma")
spbdry <- readOGR(".", "spbdry")
spsrc <- readOGR(".", "spsrc")
sproads <- readOGR(".", "sproads")


###################################################
### code chunk number 16: sppa.Rnw:408-409
###################################################
setwd("../sppa")


###################################################
### code chunk number 17: sppa.Rnw:420-436
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
plot(spbdry, axes=TRUE, lwd=0.5)
plot(sproads, add=TRUE, lwd=2, col="darkslategrey")
c_c <- (spasthma$Asthma == "case") + 1
plot(spasthma[c_c == 1,], add=TRUE, pch=4, cex=0.6, col="mediumaquamarine")
plot(spasthma[c_c == 2,], add=TRUE, pch=17, cex=0.75, col="goldenrod2")
#plot(spsrc, pch=4, add=TRUE, cex=0.9, lwd=2, col="grey40")
#plot(spsrc, pch=1, add=TRUE, cex=0.9, lwd=2, col="grey40")
plot(spsrc, pch=22, add=TRUE, cex=1.2, bg="brown4")
#plot(spsrc, pch=4, add=TRUE, cex=1)
legend("bottomright", legend=c("controls", "cases", "pollution sources"), pch=c(4, 17, 22), pt.cex=c(0.6, 0.75, 1.2), pt.bg=c(NA, NA, "brown4"), col=c("mediumaquamarine", "goldenrod2", "black"), bty="n") 
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 18: sppa.Rnw:544-545
###################################################
tic <- Sys.time()


###################################################
### code chunk number 19: sppa.Rnw:547-557 (eval = FALSE)
###################################################
## set.seed(120109)
## r <- seq(0, sqrt(2)/6, by = 0.005)
## envjap <- envelope(as(spjpines1, "ppp"), fun=Gest, r=r, nrank=2, nsim=99)
## envred <- envelope(as(spred, "ppp"), fun=Gest, r=r, nrank=2, nsim=99)
## envcells <- envelope(as(spcells, "ppp"), fun=Gest, r=r, nrank=2, nsim=99)
## Gresults <- rbind(envjap, envred, envcells) 
## Gresults <- cbind(Gresults, 
##    y=rep(c("JAPANESE", "REDWOOD", "CELLS"), each=length(r)))
## # CHANGED DATASET TO y RSB
## #	save(Gresults, envjap, envred, envcells, file="sppaGestEnv.RData")


###################################################
### code chunk number 20: sppa.Rnw:559-560
###################################################
#cat("%", difftime(Sys.time(), tic, units="secs"), "seconds\n\n")


###################################################
### code chunk number 21: sppa.Rnw:562-563
###################################################
load("Gresults.RData")


###################################################
### code chunk number 22: sppa.Rnw:595-620
###################################################
.iwidth <- 5
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())

#Alternative plot using Trellis graphics. TOO COMPICATED?
#Gresults<-data.frame(DATASET=rep(c("Japanese", "Redwood", "Cells"), each=length(r)))
#Gresults<-cbind(Gresults, G=c(Gjap, Gred, Gcells))
#Gresults<-cbind(Gresults, fbar=c(envjap$fbar, envred$fbar, envcells$fbar) )
#Gresults<-cbind(Gresults, Lenv=c(envjap$L, envred$L, envcells$L))
#Gresults<-cbind(Gresults, Uenv=c(envjap$U, envred$U, envcells$U))

#library(lattice)
# CHANGED DATASET TO y, fill to col RSB
print(xyplot(obs~theo|y , data=Gresults, type="l", 
xlab = "theoretical", ylab = "observed", # EJP
panel=function(x, y, subscripts) {
   lpolygon(c(x, rev(x)), 
   c(Gresults$lo[subscripts], rev(Gresults$hi[subscripts])),
   border="gray", col="gray"
)
llines(x, y, col="black", lwd=2)
}))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 23: sppa.Rnw:658-659
###################################################
tic <- Sys.time()


###################################################
### code chunk number 24: sppa.Rnw:661-670 (eval = FALSE)
###################################################
## set.seed(30)
## Fenvjap<-envelope(as(spjpines1, "ppp"), fun=Fest, r=r, nrank=2, nsim=99)
## Fenvred<-envelope(as(spred, "ppp"), fun=Fest, r=r, nrank=2, nsim=99)
## Fenvcells<-envelope(as(spcells, "ppp"), fun=Fest, r=r, nrank=2, nsim=99)
## Fresults<-rbind(Fenvjap, Fenvred, Fenvcells)
## Fresults<-cbind(Fresults, 
##    y=rep(c("JAPANESE", "REDWOOD", "CELLS"), each=length(r)))
## # CHANGED DATASET TO y RSB
## #	save(Fresults, Fenvjap, Fenvred, Fenvcells, file="sppaFresults.RData")


###################################################
### code chunk number 25: sppa.Rnw:672-673
###################################################
#cat("%", difftime(Sys.time(), tic, units="secs"), "seconds\n\n")


###################################################
### code chunk number 26: sppa.Rnw:675-676
###################################################
load("Fresults.RData")


###################################################
### code chunk number 27: sppa.Rnw:700-716
###################################################
.iwidth <- 5
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
# CHANGED DATASET TO y, fill to col RSB
print(xyplot(obs~theo|y , data=Fresults, type="l", 
xlab = "theoretical", ylab = "observed", # EJP
panel=function(x, y, subscripts) {
   lpolygon(c(x, rev(x)), 
   c(Fresults$lo[subscripts], rev(Fresults$hi[subscripts])),
   border="gray", col="gray"
   )
   llines(x, y, col="black", lwd=2)
}))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 28: sppa.Rnw:925-932
###################################################
set.seed(7)
x<-runif(10)
nx<-length(x)
bw<-.1

k<-density(x, bw=bw, kernel="biweight")
k$y<-k$y*nx


###################################################
### code chunk number 29: sppa.Rnw:939-951
###################################################
.iwidth <- 6
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
plot(k, ylab="Intensity", main="")
points(x, rep(0, nx), pch=20)
for(i in 1:length(x))
  lines(density(x[i], bw=bw, kernel="biweight"), lty=2)

legend(x=14, y=0.6, legend=c("Intensity", "Kernel"), lty=c(1,2))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 30: sppa.Rnw:996-997
###################################################
library(splancs)


###################################################
### code chunk number 31: sppa.Rnw:999-1008
###################################################
mserwq<-mse2d(as.points(coordinates(spred)),
 as.points(list(x=c(0,1,1,0), y=c(0,0,1,1))), 100, .15)
bwq<-mserwq$h[which.min(mserwq$mse)]
bwq

#Spatstat code
mserw<-bw.diggle(as(spred, "ppp"))
bw<-as.numeric(mserw)
bw


###################################################
### code chunk number 32: sppa.Rnw:1029-1045
###################################################
.iwidth <- 5
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
par(mfrow=c(1,2))

#Quartic kernel
plot(mserwq$h, mserwq$mse, xlab="Bandwidth", ylab="MSE", type="l", ylim=c(-2,50), main="Quartic kernel")
i<-which.min(mserwq$mse)
points(mserwq$h[i], mserwq$mse[i])

#Gaussian kernel
plot(mserw, main="Gaussian kernel", xlab="Bandwidth", ylab="MSE")
points(attr(mserw, "h")[attr(mserw, "iopt")], bw)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 33: sppa.Rnw:1074-1086
###################################################
library(splancs)
poly <- as.points(list(x = c(0, 0, 1, 1), y = c(0, 1, 1, 0)))
sG <- Sobj_SpatialGrid(spred, maxDim=100)$SG
grd <- slot(sG, "grid")
summary(grd)
k0 <- spkernel2d(spred, poly, h0=bw, grd)
k1 <- spkernel2d(spred, poly, h0=.05, grd)
k2 <- spkernel2d(spred, poly, h0=.1, grd)
k3 <- spkernel2d(spred, poly, h0=.15, grd)
df <- data.frame(k0=k0, k1=k1, k2=k2, k3=k3) 
kernels <- SpatialGridDataFrame(grd, data=df)
summary(kernels)


###################################################
### code chunk number 34: sppa.Rnw:1118-1129
###################################################
cc <- coordinates(kernels)
xy<-list(x=cc[,1], y=cc[,2])
k4<-density(as(spred, "ppp"), .5*bw, dimyx=c(100, 100), xy=xy)
kernels$k4<-as(k4, "SpatialGridDataFrame")$v
k5<-density(as(spred, "ppp"), .5*.05, dimyx=c(100, 100), xy=xy)
kernels$k5<-as(k5, "SpatialGridDataFrame")$v
k6<-density(as(spred, "ppp"), .5*.1, dimyx=c(100, 100), xy=xy)
kernels$k6<-as(k6, "SpatialGridDataFrame")$v
k7<-density(as(spred, "ppp"), .5*.15, dimyx=c(100, 100), xy=xy)
kernels$k7<-as(k7, "SpatialGridDataFrame")$v
summary(kernels)


###################################################
### code chunk number 35: sppa.Rnw:1138-1159
###################################################
.iwidth <- 6
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
library(RColorBrewer)
gp <- brewer.pal(8, "Reds")
#print(spplot(kernels, at=seq(0,2000,length.out=22),
# col.regions=colorRampPalette(gp)(21), 
#names.attr=c(paste("bw=",bw, sep="", collapse=""),
# "bw=0.05", "bw=0.1","bw=0.15", paste("bw2=",.5*bw, sep="", collapse=""),
# "bw2=0.025", "bw2=0.05","bw2=0.075") ) )
print(spplot(kernels, at=c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000),
#seq(0,4500,length.out=22),
 col.regions=colorRampPalette(gp)(15)[1:12], 
names.attr=c(paste("Q bw=",round(bw, digits=4), sep="", collapse=""),
"Q bw=0.05", "Q bw=0.1","Q bw=0.15", paste("G bw=", round(.5*bw, digits=4),
 sep="", collapse=""), "G bw=0.025", "G bw=0.05","G bw=0.075"), cex=0.7, colorkey=FALSE))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 36: sppa.Rnw:1230-1249
###################################################
#Fit parametric model
loglambda<-function(x, alpha, beta)
{
    l<-alpha+sum(beta*c(x, x*x, prod(x)))
    return(l)
}

L<-function(alphabeta, x)
{
    l<-apply(x,1,loglambda, alpha=alphabeta[1], beta=alphabeta[-1])
    l<-sum(l)
    intL<-adaptIntegrate(lowerLimit=c(0,0), upperLimit=c(1,1), fDim=1,
        tol=1e-8, f=function(x, alpha=alphabeta[1], beta=alphabeta[-1])
        {
            exp(loglambda(x, alpha, beta))
    })
    l<-l-intL$integral
    return(l)#Optim minimises
}


###################################################
### code chunk number 37: sppa.Rnw:1266-1270
###################################################
library(cubature)
data(lansing)
x<-as.points(lansing[lansing$marks=="maple",])
#x<-as.points(lansing[lansing$species=="hickory",])


###################################################
### code chunk number 38: sppa.Rnw:1272-1275 (eval = FALSE)
###################################################
## #Maximise log-likelihood
## optbeta<-optim(par=c(log(514),0,0,0,0,0), fn=L, control=list(maxit=1000, fnscale=-1), x=x)
## #	save(optbeta, file="sppaOptbeta.RData")


###################################################
### code chunk number 39: sppa.Rnw:1277-1278
###################################################
load("sppaOptbeta.RData")


###################################################
### code chunk number 40: sppa.Rnw:1309-1329
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.pwd <- 0.65
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
grd <- GridTopology(cellcentre.offset=c(0.005,0.005), cellsize=c(0.01, 0.01),
  cells.dim=c(100, 100))
lambda<-exp(apply(coordinates(grd),1, function(X, alpha, beta)
    {
        loglambda(X, alpha, beta)
    }, alpha=optbeta$par[1], beta=optbeta$par[-1]
 ))

parint<-SpatialGridDataFrame(grd, data=data.frame(intensity=lambda))

lyt<-list("sp.points", SpatialPoints(x), pch=19, col="black", cex=0.7)
print(spplot(parint, at=seq(0,1400,length.out=8),
 col.regions=colorRampPalette(gp)(7), sp.layout=list(lyt)))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 41: sppa.Rnw:1346-1348
###################################################
lmaple<-lansing[lansing$marks=="maple",]
ppm(Q=lmaple, trend=~x+y+I(x^2)+I(y^2)+I(x*y))


###################################################
### code chunk number 42: sppa.Rnw:1417-1418
###################################################
tic <- Sys.time()


###################################################
### code chunk number 43: sppa.Rnw:1420-1430
###################################################
set.seed(30)
Kenvjap<-envelope(as(spjpines1, "ppp"), fun=Kest, r=r, nrank=2, nsim=99)
Kenvred<-envelope(as(spred, "ppp"), fun=Kest, r=r, nrank=2, nsim=99)
Kenvcells<-envelope(as(spcells, "ppp"), fun=Kest, r=r, nrank=2, nsim=99)
#Merge results in a data frame
Kresults<-rbind(Kenvjap, Kenvred, Kenvcells)
Kresults<-cbind(Kresults, 
   y=rep(c("JAPANESE", "REDWOOD", "CELLS"), each=length(r)))
# CHANGED DATASET TO y RSB
#	save(Kenvjap, Kenvred, Kenvcells, Kresults, file="sppaKresults.RData")


###################################################
### code chunk number 44: sppa.Rnw:1432-1433
###################################################
#cat("%", difftime(Sys.time(), tic, units="secs"), "seconds\n\n")


###################################################
### code chunk number 45: sppa.Rnw:1468-1485
###################################################
.iwidth <- 5
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
# CHANGED DATASET TO y, fill to col RSB
print(xyplot((obs-theo)~r|y , data=Kresults, type="l",
   ylim= c(-.06, .06), ylab=expression(hat(K) (r)  - pi * r^2),
   panel=function(x, y, subscripts) {
      Ktheo<- Kresults$theo[subscripts]
      lpolygon(c(r, rev(r)),
      c(Kresults$lo[subscripts]-Ktheo, rev(Kresults$hi[subscripts]-Ktheo)),
        border="gray", col="gray"
      )
      llines(r, Kresults$obs[subscripts]-Ktheo, lty=2, lwd=1.5, col="black")
}))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 46: sppa.Rnw:1698-1706
###################################################
bwasthma<-.06

pppasthma<-as(spasthma, "ppp")
pppasthma$window<-as(spbdry, "owin")

marks(pppasthma)<-relevel(pppasthma$marks$Asthma, "control")

#bwrr<-bw.relrisk(pppasthma, hmax=.5)


###################################################
### code chunk number 47: sppa.Rnw:1716-1721
###################################################
#Here we should include some stuff on the selection of the bandwidth
#bwasthma <- .125 #.275
#VGR: This should be the equivalent bandwith for the Gaussian kernel
#used in density.ppp. In ASDAR 1st ed. we used a quartic kernel for this.
bwasthma<-.06


###################################################
### code chunk number 48: sppa.Rnw:1744-1760 (eval = FALSE)
###################################################
## #OLD CODE based on splancs
## source("cv.R")#Maybe change it to read it from the web
## 
## h<-seq(0.20, .4, by=.025)
## cvres<-cvkernratio(rbind(coordinates(cases), coordinates(controls)), ncases, 
##    ncontrols, pbdry, h, gt )
## 
## plot(h, cvres[[2]])#hopt=.275
## #cvres<-list(c(0.0110665423247158, 0.00752722499415707, 0.00555282093087163, 
## #0.00415026712356268, 0.00308783702604925, 0.00233784029819392, 
## #0.00177026933743012, 0.00139550299449857, 0.00112439898687647
## #), c(-0.00293922649041061, -0.0169882573788171, -0.0231734237606048, 
## #-0.0246841045583454, -0.0240797773260347, -0.0228125113740377, 
## #-0.0214436687245140, -0.0202705522103254, -0.0193039039468653
## #))
## 


###################################################
### code chunk number 49: sppa.Rnw:1773-1778 (eval = FALSE)
###################################################
## #OLD CODE
## library(maptools)
## sG <- Sobj_SpatialGrid(spbdry, maxDim=50)$SG
## gt <- slot(sG, "grid")
## summary(gt)


###################################################
### code chunk number 50: sppa.Rnw:1798-1800 (eval = FALSE)
###################################################
## #OLD CODE
## pbdry <- slot(slot(slot(spbdry, "polygons")[[1]], "Polygons")[[1]], "coords")


###################################################
### code chunk number 51: sppa.Rnw:1821-1830 (eval = FALSE)
###################################################
## #OLD CODE
## library(splancs)
## cases<-spasthma[spasthma$Asthma=="case",]
## ncases<-nrow(cases)
## controls<-spasthma[spasthma$Asthma=="control",]
## ncontrols<-nrow(controls)
## 
## kcases<-spkernel2d(cases, pbdry, h0=bwasthma, gt)
## kcontrols<-spkernel2d(controls, pbdry, h0=bwasthma, gt)


###################################################
### code chunk number 52: sppa.Rnw:1838-1845
###################################################
cases<-unmark(subset(pppasthma, marks(pppasthma) =="case"))
ncases<-npoints(cases)
controls<-unmark(subset(pppasthma, marks(pppasthma) =="control"))
ncontrols<-npoints(controls)

kcases<-density(cases, bwasthma)
kcontrols<-density(controls, bwasthma)


###################################################
### code chunk number 53: sppa.Rnw:1870-1877 (eval = FALSE)
###################################################
## #OLD CODE
## df0 <- data.frame(kcases=kcases, kcontrols=kcontrols)
## spkratio0 <- SpatialGridDataFrame(gt, data=df0)
## spkratio <- as(spkratio0, "SpatialPixelsDataFrame")
## spkratio$kratio <- spkratio$kcases/spkratio$kcontrols
## is.na(spkratio$kratio) <- !is.finite(spkratio$kratio)
## spkratio$logratio <- log(spkratio$kratio)-log(ncases/ncontrols)


###################################################
### code chunk number 54: sppa.Rnw:1880-1888
###################################################

spkratio0<-as(kcases, "SpatialGridDataFrame")
names(spkratio0)<-"kcases"
spkratio0$kcontrols<-as(kcontrols, "SpatialGridDataFrame")$v
spkratio<-as(spkratio0, "SpatialPixelsDataFrame")

spkratio$kratio <- spkratio$kcases/spkratio$kcontrols
spkratio$logratio <- log(spkratio$kratio)-log(ncases/ncontrols)


###################################################
### code chunk number 55: sppa.Rnw:1955-1958 (eval = FALSE)
###################################################
## #OLD CODE
## idxinbdry <- over(sG, spbdry)$x
## idxna <- !is.na(idxinbdry)


###################################################
### code chunk number 56: sppa.Rnw:1977-1983
###################################################
niter <- 99
ratio <- rep(NA, niter)
#pvaluemap <- rep(0, sum(idxna))#OLD CODE
pvaluemap <- rep(0, nrow(spkratio))
#rlabelratio <- matrix(NA, nrow=niter, ncol=sum(idxna))#OLD CODE
rlabelratio <- matrix(NA, nrow=niter, ncol=nrow(spkratio))


###################################################
### code chunk number 57: sppa.Rnw:1997-2011 (eval = FALSE)
###################################################
## #OLD CODE
## set.seed(1)
## for(i in 1:niter)
## {
## idxrel <- sample(spasthma$Asthma) == "case"
## casesrel <- spasthma[idxrel,]
## controlsrel <- spasthma[!idxrel,]
## kcasesrel <- spkernel2d(casesrel, pbdry, h0=bwasthma, gt)
## kcontrolsrel <- spkernel2d(controlsrel, pbdry, h0=bwasthma, gt)
## kratiorel <- kcasesrel[idxna]/kcontrolsrel[idxna]
##    is.na(kratiorel) <- !is.finite(kratiorel)
## rlabelratio[i,] <- kratiorel
## pvaluemap <- pvaluemap + (spkratio$kratio < kratiorel)
## }


###################################################
### code chunk number 58: sppa.Rnw:2014-2015
###################################################
tic <- Sys.time()


###################################################
### code chunk number 59: sppa.Rnw:2017-2030 (eval = FALSE)
###################################################
## set.seed(1)
## for(i in 1:niter)
## {
## pppasthma0<-rlabel(pppasthma)
## casesrel <- unmark(subset(pppasthma0, marks(pppasthma0) =="case"))
## controlsrel <- unmark(subset(pppasthma0, marks(pppasthma0) =="control"))
## 
## kcasesrel <- density(casesrel, bwasthma)
## kcontrolsrel <- density(controlsrel, bwasthma)
## kratiorel <- eval.im(kcasesrel/kcontrolsrel)
## rlabelratio[i,] <- as(as(kratiorel, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")$v
## pvaluemap <- pvaluemap + (spkratio$kratio < rlabelratio[i,])
## }


###################################################
### code chunk number 60: sppa.Rnw:2032-2033 (eval = FALSE)
###################################################
## save(rlabelratio, pvaluemap, file="loop_rlabelratio.RData")


###################################################
### code chunk number 61: sppa.Rnw:2035-2036
###################################################
#cat("%", difftime(Sys.time(), tic, units="secs"), "seconds\n\n")


###################################################
### code chunk number 62: sppa.Rnw:2063-2065
###################################################
#pre-load RSB
load("loop_rlabelratio.RData")


###################################################
### code chunk number 63: sppa.Rnw:2069-2078 (eval = FALSE)
###################################################
## #OLD CODE
## idxna2 <- apply(rlabelratio, 2, function(x) all(is.finite(x)))
## rhomean <- apply(rlabelratio[, idxna2], 2, mean)
## c <- prod(slot(gt, "cellsize"))
## ratiorho <- c*sum((spkratio$kratio[idxna2]-ncases/ncontrols)^2)
## ratio <- c*apply(rlabelratio[,idxna2], 1, 
##  function(X, rho0 ){sum((X-rho0)^2)}, rho0=ncases/ncontrols
## )
## pvaluerho <- (sum(ratio > ratiorho)+1)/(niter+1)


###################################################
### code chunk number 64: sppa.Rnw:2081-2087
###################################################
cellsize<-kcontrols$xstep*kcontrols$ystep
ratiorho <- cellsize*sum((spkratio$kratio-ncases/ncontrols)^2)
ratio <- cellsize*apply(rlabelratio, 1, 
 function(X, rho0 ){sum((X-rho0)^2)}, rho0=ncases/ncontrols
)
pvaluerho <- (sum(ratio > ratiorho)+1)/(niter+1)


###################################################
### code chunk number 65: sppa.Rnw:2125-2129
###################################################
spkratio$pvaluemap <- (pvaluemap+1)/(niter+1)
imgpvalue <- as.image.SpatialGridDataFrame(spkratio["pvaluemap"])
clpvalue <- contourLines(imgpvalue, levels=c(0,.05, .95, 1))
cl <- ContourLines2SLDF(clpvalue)


###################################################
### code chunk number 66: sppa.Rnw:2141-2170
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(RColorBrewer)
cl05 <- cl[cl$level == "0.05",]
xzx <- slot(slot(cl05, "lines")[[1]], "Lines")
cl05a <- SpatialLines(list(Lines(xzx, ID="0.05")))
lyt05 <- list("sp.lines", cl05a, lwd=2, lty=2, col="grey95")
lyt95 <- list("sp.lines", cl[cl$level == "0.95",], lwd=2, lty=1)
lytb <- list("sp.polygons", spbdry)
lytp <- list("sp.points", spsrc, cex=0.9, pch=4, col="grey95", lwd=3)
brks <- quantile(spkratio$kratio[spkratio$kratio>0], seq(0,1,1/10), na.rm=TRUE)
brks[1] <- 0
lbrks <- formatC(brks, 3, 6, "g", " ")
cols <- colorRampPalette(brewer.pal(7, "Reds"))(length(brks)-1)
# RSB quietening greys
colorkey<-list(labels=lbrks,
  at=(0:10)/10, height=.5)

print(spplot(spkratio, "kratio",
   col.regions=cols,
   do.log=TRUE, 
   colorkey=colorkey,
   at=c(0, brks[-c(1,11)], max(spkratio$kratio, na.rm=TRUE)),
   sp.layout=list(lyt05, lyt95, lytb, lytp) 
))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 67: sppa.Rnw:2233-2234
###################################################
rrbw<-bw.relrisk(pppasthma, hmax=.5)


###################################################
### code chunk number 68: sppa.Rnw:2255-2269 (eval = FALSE)
###################################################
## #OLD CODE
## source("cv.R")#Maybe change it to read it from the web
## 
## h<-seq(0.20, .4, by=.025)
## cvresbin<-cvbinreg(rbind(coordinates(cases), coordinates(controls)), ncases,
##    ncontrols, pbdry, h)
## 
## plot(h, cvresbin, type="l")#hopt=.225
## 
## #cvresbin<-c(1.56768244603478, 1.56761940865349, 1.56772491628267, 
## #1.56791042030484, 
## #1.56811754075041, 1.56831703621994, 1.56849753781631, 1.56865650893262, 
## #1.56879510590190)
## 


###################################################
### code chunk number 69: sppa.Rnw:2271-2274
###################################################
#bwasthmap <- 0.035
#bwasthmap <- 0.125 #0.225#OLD CODE
bwasthmap <- 0.06 #0.225


###################################################
### code chunk number 70: sppa.Rnw:2280-2290 (eval = FALSE)
###################################################
## #OLD CODE
## #Estimator of the probability
## lambda1 <- spkernel2d(cases, pbdry, h0=bwasthmap, gt)
## lambda0 <- spkernel2d(controls, pbdry, h0=bwasthmap, gt)
## 
## lambda1 <- lambda1[idxna]
## lambda0 <- lambda0[idxna]
## 
## spkratio$prob <- lambda1/(lambda1+lambda0)
## is.na(spkratio$prob) <- !is.finite(spkratio$prob)


###################################################
### code chunk number 71: sppa.Rnw:2293-2295
###################################################
rr<-relrisk(pppasthma, bwasthmap)
spkratio$prob<-as(as(rr, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")$v


###################################################
### code chunk number 72: sppa.Rnw:2307-2317
###################################################
.iwidth <- 5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
ats <- seq(0,max(spkratio$prob),length.out=11)
cols <- colorRampPalette(brewer.pal(8, "Reds"))(length(ats)-1)
# RSB quietening greys
print(spplot(spkratio, "prob", col.regions=cols, at=ats, sp.layout=list(lytb, lytp)))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 73: sppa.Rnw:2375-2386
###################################################
spasthma$y <- as.integer(!as.integer(spasthma$Asthma)-1)
ccasthma <- coordinates(spasthma)
spasthma$x1 <- ccasthma[,1]
spasthma$x2 <- ccasthma[,2]
spasthma$dist1 <- sqrt(spasthma$d2source1)
spasthma$dist2 <- sqrt(spasthma$d2source2)
spasthma$dist3 <- sqrt(spasthma$d2source3)
spasthma$droads <- sqrt(spasthma$roaddist2)
spasthma$smoking <- as.factor(as.numeric(spasthma$Nsmokers>0))
spasthma$Genderf<- as.factor(spasthma$Gender)
spasthma$HayFeverf<- as.factor(spasthma$HayFever)


###################################################
### code chunk number 74: sppa.Rnw:2388-2390
###################################################
library(mgcv)
gasthma<-gam(y~1+dist1+dist2+dist3+droads+Genderf+Age+HayFeverf+smoking+s(x1,x2), data=spasthma[spasthma$Gender==1 | spasthma$Gender==2, ], family=binomial)


###################################################
### code chunk number 75: sppa.Rnw:2392-2393
###################################################
summary(gasthma)


###################################################
### code chunk number 76: sppa.Rnw:2395-2398
###################################################
sumGasth <- summary(gasthma)
cpv <- sumGasth$p.pv
spv <- sumGasth$s.table[4]


###################################################
### code chunk number 77: sppa.Rnw:2483-2486
###################################################
D2_mat <- as.matrix(spasthma$dist2)
RHO <- ncases/ncontrols
expsource2<-tribble(ccflag=spasthma$y, vars=D2_mat, rho=RHO, alphas=1, betas=1)


###################################################
### code chunk number 78: sppa.Rnw:2488-2489
###################################################
print(expsource2)


###################################################
### code chunk number 79: sppa.Rnw:2491-2494
###################################################
#Hay fever
Hay_mat <- as.matrix(spasthma$HayFever)
exphay <- tribble(ccflag=spasthma$y, rho=RHO, covars=Hay_mat, thetas=1)


###################################################
### code chunk number 80: sppa.Rnw:2496-2497
###################################################
print(exphay)


###################################################
### code chunk number 81: sppa.Rnw:2524-2526
###################################################
expsource2hay<-tribble(ccflag=spasthma$y, vars=D2_mat, rho=RHO,
 alphas=1, betas=1, covars=Hay_mat, thetas=1)


###################################################
### code chunk number 82: sppa.Rnw:2626-2634
###################################################
Kdif<-function(Xppp, r, cr="border")
{
	k1<-Kest(Xppp[marks(Xppp)=="case"], r=r, correction=cr)
	k2<-Kest(Xppp[marks(Xppp)=="control"], r=r, correction=cr)

	res<-data.frame(r=r, D=k1[[cr]]-k2[[cr]])
	return(fv(res, valu="D", fname="D"))
}


###################################################
### code chunk number 83: sppa.Rnw:2642-2648 (eval = FALSE)
###################################################
## r<-seq(0, .15, by=.01)
## 
## envKdif<-envelope(pppasthma, Kdif, r=r, nsim=99, cr="iso",  nrank=2,
## savefuns=TRUE,
##    simulate=expression(rlabel(pppasthma)))
## #plot(envKdif)


###################################################
### code chunk number 84: sppa.Rnw:2651-2652
###################################################
save(file="envKdif.RData", list=c("r", "envKdif"))


###################################################
### code chunk number 85: sppa.Rnw:2659-2660
###################################################
load("envKdif.RData")


###################################################
### code chunk number 86: sppa.Rnw:2663-2669 (eval = FALSE)
###################################################
## #OLD CODE
## s<-seq(0, .15, by=.01)
## khcases<-khat(coordinates(cases), pbdry, s)
## khcontrols<-khat(coordinates(controls), pbdry, s)
## khcov<-khvmat(coordinates(cases), coordinates(controls), pbdry, s)
## T0<-sum( ((khcases-khcontrols))/sqrt(diag(khcov)))


###################################################
### code chunk number 87: sppa.Rnw:2672-2674
###################################################
khcases<-Kest(cases, r=r, correction="isotropic")
khcontrols<-Kest(controls, r=r, correction="isotropic")


###################################################
### code chunk number 88: sppa.Rnw:2687-2690
###################################################
niter<-99
T<-rep(NA, niter)
set.seed(1234)


###################################################
### code chunk number 89: sppa.Rnw:2692-2706 (eval = FALSE)
###################################################
## #OLD CODE
## khcasesrel<-matrix(NA, nrow=length(s), ncol=niter)
## khcontrolsrel<-matrix(NA, nrow=length(s), ncol=niter)
## for(i in 1:niter)
## {
##     idxrel<-sample(spasthma$Asthma)=="case"
##     casesrel<-coordinates(spasthma[idxrel,])
##     controlsrel<-coordinates(spasthma[!idxrel,])
##     khcasesrel[,i]<-khat(casesrel, pbdry, s)
##     khcontrolsrel[,i]<-khat(controlsrel, pbdry, s)
##     khdiff <- khcasesrel[,i]-khcontrolsrel[,i]
##     T[i]<-sum(khdiff/sqrt(diag(khcov)))
## }
## #	save(T, khcasesrel, khcontrolsrel, file="sppaTch58.RData")


###################################################
### code chunk number 90: sppa.Rnw:2709-2710 (eval = FALSE)
###################################################
## save(khcasesrel, khcontrolsrel, file="sppaTch84-2ed.RData")


###################################################
### code chunk number 91: sppa.Rnw:2712-2713
###################################################
load("sppaTch84-2ed.RData")


###################################################
### code chunk number 92: sppa.Rnw:2716-2728 (eval = FALSE)
###################################################
## khcasesrel<-matrix(NA, nrow=length(r), ncol=niter)
## khcontrolsrel<-matrix(NA, nrow=length(r), ncol=niter)
## 
## for(i in 1:niter)
## {
## 	pppasthma0<-rlabel(pppasthma)
## 	casesrel <- unmark(subset(pppasthma0, marks(pppasthma0) =="case"))
## 	controlsrel <- unmark(subset(pppasthma0, marks(pppasthma0) =="control"))
## 
## 	khcasesrel[,i]<-Kest(casesrel, r=r, correction="isotropic")$iso
## 	khcontrolsrel[,i]<-Kest(controlsrel, r=r, correction="isotropic")$iso
## }


###################################################
### code chunk number 93: sppa.Rnw:2742-2751
###################################################
#Compute diagonal of var-cov matrix of dif. K-functions
#khcovdiag<-apply(khcasesrel-khcontrolsrel, 1,var)
simfuns<-as.data.frame(attr(envKdif, "simfuns"))[,-1]
khcovdiag<-apply(simfuns, 1, var)
#Test statistics
T0<-sum( ((khcases$iso-khcontrols$iso)/sqrt(khcovdiag))[-1])
T<-apply(simfuns, 2, function(X){
	sum((X/sqrt(khcovdiag))[-1])
})


###################################################
### code chunk number 94: sppa.Rnw:2753-2756 (eval = FALSE)
###################################################
## #OLD CODE
## #load("sppaTch58.RData")
## load("sppaTch58-2ndEd.RData")


###################################################
### code chunk number 95: sppa.Rnw:2758-2759
###################################################
pvalue<-(sum(T>T0)+1)/(niter+1)


###################################################
### code chunk number 96: sppa.Rnw:2772-2776
###################################################
.iwidth <- 10#5
.iheight <- 5#3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())


###################################################
### code chunk number 97: sppa.Rnw:2778-2799
###################################################
#OLD CODE
#plot(s, khcases-khcontrols, type="l", 
#  ylab="D(s)", ylim=c(-.015, .015))#ylim=c(-11.5, 11.5))
#lines(s, -1.96*sqrt(diag(khcov)), lty=2)
#lines(s, +1.96*sqrt(diag(khcov)), lty=2)

#NEW CODE
#plot(r, khcases$iso-khcontrols$iso, type="l", 
#  ylab="D(s)", ylim=c(-.015, .015))#ylim=c(-11.5, 11.5))
plot(envKdif)
lines(r, -1.96*sqrt(khcovdiag), lty=2)
lines(r, +1.96*sqrt(khcovdiag), lty=2)

#Compute envelopes
#envel<-apply(khcasesrel-khcontrolsrel, 1, function(X){quantile(X, c(.025, .975))})
#lines(r, envel[1,], lty=3, lwd=2)
#lines(r, envel[2,], lty=3, lwd=2)
#
#legend("bottomleft", 
#   legend=c("Actual value", "Approx. 95% C.I.", "Sim. 95% envelopes"),
#   lty=1:3, lwd=c(1,1,2), bty="n")


###################################################
### code chunk number 98: sppa.Rnw:2801-2803
###################################################
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 99: sppa.Rnw:2920-2934 (eval = FALSE)
###################################################
## #OLD CODE
## glmasthma<-glm(y~HayFeverf, data=spasthma, family="binomial")
## prob<-fitted(glmasthma)
## weights<-exp(glmasthma$linear.predictors)
## #weights<-fitted(glmasthma)
## library(spatialkernel)
## setkernel("gaussian")
## lambda0<- lambdahat(coordinates(controls), bwasthma, coordinates(cases), 
##    pbdry, FALSE)$lambda
## lambda1<- weights[spasthma$Asthma=="case"]*lambda0
## 
## ratiocc<-ncases/ncontrols
## kihnocov<-kinhat(coordinates(cases), ratiocc*lambda0, pbdry,s)$k
## kih<-kinhat(coordinates(cases), lambda1, pbdry,s)$k


###################################################
### code chunk number 100: sppa.Rnw:2937-2938
###################################################
tic <- Sys.time()


###################################################
### code chunk number 101: sppa.Rnw:2940-2950
###################################################
glmasthma<-glm(y~HayFeverf, data=spasthma, family="binomial")
prob<-fitted(glmasthma)
weights<-exp(glmasthma$linear.predictors)
#weights<-fitted(glmasthma)
lambda0<-  interp.im (kcontrols, coords(cases)[,1], coords(cases)[,2])
lambda1<- weights[marks(pppasthma) =="case"]*lambda0

ratiocc<-ncases/ncontrols
kihnocov<-Kinhom(cases, ratiocc*lambda0, r=r)
kih<-Kinhom(cases, lambda1, r=r)


###################################################
### code chunk number 102: sppa.Rnw:2952-2953
###################################################
#cat("%", difftime(Sys.time(), tic, units="secs"), "seconds\n\n")


###################################################
### code chunk number 103: sppa.Rnw:3002-3025
###################################################
#Relabel with weights
rlabelp<-function(Xppp, ncases, prob)
{
	idxsel<-sample(1:npoints(Xppp), ncases, prob=prob)
	marks(Xppp)<-"control"
	marks(Xppp)[idxsel]<-"case"
	return(Xppp)
}

#Compute K_{I,\lambda})
KIlambda<-function(Xppp, r, cr="iso", weights, sigma)
{
    idxrel<-marks(Xppp)=="case"
    casesrel<-unmark(Xppp[idxrel])
    controlsrel<-unmark(Xppp[!idxrel])
    lambda0rel<-interp.im(density(controlsrel, sigma), coords(casesrel)[,1],
        coords(casesrel)[,2])
    lambda1rel<-weights[idxrel]*lambda0rel
    KI<-Kinhom(casesrel, lambda1rel, r=r, correction=cr)
    res<-data.frame(r=r, KI=KI[[cr]])

   return(fv(res, valu="KI", fname="K_[I,lambda]"))
}


###################################################
### code chunk number 104: sppa.Rnw:3040-3049 (eval = FALSE)
###################################################
## set.seed(4567)
## envKInocov<-envelope(pppasthma,KIlambda, r=r, cr="iso", weights=weights, 
##    sigma=bwasthma, nsim=99, nrank=2, savefuns=TRUE, 
##    simulate=expression(rlabelp(pppasthma, ncases=ncases, prob=rep(ratiocc, npoints(pppasthma)))) )
## 
## envKIcov<-envelope(pppasthma,KIlambda, r=r, cr="iso", weights=weights, 
##    sigma=bwasthma, nsim=99, nrank=2, savefuns=TRUE, 
##    simulate=expression(rlabelp(pppasthma, ncases=ncases, prob=prob)) )
## 


###################################################
### code chunk number 105: sppa.Rnw:3051-3062
###################################################
save(file="KIenv.RData", list=c("envKInocov", "envKIcov"))
#plotenv<-function(kk){
#plot(kk$r, kk$obs-kk$mmean, type="l", ylim=c(-.015, .015))
#lines(kk$r, kk$lo-kk$mmean, col="gray")
#lines(kk$r, kk$hi-kk$mmean, col="gray")
#}
#
#par(mfrow=c(1,2))
#plotenv(envKInocov)
#plotenv(envKIcov)



###################################################
### code chunk number 106: sppa.Rnw:3068-3069
###################################################
load("KIenv.RData")


###################################################
### code chunk number 107: sppa.Rnw:3072-3087 (eval = FALSE)
###################################################
## #OLD CODE
## set.seed(1234)
## niter<-99
## kinhomrelnocov<-matrix(NA, nrow=length(s), ncol=niter)
## kinhomrel<-matrix(NA, nrow=length(s), ncol=niter)
## for(i in 1:niter)
## {
##     idxrel<-sample(spasthma$Asthma, prob=prob)=="case"
##     casesrel<-coordinates(spasthma[idxrel,])
##     controlsrel<-coordinates(spasthma[!idxrel,])
##     lambda0rel<-lambdahat(controlsrel, bwasthma, casesrel, pbdry, FALSE)$lambda
##     lambda1rel<-weights[idxrel]*lambda0rel
##     kinhomrelnocov[,i]<-kinhat(casesrel, ratiocc*lambda0rel, pbdry,s)$k
##     kinhomrel[,i]<-kinhat(casesrel, lambda1rel, pbdry,s)$k
## }


###################################################
### code chunk number 108: sppa.Rnw:3091-3092
###################################################
tic <- Sys.time()


###################################################
### code chunk number 109: sppa.Rnw:3094-3110 (eval = FALSE)
###################################################
## set.seed(1234)
## niter<-99
## kinhomrelnocov<-matrix(NA, nrow=length(r), ncol=niter)
## kinhomrel<-matrix(NA, nrow=length(r), ncol=niter)
## 
## for(i in 1:niter)
## {
##     idxrel<-sample(marks(pppasthma), prob=prob)=="case"
##     casesrel<-unmark(pppasthma[idxrel])
##     controlsrel<-unmark(pppasthma[!idxrel])
##     lambda0rel<-interp.im(density(controlsrel, bwasthma), coords(casesrel)[,1],
## 	coords(casesrel)[,2])
##     lambda1rel<-weights[idxrel]*lambda0rel
##     kinhomrelnocov[,i]<-Kinhom(casesrel, ratiocc*lambda0rel, r=r)$iso
##     kinhomrel[,i]<-Kinhom(casesrel, lambda1rel, r=r)$iso
## }


###################################################
### code chunk number 110: sppa.Rnw:3118-3120
###################################################
kinhomrelnocov<-as.data.frame(attr(envKInocov, "simfuns"))[,-1]
kinhomrel<-as.data.frame(attr(envKIcov, "simfuns"))[,-1]


###################################################
### code chunk number 111: sppa.Rnw:3123-3125 (eval = FALSE)
###################################################
## #save(kinhomrelnocov, kinhomrel, file="kinhom.RData")
## save(kinhomrelnocov, kinhomrel, file="kinhomenv.RData")


###################################################
### code chunk number 112: sppa.Rnw:3127-3129
###################################################
#load("kinhom.RData")
load("kinhomenv.RData")


###################################################
### code chunk number 113: sppa.Rnw:3131-3138 (eval = FALSE)
###################################################
## #NOT USING envelope()
## kinhsdnocov<-apply(kinhomrelnocov, 1, sd)
## kihmeannocov<-apply(kinhomrelnocov, 1,mean)
## D0nocov<-sum(((kihnocov$iso-kihmeannocov)/kinhsdnocov)[-1])
## Dnocov<-apply(kinhomrelnocov, 2, 
##    function(X){ sum(((X-kihmeannocov)/kinhsdnocov)[-1])})
## pvaluenocov<-(sum(Dnocov>D0nocov)+1)/(niter+1)


###################################################
### code chunk number 114: sppa.Rnw:3140-3146
###################################################
#Using envelope()
kinhsdnocov<-apply(kinhomrelnocov, 1, sd)
D0nocov<-sum(((envKInocov$obs-envKInocov$mmean)/kinhsdnocov)[-1])
Dnocov<-apply(kinhomrelnocov, 2, 
   function(X){ sum(((X-envKInocov$mmean)/kinhsdnocov)[-1])})
pvaluenocov<-(sum(Dnocov>D0nocov)+1)/(niter+1)


###################################################
### code chunk number 115: sppa.Rnw:3148-3155 (eval = FALSE)
###################################################
## #NOT USING envelope()
## kinhsd<-apply(kinhomrel, 1, sd)
## kihmean<-apply(kinhomrel, 1,mean)
## D0<-sum( ((kih$iso-kihmean)/kinhsd)[-1])
## D<-apply(kinhomrel, 2, 
##    function(X){ sum(((X-kihmean)/kinhsd)[-1])})
## pvalue<-(sum(D>D0)+1)/(niter+1)


###################################################
### code chunk number 116: sppa.Rnw:3157-3162
###################################################
kinhsd<-apply(kinhomrel, 1, sd)
D0<-sum( ((envKIcov$obs-envKIcov$mmean)/kinhsd)[-1])
D<-apply(kinhomrel, 2, 
   function(X){ sum(((X-envKIcov$mmean)/kinhsd)[-1])})
pvalue<-(sum(D>D0)+1)/(niter+1)


###################################################
### code chunk number 117: sppa.Rnw:3167-3168
###################################################
#cat("%", difftime(Sys.time(), tic, units="secs"), "seconds\n\n")


###################################################
### code chunk number 118: sppa.Rnw:3196-3239
###################################################
.iwidth <- 5
.iheight <-3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-sppa-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
oopar <- par(mfrow=c(1,2))
#OLD CODE
#plot(s, kihnocov-kihmeannocov, type="l", 
#   ylim= c(-0.06,  0.22),
#   xlab="s", ylab= expression(hat(k)[I][","][hat(lambda)]-"E[s]"),
#    main ="No covariates" )

#plot(r, kihnocov$iso-kihmeannocov, type="l", 
#   ylim= c(-0.06,  0.22),
#   xlab="s", ylab= expression(hat(k)[I][","][hat(lambda)]-"E[s]"),
#    main ="No covariates" )
#
#envnocov<-apply(kinhomrelnocov, 1, function(X){quantile(X, c(.025, .975))})
#lines(r, envnocov[1,]-kihmeannocov, lty=2)
#lines(r, envnocov[2,]-kihmeannocov, lty=2)
plot(r, envKInocov$obs-envKInocov$mmean, type="l", 
   ylim= c(-0.06,  0.06),
   xlab="s", ylab= expression(hat(k)[I][","][hat(lambda)]-"E[s]"),
   main ="No covariates" )
lines(r, envKInocov$lo-envKInocov$mmean, lty=2)
lines(r, envKInocov$hi-envKInocov$mmean, lty=2)

#plot(r, kih$iso-kihmean, type="l", ylim=c(-0.06,  0.22), #c(-2e-4, 2e-4),
#   xlab="s", ylab= expression(hat(k)[I][","][hat(lambda)]-"E[s]"),
#   main ="Adjusting for Hay Fever"  )
#
#env<-apply(kinhomrel, 1, function(X){quantile(X, c(.025, .975))})
#lines(r, env[1,]-kihmean, lty=2)
#lines(r, env[2,]-kihmean, lty=2)
plot(r, envKIcov$obs-envKIcov$mmean, type="l", 
   ylim= c(-0.06,  0.06),
   xlab="s", ylab= expression(hat(k)[I][","][hat(lambda)]-"E[s]"),
   main ="Adjusting for Hay Fever" )
lines(r, envKIcov$lo-envKIcov$mmean, lty=2)
lines(r, envKIcov$hi-envKIcov$mmean, lty=2)

par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 119: sppa.Rnw:3301-3309
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
### code chunk number 120: sppa.Rnw:3312-3313
###################################################
options(op)


