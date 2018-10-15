###################################################
# sppa_mod.R
# packages: spatstat, sp, maptools, lattice, rgdal, splancs, cubature, mgcv,
#           spatialkernel
# datasets: north_derby_asthma.zip
# provided: cv.R


###################################################
### chunk number 1: 
###################################################
rm(list=ls())
if ((site <- Sys.getenv("ASDAR_DOWNLOAD")) != "") {
  download.file(paste(site, "north_derby_asthma.zip", sep="/"),
    "north_derby_asthma.zip", mode="wb")
}

###################################################
### chunk number 6: 
###################################################
library(spatstat)
data(japanesepines)


###################################################
### chunk number 7: 
###################################################
summary(japanesepines)


###################################################
### chunk number 8: 
###################################################
library(maptools)


###################################################
### chunk number 9: 
###################################################
spjpines <- as(japanesepines, "SpatialPoints")
summary(spjpines)


###################################################
### chunk number 10: 
###################################################
spjpines1 <- elide(spjpines, scale=TRUE, unitsq=TRUE)
summary(spjpines1)


###################################################
### chunk number 11: 
###################################################
pppjap <- as(spjpines1, "ppp")
summary(pppjap)


###################################################
### chunk number 12: 
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
### chunk number 13: 
###################################################
library(lattice)
print(xyplot(y~x|DATASET, data=dpp, pch=19, aspect=1))


###################################################
### chunk number 14: 
###################################################
#fname <- zip.file.extract(file="spasthma.dbf",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./spasthma.dbf", overwrite=TRUE)
#fname <- zip.file.extract(file="spasthma.shp",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./spasthma.shp", overwrite=TRUE)
#fname <- zip.file.extract(file="spasthma.shx",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./spasthma.shx", overwrite=TRUE)

#fname <- zip.file.extract(file="spbdry.dbf",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./spbdry.dbf", overwrite=TRUE)
#fname <- zip.file.extract(file="spbdry.shp",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./spbdry.shp", overwrite=TRUE)
#fname <- zip.file.extract(file="spbdry.shx",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./spbdry.shx", overwrite=TRUE)

#fname <- zip.file.extract(file="sproads.dbf",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./sproads.dbf", overwrite=TRUE)
#fname <- zip.file.extract(file="sproads.shp",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./sproads.shp", overwrite=TRUE)
#fname <- zip.file.extract(file="sproads.shx",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./sproads.shx", overwrite=TRUE)

#fname <- zip.file.extract(file="spsrc.dbf",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./spsrc.dbf", overwrite=TRUE)
#fname <- zip.file.extract(file="spsrc.shp",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./spsrc.shp", overwrite=TRUE)
#fname <- zip.file.extract(file="spsrc.shx",
# zipname = "north_derby_asthma.zip")
#file.copy(fname, "./spsrc.shx", overwrite=TRUE)

unzip(zipfile="north_derby_asthma.zip")

###################################################
### chunk number 15: 
###################################################
library(rgdal)
spasthma <- readOGR(".", "spasthma")
spbdry <- readOGR(".", "spbdry")
spsrc <- readOGR(".", "spsrc")
sproads <- readOGR(".", "sproads")



###################################################
### chunk number 16: 
###################################################
plot(spbdry, axes=TRUE)
plot(sproads, add=TRUE, lty=2)
plot(spasthma, add=TRUE, pch=c(4,17)[(spasthma$Asthma == "case") + 1], cex=c(0.6, 0.75)[(spasthma$Asthma == "case") + 1])
plot(spsrc, pch=22, add=TRUE, cex=1.2, bg="grey60")


###################################################
### chunk number 17: 
###################################################
set.seed(30)


###################################################
### chunk number 19:  
###################################################
# spatstat 1.45-0 enforce a tighter spacing of the r vector 2016-03-22
# r <- seq(0, sqrt(2)/6, by = 0.005)
r <- seq(0, sqrt(2)/6, by = 0.001)
envjap <- envelope(as(spjpines1, "ppp"), fun=Gest, r=r, nrank=2, nsim=99)
envred <- envelope(as(spred, "ppp"), fun=Gest, r=r, nrank=2, nsim=99)
envcells <- envelope(as(spcells, "ppp"), fun=Gest, r=r, nrank=2, nsim=99)
Gresults <- rbind(envjap, envred, envcells) 
Gresults <- cbind(Gresults, 
   y=rep(c("JAPANESE", "REDWOOD", "CELLS"), each=length(r)))


###################################################
### chunk number 20: 
###################################################
print(xyplot(obs~theo|y, data=Gresults, type="l", 
	panel=function(x, y, subscripts)
	{
		lpolygon(c(x, rev(x)), 
		   c(Gresults$lo[subscripts], rev(Gresults$hi[subscripts])),
		   border="gray", col="gray"
		)

		llines(x, y, col="black", lwd=2)
	}
))


###################################################
### chunk number 21: 
###################################################
set.seed(30)


###################################################
### chunk number 23:  
###################################################
Fenvjap<-envelope(as(spjpines1, "ppp"), fun=Fest, r=r, nrank=2, nsim=99)
Fenvred<-envelope(as(spred, "ppp"), fun=Fest, r=r, nrank=2, nsim=99)
Fenvcells<-envelope(as(spcells, "ppp"), fun=Fest, r=r, nrank=2, nsim=99)
Fresults<-rbind(Fenvjap, Fenvred, Fenvcells)
Fresults<-cbind(Fresults, 
   y=rep(c("JAPANESE", "REDWOOD", "CELLS"), each=length(r)))


###################################################
### chunk number 24: 
###################################################
print(xyplot(obs~theo|y, data=Fresults, type="l", 
	panel=function(x, y, subscripts)
	{
		lpolygon(c(x, rev(x)), 
		   c(Fresults$lo[subscripts], rev(Fresults$hi[subscripts])),
		   border="gray", col="gray"
		)

		llines(x, y, col="black", lwd=2)
	}
))


###################################################
### chunk number 25: 
###################################################
set.seed(7)
x<-runif(10)
nx<-length(x)
bw<-.1

k<-density(x, bw=bw, kernel="biweight")
k$y<-k$y*nx


###################################################
### chunk number 26: 
###################################################
plot(k, ylab="Intensity", main="")
points(x, rep(0, nx), pch=20)
for(i in 1:length(x))
	lines(density(x[i], bw=bw, kernel="biweight"), lty=2)

legend(x=14, y=0.6, legend=c("Intensity", "Kernel"), lty=c(1,2))


###################################################
### chunk number 27: 
###################################################
library(splancs)
mserw<-mse2d(as.points(coordinates(spred)),
 as.points(list(x=c(0,1,1,0), y=c(0,0,1,1))), 100, .15)
bw<-mserw$h[which.min(mserw$mse)]


###################################################
### chunk number 28: 
###################################################
plot(mserw$h, mserw$mse, xlab="Bandwidth", ylab="MSE", type="l", ylim=c(-2,50))
i<-which.min(mserw$mse)
points(mserw$h[i], mserw$mse[i])


###################################################
### chunk number 29: 
###################################################
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
### chunk number 30: 
###################################################
xy<-list(x=coordinates(kernels)[,1], y=coordinates(kernels)[,2])

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
### chunk number 31: 
###################################################
gp <- grey.colors(5, 0.9, 0.45, 2.2)
print(spplot(kernels, at=seq(0,2000,length.out=22),
 col.regions=colorRampPalette(gp)(21), 
names.attr=c(paste("Quartic bw=",bw, sep="", collapse=""),
"Quartic bw=0.05", "Quartic bw=0.1","Quartic bw=0.15", 
paste("Gaussian bw=",.5*bw, sep="", collapse=""),
 "Gaussian bw=0.025", "Gaussian bw=0.05","Gaussian bw=0.075") ) )


###################################################
### chunk number 32: 
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

#	intL<-adapt(2, c(0,0), c(1,1), functn=function(x, alpha, beta)
#		{
#			exp(loglambda(x, alpha, beta))
#		},
#		alpha=alphabeta[1], 
#		beta=alphabeta[-1])
#	l<-l-intL$value
        intL<-adaptIntegrate(lowerLimit=c(0,0), upperLimit=c(1,1), fDim=1,
                tol=1e-8,
                f=function(x, alpha=alphabeta[1], beta=alphabeta[-1])
                {
                        exp(loglambda(x, alpha, beta))
                })
        l<-l-intL$integral
#	print(l)
	return(l) #Optim minimises
}


###################################################
### chunk number 34: 
###################################################
#library(adapt)
library(cubature)
data(lansing)
x<-as.points(lansing[lansing$marks=="maple",])


###################################################
### chunk number 35:  
###################################################
#Maximise log-likelihood
optbeta<-optim(par=c(log(514),0,0,0,0,0), #method="CG", 
 fn=L, control=list(maxit=1000, fnscale=-1), x=x)


###################################################
### chunk number 36: 
###################################################
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


###################################################
### chunk number 37: 
###################################################
data(lansing)
lmaple<-lansing[lansing$marks=="maple",]
ppm(Q=lmaple, trend=~x+y+I(x^2)+I(y^2)+I(x*y))


###################################################
### chunk number 38: 
###################################################
set.seed(30)


###################################################
### chunk number 40:  
###################################################
Kenvjap<-envelope(as(spjpines1, "ppp"), fun=Kest, r=r, nrank=2, nsim=99)
Kenvred<-envelope(as(spred, "ppp"), fun=Kest, r=r, nrank=2, nsim=99)
Kenvcells<-envelope(as(spcells, "ppp"), fun=Kest, r=r, nrank=2, nsim=99)
Kresults<-rbind(Kenvjap, Kenvred, Kenvcells)
Kresults<-cbind(Kresults, 
 y=rep(c("JAPANESE", "REDWOOD", "CELLS"), each=length(r)))


###################################################
### chunk number 41: 
###################################################
print(xyplot((obs-theo)~r|y, data=Kresults, type="l", 
   ylim= c(-.06, .06), ylab=expression(hat(K) (r)  - pi * r^2),
	panel=function(x, y, subscripts)
	{
		Ktheo<- Kresults$theo[subscripts]

		lpolygon(c(r, rev(r)), 
		   c(Kresults$lo[subscripts]-Ktheo, rev(Kresults$hi[subscripts]-Ktheo)),
		   border="gray", col="gray"
		)

		llines(r, Kresults$obs[subscripts]-Ktheo, lty=2, lwd=1.5, col="black")	
	}
))

###################################################
### chunk number 44: 
###################################################
library(maptools)
sG <- Sobj_SpatialGrid(spbdry, maxDim=50)$SG
gt <- slot(sG, "grid")
summary(gt)

###################################################
### chunk number 45: 
###################################################
pbdry <- slot(slot(slot(spbdry, "polygons")[[1]], "Polygons")[[1]], "coords")


###################################################
### chunk number 42:  
###################################################
source("cv.R")#Maybe change it to read it from the web
library(spatialkernel)

h<-seq(0.20, .4, by=.025)
cases<-spasthma[spasthma$Asthma=="case",]
ncases<-nrow(cases)
controls<-spasthma[spasthma$Asthma=="control",]
ncontrols<-nrow(controls)
cvres<-cvkernratio(rbind(coordinates(cases), coordinates(controls)), ncases, 
   ncontrols, pbdry, h, gt )

plot(h, cvres[[2]])

###################################################
### chunk number 43: 
###################################################
bwasthma <- .125 #.275


###################################################
### chunk number 46: 
###################################################
kcases<-spkernel2d(cases, pbdry, h0=bwasthma, gt)
kcontrols<-spkernel2d(controls, pbdry, h0=bwasthma, gt)


###################################################
### chunk number 47: 
###################################################
df0 <- data.frame(kcases=kcases, kcontrols=kcontrols)
splancs_ge_24 <- packageDescription("splancs")$Version > "2.01-23"
if (!splancs_ge_24) idxna <- complete.cases(df0)
spkratio0 <- SpatialGridDataFrame(gt, data=df0)
spkratio <- as(spkratio0, "SpatialPixelsDataFrame")
spkratio$kratio <- spkratio$kcases/spkratio$kcontrols
is.na(spkratio$kratio) <- !is.finite(spkratio$kratio)
spkratio$logratio <- log(spkratio$kratio)-log(ncases/ncontrols)


###################################################
### chunk number 48: 
###################################################
set.seed(1234)


###################################################
### chunk number 49: 
###################################################
if (splancs_ge_24) {
if (packageVersion("sp") < "1.1.0") {
    idxinbdry <- overlay(sG, spbdry)
  } else {
    idxinbdry <- over(sG, spbdry)
  }
  idxna <- !is.na(idxinbdry)
}


###################################################
### chunk number 50: 
###################################################
niter <- 99
ratio <- rep(NA, niter)
pvaluemap <- rep(0, sum(idxna))
rlabelratio <- matrix(NA, nrow=niter, ncol=sum(idxna))


###################################################
### chunk number 52: 
###################################################
set.seed(1)


###################################################
### chunk number 53:  
###################################################
for(i in 1:niter)
{
	idxrel <- sample(spasthma$Asthma) == "case"
	casesrel <- spasthma[idxrel,]
	controlsrel <- spasthma[!idxrel,]

	kcasesrel <- spkernel2d(casesrel, pbdry, h0=bwasthma, gt)
	kcontrolsrel <- spkernel2d(controlsrel, pbdry, h0=bwasthma, gt)
	kratiorel <- kcasesrel[idxna]/kcontrolsrel[idxna]
        is.na(kratiorel) <- !is.finite(kratiorel)
	rlabelratio[i,] <- kratiorel

	pvaluemap <- pvaluemap + (spkratio$kratio < kratiorel)
}


###################################################
### chunk number 54: 
###################################################
idxna2 <- apply(rlabelratio, 2, function(x) all(is.finite(x)))
rhomean <- apply(rlabelratio[, idxna2], 2, mean)
c <- prod(slot(gt, "cellsize"))
ratiorho <- c*sum((spkratio$kratio[idxna2]-ncases/ncontrols)^2)
ratio <- c*apply(rlabelratio[,idxna2], 1, 
 function(X, rho0 ){sum((X-rho0)^2)}, rho0=ncases/ncontrols
)
pvaluerho <- (sum(ratio > ratiorho)+1)/(niter+1)


###################################################
### chunk number 55: 
###################################################
spkratio$pvaluemap <- (pvaluemap+1)/(niter+1)
imgpvalue <- as.image.SpatialGridDataFrame(spkratio["pvaluemap"])
clpvalue <- contourLines(imgpvalue, levels=c(0,.05, .95, 1))
cl <- ContourLines2SLDF(clpvalue)


###################################################
### chunk number 56: 
###################################################
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
cols <- colorRampPalette(grey.colors(5, 0.95, 0.55, 2.2))(length(brks)-1)
colorkey<-list(labels=lbrks,
  at=(0:10)/10, height=.5)

print(spplot(spkratio, "kratio",
   col.regions=cols,
   do.log=TRUE, 
   colorkey=colorkey,
   at=c(0, brks[-c(1,11)], max(spkratio$kratio, na.rm=TRUE)),
   sp.layout=list(lyt05, lyt95, lytb, lytp) 
))


###################################################
### chunk number 57:  eval=FALSE
###################################################
source("cv.R")#Maybe change it to read it from the web

h<-seq(0.20, .4, by=.025)
cvresbin<-cvbinreg(rbind(coordinates(cases), coordinates(controls)), ncases,
   ncontrols, pbdry, h)

plot(h, cvresbin, type="l")#hopt=.225


###################################################
### chunk number 58: 
###################################################
bwasthmap <- 0.125 #0.225


###################################################
### chunk number 59: 
###################################################
#Estimator of the probability
lambda1 <- spkernel2d(cases, pbdry, h0=bwasthmap, gt)
lambda0 <- spkernel2d(controls, pbdry, h0=bwasthmap, gt)

lambda1 <- lambda1[idxna]
lambda0 <- lambda0[idxna]

spkratio$prob <- lambda1/(lambda1+lambda0)
is.na(spkratio$prob) <- !is.finite(spkratio$prob)


###################################################
### chunk number 60: 
###################################################
ats <- seq(0,max(spkratio$prob),length.out=11)
cols <- colorRampPalette(grey.colors(5, 0.9, 0.5, 2.2))(length(ats)-1)
print(spplot(spkratio, "prob", col.regions=cols, at=ats, sp.layout=list(lytb, lytp)))


###################################################
### chunk number 61: 
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
### chunk number 62: 
###################################################
library(mgcv)
### change for 1.5-2
if (packageDescription("mgcv")$Version != "1.5-2") {
  gasthma<-gam(y~1+dist1+dist2+dist3+droads+Genderf+Age+HayFeverf+smoking+s(x1,x2), data=spasthma[spasthma$Gender==1 | spasthma$Gender==2, ], family=binomial)
} else {
  spasthmadf <- as(spasthma, "data.frame")
  gasthma<-gam(y~1+dist1+dist2+dist3+droads+Genderf+Age+HayFeverf+smoking+s(x1,x2), data=spasthmadf[spasthmadf$Gender==1 | spasthmadf$Gender==2, ], family=binomial)
}

###################################################
### chunk number 63: 
###################################################
summary(gasthma)


###################################################
### chunk number 64: 
###################################################
sumGasth <- summary(gasthma)
cpv <- sumGasth$p.pv
spv <- sumGasth$s.table[4]


###################################################
### chunk number 65: 
###################################################
D2_mat <- as.matrix(spasthma$dist2)
RHO <- ncases/ncontrols
expsource2<-tribble(ccflag=spasthma$y, vars=D2_mat, rho=RHO, alphas=1, betas=1) 


###################################################
### chunk number 66: 
###################################################
print(expsource2)


###################################################
### chunk number 67: 
###################################################
#Hay fever
Hay_mat <- as.matrix(spasthma$HayFever)
exphay <- tribble(ccflag=spasthma$y, rho=RHO, covars=Hay_mat, thetas=1)


###################################################
### chunk number 68: 
###################################################
print(exphay)


###################################################
### chunk number 69: 
###################################################
expsource2hay<-tribble(ccflag=spasthma$y, vars=D2_mat, rho=RHO,
 alphas=1, betas=1, covars=Hay_mat, thetas=1)



###################################################
### chunk number 71: 
###################################################

s<-seq(0, .15, by=.01)


###################################################
### chunk number 72:  
###################################################
khcases<-khat(coordinates(cases), pbdry, s)
khcontrols<-khat(coordinates(controls), pbdry, s)
khcov<-khvmat(coordinates(cases), coordinates(controls), pbdry, s)

T0<-sum( ((khcases-khcontrols))/sqrt(diag(khcov)))


###################################################
### chunk number 73: 
###################################################
#Initialise seed for next chunk
set.seed(1234)


###################################################
### chunk number 74: 
###################################################
niter<-99
T<-rep(NA, niter)


###################################################
### chunk number 76:  eval=FALSE
###################################################
khcasesrel<-matrix(NA, nrow=length(s), ncol=niter)
khcontrolsrel<-matrix(NA, nrow=length(s), ncol=niter)

for(i in 1:niter)
{

	idxrel<-sample(spasthma$Asthma)=="case"
        casesrel<-coordinates(spasthma[idxrel,])
        controlsrel<-coordinates(spasthma[!idxrel,])

	khcasesrel[,i]<-khat(casesrel, pbdry, s)
	khcontrolsrel[,i]<-khat(controlsrel, pbdry, s)
	khdiff <- khcasesrel[,i]-khcontrolsrel[,i]
	T[i]<-sum(khdiff/sqrt(diag(khcov)))
}


###################################################
### chunk number 77: 
###################################################
pvalue<-(sum(T>T0)+1)/(niter+1)


###################################################
### chunk number 78: 
###################################################
plot(s, khcases-khcontrols, type="l", 
  ylab="D(s)", ylim=c(-.015, .015))#ylim=c(-11.5, 11.5))
lines(s, -1.96*sqrt(diag(khcov)), lty=2)
lines(s, +1.96*sqrt(diag(khcov)), lty=2)

#Compute envelopes
envel<-apply(khcasesrel-khcontrolsrel, 1, function(X){quantile(X, c(.025, .975))})
lines(s, envel[1,], lty=3, lwd=2)
lines(s, envel[2,], lty=3, lwd=2)

legend("bottomleft", 
   legend=c("Actual value", "Approx. 95% C.I.", "Sim. 95% envelopes"),
   lty=1:3, lwd=c(1,1,2), bty="n")


###################################################
### chunk number 79: 
###################################################
glmasthma<-glm(y~HayFeverf, data=spasthma, family="binomial")
prob<-fitted(glmasthma)
weights<-exp(glmasthma$linear.predictors)
library(spatialkernel)
setkernel("gaussian")
lambda0<- lambdahat(coordinates(controls), bwasthma, coordinates(cases), 
   pbdry, FALSE)$lambda
lambda1<- weights[spasthma$Asthma=="case"]*lambda0

ratiocc<-ncases/ncontrols
kihnocov<-kinhat(coordinates(cases), ratiocc*lambda0, pbdry,s)$k
kih<-kinhat(coordinates(cases), lambda1, pbdry,s)$k


###################################################
### chunk number 81: 
###################################################
set.seed(1234)


###################################################
### chunk number 82:  
###################################################
niter<-99

kinhomrelnocov<-matrix(NA, nrow=length(s), ncol=niter)
kinhomrel<-matrix(NA, nrow=length(s), ncol=niter)

for(i in 1:niter)
{
        idxrel<-sample(spasthma$Asthma, prob=prob)=="case"
        casesrel<-coordinates(spasthma[idxrel,])
        controlsrel<-coordinates(spasthma[!idxrel,])

	lambda0rel<-lambdahat(controlsrel, bwasthma, casesrel, pbdry, FALSE)$lambda
	lambda1rel<-weights[idxrel]*lambda0rel

	kinhomrelnocov[,i]<-kinhat(casesrel, ratiocc*lambda0rel, pbdry,s)$k
	kinhomrel[,i]<-kinhat(casesrel, lambda1rel, pbdry,s)$k
}


###################################################
### chunk number 83: 
###################################################
kinhsdnocov<-apply(kinhomrelnocov, 1, sd)
kihmeannocov<-apply(kinhomrelnocov, 1,mean)

D0nocov<-sum((kihnocov-kihmeannocov)/kinhsdnocov)
Dnocov<-apply(kinhomrelnocov, 2, 
   function(X){ sum((X-kihmeannocov)/kinhsdnocov)})

pvaluenocov<-(sum(Dnocov>D0nocov)+1)/(niter+1)


###################################################
### chunk number 84: 
###################################################
kinhsd<-apply(kinhomrel, 1, sd)
kihmean<-apply(kinhomrel, 1,mean)

D0<-sum((kih-kihmean)/kinhsd)
D<-apply(kinhomrel, 2, 
   function(X){ sum((X-kihmean)/kinhsd)})

pvalue<-(sum(D>D0)+1)/(niter+1)


###################################################
### chunk number 85: 
###################################################
oopar <- par(mfrow=c(1,2))
plot(s, kihnocov-kihmeannocov, type="l", 
   ylim= c(-0.06,  0.22),
   xlab="s", ylab= expression(hat(k)[I][","][hat(lambda)]-"E[s]"),
    main ="No covariates" )

envnocov<-apply(kinhomrelnocov, 1, function(X){quantile(X, c(.025, .975))})
lines(s, envnocov[1,]-kihmeannocov, lty=2)
lines(s, envnocov[2,]-kihmeannocov, lty=2)
plot(s, kih-kihmean, type="l", ylim=c(-0.06,  0.22), #c(-2e-4, 2e-4),
   xlab="s", ylab= expression(hat(k)[I][","][hat(lambda)]-"E[s]"),
   main ="Adjusting for Hay Fever"  )

env<-apply(kinhomrel, 1, function(X){quantile(X, c(.025, .975))})
lines(s, env[1,]-kihmean, lty=2)
lines(s, env[2,]-kihmean, lty=2)
par(oopar)



