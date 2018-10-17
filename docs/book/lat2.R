###################################################
### chunk number 1: 
###################################################
rm(list=ls())
library(digest)
if (!exists("chkDigest")) chkDigest <- TRUE
if (!exists("online")) online <- TRUE
owidth <- getOption("width")
options("width"=70)
owarn <- options("warn")$warn
options(warn=1)
odigits <- options("digits")$digits
options(digits=4)
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
## .epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".eps", sep="")
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
## .epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".png", sep="")
## png(filename=file, width = .iwidth, height = .iheight, pointsize = .ipointsize)


###################################################
### chunk number 7: zfig_png eval=FALSE
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 8: 
###################################################
#library(rgdal)
#NY8 <- readOGR("/home/rsb/tmp/lw", "NY8_utm18")
#NY8 <- readOGR(".", "NY8_1980")
#TCE <- readOGR("/home/rsb/tmp/lw", "TCE")
#library(spdep)
#NY_nb <- read.gal("/home/rsb/tmp/lw/NY_nb.gal", region.id=row.names(as(NY8, "data.frame")))
#NY_nb <- read.gal("NY8_book.gal")
#cities <- readOGR(".", "NY8cities")
library(digest)
library(sp)
intamap <- "http://intamap.geo.uu.nl/~roger/ASDAR/data"
if (online) {
  con <- url(paste(intamap, "NY8.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("NY8.RData")
if (chkDigest && !identical(digest(NY8), "954c3040dfa89c044eefd6cdc544bb95"))
  stop("NY8.RData digest error")
if (online)  {
  con <- url(paste(intamap, "TCE.RData", sep="/"), open="rb")
  load(con)
  close(con)
} else load("TCE.RData")
if (chkDigest && !identical(digest(TCE), "bad70492f37b5176dc6ff88679f8a51b"))
  stop("TCE.RData digest error")
if (online)  {
  con <- url(paste(intamap, "NY_nb.RData", sep="/"), open="rb")
  load(con)
  close(con)
} else load("NY_nb.RData")
if (chkDigest && !identical(digest(NY_nb), "d484686ba2fa0ab43e79cdcd3f340014"))
  stop("NY_nb.RData digest error")
if (online)  {
  con <- url(paste(intamap, "cities_NY8.RData", sep="/"), open="rb")
  load(con)
  close(con)
} else load("cities_NY8.RData")
if (chkDigest && !identical(digest(cities), "8f973621536e2bb794860564e7fb5bee"))
  stop("cities_NY8.RData digest error")


###################################################
### chunk number 9: 
###################################################
library(spdep)


###################################################
### chunk number 10: 
###################################################
nylm <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8)
summary(nylm)
NY8$lmresid <- residuals(nylm)


###################################################
### chunk number 11: 
###################################################
.iwidth <- 5
.iheight <- 4.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(lattice)
trellis.par.set(canonical.theme(color = FALSE))
library(RColorBrewer)
#gry <- c(rev(brewer.pal(6, "Reds")[1:3]), brewer.pal(6, "Blues"))
#gry <- brewer.pal(9, "Greys")[-1]
# RSB quietening greys
gry <- grey.colors(9, 0.95, 0.55, 2.2)
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(NY8, "lmresid", sp.layout=list(TCEpts), col.regions=colorRampPalette(gry)(7), at=seq(-2.5,5,length.out=6))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 12: 
###################################################
NYlistw<-nb2listw(NY_nb, style = "B")
lm.morantest(nylm, NYlistw)


###################################################
### chunk number 13: 
###################################################
nysar<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME , data=NY8, listw=NYlistw)
summary(nysar)


###################################################
### chunk number 14: 
###################################################
nylam1 <- c(nysar$lambda)
nylam2 <- c(LR1.spautolm(nysar)$p.value)


###################################################
### chunk number 15: 
###################################################
nylmw <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, weights=POP8)
summary(nylmw)
NY8$lmwresid <- residuals(nylmw)


###################################################
### chunk number 16: 
###################################################
.iwidth <- 5
.iheight <- 4.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(RColorBrewer)
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(NY8, "lmwresid", sp.layout=list(TCEpts), col.regions=colorRampPalette(gry)(7), at=seq(-2.5,5,length.out=6))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 17: 
###################################################
lm.morantest(nylmw, NYlistw)


###################################################
### chunk number 18: 
###################################################
nysarw<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME , data=NY8, listw=NYlistw, weights=POP8)
summary(nysarw)


###################################################
### chunk number 19: 
###################################################
nycar<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME , data=NY8, family="CAR",
   listw=NYlistw)
summary(nycar)


###################################################
### chunk number 20: 
###################################################
nylam1 <- c(nycar$lambda)


###################################################
### chunk number 21: 
###################################################
nycarw<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="CAR",
   listw=NYlistw, weights=POP8)
summary(nycarw)


###################################################
### chunk number 22: 
###################################################
nysarwM<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SAR",
   listw=NYlistw, weights=POP8, method="Matrix")


###################################################
### chunk number 23: 
###################################################
summary(nysarwM)


###################################################
### chunk number 24: 
###################################################
1/range(eigenw(NYlistw))
nysar_ll<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SAR",
   listw=NYlistw, llprof=100)
nysarw_ll<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SAR",
   listw=NYlistw, weights=POP8, llprof=100)


###################################################
### chunk number 25: 
###################################################
.iwidth <- 6
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
ylim <- range(c(nysarw_ll$llprof$ll, nysar_ll$llprof$ll), na.rm=TRUE)
plot(nysarw_ll$llprof$lambda, nysarw_ll$llprof$ll, type="l", xlab=expression(lambda), ylab="log likelihood", ylim=ylim, lwd=2)
abline(v=nysarw_ll$lambda)
abline(h=nysarw_ll$LL)
lines(nysar_ll$llprof$lambda, nysar_ll$llprof$ll, lty=2, lwd=2)
abline(v=nysar_ll$lambda, lty=2)
abline(h=nysar_ll$LL, lty=2)
legend("bottom", legend=c("weighted SAR", "SAR"), lty=c(1,2), lwd=2, bty="n")
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 26: 
###################################################
nysmaw<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SMA",
   listw=NYlistw, weights=POP8)
summary(nysmaw)


###################################################
### chunk number 27: 
###################################################
#This example tries to fit random effects assuming that there is
#a spatial correlation between the random effects.
#
#Note that this is slightly different from kriging because we
#have the random effects PLUS the error term, whilst the kriging
#only have the spatially correlated 'error term'. (is this right?)

library(nlme)

NY8$x<-coordinates(NY8)[,1]/1000
NY8$y<-coordinates(NY8)[,2]/1000

#First value is the range value (see corExp) and second nugget effect
sp1 <- corSpatial(1, form = ~ x + y, type = "gaussian") 
scor<-Initialize(sp1, as(NY8, "data.frame")[,c("x", "y")], nugget=FALSE)


###################################################
### chunk number 28: 
###################################################
if (online) {
  con <- url(paste(intamap, "spmodel.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("spmodel.RData")
if (chkDigest && !identical(digest(spmodel), "94aa16e985115e0f5695de7c5d71585f"))
  stop("spmodel.RData digest error")


###################################################
### chunk number 29:  eval=FALSE
###################################################
## spmodel<-lme(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, random=~1|AREAKEY, data=as(NY8, "data.frame"), correlation=scor, method="ML")


###################################################
### chunk number 30: 
###################################################
summary(spmodel)

#Checking correlation structure...
#scorm<-as.matrix(scor)
#scorm[1:5,1:5]
#
#library(spatstat)
#pd<-pairdist(ppp(darea$x, darea$y))
#
#pdcor<-exp(-pd)
#pdcor[1:5, 1:5]
#
#save(spmodel, file="spmodel.RData")


###################################################
### chunk number 31: 
###################################################
library(lmtest)
bptest(nylm)


###################################################
### chunk number 32: 
###################################################
library(sandwich)
coeftest(nylm)
coeftest(nylm, vcov=vcovHC(nylm, type="HC4"))


###################################################
### chunk number 33: 
###################################################
NYlistwW <- nb2listw(NY_nb, style = "W")
res <- lm.LMtests(nylm, listw=NYlistwW, test="all")
tres <- t(sapply(res, function(x) c(x$statistic, x$parameter, x$p.value)))
colnames(tres) <- c("Statistic", "df", "p-value")
printCoefmat(tres)


###################################################
### chunk number 34: 
###################################################
nylag <- lagsarlm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW)
summary(nylag)
bptest.sarlm(nylag)


###################################################
### chunk number 35: 
###################################################
nymix <- lagsarlm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW, type="mixed")
nymix
anova(nymix, nylag)


###################################################
### chunk number 36: 
###################################################
nyerr <- errorsarlm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW)
summary(nyerr)


###################################################
### chunk number 37: 
###################################################
nystsls <- stsls(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW)
summary(nystsls)


###################################################
### chunk number 38: 
###################################################
nystslsR <- stsls(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW, robust=TRUE)
summary(nystslsR)


###################################################
### chunk number 39: 
###################################################
nyGMerr <- GMerrorsar(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW)
summary(nyGMerr)


###################################################
### chunk number 40: 
###################################################
vv <- nyGMerr$vv
lambda <- nyGMerr$lambda
s2 <- nyGMerr$s2
lseq <- seq(0, 0.4, 0.01)
s2seq <- seq(0.3, 0.55, 0.01)
lsseq <- as.matrix(expand.grid(lseq, s2seq))
res <- numeric(nrow(lsseq))
for (i in seq(along=res)) res[i] <- spdep:::.kpgm(lsseq[i,,drop=TRUE], v=vv, verbose=TRUE)
SGDF <- SpatialPixelsDataFrame(lsseq, data=data.frame(fn=res))
fullgrid(SGDF) <- TRUE


###################################################
### chunk number 41: 
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
#library(RColorBrewer)
#gry <- brewer.pal(5, "Greys")
image(SGDF, "fn", col=grey.colors(10, 1, 0.55, 2.2), axes=TRUE)
title(xlab=expression(lambda), ylab=expression(sigma^2))
contour(SGDF, "fn", add=TRUE)
points(c(lambda, nyerr$lambda), c(s2, nyerr$s2), pch=c(4, 3), lwd=2)
text(c(lambda, nyerr$lambda), c(s2, nyerr$s2), labels=c("GM", "ML"), pos=c(2, 4), offset=0.5)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 42: 
###################################################
library(mgcv)
nyGAM1 <- gam(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME+s(x,y), weights=POP8, data=NY8)
anova(nylmw, nyGAM1, test="Chisq")


###################################################
### chunk number 43: 
###################################################
nylam1 <- c(summary(nyGAM1)$edf)


###################################################
### chunk number 44: 
###################################################
nyGLMp <- glm(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8, family="poisson")


###################################################
### chunk number 45:  eval=FALSE
###################################################
## summary(nyGLMp)


###################################################
### chunk number 46: 
###################################################
xx <- capture.output(print(summary(nyGLMp)))
zz <- strwrap(paste(xx[3:4], collapse=" "))
xx[3] <- zz[1]
xx[4] <- paste("   ", zz[2])
cat(xx, sep="\n")


###################################################
### chunk number 47: 
###################################################
NY8$lmpresid <- residuals(nyGLMp, type="deviance")
lm.morantest(nyGLMp, listw=NYlistwW)


###################################################
### chunk number 48: 
###################################################
.iwidth <- 5
.iheight <- 4.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(RColorBrewer)
#gry <- c(rev(brewer.pal(6, "Reds")[1:3]), brewer.pal(6, "Blues"))
#gry <- grey.colors(9, 0.95, 0.55, 2.2)
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(NY8, "lmpresid", sp.layout=list(TCEpts), col.regions=grey.colors(9, 0.95, 0.55, 2.2), at=seq(-3,5,length.out=9))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 49: 
###################################################
nyGAMp <- gam(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8))+s(x,y), data=NY8, family="poisson")
summary(nyGAMp)


###################################################
### chunk number 50:  eval=FALSE
###################################################
## anova(nyGLMp, nyGAMp, test="Chisq")


###################################################
### chunk number 51: 
###################################################
xx <- capture.output(print(anova(nyGLMp, nyGAMp, test="Chisq")))
cat(xx[1:2], sep="\n")
cat(strwrap(xx[3], width=70, exdent=4), sep="\n")
cat(strwrap(paste(xx[4:5], collapse=" "), width=70, exdent=4), sep="\n")
cat(xx[6:8], sep="\n")


###################################################
### chunk number 52: 
###################################################
nylam1 <- c(summary(nyGAMp)$edf)


###################################################
### chunk number 53: 
###################################################
if (online) {
  con <- url(paste(intamap, "nyGLMMp.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("nyGLMMp.RData")
#if (chkDigest && !identical(digest(nyGLMMp), "3bde1bbc641250b8636a4d6442a5ebf4"))
  #stop("nyGLMMp.RData digest error")


###################################################
### chunk number 54:  eval=FALSE
###################################################
## library(MASS)
## attach(as(NY8, "data.frame"))
## nyGLMMp <- glmmPQL(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8, family=poisson, random=~1|AREAKEY, correlation=scor)
## detach("as(NY8, \"data.frame\")")


###################################################
### chunk number 55:  eval=FALSE
###################################################
## summary(nyGLMMp)


###################################################
### chunk number 56: 
###################################################
xx <- capture.output(print(summary(nyGLMMp)))
cat(xx[1:18], sep="\n")
cat(strwrap(xx[19], width=70, exdent=4), sep="\n")
cat(xx[20:36], sep="\n")


###################################################
### chunk number 57: 
###################################################
nySFE <- SpatialFiltering(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, nb=NY_nb, style="W", verbose=FALSE)
nylmSFE <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME+fitted(nySFE), data=NY8)
summary(nylmSFE)
anova(nylm, nylmSFE)


###################################################
### chunk number 58: 
###################################################
if (online) {
  con <- url(paste(intamap, "nyME.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("nyME.RData")
if (chkDigest && !identical(digest(nyME), "0bf26c87da7fa807dad3950e8694b9d2"))
  stop("nyME.RData digest error")
#save(nyME, file="nyME.RData")


###################################################
### chunk number 59:  eval=FALSE
###################################################
## nyME <- ME(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, offset=log(POP8), family="poisson", listw=NYlistwW, alpha=0.5)


###################################################
### chunk number 60: 
###################################################
nyME


###################################################
### chunk number 61: 
###################################################
NY8$eigen_24 <- fitted(nyME)[,1]
NY8$eigen_223 <- fitted(nyME)[,2]


###################################################
### chunk number 62: 
###################################################
.iwidth <- 6
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(RColorBrewer)
#gry <- brewer.pal(9, "Greys")[-1]
spplot(NY8, c("eigen_24", "eigen_223"), col.regions=grey.colors(6, 0.95, 0.55, 2.2), cuts=5)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 63: 
###################################################
nyglmME <- glm(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8))+fitted(nyME), data=NY8, family="poisson")


###################################################
### chunk number 64:  eval=FALSE
###################################################
## summary(nyglmME)


###################################################
### chunk number 65: 
###################################################
xx <- capture.output(print(summary(nyglmME)))
zz <- strwrap(paste(xx[3:4], collapse=" "), width=70, exdent=4)
xx[3] <- zz[1]
xx[4] <- zz[2]
cat(xx, sep="\n")


###################################################
### chunk number 66:  eval=FALSE
###################################################
## anova(nyGLMp, nyglmME, test="Chisq")


###################################################
### chunk number 67: 
###################################################
xx <- capture.output(print(anova(nyGLMp, nyglmME, test="Chisq")))
cat(xx[1:2], sep="\n")
cat(strwrap(xx[3], width=70, exdent=4), sep="\n")
cat(strwrap(paste(xx[4:5], collapse=" "), width=70, exdent=4), sep="\n")
cat(xx[6:8], sep="\n")


###################################################
### chunk number 68: 
###################################################
if (online) {
  con <- url(paste(intamap, "nyGWR.RData", sep="/"), open="rb")
  load(con)
  close(con)
}  else load("nyGWR.RData")
if (chkDigest && !identical(digest(gwrG), "8754940831febb7c7356c4a1cbfd8e2e"))
  stop("nyGWR.RData digest error")
#save(bwG, gwrG, gbwG, ggwrG, file="nyGWR.RData")


###################################################
### chunk number 69: 
###################################################
library(spgwr)


###################################################
### chunk number 70:  eval=FALSE
###################################################
## bwG <- gwr.sel(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, gweight=gwr.Gauss, verbose=FALSE)
## gwrG <- gwr(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, bandwidth=bwG, gweight=gwr.Gauss, hatmatrix=TRUE)


###################################################
### chunk number 71: 
###################################################
gwrG


###################################################
### chunk number 72:  eval=FALSE
###################################################
## gbwG <- ggwr.sel(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8, family="poisson", gweight=gwr.Gauss, verbose=FALSE)
## ggwrG <- ggwr(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8, family="poisson", bandwidth=gbwG, gweight=gwr.Gauss)


###################################################
### chunk number 73:  eval=FALSE
###################################################
## ggwrG


###################################################
### chunk number 74: 
###################################################
xx <- capture.output(print(ggwrG))
cat(xx[1], sep="\n")
cat(strwrap(paste(xx[2:3], collapse=" "), width=70, exdent=4), sep="\n")
cat(xx[4:11], sep="\n")


###################################################
### chunk number 75: 
###################################################
.iwidth <- 4.5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
#library(RColorBrewer)
#gry <- brewer.pal(9, "Greys")[-1]
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(ggwrG$SDF, "PEXPOSURE", sp.layout=list(TCEpts), col.regions=grey.colors(7, 0.95, 0.55, 2.2), cuts=6)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 76: 
###################################################
.iwidth <- 4.5
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-lat2-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
pairs(as(ggwrG$SDF, "data.frame")[,2:5])
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 77: 
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
sT <- capture.output(print(Sys.time()))
cat("%", sT, "\n")


###################################################
### chunk number 78: 
###################################################
options(warn=owarn)
options(width=owidth)
options(digits=odigits)


