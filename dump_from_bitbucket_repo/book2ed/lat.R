### R code from vignette source 'lat.Rnw'

###################################################
### code chunk number 1: lat.Rnw:8-14
###################################################
if (!exists("book_R_dont_trash")) rm(list=ls())
require(RColorBrewer)
pal = function(n) brewer.pal(n, "Reds")
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
### code chunk number 3: lat.Rnw:23-24
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 4: afig (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)


###################################################
### code chunk number 5: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 6: afig_png (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".png", sep="")
## png(filename=file, width = .iwidth, height = .iheight, pointsize = .ipointsize)


###################################################
### code chunk number 7: zfig_png (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 8: lat.Rnw:179-187
###################################################
library(rgdal)
setwd("../Data")
NY8 <- readOGR(".", "NY8_utm18")
TCE <- readOGR(".", "TCE")
library(spdep)
NY_nb <- read.gal("NY_nb.gal", region.id=row.names(NY8))
cities <- readOGR(".", "NY8cities")
setwd("../lat")


###################################################
### code chunk number 9: lat.Rnw:193-210
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1,2), mar=c(3,3,1,1)+0.1)
plot(NY8, border="grey60", axes=TRUE)
text(coordinates(cities), labels=as.character(cities$names), font=2, cex=0.9)
text(bbox(NY8)[1,1], bbox(NY8)[2,2], labels="a)", cex=0.8)
plot(NY8, border="grey60", axes=TRUE)
points(TCE, pch=1, cex=0.7)
points(TCE, pch=3, cex=0.7)
text(coordinates(TCE), labels=as.character(TCE$name), cex=0.7,
 font=1, pos=c(4,1,4,1,4,4,4,2,3,4,2), offset=0.3)
text(bbox(NY8)[1,1], bbox(NY8)[2,2], labels="b)", cex=0.8)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 10: lat.Rnw:281-282 (eval = FALSE)
###################################################
## vignette("nb", package="spdep")


###################################################
### code chunk number 11: lat.Rnw:310-311
###################################################
library(spdep)


###################################################
### code chunk number 12: lat.Rnw:313-316 (eval = FALSE)
###################################################
## library(rgdal)
## NY8 <- readOGR(".", "NY8_utm18")
## NY_nb <- read.gal("NY_nb.gal", region.id=row.names(NY8))


###################################################
### code chunk number 13: lat.Rnw:318-319
###################################################
summary(NY_nb)


###################################################
### code chunk number 14: lat.Rnw:321-323 (eval = FALSE)
###################################################
## plot(NY8, border="grey60")
## plot(NY_nb, coordinates(NY8), pch=19, cex=0.6, add=TRUE)


###################################################
### code chunk number 15: lat.Rnw:333-344
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 10
.pwd <- 0.75
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(3,3,1,1)+0.1)
plot(NY8, border="grey60", axes=TRUE)
plot(NY_nb, coordinates(NY8), pch=19, cex=0.6, add=TRUE)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 16: lat.Rnw:379-382
###################################################
Syracuse <- NY8[NY8$AREANAME == "Syracuse city",]
Sy0_nb <- subset(NY_nb, NY8$AREANAME == "Syracuse city")
summary(Sy0_nb)


###################################################
### code chunk number 17: lat.Rnw:408-415
###################################################
coords <- coordinates(Syracuse)
IDs <- row.names(Syracuse)
Sy8_nb <- knn2nb(knearneigh(coords, k=1), row.names=IDs)
Sy9_nb <- knn2nb(knearneigh(coords, k=2), row.names=IDs)
Sy10_nb <- knn2nb(knearneigh(coords, k=4), row.names=IDs)
dsts <- unlist(nbdists(Sy8_nb, coords))
Sy11_nb <- dnearneigh(coords, d1=0, d2=0.75*max(dsts), row.names=IDs)


###################################################
### code chunk number 18: lat.Rnw:495-499
###################################################
Sy0_lw_W <- nb2listw(Sy0_nb)
Sy0_lw_W
names(Sy0_lw_W)
names(attributes(Sy0_lw_W))


###################################################
### code chunk number 19: lat.Rnw:519-522
###################################################
1/rev(range(card(Sy0_lw_W$neighbours)))
summary(unlist(Sy0_lw_W$weights))
summary(sapply(Sy0_lw_W$weights, sum))


###################################################
### code chunk number 20: lat.Rnw:539-542
###################################################
Sy0_lw_B <- nb2listw(Sy0_nb, style="B")
summary(unlist(Sy0_lw_B$weights))
summary(sapply(Sy0_lw_B$weights, sum))


###################################################
### code chunk number 21: lat.Rnw:579-584
###################################################
dsts <- nbdists(Sy0_nb, coordinates(Syracuse))
idw <- lapply(dsts, function(x) 1/(x/1000))
Sy0_lw_idwB <- nb2listw(Sy0_nb, glist=idw, style="B")
summary(unlist(Sy0_lw_idwB$weights))
summary(sapply(Sy0_lw_idwB$weights, sum))


###################################################
### code chunk number 22: lat.Rnw:590-618
###################################################
.iwidth <- 5
.iheight <- 2.5
.ipointsize <- 10
.pwd <- 0.95
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(RColorBrewer)
pal <- brewer.pal(9, "Reds")
# RSB quietening greys
#pal <- grey.colors(9, 1, 0.5, 2.2)
oopar <- par(mfrow=c(1,3), mar=c(1,1,3,1)+0.1)
z <- t(listw2mat(Sy0_lw_W))
brks <- c(0,0.1,0.143,0.167,0.2,0.5,1)
nbr3 <- length(brks)-3
image(1:63, 1:63, z[,ncol(z):1], breaks=brks, col=pal[c(1,(9-nbr3):9)],
 main="W style", axes=FALSE)
box()
z <- t(listw2mat(Sy0_lw_B))
image(1:63, 1:63, z[,ncol(z):1], col=pal[c(1,9)], main="B style", axes=FALSE)
box()
z <- t(listw2mat(Sy0_lw_idwB))
brks <- c(0,0.35,0.73,0.93,1.2,2.6)
nbr3 <- length(brks)-3
image(1:63, 1:63, z[,ncol(z):1], breaks=brks, col=pal[c(1,(9-nbr3):9)],
 main="IDW B style", axes=FALSE)
box()
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 23: lat.Rnw:666-667 (eval = FALSE)
###################################################
## Sy0_lw_D1 <- nb2listw(Sy11_nb, style="B")


###################################################
### code chunk number 24: lat.Rnw:669-671
###################################################
t1 <- try(Sy0_lw_D1 <- nb2listw(Sy11_nb, style="B"), silent=TRUE)
cat(t1)


###################################################
### code chunk number 25: lat.Rnw:673-675
###################################################
Sy0_lw_D1 <- nb2listw(Sy11_nb, style="B", zero.policy=TRUE)
print(Sy0_lw_D1, zero.policy=TRUE)


###################################################
### code chunk number 26: lat.Rnw:744-746
###################################################
Sy14_nb <- read.gal("Sy_GeoDa1.GAL")
isTRUE(all.equal(Sy0_nb, Sy14_nb, check.attributes=FALSE))


###################################################
### code chunk number 27: lat.Rnw:775-777
###################################################
Sy16_nb <- read.gwt2nb("Sy_GeoDa4.GWT")
isTRUE(all.equal(Sy10_nb, Sy16_nb, check.attributes=FALSE))


###################################################
### code chunk number 28: lat.Rnw:882-887
###################################################
set.seed(987654)
n <- length(Sy0_nb)
uncorr_x <- rnorm(n)
rho <- 0.5
autocorr_x <- invIrW(Sy0_lw_W, rho) %*% uncorr_x


###################################################
### code chunk number 29: lat.Rnw:893-908
###################################################
.iwidth <- 5
.iheight <- 2.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1,2), mar=c(4,4,3,2)+0.1)
plot(uncorr_x, lag(Sy0_lw_W, uncorr_x), xlab="", cex.lab=0.8,
 ylab="spatial lag", main="Uncorrelated random variable", cex.main=0.8)
lines(lowess(uncorr_x, lag(Sy0_lw_W, uncorr_x)), lty=2, lwd=2)
plot(autocorr_x, lag(Sy0_lw_W, autocorr_x),
 xlab="", ylab="",
 main="Autocorrelated random variable", cex.main=0.8, cex.lab=0.8)
lines(lowess(autocorr_x, lag(Sy0_lw_W, autocorr_x)), lty=2, lwd=2)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 30: lat.Rnw:1060-1063
###################################################
moran_u <- moran.test(uncorr_x, listw=Sy0_lw_W)
moran_a <- moran.test(autocorr_x, listw=Sy0_lw_W)
moran_a1 <- moran.test(autocorr_x, listw=nb2listw(Sy9_nb, style="W"))


###################################################
### code chunk number 31: lat.Rnw:1086-1090
###################################################
et <- coords[,1] - min(coords[,1])
trend_x <- uncorr_x + 0.00025 * et
moran_t <- moran.test(trend_x, listw=Sy0_lw_W)
moran_t1 <- lm.morantest(lm(trend_x ~ et), listw=Sy0_lw_W)


###################################################
### code chunk number 32: lat.Rnw:1100-1121
###################################################
res1 <- matrix("", ncol=5, nrow=5)
rownames(res1) <- c("uncorr\\_x", "autocorr\\_x", "autocorr\\_x k=1", "trend\\_x", "lm(trend\\_x $\\sim$ et)")
colnames(res1) <- c("$I$", "$E(I)$", "$var(I)$", "st. deviate", "$p$-value")
res1[1, 1:3] <- format(moran_u$estimate, digits=3)
res1[1, 4] <- format(moran_u$statistic, digits=3)
res1[1, 5] <- format.pval(moran_u$p.value, digits=2, eps=1e-8)
res1[2, 1:3] <- format(moran_a$estimate, digits=3)
res1[2, 4] <- format(moran_a$statistic, digits=3)
res1[2, 5] <- format.pval(moran_a$p.value, digits=2, eps=1e-8)
res1[3, 1:3] <- format(moran_a1$estimate, digits=3)
res1[3, 4] <- format(moran_a1$statistic, digits=3)
res1[3, 5] <- format.pval(moran_a1$p.value, digits=2, eps=1e-8)
res1[4, 1:3] <- format(moran_t$estimate, digits=3)
res1[4, 4] <- format(moran_t$statistic, digits=3)
res1[4, 5] <- format.pval(moran_t$p.value, digits=2, eps=1e-8)
res1[5, 1:3] <- format(moran_t1$estimate, digits=3)
res1[5, 4] <- format(moran_t1$statistic, digits=3)
res1[5, 5] <- format.pval(moran_t1$p.value, digits=2, eps=1e-8)
library(xtable)
print(xtable(res1, align=c("l", rep("r", 5))), floating=FALSE,
 sanitize.text.function = function(SANITIZE) SANITIZE)


###################################################
### code chunk number 33: lat.Rnw:1239-1241
###################################################
options("width"=80, "digits"=4)
K <- moran(NY8$Cases, listw=nb2listw(NY_nb, style="B"), n=length(NY8$Cases), S0=Szero(nb2listw(NY_nb, style="B")))$K


###################################################
### code chunk number 34: lat.Rnw:1262-1263
###################################################
moran.test(NY8$Cases, listw=nb2listw(NY_nb))


###################################################
### code chunk number 35: lat.Rnw:1286-1288
###################################################
lw_B <- nb2listw(NY_nb, style="B")
moran.test(NY8$Cases, listw=lw_B)


###################################################
### code chunk number 36: lat.Rnw:1316-1317
###################################################
moran.test(NY8$Cases, listw=lw_B, randomisation=FALSE)


###################################################
### code chunk number 37: lat.Rnw:1343-1344
###################################################
lm.morantest(lm(Cases ~ 1, NY8), listw=lw_B)


###################################################
### code chunk number 38: lat.Rnw:1375-1376
###################################################
options("width"=60)


###################################################
### code chunk number 39: lat.Rnw:1378-1379
###################################################
lm.morantest.sad(lm(Cases ~ 1, NY8), listw=lw_B)


###################################################
### code chunk number 40: lat.Rnw:1381-1382
###################################################
options("width"=70)


###################################################
### code chunk number 41: lat.Rnw:1384-1385
###################################################
lm.morantest.exact(lm(Cases ~ 1, NY8), listw=lw_B)


###################################################
### code chunk number 42: lat.Rnw:1401-1404
###################################################
set.seed(1234)
bperm <- moran.mc(NY8$Cases, listw=lw_B, nsim=999)
bperm


###################################################
### code chunk number 43: lat.Rnw:1426-1433
###################################################
r <- sum(NY8$Cases)/sum(NY8$POP8)
rni <- r*NY8$POP8
CR <- function(var, mle) rpois(length(var), lambda=mle)
MoranI.pboot <- function(var, i, listw, n, S0, ...) {
  return(moran(x=var, listw=listw, n=n, S0=S0)$I)
}
set.seed(1234)


###################################################
### code chunk number 44: lat.Rnw:1435-1436
###################################################
options("width"=60)


###################################################
### code chunk number 45: lat.Rnw:1438-1440
###################################################
boot2 <- boot(NY8$Cases, statistic=MoranI.pboot, R=999, sim="parametric",
  ran.gen=CR, listw=lw_B, n=length(NY8$Cases), S0=Szero(lw_B), mle=rni)


###################################################
### code chunk number 46: lat.Rnw:1442-1443
###################################################
options("width"=70)


###################################################
### code chunk number 47: lat.Rnw:1445-1446
###################################################
pnorm((boot2$t0 - mean(boot2$t))/sd(boot2$t[,1]), lower.tail=FALSE)


###################################################
### code chunk number 48: lat.Rnw:1453-1467
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1,2))
xlim <- range(c(bperm$res, boot2$t[,1]))
hist(bperm$res[-length(bperm$res)], main="Permutation bootstrap", xlab=expression(I[std]), xlim=xlim, density=15, angle=45, ylim=c(0,260))
abline(v=bperm$statistic, lty=2)
hist(boot2$t, col=rgb(0.4,0.4,0.4), main="Parametric bootstrap", xlab=expression(I[CR]), xlim=xlim, ylim=c(0,260))
hist(bperm$res[-length(bperm$res)], density=15, angle=45, add=TRUE)
abline(v=boot2$t0, lty=2)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 49: lat.Rnw:1487-1488 (eval = FALSE)
###################################################
## rni <- fitted(glm(Cases ~ 1 + offset(log(POP8)), data=NY8, family="poisson"))


###################################################
### code chunk number 50: lat.Rnw:1521-1523
###################################################
set.seed(1234)
EBImoran.mc(n=NY8$Cases, x=NY8$POP8, listw=nb2listw(NY_nb, style="B"), nsim=999)


###################################################
### code chunk number 51: lat.Rnw:1557-1558
###################################################
cor8 <- sp.correlogram(neighbours=NY_nb, var=NY8$Cases, order=8, method="I", style="C")


###################################################
### code chunk number 52: lat.Rnw:1579-1581
###################################################
library(pgirmess)
corD <- correlog(coordinates(NY8), NY8$Cases, method="Moran")


###################################################
### code chunk number 53: lat.Rnw:1589-1599
###################################################
.iwidth <- 6
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1,2))
plot(cor8, main="Contiguity lag orders")
plot(corD, main="Distance bands")
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 54: lat.Rnw:1631-1655
###################################################
.iwidth <- 6
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1,2))
msp <- moran.plot(NY8$Cases, listw=nb2listw(NY_nb, style="C"), quiet=TRUE)
title("Moran scatterplot")
infl <- apply(msp$is.inf, 1, any)
x <- NY8$Cases
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L", "H"), include.lowest=TRUE)
wx <- lag(nb2listw(NY_nb, style="C"), NY8$Cases)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)), labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4
plot(NY8, col=brewer.pal(4, "Accent")[cols])
# RSB quietening greys
legend("topright", legend=c("None", "HL", "LH", "HH"), fill=brewer.pal(4, "Accent"), bty="n", cex=0.8, y.intersp=0.8)
title("Tracts with influence")
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 55: lat.Rnw:1683-1684 (eval = FALSE)
###################################################
## moran.plot(NY8$Cases, listw=nb2listw(NY_nb, style="C"))


###################################################
### code chunk number 56: lat.Rnw:1764-1767
###################################################
lm1 <- localmoran(NY8$Cases, listw=nb2listw(NY_nb, style="C"))
lm2 <- as.data.frame(localmoran.sad(lm(Cases ~ 1, NY8), nb=NY_nb, style="C"))
lm3 <- as.data.frame(localmoran.exact(lm(Cases ~ 1, NY8), nb=NY_nb, style="C"))


###################################################
### code chunk number 57: lat.Rnw:1777-1783
###################################################
r <- sum(NY8$Cases)/sum(NY8$POP8)
rni <- r*NY8$POP8
lw <- nb2listw(NY_nb, style="C")
sdCR <- (NY8$Cases - rni)/sqrt(rni)
wsdCR <- lag(lw, sdCR)
I_CR <- sdCR * wsdCR


###################################################
### code chunk number 58: lat.Rnw:1796-1811
###################################################
.iwidth <- 6
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
library(RColorBrewer)
gry <- c(rev(brewer.pal(8, "Reds")[1:6]), brewer.pal(6, "Blues"))
NY8$Standard <- lm1[,1]
NY8$"Constant_risk" <- I_CR
nms <- match(c("Standard", "Constant_risk"), names(NY8))
spplot(NY8, nms, at=c(-2.5,-1.4,-0.6,-0.2,0,0.2,0.6,4,7), col.regions=colorRampPalette(gry)(8))
# RSB quietening greys
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 59: lat.Rnw:1857-1871
###################################################
set.seed(1234)
nsim <- 999
N <- length(rni)
sims <- matrix(0, ncol=nsim, nrow=N)
for (i in 1:nsim) {
  y <- rpois(N, lambda=rni)
  sdCRi <- (y - rni)/sqrt(rni)
  wsdCRi <- lag(lw, sdCRi)
  sims[,i] <- sdCRi * wsdCRi 
}
xrank <- apply(cbind(I_CR, sims), 1, function(x) rank(x)[1])
diff <- nsim - xrank
diff <- ifelse(diff > 0, diff, 0)
pval <- punif((diff + 1)/(nsim + 1))


###################################################
### code chunk number 60: lat.Rnw:1877-1893
###################################################
.iwidth <- 6
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
NY8$Normal <- lm2[,3]
NY8$Randomisation <- lm1[,5]
NY8$Saddlepoint <- lm2[,5]
NY8$Exact <- lm3[,5]
NY8$Constant_risk <- pval
gry <- c(rev(brewer.pal(6, "Reds")), brewer.pal(6, "Blues"))
spplot(NY8, c("Normal", "Randomisation", "Saddlepoint", "Exact", "Constant_risk"), at=c(0,0.01,0.05,0.1,0.9,0.95,0.99,1), col.regions=colorRampPalette(gry)(7))
# RSB quietening greys
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 61: lat.Rnw:1904-1914
###################################################
.iwidth <- 6
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
spplot(NY8, c("Normal", "Exact", "Constant_risk"), xlim=c(405200, 432200), ylim=c(4652700, 4672000), at=c(0,0.01,0.05,0.1,0.9,0.95,0.99,1), col.regions=colorRampPalette(gry)(7))
# RSB quietening greys
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 62: lat.Rnw:2083-2086
###################################################
nylm <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8)
summary(nylm)
NY8$lmresid <- residuals(nylm)


###################################################
### code chunk number 63: lat.Rnw:2129-2130
###################################################
library(spdep)


###################################################
### code chunk number 64: lat.Rnw:2132-2134
###################################################
NYlistw<-nb2listw(NY_nb, style = "B")
lm.morantest(nylm, NYlistw)


###################################################
### code chunk number 65: lat.Rnw:2155-2158
###################################################
NYlistwW <- nb2listw(NY_nb, style = "W")
aple(residuals(nylm), listw=NYlistwW)
spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW)$lambda


###################################################
### code chunk number 66: lat.Rnw:2276-2278
###################################################
nysar<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistw)
summary(nysar)


###################################################
### code chunk number 67: lat.Rnw:2281-2283
###################################################
nylam1 <- c(nysar$lambda)
nylam2 <- c(LR1.spautolm(nysar)$p.value)


###################################################
### code chunk number 68: lat.Rnw:2320-2337
###################################################
.iwidth <- 5
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
NY8$sar_trend <- nysar$fit$signal_trend
NY8$sar_stochastic <- nysar$fit$signal_stochastic
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
tr_at <- seq(-1, 1.3, length.out=21)
tr_rds <- rds(sum(tr_at >= 0)*2)[-(1:(sum(tr_at >= 0)-sum(tr_at < 0)))]
tr_pl <- spplot(NY8, c("sar_trend"), at=tr_at, col="transparent", col.regions=tr_rds, main=list(label="Trend", cex=0.8))
st_at <- seq(-0.16, 0.39, length.out=21)
st_rds <- rds(sum(st_at >= 0)*2)[-(1:(sum(st_at >= 0)-sum(st_at < 0)))]
st_pl <- spplot(NY8, c("sar_stochastic"), at=st_at, col="transparent", col.regions=st_rds, main=list(label="Stochastic", cex=0.8))
plot(tr_pl, split=c(1,1,2,1), more=TRUE)
plot(st_pl, split=c(2,1,2,1), more=FALSE)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 69: lat.Rnw:2359-2362
###################################################
nylmw <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, weights=POP8)
summary(nylmw)
NY8$lmwresid <- residuals(nylmw)


###################################################
### code chunk number 70: lat.Rnw:2378-2394
###################################################
.iwidth <- 5
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(lattice)
#trellis.par.set(canonical.theme(color = FALSE))
trellis.par.set(canonical.theme("postscript", color=TRUE))
library(RColorBrewer)
gry <- c(rev(brewer.pal(6, "Reds")[1:4]), colorRampPalette(brewer.pal(5, "Blues"))(9))
#gry <- brewer.pal(9, "Greys")[-1]
# RSB quietening greys
#gry <- grey.colors(9, 0.95, 0.55, 2.2)
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(NY8, c("lmresid", "lmwresid"), sp.layout=list(TCEpts), col.regions=gry, col="transparent", lwd=0.5, at=seq(-2,4.5,0.5))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 71: lat.Rnw:2403-2404
###################################################
lm.morantest(nylmw, NYlistw)


###################################################
### code chunk number 72: lat.Rnw:2423-2425
###################################################
nysarw<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME , data=NY8, listw=NYlistw, weights=POP8)
summary(nysarw)


###################################################
### code chunk number 73: lat.Rnw:2443-2455
###################################################
.iwidth <- 5
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
NY8$sarw_trend <- nysarw$fit$signal_trend
NY8$sarw_stochastic <- nysarw$fit$signal_stochastic
tr_pl <- spplot(NY8, c("sarw_trend"), at=tr_at, col="transparent", col.regions=tr_rds, main=list(label="Trend", cex=0.8))
st_pl <- spplot(NY8, c("sarw_stochastic"), at=st_at, col="transparent", col.regions=st_rds, main=list(label="Stochastic", cex=0.8))
plot(tr_pl, split=c(1,1,2,1), more=TRUE)
plot(st_pl, split=c(2,1,2,1), more=FALSE)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 74: lat.Rnw:2510-2513
###################################################
nycar<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME , data=NY8, family="CAR",
   listw=NYlistw)
summary(nycar)


###################################################
### code chunk number 75: lat.Rnw:2516-2517
###################################################
nylam1 <- c(nycar$lambda)


###################################################
### code chunk number 76: lat.Rnw:2548-2551
###################################################
nycarw<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="CAR",
   listw=NYlistw, weights=POP8)
summary(nycarw)


###################################################
### code chunk number 77: lat.Rnw:2624-2626
###################################################
nysarwM<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SAR",
   listw=NYlistw, weights=POP8, method="Matrix")


###################################################
### code chunk number 78: lat.Rnw:2628-2629
###################################################
summary(nysarwM)


###################################################
### code chunk number 79: lat.Rnw:2659-2664
###################################################
1/range(eigenw(NYlistw))
nysar_ll<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SAR",
   listw=NYlistw, llprof=100)
nysarw_ll<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SAR",
   listw=NYlistw, weights=POP8, llprof=100)


###################################################
### code chunk number 80: lat.Rnw:2670-2684
###################################################
.iwidth <- 6
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
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
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 81: lat.Rnw:2711-2714
###################################################
nysmaw<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SMA",
   listw=NYlistw, weights=POP8)
summary(nysmaw)


###################################################
### code chunk number 82: lat.Rnw:2752-2754
###################################################
library(lmtest)
bptest(nylm)


###################################################
### code chunk number 83: lat.Rnw:2776-2779
###################################################
library(sandwich)
coeftest(nylm)
coeftest(nylm, vcov=vcovHC(nylm, type="HC4"))


###################################################
### code chunk number 84: lat.Rnw:2813-2818
###################################################
NYlistwW <- nb2listw(NY_nb, style = "W")
res <- lm.LMtests(nylm, listw=NYlistwW, test="all")
tres <- t(sapply(res, function(x) c(x$statistic, x$parameter, x$p.value)))
colnames(tres) <- c("Statistic", "df", "p-value")
printCoefmat(tres)


###################################################
### code chunk number 85: lat.Rnw:2901-2904
###################################################
nylag <- lagsarlm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW)
summary(nylag)
bptest.sarlm(nylag)


###################################################
### code chunk number 86: lat.Rnw:2942-2943
###################################################
library(McSpatial)


###################################################
### code chunk number 87: lat.Rnw:2945-2947
###################################################
McRes <- sarml(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, wmat=listw2mat(NYlistwW), eigvar=eigenw(NYlistwW), print=FALSE, data=NY8)
c(McRes$beta, rho=McRes$rho, sig2=McRes$sig2)


###################################################
### code chunk number 88: lat.Rnw:2963-2966
###################################################
nymix <- lagsarlm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW, type="mixed")
nymix
anova(nymix, nylag)


###################################################
### code chunk number 89: lat.Rnw:3048-3051
###################################################
W <- as(as_dgRMatrix_listw(NYlistwW), "CsparseMatrix")
trMat <- trW(W, type="mult")
head(trMat)


###################################################
### code chunk number 90: lat.Rnw:3088-3094
###################################################
set.seed(987654)
imps <- impacts(nymix, tr=trMat, R=1999)
imps
HPDinterval(imps, choice="direct")
HPDinterval(imps, choice="indirect")
HPDinterval(imps, choice="total")


###################################################
### code chunk number 91: lat.Rnw:3130-3133
###################################################
nyerr <- errorsarlm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW)
summary(nyerr)
LR.sarlm(nyerr, nymix)


###################################################
### code chunk number 92: lat.Rnw:3153-3155
###################################################
nyerr1 <- errorsarlm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW, etype="emixed")
coef(nyerr1)


###################################################
### code chunk number 93: lat.Rnw:3186-3188
###################################################
set.seed(987654)
resMCMC <- MCMCsamp(nyerr1, mcmc=5000, burnin=500, listw=NYlistwW)


###################################################
### code chunk number 94: lat.Rnw:3197-3213
###################################################
.iwidth <- 6
.iheight <- 3
.ipointsize <- 8
.pwd <- 0.95
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mfrow=c(1, 2), mar=c(3, 3, 3, 1))
plot(density(resMCMC[,3]), ylab="", main="Direct", xlab="", lwd=2)
lines(density(imps$sres$direct[,1]), lty=2, lwd=2)
abline(v=0)
plot(density(resMCMC[,6]), ylab="", main="Indirect", xlab="", lwd=2)
lines(density(imps$sres$indirect[,1]), lty=2, lwd=2)
abline(v=0)
legend("topright", legend=c("Error Durbin", "Durbin"), lty=1:2, bty="n", lwd=c(2, 2), cex=0.7)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 95: lat.Rnw:3243-3244
###################################################
library(sphet)


###################################################
### code chunk number 96: lat.Rnw:3246-3248
###################################################
nyGMlag <- spreg(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=NY8, listw=NYlistwW, model="lag", het=FALSE)
summary(nyGMlag)


###################################################
### code chunk number 97: lat.Rnw:3264-3266
###################################################
nyGMerr <- spreg(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW, model="error", het=FALSE)
summary(nyGMerr)


###################################################
### code chunk number 98: lat.Rnw:3287-3288
###################################################
fit <- qregspiv(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, wmat=listw2mat(NYlistwW), data=NY8, tau=.5, nboot=200)


###################################################
### code chunk number 99: lat.Rnw:3347-3352
###################################################
library(mgcv)
NY8$x<-coordinates(NY8)[,1]/1000
NY8$y<-coordinates(NY8)[,2]/1000
nyGAM1 <- gam(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME+s(x,y), weights=POP8, data=NY8)
anova(nylmw, nyGAM1, test="Chisq")


###################################################
### code chunk number 100: lat.Rnw:3355-3356
###################################################
nylam1 <- c(summary(nyGAM1)$edf)


###################################################
### code chunk number 101: lat.Rnw:3396-3398
###################################################
owarn <- options("warn")
options("warn"=-1)


###################################################
### code chunk number 102: lat.Rnw:3401-3402
###################################################
nyGLMp <- glm(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8, family="poisson")


###################################################
### code chunk number 103: lat.Rnw:3404-3405 (eval = FALSE)
###################################################
## summary(nyGLMp)


###################################################
### code chunk number 104: lat.Rnw:3407-3412
###################################################
xx <- capture.output(print(summary(nyGLMp)))
zz <- strwrap(paste(xx[3:4], collapse=" "))
xx[3] <- zz[1]
xx[4] <- paste("   ", zz[2])
cat(xx, sep="\n")


###################################################
### code chunk number 105: lat.Rnw:3415-3416
###################################################
options("warn"=owarn$warn)


###################################################
### code chunk number 106: lat.Rnw:3421-3436
###################################################
.iwidth <- 5
.iheight <- 4.5
.ipointsize <- 10
.pwd <- 0.75
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-lat-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
library(RColorBrewer)
NY8$lmpresid <- residuals(nyGLMp, type="deviance")
gry <- c(rev(brewer.pal(6, "Reds")), brewer.pal(7, "Blues"))
#gry <- grey.colors(9, 0.95, 0.55, 2.2)
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(NY8, "lmpresid", sp.layout=list(TCEpts), col.regions=gry, col="transparent", lwd=0.5, at=seq(-3,3.5,0.5))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 107: lat.Rnw:3459-3461
###################################################
owarn <- options("warn")
options("warn"=-1)


###################################################
### code chunk number 108: lat.Rnw:3464-3466
###################################################
nyGAMp <- gam(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8))+s(x,y), data=NY8, family="poisson")
summary(nyGAMp)


###################################################
### code chunk number 109: lat.Rnw:3468-3469 (eval = FALSE)
###################################################
## anova(nyGLMp, nyGAMp, test="Chisq")


###################################################
### code chunk number 110: lat.Rnw:3471-3476
###################################################
xx <- capture.output(print(anova(nyGLMp, nyGAMp, test="Chisq")))
cat(xx[1:2], sep="\n")
cat(strwrap(xx[3], width=70, exdent=4), sep="\n")
cat(strwrap(paste(xx[4:5], collapse=" "), width=70, exdent=4), sep="\n")
cat(xx[6:8], sep="\n")


###################################################
### code chunk number 111: lat.Rnw:3479-3481
###################################################
nylam1 <- c(summary(nyGAMp)$edf)
options("warn"=owarn$warn)


###################################################
### code chunk number 112: lat.Rnw:3565-3573
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
cat("\n")
ver <- system("svnversion", intern=TRUE)
cat("%SVN version", ver, "\n")
cat("\n")
sT <- capture.output(print(Sys.time()))
cat("%", sT, "\n")


###################################################
### code chunk number 113: lat.Rnw:3576-3577
###################################################
options(op)


