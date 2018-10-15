### R code from vignette source 'dismap.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: dismap.Rnw:5-10
###################################################
if (!exists("book_R_dont_trash")) rm(list=ls())
#library(digest)
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
### code chunk number 3: dismap.Rnw:19-20
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 4: afig (eval = FALSE)
###################################################
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
## lattice.options(default.theme = col.whitebg())


###################################################
### code chunk number 5: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 6: dismap.Rnw:180-181
###################################################
options("width"=60)


###################################################
### code chunk number 7: dismap.Rnw:183-191
###################################################
library(maptools)
library(spdep)
nc_file <- system.file("shapes/sids.shp", package="maptools")[1]
llCRS <- CRS("+proj=longlat +datum=NAD27")
nc <- readShapePoly(nc_file, ID="FIPSNO", proj4string=llCRS)
rn <- sapply(slot(nc, "polygons"), function(x) slot(x, "ID"))
gal_file <- system.file("etc/weights/ncCR85.gal", package="spdep")[1]
ncCR85 <- read.gal(gal_file, region.id = rn)


###################################################
### code chunk number 8: dismap.Rnw:193-194
###################################################
options("width"=70)


###################################################
### code chunk number 9: dismap.Rnw:243-244
###################################################
nc$Observed <- nc$SID74


###################################################
### code chunk number 10: dismap.Rnw:252-255
###################################################
nc$Population <- nc$BIR74
r <- sum(nc$Observed)/sum(nc$Population)
nc$Expected <- nc$Population*r


###################################################
### code chunk number 11: dismap.Rnw:263-264
###################################################
nc$SMR <- nc$Observed/nc$Expected


###################################################
### code chunk number 12: dismap.Rnw:270-291
###################################################
.iwidth <- 6
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(RColorBrewer)

#Used method proposed by Nicky Best
logSMR <- log(nc$SMR[nc$SMR>0])
nsteps <- 7
step <- (max(logSMR)-min(logSMR))/nsteps
brks <- exp(min(logSMR)+(0:nsteps)*step)
brks[1] <- 0
cols <- c(rev(brewer.pal(3, "Blues")), brewer.pal(4, "Reds"))
grps <- as.ordered(cut(nc$SMR, brks, include.lowest=TRUE))
plot(nc, col=cols[unclass(grps)], axes = FALSE)
box()
degAxis(1)
degAxis(2, at=c(34,35,36,37)) 
legend("bottomleft",legend=levels(grps), fill=cols, bty="n",cex=0.8,y.intersp=0.8) 
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 13: dismap.Rnw:318-346
###################################################
.iwidth <- 6
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(epitools)
CISMR <- pois.exact(nc$Observed, nc$Expected)
plot(1,1, type="n", xlim=c(1,100), ylim=c(0,9),
  main= "Confidence intervals of the SMR",
  xlab="County", ylab="Relative Risk", xaxt="n")
abline(h=1, lty=2)

for(i in 1:100) {
        if(CISMR$lower[i]>1 ) {
            #sig.col <- 'red'
            sig.col <- brewer.pal(4, "Reds")[4]
            col <- sig.col
            lty <- 2
            text(i, CISMR$upper[i]+.31, nc$NAME[i],
                srt=90, col=sig.col, cex=.85)
        } else {
            col <- "black"
            lty <- 1
        }
        lines(c(i,i), c(CISMR$lower[i],CISMR$upper[i]), col=col, lty=lty)
        points(i, nc$SMR[i], pch=18, col=col)
}
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 14: dismap.Rnw:425-430
###################################################
library(DCluster)
eb <- empbaysmooth(nc$Observed, nc$Expected)
nc$EBPG <- eb$smthrr
eb$nu
eb$alpha


###################################################
### code chunk number 15: dismap.Rnw:432-434
###################################################
ebnu <- eb$nu
ebalpha <- eb$alpha


###################################################
### code chunk number 16: dismap.Rnw:463-503
###################################################
.iwidth <- 6
.iheight <- 5
.ipointsize <- 10
.pwd <- 0.8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
nc$pvalpois <- ppois(nc$Observed, nc$Expected, lower.tail=FALSE)

nbparam <- calculate.mle(as(nc, "data.frame"), model="negbin")
nc$pvalnegbin <- pnbinom(nc$Observed, size=nbparam$size, prob=nbparam$prob,
  lower.tail=FALSE)

colorkeypval <- list(labels=as.character(c(0, 0.01, 0.05, 0.1, .5, 1)), 
  at=(0:5)/5, height=.5)

pvalcols <- brewer.pal(5, "Reds")
#pvalcols <- brewer.pal(5, "Greys")
# RSB quieting greys
#pvalcols <- grey.colors(5, 0.95, 0.55, 2.2)

print(spplot(nc, c("pvalpois","pvalnegbin"), col.regions=rev(pvalcols), 
  at=c(0, 0.01, 0.05, 0.1, .5, 1), axes=TRUE))#, colorkey=colorkeypval ))

#-------->Plot names of the regions
#idx <- nc@data$pvalnegbin<.05
#loc <- nc@data[idx,c("x","y")]
#txt <- as.character(nc@data[idx,"NAME"])
#
#
#kk <- spplot(nc, c("pvalpois","pvalnegbin"), col.regions=rev(pvalcols), 
#  at=c(0, 0.01, 0.05, 0.1, .5, 1), axes=TRUE, colorkey=colorkeypval)
# 
#trellis.focus("panel", 1,1) 
#grid.text(text, x, y)
#
#print(kk)
#
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 17: dismap.Rnw:558-560
###################################################
ebln <- lognormalEB(nc$Observed, nc$Expected)
nc$EBLN <- exp(ebln$smthrr)


###################################################
### code chunk number 18: dismap.Rnw:600-603
###################################################
library(spdep)
EBMarshall <- EBest(nc$Observed, nc$Expected)
nc$EBMarshall <- EBMarshall[,2]


###################################################
### code chunk number 19: dismap.Rnw:610-631
###################################################
.iwidth <- 7
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
#Display all the risk estimates
atcol <- (0:5)*max(nc$SMR)/5
colorkey <- list(labels=as.character(c(formatC(brks, format="f", dig=2))),
  at=atcol,  height=.5)

#cols <- brewer.pal(5, "Oranges")
#cols <- brewer.pal(5, "Greys")
# RSB quieting greys
#cols <- grey.colors(5, 0.95, 0.45, 2.2)
#grps <- as.ordered(cut(SMR, brks, include.lowest=TRUE))

print(spplot(nc, c("SMR","EBPG", "EBLN", "EBMarshall"), col.regions=cols, 
  at=brks, axes = TRUE))#, colorkey=colorkey))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 20: dismap.Rnw:719-720
###################################################
nc$EBMrshloc <- EBlocal(nc$Observed, nc$Expected, ncCR85)$est


###################################################
### code chunk number 21: dismap.Rnw:728-739
###################################################
.iwidth <- 6
.iheight <- 5
.ipointsize <- 10
.pwd <- 0.8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
#par(mfrow=c(1,2))
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
print(spplot(nc, c("EBMarshall", "EBMrshloc"), col.regions=cols, at=brks))#, colorkey=colorkey))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 22: dismap.Rnw:784-795
###################################################
.iwidth <- 6
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
oopar <- par(mar=c(3,7,2,1)+0.1)
boxplot(as(nc, "data.frame")[,c("SMR", "EBPG", "EBLN", "EBMarshall", "EBMrshloc")], cex.lab=.5, las=1, horizontal=TRUE)
#boxplot(as(nc, "data.frame")[,c("SMR", "EBPG", "EBLN", "EBMarshall", "EBMrshloc")], cex.lab=.5, las=3)
#boxplot(nc@data[,c("SMR", "EBPG", "EBLN", "EBMarshall", "EBMrshloc")],las=2)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 23: dismap.Rnw:971-979
###################################################
#Data loaded from a previous Windows session
#This file contains the output of MCMC simulations of the P-G model
#I have used BRugs and OpenBUGS to get the results
#
load("PG.RData")
#con <- url("http://www.bias-project.org.uk/ASDAR/PG.RData", open="rb")
#load(con)
#close(con)


###################################################
### code chunk number 24: dismap.Rnw:994-997
###################################################
library(R2WinBUGS)
N <- length(nc$Observed)
d <- list(N=N, observed=nc$Observed, expected=nc$Expected)


###################################################
### code chunk number 25: dismap.Rnw:999-1003 (eval = FALSE)
###################################################
## pgmodelfile <- paste(getwd(), "/PG-model.txt", sep="")
## wdir <- paste(getwd(), "/PG", sep="")
## if(!file.exists(wdir)){dir.create(wdir)}
## BugsDir <- "/home/rsb/.wine/dosdevices/c:/Program Files/WinBUGS14"


###################################################
### code chunk number 26: dismap.Rnw:1011-1018 (eval = FALSE)
###################################################
## MCMCres <- bugs(data = d, inits = list(list(nu = 1, alpha = 1)), 
##   working.directory = wdir,
##   parameters.to.save = c("theta", "nu", "alpha"),
##   n.chains = 1, n.iter = 20000, n.burnin = 10000, n.thin = 10,
##   model.file = pgmodelfile,
##   bugs.directory = BugsDir,
##   WINEPATH = "/usr/bin/winepath")


###################################################
### code chunk number 27: dismap.Rnw:1051-1053
###################################################
nc$PGmean <- MCMCres$mean$theta
nc$PGmedian <- MCMCres$median$theta


###################################################
### code chunk number 28: dismap.Rnw:1112-1136
###################################################
.iwidth <- 6
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
plot(1,1, type="n", xlim=c(1,100), ylim=c(0,4),
  main= "Credible intervals of the relative risks",
  xlab="County", ylab="Relative Risk", xaxt="n")
abline(h=1, lty=2)

for(i in 1:100) {
	if(MCMCres$summary[i,3]>1 ) {
		col <- sig.col #gray(.4)
		lty <- 2
		text(i, MCMCres$summary[i,7]+.31, nc$NAME[i], 
       		srt=90, col=sig.col, cex=.85)
	} else {
		col <- "black"
        lty <- 1
	}
	lines(c(i,i), c(MCMCres$summary[i,3], MCMCres$summary[i,7]), col=col, lty=lty)
	points(i, MCMCres$median$theta[i], pch=18, col=col)
}
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 29: dismap.Rnw:1158-1172
###################################################
.iwidth <- 7
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))

print(spplot(nc, c("SMR","EBPG", "PGmean", "PGmedian"),
  col.regions=cols,  at=brks, axes = TRUE))#, colorkey=colorkey))

#cat(c("alpha =",eb$alpha, "\n"))
#cat(c("nu=", eb$nu, "\n"))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 30: dismap.Rnw:1219-1221
###################################################
library(R2BayesX)
nc$AREAID <- 1:nrow(nc)


###################################################
### code chunk number 31: dismap.Rnw:1223-1226 (eval = FALSE)
###################################################
## pgbayesx <- bayesx(Observed ~ sx(AREAID, bs="re"), 
##    offset= log(nc$Expected),
##    family = "poisson", data = as(nc, "data.frame"))


###################################################
### code chunk number 32: dismap.Rnw:1231-1232
###################################################
load("pgbayesx.RData")


###################################################
### code chunk number 33: dismap.Rnw:1246-1247
###################################################
library(INLA)


###################################################
### code chunk number 34: dismap.Rnw:1249-1250 (eval = FALSE)
###################################################
## inla.version()


###################################################
### code chunk number 35: dismap.Rnw:1252-1254
###################################################
obj <- capture.output(inla.version())
cat(obj[3], "\n")


###################################################
### code chunk number 36: dismap.Rnw:1256-1260
###################################################
pginla <- inla(Observed ~ offset(log(Expected)) -1 + f(AREAID, model = "iid"),
   family = "poisson",  data = as(nc, "data.frame"),
   control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE)
)


###################################################
### code chunk number 37: dismap.Rnw:1276-1277
###################################################
library(CARBayes)


###################################################
### code chunk number 38: dismap.Rnw:1279-1284 (eval = FALSE)
###################################################
## ncdf <- as(nc, "data.frame")
## attach(ncdf)
## pgcarbayes <- poisson.independent(formula = Observed ~ offset(log(Expected)),
##    burnin=5000, n.sample=10000)
## detach(ncdf)


###################################################
### code chunk number 39: dismap.Rnw:1288-1289
###################################################
load("pgcarbayes.RData")


###################################################
### code chunk number 40: dismap.Rnw:1298-1301
###################################################
nc$PGBAYESX <- pgbayesx$fitted.values[order(pgbayesx$bayesx.setup$order),2]/nc$Expected
nc$PGINLA <- pginla$summary.fitted.values$mean/nc$Expected
nc$PGCARBAYES <- pgcarbayes$fitted.values[,1]/nc$Expected


###################################################
### code chunk number 41: dismap.Rnw:1313-1322
###################################################
.iwidth <- 7
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
trellis.par.set(canonical.theme("postscript", color=TRUE))
print(spplot(nc, c("PGmean", "PGBAYESX", "PGINLA", "PGCARBAYES"),
 col.regions=cols,  at=brks, axes = TRUE))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 42: dismap.Rnw:1456-1470
###################################################
.iwidth <- 6
.iheight <- 3.5
.ipointsize <- 10
.pwd <- 0.8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
nc$RATIO <- nc$NWBIR74/nc$BIR74
#print(spplot(nc, "RATIO", col.regions=brewer.pal(9,"Blues"), at=.8*0:9/9))
#print(spplot(nc, "RATIO", col.regions=brewer.pal(9,"Greys"), at=.8*0:9/9))
# RSB quietening greys
print(spplot(nc, "RATIO", col.regions=brewer.pal(4, "Reds"), at=.8*0:4/4))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 43: dismap.Rnw:1485-1492 (eval = FALSE)
###################################################
## idx <- match(attr(ncCR85.nb, "region.id"), nc$"CNTY_ID")
## nc.nb <- ncCR85
## nc.nb <- nc.nb[order(idx)]
## nc.nb <- lapply(nc.nb, function(X, idx){idx[X]}, idx=(idx))
## class(nc.nb) <- "nb"
## nc.nb <- nc.nb[(order(idx))]
## nc.nb <- nb2WB(nc.nb)


###################################################
### code chunk number 44: dismap.Rnw:1516-1517
###################################################
nc.nb <- nb2WB(ncCR85)


###################################################
### code chunk number 45: dismap.Rnw:1527-1528
###################################################
options("width"=60)


###################################################
### code chunk number 46: dismap.Rnw:1530-1541
###################################################
nc$nwprop <- nc$NWBIR74/nc$BIR74

d <- list(N=N, observed=nc$Observed, expected=nc$Expected,
  nonwhite=nc$nwprop,
  adj=nc.nb$adj,  weights=nc.nb$weights, num=nc.nb$num)

dwoutcov <- list(N=N, observed=nc$Observed, expected=nc$Expected,
  adj=nc.nb$adj,  weights=nc.nb$weights, num=nc.nb$num)

inits <- list(u=rep(0,N), v=rep(0,N), alpha=0, beta=0, precu=.001, precv=.001)
#inits$v[d$num==0] <- NA


###################################################
### code chunk number 47: dismap.Rnw:1543-1544
###################################################
options("width"=70)


###################################################
### code chunk number 48: dismap.Rnw:1565-1571 (eval = FALSE)
###################################################
## 
## 
## bymmodelfile <- paste(getwd(), "/BYM-model.txt", sep="")
## wdir <- paste(getwd(), "/BYM", sep="")
## if(!file.exists(wdir)){dir.create(wdir)}
## BugsDir <- "/home/rsb/.wine/dosdevices/c:/Program Files/WinBUGS14"


###################################################
### code chunk number 49: dismap.Rnw:1573-1574
###################################################
options("width"=60)


###################################################
### code chunk number 50: dismap.Rnw:1576-1583 (eval = FALSE)
###################################################
## MCMCres<- bugs(data = d, inits = list(inits),
##   working.directory = wdir,
##   parameters.to.save = c("theta", "alpha", "beta", "u", "v", "sigmau", "sigmav"),
##   n.chains = 1, n.iter = 30000, n.burnin = 20000, n.thin = 10,
##   model.file = bymmodelfile,
##   bugs.directory = BugsDir,
##   WINEPATH = "/usr/bin/winepath")


###################################################
### code chunk number 51: dismap.Rnw:1585-1586
###################################################
options("width"=70)


###################################################
### code chunk number 52: dismap.Rnw:1598-1602
###################################################
load("BYM.RData")
#con <- url("http://www.bias-project.org.uk/ASDAR/BYM.RData", open="rb")
#load(con)
#close(con)


###################################################
### code chunk number 53: dismap.Rnw:1609-1612
###################################################
nc$BYMmean <- MCMCres$mean$theta
nc$BYMumean <- MCMCres$mean$u
nc$BYMvmean <- MCMCres$mean$v


###################################################
### code chunk number 54: dismap.Rnw:1648-1660
###################################################
.iwidth <- 5
.iheight <- 6.5
.ipointsize <- 8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(coda)
#ncoutput <- read.coda("BYM/coda1.txt", "BYM/codaIndex.txt")

plot(ncoutput[,c("deviance", "alpha", "beta", "theta[94]")])

     #save(file="BYM.RData", list=c("d", "inits", "MCMCres", "ncoutput") )

dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 55: dismap.Rnw:1680-1681
###################################################
geweke.diag(ncoutput[,c("deviance", "alpha", "beta", "theta[94]")])


###################################################
### code chunk number 56: dismap.Rnw:1694-1705
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.pwd <- 0.8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
print(spplot(nc, c("SMR", "BYMmean"), at=brks, col.regions=cols,
  axes=TRUE))#, colorkey=colorkey))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 57: dismap.Rnw:1716-1742
###################################################
.iwidth <- 6
.iheight <- 5
.ipointsize <- 10
.pwd <- 0.8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
ubrks <- quantile(nc$BYMumean, seq(0,1,1/5))
ubrks[6] <- ubrks[6]*1.01
#colorkey$labels <- as.character(formatC(ubrks,3))
colorkey$labels <- as.character(trunc(1000*ubrks)/1000)
colorkey$at <- ubrks[1]+ubrks[6]*0:5/5
#RSB
ubrks <- seq(-0.35, 0.35, 0.05)
print(spplot(nc, c("BYMumean", "BYMvmean"), at= ubrks, axes=TRUE, 
col.regions=c(rev(brewer.pal(7, "Reds")), brewer.pal(7, "Blues"))))#, colorkey=colorkey))
# RSB quietening greys
#col.regions=brewer.pal(5, "Greys"), colorkey=colorkey))
#col.regions=brewer.pal(5, "BuGn"), colorkey=colorkey))
#print(spplot(nc, c("BYMumean", "BYMumedian"), axes=TRUE))

#Display a prettier map
#  at=c(-260,-215,-210,-200,-150,-210,-100,-50,0,50), 
#  col.regions=rev(brewer.pal(9, "Purples"))) )
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 58: dismap.Rnw:1787-1820
###################################################
.iwidth <- 6
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
plot(1,1, type="n", xlim=c(1,100), ylim=c(0,4.5),
  main= "Credible intervals of the relative risks",
  xlab="County", ylab="Relative Risk", xaxt="n")
abline(h=1, lty=2)

for(i in 1:100)
{
        if(MCMCres$summary[i,3]>1 )
        {
                #col <- gray(.4)
                col <- sig.col
		lty <- 2
                text(i, MCMCres$summary[i,7]+.31, nc$NAME[i],
                srt=90, col=sig.col, cex=.85)
        }
        else
        {
                col <- "black"
		lty <- 1
        }

        lines(c(i,i), c(MCMCres$summary[i,3], MCMCres$summary[i,7]), col=col, lty=lty)
#       points(i, MCMCres$mean[offset+i], pch=20, col=col)
#       legend(20,4, legend=c("Median", "Mean"), pch=c(18,20))
}

points(1:100, MCMCres$median$theta, pch=18, col=col)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 59: dismap.Rnw:1874-1879
###################################################
INLA_BYM <- inla(Observed ~ nwprop + f(FIPS, model="iid") + 
  f(AREAID, model="besag", graph=nb2mat(ncCR85, style="B")) + 
  offset(log(Expected)), 
  family="poisson", data=as(nc, "data.frame"), 
  control.predictor=list(compute=TRUE))


###################################################
### code chunk number 60: dismap.Rnw:1894-1895 (eval = FALSE)
###################################################
## ncgra <- nb2gra(ncCR85)


###################################################
### code chunk number 61: dismap.Rnw:1912-1916 (eval = FALSE)
###################################################
## bymbayesx <- bayesx(Observed~nwprop+sx(AREAID, bs="re") +
##    sx(FIPSNO,bs="spatial", map=ncgra),
##    offset= log(nc$Expected),
##    family="poisson", data=as(nc, "data.frame"))


###################################################
### code chunk number 62: dismap.Rnw:1921-1922
###################################################
load("bymbayesx.RData")


###################################################
### code chunk number 63: dismap.Rnw:1934-1940 (eval = FALSE)
###################################################
## ncdf <- as(nc, "data.frame")
## attach(ncdf)
## obj <- poisson.bymCAR(Observed ~ nwprop + offset(log(Expected)), 
##    W=nb2mat(ncCR85, style="B"), 
##    n.sample=30000, burnin=20000, thin=10)
## detach(ncdf)


###################################################
### code chunk number 64: dismap.Rnw:1944-1945
###################################################
load("obj_CARBayes.RData")


###################################################
### code chunk number 65: dismap.Rnw:1954-1957
###################################################
nc$BAYESX <- bymbayesx$fitted.values[order(bymbayesx$bayesx.setup$order),2]/nc$Expected
nc$INLA <- INLA_BYM$summary.fitted.values[,1]/nc$Expected
nc$CARBayes <- obj$fitted.values[,1]/nc$Expected


###################################################
### code chunk number 66: dismap.Rnw:1972-1982
###################################################
.iwidth <- 8
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
print(spplot(nc, c("BYMmean", "BAYESX", "INLA", "CARBayes"), 
   at=brks, axes=TRUE, col.regions=cols))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 67: dismap.Rnw:2001-2045
###################################################
.iwidth <- 7
.iheight <- 6
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
oldpar <- par(mfrow=c(2,2))
#Intercept
d1 <- density(attr(bymbayesx$fixed.effects, "sample")[,1])
plot(as.data.frame(d1[1:2]), main="Intercept", type="l")
lines(INLA_BYM$marginals.fixed[[1]], lty=2)
lines(density( obj$samples.beta[,1], bw=d1$bw), lty=3)
legend("topleft", legend=c("BayesX", "INLA", "CARBayes"), bty="n", 
   lty=1:3, cex=.65)

#Nwprop
d2 <- density(attr(bymbayesx$fixed.effects, "sample")[,2])
plot(as.data.frame(d2[1:2]), main="nwprop", type="l")
lines(INLA_BYM$marginals.fixed[[2]], lty=2)
lines(density( obj$samples.beta[,2], bw=d2$bw), lty=3)
legend("topleft", legend=c("BayesX", "INLA", "CARBayes"), bty="n",
   lty=1:3, cex=.65)

#FIXME: Add variances here
#Density of posterio means of non-spatial random effects
d3 <- density(bymbayesx$effects[["sx(AREAID)"]]$Mean)
plot(as.data.frame(d3[1:2]), main="Non-spatial r.eff.", type="l")
lines(density(INLA_BYM$summary.random$FIPS$mean, bw=d3$bw), lty=2)
lines(density(apply( obj$samples.phi,2,mean), bw=d3$bw), lty=3)
lines(density(MCMCres$mean$u, bw=d3$bw), lty=4)
legend("topleft", legend=c("BayesX", "INLA", "CARBayes", "WinBUGS"),  bty="n",
   lty=1:4, cex=.65)


#Variance of spatial random effects
#FIXME: check the prior that INLA uses for this. Seems too smoothed
d4 <- density(bymbayesx$effects[["sx(FIPSNO)"]]$Mean)
plot(as.data.frame(d4[1:2]), main="Spatial r.eff.", type="l")
lines(density(INLA_BYM$summary.random$AREAID$mean, bw=d4$bw), lty=2)
lines(density(apply( obj$samples.theta,2,mean), bw=d4$bw), lty=3)
lines(density(MCMCres$mean$v, bw=d4$bw), lty=4)
legend("topleft", legend=c("BayesX", "INLA", "CARBayes", "WinBUGS"), bty="n",
   lty=1:4, cex=.65)

par(oldpar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")


###################################################
### code chunk number 68: dismap.Rnw:2075-2078 (eval = FALSE)
###################################################
## bayesxps <- bayesx(Observed~sx(nwprop, bs="ps", knots=10),
##    offset= log(nc$Expected),
##    family="poisson", data=as(nc, "data.frame"))


###################################################
### code chunk number 69: dismap.Rnw:2090-2095 (eval = FALSE)
###################################################
## nc$long <- coordinates(nc)[,1]
## nc$lat <- coordinates(nc)[,2]
## bayesxte <- bayesx(Observed~sx(long,lat, bs="te"),
##    offset= log(nc$Expected),
##    family="poisson", data=as(nc, "data.frame"))


###################################################
### code chunk number 70: dismap.Rnw:2101-2102
###################################################
load("bayesxpste.RData")


###################################################
### code chunk number 71: dismap.Rnw:2114-2121
###################################################
.iwidth <- 7
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
plot(bayesxps)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 72: dismap.Rnw:2131-2142
###################################################
.iwidth <- 7
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
#library(RColorBrewer)
#nc$bayesxte_Mean <- bayesxte$effects[[1]]$Mean[bayesxte$bayesx.setup$order]
#print(spplot(nc, "bayesxte_Mean", at=seq(-0.92, 0.92, 0.1), col.regions=c(rev(brewer.pal(9, "Reds")), brewer.pal(9, "Blues"))))
plot(bayesxte, image=TRUE)
#plot(bayesxte)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 73: dismap.Rnw:2236-2241 (eval = FALSE)
###################################################
## #
## #I do not know why but if these objects are not removed I 
## #get 10 pages of warnings!!
## #
## rm(list=c("EBMarshall","Observed","Expected","SMR"))


###################################################
### code chunk number 74: dismap.Rnw:2245-2246
###################################################
set.seed(1)


###################################################
### code chunk number 75: dismap.Rnw:2251-2254
###################################################
chtest <- achisq.test(Observed~offset(log(Expected)), 
   as(nc, "data.frame"), "multinom", 999)
chtest


###################################################
### code chunk number 76: dismap.Rnw:2269-2270
###################################################
1- pchisq(chtest$t0, 100-1)


###################################################
### code chunk number 77: dismap.Rnw:2307-2308
###################################################
set.seed(1)


###################################################
### code chunk number 78: dismap.Rnw:2312-2314
###################################################
pwtest <- pottwhitt.test(Observed~offset(log(Expected)), 
   as(nc, "data.frame"), "multinom", 999)


###################################################
### code chunk number 79: dismap.Rnw:2323-2325
###################################################
Oplus<- sum(nc$Observed)
1- pnorm(pwtest$t0, Oplus*(Oplus-1), sqrt(2*100*Oplus*(Oplus-1)))


###################################################
### code chunk number 80: dismap.Rnw:2395-2398
###################################################
col.W <- nb2listw(ncCR85, zero.policy=TRUE)
moranI.test(Observed~offset(log(Expected)), as(nc, "data.frame"), 
   "negbin", 999, listw=col.W, n=length(ncCR85), S0=Szero(col.W))


###################################################
### code chunk number 81: dismap.Rnw:2450-2461
###################################################
data(nc.sids)
idx <- match(nc$NAME, rownames(nc.sids))
nc$x <- nc.sids$x[idx]
nc$y <- nc.sids$y[idx]
coords <- cbind(nc$x, nc$y)
dlist <- dnearneigh(coords, 0, Inf)
dlist <- include.self(dlist)
dlist.d <- nbdists(dlist, coords)
phi <- 100
col.W.tango <- nb2listw(dlist, glist=lapply(dlist.d, 
  function(x, phi) {exp(-x/phi)}, phi=phi), style="C")


###################################################
### code chunk number 82: dismap.Rnw:2475-2476
###################################################
set.seed(1)


###################################################
### code chunk number 83: dismap.Rnw:2478-2480
###################################################
tango.test(Observed~offset(log(Expected)), as(nc, "data.frame"), "negbin", 999, 
   listw=col.W.tango, zero.policy=TRUE)


###################################################
### code chunk number 84: dismap.Rnw:2533-2536
###################################################
sidsgam <- opgam(data=as(nc, "data.frame"),  radius=30, step=10, alpha=.002)
gampoints <- SpatialPoints(sidsgam[,c("x", "y")]*1000, 
   CRS("+proj=utm +zone=18 +datum=NAD27"))


###################################################
### code chunk number 85: dismap.Rnw:2539-2543
###################################################
library(rgdal)
ll <- CRS("+proj=longlat +datum=NAD27")
gampoints <- spTransform(gampoints, ll)
gam.layout <- list("sp.points", gampoints)


###################################################
### code chunk number 86: dismap.Rnw:2627-2628
###################################################
set.seed(1234)


###################################################
### code chunk number 87: dismap.Rnw:2636-2637
###################################################
options("width"=60)


###################################################
### code chunk number 88: dismap.Rnw:2639-2643
###################################################
mle <- calculate.mle(as(nc, "data.frame"), model="negbin")
thegrid <- as(nc, "data.frame")[,c("x","y")]
knresults <- opgam(data=as(nc, "data.frame"), thegrid=thegrid, alpha=.05,
   iscluster=kn.iscluster, fractpop=0.15, R=99, model="negbin", mle=mle)


###################################################
### code chunk number 89: dismap.Rnw:2645-2646
###################################################
options("width"=70)


###################################################
### code chunk number 90: dismap.Rnw:2653-2682
###################################################
.iwidth <- 6
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
#Plot all centroids and significant ones in red
#print(knresults)
#plot(nc@data$x, nc@data$y, main="Kulldorff and Nagarwalla's method",
#  xlab="Easting", ylab="Northing")

clusters <- get.knclusters(as(nc, "data.frame"), knresults)
i <- which.max(knresults$statistic)

nc$KNcluster <- "county"
nc$KNcluster[clusters[[i]]] <- "cluster"
nc$KNcluster[clusters[[i]][1]] <- "centre"
nc$KNcluster <- as.factor(nc$KNcluster)
bp0 <- brewer.pal(4, "Reds")

print(spplot(nc, "KNcluster",# main="Kulldorff's method",
  col.regions=c(bp0[4], bp0[3], bp0[1]), sp.layout=list("sp.points",
  gampoints, col=gray(.4), pch=4)))
# xlab="Easting", ylab="Northing", col.regions=c("white", "#E34A33","#FDBB84")))

#c("white", "red", "magenta") ))

dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 91: dismap.Rnw:2690-2691
###################################################
save(file="dismap.RDat", list=ls())


###################################################
### code chunk number 92: dismap.Rnw:2764-2768
###################################################
stone.stat(as(nc, "data.frame"), region=which(nc$NAME=="Anson"))
st <- stone.test(Observed~offset(log(Expected)), as(nc, "data.frame"), 
   model="negbin", 99, region=which(nc$NAME=="Anson"))
st


###################################################
### code chunk number 93: dismap.Rnw:2846-2861
###################################################
#Temporal trend
library(spacetime)
library(xts)
#Load STFDF
load("../Data/brainNM.RData")
#nmf <- slot(brainst, "sp")

obs <- xts(brainst@data$Observed,
   as.Date(as.character(brainst@data$Year), format="%Y") )
expt <- xts(brainst@data$Expected,
   as.Date(as.character(brainst@data$Year), format="%Y") )

obsy <- apply.yearly(obs, sum)
expty <- apply.yearly(expt, sum)
smry <- obsy/expty


###################################################
### code chunk number 94: dismap.Rnw:2866-2874
###################################################
.iwidth <- 7
.iheight <- 5.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
print(stplot(brainst[,,"SMR"], 1973:1991, at=c(0,.75, .9, 1, 1.1,1.25, 2, 3, 8),
   col.regions=brewer.pal(8, "Reds")) )
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 95: dismap.Rnw:2885-2892
###################################################
.iwidth <- 6
.iheight <- 2
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
plot(smry, main="Standardised Mortality Ratio by Year")
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 96: dismap.Rnw:2938-2940
###################################################
library(spdep)
neib <- poly2nb(nmf, row.names=1:length(nmf))


###################################################
### code chunk number 97: dismap.Rnw:2960-2964
###################################################
library(INLA)
hyper1 <- list(prec=list(param=c(.001, .001)))
form <- Observed~1+IDLANLre + f(Year, model="rw1", hyper = list(prec=list(param=c(.001,0.001)))) + f(ID, model="besag", graph=nb2mat(neib), hyper = hyper1)
inlares <- inla(form, family="poisson", data=slot(brainst, "data"), E=Expected, control.predictor=list(compute=TRUE), control.results=list(return.marginals.predictor=TRUE))


###################################################
### code chunk number 98: dismap.Rnw:2983-2985 (eval = FALSE)
###################################################
## nmgra <- nb2gra(neib)
## nmbayesx <- bayesx(Observed~IDLANLre+sx(Year, bs="rw1")+sx(ID, bs="spatial", map=nmgra),offset= log(brainst$Expected), family="poisson", data=as(brainst, "data.frame"))


###################################################
### code chunk number 99: dismap.Rnw:2995-2996
###################################################
load("nmbayesx.RData")


###################################################
### code chunk number 100: dismap.Rnw:3014-3030
###################################################
.iwidth <- 6
.iheight <- 2.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
#Brain cancer in New Mexico: Distance to LANL
oopar <- par(mfrow=c(1,2), mar=c(3, 4, 1, 1))
#INLA
labels.fix = names(inlares$marginals.fixed)
for (i in 1:2) {
 plot(inla.smarginal(inlares$marginals.fixed[[i]]), type="l", xlab = "", ylab = inla.nameunfix(labels.fix[i]))
 lines(density(attr(nmbayesx$fixed.effects, "sample")[,i] ), lty=2)
}
legend("topleft", legend=c("INLA", "BayesX"), lty=1:2, bty="n", cex=0.8)
par(oopar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 101: dismap.Rnw:3044-3046
###################################################
nmf$SPINLA <- inlares$summary.random$ID$mean
nmf$SPBAYESX <- nmbayesx$effects[["sx(ID)"]]$Mean


###################################################
### code chunk number 102: dismap.Rnw:3068-3089
###################################################
.iwidth <- 7
.iheight <- 3.5
.ipointsize <- 8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
oldpar<-par(mfrow=c(1,2))
#Brain cancer in New Mexico: Temporal trend
plot(1973:1991, inlares$summary.random$Year[,2], type="l", lwd=1.5,
main="Temporal trend (INLA)", ylim=c(-.25, .25), xlab="Year", ylab="")
lines(1973:1991, inlares$summary.random$Year[,4], lty=2)
lines(1973:1991, inlares$summary.random$Year[,6], lty=2)
abline(h=0, lty=3)

plot(1973:1991, nmbayesx$effects[["sx(Year)"]]$Mean, type="l", lwd=1.5,
main="Temporal trend (BayesX))", ylim=c(-.25, .25), xlab="Year", ylab="")
lines(1973:1991, nmbayesx$effects[["sx(Year)"]]$"2.5%", lty=2)
lines(1973:1991, nmbayesx$effects[["sx(Year)"]]$"97.5%", lty=2)
abline(h=0, lty=3)

par(oldpar)
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 103: dismap.Rnw:3099-3113
###################################################
.iwidth <- 7
.iheight <- 3.5
.ipointsize <- 8
.epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
lattice.options(default.theme = col.whitebg())
library(lattice)
trellis.par.set(canonical.theme("postscript", color=TRUE))
spl <- list(list("sp.points", losalamos, pch=19, col="red"),
   list("sp.text", coordinates(losalamos), "Los Alamos National Laboratory", pos=1, col="black", cex=0.7))
cols <- colorRampPalette(c(rev(brewer.pal(4, "Reds")), brewer.pal(4, "Blues")))
ats <- quantile(nmf$SPINLA, seq(0, 1, 1/9))
print(spplot(nmf, c("SPINLA", "SPBAYESX"), at=ats, col.regions=cols(length(ats)-1), 
   col="grey", sp.layout=spl, main="Spatial effects"))
dev.null <- dev.off()
system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12
.pwd <- 0.95


###################################################
### code chunk number 104: dismap.Rnw:3128-3138 (eval = FALSE)
###################################################
## #Brain cancer in New Mexico: Spatial effects (no covariates)
## form2 <- Observed~1+f(Year, model="rw1")+f(ID, model="besag", graph=nb2mat(neib))
## 
## inlares2 <- inla(form2, family="poisson", data=slot(brainst, "data"),
##    E=Expected,
## control.predictor=list(compute=TRUE),
## #   quantiles=qnts,
##    control.results=list(return.marginals.predictor=TRUE)
## )
## nmf$CAR2 <- inlares2$summary.random$ID[,2]


###################################################
### code chunk number 105: dismap.Rnw:3145-3162 (eval = FALSE)
###################################################
## .iwidth <- 7
## .iheight <- 4.5
## .ipointsize <- 8
## .epsNo <- .epsNo + 1; file <- paste("../Art/Fig-dismap-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
## lattice.options(default.theme = col.whitebg())
## library(lattice)
## trellis.par.set(canonical.theme("postscript", color=TRUE))
## spl <- list(list("sp.points", losalamos, pch=19, col="red"),
##    list("sp.text", coordinates(losalamos), "Los Alamos National Laboratory", pos=1, col="black", cex=0.8))
## cols <- colorRampPalette(c(rev(brewer.pal(4, "Reds")), brewer.pal(4, "Blues")))
## ats <- quantile(nmf$CAR, seq(0, 1, 1/9))
## CARplot <- spplot(nmf, c("CAR"), at=ats, col.regions=cols(length(ats)-1), col="grey", sp.layout=spl, main="Covariates model")
## ats <- quantile(nmf$CAR2, seq(0, 1, 1/9))
## CAR2plot <- spplot(nmf, c("CAR2"), at=ats, col.regions=cols(length(ats)-1), col="grey", sp.layout=spl, main="Null model")
## plot(CARplot, split=c(1,1,2,1), more=TRUE)
## plot(CAR2plot, split=c(2,1,2,1), more=FALSE)
## dev.null <- dev.off()
## system(paste("../scripts/mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=", .pwd, "\\textwidth]{", substr(file, 4, nchar(file)), "}\n\n", sep="")
## .iwidth <- 5
## .iheight <- 6
## .ipointsize <- 12
## .pwd <- 0.95


###################################################
### code chunk number 106: dismap.Rnw:3203-3211
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
### code chunk number 107: dismap.Rnw:3214-3215
###################################################
options(op)


