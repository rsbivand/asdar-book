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
## .epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)


###################################################
### chunk number 5: zfig eval=FALSE
###################################################
## dev.null <- dev.off()
## system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 6: 
###################################################
options("width"=60)


###################################################
### chunk number 7: 
###################################################
library(maptools)
library(spdep)

#Load data
nc_file <- system.file("shapes/sids.shp", package="maptools")[1]
llCRS <- CRS("+proj=longlat +datum=NAD27")
nc <- readShapePoly(nc_file, ID="FIPSNO", proj4string=llCRS)
#example(nc.sids, "spdep")
#nc<-nc.sids
#nc <- readShapePoly(system.file("etc/shapes/sids.shp", package="spdep")[1],
#   ID="FIPSNO", 
#proj4string=CRS("+proj=longlat +datum=NAD27"))
rn <- sapply(slot(nc, "polygons"), function(x) slot(x, "ID"))
gal_file <- system.file("etc/weights/ncCR85.gal", package="spdep")[1]
ncCR85 <- read.gal(gal_file, region.id=rn)


###################################################
### chunk number 8: 
###################################################
options("width"=70)


###################################################
### chunk number 9: 
###################################################
#names(nc)#Variables in the dataset

nc$Observed<-nc$SID74
nc$Population<-nc$BIR74#Population at risk; number of births
r<-sum(nc$Observed)/sum(nc$Population)
nc$Expected<-nc$Population*r

#Computed Standardised Mortality Ratio
nc$SMR<-nc$Observed/nc$Expected


###################################################
### chunk number 10: 
###################################################
.iwidth <- 6
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
#Display SMR 
library(RColorBrewer)

#brks<-quantile(SMR, seq(0,1,1/5))
#Used method proposed by Nicky Best
logSMR<-log(nc$SMR[nc$SMR>0])
nsteps <- 5
step<-(max(logSMR)-min(logSMR))/nsteps
brks<-exp(min(logSMR)+(0:nsteps)*step)
brks[1]<-0
#print(brks)
#cols <- brewer.pal(5, "Oranges")
#cols <- brewer.pal(5, "Greys")
# RSB quietening greys
cols <- grey.colors(nsteps, 0.95, 0.55, 2.2)
grps<-as.ordered(cut(nc$SMR, brks, include.lowest=TRUE))
plot(nc, col=cols[unclass(grps)], axes = FALSE)
box()
degAxis(1)
degAxis(2, at=c(34,35,36,37)) 
legend("bottomleft",legend=levels(grps), fill=cols, bty="n",cex=0.8,y.intersp=0.8) 
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 11: 
###################################################
.iwidth <- 6
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(epitools)

CISMR<-pois.exact(nc$Observed, nc$Expected)

plot(1,1, type="n", xlim=c(1,100), ylim=c(0,9),
  main= "Confidence intervals of the SMR",
  xlab="County", ylab="Relative Risk", xaxt="n")
abline(h=1, lty=2)

for(i in 1:100)
{
        if(CISMR$lower[i]>1 )
        {
                col<-gray(.4)
		lty<-2
                text(i, CISMR$upper[i]+.31, nc$NAME[i],
                srt=90, col=gray(.4), cex=.85)
        }
        else
        {
                col<-"black"
		lty<-1
        }

        lines(c(i,i), c(CISMR$lower[i],CISMR$upper[i]), col=col, lty=lty)
        points(i, nc$SMR[i], pch=18, col=col)
#	legend(10, 8, legend="SMR", pch=18, bty="n")
}
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 12: 
###################################################
library(DCluster)

eb<-empbaysmooth(nc$Observed, nc$Expected)

#summary(SMR)
#summary(eb$smthrr)
nc$EBPG<-eb$smthrr


###################################################
### chunk number 13: 
###################################################
ebnu <- eb$nu
ebalpha <- eb$alpha


###################################################
### chunk number 14: 
###################################################
.iwidth <- 6
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(lattice)
trellis.par.set(canonical.theme(color = FALSE))
nc$pvalpois<-ppois(nc$Observed, nc$Expected, lower.tail=FALSE)

nbparam<-calculate.mle(as(nc, "data.frame"), model="negbin")
nc$pvalnegbin<-pnbinom(nc$Observed, size=nbparam$size, prob=nbparam$prob,
  lower.tail=FALSE)

colorkeypval<-list(labels=as.character(c(0, 0.01, 0.05, 0.1, .5, 1)), 
  at=(0:5)/5, height=.5)

#pvalcols<-brewer.pal(5, "Blues")
#pvalcols<-brewer.pal(5, "Greys")
# RSB quieting greys
pvalcols <- grey.colors(5, 0.95, 0.55, 2.2)

print(spplot(nc, c("pvalpois","pvalnegbin"), col.regions=rev(pvalcols), 
  at=c(0, 0.01, 0.05, 0.1, .5, 1), axes=TRUE, colorkey=colorkeypval ))

#-------->Plot names of the regions
#idx<-nc@data$pvalnegbin<.05
#loc<-nc@data[idx,c("x","y")]
#txt<-as.character(nc@data[idx,"NAME"])
#
#
#kk<-spplot(nc, c("pvalpois","pvalnegbin"), col.regions=rev(pvalcols), 
#  at=c(0, 0.01, 0.05, 0.1, .5, 1), axes=TRUE, colorkey=colorkeypval)
# 
#trellis.focus("panel", 1,1) 
#grid.text(text, x, y)
#
#print(kk)
#
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 15: 
###################################################
#Log-normal model
ebln<-lognormalEB(nc$Observed, nc$Expected)
nc$EBLN<-exp(ebln$smthrr)


###################################################
### chunk number 16: 
###################################################
library(spdep)

#Compute Marshall risk estimator
EBMarshall<-EBest(nc$Observed, nc$Expected)
nc$EBMarshall<-EBMarshall[,2]


###################################################
### chunk number 17: 
###################################################
.iwidth <- 7
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
#Display all the risk estimates
atcol<-(0:5)*max(nc$SMR)/5
colorkey<-list(labels=as.character(c(formatC(brks, format="f", dig=2))),
  at=atcol,  height=.5)

#cols <- brewer.pal(5, "Oranges")
#cols <- brewer.pal(5, "Greys")
# RSB quieting greys
cols <- grey.colors(5, 0.95, 0.45, 2.2)
#grps<-as.ordered(cut(SMR, brks, include.lowest=TRUE))

print(spplot(nc, c("SMR","EBPG", "EBLN", "EBMarshall"), col.regions=cols, 
  at=brks, axes = TRUE, colorkey=colorkey))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 18: 
###################################################
nc$EBMrshloc<-EBlocal(nc$Observed, nc$Expected, ncCR85)$est


###################################################
### chunk number 19: 
###################################################
.iwidth <- 6
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
par(mfrow=c(1,2))
print(spplot(nc, c("EBMarshall", "EBMrshloc"), col.regions=cols, at=brks, colorkey=colorkey))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 20: 
###################################################
.iwidth <- 6
.iheight <- 3
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
oopar <- par(mar=c(3,7,2,1)+0.1)
boxplot(as(nc, "data.frame")[,c("SMR", "EBPG", "EBLN", "EBMarshall", "EBMrshloc")], cex.lab=.5, las=1, horizontal=TRUE)
#boxplot(as(nc, "data.frame")[,c("SMR", "EBPG", "EBLN", "EBMarshall", "EBMrshloc")], cex.lab=.5, las=3)
#boxplot(nc@data[,c("SMR", "EBPG", "EBLN", "EBMarshall", "EBMrshloc")],las=2)
par(oopar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 21: 
###################################################
#Data loaded from a previous Windows session
#This file contains the output of MCMC simulations of the P-G model
#I have used BRugs and OpenBUGS to get the results
#
#load("PG.RData")
con <- url("http://www.bias-project.org.uk/ASDAR/PG.RData", open="rb")
load(con)
close(con)


###################################################
### chunk number 22: 
###################################################
#Run some models
library(R2WinBUGS)

N<-length(nc$Observed)
d<-list(N=N, observed=nc$Observed, expected=nc$Expected)


###################################################
### chunk number 23:  eval=FALSE
###################################################
## pgmodelfile<-paste(getwd(), "/PG-model.txt", sep="")
## wdir<-paste(getwd(), "/PG", sep="")
## if(!file.exists(wdir)){dir.create(wdir)}
## BugsDir <- "/home/asdar/.wine/dosdevices/c:/Program Files/WinBUGS14"
## MCMCres<- bugs(data=d, inits=list(list(nu=1, alpha=1)), 
## working.directory=wdir,
## parameters.to.save=c("theta", "nu", "alpha"),
## n.chains=1, n.iter=20000, n.burnin=10000, n.thin=10,
## model.file=pgmodelfile,
## bugs.directory=BugsDir,
## WINEPATH="/usr/bin/winepath")
## 
## #     save(file="PG.RData", list=c("d", "MCMCres") )


###################################################
### chunk number 24: 
###################################################
#Add MCMC estimates
nc$PGmean<-MCMCres$mean$theta
nc$PGmedian<-MCMCres$median$theta
#MCMCres[1:2,1:6]


###################################################
### chunk number 25: 
###################################################
.iwidth <- 6
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
plot(1,1, type="n", xlim=c(1,100), ylim=c(0,4),
  main= "Credible intervals of the relative risks",
  xlab="County", ylab="Relative Risk", xaxt="n")
abline(h=1, lty=2)

for(i in 1:100)
{
	if(MCMCres$summary[i,3]>1 )
	{
		col<-gray(.4)
		lty<-2
		text(i, MCMCres$summary[i,7]+.31, nc$NAME[i], 
		srt=90, col=gray(.4), cex=.85)
	}
	else
	{
		col<-"black"
lty<-1
	}

	lines(c(i,i), c(MCMCres$summary[i,3], MCMCres$summary[i,7]), col=col, lty=lty)
	points(i, MCMCres$median$theta[i], pch=18, col=col)
#	points(i, MCMCres$mean[2+i], pch=20, col=col)
#	legend(20,4, legend=c("Median", "Mean"), pch=c(18,20))
}
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 26: 
###################################################
.iwidth <- 7
.iheight <- 4
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)

print(spplot(nc, c("SMR","EBPG", "PGmean", "PGmedian"),
  col.regions=cols,  at=brks, axes = TRUE, colorkey=colorkey))

#cat(c("alpha =",eb$alpha, "\n"))
#cat(c("nu=", eb$nu, "\n"))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 27: 
###################################################
.iwidth <- 6
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
nc$RATIO<-nc$NWBIR74/nc$BIR74
#print(spplot(nc, "RATIO", col.regions=brewer.pal(9,"Blues"), at=.8*0:9/9))
#print(spplot(nc, "RATIO", col.regions=brewer.pal(9,"Greys"), at=.8*0:9/9))
# RSB quietening greys
print(spplot(nc, "RATIO", col.regions=grey.colors(4, 0.95, 0.55, 2.2), at=.8*0:4/4))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 28: 
###################################################
#This will ensure that data and neighbour list are the same
#idx<-match(attr(ncCR85.nb, "region.id"), nc$"CNTY_ID")
#
#nc.nb<-ncCR85
#nc.nb<-nc.nb[order(idx)]
#
#nc.nb<-lapply(nc.nb, function(X, idx){idx[X]}, idx=(idx))
#class(nc.nb)<-"nb"
#nc.nb<-nc.nb[(order(idx))]
#nc.nb<-nb2WB(nc.nb)
nc.nb<-nb2WB(ncCR85)


###################################################
### chunk number 29: 
###################################################
options("width"=60)


###################################################
### chunk number 30: 
###################################################
nc$nwprop<-nc$NWBIR74/nc$BIR74

d<-list(N=N, observed=nc$Observed, expected=nc$Expected,
  nonwhite=nc$nwprop,#log(nwprop/(1-nwprop)),
  adj=nc.nb$adj,  weights=nc.nb$weights, num=nc.nb$num)

dwoutcov<-list(N=N, observed=nc$Observed, expected=nc$Expected,
  adj=nc.nb$adj,  weights=nc.nb$weights, num=nc.nb$num)

inits<-list(u=rep(0,N), v=rep(0,N), alpha=0, beta=0, precu=.001, precv=.001)
#inits$v[d$num==0]<-NA


###################################################
### chunk number 31: 
###################################################
options("width"=70)


###################################################
### chunk number 32:  eval=FALSE
###################################################
## 
## 
## bymmodelfile<-paste(getwd(), "/BYM-model.txt", sep="")
## wdir<-paste(getwd(), "/BYM", sep="")
## if(!file.exists(wdir)){dir.create(wdir)}
## BugsDir <- "/home/asdar/.wine/dosdevices/c:/Program Files/WinBUGS14"


###################################################
### chunk number 33: 
###################################################
options("width"=60)


###################################################
### chunk number 34:  eval=FALSE
###################################################
## MCMCres<- bugs(data=d, inits=list(inits),
## working.directory=wdir,
## parameters.to.save=c("theta", "alpha", "beta", "u", "v", "sigmau", "sigmav"),
## n.chains=1, n.iter=30000, n.burnin=20000, n.thin=10,
## model.file=bymmodelfile,
## bugs.directory=BugsDir,
## WINEPATH="/usr/bin/winepath")
## 
## 
## #     save(file="BYM.RData", list=c("d", "inits", "MCMCres") )


###################################################
### chunk number 35: 
###################################################
options("width"=70)


###################################################
### chunk number 36: 
###################################################
#load("BYM.RData")
con <- url("http://www.bias-project.org.uk/ASDAR/BYM.RData", open="rb")
load(con)
close(con)


###################################################
### chunk number 37: 
###################################################
#Load the data obtained by running WinBUGS in Windows
nc$BYMmean<-MCMCres$mean$theta
#nc$BYMmedian<-MCMCres$median$theta
nc$BYMumean<-MCMCres$mean$u
#nc$BYMumedian<-MCMCres$median$u
nc$BYMvmean<-MCMCres$mean$v
#nc$BYMvmedian<-MCMCres$median$v


###################################################
### chunk number 38: 
###################################################
#MCMCres[1:offset,1:6]


###################################################
### chunk number 39: 
###################################################
.iwidth <- 5
.iheight <- 6.5
.ipointsize <- 8
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
library(coda)
#ncoutput<-read.coda("BYM/coda1.txt", "BYM/codaIndex.txt")

plot(ncoutput[,c("deviance", "alpha", "beta", "theta[94]")])

     save(file="BYM.RData", list=c("d", "inits", "MCMCres", "ncoutput") )

dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 40: 
###################################################
geweke.diag(ncoutput[,c("deviance", "alpha", "beta", "theta[94]")])


###################################################
### chunk number 41: 
###################################################
.iwidth <- 5
.iheight <- 5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
print(spplot(nc, c("SMR", "BYMmean"), at=brks, col.regions=cols,
  axes=TRUE, colorkey=colorkey))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 42: 
###################################################
.iwidth <- 6
.iheight <- 4.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
ubrks<-quantile(nc$BYMumean, seq(0,1,1/5))
ubrks[6] <- ubrks[6]*1.01
#colorkey$labels<-as.character(formatC(ubrks,3))
colorkey$labels<-as.character(trunc(1000*ubrks)/1000)
colorkey$at<-ubrks[1]+ubrks[6]*0:5/5
print(spplot(nc, "BYMumean", at= ubrks, axes=TRUE, 
col.regions=grey.colors(5, 0.95, 0.55, 2.2), colorkey=colorkey))
# RSB quietening greys
#col.regions=brewer.pal(5, "Greys"), colorkey=colorkey))
#col.regions=brewer.pal(5, "BuGn"), colorkey=colorkey))
#print(spplot(nc, c("BYMumean", "BYMumedian"), axes=TRUE))

#Display a prettier map
#  at=c(-260,-215,-210,-200,-150,-210,-100,-50,0,50), 
#  col.regions=rev(brewer.pal(9, "Purples"))) )
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 43: 
###################################################
.iwidth <- 6
.iheight <- 4.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
vbrks<-quantile(nc$BYMvmean, seq(0,1,1/5))
vbrks[6] <-vbrks[6]*1.01
#vbrks[5]<-0.205#Fix to get one area inside the range
#colorkey$labels<-as.character(formatC(vbrks,3))
colorkey$labels<-as.character(trunc(1000*vbrks)/1000)
colorkey$at<-vbrks[1]+vbrks[6]*0:5/5
print(spplot(nc, "BYMvmean", at=vbrks, axes=TRUE,
   col.regions=grey.colors(5, 0.95, 0.55, 2.2), colorkey=colorkey))
# RSB quietening greys
#   col.regions=brewer.pal(5, "Greys"), colorkey=colorkey))
#   col.regions=brewer.pal(5, "BuPu"), colorkey=colorkey))
#print(spplot(nc, c("BYMvmean", "BYMvmedian"), axes=TRUE))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 44: 
###################################################
.iwidth <- 6
.iheight <- 3.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
plot(1,1, type="n", xlim=c(1,100), ylim=c(0,4.5),
  main= "Credible intervals of the relative risks",
  xlab="County", ylab="Relative Risk", xaxt="n")
abline(h=1, lty=2)

for(i in 1:100)
{
        if(MCMCres$summary[i,3]>1 )
        {
                col<-gray(.4)
		lty<-2
                text(i, MCMCres$summary[i,7]+.31, nc$NAME[i],
                srt=90, col=gray(.4), cex=.85)
        }
        else
        {
                col<-"black"
		lty<-1
        }

        lines(c(i,i), c(MCMCres$summary[i,3], MCMCres$summary[i,7]), col=col, lty=lty)
#       points(i, MCMCres$mean[offset+i], pch=20, col=col)
#       legend(20,4, legend=c("Median", "Mean"), pch=c(18,20))
}

points(1:100, MCMCres$median$theta, pch=18, col=col)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 45:  eval=FALSE
###################################################
## #
## #I do not know why but if these objects are not removed I 
## #get 10 pages of warnings!!
## #
## rm(list=c("EBMarshall","Observed","Expected","SMR"))


###################################################
### chunk number 46: 
###################################################
set.seed(1)


###################################################
### chunk number 47: 
###################################################
chtest<-achisq.test(Observed~offset(log(Expected)), as(nc, "data.frame"), "multinom", 999)
chtest


###################################################
### chunk number 48: 
###################################################
1- pchisq(chtest$t0, 100-1)


###################################################
### chunk number 49: 
###################################################
set.seed(1)


###################################################
### chunk number 50: 
###################################################
pwtest<-pottwhitt.test(Observed~offset(log(Expected)), as(nc, "data.frame"), "multinom", 999)


###################################################
### chunk number 51: 
###################################################
Oplus<- sum(nc$Observed)
1- pnorm(pwtest$t0, Oplus*(Oplus-1), sqrt(2*100*Oplus*(Oplus-1)))


###################################################
### chunk number 52: 
###################################################
col.W <- nb2listw(ncCR85, zero.policy=TRUE)


#mn<-boot(as(nc, "data.frame"), statistic=moranI.pboot, sim="parametric",
#             ran.gen=negbin.sim, R=99,  listw=col.W, n=length(ncCR85.nb),
#             S0=Szero(col.W) )
#dcluster.test(mn)

moranI.test(Observed~offset(log(Expected)), as(nc, "data.frame"), "negbin", 999, 
   listw=col.W, n=length(ncCR85), S0=Szero(col.W) )


###################################################
### chunk number 53: 
###################################################
data(nc.sids)
#Centroids
idx<-match(nc$NAME, rownames(nc.sids))

nc$x<-nc.sids$x[idx]
nc$y<-nc.sids$y[idx]
#Calculate neighbours based on distance
coords<-cbind(nc$x, nc$y)
#as.matrix(as(nc, "data.frame")[,c("x", "y")])

dlist<-dnearneigh(coords, 0, Inf)
dlist<-include.self(dlist)
dlist.d<-nbdists(dlist, coords)

#Calculate weights. They are globally standardised but it doesn't
#change significance.
phi<-100
col.W.tango<-nb2listw(dlist, glist=lapply(dlist.d, function(x, phi) {exp(-x/phi)}, 
  phi=phi), style="C")


###################################################
### chunk number 54: 
###################################################
set.seed(1)


###################################################
### chunk number 55: 
###################################################
#tn<-boot(nc@data, statistic=tango.pboot, sim="parametric",
#  ran.gen=negbin.sim, R=99, listw=col.W.tango, zero.policy=TRUE)
#dcluster.test(tn)
tango.test(Observed~offset(log(Expected)), as(nc, "data.frame"), "negbin", 999, 
   listw=col.W.tango, zero.policy=TRUE)


###################################################
### chunk number 56: 
###################################################
sidsgam<-opgam(data=as(nc, "data.frame"),  radius=30, step=10, alpha=.002)
gampoints<-SpatialPoints(sidsgam[,c("x", "y")]*1000, 
   CRS("+proj=utm +zone=18 +datum=NAD27"))


###################################################
### chunk number 57: 
###################################################
library(rgdal)
ll <- CRS("+proj=longlat +datum=NAD27")
gampoints<-spTransform(gampoints, ll)
gam.layout<-list("sp.points", gampoints)


###################################################
### chunk number 58: 
###################################################
.iwidth <- 5
.iheight <- 2.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
#Plot centroids
plot(nc, xlab="Easting", ylab="Northing")
#Plot points marked as clusters
points(gampoints, sidsgam$y, col=gray(.4), pch=4)

dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 59: 
###################################################
set.seed(1234)


###################################################
### chunk number 60: 
###################################################
options("width"=60)


###################################################
### chunk number 61: 
###################################################
mle <- calculate.mle(as(nc, "data.frame"), model="negbin")
thegrid <- as(nc, "data.frame")[,c("x","y")]
knresults<-opgam(data=as(nc, "data.frame"), thegrid=thegrid, alpha=.05,
   iscluster=kn.iscluster, fractpop=0.15, R=99, model="negbin", mle=mle)


###################################################
### chunk number 62: 
###################################################
options("width"=70)


###################################################
### chunk number 63: 
###################################################
.iwidth <- 5
.iheight <- 2.5
.ipointsize <- 10
.epsNo <- .epsNo + 1; file <- paste("Fig-dismap-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
#Plot all centroids and significant ones in red
#print(knresults)
#plot(nc@data$x, nc@data$y, main="Kulldorff and Nagarwalla's method",
#  xlab="Easting", ylab="Northing")

clusters<-get.knclusters(as(nc, "data.frame"), knresults)
i<-which.max(knresults$statistic)

nc$KNcluster<-""
nc$KNcluster[clusters[[i]]]<-"cluster"
nc$KNcluster[clusters[[i]][1]]<-"centre"
nc$KNcluster<-as.factor(nc$KNcluster)

print(spplot(nc, "KNcluster", main="Kulldorff's method",
  xlab="Easting", ylab="Northing", col.regions=c(gray(1), gray(.5), gray(.8))))
#  xlab="Easting", ylab="Northing", col.regions=c("white", "#E34A33","#FDBB84")))

#c("white", "red", "magenta") ))

dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 64: 
###################################################
save(file="dismap.RDat", list=ls())


###################################################
### chunk number 65: 
###################################################
stone.stat(as(nc, "data.frame"), region=which(nc$NAME=="Anson"))
st<-stone.test(Observed~offset(log(Expected)), as(nc, "data.frame"), model="negbin", 99, 
   region=which(nc$NAME=="Anson"))
st


###################################################
### chunk number 66: 
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
sT <- capture.output(print(Sys.time()))
cat("%", sT, "\n")


