###################################################
# dismap_WB.R
# packages: sp, maptools, spdep, epitools, DCluster, rgdal
# datasets:
#
# This script runs dismap_mod.R and then some other code
#

#SET ME FIRST!!
#
#Set these variables to suit your local installation. See the manual
#page for function 'bugs' for details
#

#WinBUGS directory
BugsDir <- "/home/asdar/.wine/dosdevices/c:/Program Files/WinBUGS14"
#Path to 'winepath'. Not used under Windows.
winepath<-"/usr/bin/winepath"

#End of SET me FIRST!!


###################################################
### chunk number 0: 
###################################################
#Download files
if(!file.exists("dismap_mod.R")){
	download.file("https://asdar-book.org/book/dismap_mod.R","dismap_mod.R")	
}

if(!file.exists("PG-model.txt")){
	download.file("https://asdar-book.org/book/PG-model.txt","PG-model.txt")	
}

if(!file.exists("BYM-model.txt")){
	download.file("https://asdar-book.org/book/BYM-model.txt","BYM-model.txt")	
}

#Run 'dismap_mod.R' to get some results and data
source("dismap_mod.R")

if(.Platform$OS.type == "windows"){
	winepath<-NULL
}


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
pgmodelfile<-paste(getwd(), "/PG-model.txt", sep="")
wdir<-paste(getwd(), "/PG", sep="")
if(!file.exists(wdir)){dir.create(wdir)}

MCMCres<- bugs(data=d, inits=list(list(nu=1, alpha=1)), 
   working.directory=wdir,
   parameters.to.save=c("theta", "nu", "alpha"),
   n.chains=1, n.iter=20000, n.burnin=10000, n.thin=10,
   model.file=pgmodelfile,
   bugs.directory=BugsDir,
   WINEPATH=winepath)
 

###################################################
### chunk number 24: 
###################################################
#Add MCMC estimates
nc$PGmean<-MCMCres$mean$theta
nc$PGmedian<-MCMCres$median$theta


###################################################
### chunk number 25: 
###################################################
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


###################################################
### chunk number 26: 
###################################################
spplot(nc, c("SMR","EBPG", "PGmean", "PGmedian"),
  col.regions=cols,  at=brks, axes = TRUE, colorkey=colorkey)

cat(c("alpha =",eb$alpha, "\n"))
cat(c("nu=", eb$nu, "\n"))


###################################################
### chunk number 27: 
###################################################
nc$RATIO<-nc$NWBIR74/nc$BIR74

spplot(nc, "RATIO", col.regions=grey.colors(4, 0.95, 0.55, 2.2), at=.8*0:4/4)


###################################################
### chunk number 28: 
###################################################
nc.nb<-nb2WB(ncCR85)


###################################################
### chunk number 30: 
###################################################
nc$nwprop<-nc$NWBIR74/nc$BIR74

d<-list(N=N, observed=nc$Observed, expected=nc$Expected,
  nonwhite=nc$nwprop,
  adj=nc.nb$adj,  weights=nc.nb$weights, num=nc.nb$num)

dwoutcov<-list(N=N, observed=nc$Observed, expected=nc$Expected,
  adj=nc.nb$adj,  weights=nc.nb$weights, num=nc.nb$num)

inits<-list(u=rep(0,N), v=rep(0,N), alpha=0, beta=0, precu=.001, precv=.001)


###################################################
### chunk number 32:  eval=FALSE
###################################################
 
 
bymmodelfile<-paste(getwd(), "/BYM-model.txt", sep="")
wdir<-paste(getwd(), "/BYM", sep="")
if(!file.exists(wdir)){dir.create(wdir)}


###################################################
### chunk number 34:  eval=FALSE
###################################################
MCMCres<- bugs(data=d, inits=list(inits),
   working.directory=wdir,
   parameters.to.save=c("theta", "alpha", "beta", "u", "v", "sigmau", "sigmav"),
   n.chains=1, n.iter=30000, n.burnin=20000, n.thin=10,
   model.file=bymmodelfile,
   bugs.directory=BugsDir,
   WINEPATH=winepath)


###################################################
### chunk number 37: 
###################################################
nc$BYMmean<-MCMCres$mean$theta
#nc$BYMmedian<-MCMCres$median$theta
nc$BYMumean<-MCMCres$mean$u
#nc$BYMumedian<-MCMCres$median$u
nc$BYMvmean<-MCMCres$mean$v
#nc$BYMvmedian<-MCMCres$median$v


###################################################
### chunk number 39: 
###################################################
library(coda)
ncoutput<-read.coda("BYM/coda1.txt", "BYM/codaIndex.txt")

plot(ncoutput[,c("deviance", "alpha", "beta", "theta[94]")])


###################################################
### chunk number 40: 
###################################################
geweke.diag(ncoutput[,c("deviance", "alpha", "beta", "theta[94]")])


###################################################
### chunk number 41: 
###################################################
spplot(nc, c("SMR", "BYMmean"), at=brks, col.regions=cols,
  axes=TRUE, colorkey=colorkey)


###################################################
### chunk number 42: 
###################################################
ubrks<-quantile(nc$BYMumean, seq(0,1,1/5))
ubrks[6] <- ubrks[6]*1.01
colorkey$labels<-as.character(trunc(1000*ubrks)/1000)
colorkey$at<-ubrks[1]+ubrks[6]*0:5/5

spplot(nc, "BYMumean", at= ubrks, axes=TRUE, 
   col.regions=grey.colors(5, 0.95, 0.55, 2.2), colorkey=colorkey)


###################################################
### chunk number 43: 
###################################################
vbrks<-quantile(nc$BYMvmean, seq(0,1,1/5))
vbrks[6] <-vbrks[6]*1.01
colorkey$labels<-as.character(trunc(1000*vbrks)/1000)
colorkey$at<-vbrks[1]+vbrks[6]*0:5/5

spplot(nc, "BYMvmean", at=vbrks, axes=TRUE,
   col.regions=grey.colors(5, 0.95, 0.55, 2.2), colorkey=colorkey)


###################################################
### chunk number 44: 
###################################################
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


