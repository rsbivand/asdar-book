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
}


