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
}


