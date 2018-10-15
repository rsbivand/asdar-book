print(xyplot((obs-theo)~r|DATASET , data=Kresults, type="l", 
   ylim= c(-.06, .06), ylab=expression(hat(K) (r)  - pi * r^2),
	panel=function(x, y, subscripts)
	{
		Ktheo<- Kresults$theo[subscripts]

		lpolygon(c(r, rev(r)), 
		   c(Kresults$lo[subscripts]-Ktheo, rev(Kresults$hi[subscripts]-Ktheo)),
		   border="gray", fill="gray"
		)

		llines(r, Kresults$obs[subscripts]-Ktheo, lty=2, lwd=1.5, col="black")	
	}
))


