print(xyplot(obs~theo|DATASET , data=Fresults, type="l", 
	panel=function(x, y, subscripts)
	{
		lpolygon(c(x, rev(x)), 
		   c(Fresults$lo[subscripts], rev(Fresults$hi[subscripts])),
		   border="gray", fill="gray"
		)

		llines(x, y, col="black", lwd=2)
	}
))


