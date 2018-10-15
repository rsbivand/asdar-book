


print(xyplot(obs~theo|DATASET , data=Gresults, type="l", 
	panel=function(x, y, subscripts)
	{
		lpolygon(c(x, rev(x)), 
		   c(Gresults$lo[subscripts], rev(Gresults$hi[subscripts])),
		   border="gray", fill="gray"
		)

		llines(x, y, col="black", lwd=2)
	}
))


