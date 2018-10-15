library(gstat)
cld <- variogram(log(zinc) ~ 1, meuse, cloud = TRUE)
svgm <- variogram(log(zinc) ~ 1, meuse)
d <- data.frame(gamma = c(cld$gamma, svgm$gamma),
	dist = c(cld$dist, svgm$dist),
	id = c(rep("cloud", nrow(cld)), rep("sample variogram", nrow(svgm)))
	)
xyplot(gamma ~ dist | id, d,
	scales = list(y = list(relation = "free", ylim = list(NULL, c(-.005,0.7)))),
	layout = c(1, 2), as.table = TRUE,
	panel = function(x,y, ...) {
		if (panel.number() == 2)
			ltext(x+10, y, svgm$np, adj = c(0,0.5)) #$
		panel.xyplot(x,y,...)
	}, 
	xlim = c(0, 1590),
	cex = .5, pch = 3
)


