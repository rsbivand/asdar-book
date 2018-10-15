v <- variogram(log(zinc) ~ 1, meuse)
print(xyplot(gamma ~ dist, v, pch = 3, type = 'b', lwd = 2,
	panel = function(x, y, ...) {
        for (i in 1:100) {
        	meuse$random = sample(meuse$zinc)
        	v = variogram(log(random) ~ 1, meuse)
        	llines(v$dist, v$gamma, col = 'grey')
		}
		panel.xyplot(x, y, ...)
	},
	ylim = c(0, 0.75), xlab = 'distance', ylab = 'semivariance'
))


