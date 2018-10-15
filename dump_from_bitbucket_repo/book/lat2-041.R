image(SGDF, "fn", col=grey.colors(10, 1, 0.55, 2.2), axes=TRUE)
title(xlab=expression(lambda), ylab=expression(sigma^2))
contour(SGDF, "fn", add=TRUE)
points(c(lambda, nyerr$lambda), c(s2, nyerr$s2), pch=c(4, 3), lwd=2)
text(c(lambda, nyerr$lambda), c(s2, nyerr$s2), labels=c("GM", "ML"), pos=c(2, 4), offset=0.5)


