def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2),1,2))
plot(meuse.sr, axes = TRUE); title("axes = TRUE")
plot(meuse.sr, axes = FALSE); title("axes added")
axis(1, at = c(178000 + 0:2 * 2000), cex.axis = .7)
axis(2, at = c(326000 + 0:3 * 4000), cex.axis = .7)
box()
par(def.par)
cat("\n")


