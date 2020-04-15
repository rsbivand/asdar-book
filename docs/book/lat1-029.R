dS <- c(0.75, 1, 1.5)*max_1nn
res <- sapply(nb_l, function(x) table(card(x)))
mx <- max(card(Sy13_nb))
res1 <- matrix(0, ncol=(mx+1), nrow=3)
rownames(res1) <- names(res)
colnames(res1) <- as.character(0:mx)
res1[1, names(res$d1)] <- res$d1
res1[2, names(res$d2)] <- res$d2
res1[3, names(res$d3)] <- res$d3
library(RColorBrewer)
pal <- grey.colors(3, 0.95, 0.55, 2.2)
barplot(res1, col=pal, beside=TRUE, legend.text=FALSE, xlab="numbers of neighbours", ylab="tracts")
legend("topright", legend=format(dS, digits=1), fill=pal, bty="n", cex=0.8, title="max. distance")


