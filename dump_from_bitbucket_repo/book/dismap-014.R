library(lattice)
trellis.par.set(canonical.theme(color = FALSE))
nc$pvalpois<-ppois(nc$Observed, nc$Expected, lower.tail=FALSE)

nbparam<-calculate.mle(as(nc, "data.frame"), model="negbin")
nc$pvalnegbin<-pnbinom(nc$Observed, size=nbparam$size, prob=nbparam$prob,
  lower.tail=FALSE)

colorkeypval<-list(labels=as.character(c(0, 0.01, 0.05, 0.1, .5, 1)), 
  at=(0:5)/5, height=.5)

pvalcols <- grey.colors(5, 0.95, 0.55, 2.2)

print(spplot(nc, c("pvalpois","pvalnegbin"), col.regions=rev(pvalcols), 
  at=c(0, 0.01, 0.05, 0.1, .5, 1), axes=TRUE, colorkey=colorkeypval ))



