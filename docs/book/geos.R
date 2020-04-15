###################################################
### chunk number 1: 
###################################################
rm(list=ls())
library(digest)
.owidth <- getOption("width")
options("width"=70)
owarn <- options("warn")$warn
options(warn=1)
.epsNo <- 0


###################################################
### chunk number 2: figreset eval=FALSE
###################################################
## .iwidth <- 5
## .iheight <- 6
## .ipointsize <- 12


###################################################
### chunk number 3: 
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 4: afig eval=FALSE
###################################################
## .epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
## postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)


###################################################
### chunk number 5: zfig eval=FALSE
###################################################
## dev.null <- dev.off()
## system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 6: afig_png eval=FALSE
###################################################
## .epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".png", sep="")
## png(filename=file, width = .iwidth, height = .iheight, pointsize = .ipointsize)


###################################################
### chunk number 7: zfig_png eval=FALSE
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### chunk number 8:  eval=FALSE
###################################################
## library(lattice)
## library(sp)
## data(meuse)
## coordinates(meuse) <- c("x", "y")
## spplot(meuse, "zinc", do.log = T)
## bubble(meuse, "zinc", do.log = T, key.space = "bottom")


###################################################
### chunk number 9: 
###################################################
library(lattice)
trellis.par.set(canonical.theme(color = FALSE))
library(sp)
data(meuse)
coordinates(meuse) <- c("x", "y")
#print(spplot(meuse, "zinc", do.log = T), split = c(1,1,2,1), more = TRUE)
#print(bubble(meuse, "zinc", do.log = T, key.space = "bottom"), split = c(2,1,2,1))


###################################################
### chunk number 10:  eval=FALSE
###################################################
## xyplot(log(zinc)~sqrt(dist), as.data.frame(meuse))
## zn.lm <- lm(log(zinc)~sqrt(dist), meuse)
## meuse$fitted.s <- predict(zn.lm, meuse) - mean(predict(zn.lm, meuse))
## meuse$residuals <- residuals(zn.lm)
## spplot(meuse, c("fitted.s", "residuals"))


###################################################
### chunk number 11: 
###################################################
.iwidth <- 6
.iheight <- 4
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
print(xyplot(log(zinc)~sqrt(dist), as.data.frame(meuse), asp = 1), split = c(1,1,2,1), more = TRUE)
zn.lm <- lm(log(zinc)~sqrt(dist), meuse)
meuse$fitted.s <- predict(zn.lm, meuse) - mean(predict(zn.lm, meuse))
meuse$residuals <- residuals(zn.lm)
# RSB quietening greys
print(spplot(meuse, c("fitted.s", "residuals"), cuts=4, col.regions=grey.colors(5, 0.95, 0.55, 2.2)), split = c(2,1,2,1))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 12: 
###################################################
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")


###################################################
### chunk number 13: 
###################################################
library(gstat)
idw.out <- idw(zinc~1, meuse, meuse.grid, idp = 2.5)
as.data.frame(idw.out)[1:5,]


###################################################
### chunk number 14: 
###################################################
zn.lm <- lm(log(zinc)~sqrt(dist), meuse)
meuse.grid$pred <- predict(zn.lm, meuse.grid)
meuse.grid$se.fit <- predict(zn.lm, meuse.grid, se.fit=TRUE)$se.fit
#$


###################################################
### chunk number 15: 
###################################################
meuse.lm <- krige(log(zinc)~sqrt(dist), meuse, meuse.grid)


###################################################
### chunk number 16: 
###################################################
meuse.tr2 <- krige(log(zinc)~1, meuse, meuse.grid, degree = 2)


###################################################
### chunk number 17:  eval=FALSE
###################################################
## lm(log(zinc)~I(x^2)+I(y^2)+I(x*y) + x + y, meuse)


###################################################
### chunk number 18:  eval=FALSE
###################################################
## lm(log(zinc)~poly(x,y,2), meuse)


###################################################
### chunk number 19:  eval=FALSE
###################################################
## hscat(log(zinc)~1,meuse,(0:9)*100)


###################################################
### chunk number 20: 
###################################################
.iwidth <- 6
.iheight <- 6
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
hscat(log(zinc)~1,meuse,(0:9)*100,pch=3,cex=.6)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 21:  eval=FALSE
###################################################
## library(gstat)
## variogram(log(zinc) ~ 1, meuse, cloud = TRUE)


###################################################
### chunk number 22: 
###################################################
.iwidth <- 6
.iheight <- 6
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
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
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 23:  eval=FALSE
###################################################
## plot(variogram(log(zinc) ~ 1, meuse))


###################################################
### chunk number 24:  eval=FALSE
###################################################
## sel <- plot(variogram(zinc ~ 1, meuse, cloud = TRUE), digitize = TRUE)
## plot(sel, meuse)


###################################################
### chunk number 25: 
###################################################
.iwidth <- 6
.iheight <- 4
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
sel <-
structure(list(x = c(145.291968730077, 266.107479142605, 320.156523274526, 
339.232656497557, 323.335878811698, 212.058435010685, 135.753902118561, 
46.7319470777507, 78.5255024494688, 142.112613192905), y = c(574649.690841889, 
581256.265954825, 627502.29174538, 822396.257577002, 1053626.38652977, 
1278249.94036961, 1255126.92747433, 792666.669568789, 634108.866858316, 
577952.978398357)), .Names = c("x", "y"))
v <- variogram(zinc ~ 1, meuse, cloud = TRUE)
v$gamma <- v$gamma/1e6
sel$y <- sel$y/1e6
p1 <- xyplot(gamma~dist, v, 
	panel = function(x, y, ...) {
		panel.xyplot(x, y, ...)
		llines(sel$x, sel$y)
	},
	pch=3, cex = .5, asp = 1, ylab = "gamma (x 1e6)")
x <-
structure(list(head = c(40, 40, 40, 54, 55, 54, 47, 80, 55, 55, 
54, 53, 54, 55, 59, 59), tail = c(41, 42, 43, 57, 57, 58, 59, 
99, 121, 122, 123, 125, 125, 125, 125, 132)), .Names = c("head", 
"tail"), row.names = as.integer(c(NA, 16)), class = c("pointPairs", 
"data.frame"))
p2 = plot(x, meuse, scales=list(draw=F), col.line = 1)
print(p1, split = c(1,1,2,1), more = TRUE)
print(p2, split = c(2,1,2,1))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 26: 
###################################################
.iwidth <- 6
.iheight <- 4
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
v <- variogram(log(zinc) ~ 1, meuse)
# INTERACTIVE mode out-commented:
#plot(v, type = 'b', pch = 3)
#fn = function(n = 100) {
#        for (i in 1:n) {
#        	meuse$random = sample(meuse$zinc)
#        	v = variogram(log(random) ~ 1, meuse)
#        	trellis.focus("panel", 1, 1, highlight = FALSE)
#        	llines(v$dist, v$gamma, col = 'grey')
#        	trellis.unfocus()
#        }
#}
#fn()
#trellis.focus("panel", 1, 1, highlight = FALSE)
#lpoints(v$dist, v$gamma, col = 'black', type = 'b', lwd = 2, pch=3)
#trellis.unfocus()
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
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 27:  eval=FALSE
###################################################
## plot(variogram(log(zinc) ~ 1, meuse))


###################################################
### chunk number 28:  eval=FALSE
###################################################
## plot(variogram(log(zinc) ~ 1, meuse, alpha = c(0, 45, 90, 135)))


###################################################
### chunk number 29:  eval=FALSE
###################################################
## plot(variogram(log(zinc) ~ 1, meuse, cutoff = 1000, width = 50))


###################################################
### chunk number 30:  eval=FALSE
###################################################
## variogram(log(zinc) ~ 1, meuse, boundaries = c(0,50,100,seq(250,1500,250)))


###################################################
### chunk number 31:  eval=FALSE
###################################################
## show.vgms()
## show.vgms(model = "Mat", kappa.range = c(.1, .2, .5, 1, 2, 5, 10), max = 10)


###################################################
### chunk number 32: 
###################################################
vgm(1, "Sph", 300)
vgm(1, "Sph", 300, 0.5)
v1 <- vgm(1, "Sph", 300, 0.5)
v2 <- vgm(0.8, "Sph", 800, add.to = v1)
v2
vgm(0.5, "Nug", 0)


###################################################
### chunk number 33: 
###################################################
vgm()


###################################################
### chunk number 34:  eval=FALSE
###################################################
## v <- variogram(log(zinc) ~ 1, meuse)
## plot(v)


###################################################
### chunk number 35: 
###################################################
fit.variogram(v, vgm(1, "Sph", 800, 1))


###################################################
### chunk number 36: 
###################################################
fit.variogram(v, vgm(1, "Sph", 10, 1))


###################################################
### chunk number 37:  eval=FALSE
###################################################
## v.fit <- fit.variogram(v, vgm(1, "Sph", 10, 1))
## if (attr(v.fit, "singular"))
## 	stop("singular fit")


###################################################
### chunk number 38: 
###################################################
.iwidth <- 6
.iheight <- 3.5
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
v <- variogram(log(zinc) ~ 1, meuse)
v.fit <- fit.variogram(v, vgm(1, "Sph", 800, 1))
plot(v, v.fit, pch = 3)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 39: 
###################################################
attr(v.fit, "SSErr")


###################################################
### chunk number 40:  eval=FALSE
###################################################
## library(geoR)
## v.eye <- eyefit(variog(as.geodata(meuse["zinc"]),max.dist=1500))
## ve.fit <- as.vgm.variomodel(v.eye[[1]])


###################################################
### chunk number 41: 
###################################################
fit.variogram(v, vgm(1, "Sph", 800, 0.06), fit.sills = c(FALSE, TRUE))


###################################################
### chunk number 42: 
###################################################
fit.variogram.reml(log(zinc)~1, meuse, model=vgm(0.6, "Sph", 800, 0.06))


###################################################
### chunk number 43: 
###################################################
v.dir <- variogram(log(zinc)~1,meuse,alpha=(0:3)*45) 
v.anis <- vgm(.6, "Sph", 1600, .05, anis=c(45, 0.3))


###################################################
### chunk number 44:  eval=FALSE
###################################################
## plot(v.dir, v.anis)


###################################################
### chunk number 45: 
###################################################
.iwidth <- 6
.iheight <- 3
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
v.dir <- variogram(log(zinc)~1,meuse,alpha=(0:3)*45) 
v.anis <- vgm(.6, "Sph", 1600, .05, anis=c(45,.3))
print(plot(v.dir, v.anis, pch=3))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 46:  eval=FALSE
###################################################
## plot(variogram(log(zinc)~1,meuse, map=TRUE, cutoff=1000, width=100))


###################################################
### chunk number 47: 
###################################################
g <- gstat(NULL, "logCd", log(cadmium)~1, meuse)
g <- gstat(g, "logCu", log(copper)~1, meuse)
g <- gstat(g, "logPb", log(lead)~1, meuse)
g <- gstat(g, "logZn", log(zinc)~1, meuse)
g 
vm <- variogram(g)
vm.fit <- fit.lmc(vm, g, vgm(1, "Sph", 800, 1))


###################################################
### chunk number 48:  eval=FALSE
###################################################
## plot(vm, vm.fit)


###################################################
### chunk number 49: 
###################################################
cor(as.data.frame(meuse)[c("cadmium", "copper", "lead", "zinc")])


###################################################
### chunk number 50: 
###################################################
.iwidth <- 7
.iheight <- 7
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
print(plot(vm, vm.fit))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 51:  eval=FALSE
###################################################
## variogram(log(zinc) ~ sqrt(dist), meuse)


###################################################
### chunk number 52: 
###################################################
f <- log(zinc) ~ sqrt(dist)
vt <- variogram(f, meuse)
vt.fit <- fit.variogram(vt, vgm(1, "Exp", 300, 1))
vt.fit
g.wls <- gstat(NULL, "log-zinc", f, meuse, model=vt.fit, set = list(gls=1))
(variogram(g.wls)$gamma - vt$gamma)/ mean(vt$gamma)
#$


###################################################
### chunk number 53: 
###################################################
lz.sk <- krige(log(zinc)~1, meuse, meuse.grid, v.fit, beta = 5.9)
lz.ok <- krige(log(zinc)~1, meuse, meuse.grid, v.fit)
lz.uk <- krige(log(zinc)~sqrt(dist), meuse, meuse.grid, vt.fit)


###################################################
### chunk number 54: 
###################################################
cok.maps <- predict(vm.fit, meuse.grid)
names(cok.maps)


###################################################
### chunk number 55:  eval=FALSE
###################################################
## spplot.vcov(cok.maps)


###################################################
### chunk number 56: 
###################################################
.iwidth <- 8
.iheight <- 8
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
print(spplot.vcov(cok.maps, cuts=8,
 col.regions=grey.colors(9, 0.9, 0.45, 2.2)))
# RSB quietening greys
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 57:  eval=FALSE
###################################################
## vm2.fit <- vm.fit
## vm2.fit$model[[3]]$range  = c(0, 900)
## predict(vm2.fit, meuse.grid)


###################################################
### chunk number 58: 
###################################################
vm2.fit <- vm.fit
vm2.fit$model[[3]]$range  = c(0, 900)


###################################################
### chunk number 59: 
###################################################
vm2.fit$set <- list(nocheck=1) #$
x <- predict(vm2.fit, meuse.grid)
names(as.data.frame(x))
any(as.data.frame(x)[c(2,4,6,8)] < 0)


###################################################
### chunk number 60: 
###################################################
g.cc <- gstat(NULL, "log.zinc", log(zinc)~1, meuse, model = v.fit)
meuse.grid$distn <- meuse.grid$dist - mean(meuse.grid$dist) + mean(log(meuse$zinc))
vd.fit <- v.fit
vov <- var(meuse.grid$distn) / var(log(meuse$zinc))
vd.fit$psill <- v.fit$psill * vov
g.cc <- gstat(g.cc, "distn", distn ~ 1, meuse.grid, nmax = 1, model=vd.fit, 
	merge = c("log.zinc","distn"))
vx.fit <- v.fit
vx.fit$psill <- sqrt(v.fit$psill * vd.fit$psill) * cor(meuse$dist, log(meuse$zinc)) #$
g.cc <- gstat(g.cc, c("log.zinc", "distn"), model = vx.fit) 
x <- predict(g.cc, meuse.grid)


###################################################
### chunk number 61: 
###################################################
.iwidth <- 6
.iheight <- 6
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
x$lz.uk <- lz.uk$var1.pred
x$lz.ok <- lz.ok$var1.pred
print(spplot(x, c("log.zinc.pred", "lz.ok", "lz.uk"),
	names.attr = c("collocated", "ordinary", "universal"),
	cuts=7, col.regions=grey.colors(8, 0.55, 0.95, 2.2)
))
# RSB quietening greys
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 62:  eval=FALSE
###################################################
## demo(pcb)


###################################################
### chunk number 63: 
###################################################
lz.ok <- krige(log(zinc)~1, meuse, meuse.grid, v.fit, block = c(40, 40))


###################################################
### chunk number 64: 
###################################################
xy <- expand.grid(x = seq(-20, 20, 4), y = seq(-20, 20, 4))
xy <- xy[(xy$x^2 + xy$y^2) <= 20^2, ]
lz.ok <- krige(log(zinc)~1, meuse, meuse.grid, v.fit, block = xy)


###################################################
### chunk number 65:  eval=FALSE
###################################################
## lz.pols <- krige(log(zinc)~1, meuse, meuse.polygons, v.fit)


###################################################
### chunk number 66:  eval=FALSE
###################################################
## spsample(polygon, n = 500, type = "regular", offset = c(.5, .5))


###################################################
### chunk number 67: 
###################################################
options("width"=55)


###################################################
### chunk number 68: 
###################################################
meuse$part.a <- idw(part.a~1, meuse.grid, meuse, nmax=1)$var1.pred


###################################################
### chunk number 69: 
###################################################
options("width"=70)


###################################################
### chunk number 70: 
###################################################
meuse$part.a <- meuse.grid$part.a[overlay(meuse.grid, meuse)]


###################################################
### chunk number 71: 
###################################################
x1 <- krige(log(zinc)~1, meuse[meuse$part.a == 0,],
meuse.grid[meuse.grid$part.a == 0,], model = vgm(.548, "Sph", 900, .0654),
nmin = 20, nmax = 40, maxdist = 1000)
x2 <- krige(log(zinc)~1, meuse[meuse$part.a == 1,],
meuse.grid[meuse.grid$part.a == 1,], model = vgm(.716, "Sph", 900),
nmin = 20, nmax = 40, maxdist = 1000)
lz.stk <- rbind(as.data.frame(x1), as.data.frame(x2))
coordinates(lz.stk) <- c("x", "y")
lz.stk <- as(x, "SpatialPixelsDataFrame")


###################################################
### chunk number 72:  eval=FALSE
###################################################
## spplot(lz.stk["var1.pred"], main = "stratified kriging predictions")


###################################################
### chunk number 73: 
###################################################
g.tr <- gstat(formula = log(zinc) ~ sqrt(dist), data = meuse, model = v.fit)
predict(g.tr, meuse[1,])
predict(g.tr, meuse[1,], BLUE = TRUE)


###################################################
### chunk number 74:  eval=FALSE
###################################################
## predict(g, meuse[1,], BLUE = TRUE, debug = 32)


###################################################
### chunk number 75: 
###################################################
meuse$Int <- rep(1, 155) #$
g.tr <- gstat(formula = log(zinc) ~ -1+Int+sqrt(dist), data = meuse, model = v.fit)
rn <- c("Intercept", "beta1")
df <- data.frame(Int = c(0,1), dist = c(1,0), row.names=rn)
spdf <- SpatialPointsDataFrame(SpatialPoints(matrix(0, 2, 2)), df)
spdf
predict(g.tr, spdf, BLUE = TRUE)


###################################################
### chunk number 76:  eval=FALSE
###################################################
## library(MASS)
## boxcox(zinc~sqrt(dist), data=as.data.frame(meuse))


###################################################
### chunk number 77: 
###################################################
meuse$zinc.ns <- qqnorm(meuse$zinc, plot.it = FALSE)$x
#$


###################################################
### chunk number 78: 
###################################################
ind.f <- I(zinc < 500) ~ 1
ind.fit <- fit.variogram(variogram(ind.f, meuse), vgm(1, "Sph", 800, 1))
ind.kr <- krige(ind.f, meuse, meuse.grid, ind.fit)
summary(ind.kr$var1.pred)
#$


###################################################
### chunk number 79:  eval=FALSE
###################################################
## meuse.dup <- rbind(as.data.frame(meuse)[1,], as.data.frame(meuse))
## coordinates(meuse.dup)=~x+y
## krige(log(zinc)~1, meuse.dup, meuse[1,], v.fit)


###################################################
### chunk number 80: 
###################################################
meuse.dup <- rbind(as.data.frame(meuse)[1,], as.data.frame(meuse))
coordinates(meuse.dup)=~x+y


###################################################
### chunk number 81: 
###################################################
zd <- zerodist(meuse.dup)
zd
krige(log(zinc)~1, meuse.dup[-zd[,1], ], meuse[1,], v.fit)


###################################################
### chunk number 82: 
###################################################
setL <- list(cn_max=1e10)
krige(log(zinc)~1, meuse.dup, meuse[1,], v.fit, set = setL)


###################################################
### chunk number 83: 
###################################################
set.seed(1357531)


###################################################
### chunk number 84: 
###################################################
sel100 <- sample(1:155, 100)
m.model <- meuse[sel100,]
m.valid <- meuse[-sel100,]
v100.fit <- fit.variogram(variogram(log(zinc)~1, m.model), vgm(1, "Sph", 800, 1))
m.valid.pr <- krige(log(zinc)~1, m.model, m.valid, v100.fit)
resid.kr <- log(m.valid$zinc) - m.valid.pr$var1.pred
summary(resid.kr)
resid.mean <- log(m.valid$zinc) - mean(log(m.valid$zinc))
R2 <- 1 - sum(resid.kr^2)/sum(resid.mean^2)
R2


###################################################
### chunk number 85: 
###################################################
m.valid.pr$res <- resid.kr #$


###################################################
### chunk number 86:  eval=FALSE
###################################################
## bubble(m.valid.pr, "res")


###################################################
### chunk number 87: 
###################################################
nfold <- 3
part <- sample(1:nfold, 155, replace = TRUE)
sel <- (part != 1)
m.model <- meuse[sel,]
m.valid <- meuse[-sel,]


###################################################
### chunk number 88: 
###################################################
v.fit <- vgm(.59, "Sph", 874, .04)
cv155 <- krige.cv(log(zinc)~1, meuse, v.fit, nfold=5)


###################################################
### chunk number 89:  eval=FALSE
###################################################
## bubble(cv155, "residual", main = "log(zinc): 5-fold CV residuals")


###################################################
### chunk number 90: 
###################################################
.iwidth <- 6
.iheight <- 6
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
print(bubble(cv155, "residual", main = "log(zinc): 5-fold CV residuals",
	maxsize = 1.5, col = c(grey(.15), grey(.7))))
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 91: 
###################################################
summary(cv155)


###################################################
### chunk number 92:  eval=FALSE
###################################################
## g.cv <- gstat.cv(g, nmax=40)


###################################################
### chunk number 93: 
###################################################
v1.fit <- vgm(0.591, "Sph", 897, .0507)
v2.fit <- vgm(0.591, "Sph", 897, add.to = vgm(0.0507, "Sph", 40))


###################################################
### chunk number 94: 
###################################################
set.seed(13331)
cv155.1 <- krige.cv(log(zinc)~1, meuse, v1.fit, nfold=5)
set.seed(13331)
cv155.2 <- krige.cv(log(zinc)~1, meuse, v2.fit, nfold=5)
summary(cv155.1$residual-cv155.2$residual)


###################################################
### chunk number 95: 
###################################################
b1 <- krige(log(zinc)~1, meuse, meuse.grid, v1.fit, block = c(40,40))$var1.var
b2 <- krige(log(zinc)~1, meuse, meuse.grid, v2.fit, block = c(40,40))$var1.var
summary((b1-b2)/b1)


###################################################
### chunk number 96: 
###################################################
lzn.sim <- krige(log(zinc)~1, meuse, meuse.grid, v.fit, nsim = 6, nmax=40)


###################################################
### chunk number 97:  eval=FALSE
###################################################
## spplot(lzn.sim)


###################################################
### chunk number 98: 
###################################################
.iwidth <- 6
.iheight <- 5
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
#print(spplot(lzn.sim, col.regions=grey(c(.1,.3,.4,.5,.7,.8,.9,.95)),cuts=6))
print(spplot(lzn.sim, col.regions=grey.colors(7, 0.55, 0.95, 2.2),cuts=6))
# RSB quietening greys
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 99: 
###################################################
`area` <-
structure(c(181130.059662363, 180917.131281033, 180769.007189672, 
180676.429632572, 180611.625342602, 180463.501251241, 180389.439205561, 
180556.078808342, 180694.945143992, 180796.780456802, 180824.553723932, 
180898.615769613, 181000.451082423, 181130.059662363, 181241.152730884, 
181268.925998014, 181194.863952334, 181130.059662363, 333718.312421053, 
333491.307789474, 333245.386105263, 332990.005894737, 332781.918315789, 
332554.913684211, 332498.162526316, 332413.035789474, 332592.747789474, 
332810.293894737, 333008.922947368, 333188.634947368, 333358.888421053, 
333491.307789474, 333566.976, 333614.268631579, 333793.980631579, 
333718.312421053), .Dim = as.integer(c(18, 2)))
area.sp = SpatialPolygons(list(Polygons(list(Polygon(area)), "area")))


###################################################
### chunk number 100: 
###################################################
nsim <- 1000
cutoff <- 500
grd <- overlay(meuse.grid, area.sp)
sel.grid <- meuse.grid[!is.na(grd),]
lzn.sim <- krige(log(zinc)~1, meuse, sel.grid, v.fit, nsim = nsim, nmax=40)
res <- apply(as.data.frame(lzn.sim)[1:nsim], 2, function(x) mean(x > log(cutoff)))


###################################################
### chunk number 101:  eval=FALSE
###################################################
## hist(res, main = paste("fraction above", cutoff), xlab = NULL, ylab = NULL)


###################################################
### chunk number 102: 
###################################################
bkr <- krige(log(zinc)~1, meuse, area.sp, v.fit)
1 - pnorm(log(cutoff), bkr$var1.pred, sqrt(bkr$var1.var))


###################################################
### chunk number 103: 
###################################################
.iwidth <- 6
.iheight <- 3
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
meuse.grid$part <- meuse.grid$part.a + meuse.grid$part.b
layout(matrix(1:2, 1, 2))
omar <- par("mar")
par(mar = rep(0,4))
image(meuse.grid["part"], col = 'gray') #$
lines(area)
par(mar = c(2,2,0.5,0)+.1)
hist(res, main = NULL, xlab = NULL, ylab = NULL)
par(mar = omar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 104:  eval=FALSE
###################################################
## cok.sims <- predict(vm.fit, meuse.grid, nsim = 1000)


###################################################
### chunk number 105: 
###################################################
table(meuse$soil) #$
s1.fit <- fit.variogram(variogram(I(soil==1)~1,meuse), vgm(1, "Sph", 800, 1))
s1.sim <- krige(I(soil==1)~1, meuse, meuse.grid, s1.fit, nsim = 6, indicators = TRUE, nmax = 40)


###################################################
### chunk number 106:  eval=FALSE
###################################################
## spplot(s1.sim)


###################################################
### chunk number 107: 
###################################################
.iwidth <- 6
.iheight <- 5
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
print(spplot(s1.sim, cuts=1, col.regions=grey.colors(2, 0.6, 0.90, 2.2)))
# RSB quietening greys
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 108:  eval=FALSE
###################################################
## m1 <- sapply(1:155, function(x) mean(krige(log(zinc)~1, meuse[-x,], 
## 	meuse.grid, v.fit)$var1.var))
## which(m1 == min(m1)) #$


###################################################
### chunk number 109:  eval=FALSE
###################################################
## plot(sort(m1))


###################################################
### chunk number 110:  eval=FALSE
###################################################
## cutoff <- 1000
## f <- function(x) {
## 	kr = krige(log(zinc)~1, meuse[-x,], meuse.grid, v.fit)
## 	mean(abs(pnorm((kr$var1.pred - log(cutoff))/sqrt(kr$var1.var)) - 0.5))
## }
## m2 <- sapply(1:155, f)
## which(m2 == max(m2))


###################################################
### chunk number 111: 
###################################################
.iwidth <- 6
.iheight <- 4
.epsNo <- .epsNo + 1; file <- paste("Fig-geos-", .epsNo, ".eps", sep="")
postscript(file=file, onefile = TRUE, paper="special", width = .iwidth, height = .iheight, pointsize = .ipointsize, horizontal=FALSE)
source("m1m2.R")
layout(matrix(1:2, 1, 2))
omar <- par("mar")
par(mar = rep(0,4))
plot(meuse)
points(meuse[m1 < quantile(m1,.1),], pch=1)
points(meuse[m1 > quantile(m1,.9),], pch=16)
plot(meuse)
points(meuse[m2 < quantile(m2,.1),], pch=16)
points(meuse[m2 > quantile(m2,.9),], pch=1)
par(mar = omar)
dev.null <- dev.off()
system(paste("./mungeEps.sh", file, "3 > outfile ; mv outfile", file, sep=" "))
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### chunk number 112: 
###################################################
library(RandomFields)
ver <- packageDescription("RandomFields")$Version


###################################################
### chunk number 113: 
###################################################
PrintModelList()


###################################################
### chunk number 114:  eval=FALSE
###################################################
## library(geoR)
## plot(variog(as.geodata(meuse["zinc"]), max.dist=1500))


###################################################
### chunk number 115: 
###################################################
sI <- capture.output(print(sessionInfo()))
cat(paste("%", sI[1:(length(sI)/2)], sep=" "), sep="\n")
sT <- capture.output(print(Sys.time()))
cat("%", sT, "\n")


