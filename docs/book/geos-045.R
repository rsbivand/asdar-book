v.dir <- variogram(log(zinc)~1,meuse,alpha=(0:3)*45) 
v.anis <- vgm(.6, "Sph", 1600, .05, anis=c(45,.3))
print(plot(v.dir, v.anis, pch=3))


