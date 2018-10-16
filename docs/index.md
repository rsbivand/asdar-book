---
layout: default
---

This web site contains scripts and datasets to reproduce all the examples in _Applied Spatial Data Analysis with R, second edition_, by Roger Bivand, Edzer Pebesma and Virgilio Gomez-Rubio.

### Data sets and scripts

Data set bundles (.zip files for all datasets occuring in a chapter of the book), for chapter
[1](/bundles2ed/hello_bundle.zip),
[2](/bundles2ed/cm_bundle.zip),
[3](/bundles2ed/vis_bundle.zip),
[4](/bundles2ed/die_bundle.zip),
[5](/bundles2ed/cm2_bundle.zip),
[6](/bundles2ed/std_bundle.zip),
[7](/bundles2ed/sppa_bundle.zip),
[8](/bundles2ed/geos_bundle.zip),
[9](/bundles2ed/lat_bundle.zip),
[10](/bundles2ed/dismap_bundle.zip).

Verbatim (unchanged) book scripts, for chapter
[1](/book2ed/hello.R),
[2](/book2ed/cm.R),
[3](/book2ed/vis.R),
[4](/book2ed/die.R),
[5](/book2ed/cm2.R),
[6](/book2ed/std.R),
[7](/book2ed/sppa.R),
[8](/book2ed/geos.R),
[9](/book2ed/lat.R),
[10](/book2ed/dismap.R).

Updated and simplifed book scripts, for chapter
[1](/book2ed/hello_mod.R),
[2](/book2ed/cm_mod.R),
[3](/book2ed/vis_mod.R),
[4](/book2ed/die_mod.R),
[5](/book2ed/cm2_mod.R),
[6](/book2ed/std_mod.R),
[7](/book2ed/sppa_mod.R),
[8](/book2ed/geos_mod.R),
[9](/book2ed/lat_mod.R),
[10](/book2ed/dismap_mod.R).

### Reproducing the whole book

A script that downloads all scripts is:
```
# ASDAR_BOOK <- "http://www.asdar-book.org/book2ed"
ASDAR_BOOK <- "http://edzer.github.io/asdar-book/book2ed"
chapters <- c("hello", "cm", "vis", "die", "cm2",
"std", "sppa", "geos", "lat", "dismap")
for (i in chapters) {
  fn <- paste(i, "mod.R", sep="_")
  download.file(paste(ASDAR_BOOK, fn, sep = "/"), fn)
}
```

To run all examples of the book, a number of packages need to be installed. Running the following script will install those that are not already present:
```
pkgs <- c("boot", "CARBayes", "classInt", "coda", "cubature",
"DCluster", "deldir", "epitools", "geoR", "ggplot2", "gstat",
"INLA", "lattice", "latticeExtra", "lmtest", "maps", "maptools",
"MASS", "McSpatial", "mgcv", "nlme", "osmar", "pgirmess", "plm",
"R2BayesX", "R2WinBUGS", "raster", "RColorBrewer", "rgdal", "rgeos",
"sandwich", "sp", "spacetime", "spatstat", "spdep", "spgwr", "sphet",
"splancs", "xts")

for (p in pkgs) {
	if (inherits(try(library(p, character.only = TRUE)), "try-error"))
		install.packages(p, character.only = TRUE)
}
```

A script that downloads all data and scripts, extracts data, and reproduces the whole book is:
```
chapters <- c("hello", "cm", "vis", "die", "cm2",
"std", "sppa", "geos", "lat", "dismap")
for (i in chapters) {
  ASDAR_BOOK <- "http://edzer.github.io/asdar-book"
  fn <- paste(i, "mod.R", sep="_")
  download.file(paste(ASDAR_BOOK, "book2ed", fn, sep = "/"), fn)
  da <- paste(i, "bundle.zip", sep = "_")
  download.file(paste(ASDAR_BOOK, "bundles2ed", da, sep = "/"), da)
  unzip(da)
  source(fn, echo = TRUE)
}
```
