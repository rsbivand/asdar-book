<?

//Here we define some 'global' variables to be used

$chapters = array("hello", "cm", "vis", "die", "cm2", "csdacm", "sppa", "geos", "lat1", "lat2", "dismap");

$nchapters = count($chapters);


$chapters2ed = array("hello", "cm", "vis", "die", "cm2", "std", "sppa", "geos", "lat", "dismap");

$nchapters2ed = count($chapters2ed);


$errchapters2ed = array("pref", "hello", "pt1", "cm", "vis", "die", "cm2", "std", "pt2", "sppa", "geos", "lat", "dismap", "after", "refs", "sind", "find", "all");

$errchapnames2ed = array("Preface", "Chapter 1", "Part 1", "Chapter 2", "Chapter 3", "Chapter 4", "Chapter 5", "Chapter 6", "Part 2", "Chapter 7", "Chapter 8", "Chapter 9", "Chapter 10", "Afterword", "References", "Subject Index", "Functions Index", "All errata");

$nerrchapters2ed = count($errchapters2ed);

$nerrschapters2ed = array(0, 1, 0, 2, 1, 7, 3, 0, 0, 2, 7, 4, 6, 1, 0, 0, 0, 0);


$fpref = array("hlo", "cm", "VIS", "die", "cm2",
   "csdacm", "sppa", "geos", "lat1", "lat2", "dismap");

//Number of figures in each chapter
//$nfigs = array(3, 9, 15, 5, 5, 3, 13, 15, 15, 8, 16);
$nfigs = array(3, 8, 15, 8, 5, 3, 13, 16, 15, 8, 18);
$offsetfigs = array(0, 3, 11, 26, 34, 39, 42, 55, 71, 86, 94);


//Vector of chunk codes associated to each figure
$chunkfig = array("hello-010.R", "hello-012.R", "hello-013.R", "cm-019.R",
"void.R", "cm-056.R", "void.R", "cm-072.R", "cm-081.R", "cm-091.R", "void.R",
"vis-012.R", "vis-014.R", "vis-016.R", "vis-018.R", "vis-020.R", "vis-022.R",
"vis-028.R", "vis-033.R", "vis-036.R", "vis-037.R", "vis-039.R", "vis-042.R",
"vis-045.R", "vis-058.R", "vis-059.R", "die-027.R", "void.R", "void.R",
"die-048.R", "die-049.R", "die-059.R", "die-061.R", "void.R", "cm2-008.R",
"cm2-011.R", "cm2-014.R", "cm2-021.R", "cm2-035.R", "csdacm-038.R",
"csdacm-044.R", "csdacm-069.R", "sppa-013.R", "sppa-016.R", "sppa-020.R",
"sppa-024.R", "sppa-026.R", "sppa-028.R", "sppa-031.R", "sppa-036.R",
"sppa-041.R", "sppa-056.R", "sppa-060.R", "sppa-078.R", "sppa-085.R",
"geos-011.R", "geos-020.R", "geos-022.R", "geos-025.R", "geos-026.R",
"geos-038.R", "geos-045.R", "geos-050.R", "void.R", "geos-056.R", "geos-061.R",
"geos-090.R", "geos-098.R", "geos-103.R", "geos-107.R", "geos-111.R",
"lat1-009.R", "lat1-014.R", "lat1-018.R", "lat1-023.R", "lat1-026.R",
"lat1-028.R", "lat1-029.R", "lat1-041.R", "lat1-049.R", "lat1-068.R",
"lat1-072.R", "lat1-079.R", "lat1-083.R", "lat1-085.R", "lat1-086.R",
"lat2-011.R", "lat2-016.R", "lat2-025.R", "lat2-041.R", "lat2-048.R",
"lat2-062.R", "lat2-075.R", "lat2-076.R", "dismap-010.R", "dismap-011.R",
"dismap-014.R", "dismap-017.R", "dismap-019.R", "dismap-020.R", "PG-model.txt",
"dismap-025.R", "dismap-026.R", "dismap-027.R", "BYM-model.txt",
"dismap-039.R", "dismap-041.R", "dismap-042.R", "dismap-043.R", "dismap-044.R",
"dismap-058.R", "dismap-063.R");


//Files associated to each figure
$figfiles=array("Fig-hlo-1.png", "Fig-hlo-2.png", "Fig-hlo-3.png",
"Fig-cm-1.png", "FignoR-cm-2.png", "Fig-cm-2.png", "FignoR-cm-4.png",
"Fig-cm-3.png", "Fig-cm-4.png", "Fig-cm-5.png", "FignoR-cm-8.png",
"Fig-VIS-1.png", "Fig-VIS-2.png", "Fig-VIS-3.png", "Fig-VIS-4.png",
"Fig-VIS-5.png", "Fig-VIS-6.png", "Fig-VIS-7.png", "Fig-VIS-8.png",
"Fig-VIS-9.png", "Fig-VIS-10.png", "Fig-VIS-11.png", "Fig-VIS-12.png",
"Fig-VIS-13.png", "Fig-VIS-14.png", "Fig-VIS-15.png", "Fig-die-1.png",
"FignoR-die-2.png", "FignoR-die-3.png", "Fig-die-2.png", "Fig-die-3.png",
"Fig-die-4.png", "Fig-die-5.png", "FignoR-die-8.png", "Fig-cm2-1.png",
"Fig-cm2-2.png", "Fig-cm2-3.png", "Fig-cm2-4.png", "Fig-cm2-5.png",
"Fig-csdacm-1.png", "Fig-csdacm-2.png", "Fig-csdacm-3.png", "Fig-sppa-1.png",
"Fig-sppa-2.png", "Fig-sppa-3.png", "Fig-sppa-4.png", "Fig-sppa-5.png",
"Fig-sppa-6.png", "Fig-sppa-7.png", "Fig-sppa-8.png", "Fig-sppa-9.png",
"Fig-sppa-10.png", "Fig-sppa-11.png", "Fig-sppa-12.png", "Fig-sppa-13.png",
"Fig-geos-1.png", "Fig-geos-2.png", "Fig-geos-3.png", "Fig-geos-4.png",
"Fig-geos-5.png", "Fig-geos-6.png", "Fig-geos-7.png", "Fig-geos-8.png",
"FignoR-geos-9.png", "Fig-geos-9.png", "Fig-geos-10.png", "Fig-geos-11.png",
"Fig-geos-12.png", "Fig-geos-13.png", "Fig-geos-14.png", "Fig-geos-15.png",
"Fig-lat1-1.png", "Fig-lat1-2.png", "Fig-lat1-3.png", "Fig-lat1-4.png",
"Fig-lat1-5.png", "Fig-lat1-6.png", "Fig-lat1-7.png", "Fig-lat1-8.png",
"Fig-lat1-9.png", "Fig-lat1-10.png", "Fig-lat1-11.png", "Fig-lat1-12.png",
"Fig-lat1-13.png", "Fig-lat1-14.png", "Fig-lat1-15.png", "Fig-lat2-1.png",
"Fig-lat2-2.png", "Fig-lat2-3.png", "Fig-lat2-4.png", "Fig-lat2-5.png",
"Fig-lat2-6.png", "Fig-lat2-7.png", "Fig-lat2-8.png", "Fig-dismap-1.png",
"Fig-dismap-2.png", "Fig-dismap-3.png", "Fig-dismap-4.png", "Fig-dismap-5.png",
"Fig-dismap-6.png", "FignoR-dismap-7.png", "Fig-dismap-7.png",
"Fig-dismap-8.png", "Fig-dismap-9.png", "FignoR-dismap-11.png",
"Fig-dismap-10.png", "Fig-dismap-11.png", "Fig-dismap-12.png",
"Fig-dismap-13.png", "Fig-dismap-14.png", "Fig-dismap-15.png",
"Fig-dismap-16.png");

$datasets=array("Auckland SRTM", "Auckland shoreline", 
"Biological cell centres", "Broad Street GRASS location", 
"Broad Street files", "California redwood trees", "Cars", "CRAN mirrors", 
"Japan shoreline", "Japanese black pine saplings", 
"Lansing Woods maple trees", "Loggerhead turtle", "Manitoulin Island", 
"Maunga Whau volcano", "Meuse bank", "New York leukemia", 
"North Carolina SIDS", "North Derbyshire asthma", "Scottish lip cancer", 
"Spearfish", "US 1999 SAT scores", "US Census 1990 Counties", 
"World volcano locations");


$datasetfiles=array("70042108.zip", "auckland_mapgen.dat", "internal", 
"snow_location.tgz", "snow_files.zip", "internal", "internal", 
"CRAN051001a.txt", "internal", "internal", "internal", "seamap105_mod.csv", 
"high.RData", "internal", "internal", "NY_data.zip", "internal", 
"north_derby_asthma.zip", "internal", "internal", "state.sat.data_mod.txt", 
"90mfips.txt", "data1964al.xy");

$datasettxt=array("auck1.txt", "auck2.txt", "cells.txt", "snow.txt", 
"snow.txt", "redwood.txt", "cars.txt", "cran.txt", "japan.txt", 
"pine.txt", "maple.txt", "turtle.txt", "gshhs.txt", "volcano.txt", 
"meuse.txt", "NY.txt", "NC-SIDS.txt", "NDa.txt", "lip.txt", 
"spearfish.txt", "sats.txt", "census.txt", "wvolcano.txt", );

$ndatasets=count($datasets);

$excodes=array("spearfish_GRASS", "snow_GRASS", "nb_GRASS", "nb_ARC", "dismap_WB");

$excodenames=array("GRASS (1): Ch. 4", "GRASS (2): Ch. 4", "GRASS: Ch. 9", "ArcGIS: Ch. 9", "WinBUGS: Ch. 11");

$nexcodes=count($excodes);

$errchapters = array("pref", "hello", "pt1", "cm", "vis", "die", "cm2", "csdacm", "pt2", "sppa", "geos", "lat1", "lat2", "dismap", "after", "refs", "sind", "find", "all");

$errchapnames = array("Preface", "Chapter 1", "Part 1", "Chapter 2", "Chapter 3", "Chapter 4", "Chapter 5", "Chapter 6", "Part 2", "Chapter 7", "Chapter 8", "Chapter 9", "Chapter 10", "Chapter 11", "Afterword", "References", "Subject Index", "Functions Index", "All errata");

$nerrchapters = count($errchapters);

$nerrschapters = array(0, 1, 0, 5, 4, 5, 5, 1, 0, 5, 3, 3, 1, 1, 2, 0, 0, 0, 0);

?>
