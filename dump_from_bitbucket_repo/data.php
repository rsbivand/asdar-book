<?

include("header.htm");

include("book-vars.php");

include("generate_menu.php");


$chapter= $_GET["chapter"];
/*$figure= $_GET["figure"];*/

echo("\n<P><h2>Chapter ". ($chapter+1) . ": Code and download bundle</h2>\n");
echo("<P>\n");

echo("\n<P>\n<P>\n\n");
$Rfile = $chapters[$chapter].".R";
$Rfile_mod = $chapters[$chapter]."_mod.R";
$Rfile_bun = $chapters[$chapter]."_bundle.zip";
echo("<h3>R Code</h3><p>\nAll verbatim chunk codes used in the Chapter are available in <a href=\"book/$Rfile\">$Rfile</a>.");
echo(" This file contains the code used in the production of the book, and may be used for reference. It will not be changed, and may fall behind current releases of R and contributed packages.");
echo(" A simplified, modified, file <a href=\"book/$Rfile_mod\">$Rfile_mod</a> is also available, and will be kept up to date with reference to current R and contributed packages; it also will include/includes coloured figure code.<br><p>It is included in the chapter script and data sets <a href=\"bundles/$Rfile_bun\">bundle</a> as a zip archive. The easiest way to reproduce the analyses in this chapter is to download and unzip the bundle, start R in the directory in which the bundle was unpacked, and say:<pre>source(\"$Rfile_mod\", echo=TRUE)</pre><P>Note that some chapters involve very time-consuming analyses!<p>\n");

echo("<P>A <a href=\"https://stat.ethz.ch/pipermail/r-sig-geo/2009-January/004751.html\">posting</a> on the R-sig-geo list includes an R code snippet for downloading all of the bundles.<P>");

echo("<P>It is possible to download the code for all chapters using this short R snippet:<pre>ASDAR_BOOK <- \"http://www.asdar-book.org/book\"</pre><pre>chapters <- c(\"hello\", \"cm\", \"vis\", \"die\", \"cm2\",</pre><pre>\"csdacm\", \"sppa\", \"geos\", \"lat1\", \"lat2\", \"dismap\")</pre><pre># setwd(<folder>) # move to download folder</pre><pre>for (i in chapters) {</pre><pre>fn <- paste(i, \"mod.R\", sep=\"_\")</pre><pre>download.file(paste(ASDAR_BOOK, fn, sep = \"/\"), fn)</pre><pre>}</pre><pre>list.files()</pre><P>");

include("footer.htm");

?>


