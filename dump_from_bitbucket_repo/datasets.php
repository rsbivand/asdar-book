<?

include("header.htm");

include("book-vars.php");

include("generate_menu.php");


$dataset= $_GET["dataset"];
/*$figure= $_GET["figure"];*/

echo("\n<P><h2>Data Sets</h2>\n");
echo("<P>\n");

echo("\n<P>\n<P>\n\n");
$DS = $datasets[$dataset];
$DSfile = $datasetfiles[$dataset];
$DStxt = $datasettxt[$dataset];

if($DSfile == "internal") echo("<h3>" . $DS . "</h3><p>\n");
else echo("<h3><a href=\"datasets/" . $DSfile . "\">Click here to download the " . $DS . " dataset</a></h3><p>\n");

include("datasets/" . $DStxt);


include("footer.htm");

?>


