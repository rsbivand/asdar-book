<?

include("header.htm");

include("book-vars.php");

include("generate_menu.php");

$excode= $_GET["excode"];


echo("\n<P><h2>Additional materials</h2>\n");
echo("<P>\n");

echo("\n<P>\n<P>\n\n");
$DS = $excodenames[$excode];
$DSfile = $excodes[$excode];
$DStxt = $excodes[$excode];

echo("<h3><a href=\"bundles/" . $DSfile . ".zip\">" . $DS . "</a></h3><p>\n");

include("datasets/" . $DStxt . ".txt");

include("footer.htm");

?>


