<?

include("header.htm");

include("book-vars.php");

include("generate_menu.php");


$chapter= $_GET["chapter"];
$figure= $_GET["figure"];

echo("\n<P><h2>Chapter ". ($chapter+1) . ": Code and Figures</h2>\n");
echo("<P>\n");

include("main-code.php");

include("footer.htm");

?>


