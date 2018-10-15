<?

include("header.htm");

include("book-vars.php");

include("generate_menu.php");


$errchapter= $_GET["errchapter"];
$errchapternm = $errchapters2ed[$errchapter];
$nerrthischapter = $nerrschapters2ed[$errchapter];

echo("\n<P><h2>Errata</h2>\n");
echo("<P>\n");

echo("<h3>" . $errchapnames2ed[$errchapter] . "</h3><p>");

if($errchapternm == "all") {
  for($i=0; $i<($nerrchapters2ed-1); $i++) {
    if($nerrschapters2ed[$i] > 0) {
      $nm = $errchapters2ed[$i];
      echo("<br>");
      include("book2ed/" . $nm . "_errata.txt");
      echo("<br>");
      echo("</p>");
    }
  }
} else {
  if($nerrthischapter == 0) include("book/no_errata.txt");
  else include("book2ed/" . $errchapternm . "_errata.txt");
}

include("footer.htm");

?>


