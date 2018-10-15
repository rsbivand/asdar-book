<?

include("header.htm");

include("book-vars.php");

include("generate_menu.php");


$errchapter= $_GET["errchapter"];
$errchapternm = $errchapters[$errchapter];
$nerrthischapter = $nerrschapters[$errchapter];

echo("\n<P><h2>Errata</h2>\n");
echo("<P>\n");

echo("<h3>" . $errchapnames[$errchapter] . "</h3><p>");

if($errchapternm == "all") {
  for($i=0; $i<($nerrchapters-1); $i++) {
    if($nerrschapters[$i] > 0) {
      $nm = $errchapters[$i];
      echo("<br>");
      include("book/" . $nm . "_errata.txt");
      echo("<br>");
      echo("</p>");
    }
  }
} else {
  if($nerrthischapter == 0) include("book/no_errata.txt");
  else include("book/" . $errchapternm . "_errata.txt");
}

include("footer.htm");

?>


