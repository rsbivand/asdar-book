
<ul id="nav">
<li><a class='menu' href="http://www.asdar-book.org">Homepage</a></li>
<li><a class='menu' href="">Code download 2nd edition</a>

<?
        echo("\t<ul>");
for($i=0; $i<$nchapters2ed; $i++)
{

//      echo("\t<ul>\n");
        $st= "\t\t<li><a class='menu' href=\"data2ed.php?chapter=$i\">Chapter " . ($i+1) . "</a>\n";
        echo("$st");

        echo("\t\t</li>\n");
}
        echo("\t</ul>");
?>


</li>
<li><a class='menu' href="">Errata 2nd edition</a>


<?
        echo("\t<ul>");
		
		$i=$nerrchapters2ed-1;
		$st= "\t\t<li><a class='menu' href=\"errata2ed.php?errchapter=$i\">" . $errchapnames2ed[$i] . "</a>\n";
		echo("$st");
						
		for($i=0; $i<$nerrchapters2ed-1; $i++)
		{
		
		//      echo("\t<ul>\n");
				$st= "\t\t<li><a class='menu' href=\"errata2ed.php?errchapter=$i\">" . $errchapnames2ed[$i] . "</a>\n";
				echo("$st");
		
				echo("\t\t</li>\n");
		}
        echo("\t</ul>");
?>


</li>
<li><a class='menu' href="">Figures 1st edition</a>

<?

//Creates the menu for the figures, etc.
//Assumes that the global variables have been LOADED.

$offset=0;


echo("<ul id=\"nav\">\n");

for($i=0; $i<$nchapters; $i++)
{

//	echo("\t<ul>\n");
 	$st= "\t\t<li><a class='menu' href=\"code.php?chapter=$i&figure=-1\">Chapter " . ($i+1) . "</a>\n"; 
        echo("$st");

        echo("\t\t\t<ul>\n");

        for($j=0; $j<$nfigs[$i]; $j++)
        {
                $st="\t\t\t<li><a class='menu' href=\"code.php?chapter=$i&figure=$j\">Figure ". ($i+1) . "." . "" . ($j+1) . " </a></li>\n";
		echo("$st");

        }

//        echo("\t\t\t</ul>\n");

        echo("\t\t</li>\n");
	echo("\t</ul>");
}


echo("</ul>");


?>

</li>
<li><a class='menu' href="">Code download 1st edition</a>

<?
        echo("\t<ul>");
for($i=0; $i<$nchapters; $i++)
{

//      echo("\t<ul>\n");
        $st= "\t\t<li><a class='menu' href=\"data.php?chapter=$i\">Chapter " . ($i+1) . "</a>\n";
        echo("$st");

        echo("\t\t</li>\n");
}
        echo("\t</ul>");
?>


</li>
<li><a class='menu' href="">Data Sets download 1st edition</a>

<?
        echo("\t<ul>");
for($i=0; $i<$ndatasets; $i++)
{

//      echo("\t<ul>\n");
        $st= "\t\t<li><a class='menu' href=\"datasets.php?dataset=$i\">" . $datasets[$i] . "</a>\n";
        echo("$st");

        echo("\t\t</li>\n");
}
        echo("\t</ul>");
?>


</li>
<li><a class='menu' href="">Additional materials 1st edition</a>

<?
        echo("\t<ul>");
		
		//$st= "\t\t<li><a class='menu' href=\"http://www.bias-project.org.uk/ASDARcourse\" target=\"_blank\">Course Material</a>\n";
		//echo("$st");
		//echo("\t\t</li>\n");
		
		$st= "\t\t<li><a class='menu' href=\"courses.php\">Courses</a>\n";
		echo("$st");
		echo("\t\t</li>\n");
	
		for($i=0; $i<$nexcodes; $i++)
		{
		
		//      echo("\t<ul>\n");
				$st= "\t\t<li><a class='menu' href=\"exercises.php?excode=$i\">" . $excodenames[$i] . "</a>\n";
				echo("$st");
		
				echo("\t\t</li>\n");
		}
        echo("\t</ul>");
?>


</li>
<li><a class='menu' href="">Errata 1st edition</a>


<?
        echo("\t<ul>");
		
		$i=$nerrchapters-1;
		$st= "\t\t<li><a class='menu' href=\"errata.php?errchapter=$i\">" . $errchapnames[$i] . "</a>\n";
		echo("$st");
						
		for($i=0; $i<$nerrchapters-1; $i++)
		{
		
		//      echo("\t<ul>\n");
				$st= "\t\t<li><a class='menu' href=\"errata.php?errchapter=$i\">" . $errchapnames[$i] . "</a>\n";
				echo("$st");
		
				echo("\t\t</li>\n");
		}
        echo("\t</ul>");
?>


</li>
</ul>
