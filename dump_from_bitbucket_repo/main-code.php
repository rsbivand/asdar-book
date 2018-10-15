<?

/*This code displays code and figures*/

/*$chapter=1; /* Index starts at ZERO!! */

/*$chapter=$var1 = $_GET["chapter"];*/

echo("\n<P>\n<P>\n\n");
$Rfile = $chapters[$chapter].".R";
echo("<h3>R Code</h3><p>\nAll verbatim chunk codes used in the Chapter are available in <a href=\"book/$Rfile\">$Rfile</a><P>\n");

echo("<h3>Figures</h3><p>");


if($figure==-1) /*Display chapter info*/
{
/*	$offset=0;
	for($i=0; $i<($chapter-1);$i++)
	{
		$offset+=$nfigs[$i];
	}
*/

	echo("Click on figure to enlarge and see associated code<p>\n");


	 echo("<table>\n");

	for($i=0; $i<$nfigs[$chapter];$i++)
	{

		if( (($i)%5) ==0)
		{
			echo("<tr>\n");
		}


//$img="Fig " . ($fpref[$chapter]) . "-". ($i+1);
//$imgfile="Fig-" . ($fpref[$chapter]) . "-". ($i+1) . ".png";

//		$img=$figfiles[$offsetfigs[$i]+$i];
		$imgfile=$figfiles[$offsetfigs[$chapter]+$i];

//		echo $offsetfigs[$chapter]+$i ."<p>\n";

		echo("<td valign=\"top\"><pre>Figure " . ($chapter+1) . "." . ($i+1) . "</pre><p>\n");
		echo("<a href=\"code.php?chapter=$chapter&figure=$i\">");
		echo("<img width=120 src=\"book/$imgfile\">\n");
		echo("</a></td>");

		if(($i)%5 ==4)
		{
			echo("</tr>\n");
		}
	}

	 echo("</table>\n");

}
else /* Display a particular figure*/
{
		$img=$chunkfig[$offsetfigs[$chapter]+$figure];
		$imgfile=$figfiles[$offsetfigs[$chapter]+$figure];

		echo("<table>\n");



		echo("<td valign=\"top\">");
		echo("<a href=\"code.php?chapter=$chapter&figure=-1\">");
		echo("<img src=\"book/$imgfile\">\n");
		echo("</a>");
		echo("</td>");


	$offset=0;
        for($i=0; $i<($chapter);$i++)
        {
                $offset+=$nfigs[$i];
        }


		echo("<b>Figure ". ($chapter+1) . "." . ($figure+1) . "</b><p>\n");

	echo("Click on figure to go back to Chapter gallery.\n");
	echo("Please see preceeding chunks for values not set in this chunk (". ($chunkfig[$offset+$figure]). ").<br><p>\n");

		echo("<td valign=\"top\">\n");
		echo("<pre>\n");
		/*echo("book/" . ($chunkfig[$offset+$figure]). "");*/
		include("book/" . ($chunkfig[$offset+$figure]). "");
		echo("</pre>\n</td>\n");
		echo("</table>\n");

}

?>
