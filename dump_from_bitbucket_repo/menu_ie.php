<title>Applied Spatial Data Analysis with R</title>
<script type="text/javascript">

window.onload=montre;
function montre(id) {
var d = document.getElementById(id);
	for (var i = 1; i<=10; i++) {
		if (document.getElementById('smenu'+i)) {document.getElementById('smenu'+i).style.display='none';}
	}
if (d) {d.style.display='block'}
}
var windowHeight = window.innerHeight;
</script>


<style type="text/css">
<!-- 
/* CSS issu des tutoriels www.alsacreations.com/articles */
body {
margin:0;
padding: 0;
background: white;
font: 80% verdana, arial, sans-serif;
list-style-type:square;
}
body li{
}
smenu5 dl, dt, dd, ul, li {
margin: 0;
padding: 0;
list-style-type: none;
}
#menu dl, dt, dd, ul, li {
margin: 0;
padding: 0;
list-style-type: none;
}
#menu {
position: absolute;
top: 0;
left: 0;
z-index:100;
width: 100%;
}
#menu dl {
float: left;
width: 16%;
margin: 0 1px;
}
#menu dl table {
float: left;
margin: 0 1px;
background:#efb32f;
color:#FFFFFF;
a:color:#FFFFFF;
}
#menu dt {
cursor: pointer;
text-align: center;
font-weight: bold;
background:#ae2f44;
color:#FFFFFF;
a:color:#FFFFFF;
border: 0px solid gray;
}
#menu dd {
border: 0px solid gray;
}
#menu li {
text-align: center;
background:#efb32f;
}
#menu li a, #menu dt a {
color: #000;
text-decoration: none;
display: block;
height: 90%;
border: 0 none;
color:#FFFFFF;
}
#menu li a:hover, #menu dt a:hover {
background: #eee;
color:#FF9966;
}
#menu_table_chapter{
font:larger;
text-decoration:underline;
}

#site {
position: absolute;
z-index: 1;
top : 70px;
left : 10px;
color: #000;
background-color: #dddd;
padding: 5px;
border: 1px solid gray; 
}

/*REmove blue border around images*/
img{  
border-style: none;
}


-->
</style>
</head>
<body>

<div id="menu">
	<dl>
		<dt onMouseOver="javascript:montre();"><a href="http://www.asdar-book.org" title="return to the homepage">homepage</a></dt>
	</dl>
	<dl>	
		<dt onMouseOver="javascript:montre('smenu3');"><a class='menu' href="">Code Download 2nd edition</a></dt>

			<dd id="smenu3" onMouseOver="javascript:montre('smenu3');" onMouseOut="javascript:montre('');">
				<ul>
					<?
						for($i=0; $i<$nchapters2ed; $i++)
						{
							$st= "<li><a href=\"data2ed.php?chapter=$i\">Chapter " . ($i+1) . "</a></li>";
							echo("$st");
						}
					?>
				</ul>
			</dd>
	</dl>
	<dl>	
		<dt onMouseOver="javascript:montre('smenu4');"><a class='menu' href="">Errata 2nd edition</a></dt>
			<dd id="smenu4" onMouseOver="javascript:montre('smenu4');" onMouseOut="javascript:montre('');">
				<ul>
					<?
						$i=$nerrchapters2ed-1;
						$st= "<li><a href=\"errata2ed.php?errchapter=$i\">" . $errchapnames2ed[$i] . "</a></li>";
						echo("$st");
						
						for($i=0; $i<($nerrchapters2ed-1); $i++)
						{
							$st= "<li><a href=\"errata2ed.php?errchapter=$i\">" . $errchapnames2ed[$i] . "</a></li>";
							echo("$st");
					}
					?>
				</ul>
			</dd>
	</dl>
	
	
	<dl>	
		<dt onMouseOver="javascript:montre('smenu2');"><a class='menu' href="">Figures 1st edition</a></dt>
			<dd id="smenu2" onMouseOver="javascript:montre('smenu2');" onMouseOut="javascript:montre('');">
				
				<ul>
					<?
						for($i=0; $i<$nchapters; $i++)
						{
							$st= "<li><a href=\"code.php?chapter=$i&figure=-1\">Chapter " . ($i+1) ."</a>"; 
							echo("$st");
							echo("</li>");
						}	
					?>
				</ul>
			</dd>
		
	</dl>
	
	<dl>	
		<dt onMouseOver="javascript:montre('smenu3');"><a class='menu' href="">Code Download 1st edition</a></dt>

			<dd id="smenu3" onMouseOver="javascript:montre('smenu3');" onMouseOut="javascript:montre('');">
				<ul>
					<?
						for($i=0; $i<$nchapters; $i++)
						{
							$st= "<li><a href=\"data.php?chapter=$i\">Chapter " . ($i+1) . "</a></li>";
							echo("$st");
						}
					?>
				</ul>
			</dd>
	</dl>
	
	<dl>			
		<dt onMouseOver="javascript:montre('smenu1');"><a class='menu' href="">Data Sets Download 1st edition</a></dt>
			<dd id="smenu1" onMouseOver="javascript:montre('smenu1');" onMouseOut="javascript:montre('');">
				<ul>
											
						<?
						for($i=0; $i<$ndatasets; $i++)
						{
						
								$st= "<li><a href=\"datasets.php?dataset=$i\">" . $datasets[$i] . "</a><li>";
								echo("$st");
						}
						?>
				</ul>
			</dd>
	</dl>
	
	<dl>	
		<dt onMouseOver="javascript:montre('smenu5');"><a class='menu' href="">Additional Materials 1st edition</a></dt>

			<dd id="smenu5" onMouseOver="javascript:montre('smenu5');" onMouseOut="javascript:montre('');">
				<ul>
					<?
						
						//$st= "<li><a href=\"http://www.bias-project.org.uk/ASDARcourse\" target=\"_blank\">Course Material</a><br></li>";
						//echo("$st");
						
						$st= "<li><a href=\"courses.php\">Courses</a><br></li>";
						echo("$st");
	
	
						for($i=0; $i<$nexcodes; $i++)
						{
						
								$st= "<li><a href=\"exercises.php?excode=$i\">" . $excodenames[$i] . "</a><li>";
								echo("$st");
						}
					?>
				</ul>
			</dd>
	</dl>
	
	<dl>	
		<dt onMouseOver="javascript:montre('smenu4');"><a class='menu' href="">Errata 1st edition</a></dt>
			<dd id="smenu4" onMouseOver="javascript:montre('smenu4');" onMouseOut="javascript:montre('');">
				<ul>
					<?
						$i=$nerrchapters-1;
						$st= "<li><a href=\"errata.php?errchapter=$i\">" . $errchapnames[$i] . "</a></li>";
						echo("$st");
						
						for($i=0; $i<($nerrchapters-1); $i++)
						{
							$st= "<li><a href=\"errata.php?errchapter=$i\">" . $errchapnames[$i] . "</a></li>";
							echo("$st");
					}
					?>
				</ul>
			</dd>
	</dl>
	
</div>
<!-- this parts move everythings on the linked pages a bit down, otherwise the menu would overlap the titles of the page-->
<table height="3em">
<tr><td> 
</td></tr>
</table>
<?
include("submenu_ie.php");
?>

</body>
</html>
