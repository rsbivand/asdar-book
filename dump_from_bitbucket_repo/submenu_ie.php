<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<title>Applied Spatial Data Analysis with R</title>

</script>


<style type="text/css">
<!-- 
/* CSS issu des tutoriels www.alsacreations.com/articles */
#link a{
background:#efb32f;
text-decoration: none;
display: inline;
height: 100%;
border: medium;
color:#FFFFFF;
}
#link a:hover{
background: #eee;
color:#FF9966;
}

-->
</style>
</head>


<body>
<table height="10em">
<tr><td> 
</td></tr>
</table>
<div id="link">
<?
//get the current chapter number
 ereg ("(code.php\?chapter=)([0-9]{1,2})", $_SERVER['REQUEST_URI'], &$regs);
 if($regs[2]<>NULL)
 {
 	//write the links for each chapter
	for($j=0; $j<$nfigs[$regs[2]]; $j++)
	{
		$st="<a href=\"code.php?chapter=$regs[2]&figure=$j\">Figure ". ($regs[2]+1) . "." . "" . ($j+1) . " </a>";
		echo("$st");
		if(($j+1)%10==0){echo("<br>");}
	}
 }
?>
</div>
</body>
</html>