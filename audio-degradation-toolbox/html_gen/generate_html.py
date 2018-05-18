import sys

args = sys.argv
wavfilestump = args[1].strip().split(".")[-2].split("/")[-1]
print wavfilestump
fileid = wavfilestump.split("_")[-1]
degid = "_".join(wavfilestump.split("_")[:-1])


html = '''<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
   "http://www.w3.org/TR/html4/loose.dtd">

<html lang="en">
<head>
	<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
	<title>untitled</title>
	<meta name="generator" content="TextMate http://macromates.com/">
	<meta name="author" content="Matthias Mauch">
	<!-- Date: 2013-08-20 -->
</head>
<body style = "background:#ddd;">
  <div>
    <div style="float:left;">
      <h3>Original ({0})</h3>
      <img src="../demoOutput/00_Original_{0}.png" width="800px"/><br>
      <div style="margin-left:auto;margin-right:auto;width:70%">
        <audio controls>
          <source src="../demoOutput/00_Original_{0}.wav" type="audio/wav">
        Your browser does not support the audio element.
        </audio>
      </div>
    </div>
    <div style="float:left;">
      <h3>{1} ({0})</h3>
      <img src="../demoOutput/{1}_{0}.png" width="800px"/><br>
      <div style="margin-left:auto;margin-right:auto;width:70%">
        <audio controls>
          <source src="../demoOutput/{1}_{0}.wav" type="audio/wav">
        Your browser does not support the audio element.
        </audio>
      </div>
    </div>
  </div>
</body>
</html>'''.format(fileid, degid)

f = open(degid+"_"+fileid+".html", "w")
print >>f, html
