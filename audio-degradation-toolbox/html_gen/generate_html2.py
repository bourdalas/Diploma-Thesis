import sys

args = sys.argv

html_top = '''<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
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
<div style="padding:20px;">
<h2>Audio Degradation Toolbox -- Examples</h2>'''

html_original = \
'''<p style="background:#ea5;">
    <audio controls>
          <source src="../demoOutput/{1}_{0}.wav" type="audio/wav">
        Your browser does not support the audio element.
    </audio>
    {1}</p>'''

html_middle = \
'''<p>
    <audio controls>
          <source src="../demoOutput/{1}_{0}.wav" type="audio/wav">
        Your browser does not support the audio element.
    </audio>
    {1}</p>'''

html_bottom = '''</div></body>
      </html>'''


fileids = ["file{0}".format(i) for i in range(1,8)]

degids = ['Unit_01_addNoise', 
'Unit_02_addSound', 
'Unit_03_applyAliasing', 
'Unit_04_applyClipping', 
'Unit_05_applyDynamicRangeCompression', 
'Unit_06_applyHarmonicDistortion', 
'Unit_07_applyMp3Compression', 
'Unit_08_applySpeedup', 
'Unit_09_applyWowResampling',
'Degr_01_liveRecording', 
 'Degr_02_strongMp3Compression', 
 'Degr_03_vinylRecording', 
 'Degr_04_radioBroadcast', 
 'Degr_05_smartPhoneRecording', 
 'Degr_06_smartPhonePlayback',]

print html_top
for fid in fileids:
    print "<h3>{0}</h3>".format(fid)
    print html_original.format(fid, "00_Original")
    for did in degids:
        print html_middle.format(fid, did)
    print '''<br>
    <br>'''
print html_bottom