source("../uifunctions.R")
initialize('con',TRUE)

shinyUI(bootstrapPage(
	head(),
	navigation(),
	titlePanel("Contact"),
	beginPage(),
	HTML(	"
<b>Analysis and design:</b><br><a href='http://www.cbs.dtu.dk/staff/show-staff.php?id=1202'>Lasse Folkersen</a>.<br><br>
				
<b>Project lead:</b><br>Anders Malarstig.<br><br>

<b>Analysis group:</b><br>
                Maria Sabater-Lleal<br>
                Rona Strawbridge<br>
                Eric Fauman<br>
                Daniel Ziemek<br>
<br><br>
<b>Web-page styling:</b><br>Stefan Delport.<br><br>

	      				"),
	footer()
))












